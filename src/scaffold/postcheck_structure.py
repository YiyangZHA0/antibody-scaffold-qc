"""Compare extracted framework rows to the uploaded structure's original model."""
from collections import defaultdict
from pathlib import Path

from .audit_structure import read_audit_structure


def _finding(chain, region, code, status, message, evidence):
    return dict(chain=chain, region=region, code=code, status=status,
                message=message, evidence=evidence)


def _match_chain(candidates, summary, rows, namespace, is_pdb):
    """Match the source instance, including auth/label collisions in old exports."""
    key = summary['source_chain_id']
    authors = {r['source_chain_auth'] for r in rows if r.get('source_chain_auth')}
    labels = {r['source_chain_label'] for r in rows if r.get('source_chain_label')}
    matches = []
    for chain in candidates:
        if authors and authors != {chain.auth_id}:
            continue
        if not is_pdb and labels and labels != {chain.label_id}:
            continue
        # Gemmi can synthesize a PDB subchain name at conversion. Match that
        # recorded label only together with the recorded native author ID.
        label_ids = labels if is_pdb and authors else {chain.label_id}
        available = {'auth': {chain.auth_id}, 'label': label_ids}
        ids = available.get(namespace, {chain.auth_id} | label_ids)
        if key in ids and chain.sequence == summary.get('full_source_sequence'):
            matches.append(chain)
    return matches[0] if len(matches) == 1 else None


def _check_framework_rows(rows, chain, out_id):
    groups = defaultdict(list)
    uncertain = set(chain.uncertain_sequence_indices)
    absent = set(chain.known_missing_positions)
    checked = 0
    total = 0
    mapped = set()
    for row in rows:
        region = row.get('region')
        if region not in ('FR1', 'FR2', 'FR3', 'FR4'):
            continue
        total += 1
        value = row.get('source_seq_index')
        site = {'source_position': value, 'aa': row.get('aa'),
                'numbering_label': row.get('numbering_label'), 'imgt_label': row.get('imgt_label')}
        try:
            index = int(value) - 1
        except (TypeError, ValueError):
            index = -1
        if index < 0 or index >= len(chain.sequence) or row.get('aa') != chain.sequence[index]:
            groups[region, 'coordinate_mapping_uncertain'].append(site)
            continue
        mapped.add(index)
        residue = chain.residues.get(index)
        if index in uncertain or (residue is None and index not in absent):
            groups[region, 'coordinate_mapping_uncertain'].append(site)
            continue
        if residue is None:
            code = 'missing_residue_coordinates'
        elif (row.get('source_auth_seq_id') not in (None, '') and
              (str(row['source_auth_seq_id']) != str(residue.auth_seq_id) or
               (row.get('source_insertion_code') or '').strip() != residue.insertion_code)):
            groups[region, 'coordinate_mapping_uncertain'].append(site)
            continue
        else:
            site.update(auth_seq_id=residue.auth_seq_id, insertion_code=residue.insertion_code,
                        missing_positive_backbone=residue.missing_backbone,
                        zero_occupancy_atoms=residue.zero_occupancy_atoms)
            if not residue.has_any_positive_atoms:
                code = 'zero_occupancy_coordinates'
            elif residue.missing_backbone:
                code = 'missing_backbone_atoms'
            else:
                checked += 1
                continue
        checked += 1
        groups[region, code].append(site)

    messages = {
        'coordinate_mapping_uncertain': 'Coordinate correspondence cannot be established at these converted positions.',
        'missing_residue_coordinates': 'The source declaration contains these framework residues, but the selected model has no corresponding coordinate records.',
        'zero_occupancy_coordinates': 'The selected model has coordinate records for these residues, all with zero or lower occupancy.',
        'missing_backbone_atoms': 'These residues lack positive-occupancy coordinates for one or more backbone atoms (N, CA, C, O).',
    }
    findings = [
        _finding(out_id, region, code, 'not_assessable' if code == 'coordinate_mapping_uncertain' else 'observed',
                 messages[code], {'sites': sites, 'sequence_basis': chain.sequence_basis})
        for (region, code), sites in groups.items()
    ]
    if not chain.declared_sequence:
        findings.append(_finding(out_id, None, 'coordinate_only_sequence', 'not_assessable',
                                 'The file supplies coordinate records only. Residues absent from those records cannot be enumerated.',
                                 {'source_chain_auth': chain.auth_id}))
    status = 'assessed' if total and checked == total and chain.declared_sequence else 'partially_assessed'
    if not checked:
        status = 'not_assessable'
    details = {'status': status, 'source_chain_auth': chain.auth_id,
               'source_chain_label': chain.label_id, 'sequence_basis': chain.sequence_basis,
               'mapping_status': chain.mapping_status, 'original_structure_warnings': chain.warnings,
               'framework_rows': total, 'framework_rows_assessed': checked,
               'declared_missing_positions_outside_assigned_frameworks': [i+1 for i in sorted(absent-mapped)]}
    return findings, details


def postcheck_coordinates(source_path, metadata, rows):
    """Report coordinate evidence for each converted FR under the chosen scheme.

    Accept conversion metadata and residue_map JSONL rows. Findings use output
    chain IDs, region names from those rows and one-based source positions.
    The caller binds the supplied file to its recorded source digest.
    """
    path = Path(source_path) if source_path else None
    findings = []
    details = {'chains': {}}
    result = {'status': 'not_assessable', 'findings': findings, 'summary': details}
    if path is None or path.suffix.lower() not in {'.pdb', '.cif', '.mmcif'}:
        details['reason'] = 'no_source_structure'
        return result
    provenance = metadata.get('structure_provenance') or {}
    index = provenance.get('model_index', metadata.get('source_model_index'))
    namespace = provenance.get('chain_namespace', metadata.get('source_chain_namespace'))
    is_pdb = path.suffix.lower() == '.pdb'
    try:
        source = read_audit_structure(path, int(index or 1), namespace='auth' if is_pdb else 'label')
    except (OSError, RuntimeError, ValueError) as exc:
        details['reason'] = 'structure_read_failed'
        findings.append(_finding(None, None, 'structure_read_failed', 'not_assessable',
                                 str(exc), {'source': str(path)}))
        return result
    details.update(model_count=source.model_count, model_provenance=source.model_provenance)
    if index is None and source.model_count > 1:
        details['reason'] = 'conversion_model_unknown'
        findings.append(_finding(None, None, 'conversion_model_unknown', 'not_assessable',
                                 'The earlier conversion did not record which structure model was used.',
                                 {'model_count': source.model_count}))
        return result
    details['model_index'] = source.model_index
    per_chain = defaultdict(list)
    for row in rows:
        per_chain[row.get('output_chain_id')].append(row)
    summaries = metadata.get('chains', [])
    for summary in summaries:
        out_id = summary['output_chain_id']
        chain_rows = per_chain[out_id]
        candidates = source.chains.values()
        chain = _match_chain(candidates, summary, chain_rows, namespace, is_pdb)
        if chain is None and is_pdb and namespace == 'label' and not any(r.get('source_chain_auth') for r in chain_rows):
            # The earlier sequence-only conversion used Gemmi's synthetic PDB
            # subchain name. Reproduce only its ID lookup, then use the original
            # declaration and occupancy-aware parser for all coordinate evidence.
            from .structure import read_structure
            key = summary['source_chain_id']
            try:
                converted_chain = read_structure(path, source.model_index, [key], 'label').chains[key]
            except (RuntimeError, ValueError, KeyError):
                converted_chain = None
            if converted_chain is not None:
                candidate = source.chains.get(converted_chain.auth_id)
                if candidate and candidate.sequence == summary.get('full_source_sequence'):
                    chain = candidate
        if chain is None:
            details['chains'][out_id] = {'status': 'not_assessable'}
            findings.append(_finding(out_id, None, 'source_chain_not_matched', 'not_assessable',
                                     'The selected source chain and its full sequence could not be uniquely verified.',
                                     {'source_chain_id': summary['source_chain_id'], 'namespace': namespace}))
            continue
        chain_findings, chain_info = _check_framework_rows(chain_rows, chain, out_id)
        findings.extend(chain_findings)
        details['chains'][out_id] = chain_info
    statuses = {c['status'] for c in details['chains'].values()}
    if statuses == {'assessed'}:
        result['status'] = 'assessed'
    elif statuses - {'not_assessable'}:
        result['status'] = 'partially_assessed'
    return result

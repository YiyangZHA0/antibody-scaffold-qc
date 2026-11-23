"""Write BoltzGen sequence-design YAML without any structure dependency.

Antibody numbering labels describe the input domain. Numeric tokens in the
BoltzGen sequence instead describe the *new length* of an entire selected CDR.
"""
from __future__ import annotations

from pathlib import Path

import yaml

from .errors import ScaffoldError
from .output import dump_json, parse_loop_lengths, ranges, write_rows


def _sequence_spec(rows, cdrs, lengths, chain_id):
    """Keep every fixed residue; replace each complete selected loop once."""
    replacements = {}
    for cdr in cdrs:
        if not cdr['designed']:
            continue
        positions = cdr['indices']
        if not positions or positions != list(range(positions[0], positions[-1] + 1)):
            raise ScaffoldError(f"{chain_id}-{cdr['name']}: CDR must be a contiguous sequence interval")
        replacements[positions[0]] = (positions[-1], str(lengths[f"{chain_id}-{cdr['name']}"]))

    parts, fixed = [], []
    i = 1
    while i <= len(rows):
        if i in replacements:
            if fixed:
                parts.append(''.join(fixed))
                fixed = []
            last, token = replacements[i]
            parts.append(token)
            i = last + 1
        else:
            fixed.append(rows[i - 1]['aa'])
            i += 1
    if fixed:
        parts.append(''.join(fixed))
    # Commas are accepted by BoltzGen and keep adjacent length tokens separate.
    return ','.join(parts)


def write_sequence_scaffold(out, metadata, annotation, records, cdrs, roles,
                            loop_specs=None, name=None):
    """Write a complete, coordinate-free BoltzGen YAML and audited mapping.

    Scope matches the structural writer: one heavy variable domain or a VH/VL
    pair. Constants, tags and an explicitly split scFv linker are not emitted.
    ``roles`` uses the same ``heavy``/``light`` names as normalize_domains.
    """
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    mapping_rows, summaries = [], []
    seen_chains = set()
    for domain in annotation['domains']:
        if domain.chain_type not in ('H', 'K', 'L'):
            raise ScaffoldError(f'Unsupported antibody chain type: {domain.chain_type}')
        chain_id = 'H' if domain.chain_type == 'H' else 'L'
        role = 'heavy' if chain_id == 'H' else 'light'
        if chain_id in seen_chains:
            raise ScaffoldError(f'Multiple domains map to output chain {chain_id}; select one antibody')
        seen_chains.add(chain_id)
        if ''.join(row['aa'] for row in domain.rows) != domain.sequence:
            raise ScaffoldError(f'{domain.record_id}: residue mapping does not match variable sequence')
        if records[domain.record_id][domain.start:domain.end] != domain.sequence:
            raise ScaffoldError(f'{domain.record_id}: variable domain does not match source sequence')

        per_chain = []
        for index, row in enumerate(domain.rows, 1):
            mapped = {
                **row,
                'source_seq_index': row['input_seq_index'],
                'output_chain_id': chain_id,
                'boltzgen_res_index': index,
                'included': True,
                'designed': row['region'] in cdrs and (roles is None or role in roles),
                # No input residue has a fixed-coordinate structure group.
                'structure_group': 0,
                'has_coordinates': False,
            }
            mapping_rows.append(mapped)
            per_chain.append(mapped)

        loop_summaries = []
        for region in ('CDR1', 'CDR2', 'CDR3'):
            region_rows = [row for row in per_chain if row['region'] == region]
            loop_summaries.append({
                'name': region,
                'sequence': ''.join(row['aa'] for row in region_rows),
                'length': len(region_rows),
                'imgt_labels': [row['imgt_label'] for row in region_rows],
                'numbering_labels': [row['numbering_label'] for row in region_rows],
                'cdr_numbering_labels': [row['cdr_label'] for row in region_rows],
                'indices': [row['boltzgen_res_index'] for row in region_rows],
                'designed': any(row['designed'] for row in region_rows),
            })
        summaries.append({
            'role': role,
            'chain_type': domain.chain_type,
            'output_chain_id': chain_id,
            'source_chain_id': domain.record_id,
            'record_id': domain.record_id,
            'domain_id': domain.domain_id,
            'full_source_sequence': records[domain.record_id],
            'full_variable_sequence': domain.sequence,
            'selected_sequence': domain.sequence,
            'length': len(domain.sequence),
            'sequence_basis': 'numbered_input_sequence',
            'trimmed_terminal_residues': 0,
            'cdrs': loop_summaries,
        })

    if not any(row['designed'] for row in mapping_rows):
        raise ScaffoldError('No design residues selected')
    lengths = parse_loop_lengths(loop_specs, summaries)
    document = {'entities': []}
    for summary in summaries:
        chain_id = summary['output_chain_id']
        chain_rows = [row for row in mapping_rows if row['output_chain_id'] == chain_id]
        sequence = _sequence_spec(chain_rows, summary['cdrs'], lengths, chain_id)
        document['entities'].append({'protein': {'id': chain_id, 'sequence': sequence, 'msa': -1}})

    display_name = ' '.join((name or out.name).splitlines())
    head = [
        f'# {display_name}',
        '# YAML kind: complete BoltzGen sequence-design task (entities/protein)',
        f"# Source molecule: {metadata['molecule_type']}",
        f"# Included domains: {'one heavy variable domain (VHH scope)' if len(summaries) == 1 else 'VH + VL (Fv scope; constant domains omitted)'}",
        f"# Scheme preset: {metadata['scheme_preset']}",
        f"# Numbering scheme: {metadata['numbering_scheme']}",
        f"# CDR definition: {metadata['cdr_definition']}",
        f"# CDR boundary numbering: {metadata['cdr_numbering_scheme']}",
        f"# CDR boundary intervals: {metadata['cdr_ranges']}",
        f"# Annotation tool/version: {metadata['numbering_engine']} {metadata['numbering_version']}",
        '# Design policy: replace each selected complete CDR; preserve FR and unselected CDR amino acid identities',
        '# Coordinate policy: no input coordinates; framework sequence is fixed, framework geometry is not fixed',
        '# Sequence syntax: amino acid letters are fixed; integer or min..max tokens specify designed lengths',
        f'# Design lengths: {lengths}',
        '# Mapping index basis: original variable-domain 1-based sequence positions BEFORE CDR replacement',
        '# Output positions shift when CDR lengths change; antibody numbering labels remain input annotations',
        '# Mapping: residue_map.csv / residue_map.jsonl; no scaffold coordinate file is required',
        '# Official check may emit checked/check.cif with zero-coordinate placeholders; this is not a predicted structure',
        f"# Input file: {metadata['input_file']}",
        f"# Input SHA256: {metadata['input_sha256']}",
    ]
    for summary in summaries:
        chain_id = summary['output_chain_id']
        head.extend([
            f"# {summary['role']} chain: source={summary['source_chain_id']}; output_label={chain_id}; type={summary['chain_type']}",
            f"# {chain_id} full source sequence: {summary['full_source_sequence']}",
            f"# {chain_id} numbered variable sequence: {summary['full_variable_sequence']}",
            f"# {chain_id} included sequence: {summary['selected_sequence']}",
        ])
        for loop in summary['cdrs']:
            head.append(f"# {chain_id}-{loop['name']}: {metadata['numbering_scheme']}={','.join(loop['numbering_labels'])}; boundary_labels={','.join(loop['cdr_numbering_labels'])}; IMGT_reference={','.join(loop['imgt_labels'])}; res_index={ranges(loop['indices'])}; sequence={loop['sequence']}; length={loop['length']}; design={loop['designed']}")
        for region in ('FR1', 'FR2', 'FR3', 'FR4'):
            region_rows = [row for row in mapping_rows if row['output_chain_id'] == chain_id and row['region'] == region]
            head.append(f"# {chain_id}-{region}: {metadata['numbering_scheme']}={','.join(row['numbering_label'] for row in region_rows)}; IMGT_reference={','.join(row['imgt_label'] for row in region_rows)}; res_index={ranges(row['boltzgen_res_index'] for row in region_rows)}; sequence={''.join(row['aa'] for row in region_rows)}; fixed_sequence=True")
    for warning in metadata.get('warnings', []):
        head.append('# Note: ' + warning.replace('\n', ' '))

    text = '\n'.join(head) + '\n\n' + yaml.safe_dump(document, sort_keys=False, width=1000)
    (out / 'scaffold.yaml').write_text(text, encoding='utf-8')
    (out / 'check.yaml').write_text(text, encoding='utf-8')
    write_rows(out, mapping_rows)
    full = {
        **metadata,
        'name': display_name,
        'chains': summaries,
        'loop_lengths': lengths,
        'design_policy': 'replace_complete_selected_cdrs',
        'scaffold_kind': 'sequence',
        'has_structure': False,
        'yaml_kind': 'boltzgen_design',
        'scaffold_output_scope': 'vhh' if len(summaries) == 1 else 'fv',
        'validation_status': 'pending_boltzgen',
    }
    dump_json(out / 'metadata.json', full)
    return full

"""Expose every domain and possible selection without silently pairing copies."""
from hashlib import sha256
from time import perf_counter
import csv
import json
from .errors import ScaffoldError
from .numbering import Numberer, read_fasta
from .output import dump_json
from .structure import read_structure, domain_mapping


def inspect_input(args, numberer=None):
    started = perf_counter()
    if args.out and args.out.exists() and (not args.out.is_dir() or any(args.out.iterdir())):
        raise ScaffoldError(f'Inspection output directory is not empty: {args.out}')
    if args.ncpu < 1 or args.model_index < 1:
        raise ScaffoldError('--ncpu and --model-index must be positive')
    path = args.input.resolve()
    if path.suffix.lower() in ('.pdb', '.cif', '.mmcif'):
        source = read_structure(path, args.model_index, namespace=args.chain_namespace)
        records = {key: c.sequence for key, c in source.chains.items()}
    elif path.suffix.lower() in ('.fasta', '.fa', '.faa'):
        source = None
        records = read_fasta(path)
    else:
        raise ScaffoldError('inspect accepts PDB, CIF/mmCIF or FASTA')
    n = numberer or Numberer(args.engine, args.ncpu, args.scheme, args.cdr_definition)
    chain_rows, domains = [], []
    for rid, sequence in records.items():
        info = {'record_id': rid, 'auth_chain': source.chains[rid].auth_id if source else None,
                'label_chain': source.chains[rid].label_id if source else None, 'sequence_length': len(sequence), 'domains': []}
        try:
            found = n.domains(rid, sequence)
            for d in found:
                domain = {'domain_id': d.domain_id, 'chain_type': d.chain_type, 'role': 'VH' if d.chain_type == 'H' else 'VL',
                          'input_start': d.start + 1, 'input_end': d.end, 'variable_length': len(d.sequence),
                          'sequence': d.sequence, 'score': d.score}
                if source:
                    mapping = domain_mapping(d, source.chains[rid])
                    domain['positions_without_coordinates'] = [d.start + i + 1 for i, p in mapping.items() if p not in source.chains[rid].residues]
                info['domains'].append(domain)
                domains.append(d)
            info['status'] = 'antibody_domains_found' if found else 'no_antibody_domain_detected'
        except ScaffoldError as exc:
            info.update(status='annotation_error', error=str(exc))
        chain_rows.append(info)
    heavy = [d for d in domains if d.chain_type == 'H']
    light = [d for d in domains if d.chain_type in ('K', 'L')]
    candidates = []
    for h in heavy:
        for l in light:
            if h.record_id == l.record_id:
                options = {'type': 'scfv', 'split_scfv': True}
                kind, distance = 'scfv_candidate', None
                if source:
                    options.update(heavy_chain=h.record_id, light_chain=l.record_id)
                else:
                    options['records'] = h.record_id
            else:
                options, kind, distance = {'type': 'fab'}, 'fab_pair_candidate', None
                if source:
                    options.update(heavy_chain=h.record_id, light_chain=l.record_id)
                    ca = []
                    for d in (h, l):
                        chain = source.chains[d.record_id]
                        mapped = domain_mapping(d, chain)
                        points = [chain.residues[p].find_atom('CA', '*') for p in mapped.values() if p in chain.residues]
                        ca.append([atom.pos for atom in points if atom is not None])
                    distance = min((a.dist(b) for a in ca[0] for b in ca[1]), default=None)
                    if distance is not None and distance > 10:
                        continue
                else:
                    options['records'] = h.record_id + ',' + l.record_id
            candidates.append({'kind': kind, 'heavy_domain': h.domain_id, 'light_domain': l.domain_id,
                               'closest_ca_angstrom': distance, 'options': options,
                               'note': 'Candidate only; geometry/chain typing does not establish the biological pair. Full conversion still validates coordinates.'})
    if not light:
        for h in heavy:
            options = {'type': 'vhh', **({'vhh_chain': h.record_id} if source else {'records': h.record_id})}
            candidates.append({'kind': 'single_vh_candidate', 'heavy_domain': h.domain_id, 'options': options,
                               'note': 'Use as VHH only when its origin is known; a single VH does not prove VHH identity.'})
    for c in candidates:
        c['options']['scheme'] = n.policy.preset
        if args.cdr_definition:
            c['options']['cdr_definition'] = args.cdr_definition
        if source:
            c['options'].update(chain_namespace=args.chain_namespace, model_index=args.model_index)
    report = {'status': 'inspected', 'input': str(path), 'input_sha256': sha256(path.read_bytes()).hexdigest(),
              'has_coordinates': source is not None, 'chain_namespace': args.chain_namespace if source else None,
              'numbering_engine': n.engine, 'numbering_version': n.version, **n.policy.metadata(),
              'chains': chain_rows, 'heavy_domains': len(heavy), 'light_domains': len(light),
              'requires_explicit_selection': len(heavy) > 1 or len(light) > 1 or len(candidates) > 1,
              'candidates': candidates, 'elapsed_seconds': perf_counter() - started,
              'note': 'Inspection does not generate a scaffold or call BoltzGen. Multiple copies are never silently selected.'}
    if args.out:
        args.out.mkdir(parents=True, exist_ok=True)
        dump_json(args.out / 'inspection.json', report)
        with (args.out / 'chains.csv').open('w', newline='', encoding='utf-8-sig') as handle:
            fields = ['record_id', 'auth_chain', 'label_chain', 'sequence_length', 'status', 'chain_types', 'variable_spans', 'error']
            writer = csv.DictWriter(handle, fieldnames=fields, extrasaction='ignore')
            writer.writeheader()
            for c in chain_rows:
                writer.writerow({**c, 'chain_types': ','.join(d['chain_type'] for d in c['domains']),
                                 'variable_spans': ';'.join(f"{d['input_start']}..{d['input_end']}" for d in c['domains'])})
    return report


def format_inspection(report):
    lines = [f"Input: {report['input']}", f"Scheme: {report['scheme_preset']}; CDR definition: {report['cdr_definition']}",
             f"Chain namespace: {report['chain_namespace'] or 'FASTA record IDs'}", '',
             'Record/chain | auth | label | full length | domains (type:start..end) | status']
    for c in report['chains']:
        spans = '; '.join(f"{d['chain_type']}:{d['input_start']}..{d['input_end']}" for d in c['domains']) or '-'
        lines.append(f"{c['record_id']} | {c['auth_chain'] or '-'} | {c['label_chain'] or '-'} | {c['sequence_length']} | {spans} | {c['status']}")
        if c.get('error'):
            lines.append('  ' + c['error'])
    lines += ['', f"Candidate selections: {len(report['candidates'])}; explicit selection required: {report['requires_explicit_selection']}"]
    for i, c in enumerate(report['candidates'], 1):
        lines.append(f"  {i}. {c['kind']}: " + json.dumps(c['options'], ensure_ascii=False))
    lines += ['', 'Use candidate options as CLI --flags or per-row manifest columns. Confirm biological pairing and VHH origin.',
              'A Fab scaffold contains VH+VL (Fv scope). Inspection is not a structure-validity certificate.']
    return '\n'.join(lines)

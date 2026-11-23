from __future__ import annotations
import argparse
from datetime import datetime
from hashlib import sha256
from pathlib import Path
from time import perf_counter
import json
import os
import sys
from uuid import uuid4
from . import __version__
from .errors import ScaffoldError
from .numbering import Numberer, classify, read_fasta
from .structure import read_structure, normalize_domains
from .resolve import resolve_fasta_structure, select_structure_domains
from .output import write_annotation, write_scaffold, write_yaml_aliases, dump_json
from .sequence import write_sequence_scaffold
from .validation import validate_scaffold
from .compare import compare_outputs
from .schemes import PRESETS, DEFINITIONS, normalize_name, resolve_policy

REPO = Path(__file__).resolve().parents[2]
EXTENSIONS = {'.pdb': 'from-pdb', '.cif': 'from-cif', '.mmcif': 'from-cif',
              '.fa': 'from-fasta', '.faa': 'from-fasta', '.fasta': 'from-fasta'}


def moldir_default():
    return Path(os.environ.get('BOLTZGEN_MOLDIR', REPO / '.cache/mols.zip'))


def add_numbering_options(p, scheme_default='imgt'):
    p.add_argument('--scheme', type=normalize_name, choices=PRESETS, default=scheme_default,
                   help='Numbering/CDR preset; contact aliases maccallum. Batch requires a scheme here or in every manifest row.')
    p.add_argument('--cdr-definition', type=normalize_name, choices=DEFINITIONS, help='Override CDR definition; AHo defaults to IMGT CDRs')
    p.add_argument('--engine', choices=['anarcii', 'anarci'], default='anarcii')
    p.add_argument('--ncpu', type=int, default=2, help='CPU threads per numbering worker (default: 2)')


def add_conversion_options(p):
    p.add_argument('--mode', choices=['auto', 'sequence', 'structure'], default='auto',
                   help='auto: FASTA uses sequence, PDB/CIF use structure; explicit FASTA structure options select structure')
    p.add_argument('--type', choices=['auto', 'fab', 'vhh', 'scfv'], default='auto')
    p.add_argument('--records', help='Comma-separated exact FASTA record IDs')
    p.add_argument('--name', help='Scaffold name for the YAML comment header')
    p.add_argument('--cdr-length', action='append', default=[], help='Repeatable: H-CDR3=8..18; complete CDR length')
    group = p.add_mutually_exclusive_group()
    group.add_argument('--heavy-chain')
    group.add_argument('--vhh-chain')
    p.add_argument('--light-chain')
    p.add_argument('--chain-namespace', choices=['auth', 'label'], default='auth')
    p.add_argument('--model-index', type=int, default=1, help='1-based structural model index')
    p.add_argument('--cdrs', default='CDR1,CDR2,CDR3')
    p.add_argument('--design-chains', default='H,L', help='H, L or H,L; VHH uses H')
    p.add_argument('--split-scfv', action='store_true', help='Explicitly remove linker and extract two-chain Fv')
    p.add_argument('--moldir', type=Path, default=moldir_default())
    p.add_argument('--structure', type=Path, help='Exact matching coordinates for a FASTA input')
    p.add_argument('--pdb-id')
    p.add_argument('--search-rcsb', action='store_true', help='Opt in to submitting selected variable-domain sequences to RCSB')
    p.add_argument('--cache', type=Path, default=REPO / '.cache/structures')
    p.add_argument('--offline', action='store_true', help='Disable remote structure lookup; sequence mode needs no structure or network')


def parser():
    p = argparse.ArgumentParser(description='Single/batch antibody scaffold conversion with selectable numbering and CDR definitions')
    p.add_argument('--version', action='version', version=f'scaffold {__version__}')
    sub = p.add_subparsers(dest='command', required=True)
    from .audit_cli import add_audit_parser
    add_audit_parser(sub)
    post = sub.add_parser('postcheck', help='Check FR1-FR4 in an existing conversion folder using its recorded scheme')
    post.add_argument('directory', type=Path)
    post.add_argument('--source', type=Path, help='Relocated original upload; must match the recorded SHA256')
    for name in ('annotate', 'from-fasta', 'from-pdb', 'from-cif', 'convert'):
        s = sub.add_parser(name)
        s.add_argument('input', type=Path)
        s.add_argument('--out', type=Path, help='New/empty output directory; default: outputs/<input>_<timestamp>_<unique>')
        add_numbering_options(s)
        if name == 'annotate':
            s.add_argument('--type', choices=['auto', 'fab', 'vhh', 'scfv'], default='auto')
            s.add_argument('--records')
        else:
            add_conversion_options(s)
    s = sub.add_parser('batch', help='Convert one/many files, directories or a CSV/JSON job manifest')
    s.add_argument('inputs', nargs='*', type=Path)
    s.add_argument('--manifest', type=Path)
    s.add_argument('--out', type=Path, help='New/empty batch root; default: outputs/batch_<timestamp>_<unique>')
    s.add_argument('--jobs', type=int, default=1, help='Worker processes; default: 1, reuses numbering model')
    s.add_argument('--recursive', action='store_true', help='Recurse into input directories')
    s.add_argument('--extensions', help='Directory/glob filter, e.g. pdb,cif,fasta (default: all supported extensions)')
    add_numbering_options(s, scheme_default=None)
    add_conversion_options(s)
    s = sub.add_parser('inspect', help='List antibody domains and candidate chain selections without generating a scaffold')
    s.add_argument('input', type=Path)
    s.add_argument('--out', type=Path, help='Optional empty directory for inspection.json and chains.csv')
    s.add_argument('--format', choices=['table', 'json'], default='table')
    s.add_argument('--chain-namespace', choices=['auth', 'label'], default='auth')
    s.add_argument('--model-index', type=int, default=1)
    add_numbering_options(s)
    s = sub.add_parser('doctor', help='Check local packages and molecule dictionary')
    s.add_argument('--moldir', type=Path, default=moldir_default())
    s = sub.add_parser('validate')
    s.add_argument('directory', type=Path)
    s.add_argument('--moldir', type=Path, default=moldir_default())
    s = sub.add_parser('compare')
    s.add_argument('directories', nargs='+', type=Path)
    return p


def conversion_mode(args):
    """Choose output semantics explicitly; FASTA auto never discovers coordinates."""
    fasta = args.command == 'from-fasta'
    external_structure = bool(args.structure or args.pdb_id or args.search_rcsb)
    if external_structure and not fasta:
        raise ScaffoldError('--structure/--pdb-id/--search-rcsb only apply to FASTA inputs')
    if args.mode == 'sequence' and external_structure:
        raise ScaffoldError('--mode sequence cannot be combined with structure-source options')
    selected = args.mode
    if selected == 'auto':
        selected = 'structure' if not fasta or external_structure else 'sequence'
    if selected == 'sequence' and fasta and any((args.heavy_chain, args.light_chain, args.vhh_chain)):
        raise ScaffoldError('FASTA sequence mode selects records with --records; --heavy-chain/--light-chain/--vhh-chain select structure chains. H/L types are assigned by antibody numbering.')
    return selected


def _conversion(args, numberer, timings):
    stamp = perf_counter()
    def mark(name):
        nonlocal stamp
        now = perf_counter()
        timings[name] = now - stamp
        stamp = now
    if args.out.exists() and (not args.out.is_dir() or any(args.out.iterdir())):
        raise ScaffoldError(f'Output directory is not empty: {args.out}. Use a new directory to preserve previous results.')
    if not args.input.is_file():
        raise ScaffoldError(f'Input does not exist: {args.input}')
    if args.ncpu < 1:
        raise ScaffoldError('--ncpu must be positive')
    if args.command == 'convert':
        if args.input.suffix.lower() not in EXTENSIONS:
            raise ScaffoldError(f'Unsupported input extension: {args.input.suffix}')
        args.command = EXTENSIONS[args.input.suffix.lower()]
    if args.command != 'annotate':
        mode = conversion_mode(args)
        if args.vhh_chain and args.light_chain:
            raise ScaffoldError('--vhh-chain cannot be combined with --light-chain')
        if args.type == 'fab' and bool(args.heavy_chain) != bool(args.light_chain):
            raise ScaffoldError('For an explicitly selected Fab, supply both --heavy-chain and --light-chain')
        if args.type == 'vhh' and args.light_chain:
            raise ScaffoldError('--type vhh cannot include --light-chain')
        if args.model_index < 1:
            raise ScaffoldError('--model-index must be positive')
    args._owns_output = True
    policy = resolve_policy(args.scheme, args.cdr_definition)
    numberer = numberer or Numberer(args.engine, args.ncpu, args.scheme, args.cdr_definition)
    if numberer.policy != policy:
        raise ScaffoldError('Supplied Numberer policy differs from CLI --scheme/--cdr-definition')
    source = None
    if args.command in ('from-pdb', 'from-cif'):
        selected = [x for x in (args.heavy_chain, args.light_chain, args.vhh_chain) if x]
        source = read_structure(args.input, args.model_index, selected or None, args.chain_namespace)
        records = {key: chain.sequence for key, chain in source.chains.items()}
    else:
        records = read_fasta(args.input)
        if args.records:
            wanted = args.records.split(',')
            if set(wanted) - set(records):
                raise ScaffoldError(f'Unknown FASTA records {set(wanted) - set(records)}')
            records = {key: records[key] for key in wanted}
    mark('read_input')
    try:
        annotation = classify(records, numberer, args.type)
    except ScaffoldError as exc:
        raise ScaffoldError(f'{exc} Use scaffold inspect to list domains and choose chains/records.') from exc
    mark('number_and_classify')
    metadata = write_annotation(args.out, annotation, records, args.input)
    if source is not None:
        metadata.update(source_model_index=args.model_index, source_chain_namespace=args.chain_namespace)
    mark('write_annotation')
    if args.command == 'annotate':
        return {'status': 'annotated', 'molecule_type': annotation['molecule_type'], **policy.metadata(),
                'domains': [(d.record_id, d.chain_type, len(d.sequence)) for d in annotation['domains']], 'output': str(args.out.resolve())}
    cdrs = args.cdrs.split(',')
    if not cdrs or set(cdrs) - {'CDR1', 'CDR2', 'CDR3'}:
        raise ScaffoldError('--cdrs must be comma-separated CDR1/CDR2/CDR3')
    chains = args.design_chains.split(',')
    if not chains or set(chains) - {'H', 'L'}:
        raise ScaffoldError('--design-chains must be H, L or H,L')
    roles = [{'H': 'heavy', 'L': 'light'}[x] for x in chains]
    if annotation['molecule_type'] == 'scfv' and not args.split_scfv:
        raise ScaffoldError('scFv annotation written. Use --split-scfv to explicitly extract its VH/VL into a two-chain Fv scaffold.')
    chain_map = {k: v for k, v in {'heavy': args.heavy_chain or args.vhh_chain, 'light': args.light_chain}.items() if v}
    if mode == 'sequence':
        if source is not None:
            # Validate the requested heavy/light chain roles even without keeping coordinates.
            select_structure_domains(source, annotation['domains'], chain_map)
        if annotation['molecule_type'] == 'scfv':
            metadata['warnings'].append('Explicit --split-scfv: linker removed; output is two-chain Fv, not a continuous scFv')
        write_sequence_scaffold(args.out, metadata, annotation, records, cdrs, roles, args.cdr_length, args.name)
        write_yaml_aliases(args.out)
        mark('write_sequence_scaffold')
        verification = validate_scaffold(args.out, args.moldir)
        mark('official_validation')
        return {'status': 'scaffold_validated', 'molecule_type': annotation['molecule_type'], 'scaffold_kind': 'sequence',
                **policy.metadata(), 'output': str(args.out.resolve()), 'chains': verification['chains']}
    if source is None:
        source, selections, provenance = resolve_fasta_structure(annotation, records, args.cache, args.structure, args.pdb_id,
            args.search_rcsb, chain_map, args.model_index, args.chain_namespace, input_path=args.input, offline=args.offline)
    else:
        selections = select_structure_domains(source, annotation['domains'], chain_map)
        provenance = {'resolution_method': 'coordinates from input PDB/mmCIF'}
    provenance.update({'structure_file': str(source.path), 'structure_sha256': sha256(source.path.read_bytes()).hexdigest(),
                       'model_index': source.model_index, 'chain_namespace': args.chain_namespace, 'source_kind': args.command})
    mark('resolve_structure')
    if annotation['molecule_type'] == 'scfv':
        annotation['warnings'].append('Explicit --split-scfv: linker removed; output is two-chain Fv, not a continuous scFv')
    rows, summaries, warnings = normalize_domains(source, selections, args.out / 'scaffold.cif', cdrs, roles)
    mark('normalize_structure')
    write_scaffold(args.out, metadata, rows, summaries, warnings, provenance, args.cdr_length, args.name)
    write_yaml_aliases(args.out)
    mark('write_scaffold')
    verification = validate_scaffold(args.out, args.moldir)
    mark('official_validation')
    return {'status': 'scaffold_validated', 'molecule_type': annotation['molecule_type'], 'scaffold_kind': 'structure',
            **policy.metadata(), 'output': str(args.out.resolve()), 'chains': verification['chains']}


def execute(args, numberer=None):
    if args.command == 'postcheck':
        from .postcheck import write_postcheck
        return write_postcheck(args.directory, args.source)
    if args.command == 'audit':
        from .audit_cli import run_audit
        return run_audit(args, numberer)
    if args.command in ('batch', 'annotate', 'convert', 'from-fasta', 'from-pdb', 'from-cif') and args.out is None:
        label = 'batch' if args.command == 'batch' else args.input.stem
        # Avoid reserved Windows names and retain the user source name in metadata.
        import re
        label = re.sub(r'[^A-Za-z0-9_-]', '_', label)[:50] or 'input'
        args.out = Path.cwd() / 'outputs' / f'{label}_{datetime.now():%Y%m%d_%H%M%S}_{uuid4().hex[:8]}'
    if args.command == 'batch':
        from .batch import run_batch
        return run_batch(args)
    if args.command == 'inspect':
        from .inspection import inspect_input
        return inspect_input(args, numberer)
    if args.command == 'doctor':
        from .diagnostics import doctor
        return doctor(args.moldir)
    if args.command == 'compare':
        return compare_outputs(args.directories)
    if args.command == 'validate':
        return validate_scaffold(args.directory, args.moldir)
    started = perf_counter()
    timings = {}
    args._owns_output = False
    try:
        result = _conversion(args, numberer, timings)
        if result['status'] == 'scaffold_validated':
            from .postcheck import write_postcheck
            check_started = perf_counter()
            result['framework_check'] = write_postcheck(args.out)
            timings['framework_check'] = perf_counter() - check_started
    except Exception as exc:
        if args._owns_output and args.out.is_dir():
            dump_json(args.out / 'timing.json', {'status': 'failed', 'elapsed_seconds': perf_counter() - started, 'stages_seconds': timings})
            if not (args.out / 'validation.json').exists():
                dump_json(args.out / 'failure.json', {'status': 'failed', 'error': str(exc), 'scaffold_is_validated': False})
        raise
    result.update(elapsed_seconds=perf_counter() - started, stages_seconds=timings)
    dump_json(args.out / 'timing.json', {'status': result['status'], 'elapsed_seconds': result['elapsed_seconds'], 'stages_seconds': timings})
    return result


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        result = execute(args)
        if args.command == 'inspect' and args.format == 'table':
            from .inspection import format_inspection
            print(format_inspection(result))
        else:
            print(json.dumps(result, ensure_ascii=False, indent=2))
        return 0 if result.get('passed', True) else 1
    except KeyboardInterrupt:
        print('Interrupted; completed batch outputs and progress.json are preserved.', file=sys.stderr)
        return 130
    except (ScaffoldError, ValueError, OSError, ImportError, RuntimeError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())

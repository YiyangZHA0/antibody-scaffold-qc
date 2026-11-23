"""Deterministic job planning, isolated outputs and reusable process workers."""
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter
import csv
import glob
import json
import multiprocessing
import os
import re
import sys

from .errors import ScaffoldError
from .output import dump_json

OPTIONS = ('mode', 'scheme', 'cdr_definition', 'engine', 'ncpu', 'type', 'records', 'name', 'cdr_length',
           'heavy_chain', 'light_chain', 'vhh_chain', 'chain_namespace', 'model_index', 'cdrs',
           'design_chains', 'split_scfv', 'moldir', 'structure', 'pdb_id', 'search_rcsb', 'cache', 'offline')
BOOLEANS = {'split_scfv', 'search_rcsb', 'offline'}
PATHS = {'structure', 'moldir', 'cache'}
NUMBERS = {'ncpu', 'model_index'}
_NUMBERERS = {}


def boolean(value):
    if isinstance(value, bool):
        return value
    if str(value).lower() in ('true', '1', 'yes'):
        return True
    if str(value).lower() in ('false', '0', 'no'):
        return False
    raise ScaffoldError(f'Invalid boolean {value!r}; use true/false')


def job_id(value):
    if not isinstance(value, str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,79}', value) or value.endswith('.'):
        raise ScaffoldError(f'Invalid job id {value!r}; use 1..80 ASCII letters/digits/dot/underscore/hyphen, starting with a letter/digit')
    if re.fullmatch(r'CON|PRN|AUX|NUL|COM[0-9]|LPT[0-9]', value.split('.')[0], re.I):
        raise ScaffoldError(f'Reserved Windows job id: {value}')
    return value


def read_manifest(path):
    path = Path(path).resolve()
    if path.suffix.lower() == '.csv':
        with path.open(encoding='utf-8-sig', newline='') as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None or len(reader.fieldnames) != len(set(reader.fieldnames)):
                raise ScaffoldError('CSV manifest needs a header with unique column names')
            rows = list(reader)
    elif path.suffix.lower() == '.json':
        obj = json.loads(path.read_text(encoding='utf-8-sig'))
        if isinstance(obj, dict):
            if set(obj) != {'jobs'}:
                raise ScaffoldError('JSON manifest object must contain only a jobs list')
            rows = obj['jobs']
        else:
            rows = obj
    else:
        raise ScaffoldError('Manifest must be .csv or .json')
    if not isinstance(rows, list) or not rows:
        raise ScaffoldError('Manifest contains no jobs')
    for i, row in enumerate(rows, 1):
        if not isinstance(row, dict) or not isinstance(row.get('input'), str) or not row['input'].strip():
            raise ScaffoldError(f'Manifest row {i} needs a nonempty input path string')
        if row.get('id') is not None and not isinstance(row['id'], str):
            raise ScaffoldError(f'Manifest row {i}: id must be a string')
        unknown = set(row) - {'id', 'input', *OPTIONS}
        if unknown:
            raise ScaffoldError(f'Manifest row {i}: unknown columns {sorted(map(str, unknown))}')
    return rows


def expand_inputs(inputs, recursive=False, extensions=None):
    from .cli import EXTENSIONS
    allowed = {'.' + x.strip().lower().lstrip('.') for x in extensions.split(',')} if extensions else set(EXTENSIONS)
    if not allowed or allowed - set(EXTENSIONS):
        raise ScaffoldError(f'Invalid extension filter: {extensions}')
    found = []
    for raw in inputs:
        path = Path(raw)
        if path.is_dir():
            found.extend(p for p in (path.rglob('*') if recursive else path.iterdir()) if p.is_file() and p.suffix.lower() in allowed)
        elif path.is_file():
            if path.suffix.lower() not in EXTENSIONS:
                raise ScaffoldError(f'Unsupported input: {path}')
            found.append(path)
        elif glob.has_magic(str(path)):
            found.extend(Path(p) for p in glob.glob(str(path), recursive=recursive) if Path(p).is_file() and Path(p).suffix.lower() in allowed)
        else:
            raise ScaffoldError(f'Input does not exist: {path}')
    unique = {os.path.normcase(str(p.resolve())): p.resolve() for p in found}
    if not unique:
        raise ScaffoldError('No supported input files found')
    return sorted(unique.values(), key=lambda p: str(p).casefold())


def conversion_argv(job):
    argv = ['convert', job['input'], '--out', job['output']]
    for name, value in job['options'].items():
        if value is None or value == '':
            continue
        option = '--' + name.replace('_', '-')
        if name in BOOLEANS:
            if value:
                argv.append(option)
        elif name == 'cdr_length':
            for item in value:
                argv += [option, str(item)]
        else:
            argv += [option, str(value)]
    return argv


def plan_jobs(args):
    from .cli import parser, EXTENSIONS
    if args.jobs < 1 or args.ncpu < 1:
        raise ScaffoldError('--jobs and --ncpu must be positive')
    if bool(args.inputs) == bool(args.manifest):
        raise ScaffoldError('Supply input files/directories OR --manifest, but not both')
    if args.manifest and (args.recursive or args.extensions):
        raise ScaffoldError('--recursive/--extensions apply to directory/glob inputs, not manifests')
    root = args.out.resolve()
    if root.exists() and (not root.is_dir() or any(root.iterdir())):
        raise ScaffoldError(f'Output directory is not empty: {root}; choose a new batch directory')
    defaults = {name: getattr(args, name) for name in OPTIONS}
    defaults = {name: str(Path(value).resolve()) if name in PATHS and value else value for name, value in defaults.items()}
    rows = read_manifest(args.manifest) if args.manifest else [{'input': str(p)} for p in expand_inputs(args.inputs, args.recursive, args.extensions)]
    base = args.manifest.resolve().parent if args.manifest else Path.cwd()
    jobs, ids = [], set()
    for i, row in enumerate(rows, 1):
        path = Path(row['input']).expanduser()
        path = (base / path).resolve() if not path.is_absolute() else path.resolve()
        if path.suffix.lower() not in EXTENSIONS:
            raise ScaffoldError(f'Job {i}: unsupported input extension {path.suffix}')
        options = dict(defaults)
        for name in OPTIONS:
            value = row.get(name)
            if value is None or value == '':
                continue  # Empty CSV fields inherit command-line defaults.
            if name in BOOLEANS:
                value = boolean(value)
            elif name in PATHS:
                p = Path(value).expanduser()
                value = str((base / p).resolve() if not p.is_absolute() else p.resolve())
            elif name in NUMBERS:
                try:
                    integer = int(value)
                except (TypeError, ValueError) as exc:
                    raise ScaffoldError(f'Job {i}: {name} must be a positive integer') from exc
                if integer < 1 or str(integer) != str(value):
                    raise ScaffoldError(f'Job {i}: {name} must be a positive integer')
                value = integer
            elif name == 'cdr_length':
                value = value.split(';') if isinstance(value, str) else value
                if not isinstance(value, list) or not all(isinstance(x, str) for x in value):
                    raise ScaffoldError('cdr_length must be a JSON list or semicolon-separated CSV value')
            elif not isinstance(value, str):
                raise ScaffoldError(f'Job {i}: {name} must be a string')
            options[name] = value
        if not options['scheme']:
            raise ScaffoldError(f'Job {i}: specify --scheme or a scheme in the manifest')
        stem = re.sub(r'[^A-Za-z0-9_-]', '_', path.stem)[:45] or 'input'
        identifier = job_id(row.get('id') or f'{i:04d}_{stem}_{path.suffix[1:]}')
        if identifier.casefold() in ids:
            raise ScaffoldError(f'Duplicate job id (case insensitive): {identifier}')
        ids.add(identifier.casefold())
        output = (root / identifier).resolve()
        if output.parent != root:
            raise ScaffoldError('Job output escapes batch directory')
        job = {'id': identifier, 'input': str(path), 'output': str(output), 'options': options}
        try:
            checked = parser().parse_args(conversion_argv(job))
        except SystemExit as exc:
            raise ScaffoldError(f'Invalid conversion options in job {identifier}') from exc
        options['scheme'], options['cdr_definition'] = checked.scheme, checked.cdr_definition
        jobs.append(job)
    return jobs


def _run_job(job):
    from .cli import execute, parser
    from .numbering import Numberer
    from .schemes import resolve_policy
    started = perf_counter()
    row = {'id': job['id'], 'input': job['input'], 'output': job['output'], 'worker_pid': os.getpid(),
           'format': Path(job['input']).suffix.lower().lstrip('.'), 'scheme': job['options']['scheme']}
    try:
        args = parser().parse_args(conversion_argv(job))
        policy = resolve_policy(args.scheme, args.cdr_definition)
        key = (args.engine, args.ncpu, policy)
        if key not in _NUMBERERS:
            _NUMBERERS[key] = Numberer(args.engine, args.ncpu, args.scheme, args.cdr_definition)
        # PyTorch thread count is process-global, whereas models are cached per
        # policy/ncpu. Restore it when switching back to an already loaded model.
        cached = _NUMBERERS[key]
        if args.engine == 'anarcii' and cached._model is not None:
            import torch
            if torch.get_num_threads() != cached._model.ncpu:
                torch.set_num_threads(cached._model.ncpu)
        result = execute(args, cached)
        row.update(passed=True, **result)
    except Exception as exc:
        row.update(passed=False, status='failed', error=str(exc), error_type=type(exc).__name__)
        out = Path(job['output'])
        # Do not write into a directory independently created by another actor.
        if not out.exists():
            out.mkdir(parents=True)
            dump_json(out / 'failure.json', {'status': 'failed', 'error': str(exc), 'scaffold_is_validated': False})
    row['elapsed_seconds'] = perf_counter() - started
    return row


def _atomic_json(path, obj):
    temporary = path.with_suffix(path.suffix + '.tmp')
    dump_json(temporary, obj)
    temporary.replace(path)


def run_batch(args):
    started = perf_counter()
    jobs = plan_jobs(args)
    root = args.out.resolve()
    root.mkdir(parents=True, exist_ok=True)
    workers = min(args.jobs, len(jobs))
    dump_json(root / 'plan.json', {'created_utc': datetime.now(timezone.utc).isoformat(), 'jobs': jobs,
                                 'requested_workers': args.jobs, 'effective_workers': workers})
    results = {}
    def record(row):
        results[row['id']] = row
        ordered = [results[j['id']] for j in jobs if j['id'] in results]
        _atomic_json(root / 'progress.json', {'status': 'running', 'total': len(jobs), 'completed': len(results), 'results': ordered})
        print(f"[{len(results)}/{len(jobs)}] {'PASS' if row['passed'] else 'FAIL'} {row['id']} ({row['elapsed_seconds']:.2f}s)", file=sys.stderr, flush=True)
    _atomic_json(root / 'progress.json', {'status': 'running', 'total': len(jobs), 'completed': 0, 'results': []})
    if workers == 1:
        for job in jobs:
            record(_run_job(job))
    else:
        with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context('spawn')) as pool:
            futures = {pool.submit(_run_job, job): job for job in jobs}
            for future in as_completed(futures):
                job = futures[future]
                try:
                    row = future.result()
                except Exception as exc:
                    row = {'id': job['id'], 'input': job['input'], 'output': job['output'], 'passed': False,
                           'scheme': job['options']['scheme'], 'status': 'failed', 'error': str(exc), 'elapsed_seconds': 0}
                record(row)
    ordered = [results[j['id']] for j in jobs]
    elapsed = perf_counter() - started
    successful = sum(row['passed'] for row in ordered)
    summary = {'status': 'PASS' if successful == len(jobs) else 'PARTIAL_FAILURE', 'passed': successful == len(jobs),
               'total': len(jobs), 'succeeded': successful, 'failed': len(jobs) - successful, 'workers': workers,
               'elapsed_seconds': elapsed, 'throughput_files_per_minute': successful * 60 / elapsed,
               'mean_job_seconds': sum(r['elapsed_seconds'] for r in ordered) / len(jobs),
               'timing_note': 'Wall clock includes worker/model startup and official validation; parallel job durations overlap.',
               'results': ordered}
    dump_json(root / 'summary.json', summary)
    _atomic_json(root / 'progress.json', {'status': summary['status'], 'total': len(jobs), 'completed': len(jobs), 'results': ordered})
    fields = ['id', 'input', 'format', 'scheme', 'cdr_definition', 'passed', 'elapsed_seconds', 'worker_pid', 'output', 'error']
    with (root / 'summary.csv').open('w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(ordered)
    failed = [{'id': j['id'], 'input': j['input'], **j['options']} for j in jobs if not results[j['id']]['passed']]
    dump_json(root / 'failed_jobs.json', {'jobs': failed})
    return {k: v for k, v in summary.items() if k != 'results'} | {'output': str(root), 'summary_json': str(root / 'summary.json'),
            'summary_csv': str(root / 'summary.csv'), 'failed_manifest': str(root / 'failed_jobs.json')}

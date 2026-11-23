from pathlib import Path
import csv
import json
import pytest
from scaffold.batch import plan_jobs, read_manifest, expand_inputs, job_id, run_batch
from scaffold.cli import parser, execute, main
from scaffold.compare import compare_outputs
from scaffold.errors import ScaffoldError
from scaffold.numbering import Numberer

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'tests/data'


def test_directory_glob_recursion_extension_and_dedup(tmp_path):
    (tmp_path / 'a.PDB').write_text('fixture')
    (tmp_path / 'a.fasta').write_text('fixture')
    (tmp_path / 'notes.txt').write_text('not input')
    (tmp_path / 'sub').mkdir()
    (tmp_path / 'sub/b.cif').write_text('fixture')
    assert len(expand_inputs([tmp_path, tmp_path / 'a.PDB'])) == 2
    assert len(expand_inputs([tmp_path], recursive=True)) == 3
    assert expand_inputs([tmp_path], extensions='pdb') == [(tmp_path / 'a.PDB').resolve()]
    assert expand_inputs([tmp_path / '*.fasta']) == [(tmp_path / 'a.fasta').resolve()]


@pytest.mark.parametrize('identifier', ['../escape', '..', '/abs', 'A/B', 'C:\\x', 'CON', 'LPT1.txt', 'a.', ''])
def test_unsafe_ids_rejected(identifier):
    with pytest.raises(ScaffoldError):
        job_id(identifier)


def test_manifest_relative_paths_per_row_schemes_and_csv_types(tmp_path):
    path = tmp_path / 'jobs.csv'
    with path.open('w', newline='', encoding='utf-8-sig') as handle:
        writer = csv.DictWriter(handle, fieldnames=['id', 'input', 'scheme', 'structure', 'ncpu', 'offline', 'cdr_length'])
        writer.writeheader()
        writer.writerow({'id': 'first', 'input': 'a.fasta', 'scheme': 'Chothia', 'structure': 'a.cif', 'ncpu': 1, 'offline': 'false', 'cdr_length': 'H-CDR1=6;H-CDR3=8..18'})
        writer.writerow({'id': 'second', 'input': 'b.pdb', 'scheme': 'Contact', 'offline': 'true'})
    args = parser().parse_args(['batch', '--manifest', str(path), '--out', str(tmp_path / 'out'), '--offline'])
    jobs = plan_jobs(args)
    assert [j['options']['scheme'] for j in jobs] == ['chothia', 'maccallum']
    assert jobs[0]['input'] == str((tmp_path / 'a.fasta').resolve())
    assert jobs[0]['options']['structure'] == str((tmp_path / 'a.cif').resolve())
    assert jobs[0]['options']['offline'] is False
    assert jobs[0]['options']['ncpu'] == 1
    assert jobs[0]['options']['cdr_length'] == ['H-CDR1=6', 'H-CDR3=8..18']
    assert jobs[1]['options']['offline'] is True


@pytest.mark.parametrize('rows,match', [
    ([{'input': 'x.pdb'}], 'scheme'),
    ([{'input': 'x.pdb', 'scheme': 'imgt', 'typo': 1}], 'unknown'),
    ([{'id': 'A', 'input': 'x.pdb', 'scheme': 'imgt'}, {'id': 'a', 'input': 'x.pdb', 'scheme': 'imgt'}], 'Duplicate'),
    ([{'input': 'x.pdb', 'scheme': 'imgt', 'offline': 'maybe'}], 'boolean'),
    ([{'input': 'x.pdb', 'scheme': 'imgt', 'ncpu': '1.5'}], 'integer'),
])
def test_bad_manifest_fails_before_output_creation(tmp_path, rows, match):
    manifest = tmp_path / 'jobs.json'
    manifest.write_text(json.dumps({'jobs': rows}))
    out = tmp_path / 'out'
    args = parser().parse_args(['batch', '--manifest', str(manifest), '--out', str(out)])
    with pytest.raises(ScaffoldError, match=match):
        run_batch(args)
    assert not out.exists()


def test_one_file_is_a_batch_and_scheme_is_explicit(tmp_path):
    input_path = tmp_path / 'one.pdb'
    input_path.touch()
    args = parser().parse_args(['batch', str(input_path), '--out', str(tmp_path / 'out')])
    with pytest.raises(ScaffoldError, match='scheme'):
        plan_jobs(args)
    args.scheme = 'imgt'
    jobs = plan_jobs(args)
    assert len(jobs) == 1
    assert Path(jobs[0]['output']).parent == (tmp_path / 'out').resolve()


@pytest.mark.parametrize('row', [{'input': 123}, {'input': ['a.pdb']}, {'input': '  '}, {'input': 'a.pdb', 'id': 0}])
def test_manifest_path_and_id_types(tmp_path, row):
    path = tmp_path / 'invalid.json'
    path.write_text(json.dumps([row]))
    with pytest.raises(ScaffoldError, match='string'):
        read_manifest(path)


def test_cached_model_restores_its_cpu_thread_count(tmp_path, monkeypatch):
    import torch
    from types import SimpleNamespace
    import scaffold.batch as module
    import scaffold.cli as cli
    from scaffold.schemes import resolve_policy
    previous = torch.get_num_threads()
    model = SimpleNamespace(_model=SimpleNamespace(ncpu=2))
    monkeypatch.setattr(module, '_NUMBERERS', {('anarcii', 2, resolve_policy('imgt')): model})
    observed = []
    def fake_execute(args, numberer):
        observed.append(torch.get_num_threads())
        return {'status': 'scaffold_validated'}
    monkeypatch.setattr(cli, 'execute', fake_execute)
    job = {'id': 'cached', 'input': str(tmp_path / 'x.pdb'), 'output': str(tmp_path / 'out'), 'options': {'scheme': 'imgt', 'ncpu': 2}}
    try:
        torch.set_num_threads(1)
        result = module._run_job(job)
        assert result['passed'] and observed == [2]
    finally:
        torch.set_num_threads(previous)


def test_manifest_and_paths_cannot_be_mixed(tmp_path):
    args = parser().parse_args(['batch', 'x.pdb', '--manifest', 'a.csv', '--out', str(tmp_path)])
    with pytest.raises(ScaffoldError, match='not both'):
        plan_jobs(args)


def test_partial_failure_summary_order_retry_manifest_and_exit_code(tmp_path, monkeypatch):
    import scaffold.batch as module
    manifest = tmp_path / 'jobs.json'
    manifest.write_text(json.dumps([{'id': str(i), 'input': f'{i}.pdb', 'scheme': 'imgt'} for i in range(3)]))
    def fake(job):
        return {'id': job['id'], 'input': job['input'], 'output': job['output'], 'scheme': 'imgt', 'passed': job['id'] != '1', 'elapsed_seconds': .1}
    monkeypatch.setattr(module, '_run_job', fake)
    out = tmp_path / 'out'
    assert main(['batch', '--manifest', str(manifest), '--out', str(out)]) == 1
    summary = json.loads((out / 'summary.json').read_text())
    assert (summary['total'], summary['succeeded'], summary['failed']) == (3, 2, 1)
    assert [r['id'] for r in summary['results']] == ['0', '1', '2']
    failed = read_manifest(out / 'failed_jobs.json')
    assert len(failed) == 1 and failed[0]['id'] == '1'
    retry = parser().parse_args(['batch', '--manifest', str(out / 'failed_jobs.json'), '--out', str(tmp_path / 'retry')])
    assert len(plan_jobs(retry)) == 1
    original = (out / 'summary.json').read_bytes()
    with pytest.raises(ScaffoldError, match='not empty'):
        run_batch(parser().parse_args(['batch', '--manifest', str(manifest), '--out', str(out)]))
    assert (out / 'summary.json').read_bytes() == original


@pytest.fixture(scope='module')
def numberer():
    return Numberer()


def test_inspect_multiple_vhh_copies_and_namespace(tmp_path, numberer):
    out = tmp_path / 'inspect'
    report = execute(parser().parse_args(['inspect', str(DATA / '7xl0.pdb'), '--out', str(out)]), numberer)
    assert report['heavy_domains'] == 2 and report['light_domains'] == 0
    assert report['requires_explicit_selection']
    assert {c['options']['vhh_chain'] for c in report['candidates']} == {'A', 'B'}
    assert (out / 'chains.csv').exists() and not (out / 'scaffold.yaml').exists()
    assert any(c['auth_chain'] != c['label_chain'] for c in report['chains'])
    label_report = execute(parser().parse_args(['inspect', str(DATA / '7xl0.pdb'), '--chain-namespace', 'label']), numberer)
    assert {c['options']['vhh_chain'] for c in label_report['candidates']} == {c['label_chain'] for c in report['chains']}


def test_inspect_fab_with_antigen_and_missing_coordinates(numberer):
    report = execute(parser().parse_args(['inspect', str(DATA / '5j13.cif')]), numberer)
    assert report['heavy_domains'] == report['light_domains'] == 1
    assert any(c['status'] == 'no_antibody_domain_detected' for c in report['chains'])
    opts = report['candidates'][0]['options']
    assert opts['heavy_chain'] == 'C' and opts['light_chain'] == 'B'
    incomplete = execute(parser().parse_args(['inspect', str(DATA / '6wgb.cif')]), numberer)
    assert incomplete['requires_explicit_selection']
    assert any(d.get('positions_without_coordinates') for c in incomplete['chains'] for d in c['domains'])


@pytest.mark.parametrize('flags', [['--type', 'fab', '--heavy-chain', 'H'], ['--vhh-chain', 'B', '--light-chain', 'A']])
def test_conflicting_or_incomplete_chain_selection_is_clear(tmp_path, flags):
    with pytest.raises(ScaffoldError, match='both|cannot'):
        execute(parser().parse_args(['convert', str(DATA / '7eow.pdb'), '--out', str(tmp_path / 'out'), *flags]))


@pytest.mark.integration
def test_real_parallel_batch_all_schemes_formats_and_failure_isolation(tmp_path):
    rows = []
    for i, scheme in enumerate(['imgt', 'kabat', 'chothia', 'abm', 'maccallum', 'aho']):
        pid = '7eow' if i % 2 == 0 else '6cr1'
        fmt = ('pdb', 'cif', 'fasta')[i % 3]
        select = {'type': 'vhh', 'vhh_chain': 'B'} if pid == '7eow' else {'type': 'fab', 'heavy_chain': 'H', 'light_chain': 'L'}
        rows.append({'id': scheme, 'input': str(DATA / f'{pid}.{fmt}'), 'scheme': scheme, 'mode': 'structure', **select})
    rows.insert(1, {'id': 'missing', 'input': str(tmp_path / 'missing.pdb'), 'scheme': 'imgt'})
    manifest = tmp_path / 'mixed.json'
    manifest.write_text(json.dumps(rows))
    out = tmp_path / 'parallel'
    result = execute(parser().parse_args(['batch', '--manifest', str(manifest), '--out', str(out), '--jobs', '2', '--offline']))
    assert (result['succeeded'], result['failed'], result['workers']) == (6, 1, 2)
    assert not result['passed']
    summary = json.loads((out / 'summary.json').read_text())
    assert len({r['worker_pid'] for r in summary['results']}) == 2
    assert read_manifest(out / 'failed_jobs.json')[0]['id'] == 'missing'
    serial_out = tmp_path / 'serial'
    serial = execute(parser().parse_args(['batch', '--manifest', str(manifest), '--out', str(serial_out), '--jobs', '1', '--offline']))
    assert (serial['succeeded'], serial['failed']) == (6, 1)
    for i, scheme in enumerate(['imgt', 'kabat', 'chothia', 'abm', 'maccallum', 'aho']):
        pid = '7eow' if i % 2 == 0 else '6cr1'
        fmt = ('pdb', 'cif', 'fasta')[i % 3]
        baseline = serial_out / scheme
        assert compare_outputs([out / scheme, baseline])['passed']
        timing = json.loads((out / scheme / 'timing.json').read_text())
        assert timing['elapsed_seconds'] >= timing['stages_seconds']['official_validation'] > 0

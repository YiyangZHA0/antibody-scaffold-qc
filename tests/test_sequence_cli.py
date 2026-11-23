"""Standalone inputs, output routing and mixed-format CLI regressions."""
from pathlib import Path
from unittest.mock import patch
import json
import shutil

import pytest

from scaffold.cli import conversion_mode, execute, parser
from scaffold.batch import plan_jobs
from scaffold.compare import compare_outputs
from scaffold.errors import ScaffoldError
from scaffold.numbering import Numberer, read_fasta

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'tests/data'


@pytest.mark.parametrize('command,flags,expected', [
    ('from-fasta', [], 'sequence'),
    ('from-fasta', ['--offline'], 'sequence'),
    ('from-fasta', ['--structure', 'template.cif'], 'structure'),
    ('from-fasta', ['--mode', 'structure'], 'structure'),
    ('from-pdb', [], 'structure'),
    ('from-cif', [], 'structure'),
    ('from-cif', ['--mode', 'sequence'], 'sequence'),
])
def test_explicit_and_default_modes(command, flags, expected):
    args = parser().parse_args([command, 'input', *flags])
    assert conversion_mode(args) == expected


@pytest.mark.parametrize('command,flags', [
    ('from-fasta', ['--mode', 'sequence', '--structure', 'a.cif']),
    ('from-fasta', ['--heavy-chain', 'A']),
    ('from-pdb', ['--structure', 'a.cif']),
])
def test_conflicting_sequence_options_rejected(command, flags):
    with pytest.raises(ScaffoldError):
        conversion_mode(parser().parse_args([command, 'input', *flags]))


def test_batch_modes_can_be_selected_per_row(tmp_path):
    manifest = tmp_path / 'jobs.json'
    manifest.write_text(json.dumps([
        {'input': 'a.fasta', 'scheme': 'imgt'},
        {'input': 'b.pdb', 'scheme': 'chothia', 'mode': 'sequence'},
    ]))
    args = parser().parse_args(['batch', '--manifest', str(manifest), '--out', str(tmp_path / 'out')])
    jobs = plan_jobs(args)
    assert [j['options']['mode'] for j in jobs] == ['auto', 'sequence']


@pytest.fixture(scope='module')
def numberer():
    return Numberer()


@pytest.mark.integration
def test_fasta_alone_no_network_and_complete_cdr_replacement(tmp_path, numberer):
    fasta = tmp_path / 'uploaded.fasta'
    records = read_fasta(DATA / '5y9k.fasta')
    fasta.write_text(''.join(f'>custom_{i}\n{seq}\n' for i, seq in enumerate(records.values())))
    out = tmp_path / 'converted'
    with patch('requests.sessions.Session.request', side_effect=AssertionError('FASTA must not request coordinates')):
        result = execute(parser().parse_args([
            'convert', str(fasta), '--scheme', 'imgt', '--type', 'fab',
            '--cdrs', 'CDR3', '--design-chains', 'H', '--cdr-length', 'H-CDR3=12',
            '--out', str(out),
        ]), numberer)
    assert result['scaffold_kind'] == 'sequence'
    assert not (out / 'scaffold.cif').exists()
    assert (out / 'scaffold.yaml').read_bytes() == (out / 'design.yaml').read_bytes()
    assert [c['designed_residues'] for c in result['chains']] == [12, 0]
    assert all(c['fixed_sequence_matches'] for c in result['chains'])
    assert json.loads((out / 'validation.json').read_text(encoding='utf-8'))['cli_check_passed']


def test_automatic_output_directory_is_unique(tmp_path, numberer, monkeypatch):
    fasta = tmp_path / 'input.fasta'
    shutil.copy2(DATA / '5y9k.fasta', fasta)
    monkeypatch.chdir(tmp_path)
    first = execute(parser().parse_args(['annotate', str(fasta)]), numberer)
    second = execute(parser().parse_args(['annotate', str(fasta)]), numberer)
    assert first['output'] != second['output']
    for result in (first, second):
        out = Path(result['output'])
        assert out.parent == tmp_path / 'outputs'
        assert (out / 'annotation.json').is_file()


@pytest.mark.integration
def test_mixed_directory_batch_works_without_output_argument(tmp_path, monkeypatch):
    incoming = tmp_path / 'incoming'
    incoming.mkdir()
    for extension in ('pdb', 'cif', 'fasta'):
        shutil.copy2(DATA / f'5y9k.{extension}', incoming / f'upload_{extension}.{extension}')
    monkeypatch.chdir(tmp_path)
    result = execute(parser().parse_args(['batch', str(incoming), '--scheme', 'chothia', '--jobs', '2']))
    assert result['passed'] and result['succeeded'] == 3 and result['failed'] == 0
    # The returned root contains a separate directory for each independent input.
    batch_roots = list((tmp_path / 'outputs').iterdir())
    assert len(batch_roots) == 1
    summary = json.loads((batch_roots[0] / 'summary.json').read_text(encoding='utf-8'))
    rows = summary['results']
    assert sorted(r['scaffold_kind'] for r in rows) == ['sequence', 'structure', 'structure']
    assert all((Path(r['output']) / 'validation.json').is_file() for r in rows)
    for job in rows:
        directory = Path(job['output'])
        audit = json.loads((directory / 'framework_check.json').read_text(encoding='utf-8'))
        assert audit['cdr_definition'] == 'chothia' and len(audit['frameworks']) == 8
        assert audit['original_verified'] and (directory / 'framework_check.md').is_file()
        if job['scaffold_kind'] == 'sequence':
            assert audit['coordinates']['status'] == 'not_assessable'
        else:
            assert any(f['code'] == 'missing_residue_coordinates' for f in audit['findings'])
    comparison = compare_outputs([r['output'] for r in rows])
    # 5Y9K has three unresolved terminal FR residues. Sequence mode correctly
    # retains them; structural mode trims them. CDR design is still identical.
    assert comparison['checks']['cdr_design']
    assert not comparison['checks']['sequence_design']
    assert comparison['comparison_scope'] == 'sequence_design'
    assert comparison['coordinates_compared'] is False

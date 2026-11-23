"""Regressions for tracing converted framework regions to original coordinates."""
from pathlib import Path

import pytest

from scaffold.postcheck_structure import postcheck_coordinates
from scaffold.structure import normalize_domains, read_structure


def _pdb_case(tmp_path, declaration, observations):
    """Place CA-only coordinate records under unequivocal native author IDs."""
    lines = []
    if declaration:
        lines.append(f"SEQRES   1 H {len(declaration):4}  {' '.join(declaration)}")
    for atom_number, (aa, auth_pos, occ) in enumerate(observations, start=1):
        lines.append(f'ATOM  {atom_number:5d}  CA  {aa:3s} H{auth_pos:4d}    '
                     f'{float(atom_number):8.3f}{0.:8.3f}{0.:8.3f}{occ:6.2f}{20.:6.2f}           C')
    path = tmp_path / 'original.pdb'
    path.write_text('\n'.join(lines + ['END']) + '\n', encoding='ascii')
    return path


def _check_rows(path, sequence, regions):
    rows = [dict(source_seq_index=p, aa=aa, output_chain_id='H', region=r)
            for p, (aa, r) in enumerate(zip(sequence, regions), start=1)]
    metadata = {'chains': [dict(output_chain_id='H', source_chain_id='H', full_source_sequence=sequence)]}
    return postcheck_coordinates(path, metadata, rows)


def test_ambiguous_pdb_keeps_unknowns_distinct_from_zero_occupancy(tmp_path):
    path = _pdb_case(tmp_path, ['ALA', 'ALA', 'ALA', 'VAL', 'CYS', 'ASP'],
                     [('ALA', 101, 1.), ('ALA', 104, 1.), ('CYS', 105, 0.), ('ASP', 106, 1.)])
    report = _check_rows(path, 'AAAVCD', ['FR1']*3 + ['FR2', 'FR3', 'FR4'])
    findings = {f['code']: f for f in report['findings']}
    assert report['status'] == 'partially_assessed'
    expected = {
        'coordinate_mapping_uncertain': ('FR1', [1, 2, 3]),
        'missing_residue_coordinates': ('FR2', [4]),
        'zero_occupancy_coordinates': ('FR3', [5]),
        'missing_backbone_atoms': ('FR4', [6]),
    }
    for code, (region, positions) in expected.items():
        assert findings[code]['region'] == region
        assert [site['source_position'] for site in findings[code]['evidence']['sites']] == positions
    assert findings['missing_backbone_atoms']['evidence']['sites'][0]['missing_positive_backbone'] == ['N', 'C', 'O']


def test_atom_only_is_not_terminal_deletion_evidence(tmp_path):
    path = _pdb_case(tmp_path, None, [('SER', 98, 1.)])
    report = _check_rows(path, 'S', ['FR4'])
    codes = {f['code'] for f in report['findings']}
    assert report['status'] == 'partially_assessed'
    assert 'coordinate_only_sequence' in codes
    assert 'missing_residue_coordinates' not in codes


@pytest.mark.integration
def test_original_cif_label_selection_from_real_chothia_conversion(tmp_path):
    from scaffold.numbering import Numberer

    path = Path(__file__).parent / 'data' / '5y9k.cif'
    source = read_structure(path, selected_chains=['B'], namespace='label')
    heavy = source.chains['B']
    assert heavy.auth_id == 'H'
    domain, = Numberer(scheme='chothia').domains('B', heavy.sequence)
    rows, summaries, _ = normalize_domains(source, [(domain, 'B', 'H')], tmp_path / 'normalized.cif')
    metadata = {'chains': summaries, 'structure_provenance': {'model_index': 1, 'chain_namespace': 'label'}}
    report = postcheck_coordinates(path, metadata, rows)
    assert report['summary']['chains']['H']['source_chain_label'] == 'B'
    assert report['summary']['chains']['H']['source_chain_auth'] == 'H'
    assert report['status'] == 'assessed'
    finding = next(f for f in report['findings'] if f['code'] == 'missing_residue_coordinates')
    assert finding['region'] == 'FR1'
    assert finding['evidence']['sites'][0]['aa'] == 'Q'
    assert not rows[0]['included']


def test_source_parse_error_is_a_finding(tmp_path):
    report = _check_rows(tmp_path / 'missing.pdb', 'S', ['FR4'])
    assert report['findings'][0]['code'] == 'structure_read_failed'
    assert report['status'] == 'not_assessable'

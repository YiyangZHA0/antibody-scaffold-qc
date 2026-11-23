import json
from pathlib import Path

import pytest

from scaffold.errors import ScaffoldError
from scaffold.numbering import region_for
from scaffold.postcheck import _verified_source, _check_sequence_payload, _fasta_endpoints
from scaffold.schemes import region_for_label
from scaffold.terminal_references import compare_terminals


def _reference_scaffold():
    path = Path(__file__).parents[1] / 'src/scaffold/static/terminal_references.json'
    entries = json.loads(path.read_text(encoding='utf-8'))['entries']
    return next(r['sequence'] for r in entries if r.get('pdb') == '3EAK')


def _check(query, source=None):
    return compare_terminals({'chains': [dict(output_chain_id='H', chain_type='H',
                            full_variable_sequence=query, full_source_sequence=source or query)]})


def test_two_published_endpoint_versions():
    full = _reference_scaffold()
    assert _check(full)['findings'] == []
    result = _check(full[:-1])
    issue, = result['findings']
    assert (issue['code'], issue['region'], issue['status']) == ('known_reference_terminal_difference', 'FR4', 'reference_relative')
    assert issue['evidence']['residues'] == 'S'
    assert 'RFantibody' in result['comparisons'][0]['exact_versions'][0]
    assert _check(full[:-1], full)['findings'][0]['code'] == 'reference_residues_outside_numbered_domain'


def test_unrecognized_heavy_domain_with_similar_tail_only_gets_a_hint():
    full = _reference_scaffold()
    query = 'A' + full[1:-1]
    issue, = _check(query)['findings']
    assert (issue['code'], issue['status']) == ('terminal_motif_review', 'suspected')
    assert _check(query)['comparisons'][0]['status'] == 'no_domain_match'
    other_tail = query[:-4] + 'RRRR'
    assert not _check(other_tail)['findings']


def test_earlier_validation_cannot_silently_substitute_a_different_upload(tmp_path):
    path = tmp_path / 'wrong.fasta'
    path.write_text('>H\nAA\n')
    meta = dict(input_file=str(path), input_sha256='not-the-content-digest')
    findings = []
    assert _verified_source(meta, None, findings) is None
    assert findings[0]['code'] == 'original_unverified'
    with pytest.raises(ScaffoldError):
        _verified_source(meta, path, [])


def test_sequence_yaml_body_edit_detected_even_with_correct_comment(tmp_path):
    path = tmp_path / 'scaffold.yaml'
    path.write_text('# H-FR4: sequence=SS\nentities:\n- protein: {id: H, sequence: S}\n')
    meta = dict(scaffold_kind='sequence', chains=[{'output_chain_id': 'H'}])
    rows = [dict(output_chain_id='H', included=True, designed=False, region='FR4', aa='S')]*2
    assert _check_sequence_payload(path, meta, rows)[0]['code'] == 'sequence_yaml_payload_mismatch'


def test_imgt_fr4_last_possible_label_is_retained_by_conversion():
    assert region_for(129) == region_for_label(129, '', 'H', 'imgt') == 'FR4'
    with pytest.raises(ScaffoldError):
        region_for_label(130, '', 'H', 'imgt')


def test_endpoint_panel_never_uses_an_unverified_upload():
    result = _fasta_endpoints({'input_file': 'older.fasta'}, [], None, {}, [])
    assert result['reason'] == 'original_fasta_unavailable'
    assert _fasta_endpoints({'input_file': 'model.cif'}, [], None, {}, [])['reason'] == 'not_a_fasta_upload'


def test_endpoint_panel_rejects_wrong_record_even_with_matching_substring():
    chain = dict(record_id='chainA', domain_id='a', output_chain_id='H', full_variable_sequence='AAA')
    meta = dict(chains=[chain], full_input_sequences={'chainA': 'AAAQ'},
                domains=[dict(domain_id='a', start=0, end=3)])
    result = _fasta_endpoints(meta, [], {'chainA': 'AAAW'}, {}, [])
    assert [(a['side'], a['state']) for a in result['assessments']] == [
        ('N', 'indeterminate'), ('C', 'indeterminate')]


def test_endpoint_panel_verifies_numbered_copy_of_repeated_domain():
    chain = dict(record_id='a', domain_id='a_2', output_chain_id='H', full_variable_sequence='AAA')
    raw = {'a': 'AAAWAAA'}
    meta = dict(chains=[chain], full_input_sequences=raw, domains=[dict(domain_id='a_2', start=4, end=7)])
    rows = [dict(output_chain_id='H', input_seq_index=p) for p in (1, 2, 3)]
    result = _fasta_endpoints(meta, rows, raw, {}, [])
    assert result['status'] == 'not_assessed'
    assert all(a['evidence']['reason'] == 'source_domain_mapping_unverified' for a in result['assessments'])

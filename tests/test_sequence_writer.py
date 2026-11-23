"""Writer semantics on explicit annotated inputs, without numbering inference."""
import json

import pytest
import yaml

from scaffold.errors import ScaffoldError
from scaffold.numbering import Domain
from scaffold.sequence import write_sequence_scaffold


def annotated_domain(record='source_H', chain_type='H'):
    # Each region has distinct residues so boundary shifts are observable.
    fragments = [('FR1', 'AC'), ('CDR1', 'DE'), ('FR2', 'FG'),
                 ('CDR2', 'HI'), ('FR3', 'KL'), ('CDR3', 'MNP'), ('FR4', 'QR')]
    sequence = ''.join(fragment for _, fragment in fragments)
    rows = []
    for region, fragment in fragments:
        for aa in fragment:
            index = len(rows) + 1
            rows.append({'record_id': record, 'input_seq_index': index + 2,
                         'aa': aa, 'region': region, 'numbering_label': str(index),
                         'cdr_label': str(index), 'imgt_label': str(index)})
    return Domain(record, record + ':domain1', chain_type, 2, 2 + len(sequence), sequence, rows, None)


def write_fixture(out, domains=None, cdrs=('CDR1', 'CDR2', 'CDR3'), roles=None, lengths=None):
    domains = domains or [annotated_domain()]
    metadata = {'molecule_type': 'vhh' if len(domains) == 1 else 'fab',
                'scheme_preset': 'imgt', 'numbering_scheme': 'imgt',
                'cdr_definition': 'imgt', 'cdr_numbering_scheme': 'imgt',
                'cdr_ranges': {}, 'numbering_engine': 'explicit_test_fixture',
                'numbering_version': 'none', 'input_file': 'fixture.fasta',
                'input_sha256': 'fixture', 'warnings': []}
    records = {domain.record_id: 'SS' + domain.sequence + 'TT' for domain in domains}
    result = write_sequence_scaffold(out, metadata, {'domains': domains}, records,
                                    cdrs, roles, loop_specs=lengths)
    document = yaml.safe_load((out / 'scaffold.yaml').read_text())
    rows = [json.loads(line) for line in (out / 'residue_map.jsonl').read_text().splitlines()]
    return result, document, rows


def test_complete_loop_replacement_and_original_mapping(tmp_path):
    result, document, rows = write_fixture(tmp_path, lengths=['H-CDR3=5'])
    assert document == {'entities': [{'protein': {'id': 'H', 'sequence': 'AC,2,FG,2,KL,5,QR', 'msa': -1}}]}
    assert result['loop_lengths'] == {'H-CDR1': 2, 'H-CDR2': 2, 'H-CDR3': 5}
    assert [row['boltzgen_res_index'] for row in rows] == list(range(1, 16))
    # Length replacement must not change the input annotation or downstream FR4.
    assert result['chains'][0]['cdrs'][2]['indices'] == [11, 12, 13]
    assert [row['aa'] for row in rows if row['region'] == 'FR4'] == ['Q', 'R']
    assert all(row['input_seq_index'] == row['boltzgen_res_index'] + 2 for row in rows)
    assert sum(row['designed'] for row in rows) == 7


def test_unselected_loops_and_light_chain_remain_fixed(tmp_path):
    domains = [annotated_domain(), annotated_domain('source_L', 'K')]
    result, document, rows = write_fixture(tmp_path, domains, cdrs=['CDR3'], roles=['heavy'], lengths=['H-CDR3=4..8'])
    proteins = [entity['protein'] for entity in document['entities']]
    assert proteins[0]['sequence'] == 'ACDEFGHIKL,4..8,QR'
    assert proteins[1]['sequence'] == 'ACDEFGHIKLMNPQR'
    assert result['scaffold_output_scope'] == 'fv'
    assert [chain['role'] for chain in result['chains']] == ['heavy', 'light']
    assert all(not row['designed'] for row in rows if row['output_chain_id'] == 'L')


def test_sequence_yaml_is_self_contained_and_coordinate_free(tmp_path):
    result, document, rows = write_fixture(tmp_path)
    assert not (tmp_path / 'scaffold.cif').exists()
    assert (tmp_path / 'check.yaml').read_bytes() == (tmp_path / 'scaffold.yaml').read_bytes()
    assert result['has_structure'] is False and result['scaffold_kind'] == 'sequence'
    assert result['yaml_kind'] == 'boltzgen_design'
    assert all(row['included'] and not row['has_coordinates'] and row['structure_group'] == 0 for row in rows)
    assert all('file' not in entity for entity in document['entities'])
    assert 'framework geometry is not fixed' in (tmp_path / 'scaffold.yaml').read_text()


def test_cannot_override_an_unselected_loop(tmp_path):
    with pytest.raises(ScaffoldError, match='selected CDR'):
        write_fixture(tmp_path, cdrs=['CDR3'], lengths=['H-CDR1=8'])


def test_adjacent_loop_length_tokens_stay_separate(tmp_path):
    domain = annotated_domain()
    # Artificial stress case for the token serializer, not an antibody example.
    for row in domain.rows:
        if row['region'] == 'FR2':
            row['region'] = 'CDR2'
    _, document, _ = write_fixture(tmp_path, [domain])
    assert document['entities'][0]['protein']['sequence'] == 'AC,2,4,KL,3,QR'

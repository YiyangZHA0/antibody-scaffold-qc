"""Scheme changes must move actual residue masks, preserve AA mapping and stay valid."""
from pathlib import Path
from unittest.mock import patch
import json
import pytest
import yaml

from scaffold.cli import execute, parser
from scaffold.compare import compare_outputs
from scaffold.errors import ScaffoldError
from scaffold.numbering import Numberer, classify, read_fasta
from scaffold.schemes import PRESETS, resolve_policy, region_for_label, project_numbering

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'tests/data'


@pytest.mark.parametrize('scheme', ['chothia', 'abm', 'maccallum', 'chothia-consensus'])
def test_boundaries_against_independent_primary_source_fixture(scheme):
    # Extracted from Oxford's literal table and UCL HTML, not CDR_BOUNDS.
    reference = json.loads((ROOT / 'tests/reference/cdr_boundaries.json').read_text())
    for loop, (lo, hi) in reference['boundaries'][scheme].items():
        chain, i = loop[0], int(loop[1])
        for role in (['H'] if chain == 'H' else ['K', 'L']):
            assert region_for_label(lo - 1, '', role, scheme) == f'FR{i}'
            assert region_for_label(lo, '', role, scheme) == f'CDR{i}'
            assert region_for_label(hi, 'A', role, scheme) == f'CDR{i}'
            assert region_for_label(hi + 1, '', role, scheme) == f'FR{i + 1}'


def test_kabat_h1_exception_and_light_insertions():
    assert [region_for_label(35, ins, 'H', 'kabat') for ins in ('', 'A', 'B', 'C', 'D')] == ['CDR1'] * 3 + ['FR2'] * 2
    assert region_for_label(30, '', 'H', 'kabat') == 'FR1'
    assert region_for_label(31, '', 'H', 'kabat') == 'CDR1'
    assert region_for_label(34, 'C', 'K', 'kabat') == 'CDR1'
    assert region_for_label(100, 'K', 'H', 'kabat') == 'CDR3'


@pytest.mark.parametrize('command', ['convert', 'from-pdb', 'from-cif', 'from-fasta', 'annotate'])
def test_every_input_command_exposes_scheme_and_definition(command):
    args = parser().parse_args([command, 'input.fasta', '--out', 'out', '--scheme', 'AHo', '--cdr-definition', 'Contact'])
    assert args.scheme == 'aho' and args.cdr_definition == 'maccallum'


def test_policy_defaults_aliases_and_errors():
    assert resolve_policy().numbering_scheme == resolve_policy().cdr_definition == 'imgt'
    assert resolve_policy('AbM').numbering_scheme == 'chothia'
    assert resolve_policy('MacCallum') == resolve_policy('contact')
    assert resolve_policy('AHo').cdr_definition == 'imgt'
    assert resolve_policy('aho', 'abm').cdr_numbering_scheme == 'chothia'
    for scheme, definition in [('unknown', None), ('imgt', 'aho'), ('aho', 'unknown')]:
        with pytest.raises(ScaffoldError):
            resolve_policy(scheme, definition)


@pytest.fixture(scope='module')
def annotations():
    raw = Numberer()
    result = {}
    for preset in PRESETS:
        n = Numberer(scheme=preset)
        # Only these unit fixtures reuse raw IMGT model inference. Full matrix
        # conversions deliberately use independent models for each file format.
        n._number_one = raw._number_one
        result[preset] = {pid: classify(read_fasta(DATA / f'{pid}.fasta'), n, kind)
                          for pid, kind in [('6cr1', 'fab'), ('7eow', 'vhh')]}
    return result


def loops(domain):
    return [''.join(r['aa'] for r in domain.rows if r['region'] == f'CDR{i}') for i in (1, 2, 3)]


def test_known_adalimumab_and_caplacizumab_loop_sequences(annotations):
    # Literal loop sequences checked against native scheme labels on the public
    # 6CR1 and 7EOW sequences; these catch shifts at both sides of a loop.
    assert loops(annotations['kabat']['6cr1']['domains'][0]) == ['DYAMH', 'AITWNSGHIDYADSVEG', 'VSYLSTASSLDY']
    assert loops(annotations['chothia']['6cr1']['domains'][0]) == ['GFTFDDY', 'TWNSGH', 'VSYLSTASSLDY']
    assert loops(annotations['abm']['6cr1']['domains'][0]) == ['GFTFDDYAMH', 'AITWNSGHID', 'VSYLSTASSLDY']
    assert loops(annotations['maccallum']['6cr1']['domains'][0]) == ['DDYAMH', 'WVSAITWNSGHID', 'AKVSYLSTASSLD']
    assert loops(annotations['kabat']['7eow']['domains'][0]) == ['YNPMG', 'AISRTGGSTYYPDSVEG', 'AGVRAEDGRVRTLPSEYTF']


@pytest.mark.parametrize('scheme', PRESETS)
def test_real_projection_preserves_sequence_offsets_and_unique_labels(annotations, scheme):
    for pid, result in annotations[scheme].items():
        reference = annotations['imgt'][pid]
        assert [d.chain_type for d in result['domains']] == [d.chain_type for d in reference['domains']]
        for d, old in zip(result['domains'], reference['domains']):
            assert (d.start, d.end, d.sequence) == (old.start, old.end, old.sequence)
            assert [(r['aa'], r['input_seq_index'], r['imgt_label']) for r in d.rows] == [(r['aa'], r['input_seq_index'], r['imgt_label']) for r in old.rows]
            assert len(set(r['numbering_label'] for r in d.rows)) == len(d.rows)
            assert ''.join(r['aa'] for r in d.rows) == d.sequence
            if scheme == 'aho':
                assert [r['region'] for r in d.rows] == [r['region'] for r in old.rows]
                assert any(r['numbering_label'] != r['imgt_label'] for r in d.rows)
                assert any(r['numbering_position'] > 128 for r in d.rows)


def test_override_cdr_definition_is_projected_in_its_own_numbering_basis(annotations):
    seq = annotations['imgt']['6cr1']['domains'][0].sequence
    base = Numberer()
    for scheme, definition in [('kabat', 'chothia'), ('aho', 'abm')]:
        n = Numberer(scheme=scheme, cdr_definition=definition)
        n._number_one = base._number_one
        domain = n.domains('override', seq)[0]
        expected = annotations[definition]['6cr1']['domains'][0]
        assert loops(domain) == loops(expected)
        assert all(r['numbering_scheme'] == scheme for r in domain.rows)
        assert all(r['cdr_numbering_scheme'] == 'chothia' for r in domain.rows)


def test_projection_rejects_sequence_loss_and_duplicate_labels():
    positions = [((1, ' '), 'A'), ((2, ' '), 'C')]
    policy = resolve_policy('kabat')
    for bad, message in [([((1, ' '), 'A')], 'changed'), ([((1, ' '), 'A'), ((1, ' '), 'C')], 'duplicate')]:
        with patch('anarcii.output_data_processing.schemes.convert_number_scheme', return_value={'query': {'numbering': bad}}):
            with pytest.raises(ScaffoldError, match=message):
                project_numbering(positions, 'H', 1.0, policy)


def test_cli_rejects_injected_numberer_with_different_policy(tmp_path):
    args = parser().parse_args(['convert', str(DATA / '7eow.pdb'), '--out', str(tmp_path / 'out'), '--scheme', 'aho'])
    with pytest.raises(ScaffoldError, match='policy differs'):
        execute(args, Numberer())


@pytest.mark.integration
def test_aho_override_all_three_formats_and_compare_policy(tmp_path):
    outputs = []
    for fmt in ('pdb', 'cif', 'fasta'):
        out = tmp_path / fmt
        result = execute(parser().parse_args(['convert', str(DATA / f'7eow.{fmt}'), '--out', str(out),
                         '--mode', 'structure', '--type', 'vhh', '--vhh-chain', 'B', '--scheme', 'aho', '--cdr-definition', 'abm']))
        assert result['cdr_definition'] == 'abm' and result['numbering_scheme'] == 'aho'
        outputs.append(out)
    assert compare_outputs(outputs)['passed']
    text = (outputs[0] / 'scaffold.yaml').read_text()
    assert '# Numbering scheme: aho' in text and '# CDR definition: abm' in text
    assert yaml.safe_load(text)['path'] == 'scaffold.cif'
    other = tmp_path / 'abm'
    execute(parser().parse_args(['convert', str(DATA / '7eow.cif'), '--out', str(other), '--type', 'vhh',
                                '--vhh-chain', 'B', '--scheme', 'abm']))
    # Same CDR design body, different number-label maps and scheme policy.
    compared = compare_outputs([outputs[1], other])
    assert compared['checks']['yaml_body_text']
    assert not compared['checks']['selected_numbering_map']
    assert not compared['checks']['scheme_policy']
    assert not compared['passed']

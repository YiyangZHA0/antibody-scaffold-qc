"""Controlled endpoint decision tests, independent of the release/reference benchmark."""
import pytest
import json
from importlib.resources import files

from scaffold.germline_endpoints import assess_endpoints
from scaffold.terminal_references import compare_terminals


def _v():
    # Fixed alignment fixture. C23 anchors the held-back N-terminal segment.
    seq = list(('ACDEFGHIKLMNPQRSTVWY' * 6)[:104])
    seq[22] = 'C'
    seq[40] = 'W'
    seq[103] = 'C'
    return ''.join(seq)


def _entry(role, seg, seq, gene='test'):
    return dict(role=role, segment=seg, aligned=seq, species='fixture', genes=[gene])


def _run(*, tail='WGLGTSVTVSS', ncut=0, flank_n='', flank_c='', role='H',
         jrefs=('WGLGTSVTVSS',), bad=False, exact=False, extra_refs=(), context_only=False,
         substitutions=(), return_report=False):
    v = _v()
    query = list('Y' * len(v) if bad else v)
    for index, aa in substitutions:
        query[index-1] = aa
    observed = ''.join(query[ncut:]) + tail
    raw = flank_n + observed + flank_c
    meta = {'chains': [dict(output_chain_id='X', chain_type=role, full_variable_sequence=observed,
                            record_id='upload', full_source_sequence=raw)]}
    labels = list(range(1 + ncut, 105)) + list(range(118, 118 + len(tail)))
    rows = [dict(output_chain_id='X', input_seq_index=idx + len(flank_n), aa=aa, imgt_position=p)
            for idx, (aa, p) in enumerate(zip(observed, labels), 1)]
    refs = [_entry(role, 'V', v)] + [_entry(role, 'J', ref) for ref in jrefs] + list(extra_refs)
    known = dict(comparisons=[dict(chain='X', exact_versions=['exact fixture'])]) if exact else None
    report = assess_endpoints(meta, rows, catalog={'version':'fixture', 'entries':refs}, known_versions=known,
                              use_germline=not context_only)
    if return_report:
        return report
    return {r['side']:r for r in report['assessments']}


@pytest.mark.parametrize('ncut', [1, 2, 3])
def test_missing_n_residues_relative_to_independently_selected_v(ncut):
    n = _run(ncut=ncut)['N']
    assert n['state'] == 'suspected_short'
    assert n['evidence']['candidate_extensions'] == [_v()[:ncut]]


def test_internal_fr1_gap_cannot_be_called_an_n_terminal_deletion():
    v = _v()
    # Remove IMGT 12 while the true first residue remains in the upload.
    query = v[:11] + v[12:]
    labels = list(range(1, 12)) + list(range(13, 105))
    meta = {'chains': [dict(output_chain_id='H', chain_type='H', record_id='q',
                            full_variable_sequence=query, full_source_sequence=query)]}
    rows = [dict(output_chain_id='H', aa=aa, input_seq_index=i, imgt_position=label)
            for i, (aa, label) in enumerate(zip(query, labels), 1)]
    r = assess_endpoints(meta, rows, catalog={'entries': [_entry('H', 'V', v)]})
    n = r['assessments'][0]
    assert n['state'] == 'indeterminate' and n['evidence']['reason'] == 'terminal_register_unresolved'


def test_missing_residue_present_in_raw_flank():
    n = _run(ncut=1, flank_n=_v()[0])['N']
    assert n['state'] == 'boundary_review' and n['evidence']['matched_adjacent'] == [_v()[0]]
    assert n['evidence']['raw_domain_span'][0] == 2


def test_alternative_v_endpoints_win_over_duplicated_allele_votes():
    full = _entry('H', 'V', _v())
    short = _entry('H', 'V', '.' + _v()[1:])
    n = _run(ncut=1, extra_refs=[full]*12 + [short])['N']
    assert n['state'] == 'indeterminate'
    assert n['evidence']['candidate_extensions'] == ['', _v()[0]]
    assert n['evidence']['distinct_reference_proteins'] == 2


def test_j_consensus_can_raise_hint_but_a_vhh_single_s_cannot_establish_a_deletion():
    assert _run(tail='WGLGTSVTV')['C']['state'] == 'suspected_short'
    c = _run(tail='WGLGTSVTVS')['C']
    assert c['state'] == 'indeterminate' and c['evidence']['reason'] == 'known_single_and_double_serine_endpoints'
    c = _run(tail='WGLGTSVTVS', flank_c='S')['C']
    assert c['state'] == 'boundary_review'


def test_j_locus_and_ambiguity_handling():
    kappaj = 'FGQGTKVEIK'
    assert _run(role='K', tail=kappaj, jrefs=[kappaj], extra_refs=[_entry('H', 'J', kappaj + 'S')])['C']['state'] == 'compatible'
    light = _run(role='L', tail=kappaj[:-1], jrefs=[kappaj[:-1], kappaj])['C']
    assert light['state'] == 'indeterminate'


def test_weak_framework_support_abstains_even_with_exact_j_motif():
    r = _run(tail='WGLGTSVT', bad=True)
    for side in ('N','C'):
        assert r[side]['state'] == 'indeterminate'
        assert r[side]['evidence']['reason'] == 'insufficient_internal_framework_support'


def test_tag_and_constant_boundary_use_only_observed_context():
    assert _run(tail='WGLGTSVT', flank_c='HHHHHH')['C']['state'] == 'indeterminate'
    c = _run(tail='WGLGTSVT', flank_c='ASTKGPSV')['C']
    assert c['state'] == 'suspected_short' and c['evidence']['flank_type'] == 'constant_start_motif'
    assert _run(flank_n='M', context_only=True)['N']['state'] == 'indeterminate'
    assert _run(flank_c='S', context_only=True)['C']['state'] == 'boundary_review'


@pytest.mark.parametrize('role,j,appendage', [
    ('H', 'WGLGTSVTVSS', 'HHHHHH'), ('H', 'WGLGTSVTVSS', 'ASTKGPSV'),
    ('H', 'WGLGTSVTVSS', 'GGGGS' * 2), ('H', 'WGLGTSVTVSS', 'GGGGS' * 3),
    ('K', 'FGQGTKVEIK', 'RTVAAPSV'), ('L', 'FGGGTKLTVL', 'GQPKAAPSV'),
])
def test_clean_junction_vs_one_appendage_residue_included(role, j, appendage):
    cfg = dict(role=role, jrefs=[j])
    complete = _run(tail=j, flank_c=appendage, **cfg)['C']
    assert complete['state'] == 'compatible' and 'downstream_motif_starts' not in complete['evidence']
    borrowed = j[:-1] + appendage[0]
    report = _run(tail=borrowed, flank_c=appendage[1:], **cfg, return_report=True)
    c = report['assessments'][1]
    assert c['state'] == 'boundary_review'
    assert c['evidence']['reason'] == 'downstream_motif_overlaps_numbering'
    assert c['evidence']['downstream_motif_starts'][0]['source_start'] == 104 + len(j)
    assert c['evidence']['downstream_motif_starts'][0]['numbered_overlap'] == 1
    assert [f['code'] for f in report['findings']] == ['downstream_motif_overlaps_numbering']
    assert _run(tail=borrowed, **cfg)['C']['state'] == 'compatible'  # Unfused point variant.
    tagged_exact = _run(tail=borrowed, flank_c=appendage[1:], exact=True, **cfg)['C']
    assert tagged_exact['state'] == 'compatible' and tagged_exact['evidence']['downstream_motif_starts']


def test_shared_native_terminal_h_and_his_tag_retains_alternative_junctions():
    j = 'WGLGTSVTVSH'
    c = _run(tail=j, jrefs=[j], flank_c='HHHHHH')['C']
    assert c['state'] == 'indeterminate'
    assert c['evidence']['motif_also_starts_outside']
    assert {h['numbered_overlap'] for h in c['evidence']['downstream_motif_starts']} == {0,1}
    assert c['evidence']['native_residue_reference_endings'] == [j]
    c = _run(tail=j, jrefs=[j], flank_c='HHHHH')['C']
    assert c['state'] == 'indeterminate' and not c['evidence']['motif_also_starts_outside']
    assert c['evidence']['reason'] == 'ambiguous_downstream_junction'


def test_tandem_linker_requires_two_complete_repeats_and_respects_native_terminal_g():
    tandem = 'GGGGS' * 2
    j = 'WGLGTSVTVSG'
    native = _run(tail=j, jrefs=[j], flank_c=tandem)['C']
    assert native['state'] == 'compatible'
    assert 'downstream_motif_starts' not in native['evidence']
    ambiguous = _run(tail=j, jrefs=[j], flank_c=tandem[1:])['C']
    assert ambiguous['state'] == 'indeterminate'
    assert ambiguous['evidence']['native_residue_reference_endings'] == [j]
    assert ambiguous['evidence']['downstream_motif_starts'][0]['motif'] == tandem
    for partial in ('GGGS', 'GGGSGGGG', 'GGGSGGGGA'):
        c = _run(tail=j, flank_c=partial)['C']
        assert c['state'] == 'indeterminate' and 'downstream_motif_starts' not in c['evidence']


def test_junction_review_requires_an_actual_motif_and_honors_j_register_gates():
    j = 'WGLGTSVTVSH'
    for tail, flank in [(j, 'HHHH'), ('WGLGTSVTVSR', 'TVAAPS')]:  # Incomplete His and a K motif in H.
        c = _run(tail=tail, flank_c=flank)['C']
        assert c['state'] == 'indeterminate'
        assert c['evidence']['reason'] == 'same_length_terminal_residue_disagrees'
        assert 'downstream_motif_starts' not in c['evidence']
    c = _run(tail=j, flank_c='HHHHH', bad=True)['C']
    assert c['state'] == 'indeterminate' and c['evidence']['reason'] == 'insufficient_internal_framework_support'
    c = _run(tail=j, flank_c='HHHHH', substitutions=[(104,'A')])['C']
    assert c['state'] == 'indeterminate' and c['evidence']['reason'] == 'missing_conserved_cys104'


def test_same_register_terminal_substitution_remains_distinct_from_shortness():
    changed = 'WGLGTSVTVSA'
    refs = ['WGLGTSVTVSS']
    bare = _run(tail=changed, jrefs=refs)['C']
    assert bare['evidence']['reason'] == 'endpoint_length_only_residue_differs'
    fused = _run(tail=changed, jrefs=refs, flank_c='PQYT')['C']
    assert fused['state'] == 'indeterminate' and fused['evidence']['reference_endpoint_residues'] == ['S']
    supported = _run(tail=changed, jrefs=[*refs, changed], flank_c='PQYT')['C']
    assert supported['state'] == 'compatible'
    start_variant = [(1, 'T')]
    n = _run(substitutions=start_variant)['N']
    assert n['state'] == 'indeterminate' and n['evidence']['observed_endpoint_residue'] == 'T'
    n = _run(substitutions=start_variant, extra_refs=[_entry('H', 'V', 'T' + _v()[1:])])['N']
    assert n['state'] == 'compatible' and n['evidence']['reference_endpoint_residues'] == ['A','T']


def test_exact_construct_takes_precedence():
    c = _run(tail='WGLGTSVT', exact=True)['C']
    assert c['state'] == 'compatible' and c['evidence']['germline_only_state'] == 'suspected_short'


def test_noncanonical_anchors_and_numbered_start_methionine():
    r = _run(substitutions=[(104, 'A'), (1, 'M')], tail='WGLGTSVT')
    assert r['N']['state'] == 'boundary_review' and r['N']['evidence']['reason'] == 'ambiguous_numbered_start_methionine'
    assert r['C']['state'] == 'indeterminate' and r['C']['evidence']['reason'] == 'missing_conserved_cys104'


def test_identical_m_starts_from_substitution_or_initiator_yield_same_molecular_uncertainty():
    v = _v()
    engineered = 'M' + v[1:]
    clipped_plus_initiator = 'M' + v.removeprefix(v[0])
    assert engineered == clipped_plus_initiator
    for observed in (engineered, clipped_plus_initiator):
        r = _run(substitutions=list(enumerate(observed, 1)), return_report=True)
        finding, = r['findings']
        assert finding['code'] == 'ambiguous_numbered_start_methionine'
        n = r['assessments'][0]
        assert n['state'] == 'boundary_review' and len(n['evidence']['possible_explanations']) == 2
    cfg = dict(substitutions=[(1,'M')])
    assert _run(**cfg, extra_refs=[_entry('H','V',engineered)])['N']['state'] == 'indeterminate'
    assert _run(**cfg, exact=True)['N']['state'] == 'compatible'
    assert _run(**cfg, bad=True)['N']['evidence']['reason'] == 'insufficient_internal_framework_support'
    assert _run(substitutions=[(1,'M'),(23,'A')])['N']['evidence']['reason'] == 'missing_v_anchor'
    assert _run(flank_n='M')['N']['state'] == 'compatible'  # Numbered domain starts after M.


def test_raw_terminal_s_agrees_with_longer_known_version_even_when_127_domain_matches():
    refs = json.loads(files('scaffold').joinpath('static/terminal_references.json').read_text(encoding='utf-8'))
    short = next(r['sequence'] for r in refs['entries'] if 'RFantibody' in r['name'])
    raw = short + 'S'
    meta = {'chains': [dict(output_chain_id='H', chain_type='H', record_id='query',
                            full_variable_sequence=short, full_source_sequence=raw)]}
    # No germline evidence is needed to resolve this construct-version conflict.
    rows = [dict(output_chain_id='H', aa=aa, imgt_position=i, input_seq_index=i)
            for i, aa in enumerate(short, 1)]
    known = compare_terminals(meta)
    result = assess_endpoints(meta, rows, catalog={'entries': []}, known_versions=known)
    n, c = result['assessments']
    assert n['state'] == 'compatible'
    assert c['state'] == 'boundary_review' and c['evidence']['reason'] == 'source_supports_longer_catalog_boundary'
    assert c['evidence']['longer_versions_supported_by_flank'][0]['residues'] == 'S'
    assert c['evidence']['exact_catalog_versions']
    assert not result['findings']  # compare_terminals already provides the specific finding.

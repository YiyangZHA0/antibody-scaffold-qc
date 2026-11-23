"""Screen antibody variable-domain endpoints against germline *protein* references.

IMGT coordinates here align the evidence across numbering schemes. They do not
replace the extraction's chosen region definition or specify a required length.
"""
from collections import defaultdict
from functools import lru_cache
from importlib.resources import files
import json


# Fixed before evaluation. Identity thresholds are screening rules, not probabilities.
_FRAMEWORK = (*range(11, 27), *range(39, 56), *range(66, 105))
_MIN_FRAMEWORK_OVERLAP = 48
_MIN_FRAMEWORK_IDENTITY = 0.72
_NEAR_IDENTITY = 0.03
_CONSTANT_PREFIXES = {
    'H': ('ASTKG', 'ASTTT', 'AKTTP', 'ASTTS'),
    'K': ('RTVAA', 'RTVAP'),
    'L': ('GQPKA', 'GQPKS'),
}


@lru_cache(maxsize=1)
def _catalog():
    return json.loads(files('scaffold').joinpath('static/germline_endpoints.json').read_text(encoding='utf-8'))


def _flank_type(side, flank, role):
    if not flank:
        return 'record_edge'
    if side == 'N' and flank == 'M':
        return 'possible_initiator_M'
    if side == 'C' and flank.startswith(_CONSTANT_PREFIXES.get(role, ())):
        return 'constant_start_motif'
    if flank.startswith(('HHHHH', 'GGGGS')) or flank.endswith(('HHHHH', 'GGGGS')):
        return 'tag_or_linker_motif'
    if len(flank) <= 3:
        return 'short_unassigned_flank'
    return 'unassigned_flank'


def _downstream_crossing(seq, flank, role, end):
    """Look for a complete downstream motif starting just inside the numbered end."""
    if not flank:
        return None
    motifs = ('HHHHHH', 'GGGGS' * 2, *_CONSTANT_PREFIXES.get(role, ()))
    adjacent = flank[:max(map(len, motifs))]
    starts = []
    for overlap in range(min(3, len(seq)) + 1):
        suffix = seq[len(seq) - overlap:] + adjacent
        for motif in motifs:
            if suffix.startswith(motif):
                starts.append(dict(motif=motif, numbered_overlap=overlap, source_start=end + 1 - overlap))
    if any(hit['numbered_overlap'] for hit in starts):
        return starts
    return None


def _review_crossing(hits, seq, positions, opts, state, evidence, use_germline):
    # A native terminal H followed by His6 can share the same observations as a
    # numbered tag residue. Keep both possible junctions and check independent J alleles.
    at_junction = any(hit['numbered_overlap'] == 0 for hit in hits)
    covered = {hit['numbered_overlap'] for hit in hits if hit['numbered_overlap']}
    native_support = []
    if opts:
        anchor = positions.index((118, ''))
        tail = seq[anchor:]
        native_support = sorted({r['aligned'] for _, r in opts
                                 if len(r['aligned']) == len(tail)
                                 and any(r['aligned'][-n:] == tail[-n:] for n in covered)})
    evidence.update(downstream_motif_starts=hits, motif_also_starts_outside=at_junction,
                    native_residue_reference_endings=native_support,
                    pre_junction_state=state, pre_junction_reason=evidence['reason'])
    if use_germline and not opts:
        # Preserve the applicability/anchor gate even if a fusion motif is visible.
        return state
    if at_junction or native_support:
        evidence['reason'] = 'ambiguous_downstream_junction'
        return 'indeterminate'
    evidence['reason'] = 'downstream_motif_overlaps_numbering'
    return 'boundary_review'


def _domain_context(summary, rows, records):
    seq = summary.get('full_variable_sequence') or ''
    raw = records.get(summary.get('record_id'), summary.get('full_source_sequence', ''))
    if not rows or not seq or not raw or seq != ''.join(r.get('aa', '') for r in rows):
        return None
    try:
        indices = [int(r['input_seq_index']) for r in rows]
        positions = [(int(r['imgt_position']), (r.get('imgt_insertion') or '').strip()) for r in rows]
    except (KeyError, ValueError, TypeError):
        return None
    start, end = indices[0] - 1, indices[-1]
    if start < 0 or indices != list(range(start + 1, end + 1)) or raw[start:end] != seq:
        return None
    if len(set(positions)) != len(positions):
        return None
    return seq, raw, start, end, positions


def _select_v(column, refs):
    """Rank V candidates without the first ten positions, CDRs, or J residues."""
    scores = []
    for ref in refs:
        aligned = ref['aligned']
        if len(aligned) < 104:
            continue
        pairs = [(column[p], aligned[p - 1]) for p in _FRAMEWORK
                 if p in column and aligned[p - 1] not in '.-X']
        if len(pairs) >= _MIN_FRAMEWORK_OVERLAP:
            matches = sum(q == r for q, r in pairs)
            scores.append((matches / len(pairs), len(pairs), ref))
    if not scores:
        return [], {'internal_identity': None, 'internal_overlap': 0}
    best = max(scores, key=lambda s: (s[0], s[1]))
    close = [ref for identity, _, ref in scores
             if identity >= max(_MIN_FRAMEWORK_IDENTITY, best[0] - _NEAR_IDENTITY)]
    return close, {'internal_identity': round(best[0], 3), 'internal_overlap': best[1],
                   'internal_positions': 'IMGT 11-26,39-55,66-104'}


def _n_options(seq, positions, vrefs):
    # The shared internal anchor fixes the register. A missing position inside
    # this window, after the query has started, is not a missing N terminus.
    indices = {p: i for i, (p, ins) in enumerate(positions) if not ins}
    anchor = 23
    if anchor not in indices or seq[indices[anchor]] != 'C':
        return [], 'missing_v_anchor'
    first = positions[0][0]
    if first < 1 or first > 9 or positions[0][1]:
        return [], 'start_register_uncertain'
    stop = indices[anchor] + 1
    observed = seq[:stop]
    candidates = []
    for ref in vrefs:
        aligned = ref['aligned']
        if aligned[anchor - 1] != seq[stop - 1]:
            continue
        # Examine every retained column, including its gap pattern. This prevents
        # mislabelling internal V-gene indels or FR1 substitutions as deletions.
        reference = ''.join(c for c in aligned[first - 1:anchor] if c not in '.-')
        if len(reference) != len(observed):
            continue
        if len(observed) < 12 or sum(a == b for a, b in zip(reference, observed)) / len(observed) < 0.70:
            continue
        extension = ''.join(c for c in aligned[:first - 1] if c not in '.-')
        candidates.append((extension, ref))
    return candidates, 'terminal_register_unresolved'


def _j_options(seq, positions, jrefs):
    cys = next((i for i, pos in enumerate(positions) if pos == (104, '')), None)
    if cys is None or seq[cys] != 'C':
        return [], 'missing_conserved_cys104'
    start = next((i for i, pos in enumerate(positions) if pos == (118, '')), None)
    if start is None or seq[start] not in 'FW':
        return [], 'missing_j_anchor'
    tail = seq[start:]
    if positions[start:] != [(p, '') for p in range(118, 118 + len(tail))]:
        return [], 'j_register_uncertain'
    # Choose J candidates by a fixed common prefix before the tested endpoint.
    # The terminal query residue does not participate in choosing the J set.
    n = min(8, len(tail) - 1)
    if n < 6:
        return [], 'short_j_overlap'
    scored = [(sum(a != b for a, b in zip(tail[:n], r['aligned'][:n])), r)
              for r in jrefs if len(r['aligned']) >= n and r['aligned'][0] == tail[0]]
    if not scored:
        return [], 'no_compatible_j'
    best = min(distance for distance, _ in scored)
    if best > 1:
        return [], 'weak_j_match'
    candidates = []
    for distance, r in scored:
        if distance > min(1, best + 1):
            continue
        ref = r['aligned']
        # Internal substitutions near the last residue remain possible. A shorter
        # alternative remains evidence against a unique truncation interpretation.
        overlap = min(len(ref), len(tail))
        if sum(a != b for a, b in zip(tail[:overlap], ref[:overlap])) > 1:
            continue
        extension = ref[len(tail):] if len(ref) >= len(tail) else None
        candidates.append((extension, r))
    return candidates, 'j_endpoint_unresolved'


def _interpret(options, side, terminal_aa, flank, kind, reason):
    if not options:
        return 'indeterminate', {'reason': reason}
    extensions = sorted({ext for ext, _ in options if ext is not None})
    shorter = any(ext is None for ext, _ in options)
    evidence = {
        'candidate_extensions': extensions,
        'candidate_extension_lengths': sorted({len(ext) for ext in extensions}),
        'shorter_reference_present': shorter,
        'distinct_reference_proteins': len({r['aligned'] for _, r in options}),
        'candidate_species': sorted({r['species'] for _, r in options}),
        'reference_gene_examples': sorted({name for _, r in options for name in r['genes']})[:8],
    }
    if extensions == [''] and not shorter:
        # Register/length alone cannot distinguish a native substitution from
        # a residue borrowed at an unrecognized fusion or expression junction.
        endpoint_aas = sorted({r['aligned'].lstrip('.-')[0] if side == 'N' else r['aligned'][-1]
                               for _, r in options})
        evidence.update(reference_endpoint_residues=endpoint_aas, observed_endpoint_residue=terminal_aa)
        if terminal_aa not in endpoint_aas:
            if side == 'N' or flank:
                return 'indeterminate', {**evidence, 'reason': 'same_length_terminal_residue_disagrees'}
            return 'compatible', {**evidence, 'reason': 'endpoint_length_only_residue_differs'}
        return 'compatible', {**evidence, 'reason': 'matches_candidate_endpoint'}
    if not extensions or '' in extensions or shorter:
        return 'indeterminate', {**evidence, 'reason': 'alternative_endpoints_or_longer_query'}
    if any(len(ext) > 3 for ext in extensions):
        return 'indeterminate', {**evidence, 'reason': 'outside_short_deletion_scope'}
    if flank:
        adjacent = [ext for ext in extensions if (flank.endswith(ext) if side == 'N' else flank.startswith(ext))]
        if adjacent:
            return 'boundary_review', {**evidence, 'reason': 'candidate_residues_already_in_source',
                                       'matched_adjacent': adjacent}
        if kind != 'constant_start_motif':
            return 'indeterminate', {**evidence, 'reason': 'unresolved_flank_origin'}
    return 'suspected_short', {**evidence, 'reason': 'retained_candidates_extend_beyond_endpoint'}


def assess_endpoints(metadata, rows, *, catalog=None, use_germline=True, known_versions=None):
    """Return per-side evidence for the original numbered domain, without modifying it.

    If use_germline=False only flanking sequence context is considered. Callers
    may provide the result of compare_terminals as known_versions. A verified
    catalog endpoint takes precedence over a generic germline suggestion.
    """
    db = (catalog if catalog is not None else _catalog()) if use_germline else {}
    refs = defaultdict(list)
    for r in db.get('entries', []):
        refs[r['role'], r['segment']].append(r)
    versions = {v['chain']: v for v in (known_versions or {}).get('comparisons', [])}
    grouped = defaultdict(list)
    for r in rows:
        grouped[r.get('output_chain_id')].append(r)
    assessments, findings = [], []
    for summary in metadata['chains']:
        chain, role = summary['output_chain_id'], summary.get('chain_type')
        context = _domain_context(summary, grouped[chain], metadata.get('full_input_sequences') or {})
        options = {}
        gate = {}
        if context:
            seq, raw, start, end, positions = context
            if use_germline and role in {'H', 'K', 'L'}:
                column = {p: seq[i] for i, (p, ins) in enumerate(positions) if not ins}
                vrefs, gate = _select_v(column, refs[role, 'V'])
                if vrefs:
                    options = {'N': _n_options(seq, positions, vrefs), 'C': _j_options(seq, positions, refs[role, 'J'])}
        for side, region in [('N', 'FR1'), ('C', 'FR4')]:
            state = 'indeterminate'
            evidence = {'reason': 'source_mapping_unverified', 'chain_type': role,
                        'sequence_basis': summary.get('sequence_basis')}
            if context:
                flank = raw[:start] if side == 'N' else raw[end:]
                kind = _flank_type(side, flank, role)
                evidence.update(raw_domain_span=[start + 1, end], adjacent_length=len(flank),
                                adjacent_sequence=flank[-12:] if side == 'N' else flank[:12],
                                flank_type=kind, **gate)
                if not use_germline:
                    evidence['reason'] = 'context_only'
                    if kind == 'short_unassigned_flank':
                        state = 'boundary_review'
                else:
                    opts, reason = options.get(side, ([], 'insufficient_internal_framework_support'))
                    terminal_aa = seq[0] if side == 'N' else seq[-1]
                    state, result = _interpret(opts, side, terminal_aa, flank, kind, reason)
                    evidence.update(result)
                    # A single-S heavy endpoint exists in independently documented
                    # engineered constructs. Mark this version ambiguity explicitly
                    # even when a germline snapshot happens to contain only SS.
                    if (side == 'C' and role == 'H' and state == 'suspected_short' and seq.endswith('VTVS')
                            and 'S' in evidence.get('candidate_extensions', [])):
                        state, evidence['reason'] = 'indeterminate', 'known_single_and_double_serine_endpoints'
                    if kind == 'possible_initiator_M' and state == 'boundary_review':
                        state, evidence['reason'] = 'indeterminate', 'initiator_or_domain_residue'
                    if side == 'N' and start == 0 and seq.startswith('M') and opts:
                        if positions[0] == (1, '') and evidence['reason'] == 'same_length_terminal_residue_disagrees':
                            state, evidence['reason'] = 'boundary_review', 'ambiguous_numbered_start_methionine'
                            evidence['possible_explanations'] = ['intended_first_residue_is_methionine',
                                                                 'initiator_methionine_replaces_original_first_residue']
                        else:
                            state, evidence['reason'] = 'indeterminate', 'numbered_initial_methionine'
                if side == 'C':
                    hits = _downstream_crossing(seq, flank, role, end)
                    if hits:
                        opts, _ = options.get(side, ([], ''))
                        state = _review_crossing(hits, seq, positions, opts, state, evidence, use_germline)
            known = versions.get(chain, {})
            known_neighbors = []
            if context:
                for variant in known.get('longer_versions', []):
                    extra = variant.get('residues', '')
                    if (variant.get('side') == side and extra
                            and (flank.endswith(extra) if side == 'N' else flank.startswith(extra))
                            and variant.get('adjacent_source_residues') == extra):
                        known_neighbors.append({'reference': variant['reference'], 'residues': extra})
            if known_neighbors:
                evidence.update(germline_only_state=state, germline_only_reason=evidence['reason'],
                                exact_catalog_versions=known.get('exact_versions', []),
                                longer_versions_supported_by_flank=known_neighbors,
                                reason='source_supports_longer_catalog_boundary')
                state = 'boundary_review'
            elif known.get('exact_versions') and context:
                evidence.update(germline_only_state=state, germline_only_reason=evidence['reason'],
                                exact_catalog_versions=known['exact_versions'], reason='exact_construct_version_available')
                state = 'compatible'
            assessments.append(dict(chain=chain, side=side, state=state, evidence=evidence))
            # The construct comparison already reports a matched flank with its source URL.
            if use_germline and not known_neighbors and state in ('suspected_short', 'boundary_review'):
                code = 'germline_endpoint_short' if state == 'suspected_short' else 'germline_endpoint_boundary'
                message = ('Mapped candidates extend past this endpoint; verify the intended construct.'
                           if state == 'suspected_short' else 'Candidate residues occur adjacent to the numbered domain in the source.')
                if evidence['reason'] == 'downstream_motif_overlaps_numbering':
                    code = 'downstream_motif_overlaps_numbering'
                    message = 'A possible tag, tandem G4S linker or constant-region start overlaps the numbered domain end. Check the construct junction.'
                elif evidence['reason'] == 'ambiguous_numbered_start_methionine':
                    code = 'ambiguous_numbered_start_methionine'
                    message = ('The initial M at IMGT 1 differs from the retained V starts. An intended M and an initiator M '
                               'at the original first residue position produce the same sequence. Verify the construct record.')
                findings.append(dict(chain=chain, region=region, code=code, status='suspected',
                                     message=message, evidence=evidence))
    return dict(version=db.get('version'), assessments=assessments, findings=findings)

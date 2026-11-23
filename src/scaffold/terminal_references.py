"""Conservative terminal comparisons with versioned, locally bundled variable domains."""
import json
from importlib.resources import files


def compare_terminals(metadata):
    catalog = json.loads(files('scaffold').joinpath('static/terminal_references.json').read_text(encoding='utf-8'))
    comparisons, findings = [], []
    for chain in metadata['chains']:
        chain_id = chain['output_chain_id']
        query = chain['full_variable_sequence']
        source = chain['full_source_sequence']
        role = 'H' if chain['chain_type'] == 'H' else 'L'
        candidates = [r for r in catalog['entries'] if r['role'] == role]
        exact = [r['name'] for r in candidates if r['sequence'] == query]
        longer = []
        for ref in candidates:
            ref_seq = ref['sequence']
            if not (1 <= len(ref_seq) - len(query) <= 3 and ref_seq.count(query) == 1):
                continue
            start = ref_seq.index(query)
            end = start + len(query)
            for side, region, extra in [('N', 'FR1', ref_seq[:start]), ('C', 'FR4', ref_seq[end:])]:
                if not extra:
                    continue
                # Compare the adjacent *source* residues as well. Numbering may have stopped early.
                domain_at = source.find(query)
                adjacent = None
                if source.count(query) == 1:
                    adjacent = (source[max(0, domain_at-len(extra)):domain_at] if side == 'N' else
                                source[domain_at+len(query):domain_at+len(query)+len(extra)])
                evidence = dict(reference=ref['name'], source_url=ref['source'], side=side, residues=extra,
                                numbered_domain_length=len(query), reference_length=len(ref_seq),
                                adjacent_source_residues=adjacent, exact_matching_versions=exact)
                longer.append(evidence)
                if adjacent == extra:
                    message = f'{extra} occurs next to the numbered domain in the source and fits the longer reference. Review the extracted boundary.'
                    code = 'reference_residues_outside_numbered_domain'
                else:
                    message = f'The numbered domain is shorter at the {side} terminus by {extra} relative to {ref["name"]}. Confirm the intended version.'
                    code = 'known_reference_terminal_difference'
                findings.append(dict(chain=chain_id, region=region, code=code, status='reference_relative',
                                     message=message, evidence=evidence))
        # A short terminal motif is weak evidence for an anonymous heavy domain. Require
        # the raw sequence to end there. Never assign a construct or a missing residue from it.
        if not exact and not longer and role == 'H' and len(query) >= 9 and source.endswith(query):
            ending = query[-9:]
            motifs = [r['name'] for r in candidates if r['sequence'][-10:-1] == ending]
            if motifs:
                findings.append(dict(chain=chain_id, region='FR4', code='terminal_motif_review',
                                     status='suspected', message='The input ends at a motif that continues by one residue in a catalogued heavy domain. Check the recorded construct endpoint.',
                                     evidence=dict(terminal_motif=ending, comparison_examples=motifs)))
        state = 'compare_versions' if longer else ('exact_catalog_domain' if exact else 'no_domain_match')
        comparisons.append(dict(chain=chain_id, status=state, exact_versions=exact, longer_versions=longer))
    return dict(catalog_version=catalog['version'], catalog_size=len(catalog['entries']),
                comparisons=comparisons, findings=findings)

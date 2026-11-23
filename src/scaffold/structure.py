"""Read author-chain structure sequences and write resolved, reindexed variable domains."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import gemmi
from .errors import ScaffoldError
from .numbering import AA

def one_letter(name):
    letter = gemmi.find_tabulated_residue(name).one_letter_code
    if letter not in AA:
        raise ScaffoldError(f'Unsupported or nonstandard residue {name}; normalize explicitly before conversion')
    return letter

@dataclass
class SourceChain:
    auth_id: str
    label_id: str
    sequence: str
    residues: dict[int, object] # full sequence zero-based index -> gemmi.Residue
    sequence_basis: str

@dataclass
class SourceStructure:
    path: Path
    structure: object
    chains: dict[str,SourceChain]
    model_index: int

def read_structure(path, model_index=1, selected_chains=None, namespace='auth'):
    path = Path(path).resolve()
    st = gemmi.read_structure(str(path))
    if not 1 <= model_index <= len(st):
        raise ScaffoldError(f'Model index {model_index} is outside 1..{len(st)}')
    st.setup_entities()
    st.assign_label_seq_id()
    chains = {}
    for chain in st[model_index-1]:
        polymer = chain.get_polymer()
        if not polymer:
            continue
        if not gemmi.find_tabulated_residue(polymer[0].name).is_amino_acid():
            continue
        label = polymer[0].subchain
        key = chain.name if namespace=='auth' else label
        if selected_chains and key not in selected_chains:
            continue
        if key in chains:
            raise ScaffoldError(f'Duplicate chain ID {key}; select a unique model/chain instance')
        entity = st.get_entity_of(polymer)
        has_full = bool(entity and entity.full_sequence)
        names = list(entity.full_sequence) if has_full else [r.name for r in polymer]
        sequence = ''.join(one_letter(n) for n in names)
        residues = {}
        for i,res in enumerate(polymer):
            if has_full and res.label_seq is None:
                raise ScaffoldError(f'{key}: coordinate residue {res.seqid} could not be mapped to the full sequence')
            position = res.label_seq-1 if has_full and res.label_seq is not None else i
            if position in residues:
                raise ScaffoldError(f'{key}: microheterogeneous/duplicate residue at sequence position {position+1}')
            if position < 0 or position >= len(sequence) or sequence[position] != one_letter(res.name):
                raise ScaffoldError(f'{key}: full-sequence/coordinate mismatch at position {position+1}')
            residues[position] = res
        chains[key] = SourceChain(chain.name,label,sequence,residues,'SEQRES/entity_poly_seq' if has_full else 'coordinates_only')
    if selected_chains and set(selected_chains)-set(chains):
        raise ScaffoldError(f'Unknown/non-protein chains: {sorted(set(selected_chains)-set(chains))}; namespace={namespace}')
    if not chains:
        raise ScaffoldError('No selected protein chains')
    return SourceStructure(path,st,chains,model_index)

def domain_mapping(domain, chain: SourceChain):
    """Only accept an exact complete variable-sequence match, never a homolog mutation."""
    starts = []
    cursor = 0
    while (cursor := chain.sequence.find(domain.sequence,cursor)) >= 0:
        starts.append(cursor)
        cursor += 1
    if len(starts) != 1:
        raise ScaffoldError(f'{domain.domain_id}: expected one exact variable sequence match in chain {chain.auth_id}, found {len(starts)}')
    return {i: starts[0]+i for i in range(len(domain.sequence))}

def _clone_residue(res, seq_index, chain_id):
    clone = gemmi.Residue()
    clone.name = res.name
    clone.seqid = gemmi.SeqId(seq_index,' ')
    clone.label_seq = seq_index
    clone.subchain = chain_id
    clone.entity_type = gemmi.EntityType.Polymer
    clone.het_flag = 'A'
    best = {}
    for atom in res:
        if atom.element.is_hydrogen:
            continue
        prior = best.get(atom.name)
        if prior is None or (atom.occ,atom.altloc in ('\x00','A')) > (prior.occ,prior.altloc in ('\x00','A')):
            best[atom.name] = atom
    for atom in best.values():
        copied = atom.clone()
        copied.altloc = '\x00'
        clone.add_atom(copied)
    if not {'N','CA','C'} <= set(best):
        raise ScaffoldError(f'{chain_id}:{seq_index} missing N/CA/C backbone atoms')
    return clone

def normalize_domains(source, selections, output_path, selected_cdrs=('CDR1','CDR2','CDR3'), design_roles=None):
    """selections: (Domain, source chain key, output chain ID), one per domain."""
    new = gemmi.Structure()
    new.name = 'scaffold'
    new.cell = source.structure.cell
    new.spacegroup_hm = source.structure.spacegroup_hm or 'P 1'
    model = gemmi.Model('1')
    mapping_rows, summaries, warnings = [], [], []
    for domain,key,outid in selections:
        chain = source.chains[key]
        positions = domain_mapping(domain,chain)
        present = [i for i in range(len(domain.sequence)) if positions[i] in chain.residues]
        if not present:
            raise ScaffoldError(f'{key}: no coordinates for variable domain')
        first,last = min(present),max(present)
        holes = [i for i in range(first,last+1) if positions[i] not in chain.residues]
        if holes:
            raise ScaffoldError(f'{key}: internal missing variable-domain residues at input positions {[domain.start+i+1 for i in holes]}')
        trimmed = list(range(first)) + list(range(last+1,len(domain.sequence)))
        if any(domain.rows[i]['region'].startswith('CDR') for i in trimmed):
            raise ScaffoldError(f'{key}: trimming would discard CDR residues')
        if trimmed:
            warnings.append(f'{key}: trimmed {len(trimmed)} unresolved terminal framework residues; full mapping retained')
        outchain = gemmi.Chain(outid)
        emitted = 0
        per_chain_rows = []
        previous = None
        for i,row in enumerate(domain.rows):
            original = chain.residues.get(positions[i])
            included = first <= i <= last
            if included:
                emitted += 1
                clone = _clone_residue(original,emitted,outid)
                if previous is not None:
                    distance = previous.find_atom('C','*').pos.dist(clone.find_atom('N','*').pos)
                    if distance > 2.2:
                        raise ScaffoldError(f'{key}: backbone break ({distance:.2f} A) before source residue {original.seqid}; cannot silently concatenate')
                outchain.add_residue(clone)
                previous = clone
            role = 'heavy' if domain.chain_type=='H' else 'light'
            designed = included and row['region'] in selected_cdrs and (design_roles is None or role in design_roles)
            mapped = {**row,'source_chain_auth':chain.auth_id,'source_chain_label':chain.label_id,'source_seq_index':positions[i]+1,'source_auth_seq_id':original.seqid.num if original else None,'source_insertion_code':original.seqid.icode.strip() if original else None,'source_label_seq_id':original.label_seq if original else None,'output_chain_id':outid,'boltzgen_res_index':emitted if included else None,'included':included,'designed':designed,'structure_group':(0 if designed else 2) if included else None,'has_coordinates':original is not None}
            mapping_rows.append(mapped)
            per_chain_rows.append(mapped)
        model.add_chain(outchain)
        cdrs = []
        for region in ('CDR1','CDR2','CDR3'):
            rows = [r for r in per_chain_rows if r['region']==region]
            cdrs.append({'name':region,'sequence':''.join(r['aa'] for r in rows),'length':len(rows),'imgt_labels':[r['imgt_label'] for r in rows],'numbering_labels':[r['numbering_label'] for r in rows],'cdr_numbering_labels':[r['cdr_label'] for r in rows],'indices':[r['boltzgen_res_index'] for r in rows if r['included']],'designed':any(r['designed'] for r in rows)})
        summaries.append({'role':role,'chain_type':domain.chain_type,'output_chain_id':outid,'source_chain_id':key,'record_id':domain.record_id,'domain_id':domain.domain_id,'full_source_sequence':chain.sequence,'full_variable_sequence':domain.sequence,'selected_sequence':domain.sequence[first:last+1],'length':emitted,'sequence_basis':chain.sequence_basis,'trimmed_terminal_residues':len(trimmed),'cdrs':cdrs})
    if len(model)==2:
        # Reject obviously unrelated chain copies; this does not infer biological assembly.
        ca = [[r.find_atom('CA','*').pos for r in c] for c in model]
        minimum = min(a.dist(b) for a in ca[0] for b in ca[1])
        if minimum > 10.0:
            raise ScaffoldError(f'VH/VL have no plausible interface (closest CA distance {minimum:.2f} A). Select the correct chain pair/assembly.')
    new.add_model(model)
    new.setup_entities()
    for chain in new[0]:
        polymer = chain.get_polymer()
        entity = new.get_entity_of(polymer)
        entity.full_sequence = [r.name for r in polymer]
    new.make_mmcif_document().write_file(str(output_path))
    # Readback ensures serialized label chain IDs and sequence order are preserved.
    check = read_structure(output_path,namespace='label')
    for item in summaries:
        if check.chains[item['output_chain_id']].sequence != item['selected_sequence']:
            raise ScaffoldError('Normalized mmCIF sequence readback mismatch')
    return mapping_rows,summaries,warnings

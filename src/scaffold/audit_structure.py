"""Conservative structure input for audits, without trimming or filling residues.

Unlike conversion, an audit must retain incomplete and conflicting observations.
Indices in this module are zero-based *declared sequence* positions when a
declaration exists. Author residue numbers are identifiers, never gap evidence.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
import re

import gemmi

from .errors import ScaffoldError

BACKBONE = ("N", "CA", "C", "O")
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
DECLARED_INDEX_MARKERS = {
    'REMARK 900 RESNUM = 1-BASED SEQUENCE POSITION, NOT CHOTHIA NUMBER.',
    'REMARK 900 RESNUM IS CORRECTED SEQUENCE POSITION; H=VH/VHH L=VL.',
    'REMARK 900 RESIDUE NUMBERS ARE 1-BASED POSITIONS IN THE COMPLETE VARIABLE SEQUENCE',
    'REMARK 900 RESIDUE NUMBERS ARE 1-BASED POSITIONS IN THE COMPLETE VARIABLE SEQUENCE.',
}


@dataclass
class AuditResidue:
    aa: str
    resname: str
    auth_seq_id: int
    insertion_code: str
    label_seq_id: int | None
    auth_chain_id: str
    label_chain_id: str
    atom_names: list[str]
    positive_atom_names: list[str]
    missing_backbone: list[str]
    zero_occupancy_atoms: list[str]
    has_positive_ca: bool
    has_any_positive_atoms: bool
    nonstandard: bool
    atom_occupancy: dict[str, float]
    seq_index: int | None = None


@dataclass
class AuditChain:
    auth_id: str
    label_id: str
    sequence: str
    declared_sequence: str | None
    observed_sequence: str
    sequence_basis: str
    residues: dict[int, AuditResidue]
    mapping_status: str
    warnings: list[dict] = field(default_factory=list)
    unmapped_residues: list[AuditResidue] = field(default_factory=list)
    uncertain_sequence_indices: list[int] = field(default_factory=list)
    known_missing_positions: list[int] = field(default_factory=list)
    entity_id: str | None = None
    sequence_index_convention: str | None = None

    @property
    def mapping_complete(self):
        return self.mapping_status not in {"ambiguous", "conflict"}


@dataclass
class AuditStructure:
    path: Path
    chains: dict[str, AuditChain]
    model_index: int
    model_count: int
    warnings: list[dict]
    model_provenance: list[str]


def _aa(name):
    tab = gemmi.find_tabulated_residue(name)
    aa = tab.one_letter_code.upper()
    return aa if aa in STANDARD_AA else "X"


def _residue(residue, auth, label):
    occupancy = {}
    for atom in residue:
        occupancy[atom.name] = max(float(atom.occ), occupancy.get(atom.name, float('-inf')))
    names = sorted(occupancy)
    positive = sorted(n for n, occ in occupancy.items() if occ > 0)
    return AuditResidue(
        aa=_aa(residue.name), resname=residue.name,
        auth_seq_id=residue.seqid.num, insertion_code=residue.seqid.icode.strip().strip('\x00'),
        label_seq_id=residue.label_seq, auth_chain_id=auth, label_chain_id=label,
        atom_names=names, positive_atom_names=positive,
        missing_backbone=[n for n in BACKBONE if n not in positive],
        zero_occupancy_atoms=sorted(n for n, occ in occupancy.items() if occ <= 0),
        has_positive_ca='CA' in positive, has_any_positive_atoms=bool(positive),
        nonstandard=residue.name not in {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'},
        atom_occupancy=occupancy,
    )


def exact_subsequence_mapping(declared, observed):
    """Return all certain correspondences and ambiguity, with no substitutions.

    Every observed residue must occur in order in the declaration. The earliest
    and latest embeddings bound every valid embedding. Only positions shared by
    both bounds are asserted. Ambiguous possible coordinates are not called
    absent. A failed embedding does not authorize forced alignment or mutation.
    """
    if not observed:
        return {}, [], list(range(len(declared))), 'exact_unique'
    earliest, cursor = [], 0
    for aa in observed:
        index = declared.find(aa, cursor)
        if index < 0:
            return {}, list(range(len(declared))), [], 'conflict'
        earliest.append(index)
        cursor = index + 1
    latest, cursor = [], len(declared)
    for aa in reversed(observed):
        index = declared.rfind(aa, 0, cursor)
        if index < 0:
            return {}, list(range(len(declared))), [], 'conflict'
        latest.append(index)
        cursor = index
    latest.reverse()
    certain, possible = {}, set()
    for i, (left, right) in enumerate(zip(earliest, latest)):
        if left == right:
            certain[i] = left
            possible.add(left)
        else:
            possible.update(p for p in range(left, right + 1) if declared[p] == observed[i])
    uncertain = sorted(possible - set(certain.values()))
    missing = sorted(set(range(len(declared))) - possible)
    return certain, uncertain, missing, 'ambiguous' if uncertain else 'exact_unique'


def _map_chain(auth, label, declared, records, basis, use_labels=False, entity=None, position_convention=None):
    observed = ''.join(r.aa for r in records)
    warnings = []
    if declared is None:
        mapping = {i: i for i in range(len(records))}
        sequence, uncertain, missing, status = observed, [], [], 'coordinates_only'
        warnings.append({'code':'DECLARED_SEQUENCE_UNAVAILABLE', 'message':'Only coordinate sequence is available; unobserved terminal or internal residues cannot be enumerated from author-number gaps.'})
    elif records and (position_convention or (use_labels and all(r.label_seq_id is not None for r in records))):
        sequence, mapping = declared, {}
        uncertainty = set()
        groups = defaultdict(list)
        for i, r in enumerate(records):
            p = (r.auth_seq_id if position_convention else r.label_seq_id) - 1
            groups[p].append(i)
            if position_convention and r.insertion_code:
                uncertainty.update(range(len(declared)))
                warnings.append({'code':'DECLARED_INDEX_CONVENTION_CONFLICT','position':p+1,'message':'Explicit continuous sequence-index convention conflicts with an author insertion code; no complete mapping is asserted.'})
        for p, indices in groups.items():
            if len(indices) != 1 or p < 0 or p >= len(declared):
                uncertainty.add(p)
                warnings.append({'code':'DECLARED_INDEX_CONVENTION_CONFLICT' if position_convention else 'LABEL_SEQUENCE_CONFLICT','position':p+1,'message':'Duplicate or out-of-range declared sequence index; not treated as a known missing residue.'})
                continue
            i = indices[0]
            if declared[p] != records[i].aa:
                uncertainty.add(p)
                warnings.append({'code':'DECLARATION_COORDINATE_MISMATCH','position':p+1,'declared_aa':declared[p],'observed_aa':records[i].aa})
            else:
                mapping[i] = p
        status = 'conflict' if warnings else ('exact_declared_sequence_position' if position_convention else 'exact_label_seq_id')
        # Out-of-range coordinates invalidate completeness claims, not just one site.
        if any(p < 0 or p >= len(declared) for p in uncertainty):
            uncertainty.update(range(len(declared)))
        uncertain = sorted(p for p in uncertainty if 0 <= p < len(declared))
        missing = sorted(set(range(len(declared))) - set(mapping.values()) - set(uncertain))
    else:
        sequence = declared
        mapping, uncertain, missing, status = exact_subsequence_mapping(declared, observed)
        if use_labels and records:
            warnings.append({'code':'LABEL_SEQ_ID_UNAVAILABLE','message':'One or more coordinate label_seq_id values are absent; used conservative sequence embedding.'})
        if status == 'ambiguous':
            warnings.append({'code':'AMBIGUOUS_SEQUENCE_MAPPING','message':'Multiple exact SEQRES/coordinate embeddings exist. Only common positions are mapped; potential positions are marked uncertain.'})
        elif status == 'conflict':
            warnings.append({'code':'DECLARATION_COORDINATE_CONFLICT','message':'The coordinate sequence is not an exact subsequence of the declaration. No mismatch was silently corrected and no missing position was inferred.'})
    mapped = {}
    for i, p in mapping.items():
        records[i].seq_index = p
        mapped[p] = records[i]
    unmapped = [r for i, r in enumerate(records) if i not in mapping]
    if any(r.nonstandard for r in records) or 'X' in sequence:
        warnings.append({'code':'NONSTANDARD_OR_UNKNOWN_RESIDUE','message':'Modified/unknown residue identities are retained; canonical one-letter forms or X are used only for sequence recognition.'})
    seen = CounterKey((r.auth_seq_id, r.insertion_code) for r in records)
    if any(n > 1 for n in seen.values()):
        warnings.append({'code':'DUPLICATE_AUTHOR_RESIDUE_ID','message':'Multiple residue records share an author ID; inspect microheterogeneity before use.'})
    return AuditChain(auth, label, sequence, declared, observed, basis, mapped, status,
                      warnings, unmapped, uncertain, missing, entity, position_convention)


def CounterKey(items):
    out = defaultdict(int)
    for item in items:
        out[item] += 1
    return out


def _category(block, prefix):
    table = block.find_mmcif_category(prefix)
    columns = [x.split('.',1)[1] for x in table.tags]
    return [dict(zip(columns, (gemmi.cif.as_string(str(value)) if str(value) not in {'.','?'} else str(value) for value in row))) for row in table]


def _is_protein_name(name):
    return name == 'UNK' or gemmi.find_tabulated_residue(name).is_amino_acid()


def _pdb_chains(st, model, text):
    declared_names = defaultdict(list)
    declared_counts = {}
    position_convention = next((line.strip() for line in text.splitlines() if line.strip() in DECLARED_INDEX_MARKERS), None)
    for line in text.splitlines():
        if line.startswith('SEQRES'):
            declared_names[line[11]].extend(line[19:70].split())
            try:
                declared_counts[line[11]] = int(line[13:17])
            except ValueError:
                pass
    chains = []
    observed = defaultdict(list)
    for c in model:
        polymer = list(c.get_polymer())
        if not polymer:
            polymer = [r for r in c if r.het_flag == 'A' and _is_protein_name(r.name)]
        if polymer and any(_is_protein_name(r.name) for r in polymer):
            observed[c.name].extend(_residue(r, c.name, '') for r in polymer if _is_protein_name(r.name))
    for auth in dict.fromkeys([*declared_names, *observed]):
        names = declared_names.get(auth, [])
        if names and not any(_is_protein_name(n) for n in names):
            continue
        declared = ''.join(_aa(n) for n in names) if names else None
        chain = _map_chain(auth, '', declared, observed.get(auth, []), 'SEQRES' if names else 'coordinates_only', position_convention=position_convention)
        if auth in declared_counts and declared_counts[auth] != len(names):
            chain.warnings.append({'code':'SEQRES_LENGTH_INCONSISTENT','declared_count':declared_counts[auth],'parsed_count':len(names),'message':'SEQRES record count does not match its declared length; declaration may itself be incomplete.'})
        chains.append(chain)
    return chains


def _cif_chains(block, st, model):
    entity_rows = {r['entity_id']:r for r in _category(block, '_entity_poly.')}
    entity_types = {key:row.get('type','') for key,row in entity_rows.items()}
    names = defaultdict(dict)
    for row in _category(block, '_entity_poly_seq.'):
        names[row['entity_id']][int(row['num'])] = row['mon_id']
    asym = {r['id']:r['entity_id'] for r in _category(block,'_struct_asym.')}
    auth_from_scheme = defaultdict(set)
    for r in _category(block, '_pdbx_poly_seq_scheme.'):
        if r.get('pdb_strand_id') not in (None,'.','?'):
            auth_from_scheme[r['asym_id']].add(r['pdb_strand_id'])
    observations = defaultdict(list)
    for chain in model:
        for residue in chain:
            label = residue.subchain
            entity = asym.get(label)
            protein = entity_types.get(entity, '').startswith('polypeptide')
            if protein or (not entity_types.get(entity) and _is_protein_name(residue.name)):
                if not _is_protein_name(residue.name):
                    continue
                observations[(chain.name,label)].append(_residue(residue,chain.name,label))
    pairs = list(observations)
    seen_labels = {p[1] for p in pairs}
    for label, entity in asym.items():
        if label not in seen_labels and entity_types.get(entity,'').startswith('polypeptide'):
            ids=auth_from_scheme.get(label,set())
            pairs.append((next(iter(ids)) if len(ids)==1 else label,label))
    chains=[]
    for auth,label in pairs:
        entity=asym.get(label)
        nmap=names.get(entity,{})
        declared=None
        basis='entity_poly_seq'
        if nmap:
            declared=''.join(_aa(nmap.get(i,'UNK')) for i in range(1,max(nmap)+1))
        else:
            canonical=re.sub(r'\s+','',entity_rows.get(entity,{}).get('pdbx_seq_one_letter_code_can',''))
            if canonical and all(aa in STANDARD_AA | {'X'} for aa in canonical):
                declared=canonical
                basis='entity_poly.pdbx_seq_one_letter_code_can'
        chain=_map_chain(auth,label,declared,observations.get((auth,label),[]),
                         basis if declared is not None else 'coordinates_only',True,entity)
        if nmap and set(nmap)!=set(range(1,max(nmap)+1)):
            chain.warnings.append({'code':'ENTITY_SEQUENCE_NUMBER_GAPS','message':'Missing entity_poly_seq entries represented by X, not fabricated amino acids.'})
        chains.append(chain)
    return chains


def read_audit_structure(path, model_index=1, selected_chains=None, namespace='auth'):
    """Read a single model with explicit evidence and conservative mapping.

    ``known_missing_positions`` are reliable absent-coordinate sequence indices;
    ``uncertain_sequence_indices`` must not be treated as confirmed missing.
    Residues with zero occupancy still exist as records but fail positive-atom
    coverage. No selection guesses between repeated antibody copies are made.
    """
    path=Path(path).resolve()
    if namespace not in {'auth','label'}:
        raise ScaffoldError('Chain namespace must be auth or label')
    try:
        text=path.read_text(encoding='utf-8-sig')
        st=gemmi.read_structure(str(path))
    except (OSError,RuntimeError,ValueError) as exc:
        raise ScaffoldError(f'Cannot read structure {path}: {exc}') from exc
    if not 1 <= model_index <= len(st):
        raise ScaffoldError(f'Model index {model_index} is outside 1..{len(st)}')
    st.setup_entities()
    model=st[model_index-1]
    is_cif=path.suffix.lower() in {'.cif','.mmcif'} or text.lstrip().startswith('data_')
    provenance=[]
    if is_cif:
        block=gemmi.cif.read_file(str(path)).sole_block()
        candidates=_cif_chains(block,st,model)
        for tag in ('_exptl.method','_ma_model_list.model_type','_ma_model_list.model_name'):
            for value in block.find_values(tag):
                provenance.append(f'{tag}: {gemmi.cif.as_string(value)}')
    else:
        if namespace=='label':
            raise ScaffoldError('PDB has no native label_asym_id namespace; select author chain IDs')
        candidates=_pdb_chains(st,model,text)
        provenance=[line.rstrip() for line in text.splitlines() if line.startswith(('REMARK','EXPDTA')) and re.search(r'COMPUT|MODEL|PREDICT|PROVENANCE|PLACEHOLDER|EXPERIMENT|RESNUM|SEQUENCE POSITION',line,re.I)]
    chains={}
    selected=set(selected_chains or [])
    for chain in candidates:
        key=chain.auth_id if namespace=='auth' else chain.label_id
        if selected and key not in selected:
            continue
        if key in chains:
            raise ScaffoldError(f'Duplicate chain ID {key} in {namespace} namespace; use label chain IDs for distinct polymer instances')
        chains[key]=chain
    if selected-set(chains):
        raise ScaffoldError(f'Unknown/non-protein chains: {sorted(selected-set(chains))}; namespace={namespace}')
    if not chains:
        raise ScaffoldError('No selected protein chains')
    warnings=[]
    if len(st)>1:
        warnings.append({'code':'MULTIPLE_MODELS','message':f'Only model {model_index} of {len(st)} was audited.'})
    return AuditStructure(path,chains,model_index,len(st),warnings,provenance)

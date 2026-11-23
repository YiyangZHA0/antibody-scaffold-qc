"""Audit input evidence must not turn unknown mapping into confirmed loss."""
from pathlib import Path

import gemmi
import pytest

from scaffold.audit_structure import exact_subsequence_mapping, read_audit_structure
from scaffold.errors import ScaffoldError


def make_structure(records, declared=None):
    """records: (resname, author number, insertion, label position, occupancy)."""
    st=gemmi.Structure()
    model=gemmi.Model('1')
    chain=gemmi.Chain('H')
    for resname, number, insertion, label, occupancy in records:
        residue=gemmi.Residue()
        residue.name=resname
        residue.seqid=gemmi.SeqId(number,insertion)
        residue.label_seq=label
        residue.subchain='A'
        residue.entity_type=gemmi.EntityType.Polymer
        residue.het_flag='A'
        for j,name in enumerate(('N','CA','C','O')):
            atom=gemmi.Atom()
            atom.name=name
            atom.element=gemmi.Element(name[0])
            atom.pos=gemmi.Position(number*4+j,0,0)
            atom.occ=occupancy
            residue.add_atom(atom)
        chain.add_residue(residue)
    model.add_chain(chain)
    st.add_model(model)
    entity=gemmi.Entity('1')
    entity.entity_type=gemmi.EntityType.Polymer
    entity.polymer_type=gemmi.PolymerType.PeptideL
    entity.subchains=['A']
    if declared:
        entity.full_sequence=declared
    st.entities.append(entity)
    return st


def pdb_file(tmp_path, records, declared=None):
    st=make_structure(records)
    header=''
    if declared:
        header=''.join(f"SEQRES {i//13+1:3d} H {len(declared):4d}  {' '.join(declared[i:i+13])}\n" for i in range(0,len(declared),13))
    path=tmp_path/'input.pdb'
    path.write_text(header+st.make_pdb_string(),encoding='ascii')
    return path


def cif_file(tmp_path, records, declared):
    st=make_structure(records,declared)
    doc=st.make_mmcif_document()
    doc.sole_block().set_pair('_exptl.method', "'X-RAY DIFFRACTION'")
    path=tmp_path/'input.cif'
    doc.write_file(str(path))
    return path


def test_exact_embedding_retains_all_ambiguity():
    mapping, uncertain, missing, status=exact_subsequence_mapping('AAACD','AACD')
    assert status=='ambiguous'
    assert mapping=={2:3,3:4}
    assert uncertain==[0,1,2]
    assert missing==[]  # one A is absent, but its exact identity is unknowable.


def test_missing_positions_require_declaration_evidence(tmp_path):
    records=[('CYS',10,' ',None,1),('ASP',20,' ',None,1),('PHE',30,' ',None,1)]
    path=pdb_file(tmp_path,records,['ALA','CYS','ASP','GLU','PHE'])
    chain=read_audit_structure(path).chains['H']
    assert chain.sequence=='ACDEF'
    assert chain.observed_sequence=='CDF'
    assert chain.mapping_status=='exact_unique'
    assert chain.known_missing_positions==[0,3]
    assert list(chain.residues)==[1,2,4]
    path=pdb_file(tmp_path,records)
    coordinate_only=read_audit_structure(path).chains['H']
    assert coordinate_only.sequence=='CDF'
    assert coordinate_only.sequence_basis=='coordinates_only'
    assert coordinate_only.known_missing_positions==[]
    assert len(coordinate_only.residues)==3  # author 10/20/30 does not imply 27 losses.


def test_mismatch_is_reported_without_forcing_an_alignment(tmp_path):
    records=[('ALA',1,' ',None,1),('GLU',2,' ',None,1),('ASP',3,' ',None,1)]
    path=pdb_file(tmp_path,records,['ALA','CYS','ASP'])
    chain=read_audit_structure(path).chains['H']
    assert chain.mapping_status=='conflict'
    assert chain.residues=={}
    assert chain.known_missing_positions==[]
    assert chain.uncertain_sequence_indices==[0,1,2]
    assert len(chain.unmapped_residues)==3


def test_native_cif_label_coordinates_and_author_insertion_retained(tmp_path):
    records=[('ALA',101,'A',1,1),('ASP',950,' ',3,1)]
    path=cif_file(tmp_path,records,['ALA','CYS','ASP','GLU'])
    source=read_audit_structure(path,selected_chains=['A'],namespace='label')
    chain=source.chains['A']
    assert chain.auth_id=='H' and chain.label_id=='A'
    assert chain.sequence_basis=='entity_poly_seq'
    assert chain.mapping_status=='exact_label_seq_id'
    assert chain.known_missing_positions==[1,3]
    assert chain.residues[0].auth_seq_id==101
    assert chain.residues[0].insertion_code=='A'
    assert chain.residues[0].label_seq_id==1
    assert source.model_provenance==['_exptl.method: X-RAY DIFFRACTION']


def test_cif_identity_conflict_is_not_marked_missing(tmp_path):
    records=[('ALA',1,' ',1,1),('GLU',2,' ',2,1)]
    path=cif_file(tmp_path,records,['ALA','CYS','ASP'])
    chain=read_audit_structure(path).chains['H']
    assert chain.mapping_status=='conflict'
    assert chain.known_missing_positions==[2]
    assert chain.uncertain_sequence_indices==[1]
    assert chain.unmapped_residues[0].aa=='E'


def test_zero_occupancy_backbone_and_partial_atom_coverage(tmp_path):
    records=[('ALA',1,' ',1,0),('CYS',2,' ',2,1)]
    st=make_structure(records,['ALA','CYS'])
    st[0]['H'][1].remove_atom('O','\x00',gemmi.Element('O'))
    path=tmp_path/'zero.cif'
    st.make_mmcif_document().write_file(str(path))
    chain=read_audit_structure(path).chains['H']
    assert chain.known_missing_positions==[]
    assert chain.residues[0].has_positive_ca is False
    assert chain.residues[0].has_any_positive_atoms is False
    assert chain.residues[0].missing_backbone==['N','CA','C','O']
    assert chain.residues[0].zero_occupancy_atoms==['C','CA','N','O']
    assert chain.residues[1].has_positive_ca
    assert chain.residues[1].missing_backbone==['O']


def test_native_canonical_entity_sequence_fallback(tmp_path):
    path=cif_file(tmp_path,[('ALA',1,' ',1,1),('ASP',3,' ',3,1)],['ALA','CYS','ASP'])
    doc=gemmi.cif.read_file(str(path))
    block=doc.sole_block()
    block.find_mmcif_category('_entity_poly_seq.').erase()
    block.set_mmcif_category('_entity_poly.', {'entity_id':['1'], 'type':['polypeptide(L)'], 'pdbx_seq_one_letter_code_can':['ACD']})
    doc.write_file(str(path))
    chain=read_audit_structure(path).chains['H']
    assert chain.sequence=='ACD'
    assert chain.sequence_basis=='entity_poly.pdbx_seq_one_letter_code_can'
    assert chain.known_missing_positions==[1]


def test_unknown_chain_and_invalid_model_are_explicit_errors(tmp_path):
    path=pdb_file(tmp_path,[('ALA',1,' ',None,1)])
    with pytest.raises(ScaffoldError,match='Unknown'):
        read_audit_structure(path,selected_chains=['Z'])
    with pytest.raises(ScaffoldError,match='Model index'):
        read_audit_structure(path,model_index=2)
    with pytest.raises(ScaffoldError,match='no native label'):
        read_audit_structure(path,namespace='label')


def test_computational_provenance_is_not_inferred_from_bfactor(tmp_path):
    path=pdb_file(tmp_path,[('ALA',1,' ',None,1)])
    source=read_audit_structure(path)
    assert not any('COMPUTATIONALLY' in r for r in source.model_provenance)
    path.write_text('REMARK 900 COMPUTATIONALLY COMPLETED VARIABLE-DOMAIN SCAFFOLD.\n'+path.read_text(),encoding='ascii')
    assert any('COMPUTATIONALLY' in r for r in read_audit_structure(path).model_provenance)


def test_real_cif_framework_gap_is_available_without_conversion(tmp_path):
    path=Path(__file__).parent/'data/6wgb.cif'
    source=read_audit_structure(path,selected_chains=['A','B'])
    assert 61 in source.chains['A'].known_missing_positions
    assert set(range(100,110)) <= set(source.chains['A'].known_missing_positions)
    assert {31,32} <= set(source.chains['B'].known_missing_positions)


def test_uncertain_coordinate_positions_are_not_confirmed_missing(tmp_path):
    from scaffold.audit import AuditDomain, _assess
    path=pdb_file(tmp_path,[('ALA',1,' ',None,1),('ALA',2,' ',None,1),('CYS',3,' ',None,1),('ASP',4,' ',None,1)],['ALA','ALA','ALA','CYS','ASP'])
    chain=read_audit_structure(path).chains['H']
    rows=[{'input_seq_index':i+1,'domain_seq_index':i+1,'aa':aa,'imgt_position':i+1,'imgt_label':str(i+1),'region':'FR1'} for i,aa in enumerate(chain.sequence)]
    domain=AuditDomain('H','H:domain1','H',0,len(chain.sequence),chain.sequence,rows,1)
    assessed=_assess(domain,chain.sequence,coordinate_chain=chain)
    assert [r['coordinate_state'] for r in assessed['rows'][:3]]==['mapping_uncertain']*3
    assert not any(r['code']=='missing_residue_coordinates' for r in assessed['risks'])
    assert assessed['regions']['FR1']['coordinate_mapping_uncertain_count']==3


def test_paired_multidomain_record_reports_its_unmatched_domain(tmp_path):
    from scaffold.audit import audit_file
    heavy='A'*30+'C'*5
    light='D'*30+'E'*5

    class FakeNumberer:
        def _number_one(self,seq):
            for domain,kind in ((heavy,'H'),(light,'K')):
                if domain in seq:
                    return [([((i+1,''),aa) for i,aa in enumerate(domain)],kind,1,seq.index(domain))]
            return []

    records=[('ALA' if aa=='A' else 'CYS',i+1,' ',None,1) for i,aa in enumerate(heavy)]
    path=pdb_file(tmp_path,records,['ALA']*30+['CYS']*5)
    paired=tmp_path/'fusion.fasta'
    paired.write_text('>fusion\n'+heavy+'GGGGS'*3+light+'\n',encoding='ascii')
    report=audit_file(path,FakeNumberer(),sequence_path=paired)
    assert report['domain_count']==1
    assert report['risk_counts'].get('unmatched_fasta_domain')==1


def test_explicit_normalized_sequence_indices_disambiguate_repeated_residues(tmp_path):
    path=pdb_file(tmp_path,[('GLN',3,' ',None,1),('LEU',4,' ',None,1)],['GLN','VAL','GLN','LEU'])
    ordinary=read_audit_structure(path).chains['H']
    assert ordinary.mapping_status=='ambiguous'
    assert ordinary.known_missing_positions==[1]
    marker='REMARK 900 RESNUM = 1-BASED SEQUENCE POSITION, NOT CHOTHIA NUMBER.'
    path.write_text(marker+'\n'+path.read_text(),encoding='ascii')
    normalized=read_audit_structure(path).chains['H']
    assert normalized.mapping_status=='exact_declared_sequence_position'
    assert normalized.known_missing_positions==[0,1]
    assert normalized.sequence_index_convention==marker
    assert normalized.residues[2].label_seq_id is None  # Do not invent a native CIF label for PDB.


def test_corrupt_explicit_sequence_index_convention_cannot_force_mapping(tmp_path):
    path=pdb_file(tmp_path,[('GLN',42,' ',None,1),('LEU',43,' ',None,1)],['GLN','VAL','GLN','LEU'])
    marker='REMARK 900 RESNUM IS CORRECTED SEQUENCE POSITION; H=VH/VHH L=VL.'
    path.write_text(marker+'\n'+path.read_text(),encoding='ascii')
    chain=read_audit_structure(path).chains['H']
    assert chain.mapping_status=='conflict'
    assert not chain.residues
    assert not chain.known_missing_positions
    assert chain.uncertain_sequence_indices==[0,1,2,3]
    assert any(w['code']=='DECLARED_INDEX_CONVENTION_CONFLICT' for w in chain.warnings)

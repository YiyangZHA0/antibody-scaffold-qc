from pathlib import Path
import json
import gemmi
import pytest
import yaml
from scaffold.cli import execute, parser, main
from scaffold.errors import ScaffoldError
from scaffold.numbering import Numberer, classify, read_fasta, region_for
from scaffold.structure import read_structure, normalize_domains, domain_mapping
from scaffold.resolve import select_structure_domains
from scaffold.validation import validate_scaffold

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'tests/data'

@pytest.fixture(scope='session')
def numberer():
    return Numberer()

@pytest.fixture(scope='session')
def fab(numberer):
    return classify(read_fasta(DATA/'5y9k.fasta'),numberer,'fab')

def source_and_selection(fab):
    source=read_structure(DATA/'5y9k.pdb')
    return source,select_structure_domains(source,fab['domains'])

def test_fab_labels_and_full_chain_offsets(fab):
    records=read_fasta(DATA/'5y9k.fasta')
    assert [d.chain_type for d in fab['domains']]==['H','L']
    assert [len(d.sequence) for d in fab['domains']]==[123,108]
    for domain in fab['domains']:
        assert len(records[domain.record_id])>len(domain.sequence)
        for row in domain.rows:
            assert row['aa']==records[row['record_id']][row['input_seq_index']-1]

def test_real_scfv_reverse_orientation_and_non_g4s_linker(numberer):
    records=read_fasta(DATA/'1x9q.fasta')
    assert len(records)==1 and len(next(iter(records.values())))<300
    result=classify(records,numberer)
    assert result['molecule_type']=='scfv'
    assert result['linkers'][0]['orientation']=='K-H'
    assert result['linkers'][0]['matched_pattern'] is None
    assert result['domains'][0].start>result['domains'][1].start
    for d in result['domains']:
        assert ''.join(r['aa'] for r in d.rows)==records[d.record_id][d.start:d.end]

@pytest.mark.parametrize('reverse',[False,True])
def test_derived_g4s_scfv_both_orientations(fab,numberer,reverse):
    # Synthetic construct from real 5Y9K VH/VL, explicitly only a regression fixture.
    domains=fab['domains'][::-1] if reverse else fab['domains']
    seq=domains[0].sequence+'GGGGS'*3+domains[1].sequence
    result=classify({'derived':seq},numberer)
    assert result['molecule_type']=='scfv'
    assert result['linkers'][0]['sequence']=='GGGGS'*3
    assert result['linkers'][0]['matched_pattern']
    assert [d.chain_type for d in result['domains']]==['H','L']

def test_vhh_antigen_ignored_and_insertions_preserved(numberer):
    result=classify(read_fasta(DATA/'7eow.fasta'),numberer,'vhh')
    assert result['ignored_records']
    assert len(result['domains'])==1
    rows=result['domains'][0].rows
    inserted=[r for r in rows if r['imgt_insertion']]
    assert inserted and all(r['region']=='CDR3' for r in inserted)
    assert len({r['imgt_label'] for r in rows})==len(rows)

@pytest.mark.parametrize('position,expected',[(26,'FR1'),(27,'CDR1'),(38,'CDR1'),(39,'FR2'),(55,'FR2'),(56,'CDR2'),(65,'CDR2'),(66,'FR3'),(104,'FR3'),(105,'CDR3'),(117,'CDR3'),(118,'FR4')])
def test_imgt_boundaries(position,expected):
    assert region_for(position)==expected

@pytest.mark.parametrize('content', ['>same\nAAAA\n>same\nCCCC\n','>unknown\nEVQLX\n',''])
def test_invalid_fasta_rejected(tmp_path,content):
    path=tmp_path/'input.fasta';path.write_text(content)
    with pytest.raises(ScaffoldError):
        read_fasta(path)

def test_multiple_heavy_domains_rejected(fab,numberer):
    seq=fab['domains'][0].sequence
    with pytest.raises(ScaffoldError,match='unambiguously'):
        classify({'H1':seq,'H2':seq},numberer)

def test_internal_missing_residue_rejected(tmp_path,fab):
    source,selections=source_and_selection(fab)
    d,key,_=selections[0]
    idx=next(i for i,r in enumerate(d.rows) if r['region']=='CDR2')
    del source.chains[key].residues[domain_mapping(d,source.chains[key])[idx]]
    with pytest.raises(ScaffoldError,match='internal missing'):
        normalize_domains(source,selections,tmp_path/'missing.cif')

def test_auth_insertion_and_terminal_missing_mapping(tmp_path,fab):
    source,selections=source_and_selection(fab)
    d,key,_=selections[0]
    chain=source.chains[key]
    pos=domain_mapping(d,chain)[30]
    chain.residues[pos].seqid=gemmi.SeqId(777,'A')
    rows,summary,warnings=normalize_domains(source,selections,tmp_path/'normal.cif')
    changed=[r for r in rows if r['source_auth_seq_id']==777]
    assert changed[0]['source_insertion_code']=='A'
    assert changed[0]['boltzgen_res_index']!=777
    assert sum(not r['included'] for r in rows)==3
    assert all(not r['region'].startswith('CDR') for r in rows if not r['included'])
    assert len(warnings)==2

def test_distant_chain_pair_rejected(tmp_path,fab):
    source,selections=source_and_selection(fab)
    for residue in source.chains[selections[1][1]].residues.values():
        for atom in residue:
            atom.pos.x+=1000
    with pytest.raises(ScaffoldError,match='no plausible interface'):
        normalize_domains(source,selections,tmp_path/'distant.cif')

def test_similar_but_nonidentical_structure_rejected(fab):
    source,_=source_and_selection(fab)
    for chain in source.chains.values():
        chain.sequence=chain.sequence.replace('C','A')
    with pytest.raises(ScaffoldError,match='No exact'):
        select_structure_domains(source,fab['domains'])

def test_preexisting_output_untouched(tmp_path):
    out=tmp_path/'keep';out.mkdir();(out/'mine.txt').write_text('important')
    assert main(['from-pdb',str(DATA/'5y9k.pdb'),'--out',str(out)])==2
    assert [f.name for f in out.iterdir()]==['mine.txt']
    assert (out/'mine.txt').read_text()=='important'

def test_fasta_without_structure_leaves_annotation_only(tmp_path,fab,numberer):
    path=tmp_path/'custom.fasta'
    path.write_text(''.join(f'>custom_{d.chain_type}\n{d.sequence}\n' for d in fab['domains']))
    out=tmp_path/'annotation_only'
    with pytest.raises(ScaffoldError,match='needs coordinates'):
        execute(parser().parse_args(['from-fasta',str(path),'--mode','structure','--out',str(out)]),numberer)
    assert (out/'numbering.csv').exists() and (out/'annotation.json').exists()
    assert not (out/'scaffold.yaml').exists() and not (out/'scaffold.cif').exists()

@pytest.mark.integration
def test_official_parser_detects_mask_tampering_and_invalidates_status(tmp_path,numberer):
    out=tmp_path/'official'
    result=execute(parser().parse_args(['from-pdb',str(DATA/'5y9k.pdb'),'--out',str(out)]),numberer)
    assert result['status']=='scaffold_validated'
    document=yaml.safe_load((out/'scaffold.yaml').read_text())
    document['design'][0]['chain']['res_index']='1'
    (out/'scaffold.yaml').write_text(yaml.safe_dump(document,sort_keys=False))
    with pytest.raises(ScaffoldError,match='design mask disagrees|structure groups disagree'):
        validate_scaffold(out,ROOT/'.cache/mols.zip',run_cli=False)
    assert not json.loads((out/'validation.json').read_text())['official_parser_passed']
    assert json.loads((out/'metadata.json').read_text())['validation_status']=='failed'

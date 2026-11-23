"""Real input regressions for exact matching, structural edits and validation."""
from pathlib import Path
from unittest.mock import patch
import json
import shutil
import numpy as np
import pytest
import yaml
from scaffold.cli import execute,parser
from scaffold.compare import compare_outputs,compare_coordinates
from scaffold.errors import ScaffoldError
from scaffold.numbering import Numberer,read_fasta,classify
from scaffold.output import parse_loop_lengths
from scaffold.resolve import resolve_fasta_structure
from scaffold.validation import validate_scaffold

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'tests/data'


@pytest.fixture(scope='module')
def numberer():
    return Numberer()


@pytest.fixture(scope='module')
def native(tmp_path_factory,numberer):
    out=tmp_path_factory.mktemp('real_native')
    execute(parser().parse_args(['convert',str(DATA/'7eow.pdb'),'--vhh-chain','B','--type','vhh','--name','caplacizumab','--out',str(out)]),numberer)
    return out


def test_requested_schema_and_imgt_header(native):
    text=(native/'scaffold.yaml').read_text()
    doc=yaml.safe_load(text)
    assert set(doc)=={'path','include','design','structure_groups','exclude','design_insertions','reset_res_index'}
    assert '# H-FR1:' in text and '# H-CDR3:' in text and '111A' in text
    assert doc['design']==[{'chain':{'id':'H','res_index':'26..33,51..58,97..117'}}]
    assert [x['insertion']['num_residues'] for x in doc['design_insertions']]==[8,8,21]
    assert (native/'scaffold.yml').read_bytes()==(native/'scaffold.yaml').read_bytes()
    parsed=[json.loads(line) for line in (native/'parsed_residue_map.jsonl').read_text().splitlines()]
    assert sum(r['designed'] for r in parsed)==37
    assert all(r['source_res_index'] is None for r in parsed if r['designed'])


@pytest.mark.parametrize('field',['include','design','structure_groups','exclude','design_insertions','reset_res_index'])
def test_official_parser_rejects_wrong_semantic_edits(native,tmp_path,field):
    out=tmp_path/'tampered';out.mkdir()
    for src in native.iterdir():
        if src.is_file():
            shutil.copy2(src,out/src.name)
    doc=yaml.safe_load((out/'scaffold.yaml').read_text())
    if field=='include':
        doc[field][0]['chain']['res_index']='2..128'
    elif field=='design':
        doc[field][0]['chain']['res_index']='1'
    elif field=='structure_groups':
        doc[field][0]['group']['visibility']=1
    elif field=='exclude':
        doc[field][0]['chain']['res_index']='26..32'
    elif field=='design_insertions':
        doc[field][0]['insertion']['res_index']=1
    else:
        doc[field]=[]
    (out/'scaffold.yaml').write_text(yaml.safe_dump(doc,sort_keys=False))
    with pytest.raises(ScaffoldError):
        validate_scaffold(out,ROOT/'.cache/mols.zip',run_cli=False)
    assert json.loads((out/'metadata.json').read_text())['validation_status']=='failed'


@pytest.mark.parametrize('mutation',['widened_range','missing','extra','position','chain','duplicate'])
def test_structural_insertions_checked_before_sampling(native,tmp_path,mutation):
    out=tmp_path/'changed_insertion';out.mkdir()
    for src in native.iterdir():
        if src.is_file():
            shutil.copy2(src,out/src.name)
    path=out/'scaffold.yaml'
    doc=yaml.safe_load(path.read_text())
    insertion=doc['design_insertions'][0]['insertion']
    if mutation=='widened_range':
        # Previously an 8..9 specification passed when the parser sampled 8.
        length=insertion['num_residues']
        insertion['num_residues']=f'{length}..{length+1}'
    elif mutation=='missing':
        doc['design_insertions'].pop()
    elif mutation=='extra':
        doc['design_insertions'].append({'insertion':{'id':'H','res_index':1,'num_residues':1}})
    elif mutation=='position':
        insertion['res_index']+=1
    elif mutation=='chain':
        insertion['id']='L'
    else:
        doc['design_insertions'].append({'insertion':dict(insertion)})
    path.write_text(yaml.safe_dump(doc,sort_keys=False))
    with patch('numpy.random.choice',side_effect=lambda a,*args,**kwargs:np.min(a)):
        with pytest.raises(ScaffoldError,match='insertion'):
            validate_scaffold(out,ROOT/'.cache/mols.zip',run_cli=False)
    assert json.loads((out/'metadata.json').read_text())['validation_status']=='failed'


def test_custom_fasta_companion_offline_and_exact_match(tmp_path,numberer):
    path=tmp_path/'custom.fasta'
    records=read_fasta(DATA/'7eow.fasta')
    records={f'custom_{i}':seq for i,seq in enumerate(records.values())}
    path.write_text(''.join(f'>{rid}\n{seq}\n' for rid,seq in records.items()))
    shutil.copy2(DATA/'7eow.cif',path.with_suffix('.cif'))
    annotation=classify(records,numberer,'vhh')
    with patch('requests.sessions.Session.request',side_effect=AssertionError('Must use local companion')):
        source,selection,provenance=resolve_fasta_structure(annotation,records,tmp_path/'cache',chain_map={'heavy':'B'},input_path=path)
    assert len(selection)==1 and source.path==path.with_suffix('.cif')
    assert 'exact' in provenance['resolution_method']
    # Even one mismatch must not be accepted as a structural template for fixed FR.
    annotation['domains'][0].sequence='A'+annotation['domains'][0].sequence[1:]
    with pytest.raises(ScaffoldError,match='No exact'):
        resolve_fasta_structure(annotation,records,tmp_path/'cache',chain_map={'heavy':'B'},input_path=path)


def test_duplicate_structure_copies_require_chain_selection(numberer,tmp_path):
    records=read_fasta(DATA/'7xl0.fasta')
    ann=classify(records,numberer,'vhh')
    with pytest.raises(ScaffoldError,match='Ambiguous'):
        resolve_fasta_structure(ann,records,tmp_path/'cache',structure_path=DATA/'7xl0.cif')


@pytest.mark.parametrize('fmt',['pdb','cif','fasta'])
def test_real_missing_internal_framework_rejected(tmp_path,numberer,fmt):
    out=tmp_path/f'incomplete_{fmt}'
    with pytest.raises(ScaffoldError,match='internal missing'):
        execute(parser().parse_args(['convert',str(DATA/f'6wgb.{fmt}'),'--mode','structure','--heavy-chain','A','--light-chain','B','--out',str(out)]),numberer)
    assert not (out/'scaffold.yaml').exists()


def test_variable_loop_lengths_official_parser_min_and_max(tmp_path,numberer):
    out=tmp_path/'variable'
    execute(parser().parse_args(['convert',str(DATA/'7eow.cif'),'--vhh-chain','B','--type','vhh','--cdr-length','H-CDR1=6..10','--cdr-length','H-CDR3=8..18','--out',str(out)]),numberer)
    for select,expected in [(np.min,[6,8,8]),(np.max,[10,8,18])]:
        # Control sampling inputs to the unmodified official parser, to verify
        # both extreme lengths and downstream FR positions/CA coordinates.
        with patch('numpy.random.choice',side_effect=lambda a,*args,**kwargs:select(a)):
            result=validate_scaffold(out,ROOT/'.cache/mols.zip',run_cli=False)
        assert [x['sampled_length'] for x in result['chains'][0]['loops']]==expected
        assert result['chains'][0]['fixed_residues']==91


@pytest.mark.parametrize('spec',['H-CDR3=20..5','H-CDR4=5','L-CDR3=5','H-CDR3=0','H-CDR3=101'])
def test_invalid_loop_specs(native,spec):
    summaries=json.loads((native/'metadata.json').read_text())['chains']
    with pytest.raises(ScaffoldError):
        parse_loop_lengths([spec],summaries)


def test_compare_detects_changed_artifact(native,tmp_path):
    out=tmp_path/'changed';out.mkdir()
    for src in native.iterdir():
        if src.is_file():
            shutil.copy2(src,out/src.name)
    assert compare_outputs([native,out])['passed']
    with (out/'scaffold.yaml').open('a') as handle:
        handle.write('# edited after check\n')
    with pytest.raises(ScaffoldError,match='changed since validation'):
        compare_outputs([native,out])


def test_pdb_precision_and_real_coordinate_mismatch():
    cif=[('H',11,'LEU','CG',-41.46350,79.74467,47.41497)]
    pdb=[('H',11,'LEU','CG',-41.463,79.745,47.415)]
    assert compare_coordinates(cif,pdb)['passed']
    changed=[('H',11,'LEU','CG',-41.46,79.745,47.415)]
    assert not compare_coordinates(cif,changed)['passed']
    renamed=[('L',11,'LEU','CG',-41.463,79.745,47.415)]
    assert not compare_coordinates(cif,renamed)['passed']


@pytest.mark.parametrize('mode,lengths',[('structure',[117,110]),('sequence',[117,112])])
def test_previous_scfv_still_converts_when_explicitly_split(tmp_path,numberer,mode,lengths):
    out=tmp_path/'scfv'
    result=execute(parser().parse_args(['convert',str(ROOT/'tests/data/1x9q.fasta'),'--mode',mode,'--split-scfv','--out',str(out)]),numberer)
    assert [c['chain'] for c in result['chains']]==['H','L']
    # FASTA keeps two unresolved terminal VL residues which structure mode trims.
    assert [c['residues'] for c in result['chains']]==lengths

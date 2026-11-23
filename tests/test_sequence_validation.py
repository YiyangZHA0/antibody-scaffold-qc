"""Independent sequence fixtures exercise the unmodified official parser."""
import json
import os
from pathlib import Path
from unittest.mock import patch

import pytest
import yaml

from scaffold.errors import ScaffoldError
from scaffold.validation import validate_scaffold


ROOT=Path(__file__).resolve().parents[1]


@pytest.fixture(scope='module')
def moldir():
    candidates=[os.environ.get('SCAFFOLD_MOLDIR'), ROOT/'.cache/mols.zip']
    for candidate in candidates:
        if candidate and Path(candidate).exists():
            return Path(candidate)
    pytest.skip('Official BoltzGen molecule dictionary is unavailable')


@pytest.fixture
def sequence_output(tmp_path):
    # Explicit fixture answers; no numbering or production YAML writer is used.
    chunks=[('FR1','ACD',False),('CDR1','YY',True),('FR2','EF',False),
            ('CDR2','WWW',True),('FR3','GH',False),('CDR3','PQR',True),
            ('FR4','IKL',False)]
    rows=[]
    for region,sequence,designed in chunks:
        for aa in sequence:
            rows.append({'output_chain_id':'H','boltzgen_res_index':len(rows)+1,
                         'aa':aa,'region':region,'included':True,'designed':designed})
    # Include an entirely fixed chain to verify no-design entities too.
    rows.extend({'output_chain_id':'L','boltzgen_res_index':i+1,'aa':aa,
                 'region':'FR1','included':True,'designed':False}
                for i,aa in enumerate('MNPQ'))
    metadata={'scaffold_kind':'sequence','chains':[{'output_chain_id':'H'},
                  {'output_chain_id':'L'}],'loop_lengths':{'H-CDR1':2,
                  'H-CDR2':3,'H-CDR3':'2..5'},'numbering_scheme':'imgt',
                  'cdr_definition':'imgt'}
    document={'entities':[{'protein':{'id':'H','sequence':'ACD,2,EF,3,GH,2..5,IKL','msa':-1}},
                           {'protein':{'id':'L','sequence':'MNPQ','msa':-1}}]}
    (tmp_path/'metadata.json').write_text(json.dumps(metadata),encoding='utf-8')
    (tmp_path/'residue_map.jsonl').write_text(''.join(json.dumps(row)+'\n' for row in rows),encoding='utf-8')
    (tmp_path/'scaffold.yaml').write_text(yaml.safe_dump(document),encoding='utf-8')
    return tmp_path


@pytest.mark.parametrize('choose_max',[False,True])
def test_sequence_parser_length_endpoints_and_fixed_chain(sequence_output,moldir,choose_max):
    with patch('numpy.random.randint',side_effect=lambda lo,hi:hi-1 if choose_max else lo):
        result=validate_scaffold(sequence_output,moldir,run_cli=False)
    heavy,light=result['chains']
    assert [loop['sampled_length'] for loop in heavy['loops']]==[2,3,5 if choose_max else 2]
    assert heavy['fixed_residues']==10
    assert light['designed_residues']==0 and light['fixed_residues']==4
    assert not result['has_input_coordinates']
    assert all(not chain['has_input_coordinates'] for chain in result['chains'])
    assert 'scaffold.cif' not in result['artifact_sha256']
    assert not (sequence_output/'scaffold.cif').exists()
    assert (sequence_output/'check.yaml').read_bytes()==(sequence_output/'scaffold.yaml').read_bytes()
    assert (sequence_output/'design.yaml').read_bytes()==(sequence_output/'scaffold.yaml').read_bytes()
    assert result['official_parser_passed']


@pytest.mark.parametrize('mutation',['fixed_aa','cdr_length','range_bound','chain_id','extra_chain'])
def test_sequence_semantic_tampering_rejected(sequence_output,moldir,mutation):
    path=sequence_output/'scaffold.yaml'
    doc=yaml.safe_load(path.read_text(encoding='utf-8'))
    entity=doc['entities'][0]['protein']
    if mutation=='fixed_aa':
        entity['sequence']=entity['sequence'].replace('ACD','VCD')
    elif mutation=='cdr_length':
        entity['sequence']=entity['sequence'].replace('ACD,2,','ACD,4,')
    elif mutation=='range_bound':
        entity['sequence']=entity['sequence'].replace('2..5','2..6')
    elif mutation=='chain_id':
        entity['id']='X'
    else:
        doc['entities'].append({'protein':{'id':'X','sequence':'AV','msa':-1}})
    path.write_text(yaml.safe_dump(doc),encoding='utf-8')
    # Even an in-range random sample must not hide a changed range endpoint.
    with patch('numpy.random.randint',side_effect=lambda lo,hi:lo):
        with pytest.raises(ScaffoldError):
            validate_scaffold(sequence_output,moldir,run_cli=False)
    assert json.loads((sequence_output/'metadata.json').read_text())['validation_status']=='failed'


@pytest.mark.integration
def test_sequence_actual_boltzgen_check(sequence_output,moldir):
    result=validate_scaffold(sequence_output,moldir,run_cli=True)
    assert result['official_parser_passed'] and result['cli_check_passed']
    assert result['cli_exit_code']==0
    assert (sequence_output/'boltzgen_check.log').exists()
    assert (sequence_output/'checked/check.cif').exists()
    assert not (sequence_output/'scaffold.cif').exists()
    assert 'not a predicted structure' in result['coordinate_note']

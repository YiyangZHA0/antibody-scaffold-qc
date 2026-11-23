"""Execute the installed, unmodified BoltzGen parser and its real check command."""
from __future__ import annotations
from importlib.metadata import version
from hashlib import sha256
from pathlib import Path
from itertools import groupby
import json
import re
import subprocess
import sys
import numpy as np
import yaml
import gemmi
from .errors import ScaffoldError
from .output import dump_json, write_yaml_aliases
from .structure import one_letter

def validate_scaffold(directory, moldir, run_cli=True):
    directory=Path(directory).resolve()
    metadata_path=directory/'metadata.json'
    if not metadata_path.exists():
        raise ScaffoldError(f'No scaffold metadata in {directory}')
    metadata=json.loads(metadata_path.read_text(encoding='utf-8'))
    metadata['validation_status']='pending_boltzgen'
    dump_json(metadata_path,metadata)
    try:
        if metadata.get('scaffold_kind') == 'sequence':
            # A sequence scaffold is already a complete task, not a file fragment.
            (directory/'check.yaml').write_bytes((directory/'scaffold.yaml').read_bytes())
        else:
            (directory/'check.yaml').write_text(yaml.safe_dump({'entities':[{'file':{'path':'scaffold.yaml'}}]},sort_keys=False),encoding='utf-8')
        write_yaml_aliases(directory)
        return _validate_scaffold(directory,moldir,run_cli)
    except Exception as exc:
        metadata['validation_status']='failed'
        dump_json(metadata_path,metadata)
        dump_json(directory/'validation.json',{'official_parser_passed':False,'cli_check_passed':False,'error':str(exc)})
        raise


def _validate_scaffold(directory, moldir, run_cli=True):
    directory=Path(directory).resolve()
    moldir=Path(moldir).resolve()
    if not moldir.exists():
        raise ScaffoldError(f'BoltzGen molecule dictionary not found: {moldir}. Run scripts/fetch_boltzgen_mols.py or specify --moldir.')
    from boltzgen.data.parse.schema import YamlDesignParser
    from boltzgen.data.mol import load_canonicals
    import boltzgen.data.parse.schema as schema_module
    meta=json.loads((directory/'metadata.json').read_text(encoding='utf-8'))
    if meta.get('scaffold_kind') == 'structure':
        _validate_structure_insertions(directory,meta)
    mols=load_canonicals(moldir)
    parser=YamlDesignParser(moldir)
    parsed=parser.parse_yaml(directory/'check.yaml',mols,moldir)
    if meta.get('scaffold_kind') == 'sequence':
        return _validate_sequence_scaffold(directory,moldir,parsed,meta,schema_module,run_cli)
    if meta.get('scaffold_kind')!='structure':
        raise ScaffoldError('This validator accepts structural scaffold fragments only')
    expected_rows=[json.loads(line) for line in (directory/'residue_map.jsonl').read_text(encoding='utf-8').splitlines()]
    chain_results=[]
    parsed_rows=[]
    coordinate_structure=gemmi.read_structure(str(directory/'scaffold.cif'))
    actual_names=set(parsed.structure.chains['name'])
    if actual_names != {c['output_chain_id'] for c in meta['chains']}:
        raise ScaffoldError(f'BoltzGen chain identity changed unexpectedly: {actual_names}')
    for chain in parsed.structure.chains:
        name=str(chain['name']); start=int(chain['res_idx']); end=start+int(chain['res_num'])
        rows=[r for r in expected_rows if r['output_chain_id']==name and r['included']]
        actual_mask=parsed.design_info.res_design_mask[start:end].astype(bool)
        actual_groups=parsed.design_info.res_structure_groups[start:end]
        if not np.array_equal(actual_groups,np.where(actual_mask,0,2)):
            raise ScaffoldError(f'{name}: official BoltzGen structure groups disagree with intended mask')
        residues=parsed.structure.residues[start:end]
        expected_runs=[(flag,list(items)) for flag,items in groupby(rows,key=lambda r:r['designed'])]
        actual_runs=[(flag,list(indices)) for flag,indices in groupby(range(len(residues)),key=lambda i:bool(actual_mask[i]))]
        if len(expected_runs)!=len(actual_runs):
            raise ScaffoldError(f'{name}: official BoltzGen design mask disagrees with CDR map (run count)')
        loop_results=[]
        for (expected_design,original),(actual_design,indices) in zip(expected_runs,actual_runs):
            if expected_design!=actual_design:
                raise ScaffoldError(f'{name}: official BoltzGen design mask disagrees with CDR map')
            if expected_design:
                region=original[0]['region']
                spec=meta['loop_lengths'][f'{name}-{region}']
                bounds=[int(x) for x in str(spec).split('..')]
                if not bounds[0]<=len(indices)<=bounds[-1]:
                    raise ScaffoldError(f'{name}-{region}: parsed length {len(indices)} outside {spec}')
                loop_results.append({'cdr':region,'length_spec':spec,'sampled_length':len(indices)})
                for i in indices:
                    parsed_rows.append({'chain':name,'final_res_index':i+1,'region':region,'designed':True,'source_res_index':None,'aa':None})
            else:
                if len(indices)!=len(original):
                    raise ScaffoldError(f'{name}: official BoltzGen design mask disagrees with CDR map (fixed run length)')
                for i,row in zip(indices,original):
                    residue=residues[i]
                    if one_letter(str(residue['name']))!=row['aa']:
                        raise ScaffoldError(f'{name}: BoltzGen changed fixed framework sequence')
                    if not residue['is_present']:
                        raise ScaffoldError(f'{name}: fixed framework residue lacks coordinates')
                    original_ca=coordinate_structure[0][name][row['boltzgen_res_index']-1].find_atom('CA','*').pos
                    atoms=parsed.structure.atoms[int(residue['atom_idx']):int(residue['atom_idx'])+int(residue['atom_num'])]
                    parsed_ca=atoms[atoms['name']=='CA']['coords']
                    if parsed_ca.shape!=(1,3) or not np.allclose(parsed_ca[0],[original_ca.x,original_ca.y,original_ca.z],atol=0.001,rtol=0):
                        raise ScaffoldError(f'{name}:{i+1}: fixed framework CA coordinates changed')
                    parsed_rows.append({'chain':name,'final_res_index':i+1,'region':row['region'],'designed':False,'source_res_index':row['boltzgen_res_index'],'aa':row['aa']})
        if not np.array_equal(residues['res_idx'],np.arange(len(residues))):
            raise ScaffoldError(f'{name}: reset_res_index did not produce consecutive positions')
        chain_results.append({'chain':name,'residues':len(residues),'designed_residues':int(actual_mask.sum()),'fixed_residues':int((~actual_mask).sum()),'mask_matches':True,'groups_match':True,'fixed_sequence_matches':True,'fixed_ca_coordinates_match':True,'reset_indices_match':True,'has_input_coordinates':True,'loops':loop_results})
    report={'boltzgen_version':version('boltzgen'),'parser_source_sha256':sha256(Path(schema_module.__file__).read_bytes()).hexdigest(),'official_parser_passed':True,'chains':chain_results,'cli_check_passed':False,'scheme_preset':meta.get('scheme_preset','imgt'),'numbering_scheme':meta['numbering_scheme'],'cdr_definition':meta['cdr_definition'],'cdr_numbering_scheme':meta.get('cdr_numbering_scheme','imgt')}
    report['scaffold_kind']='structure'
    report['coordinate_note']='Framework input coordinates verified. Inserted design residues have no experimentally observed coordinates.'
    (directory/'parsed_residue_map.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in parsed_rows),encoding='utf-8')
    names=['scaffold.yaml','scaffold.yml','design.yaml','residue_map.jsonl','scaffold.cif','parsed_residue_map.jsonl']
    report['artifact_sha256']={name:sha256((directory/name).read_bytes()).hexdigest() for name in names}
    if run_cli:
        _run_official_check(directory,moldir,report)
    dump_json(directory/'validation.json',report)
    meta['validation_status']='boltzgen_cli_and_masks_passed' if run_cli else 'boltzgen_parser_and_masks_passed'
    dump_json(directory/'metadata.json',meta)
    return report


def _length_bounds(value):
    """Canonical semantic length bounds, without sampling a designed loop."""
    if isinstance(value,bool) or not re.fullmatch(r'[1-9][0-9]*(?:\.\.[1-9][0-9]*)?',str(value)):
        raise ScaffoldError(f'Invalid design insertion length specification: {value!r}')
    bounds=[int(part) for part in str(value).split('..')]
    if bounds[0]>bounds[-1]:
        raise ScaffoldError(f'Invalid design insertion length specification: {value!r}')
    return bounds[0],bounds[-1]


def _validate_structure_insertions(directory,meta):
    """Check complete length bounds, chain, position and multiplicity before sampling."""
    expected={}
    for chain in meta['chains']:
        name=chain['output_chain_id']
        for cdr in chain['cdrs']:
            if not cdr['designed']:
                continue
            indices=cdr['indices']
            if not indices:
                raise ScaffoldError(f"{name}-{cdr['name']}: designed CDR has no original residue indices")
            key=(name,min(indices))
            if key in expected:
                raise ScaffoldError(f'{name}: multiple designed CDRs share one insertion position')
            expected[key]=_length_bounds(meta['loop_lengths'][f"{name}-{cdr['name']}"])
    doc=yaml.safe_load((directory/'scaffold.yaml').read_text(encoding='utf-8'))
    items=doc.get('design_insertions',[])
    if not isinstance(items,list):
        raise ScaffoldError('Structural design_insertions must be a list')
    actual={}
    for item in items:
        if not isinstance(item,dict) or set(item)!={'insertion'} or not isinstance(item['insertion'],dict):
            raise ScaffoldError('Malformed structural design insertion')
        insertion=item['insertion']
        name=insertion.get('id')
        index=insertion.get('res_index')
        if not isinstance(name,str) or isinstance(index,bool) or not re.fullmatch(r'[1-9][0-9]*',str(index)):
            raise ScaffoldError('Structural design insertion needs one chain id and positive residue index')
        key=(name,int(index))
        if key in actual:
            raise ScaffoldError(f'{name}:{index}: duplicate design insertion')
        actual[key]=_length_bounds(insertion.get('num_residues'))
    if actual!=expected:
        raise ScaffoldError('Structural design insertions differ from expected CDR chains, positions or length specifications')


def _run_official_check(directory,moldir,report):
    """Run the real CLI identically for structural and sequence scaffolds."""
    previous=directory/'checked/check.cif'
    if previous.exists():
        previous.unlink()
    command=[sys.executable,'-m','boltzgen.cli.boltzgen','check',str(directory/'check.yaml'),'--moldir',str(moldir),'--output',str(directory/'checked')]
    process=subprocess.run(command,capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=180)
    (directory/'boltzgen_check.log').write_text(process.stdout+'\n'+process.stderr,encoding='utf-8')
    report['cli_command']=command
    report['cli_exit_code']=process.returncode
    if process.returncode!=0 or not previous.exists():
        dump_json(directory/'validation.json',report)
        raise ScaffoldError(f'Official boltzgen check failed; see {directory / "boltzgen_check.log"}')
    report['cli_check_passed']=True


def _sequence_tokens(sequence):
    """Read the explicit amino-acid/length contract without sampling ranges."""
    compact=re.sub(r'[,\s]', '', str(sequence))
    tokens=re.findall(r'\d+\.\.\d+|\d+|[A-Z]',compact)
    if ''.join(tokens)!=compact:
        raise ScaffoldError('Sequence scaffold contains unsupported sequence syntax')
    return tokens


def _validate_sequence_scaffold(directory,moldir,parsed,meta,schema_module,run_cli):
    """Check unmodified BoltzGen output against the original residue/loop map.

    BoltzGen constructs zero-coordinate placeholders for sequence entities. They
    are neither experimental coordinates nor a structure prediction.
    """
    expected_rows=[json.loads(line) for line in (directory/'residue_map.jsonl').read_text(encoding='utf-8').splitlines()]
    document=yaml.safe_load((directory/'scaffold.yaml').read_text(encoding='utf-8'))
    entities=document.get('entities',[])
    if not entities or any(set(entity)!={'protein'} for entity in entities):
        raise ScaffoldError('Sequence scaffold must contain protein entities only')
    expected_names={c['output_chain_id'] for c in meta['chains']}
    entity_names=[entity['protein'].get('id') for entity in entities]
    if any(not isinstance(name,str) for name in entity_names) or len(entity_names)!=len(set(entity_names)) or set(entity_names)!=expected_names:
        raise ScaffoldError('Sequence scaffold chain identity changed unexpectedly')
    entity_by_name={entity['protein']['id']:entity['protein'] for entity in entities}
    actual_names=[str(chain['name']) for chain in parsed.structure.chains]
    if set(actual_names)!=expected_names or len(actual_names)!=len(expected_names):
        raise ScaffoldError(f'BoltzGen chain identity changed unexpectedly: {actual_names}')
    chain_results=[]
    parsed_rows=[]
    for chain in parsed.structure.chains:
        name=str(chain['name']); start=int(chain['res_idx']); end=start+int(chain['res_num'])
        rows=[r for r in expected_rows if r['output_chain_id']==name and r['included']]
        expected_runs=[(flag,list(items)) for flag,items in groupby(rows,key=lambda r:r['designed'])]
        expected_tokens=[]
        for designed,original in expected_runs:
            if designed:
                if len({row['region'] for row in original})!=1:
                    raise ScaffoldError(f'{name}: consecutive designed regions are ambiguous')
                expected_tokens.append(str(meta['loop_lengths'][f"{name}-{original[0]['region']}"]))
            else:
                expected_tokens.extend(row['aa'] for row in original)
        if _sequence_tokens(entity_by_name[name].get('sequence',''))!=expected_tokens:
            raise ScaffoldError(f'{name}: sequence YAML differs from expected fixed sequence or CDR length specification')
        actual_mask=parsed.design_info.res_design_mask[start:end].astype(bool)
        actual_groups=parsed.design_info.res_structure_groups[start:end]
        if not np.array_equal(actual_groups,np.zeros(len(actual_groups),dtype=actual_groups.dtype)):
            raise ScaffoldError(f'{name}: sequence scaffold must not constrain structure groups')
        residues=parsed.structure.residues[start:end]
        actual_runs=[(flag,list(indices)) for flag,indices in groupby(range(len(residues)),key=lambda i:bool(actual_mask[i]))]
        if len(expected_runs)!=len(actual_runs):
            raise ScaffoldError(f'{name}: official BoltzGen design mask disagrees with CDR map (run count)')
        loop_results=[]
        for (expected_design,original),(actual_design,indices) in zip(expected_runs,actual_runs):
            if expected_design!=actual_design:
                raise ScaffoldError(f'{name}: official BoltzGen design mask disagrees with CDR map')
            if expected_design:
                region=original[0]['region']
                spec=meta['loop_lengths'][f'{name}-{region}']
                bounds=[int(x) for x in str(spec).split('..')]
                if not bounds[0]<=len(indices)<=bounds[-1]:
                    raise ScaffoldError(f'{name}-{region}: parsed length {len(indices)} outside {spec}')
                loop_results.append({'cdr':region,'length_spec':spec,'sampled_length':len(indices)})
                parsed_rows.extend({'chain':name,'final_res_index':i+1,'region':region,'designed':True,'source_res_index':None,'aa':None} for i in indices)
            else:
                if len(indices)!=len(original):
                    raise ScaffoldError(f'{name}: official BoltzGen design mask disagrees with CDR map (fixed run length)')
                for i,row in zip(indices,original):
                    if one_letter(str(residues[i]['name']))!=row['aa']:
                        raise ScaffoldError(f'{name}: BoltzGen changed fixed framework sequence')
                    parsed_rows.append({'chain':name,'final_res_index':i+1,'region':row['region'],'designed':False,'source_res_index':row['boltzgen_res_index'],'aa':row['aa']})
        if not np.array_equal(residues['res_idx'],np.arange(len(residues))):
            raise ScaffoldError(f'{name}: sequence positions are not consecutive')
        # The official sequence parser currently emits all-zero coordinates.
        for residue in residues:
            atoms=parsed.structure.atoms[int(residue['atom_idx']):int(residue['atom_idx'])+int(residue['atom_num'])]
            if not np.all(atoms['coords']==0):
                raise ScaffoldError(f'{name}: sequence entity unexpectedly contains non-placeholder coordinates')
        chain_results.append({'chain':name,'residues':len(residues),'designed_residues':int(actual_mask.sum()),'fixed_residues':int((~actual_mask).sum()),'mask_matches':True,'groups_match':True,'fixed_sequence_matches':True,'length_specifications_match':True,'reset_indices_match':True,'has_input_coordinates':False,'loops':loop_results})
    report={'boltzgen_version':version('boltzgen'),'parser_source_sha256':sha256(Path(schema_module.__file__).read_bytes()).hexdigest(),'official_parser_passed':True,'chains':chain_results,'cli_check_passed':False,'scheme_preset':meta.get('scheme_preset','imgt'),'numbering_scheme':meta['numbering_scheme'],'cdr_definition':meta['cdr_definition'],'cdr_numbering_scheme':meta.get('cdr_numbering_scheme','imgt'),'scaffold_kind':'sequence','has_input_coordinates':False,'coordinate_note':'No input three-dimensional coordinates. BoltzGen check produces zero-coordinate placeholders, not a predicted structure.'}
    (directory/'parsed_residue_map.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in parsed_rows),encoding='utf-8')
    names=['scaffold.yaml','scaffold.yml','design.yaml','check.yaml','residue_map.jsonl','parsed_residue_map.jsonl']
    report['artifact_sha256']={name:sha256((directory/name).read_bytes()).hexdigest() for name in names}
    if run_cli:
        _run_official_check(directory,moldir,report)
    dump_json(directory/'validation.json',report)
    meta['validation_status']='boltzgen_cli_and_masks_passed' if run_cli else 'boltzgen_parser_and_masks_passed'
    dump_json(directory/'metadata.json',meta)
    return report

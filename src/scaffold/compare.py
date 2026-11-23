"""Compare coordinates, numbering policies and residue masks as well as YAML."""
from hashlib import sha256
from pathlib import Path
import json
import gemmi
import numpy as np
import yaml
from .errors import ScaffoldError

MAP_FIELDS=('output_chain_id','boltzgen_res_index','aa','imgt_label','region','included','designed','structure_group')
SCHEME_FIELDS=('numbering_scheme','numbering_label','cdr_definition','cdr_numbering_scheme','cdr_label')


def signature(directory):
    directory=Path(directory)
    text=(directory/'scaffold.yaml').read_text(encoding='utf-8')
    doc=yaml.safe_load(text)
    body='\n'.join(line for line in text.splitlines() if line.strip() and not line.startswith('#'))
    rows=[json.loads(line) for line in (directory/'residue_map.jsonl').read_text(encoding='utf-8').splitlines()]
    metadata=json.loads((directory/'metadata.json').read_text(encoding='utf-8'))
    kind=metadata.get('scaffold_kind','structure')
    mapped=[{key:row[key] for key in MAP_FIELDS} for row in rows]
    selected_mapping=[{key:row.get(key) for key in SCHEME_FIELDS} for row in rows]
    atoms=[]
    if kind=='structure':
        for chain in gemmi.read_structure(str(directory/'scaffold.cif'))[0]:
            for residue in chain:
                for atom in residue:
                    atoms.append((chain.name,residue.seqid.num,residue.name,atom.name,atom.pos.x,atom.pos.y,atom.pos.z))
    validation=json.loads((directory/'validation.json').read_text(encoding='utf-8'))
    if not validation.get('official_parser_passed') or not validation.get('cli_check_passed'):
        raise ScaffoldError(f'{directory}: official parser/check has not passed')
    for name,digest in validation['artifact_sha256'].items():
        if sha256((directory/name).read_bytes()).hexdigest()!=digest:
            raise ScaffoldError(f'{directory/name}: changed since validation')
    policy={key:validation.get(key) for key in ('scheme_preset','numbering_scheme','cdr_definition','cdr_numbering_scheme')}
    semantic_fields=('output_chain_id','aa','imgt_label','numbering_label','cdr_label','region','designed')
    sequence_design={'residues':[{key:r.get(key) for key in semantic_fields} for r in rows if r['included']],
                     'loop_lengths':metadata['loop_lengths']}
    return {'scaffold_kind':kind,'sequence_design':sequence_design,'yaml':doc,'yaml_body_text':body,'imgt_map':mapped,'selected_numbering_map':selected_mapping,'scheme_policy':policy,'atoms':sorted(atoms),'parser_chains':validation['chains']}


def compare_coordinates(reference, other):
    identities_match=[r[:4] for r in reference]==[r[:4] for r in other]
    if not identities_match or not reference:
        return {'passed':False,'atom_identities_match':identities_match}
    delta=np.abs(np.array([r[4:] for r in reference])-np.array([r[4:] for r in other]))
    # PDB XYZ fields retain 3 decimal places; raw RCSB CIF can retain 5.
    # Keep original coordinates and report the measured error instead of
    # rounding both sides (ties may round in opposite directions).
    maximum=float(delta.max())
    tolerance=0.00050001
    return {'passed':maximum<=tolerance,'atom_identities_match':True,'atoms_compared':len(reference),'max_axis_difference_angstrom':maximum,'max_atom_distance_angstrom':float(np.linalg.norm(delta,axis=1).max()),'axis_tolerance_angstrom':tolerance}


def compare_outputs(directories):
    if len(directories)<2:
        raise ScaffoldError('Compare requires at least two conversion output directories')
    snapshots=[signature(path) for path in directories]
    kinds={s['scaffold_kind'] for s in snapshots}
    if 'sequence' in kinds:
        # Different YAML schemas can express the same amino-acid design. There
        # is no coordinate equivalence claim when any input has no structure.
        fields=('sequence_design','scheme_policy')
        if kinds=={'sequence'}:
            fields+=('yaml','yaml_body_text')
        checks={key:all(s[key]==snapshots[0][key] for s in snapshots[1:]) for key in fields}
        def cdr_design(snapshot):
            spec = snapshot['sequence_design']
            return {'residues': [r for r in spec['residues'] if r['designed']], 'loop_lengths': spec['loop_lengths']}
        checks['cdr_design']=all(cdr_design(s)==cdr_design(snapshots[0]) for s in snapshots[1:])
        return {'passed':all(checks.values()),'directories':[str(Path(p).resolve()) for p in directories],
                'comparison_scope':'sequence_design','checks':checks,'coordinates_compared':False,
                'scaffold_kinds':[s['scaffold_kind'] for s in snapshots],
                'comparison_note':'Compare retained residue identities, numbering/CDR membership and designed lengths. Sequence input has no framework coordinates; cross-mode YAML bodies and structure constraints are not equivalent.'}
    fields=('yaml','yaml_body_text','imgt_map','selected_numbering_map','scheme_policy')
    checks={key:all(s[key]==snapshots[0][key] for s in snapshots[1:]) for key in fields}
    coordinates=[compare_coordinates(snapshots[0]['atoms'],s['atoms']) for s in snapshots[1:]]
    checks['coordinates_at_pdb_precision']=all(c['passed'] for c in coordinates)
    # With length ranges, the parser samples new lengths per invocation. Compare
    # constraints and input masks, not equality of independent random samples.
    fixed=all(isinstance(v['insertion']['num_residues'],int) for v in snapshots[0]['yaml']['design_insertions'])
    if fixed:
        checks['parser_chains']=all(s['parser_chains']==snapshots[0]['parser_chains'] for s in snapshots[1:])
    return {'passed':all(checks.values()),'directories':[str(Path(p).resolve()) for p in directories],'checks':checks,'coordinate_comparisons':coordinates,'yaml_body_sha256':[sha256(s['yaml_body_text'].encode()).hexdigest() for s in snapshots],'comparison_note':'Input provenance comments differ intentionally. YAML body, selected numbering/CDR policy and mappings must agree exactly; coordinates must agree within PDB decimal precision. Fixed-length parser masks and lengths must also agree.'}

from __future__ import annotations
from datetime import datetime,timezone
from hashlib import sha256
from pathlib import Path
import csv
import json
import re
import yaml
from .errors import ScaffoldError

def ranges(indices):
    values=sorted(set(indices))
    groups=[]
    for i in values:
        if groups and i==groups[-1][-1]+1:
            groups[-1].append(i)
        else:
            groups.append([i])
    return ','.join(str(g[0]) if len(g)==1 else f'{g[0]}..{g[-1]}' for g in groups)

def dump_json(path,obj):
    Path(path).write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')

def write_rows(directory, rows, stem='residue_map'):
    directory=Path(directory)
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with (directory/f'{stem}.csv').open('w',newline='',encoding='utf-8-sig') as out:
        writer=csv.DictWriter(out,fieldnames=fields)
        writer.writeheader();writer.writerows(rows)
    with (directory/f'{stem}.jsonl').open('w',encoding='utf-8') as out:
        for row in rows:
            out.write(json.dumps(row,ensure_ascii=False)+'\n')

def write_annotation(out,annotation,records,input_path):
    out=Path(out)
    out.mkdir(parents=True,exist_ok=True)
    rows=[row for d in annotation['domains'] for row in d.rows]
    write_rows(out,rows,'numbering')
    with (out/'variable_domains.fasta').open('w',encoding='utf-8') as handle:
        for d in annotation['domains']:
            role='VH' if d.chain_type=='H' else 'VL'
            handle.write(f'>{d.record_id}__{role} chain_type={d.chain_type} original={d.start+1}..{d.end}\n{d.sequence}\n')
    metadata={k:v for k,v in annotation.items() if k!='domains'}
    metadata['domains']=[d.to_dict() for d in annotation['domains']]
    metadata['input_file']=str(Path(input_path).resolve())
    metadata['input_sha256']=sha256(Path(input_path).read_bytes()).hexdigest()
    metadata['full_input_sequences']=records
    metadata['created_utc']=datetime.now(timezone.utc).isoformat()
    dump_json(out/'annotation.json',metadata)
    return metadata

def parse_loop_lengths(specs, summaries):
    """Length values are design lengths, never IMGT residue labels."""
    lengths={f"{s['output_chain_id']}-{c['name']}":c['length'] for s in summaries for c in s['cdrs'] if c['designed']}
    seen=set()
    for spec in specs or []:
        match=re.fullmatch(r'(H|L)-(CDR[123])=([1-9][0-9]*)(?:\.\.([1-9][0-9]*))?',spec)
        if not match:
            raise ScaffoldError(f'Invalid --cdr-length {spec!r}; use H-CDR3=8..18 or L-CDR1=6')
        chain,cdr,lower,upper=match.groups()
        key=f'{chain}-{cdr}'
        if key not in lengths or key in seen:
            raise ScaffoldError(f'{key}: length override must name a selected CDR exactly once')
        lower=int(lower); upper=int(upper or lower)
        if not 1 <= lower <= upper <= 100:
            raise ScaffoldError(f'{key}: require 1 <= minimum <= maximum <= 100')
        lengths[key]=lower if lower==upper else f'{lower}..{upper}'
        seen.add(key)
    return lengths


def write_scaffold(out,metadata,rows,summaries,warnings,provenance,loop_specs=None,name=None):
    out=Path(out)
    lengths=parse_loop_lengths(loop_specs,summaries)
    document={'path':'scaffold.cif','include':[],'design':[],'structure_groups':[],'exclude':[],'design_insertions':[],'reset_res_index':[]}
    for s in summaries:
        cid=s['output_chain_id']
        document['include'].append({'chain':{'id':cid,'res_index':f"1..{s['length']}"}})
        document['structure_groups'].append({'group':{'id':cid,'visibility':2}})
    for s in summaries:
        cid=s['output_chain_id']
        designed=ranges(r['boltzgen_res_index'] for r in rows if r['output_chain_id']==cid and r['designed'])
        if designed:
            document['design'].append({'chain':{'id':cid,'res_index':designed}})
            document['structure_groups'].append({'group':{'id':cid,'res_index':designed,'visibility':0}})
        # Indices refer to the original normalized CIF, before any deletion/insertion.
        # BoltzGen applies all insertions before its exclusion mask; do not subtract
        # preceding excluded residues from these positions.
        for cdr in s['cdrs']:
            if not cdr['designed']:
                continue
            document['exclude'].append({'chain':{'id':cid,'res_index':ranges(cdr['indices'])}})
            document['design_insertions'].append({'insertion':{'id':cid,'res_index':min(cdr['indices']),'num_residues':lengths[f"{cid}-{cdr['name']}"]}})
        document['reset_res_index'].append({'chain':{'id':cid}})
    if not document['design']:
        raise ScaffoldError('No design residues selected')
    display_name=' '.join((name or out.name).splitlines())
    head=[
        f'# {display_name}',
        f"# Source molecule: {metadata['molecule_type']}",
        f"# Included domains: {'one VHH domain' if len(summaries)==1 else 'VH + VL (Fv scope; constant domains omitted)'}",
        f"# Scheme preset: {metadata['scheme_preset']}",
        f"# Numbering scheme: {metadata['numbering_scheme']}",
        f"# CDR definition: {metadata['cdr_definition']}",
        f"# CDR boundary numbering: {metadata['cdr_numbering_scheme']}",
        f"# CDR boundary intervals: {metadata['cdr_ranges']}",
        f"# Annotation tool/version: {metadata['numbering_engine']} {metadata['numbering_version']}",
        '# Design policy: replace each selected complete CDR; preserve framework sequence and structure',
        f'# Design lengths: {lengths}',
        '# Index basis: normalized scaffold.cif parsed-chain 1-based positions; NOT antibody numbering labels',
        '# All include/design/exclude/insertion indices refer to CIF before edits; reset afterward',
        '# Mapping: residue_map.csv / residue_map.jsonl; IMGT labels retained as reference',
        f"# Input file: {metadata['input_file']}",
        f"# Input SHA256: {metadata['input_sha256']}",
        f"# Source structure: {provenance['structure_file']}",
        f"# Source structure SHA256: {provenance['structure_sha256']}",
        f"# Structure resolution: {provenance['resolution_method']}",
    ]
    for s in summaries:
        cid=s['output_chain_id']
        head.extend([f"# {s['role']} chain: source={s['source_chain_id']}; output_label={cid}; type={s['chain_type']}",f"# {cid} full source sequence: {s['full_source_sequence']}",f"# {cid} numbered variable sequence: {s['full_variable_sequence']}",f"# {cid} included sequence: {s['selected_sequence']}"])
        for cdr in s['cdrs']:
            head.append(f"# {cid}-{cdr['name']}: {metadata['numbering_scheme']}={','.join(cdr['numbering_labels'])}; boundary_labels={','.join(cdr['cdr_numbering_labels'])}; IMGT_reference={','.join(cdr['imgt_labels'])}; res_index={ranges(cdr['indices'])}; sequence={cdr['sequence']}; length={cdr['length']}; design={cdr['designed']}")
        for region in ('FR1','FR2','FR3','FR4'):
            region_rows=[r for r in rows if r['output_chain_id']==cid and r['region']==region]
            head.append(f"# {cid}-{region}: {metadata['numbering_scheme']}={','.join(r['numbering_label'] for r in region_rows)}; IMGT_reference={','.join(r['imgt_label'] for r in region_rows)}; res_index={ranges(r['boltzgen_res_index'] for r in region_rows if r['included'])}; sequence={''.join(r['aa'] for r in region_rows)}; fixed=True")
    for warning in metadata.get('warnings',[])+warnings:
        head.append('# Note: '+warning.replace('\n',' '))
    (out/'scaffold.yaml').write_text('\n'.join(head)+'\n\n'+yaml.safe_dump(document,sort_keys=False,width=1000),encoding='utf-8')
    write_rows(out,rows)
    full={**metadata,'name':display_name,'chains':summaries,'loop_lengths':lengths,'design_policy':'replace_complete_selected_cdrs','warnings':metadata.get('warnings',[])+warnings,'structure_provenance':provenance,'scaffold_kind':'structure','has_structure':True,'yaml_kind':'boltzgen_file_fragment','scaffold_output_scope':'vhh' if len(summaries)==1 else 'fv','validation_status':'pending_boltzgen'}
    dump_json(out/'metadata.json',full)
    # Scaffold is not a main task; this wrapper is only a parser/check entry point.
    (out/'check.yaml').write_text(yaml.safe_dump({'entities':[{'file':{'path':'scaffold.yaml'}}]},sort_keys=False),encoding='utf-8')
    return full


def write_yaml_aliases(out):
    out=Path(out)
    (out/'scaffold.yml').write_bytes((out/'scaffold.yaml').read_bytes())
    # Uniform main-task entry point for either scaffold kind.
    (out/'design.yaml').write_bytes((out/'check.yaml').read_bytes())

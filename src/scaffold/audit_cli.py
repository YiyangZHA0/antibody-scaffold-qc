"""CLI orchestration for read-only single/paired/mixed-batch IMGT audits."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
from time import perf_counter
from uuid import uuid4
import json
import re
from .audit import audit_file, write_audit
from .errors import ScaffoldError
from .numbering import Numberer


def add_audit_parser(sub):
    p=sub.add_parser('audit',help='Audit IMGT FR sequence/coordinate risks without modifying inputs or generating YAML')
    p.add_argument('inputs',nargs='+',type=Path,help='One/many FASTA/PDB/CIF files, directories or globs')
    p.add_argument('--out',type=Path,help='New/empty report folder; default outputs/audit_<UTC timestamp>_<unique>')
    p.add_argument('--recursive',action='store_true')
    p.add_argument('--extensions',help='Directory/glob filter, e.g. pdb,cif,fasta')
    p.add_argument('--sequence',type=Path,help='Full construct FASTA to pair with ONE structure input')
    p.add_argument('--reference',type=Path,help='Intended complete construct reference FASTA for ONE input; differences are reference-relative')
    p.add_argument('--records',help='Comma-separated FASTA record IDs (primary FASTA or --sequence)')
    p.add_argument('--chain',action='append',default=[],help='Repeatable exact structure chain ID; default audits every recognized antibody domain')
    p.add_argument('--heavy-chain',help='Select and verify a heavy structure chain')
    p.add_argument('--light-chain',help='Select and verify a light structure chain')
    p.add_argument('--vhh-chain',help='Select a heavy variable-domain structure chain; origin is user-specified')
    p.add_argument('--chain-namespace',choices=['auth','label'],default='auth')
    p.add_argument('--model-index',type=int,default=1)
    p.add_argument('--scheme',choices=['imgt'],default='imgt',help='Integrity rules currently defined for IMGT only')
    p.add_argument('--engine',choices=['anarcii','anarci'],default='anarcii')
    p.add_argument('--ncpu',type=int,default=2)


def run_audit(args, numberer=None):
    from .batch import expand_inputs
    paths=expand_inputs(args.inputs,args.recursive,args.extensions)
    if args.ncpu<1 or args.model_index<1:
        raise ScaffoldError('--ncpu and --model-index must be positive')
    if (args.sequence or args.reference) and len(paths)!=1:
        raise ScaffoldError('--sequence/--reference apply to one input only; audit each intended construct/reference pair separately')
    if args.sequence and paths[0].suffix.lower() not in {'.pdb','.cif','.mmcif'}:
        raise ScaffoldError('--sequence requires a PDB/CIF primary input')
    if args.heavy_chain and args.vhh_chain:
        raise ScaffoldError('Choose --heavy-chain or --vhh-chain')
    if args.vhh_chain and args.light_chain:
        raise ScaffoldError('--vhh-chain cannot be paired with --light-chain')
    if (args.heavy_chain or args.vhh_chain) and (args.heavy_chain or args.vhh_chain)==args.light_chain:
        raise ScaffoldError('Heavy and light chain IDs must differ')
    selected=list(dict.fromkeys(args.chain+[x for x in (args.heavy_chain,args.light_chain,args.vhh_chain) if x])) or None
    records=args.records.split(',') if args.records else None
    if selected and any(p.suffix.lower() in {'.fa','.faa','.fasta'} for p in paths):
        raise ScaffoldError('Chain options apply to structure files; use --records for FASTA')
    if any(p.suffix.lower() in {'.pdb','.cif','.mmcif'} for p in paths) and records and not args.sequence:
        raise ScaffoldError('--records applies to FASTA input or the paired --sequence FASTA')
    root=(args.out or Path.cwd()/'outputs'/f'audit_{datetime.now(timezone.utc):%Y%m%d_%H%M%S}UTC_{uuid4().hex[:8]}').resolve()
    if root.exists() and (not root.is_dir() or any(root.iterdir())):
        raise ScaffoldError(f'Output directory is not empty: {root}')
    root.mkdir(parents=True,exist_ok=True)
    numberer=numberer or Numberer(args.engine,args.ncpu,'imgt')
    start=perf_counter();results=[]
    for i,path in enumerate(paths,1):
        label=re.sub(r'[^A-Za-z0-9_.-]','_',path.stem)[:45] or 'input'
        folder=root if len(paths)==1 else root/f'{i:03d}_{label}_{path.suffix[1:]}'
        before=sha256(path.read_bytes()).hexdigest()
        try:
            report=audit_file(path,numberer,args.sequence,args.reference,selected,args.chain_namespace,args.model_index,records)
            for expected,key in [('H',args.heavy_chain or args.vhh_chain),('L',args.light_chain)]:
                if not key:
                    continue
                selected_domains=[d for d in report['domains'] if d['structure_chain_'+args.chain_namespace]==key]
                if not selected_domains or any((d['chain_type']=='H')!=(expected=='H') for d in selected_domains):
                    raise ScaffoldError(f'Chain {key} does not unambiguously contain the requested {expected} domain type')
            report['input_unchanged']=sha256(path.read_bytes()).hexdigest()==before
            assert report['input_unchanged']
            write_audit(report,folder)
            results.append({'input':str(path),'status':report['status'],'domains':report['domain_count'],
                            'risk_counts':report['risk_counts'],'output':str(folder)})
        except (ScaffoldError,ValueError,OSError,RuntimeError) as exc:
            folder.mkdir(parents=True,exist_ok=True)
            failure={'input':str(path),'status':'failed','error':str(exc),'output':str(folder)}
            (folder/'failure.json').write_text(json.dumps(failure,ensure_ascii=False,indent=2),encoding='utf-8')
            results.append(failure)
        (root/'progress.json').write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding='utf-8')
    result={'status':'audit_completed','files':len(paths),'failed_files':sum(r['status']=='failed' for r in results),
            'unassessable_files':sum(r['status']=='not_assessable' for r in results),'elapsed_seconds':perf_counter()-start,
            'output':str(root),'results':results,'passed':not any(r['status']=='failed' for r in results),
            'passed_meaning':'Execution completed; this is not a declaration that the input is complete or risk-free.'}
    (root/'summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    if len(paths)>1:
        text=['# IMGT audit batch','',result['passed_meaning'],'','| Input | Status | Domains | Report |','|---|---|---:|---|']
        for item in results:
            relative=Path(item['output']).relative_to(root).as_posix()
            text.append(f"| {Path(item['input']).name} | {item['status']} | {item.get('domains',0)} | [{relative}]({relative}) |")
        (root/'README_audit_batch.md').write_text('\n'.join(text)+'\n',encoding='utf-8-sig')
    return result

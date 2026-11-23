"""Run real/synthetic FR-audit regressions with externally specified expectations.

This script uses the previously curated local library as test data, and writes
input snapshots, per-case reports, timing and explicit assertions. Synthetic
deletions have known construction provenance; expected missing coordinates are
from original depositions/prior independent structure checks, not from the new
auditor. It is not an independent numbering gold-standard benchmark.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter

from scaffold.audit import audit_file, write_audit
from scaffold.numbering import Numberer


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def dump(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')


def fasta_read(path):
    records={}; current=None
    for line in path.read_text(encoding='utf-8-sig').splitlines():
        if line.startswith('>'):
            current=line[1:].split()[0];records[current]=''
        elif line.strip():
            records[current]+=line.strip()
    return records


def all_events(report):
    return report['risks']+[r for d in report['domains'] for r in d['risks']]


def selected(report, code):
    return [r for r in all_events(report) if r['code']==code]


def check(report, expected):
    checks=[]
    def add(name, actual, wanted):
        checks.append({'check':name,'actual':actual,'expected':wanted,'passed':actual==wanted})
    add('domain_count',report['domain_count'],expected['domains'])
    if 'source_auth_chains' in expected:
        add('source_auth_chains',sorted(d['structure_chain_auth'] for d in report['domains']),sorted(expected['source_auth_chains']))
    if 'source_label_chains' in expected:
        add('source_label_chains',sorted(d['structure_chain_label'] for d in report['domains']),sorted(expected['source_label_chains']))
    if 'FR1_missing_CA' in expected:
        add('FR1_missing_CA',sum(d['regions']['FR1']['missing_positive_CA_count'] or 0 for d in report['domains']),expected['FR1_missing_CA'])
    if 'all_rows_missing_CA' in expected:
        add('all_rows_missing_CA',sum(r.get('has_positive_ca') is False for d in report['domains'] for r in d['rows']),expected['all_rows_missing_CA'])
    if 'zero_occupancy_rows' in expected:
        add('zero_occupancy_rows',sum(r['coordinate_state']=='zero_occupancy' for d in report['domains'] for r in d['rows']),expected['zero_occupancy_rows'])
    if 'uncertain_mapping_rows' in expected:
        add('uncertain_mapping_rows',sum(r['coordinate_state']=='mapping_uncertain' for d in report['domains'] for r in d['rows']),expected['uncertain_mapping_rows'])
    if 'FR_missing_CA' in expected:
        add('FR_missing_CA',sum(fr['missing_positive_CA_count'] or 0 for d in report['domains'] for fr in d['regions'].values()),expected['FR_missing_CA'])
    for code in expected.get('required_codes',[]):
        add('code_present:'+code,bool(selected(report,code)),True)
    for code in expected.get('absent_codes',[]):
        add('code_absent:'+code,bool(selected(report,code)),False)
    if expected.get('coordinate_unassessed'):
        add('coordinate_unassessed',all(not fr['coordinate_assessed'] and fr['missing_positive_CA_count'] is None for d in report['domains'] for fr in d['regions'].values()),True)
    if 'reference_deleted_residues' in expected:
        changes=selected(report,'reference_terminal_deletion')+selected(report,'reference_internal_deletion')
        add('reference_deleted_residues',sum(len(r['evidence']['reference_residues']) for r in changes),expected['reference_deleted_residues'])
        add('reference_certainty',all(r['certainty']=='reference_relative' for r in changes),True)
    if 'reference_deleted_regions' in expected:
        changes=selected(report,'reference_terminal_deletion')+selected(report,'reference_internal_deletion')
        add('reference_deleted_regions',sorted({r['region'] for r in changes}),sorted(expected['reference_deleted_regions']))
    if 'minimum_noncanonical' in expected:
        add('minimum_noncanonical',len(selected(report,'noncanonical_anchor'))>=expected['minimum_noncanonical'],True)
    return checks


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--library',type=Path,default=Path(r'D:\Downloads\VHH_Fab预制骨架结构汇总'))
    parser.add_argument('--curated',type=Path,default=Path(r'D:\Downloads\VHH_Fab骨架完整交付_20261119'))
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--ncpu',type=int,default=2)
    args=parser.parse_args()
    out=args.out.resolve()
    if out.exists() and any(out.iterdir()):
        raise SystemExit('Choose a new/empty output folder')
    out.mkdir(parents=True,exist_ok=True)
    inputs=out/'inputs';inputs.mkdir()
    source_manifest=[]

    def bundle(tid):
        hits=[p for p in args.curated.glob('*/*') if p.is_dir() and p.name.startswith(tid+'_')]
        assert len(hits)==1,(tid,hits)
        return hits[0]

    def archive(path, name):
        target=inputs/name
        shutil.copy2(path,target)
        source_manifest.append({'file':str(target.relative_to(out)),'source':str(path.resolve()),'sha256':digest(target),'synthetic':False})
        return target

    def generated(name, seq, source, change):
        target=inputs/name
        target.write_text('>H\n'+seq+'\n',encoding='ascii')
        source_manifest.append({'file':str(target.relative_to(out)),'source':str(source),'sha256':digest(target),'synthetic':True,'exact_construction':change})
        return target

    rf_source=bundle('08')/'sequence.fasta'
    rf=fasta_read(rf_source)['H']
    assert len(rf)==128 and rf.endswith('SS')
    full=archive(rf_source,'RF_full_128.fasta')
    short=generated('RF_minus_last_S.fasta',rf[:-1],rf_source,'Delete final sequence position128 S only')
    cut_n=generated('RF_minus_N5.fasta',rf[5:],rf_source,'Delete first5 residues only')
    cut_fr4=generated('RF_minus_FR4.fasta',rf[:-11],rf_source,'Delete exact final11 aa WGQGTLVTVSS only')
    fragment=generated('RF_first60.fasta',rf[:60],rf_source,'Keep first60 amino acids only')
    bad=generated('RF_unknown_X.fasta',rf[:20]+'X'+rf[21:],rf_source,'Replace linear position21 with X; no sequence deletion')
    both=inputs/'valid_and_invalid.fasta'
    both.write_text('>valid\n'+rf+'\n>invalid\n'+rf[:20]+'X'+rf[21:]+'\n',encoding='ascii')
    source_manifest.append({'file':str(both.relative_to(out)),'sha256':digest(both),'synthetic':True,'source':str(rf_source),'exact_construction':'One unchanged RF128 record plus a copy with position21 X'})
    original={tid:archive(next(args.library.rglob(tid+'_*.pdb')),f'{tid}_original.pdb') for tid in ('04','08','10','13','15','22')}
    cif_10=archive(next((bundle('10')/'sources/upstream').glob('*.cif')),'10_original.cif')
    cif_22=archive(next((bundle('22')/'sources/upstream').glob('*.cif')),'22_original.cif')
    full_10=archive(bundle('10')/'sequence.fasta','10_full_sequence.fasta')
    fab_20=archive(bundle('20')/'sequence.fasta','20_complete_Fab.fasta')
    atom_only=inputs/'04_atom_only_renumbered_41.pdb'
    lines=[]
    for line in original['04'].read_text(encoding='ascii').splitlines():
        if line.startswith(('SEQRES','REMARK')):
            continue
        if line.startswith(('ATOM  ','HETATM')):
            line=line[:22]+f'{int(line[22:26])+40:4d}'+line[26:]
        if line.startswith('TER'):
            continue
        lines.append(line)
    atom_only.write_text('\n'.join(lines)+'\n',encoding='ascii')
    source_manifest.append({'file':str(atom_only.relative_to(out)),'sha256':digest(atom_only),'source':str(original['04']),'synthetic':True,'exact_construction':'Remove SEQRES/REMARK/TER; increase author residue numbers by40; keep all atoms/identities'})
    no_marker=inputs/'10_without_explicit_index_marker.pdb'
    no_marker.write_text('\n'.join(line for line in original['10'].read_text(encoding='ascii').splitlines()
                                   if not line.startswith('REMARK'))+'\n',encoding='ascii')
    source_manifest.append({'file':str(no_marker.relative_to(out)),'sha256':digest(no_marker),'source':str(original['10']),'synthetic':True,
                            'exact_construction':'Remove REMARK provenance/index convention only; retain SEQRES and all coordinates. Q1/Q3 mapping must then stay uncertain.'})
    native_22=archive(Path(__file__).resolve().parents[1]/'tests/data/6wgb.pdb','22_native_archive.pdb')
    atezo=archive(next((bundle('36')/'sources/upstream').glob('*.cif')),'36_atezolizumab.cif')

    no_deletion=['reference_terminal_deletion','reference_internal_deletion']
    cases=[
        ('01_complete_FASTA',full,{}, {'domains':1,'coordinate_unassessed':True,'absent_codes':no_deletion+['possible_c_terminal_truncation']},'Complete curated RF128 sequence; FASTA has no coordinate information.'),
        ('02_RF127_no_reference',short,{}, {'domains':1,'coordinate_unassessed':True,'absent_codes':no_deletion+['possible_c_terminal_truncation']},'Known constructed minus1S, but absent reference a10-aa FR4 cannot prove truncation.'),
        ('03_RF127_with_reference',short,{'reference_path':full},{'domains':1,'reference_deleted_residues':1,'reference_deleted_regions':['FR4']},'Synthetic exact terminal1S deletion; expected relative difference is1.'),
        ('04_Nminus5_with_reference',cut_n,{'reference_path':full},{'domains':1,'reference_deleted_residues':5,'reference_deleted_regions':['FR1']},'Synthetic exact N-terminal5aa deletion.'),
        ('05_entire_FR4_removed',cut_fr4,{'reference_path':full},{'domains':1,'reference_deleted_residues':11,'reference_deleted_regions':['FR4'],'required_codes':['incomplete_domain_numbering']},'Synthetic exact11aa FR4 deletion; existing unnumberedCDR3 tail must not inflate missing-reference count.'),
        ('06_partial60',fragment,{}, {'domains':1,'required_codes':['incomplete_domain_numbering','possible_c_terminal_truncation']},'Synthetic partial domain shorter than conversion threshold80; auditor must still report it.'),
        ('07_unknown_X',bad,{}, {'domains':0,'required_codes':['invalid_sequence_symbols']},'Known invalid symbol; no completeness pass.'),
        ('08_mixed_valid_invalid_FASTA',both,{}, {'domains':1,'required_codes':['invalid_sequence_symbols']},'One invalid record must not block the valid record.'),
        ('09_5M2M_PDB_N_gap',original['10'],{}, {'domains':1,'FR1_missing_CA':2,'all_rows_missing_CA':2,'required_codes':['missing_residue_coordinates']},'Original deposited/local sequence contains Q1/V2; experimental coordinates lack them.'),
        ('10_5M2M_CIF_N_gap',cif_10,{'selected_chains':['D']}, {'domains':1,'FR1_missing_CA':2,'all_rows_missing_CA':2,'required_codes':['missing_residue_coordinates']},'Same selected5M2M source chain D in mmCIF; declared sequence/atom evidence.'),
        ('11_5M2M_PDB_plus_FASTA',original['10'],{'sequence_path':full_10}, {'domains':1,'FR1_missing_CA':2,'all_rows_missing_CA':2},'Independent completeFASTA paired with original structure; missing coordinates remain2.'),
        ('12_cAbLys3_zero_occupancy',original['13'],{}, {'domains':1,'all_rows_missing_CA':6,'zero_occupancy_rows':5,'required_codes':['zero_occupancy_coordinates']},'1MEL source: N-terminalD absent; five original residues73–77 have zero occupancy.'),
        ('13_Dupilumab_PDB',original['22'],{}, {'domains':2,'all_rows_missing_CA':13,'FR_missing_CA':1},'Prior independent coordinate audit: H62,H101–110,L32–33; only H62 is IMGT FR.'),
        ('14_Dupilumab_CIF',cif_22,{'selected_chains':['A','B']}, {'domains':2,'all_rows_missing_CA':13,'FR_missing_CA':1},'Same source declaration and selectedA/B copies in mmCIF; constant domains do not enter FR.'),
        ('15_ALB8_noncanonical',original['04'],{}, {'domains':1,'all_rows_missing_CA':0,'minimum_noncanonical':2,'absent_codes':no_deletion},'Archived source supports IMGT118/119 S/S; retain sequence.'),
        ('16_antiTNF_noncanonical',original['15'],{}, {'domains':1,'all_rows_missing_CA':0,'minimum_noncanonical':1,'absent_codes':no_deletion},'Archived source supports IMGT118 R; retain sequence.'),
        ('17_light_ends127',fab_20,{}, {'domains':2,'coordinate_unassessed':True,'absent_codes':no_deletion+['possible_c_terminal_truncation']},'CompleteBelimumabVL has10-aa FR4 endingIMGT127; Thera sequence verified previously.'),
        ('18_atom_only_author_starts41',atom_only,{}, {'domains':1,'required_codes':['coordinate_only_sequence'],'absent_codes':no_deletion+['missing_residue_coordinates']},'Author numbering offset is not evidence of40 missing sequence residues.'),
        ('19_RF127_PDB_plus128_FASTA',original['08'],{'sequence_path':full}, {'domains':1,'FR_missing_CA':1,'all_rows_missing_CA':1,'required_codes':['sequence_structure_disagreement','missing_residue_coordinates']},'Full intended128FASTA paired with127author template; terminalS lacks coordinates relativeFASTA.'),
        ('20_PDB_Q1_Q3_mapping_ambiguous',no_marker,{}, {'domains':1,'FR1_missing_CA':1,'all_rows_missing_CA':1,'uncertain_mapping_rows':2,'required_codes':['ambiguous_coordinate_mapping']},'Without explicit sequence-index convention, Q1/Q3 residue assignment is ambiguous; onlyV2 is definitely absent, whileCIF provides exactlabels.'),
        ('21_native_PDB_repeats_ambiguous',native_22,{'selected_chains':['A','B']}, {'domains':2,'all_rows_missing_CA':11,'FR_missing_CA':1,'uncertain_mapping_rows':3,'required_codes':['ambiguous_coordinate_mapping']},'NativePDB lacks label_seq_id; R100/R106/R108 assignments remain ambiguous. Definite missing positions11 is a lower bound;CIF establishes13.'),
        ('22_CIF_auth_chain_namespace',atezo,{'selected_chains':['F','K'],'namespace':'auth'}, {'domains':2,'all_rows_missing_CA':0,'source_auth_chains':['F','K'],'source_label_chains':['C','B']},'5X8L authF/K correspond to selectedVH/VL; both domains coordinate-complete.'),
        ('23_CIF_label_chain_namespace',atezo,{'selected_chains':['C','B'],'namespace':'label'}, {'domains':2,'all_rows_missing_CA':0,'source_auth_chains':['F','K'],'source_label_chains':['C','B']},'5X8L labelC/B must select the same domains as authF/K.'),
    ]
    dump(out/'input_provenance.json',source_manifest)
    numberer=Numberer(engine='anarcii',ncpu=args.ncpu,scheme='imgt')
    begin=perf_counter();results=[]
    for name,path,kwargs,expected,basis in cases:
        before={str(p):digest(p) for p in [path,*[p for p in kwargs.values() if isinstance(p,Path)]]}
        started=perf_counter()
        try:
            report=audit_file(path,numberer=numberer,**kwargs)
            write_audit(report,out/'cases'/name)
            checks=check(report,expected)
            unchanged=all(digest(Path(p))==h for p,h in before.items())
            checks.append({'check':'inputs_unchanged','actual':unchanged,'expected':True,'passed':unchanged})
            row={'case':name,'passed':all(c['passed'] for c in checks),'seconds':perf_counter()-started,
                 'expectation_basis':basis,'checks':checks,'risk_counts':report['risk_counts']}
        except Exception as exc:
            row={'case':name,'passed':False,'seconds':perf_counter()-started,'expectation_basis':basis,'error':repr(exc)}
        results.append(row)
        print(json.dumps({'case':name,'passed':row['passed'],'seconds':round(row['seconds'],3)},ensure_ascii=True),flush=True)
        dump(out/'progress.json',results)
    summary={'cases':len(results),'passed':sum(r['passed'] for r in results),'failed':sum(not r['passed'] for r in results),
             'total_seconds':perf_counter()-begin,'ncpu':args.ncpu,'cpu_only':True,
             'timing_scope':'One process reuses one model/cache. First case includes model load; later cases may reuse sequence numbering.',
             'scope':'Specified coordinate/source/synthetic-construction facts and conservative-reporting regressions; not an independent biological boundary gold standard.',
             'results':results}
    dump(out/'summary.json',summary)
    lines=['# FR审计真实回归','',f"共{summary['cases']}例，通过{summary['passed']}，失败{summary['failed']}。总耗时{summary['total_seconds']:.2f}秒。",'',
           '同进程复用模型；首例包含模型加载，后续可能使用序列缓存，不能把该时间当成每次独立启动耗时。','',
           '预期依据为原结构声明/坐标、此前独立核对以及明确构造的删改。合成样本不冒充实验构建体；本测试也不构成全部FR边界的独立生物学金标准。','',
           '| 用例 | 通过 | 秒 | 预期依据 |','|---|---|---:|---|']
    lines += [f"| {r['case']} | {r['passed']} | {r['seconds']:.3f} | {r['expectation_basis']} |" for r in results]
    (out/'README_回归结果.md').write_text('\n'.join(lines)+'\n',encoding='utf-8-sig')
    print(json.dumps({k:summary[k] for k in ('cases','passed','failed','total_seconds')},ensure_ascii=True))
    if summary['failed']:
        raise SystemExit(1)


if __name__=='__main__':
    main()

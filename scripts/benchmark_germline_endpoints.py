"""Frozen independent endpoint challenge. Run prepare, number, then score (see --help)."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
from hashlib import sha256
from importlib import metadata, resources
import json
from pathlib import Path

from Bio.Align import PairwiseAligner
from scaffold.numbering import AA, Numberer, classify
from scaffold.terminal_references import compare_terminals


SOURCE = 'https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/static/downloads/TheraSAbDab_SeqStruc_OnlineDownload.csv'
SOURCE_SHA = 'cc9402788297ed7b6ab3d676daa21953a6d28ebbd9cc78ec7248b19f74d49d6d'
INN_ALCESTOBART = 'https://cdn.who.int/media/docs/default-source/international-nonproprietary-names-%28inn%29/rl93.pdf'
GENETICS = 'Genetics (Bispecifics delimited with semicolon)'
GOALS = {'VHH': 8, 'VH': 10, 'K': 10, 'L': 10, 'out_of_panel': 8}
JUNCTION_GOALS = {'VHH': 3, 'VH': 3, 'K': 3, 'L': 3, 'out_of_panel': 2}
LINKER_GOALS = {'VHH': 2, 'VH': 2, 'K': 2, 'L': 2, 'out_of_panel': 2}
# These are synthetic V+C controls using published human CH1, kappa and lambda
# constant N-termini. They do not assert the actual isotype of that therapeutic drug.
CONSTANT = {'H': 'ASTKGPSVFPLAPSSKSTSGGTAALG', 'K': 'RTVAAPSVFIFPPSDEQLKSGTAS',
            'L': 'GQPKAAPSVTLFPPSDEQLKSGTAS'}


def checksum(text):
    return sha256(text.encode('ascii')).hexdigest()


def save(path, data):
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')


def load(path):
    return json.loads(path.read_text(encoding='utf-8'))


def eligible(seq):
    return 85 <= len(seq) <= 155 and set(seq) <= AA


def candidates(row):
    fmt, name, genetics, notes = row['Format'], row['Therapeutic'], row[GENETICS], row['Notes']
    heavy, light = row['HeavySequence'], row['LightSequence']
    if not eligible(heavy):
        return
    item = dict(name=name, genetics=genetics, format=fmt, notes=notes)
    if 'Nanobody' in fmt or (('Single Domain' in fmt or 'sc-5-VH' in fmt)
            and any(k in notes.lower() for k in ['llama','vicpac','alpaca','camelid'])):
        yield {**item, 'group':'VHH', 'role':'H', 'seq':heavy, 'partner':''}
    if fmt not in ('Whole mAb', 'Fab') or not eligible(light) or row['VD LC'] not in ('Kappa', 'Lambda'):
        return
    if genetics not in ('Humanised', 'Genetically human', 'Murine'):
        return
    light_type = {'Kappa':'K', 'Lambda':'L'}[row['VD LC']]
    for role, seq, partner in [('H',heavy,light),(light_type,light,heavy)]:
        group = ('out_of_panel' if genetics == 'Murine' else ('VH' if role == 'H' else role))
        yield {**item, 'group':group, 'role':role, 'seq':seq, 'partner':partner}


def select(rows, prior, quotas):
    refs = json.loads(resources.files('scaffold').joinpath('static/terminal_references.json').read_text(encoding='utf-8'))['entries']
    align = PairwiseAligner(mode='global', match_score=0, mismatch_score=-1, open_gap_score=-1, extend_gap_score=-1)
    def too_close(s, t):
        return abs(len(s)-len(t)) <= max(len(s),len(t))*.10 and -align.score(s,t) <= max(len(s),len(t))*.10
    pool = [item for row in rows for item in candidates(row)]
    pool.sort(key=lambda d: checksum(d['seq'] + d['name']))
    selected, excluded = [], Counter()
    for group, quota in quotas.items():
        accepted = 0
        for item in (i for i in pool if i['group'] == group):
            if any(i['name'].rstrip('1234567890') == item['name'].rstrip('1234567890') for i in [*prior, *selected]):
                excluded['same_therapeutic'] += 1
                continue
            if any(too_close(item['seq'],r['sequence']) for r in refs):
                excluded['close_to_catalog'] += 1
                continue
            if any(too_close(item['seq'],i['seq']) for i in [*prior, *selected]):
                excluded['close_to_selected_representative'] += 1
                continue
            selected.append(item)
            accepted += 1
            if accepted == quota:
                break
    return selected, excluded


def perturb(item, junction=False, linker_only_controls=False):
    seq = item['seq']
    yield 'clean', seq, 'none', 0
    for k in (1,2,3):
        yield f'N_minus_{k}', seq[k:], 'N', k
        yield f'C_minus_{k}', seq[:-k], 'C', k
    yield 'start_M_added', 'M'+seq, 'none', 0
    yield 'start_M_replaces_first', 'M'+seq[1:], 'none', 0
    yield 'C_last_substitution', seq[:-1]+('A' if seq[-1]!='A' else 'G'), 'none', 0
    yield 'intact_His6', seq+'HHHHHH', 'none', 0
    yield 'C_minus_1_His6', seq[:-1]+'HHHHHH', 'C', 1
    if item['partner']:
        constant = CONSTANT[item['role']]
        yield 'intact_fused_constant', seq+constant, 'none', 0
        yield 'C_minus_1_fused_constant', seq[:-1]+constant, 'C', 1
        if junction:
            yield 'intact_last_equals_constant_start',seq[:-1]+constant[0]+constant,'none',0
    if junction:
        linker = 'GGGGS'*3
        yield 'intact_last_H_His6',seq[:-1]+'H'+'HHHHHH','none',0
        # Identical uploaded FASTA to C_minus_1_His6 with a different intended boundary.
        yield 'intact_last_H_His5_ambiguous',seq[:-1]+'H'+'HHHHH','none',0
        yield 'intact_GS_linker', seq+linker,'none',0
        yield 'C_minus_1_GS_linker', seq[:-1]+linker,'C',1
        yield 'intact_last_G_GS_linker',seq[:-1]+'G'+linker,'none',0
        if linker_only_controls:
            yield 'intact_last_G_shortened_linker_ambiguous',seq[:-1]+'G'+linker[1:],'none',0
            yield 'intact_GS_linker_2x',seq+'GGGGS'*2,'none',0
            yield 'C_minus_1_GS_linker_2x',seq[:-1]+'GGGGS'*2,'C',1
            yield 'C_minus_2_His6',seq[:-2]+'HHHHHH','C',2
            yield 'intact_EPEA',seq+'EPEA','none',0
            yield 'C_minus_1_EPEA',seq[:-1]+'EPEA','C',1
            yield 'intact_initial_A_substitution','A'+seq[1:],'none',0
            yield 'N_minus_1_add_A','A'+seq[1:],'N',1


def prepare(csv_path, out, prior_panels=None):
    if (out/'panel.json').exists():
        raise ValueError('Panel is already frozen; use its existing number/score stages')
    if sha256(csv_path.read_bytes()).hexdigest() != SOURCE_SHA:
        raise ValueError('Expected the public PL132 snapshot specified in this script')
    prior_panels=prior_panels or []
    prior=[parent for path in prior_panels for parent in load(path)['parents']]
    quotas=(LINKER_GOALS if len(prior_panels)>1 else JUNCTION_GOALS) if prior else GOALS
    with csv_path.open(encoding='utf-8-sig', newline='') as handle:
        selected, exclusions = select(csv.DictReader(handle),prior,quotas)
    out.mkdir(exist_ok=True,parents=True)
    save(out/'panel.json', dict(source=SOURCE, csv_sha256=SOURCE_SHA, snapshot='PL132',
         selection='One deterministic representative per therapeutic and 90% edit-similarity neighbourhood. No parent selected on numbering, endpoint result or germline prediction. No model is fitted.',
         excluded_while_filling=dict(exclusions), requested=quotas, parents=selected,
         **({'prior_panel_sha256':sha256(prior_panels[0].read_bytes()).hexdigest()} if len(prior_panels)==1 else
            {'prior_panels_sha256':[sha256(p.read_bytes()).hexdigest() for p in prior_panels]} if prior else {})))
    cases=[]
    for item in selected:
        for treatment, query, side, k in perturb(item, junction=bool(prior), linker_only_controls=len(prior_panels)>1):
            cases.append(dict(id=f"{item['name']}-{item['role']}-{treatment}", name=item['name'], group=item['group'],
                              role=item['role'], treatment=treatment, query=query, partner=item['partner'], expected_side=side, removed=k))
    save(out/'cases.json',cases)
    print('Frozen:',dict(Counter(i['group'] for i in selected)),'cases:',len(cases),flush=True)


def number(out, ncpu):
    from anarcii import Anarcii
    cases=load(out/'cases.json')
    seqs=sorted({s for c in cases for s in (c['query'],c['partner']) if s})
    model=Anarcii(seq_type='antibody', mode='accuracy', cpu=True, ncpu=ncpu, verbose=False)
    numberer=Numberer(ncpu=ncpu)
    numberer._model=model
    for i in range(0,len(seqs),64):
        batch=seqs[i:i+64]
        results=model.number([(str(j),s) for j,s in enumerate(batch)])
        for j,s in enumerate(batch):
            result=results[str(j)]
            numberer._cache[s] = ([(result['numbering'],result['chain_type'],result['score'],result['query_start'])]
                                  if result.get('numbering') and not result.get('error') else [])
        print('Fresh numbering:', min(i+64,len(seqs)), '/', len(seqs),flush=True)
    # Run production domain recovery, chain typing and selection on actual variants.
    # Cache one compact record per input to avoid duplicating partner chains in the evidence.
    with gzip.open(out/'numbered.jsonl.gz','wt',encoding='utf-8') as handle:
        for c in cases:
            records={'test':c['query']}
            if c['partner']:
                records['partner']=c['partner']
            try:
                annotation=classify(records,numberer,'fab' if c['partner'] else 'vhh')
                domain,=[d for d in annotation['domains'] if d.record_id=='test']
                if domain.chain_type!=c['role']:
                    raise ValueError(f'Expected {c["role"]}, got {domain.chain_type}')
                chain='H' if domain.chain_type=='H' else 'L'
                meta=dict(chains=[dict(output_chain_id=chain,chain_type=domain.chain_type,domain_id=domain.domain_id,
                              record_id='test',full_variable_sequence=domain.sequence,full_source_sequence=c['query'])],
                          full_input_sequences=records, domains=[domain.to_dict()], numbering_scheme='imgt', cdr_definition='imgt')
                rows=[{**r,'output_chain_id':chain} for r in domain.rows]
                record=dict(id=c['id'],meta=meta,rows=rows)
            except (ValueError,RuntimeError) as exc:
                record=dict(id=c['id'],error=str(exc))
            handle.write(json.dumps(record,separators=(',',':'))+'\n')
    save(out/'numbering_provenance.json',dict(engine='anarcii',version=metadata.version('anarcii'),cpu=True,mode='accuracy',
         inputs_sha256=sha256((out/'cases.json').read_bytes()).hexdigest(), sequences_numbered=len(seqs)))


def preliminary(meta, rows, versions):
    """Obtain baseline length/anchor and numbering evidence from the production rules."""
    from scaffold.postcheck_sequence import _terminal_checks, _check_framework
    warnings = dict(N=[], C=[])
    chain = meta['chains'][0]['output_chain_id']
    def collect(region, code, status, message, **evidence):
        if code in {'possible_n_terminal_truncation','possible_c_terminal_truncation',
                    'atypical_imgt_length','noncanonical_anchor','terminal_motif_review'}:
            warnings[{'FR1':'N','FR4':'C'}[region]].append(code)
    _terminal_checks(meta['domains'][0], meta['full_input_sequences']['test'], collect)
    for side, region in [('N','FR1'), ('C','FR4')]:
        segment = [r for r in rows if r['region'] == region]
        _check_framework(chain, region, segment, None, [], meta, collect)
    for f in versions['findings']:
        if f['code']=='terminal_motif_review':
            collect(f['region'],f['code'],f['status'],f['message'])
    return warnings


def measure(group):
    return dict(total=len(group), failed=sum(r['state']=='extraction_failed' for r in group),
                any_review=sum(r['alert'] for r in group), compatible=sum(r['state']=='compatible' for r in group),
                suspected_short=sum(r['state']=='suspected_short' for r in group),
                boundary_review=sum(r['state']=='boundary_review' for r in group),
                abstain=sum(r['state'] in ('indeterminate','extraction_failed') for r in group))


def score(out, fixtures=None, regression=False):
    from scaffold.germline_endpoints import assess_endpoints
    from scaffold.postcheck import _fasta_endpoints
    fixtures=fixtures or out
    out.mkdir(parents=True,exist_ok=True)
    if (out/'predictions.json').exists():
        raise ValueError('Scored output already exists. Use --out NEW_FOLDER --fixtures FROZEN_FIXTURES to rescore.')
    cases = {c['id']:c for c in load(fixtures/'cases.json')}
    if load(fixtures/'numbering_provenance.json')['inputs_sha256'] != sha256((fixtures/'cases.json').read_bytes()).hexdigest():
        raise ValueError('Numbering cache does not match frozen cases')
    results=[]
    stages=['existing_rules','plus_context','germline_deployed']
    with gzip.open(fixtures/'numbered.jsonl.gz','rt',encoding='utf-8') as handle:
        for line in handle:
            example=json.loads(line)
            c=cases[example['id']]
            if 'error' not in example:
                meta,rows=example['meta'],example['rows']
                versions=compare_terminals(meta)
                warnings=preliminary(meta,rows,versions)
                context={a['side']:a for a in assess_endpoints(meta,rows,use_germline=False)['assessments']}
                current_findings=[dict(chain=meta['chains'][0]['output_chain_id'],region=region,code=code)
                                  for side,region in [('N','FR1'),('C','FR4')] for code in warnings[side]]
                current=_fasta_endpoints(meta, rows, meta['full_input_sequences'], versions, current_findings)
                germline={a['side']:a for a in current['assessments']}
            for side,region in [('N','FR1'),('C','FR4')]:
                for stage in stages:
                    reason=None
                    if 'error' in example:
                        state='extraction_failed'
                        reason=example['error']
                        codes=[]
                    elif stage=='germline_deployed':
                        state=germline[side]['state']
                        reason=germline[side]['evidence']['reason']
                        codes=[f['code'] for f in current_findings if f['region']==region]
                    else:
                        codes=list(warnings[side])
                        state='legacy_review' if codes else 'indeterminate'
                        if stage=='plus_context' and context[side]['state']=='boundary_review':
                            state='boundary_review'
                            codes.append('context_boundary_review')
                    results.append(dict(id=c['id'],name=c['name'],group=c['group'],role=c['role'],
                        treatment=c['treatment'],side=side,truth_missing=c['expected_side']==side,
                        removed=c['removed'] if c['expected_side']==side else 0, stage=stage, state=state,
                        alert=bool(codes),codes=codes,reason=reason,
                        fr4_length=None if 'error' in example else sum(r['region']=='FR4' for r in rows),
                        input_ends_VTVS=c['query'].endswith('VTVS')))
    save(out/'predictions.json',results)
    cells=[]
    for stage in stages:
        for group in ['all', *GOALS]:
            for side in ('N','C'):
                selected=[r for r in results if r['stage']==stage and r['side']==side and (group=='all' or r['group']==group)]
                for treatment in sorted({r['treatment'] for r in selected}):
                    sample=[r for r in selected if r['treatment']==treatment]
                    cells.append(dict(stage=stage,group=group,side=side,treatment=treatment,**measure(sample)))
    panel=load(fixtures/'panel.json')
    path=out/'metrics.json'
    save(path,dict(parents=len(panel['parents']),cases=len(cases),parent_group_counts=dict(Counter(p['group'] for p in panel['parents'])),
        unscored_cases=sum(r['state']=='extraction_failed' and r['stage']==stages[-1] and r['side']=='C' for r in results),
        unscored_by_reason=dict(Counter(r['reason'] for r in results if r['state']=='extraction_failed' and r['stage']==stages[-1] and r['side']=='C')),
        catalogue_version=json.loads(resources.files('scaffold').joinpath('static/germline_endpoints.json').read_text(encoding='utf-8'))['version'],
        definitions=dict(any_review='Any side-specific anomaly, including noncanonical anchors or an atypical FR length. Not a claim of physical deletion.',
        suspected_short='Germline endpoint-specific, references predict a longer end. Boundary review and noncanonical anchor do not count.',
        abstain='No endpoint decision or extraction failed. The baseline cannot certify intactness from absent alerts.'),
        system_sha256={f:sha256(resources.files('scaffold').joinpath(f).read_bytes()).hexdigest()
                       for f in ('germline_endpoints.py','postcheck.py','postcheck_sequence.py','terminal_references.py',
                                 'static/germline_endpoints.json','static/terminal_references.json')},
        inputs_dir=str(fixtures.resolve()), panel_sha256=sha256((fixtures/'panel.json').read_bytes()).hexdigest(), rows=cells))
    report(out, panel, cases, results, regression)
    print('Scored',len(cases),'queries;',path,flush=True)


def fraction(value, n):
    return f'{value}/{n} ({100*value/n:.1f}%)' if n else '无样本'


def report(out, panel, cases, results, regression):
    stats=load(out/'metrics.json')
    source_conflicts = stats['unscored_by_reason']
    tested=[r for r in results if r['stage']=='germline_deployed']
    is_new_panel=bool(panel.get('prior_panel_sha256') or panel.get('prior_panels_sha256'))
    title=f"原 {len(cases)} 项回归" if regression else (
          '独立连接肽留出测试' if panel.get('prior_panels_sha256') else
          '独立衔接区留出测试' if is_new_panel else '库外可变域端点盲裁测试')
    lines=['# '+title, '',
       f"输入来自固定版本的 Thera-SAbDab PL132。{len(panel['parents'])} 个母本先去重再构造 {len(cases)} 个测试输入。"
       f"分组为 {stats['parent_group_counts']}。"+('本集合用于重查已发现的问题。新的留出集单独报告。' if regression
           else '检测规则冻结后直接评分，未用本集合拟合阈值。'),
       '', '每个变体使用实跑 ANARCII（IMGT、CPU、accuracy），随后执行生产代码的域识别和链型检查。'
       '对 Fab，伴侣链保留原序列；评估单位为被扰动链的一个端点。每个母本的多个裁短样例相关。'
       f"未计分 {stats['unscored_cases']} 份，保留在全量分母及弃判中。测试对象是编号和后检，不计 BoltzGen 运行。",
       '', (f"Alcestobart 的 14 个变体均得到 κ，Thera-SAbDab PL132 标成 λ，源端点结尾为 TVL；WHO 列出 κ、V 区至第 108 位，工程化末端 TVLG。因此其 14 个输入没有合格的参考答案。保留原清单和失败原因，未把该冲突计为算法分类错误。[WHO 推荐国际非专利药名清单 93]({INN_ALCESTOBART})"
            if 'Alcestobart-L-clean' in cases else f"无法评分的原因及数量为 {source_conflicts or '无'}。"),
       '', '## 端点提示的消融比较', '',
       '表格前 3 列计入该端任何需复核事项，包含长度及锚点异常；末两列分别统计明确的胚系“疑似短端”及无法判定。'
       '完整序列行对应误报比例，裁短行对应同端检出比例。分子/分母均为域的数目。', '',
       '| 端点及处理 | 既有规则提示 | 加上下文提示 | 加胚系后的实际提示 | 其中疑似短端 | 胚系弃判 |',
       '|---|---:|---:|---:|---:|---:|']
    def cell(group,stage,side,condition):
        return next(c for c in stats['rows'] if c['group']==group and c['stage']==stage and c['side']==side and c['treatment']==condition)
    for side in ('N','C'):
        for condition in ['clean', *(f'{side}_minus_{k}' for k in (1,2,3))]:
            a,b,c=[cell('all',stage,side,condition) for stage in ('existing_rules','plus_context','germline_deployed')]
            lines.append(f'| {side} {condition} | '+' | '.join(fraction(x['any_review'],x['total']) for x in [a,b,c])+
                         f" | {fraction(c['suspected_short'],c['total'])} | {fraction(c['abstain'],c['total'])} |")
    lines.extend(['','## 对照与局限','',
         '| 已知处理 | 端点 | 样本 | 任何提示 | 疑似短端 | 边界复核 | 与参考长度兼容 | 弃判 |','|---|---|---:|---:|---:|---:|---:|---:|'])
    challenge_names=['intact_last_H_His6','intact_last_H_His5_ambiguous','intact_GS_linker','C_minus_1_GS_linker',
                     'intact_last_G_GS_linker','intact_last_equals_constant_start',
                     'intact_last_G_shortened_linker_ambiguous','intact_GS_linker_2x','C_minus_1_GS_linker_2x',
                     'C_minus_2_His6','intact_EPEA','C_minus_1_EPEA']
    extra=[(name,'C') for name in challenge_names if any(c['treatment']==name for c in cases.values())]
    extra += [(name,'N') for name in ('intact_initial_A_substitution','N_minus_1_add_A')
              if any(c['treatment']==name for c in cases.values())]
    for condition,side in [('start_M_added','N'),('start_M_replaces_first','N'),('C_last_substitution','C'),
                           ('intact_His6','C'),('C_minus_1_His6','C'),('intact_fused_constant','C'),('C_minus_1_fused_constant','C'),*extra]:
        c=cell('all','germline_deployed',side,condition)
        lines.append(f"| {condition} | {side} | {c['total']} | {c['any_review']} | {c['suspected_short']} | {c['boundary_review']} | {c['compatible']} | {c['abstain']} |")
    if extra:
        lines+=['','每个母本都有一对完全相同的上传字符串。其可解释为域缺一位再接 His6，也可解释为有意把域末残基换为 H 再接 His5。该配对只检验判定的一致性。']
    lines.extend(['','添加起始 M、替换第一个氨基酸和替换最后一个氨基酸的样本分别列出。这些原始字符串可能有多个生物学解释，均不纳入完整母本误报率。',
            '', '带标签及恒定区的裁短序列可能将下游首个残基编号进 V 域。与胚系尾长兼容仍可能漏报，尤其是端点带有允许的残基替换时。',
            '', '## 分组和已发表的短尾构建体','',
            '| 样本 | 端点与处理 | 数量 | 任何提示 | 疑似短端 | 边界复核 | 弃判 |',
            '|---|---|---:|---:|---:|---:|---:|'])
    def show(group, title, sample):
        v=measure(sample)
        lines.append(f"| {group} | {title} | {v['total']} | {v['any_review']} | {v['suspected_short']} | {v['boundary_review']} | {v['abstain']} |")
    deployed=[r for r in results if r['stage']=='germline_deployed']
    for group in GOALS:
        for side in ('N','C'):
            for condition in ('clean', f'{side}_minus_1'):
                show(group, side+' '+condition, [r for r in deployed if r['group']==group and r['side']==side and r['treatment']==condition])
    short_clean=[r for r in deployed if r['side']=='C' and r['treatment']=='clean' and r['fr4_length'] is not None and r['fr4_length']<=10]
    show('原始 FR4 ≤10 aa', 'C clean', short_clean)
    show('H 裁后以 VTVS 结尾','C C_minus_1',[r for r in deployed if r['side']=='C' and r['role']=='H' and r['treatment']=='C_minus_1' and r['input_ends_VTVS']])
    for stage in ('existing_rules','plus_context','germline_deployed'):
        selected=[r for r in results if r['stage']==stage]
        show(stage,'两端，全部输入',selected)
    lines+=['','## 有可用参考答案时的条件指标','',
            '| 端点 | 处理 | 全量 | 可计分 | 疑似短端/可计分 | 边界复核/可计分 | 弃判/可计分 |',
            '|---|---|---:|---:|---:|---:|---:|']
    for side in ('N','C'):
        for treatment in ('clean',*(f'{side}_minus_{k}' for k in (1,2,3))):
            x=[r for r in tested if r['side']==side and r['treatment']==treatment]
            valid=[r for r in x if r['state']!='extraction_failed']
            m=measure(valid);n=m['total']
            lines.append(f"| {side} | {treatment} | {len(x)} | {n} | {fraction(m['suspected_short'],n)} | {fraction(m['boundary_review'],n)} | {fraction(m['abstain'],n)} |")
    lines.extend(['','已发表的原始短尾按照完整母本重新编号后 FR4 ≤10 aa 分组。此分组描述可变域的端点长度，治疗构建体包含天然及工程序列。潜在阳性单独列出所有重链从双 S 裁成 VTVS 的样本。',
            '', '相似性筛选采用全长编辑距离，排除与项目中 17 个已知版本或任一入选母本达到 90% 相似度的候选。每种药物最多贡献一个母本；未进行阈值调参。'
            '这种分组防止近重复序列泄漏，不能据此断言各组的 V/J 基因家族互不重叠。',
            '', '参考端点来自外部整理的命名治疗构建体。鼠源完整域是人/羊驼覆盖范围外的阴性。标签和恒定区为人工拼接，真实应用中其起点可不同。'
            '同一原始 FASTA 可能来自有意或无意的截短，因此本验证只针对有记录的端点和受控差异。'
            '数据集富集治疗抗体，覆盖外、短尾和各链型实际比例均不代表客户输入。请勿由此估计商业部署精度或阳性预测值。',
            '', '文件：`panel.json` 保留母本、分组、来源，`cases.json` 保留独立构造答案，'
            '`predictions.json` 保留每端逐例判定，`metrics.json` 含各链型、处理、端点和阶段的分母。',
            '', f'[Thera-SAbDab 数据源]({SOURCE})。下载文件 SHA256 为 `{SOURCE_SHA}`。数据库后续版本可能改变；冻结的母本和变体均保存在本目录。',
            '', '在项目仓库已激活的 scaffold 环境运行下列命令。已有冻结清单可从 number 开始，已有带编号的压缩结果可从 score 开始。','',
            '```powershell',
            'python scripts/benchmark_germline_endpoints.py prepare --csv TheraSAbDab_SeqStruc_OnlineDownload.csv --out benchmark',
            'python scripts/benchmark_germline_endpoints.py number --ncpu 4 --out benchmark',
            'python scripts/benchmark_germline_endpoints.py score --out benchmark','```'])
    (out/'README_盲测结果.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['prepare','number','score'])
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--fixtures', type=Path, help='Score immutable artifacts from another directory')
    parser.add_argument('--regression', action='store_true', help='Label an evaluation of previously seen cases as a regression')
    parser.add_argument('--after-panel', action='append', type=Path, help='Exclude all parents of this previous panel (repeatable)')
    parser.add_argument('--csv', type=Path, help='Frozen Thera-SAbDab PL132 table for prepare')
    parser.add_argument('--ncpu', type=int, default=2)
    args=parser.parse_args()
    if args.stage=='prepare':
        prepare(args.csv,args.out,args.after_panel)
    elif args.stage=='number':
        number(args.out,args.ncpu)
    else:
        score(args.out,args.fixtures,args.regression)


if __name__ == '__main__':
    main()

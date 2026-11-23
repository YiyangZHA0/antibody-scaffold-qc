"""Inspect the artifacts of one completed extraction using its saved region definition."""
from collections import Counter
from hashlib import sha256
from itertools import groupby
import json
from pathlib import Path
from uuid import uuid4

import yaml

from .errors import ScaffoldError
from .numbering import read_fasta
from .postcheck_sequence import postcheck_sequences
from .postcheck_structure import postcheck_coordinates
from .terminal_references import compare_terminals
from .validation import _sequence_tokens


def _finding(code, status, message, **evidence):
    return dict(chain=None, region=None, code=code, status=status, message=message, evidence=evidence)


def _verified_source(meta, relocated, findings):
    original = Path(relocated or meta.get('input_file') or '')
    expected = meta.get('input_sha256')
    if original.is_file() and expected and sha256(original.read_bytes()).hexdigest() == expected:
        return original
    if relocated:
        raise ScaffoldError(f'The specified source file is missing or differs from the recorded input digest: {original}')
    findings.append(_finding('original_unverified', 'not_assessable',
                             'Original upload unavailable or changed. Comparisons use saved conversion records.',
                             recorded_source=str(original)))
    return None


def _validation_evidence(directory, findings):
    validation = directory / 'validation.json'
    if not validation.is_file():
        return 'unavailable'
    saved = json.loads(validation.read_text(encoding='utf-8-sig'))
    digests = saved.get('artifact_sha256') or {}
    changed = []
    for name in ('scaffold.yaml', 'scaffold.cif', 'residue_map.jsonl', 'design.yaml', 'scaffold.yml'):
        if name in digests:
            path = directory / name
            if not path.is_file() or sha256(path.read_bytes()).hexdigest() != digests[name]:
                changed.append(name)
    if changed:
        findings.append(_finding('changed_since_validation', 'observed',
                                 'Files differ from the last saved BoltzGen validation.', files=changed))
        return 'changed'
    if saved.get('official_parser_passed') and digests.get('scaffold.yaml') and digests.get('residue_map.jsonl'):
        return 'recorded_parser_pass_hashes_match'
    return 'unverified'


def _check_sequence_payload(path, meta, rows):
    """Compare the fixed letters and loop tokens read by BoltzGen with the saved map."""
    if meta['scaffold_kind'] != 'sequence':
        return []
    findings = []
    try:
        entities = yaml.safe_load(path.read_text(encoding='utf-8-sig'))['entities']
        proteins = [entity['protein'] for entity in entities]
        by_id = {p['id']: _sequence_tokens(p['sequence']) for p in proteins}
        expected_ids = {chain['output_chain_id'] for chain in meta['chains']}
        if len(proteins) != len(by_id) or set(by_id) != expected_ids:
            raise ValueError('protein chain IDs differ from the saved extraction')
        for chain_id, actual in by_id.items():
            expected = []
            selected = (r for r in rows if r['output_chain_id'] == chain_id and r['included'])
            for (designed, region), run in groupby(selected, key=lambda r: (r['designed'], r['region'])):
                if designed:
                    expected.append(str(meta['loop_lengths'][f'{chain_id}-{region}']))
                else:
                    expected.extend(r['aa'] for r in run)
            if actual != expected:
                issue = _finding('sequence_yaml_payload_mismatch', 'observed',
                                 'Fixed letters or loop lengths in entities/protein differ from the residue map.',
                                 expected=''.join(expected), actual=''.join(actual))
                issue['chain'] = chain_id
                findings.append(issue)
    except (KeyError, TypeError, ValueError, yaml.YAMLError, ScaffoldError) as exc:
        findings.append(_finding('sequence_yaml_payload_invalid', 'observed',
                                 'Cannot compare the exported protein specification.', reason=str(exc)))
    return findings


def _fasta_endpoints(meta, rows, raw, versions, findings):
    """Original-record evidence is only available for a verified FASTA upload."""
    if raw is None:
        reason = ('original_fasta_unavailable' if Path(meta.get('input_file', '')).suffix.lower()
                  in {'.fa', '.faa', '.fasta'} else 'not_a_fasta_upload')
        return {'status': 'not_assessed', 'reason': reason, 'assessments': []}
    from .germline_endpoints import assess_endpoints
    domains = {d['domain_id']: d for d in meta.get('domains', [])}
    mapped_spans = {}
    for row in rows:
        chain_id = row.get('output_chain_id')
        position = row.get('input_seq_index')
        if chain_id not in mapped_spans:
            mapped_spans[chain_id] = [position, position]
        else:
            mapped_spans[chain_id][1] = position
    eligible, excluded = [], []
    for chain in meta['chains']:
        domain = domains.get(chain.get('domain_id'), {})
        record = raw.get(chain.get('record_id'))
        start, end = domain.get('start'), domain.get('end')
        if (record and isinstance(start, int) and isinstance(end, int) and 0 <= start < end <= len(record)
                and mapped_spans.get(chain['output_chain_id']) == [start + 1, end]
                and record[start:end] == chain.get('full_variable_sequence')
                and record == (meta.get('full_input_sequences') or {}).get(chain.get('record_id'))):
            eligible.append(chain)
        else:
            excluded.append(chain['output_chain_id'])
    result = assess_endpoints({**meta, 'chains': eligible}, rows, known_versions=versions) if eligible else {'assessments': [], 'findings': []}
    by_end = {(a['chain'], a['side']): a for a in result['assessments']}
    # A single ending motif or a numbering gap must not overrule a comparison
    # that explicitly found plausible alternative endpoints. Preserve the evidence.
    old_hints = {'possible_n_terminal_truncation', 'possible_c_terminal_truncation',
                 'short_unnumbered_terminal_flank', 'terminal_motif_review'}
    retained = []
    for item in findings:
        key = (item['chain'], {'FR1': 'N', 'FR4': 'C'}.get(item['region']))
        if item['code'] in old_hints and key in by_end:
            by_end[key].setdefault('preliminary_hints', []).append(item)
        else:
            retained.append(item)
    findings[:] = retained + result.pop('findings')
    for chain_id in excluded:
        for side in ('N', 'C'):
            result['assessments'].append(dict(chain=chain_id, side=side, state='indeterminate',
                                              evidence={'reason': 'source_domain_mapping_unverified'}))
    result['status'] = 'assessed' if not excluded else ('partially_assessed' if eligible else 'not_assessed')
    return result


_ENDPOINT_REASONS_ZH = {
    'source_domain_mapping_unverified': '无法核实原文件与可变域的对应关系',
    'source_mapping_unverified': '逐位映射无法与上传序列对齐',
    'exact_construct_version_available': '有全域序列精确匹配的已知版本，详见上表',
    'source_supports_longer_catalog_boundary': '已提取片段匹配短版本，紧邻的原始残基同时支持长版本的边界',
    'downstream_motif_overlaps_numbering': '标签、串联 G4S 接头或恒定区的候选起点与已编号的 FR4 重叠',
    'ambiguous_downstream_junction': '域末残基和下游片段两种归属均有依据',
    'same_length_terminal_residue_disagrees': '参考长度相容，但端点残基与所有同长候选不符，无法区分替换或边界偏移',
    'endpoint_length_only_residue_differs': '仅长度相容，末位氨基酸与所有同长候选不同，请核对构建版本',
    'matches_candidate_endpoint': '保留的候选参考与该端点长度相容',
    'insufficient_internal_framework_support': '可比较的内部框架长度或一致度不足',
    'missing_v_anchor': 'N 端定位缺少可靠的内部编号锚点',
    'missing_j_anchor': 'C 端缺少可靠的 J 区起始锚点',
    'missing_conserved_cys104': 'IMGT 104 位缺少典型 C，C 端的常规判断不适用',
    'numbered_initial_methionine': '原文件和编号域都从 M 开始，其归属需要构建记录',
    'ambiguous_numbered_start_methionine': '编号首位 M 与候选 V 起点不同，工程替换或起始 M 占据原首位均可解释',
    'start_register_uncertain': 'N 端超出本次短截检测范围',
    'terminal_register_unresolved': '参考与 N 端无法可靠对齐',
    'j_register_uncertain': 'J 区编号不连续或存在插入',
    'short_j_overlap': 'J 区重叠太短',
    'no_compatible_j': '没有相容的同链型 J 参考',
    'weak_j_match': 'J 区匹配不足',
    'j_endpoint_unresolved': '未找到可靠的相似 J 端点',
    'alternative_endpoints_or_longer_query': '候选支持不同终点，或输入长于候选',
    'candidate_residues_already_in_source': '候选延伸残基已在原序列紧邻域边界处出现',
    'unresolved_flank_origin': '原序列有域外残基，其来源无法确定',
    'initiator_or_domain_residue': '起始 M 与可变域起点的归属有歧义',
    'known_single_and_double_serine_endpoints': '已有可信单 S 和双 S 端点',
    'retained_candidates_extend_beyond_endpoint': '可比较的候选参考均延伸到该域端点之外',
    'outside_short_deletion_scope': '差异超过本方法评估的 1–3 位',
}


_DESCRIPTIONS_ZH = {
    'missing_residue_coordinates': '结构声明了这些残基，但在该模型中找不到对应坐标。核对结构版骨架实际纳入的范围。',
    'zero_occupancy_coordinates': '结构中的这些残基全部原子的占有率均不大于零，请复核原模型。',
    'missing_backbone_atoms': '这些残基缺少具有正占有率的主链原子，请核对原模型。',
    'coordinate_mapping_uncertain': '这些位置与坐标的对应关系有歧义。',
    'coordinate_only_sequence': '文件仅含坐标记录，未记录的末端残基无法计数。',
    'framework_header_mismatch': '导出的 FR 注释与逐残基映射有差异，请核对导出内容。',
    'framework_annotation_mismatch': 'FR 映射与提取时的编号记录不同，请核对序列、顺序和边界。',
    'sequence_yaml_payload_mismatch': 'YAML 实体中的固定氨基酸或 CDR 长度与提取记录不一致。',
    'changed_since_validation': '文件与上次 BoltzGen 校验时的内容不同。',
    'atypical_imgt_length': '该片段长度超出典型 IMGT 范围，请核实所用构建体。',
    'noncanonical_anchor': '保守位置出现非典型氨基酸，请查原始序列及工程化记录。',
    'short_unnumbered_terminal_flank': '原始序列与编号域之间留有一至两个残基，请复核提取端点。',
    'possible_n_terminal_truncation': '输入起点的 IMGT 编号晚于 1，请核对 N 端。',
    'possible_c_terminal_truncation': '原始序列到达末尾时仍未进入典型 IMGT FR4，请核对 C 端。',
    'terminal_motif_review': '末端 9 个残基与库中某些重链末位之前的片段相同。仅凭该短片段无法识别构建体，请核对预期终点。',
    'germline_endpoint_short': '提取域在此端比相容的胚系候选短。请查该构建版本的完整序列，再决定是否需要处理。',
    'germline_endpoint_boundary': '原始 FASTA 的相邻序列支持另一种域边界，请核对编号片段与原始记录的交界。',
    'downstream_motif_overlaps_numbering': '标签、串联 G4S 接头或恒定区的候选起点落在已编号的 FR4 内，请对照原始构建图核实交界。',
    'ambiguous_numbered_start_methionine': '原文件首位 M 被编号进 FR1，与候选 V 起点不同。有意引入 M 和缺少原首位后接起始 M 均可产生此序列，请核对表达构建记录。',
    'original_unverified': '原始文件缺失或已变更，只能核对保存的转换记录。',
    'original_numbering_unavailable': '缺少原始编号记录，无法交叉核对提取结果。',
    'source_record_unavailable': '未找到对应的原始序列记录。',
    'nonconsecutive_source_positions': '映射中的原始序列位置不连续，请核对提取范围。',
    'source_sequence_mismatch': '映射中的氨基酸与原始序列不一致。',
    'numbering_mapping_inconsistent': '残基映射中的编号记录不一致。',
    'region_definition_mismatch': '映射中的区域归属与所选划分规则不符。',
    'framework_header_count': '该链此 FR 的 YAML 序列注释缺失或重复。',
    'framework_empty': '该 FR 在导出的残基映射中为空。',
    'nonstandard_framework_residue': '该 FR 含未知或非标准氨基酸，请核对原始记录。',
    'unexpected_output_chains': '残基映射包含提取记录以外的链。',
    'structure_read_failed': '原结构文件读取失败，无法核对坐标。',
    'conversion_model_unknown': '无法确认转换使用的结构模型。',
    'source_chain_not_matched': '无法唯一确定对应的原结构链。',
    'sequence_yaml_payload_invalid': '无法解析导出的蛋白质定义以核对固定序列和设计区。',
}

def _localize_finding(f):
    evidence = f['evidence']
    code = f['code']
    if code == 'known_reference_terminal_difference':
        zh = f"相对 {evidence['reference']}，提取域的 {evidence['side']} 端少 {evidence['residues']}。请核实所需版本。"
    elif code == 'reference_residues_outside_numbered_domain':
        zh = f"参考版本的 {evidence['side']} 端残基 {evidence['residues']} 位于上传序列中、提取域之外。请核对域边界。"
    else:
        zh = _DESCRIPTIONS_ZH[code]
    en = f['message']
    if code == 'structure_read_failed':
        en = f"Could not read the source structure to check coordinates. {en}"
    if code == 'germline_endpoint_short':
        extra = '/'.join(x for x in evidence.get('candidate_extensions', []) if x)
        if extra:
            zh += f" 候选延伸残基：{extra}。"
            en += f" Candidate extension: {extra}."
    if code == 'atypical_imgt_length':
        lo, hi = evidence['typical']
        zh += f" 实际 {evidence['length']} aa，典型范围 {lo}–{hi} aa。"
        en += f" Observed {evidence['length']} aa, typical range {lo}–{hi} aa."
    if code == 'downstream_motif_overlaps_numbering':
        sites = ', '.join(f"{h['motif']}@{h['source_start']}" for h in evidence.get('downstream_motif_starts', []))
        if sites:
            zh += f" 候选片段及原序列起点：{sites}。"
            en += f" Candidate motifs and source start positions: {sites}."
    if evidence.get('sites'):
        positions = ', '.join(f"{s.get('aa', '?')}{s.get('source_position', '?')}" for s in evidence['sites'])
        zh += f" 原始序列位置：{positions}。"
        en += f" Source positions: {positions}."
    if code == 'noncanonical_anchor':
        detail = f"IMGT {evidence['imgt_label']} ({evidence['aa']}, {'/'.join(evidence['typical'])})"
        zh += f" 位点（实测，典型）：{detail}。"
        en += f" Position (observed, typical): {detail}."
    f.update(message_zh=zh, message_en=en)


def _risk_summary(actionable, incomplete):
    if actionable:
        status = 'review'
        zh = '检查发现需复核的项目，请按下列位置和证据核对原始文件或构建记录。'
        en = 'The checks flagged items for review. Verify the listed locations and evidence against the original file or construct record.'
    elif incomplete:
        status = 'insufficient_evidence'
        zh = '已完成的检查未发现明确异常，部分位置证据不足。请自行核对对应位置及构建记录。'
        en = 'The completed checks found no specific anomaly. Evidence is insufficient at some locations. Please verify them against the construct record.'
    else:
        status = 'no_flag'
        zh = '基于内置参考库检索及本次可执行的检查，未发现骨架异常。请自行对照原始构建记录进一步核对。'
        en = 'The built-in reference search and applicable checks found no framework anomaly. Please independently verify the result against the original construct record.'
    return dict(status=status, message_zh=zh, message_en=en)


def build_postcheck(directory, source=None):
    directory = Path(directory)
    for name in ('metadata.json', 'scaffold.yaml', 'residue_map.jsonl'):
        if not (directory / name).is_file():
            raise ScaffoldError(f'{directory}: missing conversion artifact {name}')
    meta = json.loads((directory / 'metadata.json').read_text(encoding='utf-8-sig'))
    required = ('chains', 'cdr_definition', 'numbering_scheme', 'scaffold_kind')
    if not isinstance(meta, dict) or not all(meta.get(k) for k in required):
        raise ScaffoldError('Conversion metadata lacks chains, scaffold kind or region definition')
    rows = [json.loads(line) for line in (directory / 'residue_map.jsonl').read_text(encoding='utf-8-sig').splitlines() if line.strip()]
    findings = []
    original = _verified_source(meta, source, findings)
    raw = read_fasta(original) if original and original.suffix.lower() in {'.fasta', '.fa', '.faa'} else None
    yaml_path = directory / 'scaffold.yaml'
    sequence = postcheck_sequences(meta, rows, yaml_path, raw)
    structure = postcheck_coordinates(original, meta, rows)
    versions = compare_terminals(meta)
    findings.extend(sequence['findings'] + structure.pop('findings') + versions.pop('findings'))
    endpoints = _fasta_endpoints(meta, rows, raw, versions, findings)
    findings.extend(_check_sequence_payload(yaml_path, meta, rows))
    previous_validation = _validation_evidence(directory, findings)
    counts = dict(Counter(f['status'] for f in findings))
    actionable = any(key in counts for key in ('observed', 'suspected', 'reference_relative'))
    for assessment in endpoints['assessments']:
        reason = assessment['evidence']['reason']
        assessment['reason_zh'] = _ENDPOINT_REASONS_ZH.get(reason, reason)
        assessment['reason_en'] = reason.replace('_', ' ').replace('imgt', 'IMGT')
    for item in findings:
        _localize_finding(item)
    incomplete = counts.get('not_assessable', 0) or any(a['state'] == 'indeterminate' for a in endpoints['assessments'])
    return dict(schema_version=2, risk_summary=_risk_summary(actionable, incomplete), conclusion='review' if actionable else 'no_flag_in_checks_performed',
                input=meta.get('input_file'), original_verified=original is not None,
                scheme_preset=meta.get('scheme_preset'), numbering_scheme=meta['numbering_scheme'],
                cdr_definition=meta['cdr_definition'], scaffold_kind=meta['scaffold_kind'],
                previous_validation=previous_validation, frameworks=sequence['frameworks'], coordinates=structure,
                terminal_comparison=versions, germline_endpoints=endpoints, counts=counts, findings=findings)


def _markdown(report):
    status = {'observed': '文件中已观察到', 'reference_relative': '相对参考有差异',
              'suspected': '建议复核', 'not_assessable': '证据不足'}
    coordinate_status = {'assessed': '已核对', 'partially_assessed': '部分可核对', 'not_assessable': '无法评估'}
    def cell(text):
        return str(text).replace('|', '&#124;').replace('\n', ' ')
    coord = report['coordinates']
    coordinate_text = coordinate_status[coord['status']]
    if coord['summary'].get('reason') == 'no_source_structure' and report['original_verified']:
        coordinate_text = '本次上传的是序列，没有可核对的三维坐标'
    lines = ['# 骨架序列检查', '',
             report['risk_summary']['message_zh'], '', report['risk_summary']['message_en'], '',
             f"输入文件 `{report['input']}`  ",
             f"编号 `{report['numbering_scheme']}`，区域定义 `{report['cdr_definition']}`。原始文件摘要核对{'成功' if report['original_verified'] else '不可用'}。  ",
             f"原结构坐标覆盖：{coordinate_text}。结论仅覆盖成功提取的可变域。",
             '', '| 输出链 | 区域 | 从来源提取的序列 | 提取长度 | 实际纳入骨架长度 |',
             '|---|---|---|---:|---:|']
    for fr in report['frameworks']:
        lines.append(f"| {cell(fr['chain'])} | {fr['region']} | `{cell(fr['sequence'])}` | {fr['length']} | {fr['included_length']} |")
    lines.extend(['', '## 需复核的位置', ''])
    if not report['findings']:
        uncertain = [a for a in report['germline_endpoints']['assessments'] if a['state'] == 'indeterminate']
        lines.append('本次检查未发现具体异常。部分末端证据不足，详见下方端点表。' if uncertain
                     else '本次可执行的检查没有发现异常。')
    for f in report['findings']:
        evidence = f['evidence']
        message = f['message_zh']
        lines.append(f"- **{cell(f['chain'] or '整体')} {f['region'] or ''} · {status.get(f['status'], f['status'])}** {message}")
        lines.append(f"  **{cell(f['chain'] or 'Overall')} {f['region'] or ''} · English** {cell(f['message_en'])}")
        if 'reference' in evidence:
            lines[-1] += f" [参考来源]({evidence['source_url']})"
        elif evidence:
            detail = cell(json.dumps(evidence, ensure_ascii=False))
            if len(detail) <= 240:
                lines[-1] += ' `' + detail + '`'
            else:
                lines[-1] += f" 逐位差异见 JSON 中的 `{f['code']}`。"
    lines.extend(['', '## 末端参考', ''])
    for comparison in report['terminal_comparison']['comparisons']:
        exact = '、'.join(comparison['exact_versions']) or '无完整域精确匹配'
        lines.append(f"- {cell(comparison['chain'])}：{exact}。")
    if report['germline_endpoints']['assessments']:
        labels = {'suspected_short': '比候选参考短，建议核对', 'boundary_review': '建议复核域边界',
                  'compatible': '与参考终点相容', 'indeterminate': '无法区分或证据不足'}
        lines.extend(['', '## FASTA 端点检查', '',
                      '相容仅表示参考支持这个长度。参考库只覆盖部分物种与构建版本，候选相似度不能确定物种。',
                      '', '| 链 | 末端 | 结果 | 判断依据 | 原始 FASTA 域范围（从 1 起） | 域外相邻片段 | 候选延伸 |',
                      '|---|---|---|---|---|---|---|'])
        for a in report['germline_endpoints']['assessments']:
            evidence = a['evidence']
            reason = evidence.get('reason', '')
            explanation = a['reason_zh']
            junctions = evidence.get('downstream_motif_starts', [])
            if junctions:
                starts = '、'.join(f"{hit['motif']}@{hit['source_start']}" for hit in junctions)
                explanation += f'。原文件中的候选起点：{starts}'
            if reason in {'same_length_terminal_residue_disagrees', 'endpoint_length_only_residue_differs'}:
                explanation += f"（输入 {evidence['observed_endpoint_residue']}，候选 {'/'.join(evidence['reference_endpoint_residues'])}）"
            span = '–'.join(map(str, evidence.get('raw_domain_span', []))) or '未知'
            flank = evidence.get('adjacent_sequence') or ('已到文件边界' if evidence.get('adjacent_length') == 0 else '未知')
            extensions = [value or '原位终止' for value in evidence.get('candidate_extensions', [])]
            if evidence.get('shorter_reference_present'):
                extensions.append('存在更短参考')
            lines.append(f"| {cell(a['chain'])} | {a['side']} | {labels[a['state']]} / {a['state'].replace('_', ' ')} | {cell(explanation)} / {cell(a['reason_en'])} | {span} | {cell(flank)} | {cell(' / '.join(extensions) or '—')} |")
    lines.append('\n库外序列的末端完整性仍待核实。保守位点和长度提示用于复核序列，'
                 '表达、稳定性、结合活性需要另外评估。逐项证据见同目录 framework_check.json。')
    return '\n'.join(lines) + '\n'


def _replace(path, text):
    temp = path.with_name(f'.{path.name}.{uuid4().hex}.tmp')
    try:
        temp.write_text(text, encoding='utf-8')
        temp.replace(path)
    finally:
        temp.unlink(missing_ok=True)


def write_postcheck(directory, source=None):
    directory = Path(directory).resolve()
    report = build_postcheck(directory, source)
    markdown = _markdown(report)
    _replace(directory / 'framework_check.json', json.dumps(report, ensure_ascii=False, indent=2) + '\n')
    _replace(directory / 'framework_check.md', markdown)
    return dict(conclusion=report['conclusion'], risk_summary=report['risk_summary'], findings_by_evidence=report['counts'],
                fasta_endpoint_states=dict(Counter(a['state'] for a in report['germline_endpoints']['assessments'])),
                coordinate_assessment=report['coordinates']['status'], report=str(directory / 'framework_check.md'))

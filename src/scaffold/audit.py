"""Read-only IMGT framework integrity audit; never repairs or trims an input."""
from __future__ import annotations
from collections import Counter, defaultdict
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
import csv
import json
import re
from Bio import Align, SeqIO
from .errors import ScaffoldError
from .numbering import Numberer

AA = set('ACDEFGHIKLMNPQRSTVWY')
BOUNDS = ((1,26,'FR1'),(27,38,'CDR1'),(39,55,'FR2'),(56,65,'CDR2'),(66,104,'FR3'),(105,117,'CDR3'),(118,129,'FR4'))
TYPICAL = {'FR1':(25,26),'FR2':(16,17),'FR3':(36,39),'FR4':(10,12)}
IMGT_SOURCE = 'https://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html'


def risk(code, message, *, severity='warning', certainty='suspected', region=None,
         positions=None, labels=None, evidence=None, recommendation=None):
    return {'code':code,'severity':severity,'certainty':certainty,'region':region,
            'sequence_positions':positions or [],'imgt_labels':labels or [],
            'message':message,'evidence':evidence or {},'recommendation':recommendation or 'Review the stated evidence; do not automatically alter the sequence.'}


def region_for(position):
    return next((name for lo,hi,name in BOUNDS if lo <= position <= hi), 'UNASSIGNED')


@dataclass
class AuditDomain:
    record_id: str
    domain_id: str
    chain_type: str
    start: int
    end: int
    sequence: str
    rows: list[dict]
    score: float | None


def number_domains(record_id, sequence, numberer):
    """Use raw numbering so incomplete domains are not rejected by conversion rules."""
    pending=[(0,len(sequence))]
    output=[]
    while pending:
        left,right=pending.pop(0)
        if right-left < 30:
            continue
        fragment=sequence[left:right]
        located=[]
        for positions,chain_type,score,hint in numberer._number_one(fragment):
            real=[(p,aa) for p,aa in positions if aa != '-']
            text=''.join(aa for _,aa in real)
            if not text:
                continue
            starts=[m.start() for m in re.finditer(f'(?={re.escape(text)})',fragment)]
            if hint is not None and int(hint) in starts:
                start=int(hint)
            elif len(starts)==1:
                start=starts[0]
            else:
                raise ScaffoldError(f'{record_id}: numbered domain cannot be located unambiguously')
            start+=left
            end=start+len(text)
            rows=[]
            seen=set()
            for offset,((pos,ins),aa) in enumerate(real):
                pos=int(pos); ins=str(ins).strip(); label=f'{pos}{ins}'
                if label in seen:
                    raise ScaffoldError(f'{record_id}: duplicate IMGT label {label}')
                seen.add(label)
                rows.append({'input_seq_index':start+offset+1,'domain_seq_index':offset+1,
                             'aa':aa,'imgt_position':pos,'imgt_label':label,'region':region_for(pos)})
            output.append(AuditDomain(record_id,'',chain_type,start,end,text,rows,float(score) if score is not None else None))
            located.append((start,end))
        if located:
            cursor=left
            for start,end in sorted(located):
                if start < cursor:
                    raise ScaffoldError(f'{record_id}: overlapping domain predictions')
                if start-cursor>=30:
                    pending.append((cursor,start))
                cursor=end
            if right-cursor>=30:
                pending.append((cursor,right))
    output.sort(key=lambda d:d.start)
    for i,d in enumerate(output,1):
        d.domain_id=f'{record_id}:domain{i}'
    return output


def read_audit_fasta(path, selected=None):
    records={}; problems=[]
    for item in SeqIO.parse(str(path),'fasta'):
        if selected and item.id not in selected:
            continue
        if item.id in records:
            raise ScaffoldError(f'Duplicate FASTA record ID: {item.id}')
        seq=str(item.seq).upper()
        records[item.id]=seq
        invalid=[{'position':i+1,'symbol':aa} for i,aa in enumerate(seq) if aa not in AA]
        if invalid or not seq:
            problems.append(risk('invalid_sequence_symbols',f'{item.id}: empty or non-standard amino-acid sequence; this record is not numbered.',
                                 severity='high',certainty='observed',evidence={'record_id':item.id,'invalid':invalid},
                                 recommendation='Check unknown residues, alignment gaps, stop symbols, and sequence type. Other valid records are still audited.'))
    if selected and set(selected)-set(records):
        raise ScaffoldError(f'Unknown FASTA records: {sorted(set(selected)-set(records))}')
    if not records:
        raise ScaffoldError('No selected FASTA records')
    return records,problems


def align_pair(reference, query):
    """Global affine-gap comparison of two variable domains, with ambiguity control.

    Exact contiguous matches use a documented terminal-deletion convention.
    Otherwise every equally optimal mapping must agree (up to a bounded check).
    """
    starts=[m.start() for m in re.finditer(f'(?={re.escape(query)})',reference)] if query else []
    if len(starts)==1:
        mapping={i:starts[0]+i for i in range(len(query))}
        gaps=[starts[0],len(reference)-starts[0]-len(query)]
        score=2*len(query)+sum(-6-0.5*(length-1) for length in gaps if length)
        return {'mapping':mapping,'ambiguous':False,'identity':1.0,'score':score,
                'method':'unique_exact_contiguous_match; terminal deletions placed outside this match',
                'missing_reference':[i for i in range(len(reference)) if i not in mapping.values()], 'insertions':[], 'mismatches':[]}
    aligner=Align.PairwiseAligner(mode='global',match_score=2,mismatch_score=-3,open_gap_score=-6,extend_gap_score=-0.5)
    alignments=aligner.align(reference,query)
    mappings=[]
    score=None
    for i,alignment in enumerate(alignments):
        if i>=16:
            break
        score=float(alignment.score)
        mapping={}
        for (r0,r1),(q0,q1) in zip(*alignment.aligned):
            mapping.update({q:r for q,r in zip(range(int(q0),int(q1)),range(int(r0),int(r1)))})
        mappings.append(mapping)
    if not mappings:
        return {'mapping':{},'ambiguous':True,'identity':0,'score':float('-inf'),'method':'no_alignment'}
    first=mappings[0]
    ambiguous=len(mappings)>=16 or any(m!=first for m in mappings[1:])
    mismatches=[{'reference_pos1':r+1,'query_pos1':q+1,'reference_aa':reference[r],'query_aa':query[q]} for q,r in first.items() if reference[r]!=query[q]]
    return {'mapping':first,'ambiguous':ambiguous,'identity':(len(first)-len(mismatches))/max(1,len(first)),
            'score':score,'method':'global_affine_gap_alignment',
            'missing_reference':[i for i in range(len(reference)) if i not in first.values()],
            'insertions':[i for i in range(len(query)) if i not in first], 'mismatches':mismatches}


def choose_match(query, candidates):
    scored=[]
    for candidate in candidates:
        if candidate.chain_type!=query.chain_type:
            continue
        a=align_pair(candidate.sequence,query.sequence)
        coverage=len(a['mapping'])/max(1,len(query.sequence))
        if a['identity']>=0.9 and coverage>=0.8:
            scored.append((a['score'],candidate,a))
    if not scored:
        return None,'no_sufficient_match'
    scored.sort(key=lambda item:item[0],reverse=True)
    ties=[x for x in scored if abs(x[0]-scored[0][0])<1e-6]
    if len({x[1].sequence for x in ties})>1:
        return None,'multiple_equally_supported_references'
    best=ties[0]
    if best[2]['ambiguous']:
        return None,'multiple_optimal_residue_mappings'
    return (best[1],best[2]),None


def _reference_risks(domain, references, input_sequence):
    # For a short record, compare the actual supplied sequence, not just the
    # model's extracted fragment. Truncated domains can leave real AA unnumbered.
    short_refs=[r for r in references if r.chain_type==domain.chain_type and len(input_sequence)<=len(r.sequence)+25]
    comparison=domain
    scope='numbered_variable_domain_only'
    if short_refs:
        comparison=AuditDomain(domain.record_id,domain.domain_id,domain.chain_type,0,len(input_sequence),input_sequence,domain.rows,domain.score)
        scope='complete_supplied_record_including_unnumbered_flanks'
    match,reason=choose_match(comparison,short_refs or references)
    if match is None:
        return [risk('ambiguous_reference_mapping' if reason!='no_sufficient_match' else 'reference_not_matched',
                     'The supplied reference cannot be assigned unambiguously to this domain.',certainty='not_assessable',
                     evidence={'reason':reason})], {'status':'not_assessable','reason':reason}
    ref,a=match
    result=[]
    prefix=[];suffix=[];internal=[]
    used=list(a['mapping'].values())
    lo,hi=min(used),max(used)
    for index in a['missing_reference']:
        (prefix if index<lo else suffix if index>hi else internal).append(index)
    for side,indices,code in [('N',prefix,'reference_terminal_deletion'),('C',suffix,'reference_terminal_deletion'),('internal',internal,'reference_internal_deletion')]:
        if not indices:
            continue
        by_region=defaultdict(list)
        for index in indices:
            by_region[ref.rows[index]['region']].append(index)
        for region,subset in by_region.items():
            result.append(risk(code,f'{len(subset)} residue(s) absent relative to the supplied reference ({side} segment). This may be an intentional construct difference.',
                               severity='high' if region.startswith('FR') else 'warning',certainty='reference_relative',region=region,
                               labels=[ref.rows[i]['imgt_label'] for i in subset],
                               evidence={'reference_id':ref.domain_id,'reference_positions':[ref.start+i+1 for i in subset],
                                                      'reference_residues':''.join(ref.sequence[i] for i in subset),'side':side,'alignment_method':a['method'],'comparison_scope':scope},
                               recommendation='Confirm that this reference is the intended exact construct before restoring any residues; the input is not modified.'))
    if a['mismatches'] or a['insertions']:
        result.append(risk('reference_mismatch','Sequence substitutions/insertions differ from the supplied reference; these are not automatically errors.',
                           certainty='reference_relative',evidence={'reference_id':ref.domain_id,'mismatches':a['mismatches'],'query_insertion_positions':[comparison.start+i+1 for i in a['insertions']],'comparison_scope':scope}))
    return result,{'status':'exact_variable_domain_match' if comparison.sequence==ref.sequence else 'differences_relative_to_supplied_reference',
                   'reference_id':ref.domain_id,'identity_on_aligned_residues':a['identity'],'alignment_method':a['method'],'comparison_scope':scope}


def _assess(domain, input_sequence, coordinate_chain=None, coordinate_map=None, coordinate_mapping_reliable=True, references=None, paired_evidence=None, coordinate_uncertain=None):
    risks=[]
    rows=[]
    for original in domain.rows:
        row=dict(original)
        index=row['input_seq_index']-1
        residue=(coordinate_chain.residues.get(index) if coordinate_map is None else coordinate_map.get(row['domain_seq_index']-1)) if coordinate_chain else None
        uncertain=(index in coordinate_chain.uncertain_sequence_indices if coordinate_map is None else row['domain_seq_index']-1 in (coordinate_uncertain or set())) if coordinate_chain else False
        assessed=coordinate_chain is not None and coordinate_mapping_reliable and not uncertain
        row.update(has_positive_ca=bool(residue and residue.has_positive_ca) if assessed else None,
                   source_auth_seq_id=residue.auth_seq_id if residue else None,
                   source_insertion_code=residue.insertion_code if residue else None,
                   source_label_seq_id=residue.label_seq_id if residue else None,
                   coordinate_state='not_assessable_no_structure' if coordinate_chain is None else 'mapping_uncertain' if not assessed else
                       'missing_residue' if residue is None else 'zero_occupancy' if not residue.has_any_positive_atoms else
                       'missing_backbone_atoms' if residue.missing_backbone else 'backbone_present')
        if residue:
            row['missing_backbone_atoms']=list(residue.missing_backbone)
            row['zero_occupancy_atoms']=list(residue.zero_occupancy_atoms)
            row['source_resname']=residue.resname
        rows.append(row)
    regions={}
    for name,(minimum,maximum) in TYPICAL.items():
        selected=[r for r in rows if r['region']==name]
        coords=coordinate_chain is not None and coordinate_mapping_reliable
        regions[name]={'sequence':''.join(r['aa'] for r in selected),'length':len(selected),
                       'input_positions':[r['input_seq_index'] for r in selected], 'imgt_labels':[r['imgt_label'] for r in selected],
                       'coordinate_assessed':coords,'coordinate_mapping_uncertain_count':sum(r['coordinate_state']=='mapping_uncertain' for r in selected),
                       'missing_positive_CA_count':sum(r['has_positive_ca'] is False for r in selected) if coords else None,
                       'sequence_completeness':'not_established_without_construct_reference'}
        if not selected:
            risks.append(risk('incomplete_domain_numbering',f'{name} has no assigned residues; the domain may be partial or numbering may have failed.',
                              severity='high',certainty='suspected',region=name))
        elif len(selected)<minimum or len(selected)>maximum:
            risks.append(risk('outside_typical_fr_length',f'{name} length {len(selected)} is outside the typical IMGT range {minimum}..{maximum}; this alone does not establish a deletion.',
                              region=name,positions=regions[name]['input_positions'],evidence={'reference':IMGT_SOURCE}))
        if coords:
            for state,code in [('missing_residue','missing_residue_coordinates'),('zero_occupancy','zero_occupancy_coordinates'),('missing_backbone_atoms','missing_backbone_atoms')]:
                bad=[r for r in selected if r['coordinate_state']==state]
                if bad:
                    risks.append(risk(code,f'{name}: {len(bad)} residue(s) have unavailable/incomplete backbone coordinates relative to the sequence used for this audit.',
                                      severity='high',certainty='observed',region=name,
                                      positions=[r['input_seq_index'] for r in bad],labels=[r['imgt_label'] for r in bad],
                                      evidence={'sequence_basis':coordinate_chain.sequence_basis,'paired_sequence':paired_evidence,
                                                'details':[{'position':r['input_seq_index'],'missing_backbone':r.get('missing_backbone_atoms',[])} for r in bad]},
                                      recommendation='Retain the sequence. Obtain a matching complete structure or explicitly model and separately validate these coordinates.'))
    first,last=domain.rows[0],domain.rows[-1]
    if domain.start==0 and first['imgt_position']>1:
        risks.append(risk('possible_n_terminal_truncation','The input starts after IMGT position 1; this is a screening hint, not proof that the intended construct is truncated.',
                          region='FR1',evidence={'first_imgt_label':first['imgt_label']},recommendation='Compare with the exact construct sequence, SEQRES/entity sequence, or a user-supplied complete reference.'))
    if domain.end==len(input_sequence) and (last['imgt_position']<118 or regions['FR4']['length']<10):
        risks.append(risk('possible_c_terminal_truncation','The input ends before a typical complete FR4; terminal truncation or a numbering problem is possible.',
                          region='FR4',evidence={'last_imgt_label':last['imgt_label'],'FR4_length':regions['FR4']['length']}))
    if domain.start>0 or domain.end<len(input_sequence):
        risks.append(risk('unnumbered_flanks','Residues outside this variable domain were retained in the input but excluded from FR extraction; they may be tags, signal peptides, constant domains or linkers.',
                          severity='info',certainty='observed',evidence={'before_domain':domain.start,'after_domain':len(input_sequence)-domain.end}))
    by_label={r['imgt_label']:r for r in rows}
    for label,allowed in [('23','C'),('41','W'),('89','ILVFM'),('104','C'),('118','FW'),('119','G'),('121','G')]:
        r=by_label.get(label)
        if r and r['aa'] not in allowed:
            risks.append(risk('noncanonical_anchor',f"IMGT {label} is {r['aa']} (typical: {allowed}); known engineered scaffolds can legitimately differ.",
                              region=r['region'],positions=[r['input_seq_index']],labels=[label],evidence={'observed_aa':r['aa'],'typical':allowed},
                              recommendation='Verify construct identity; do not replace this residue solely to match a consensus.'))
    reference={'status':'not_provided','terminal_completeness':'not_determined'}
    if references is not None:
        found,reference=_reference_risks(domain,references,input_sequence)
        risks.extend(found)
        for info in regions.values():
            info['sequence_completeness']=reference['status']
    if coordinate_chain is None:
        risks.append(risk('coordinates_not_supplied','FASTA contains no coordinates; coordinate coverage was not assessed.',severity='info',certainty='not_assessable'))
    elif not coordinate_mapping_reliable or any(r['coordinate_state']=='mapping_uncertain' for r in rows):
        risks.append(risk('ambiguous_coordinate_mapping','Sequence-to-coordinate mapping is incomplete or ambiguous; precise missing-coordinate positions are withheld.',severity='high',certainty='not_assessable'))
    elif coordinate_chain.sequence_basis=='coordinates_only' and not paired_evidence:
        risks.append(risk('coordinate_only_sequence','No declared full sequence is available. Missing residues absent from ATOM records may be invisible; author numbering gaps alone are not proof of deletions.',severity='warning',certainty='not_assessable'))
    return {'domain_id':domain.domain_id,'record_id':domain.record_id,'chain_type':domain.chain_type,
            'identity_note':'H denotes VH/VHH; sequence numbering alone does not establish VHH origin.' if domain.chain_type=='H' else 'K/L denotes the light-chain subtype.',
            'input_start':domain.start+1,'input_end':domain.end,'sequence':domain.sequence,'numbering_score':domain.score,
            'structure_chain_auth':coordinate_chain.auth_id if coordinate_chain else None,
            'structure_chain_label':coordinate_chain.label_id if coordinate_chain else None,
            'structure_mapping':{'status':coordinate_chain.mapping_status,'sequence_basis':coordinate_chain.sequence_basis,'warnings':coordinate_chain.warnings} if coordinate_chain else None,
            'coordinate_summary':{'numbered_residues':len(rows),
                                  'missing_positive_CA_count':sum(r['has_positive_ca'] is False for r in rows) if coordinate_chain and coordinate_mapping_reliable else None,
                                  'uncertain_positions':[r['input_seq_index'] for r in rows if r['coordinate_state']=='mapping_uncertain'],
                                  'scope':'All numbered FR/CDR residues; region summaries and risk rules below focus on FR1-FR4.'},
            'regions':regions,'rows':rows,'risks':risks,'reference_comparison':reference}


def _number_records(records, numberer, problems):
    output=[]
    for rid,seq in records.items():
        if not seq or set(seq)-AA:
            continue
        try:
            ds=number_domains(rid,seq,numberer)
        except (ScaffoldError,ValueError,RuntimeError) as exc:
            problems.append(risk('numbering_failed',f'{rid}: {exc}',severity='high',certainty='not_assessable',evidence={'record_id':rid}))
            continue
        if not ds:
            problems.append(risk('antibody_domain_not_detected',f'{rid}: no antibody variable domain could be assessed; this is not a completeness pass.',certainty='not_assessable',evidence={'record_id':rid}))
        output.extend(ds)
    return output


def audit_file(path, numberer=None, sequence_path=None, reference_path=None, selected_chains=None,
               namespace='auth', model_index=1, records=None):
    path=Path(path)
    if not path.is_file():
        raise ScaffoldError(f'Input not found: {path}')
    numberer=numberer or Numberer('anarcii',2,'imgt')
    problems=[]; domains=[]; source=None
    reference_domains=None
    if reference_path:
        refs,ref_problems=read_audit_fasta(reference_path)
        problems.extend(ref_problems)
        reference_domains=_number_records(refs,numberer,problems)
        incomplete=[d for d in reference_domains if not set(TYPICAL)<={r['region'] for r in d.rows}]
        for d in incomplete:
            problems.append(risk('incomplete_reference_domain',f'{d.domain_id}: the proposed complete reference lacks one or more numbered FR regions; it was excluded from precise reference comparisons.',severity='high',certainty='not_assessable'))
        reference_domains=[d for d in reference_domains if d not in incomplete]
    is_structure=path.suffix.lower() in {'.pdb','.cif','.mmcif'}
    if not is_structure and path.suffix.lower() not in {'.fa','.faa','.fasta'}:
        raise ScaffoldError(f'Unsupported audit input: {path.suffix}')
    if sequence_path and not is_structure:
        raise ScaffoldError('--sequence pairs a FASTA with one PDB/CIF input')
    if selected_chains and not is_structure:
        raise ScaffoldError('Use --records for FASTA; chain selection applies to PDB/CIF')
    if is_structure:
        from .audit_structure import read_audit_structure
        source=read_audit_structure(path,model_index,selected_chains,namespace)
        input_records={key:c.sequence for key,c in source.chains.items()}
        for key,chain in source.chains.items():
            for warning in chain.warnings:
                problems.append(risk('structure_sequence_mapping_warning',f"{key}: {warning.get('message',warning['code'])}",
                                     certainty='observed',evidence={'chain':key,**warning}))
            if not chain.sequence or set(chain.sequence)-AA:
                problems.append(risk('invalid_sequence_symbols',f'{key}: structure sequence contains non-standard or unknown residues; numbering was not attempted.',severity='high',certainty='observed',evidence={'record_id':key}))
        structure_domains=_number_records(input_records,numberer,problems)
        if sequence_path:
            paired_records,pair_problems=read_audit_fasta(sequence_path,records)
            problems.extend(pair_problems)
            paired_domains=_number_records(paired_records,numberer,problems)
            used_domains=set()
            for structural in structure_domains:
                chain=source.chains[structural.record_id]
                # Match short complete records, including any residues that a
                # partial-domain prediction left unnumbered. Keep IMGT rows from
                # the original prediction and translate the coordinate offsets.
                structural_compare=structural
                if sum(d.record_id==structural.record_id for d in structure_domains)==1 and any(len(chain.sequence)<=len(d.sequence)+25 for d in paired_domains if d.chain_type==structural.chain_type):
                    structural_compare=AuditDomain(structural.record_id,structural.domain_id,structural.chain_type,0,len(chain.sequence),chain.sequence,structural.rows,structural.score)
                pair_candidates=[]
                for candidate in paired_domains:
                    seq=paired_records[candidate.record_id]
                    if sum(d.record_id==candidate.record_id for d in paired_domains)==1 and len(seq)<=len(structural_compare.sequence)+25:
                        candidate=AuditDomain(candidate.record_id,candidate.domain_id,candidate.chain_type,0,len(seq),seq,candidate.rows,candidate.score)
                    pair_candidates.append(candidate)
                match,reason=choose_match(structural_compare,pair_candidates)
                if match is None:
                    problems.append(risk('paired_sequence_not_matched','Cannot assign the supplied FASTA unambiguously to this structural domain.',severity='high',certainty='not_assessable',evidence={'structure_domain':structural.domain_id,'reason':reason}))
                    continue
                comparison,a=match
                main=next(d for d in paired_domains if d.domain_id==comparison.domain_id)
                used_domains.add(main.domain_id)
                coord_map={comparison.start+ref_index-main.start:chain.residues.get(structural_compare.start+query_index) for query_index,ref_index in a['mapping'].items()}
                reliable=not a['mismatches']
                uncertain={comparison.start+ref_index-main.start for query_index,ref_index in a['mapping'].items() if structural_compare.start+query_index in chain.uncertain_sequence_indices}
                assessed=_assess(main,paired_records[main.record_id],chain,coord_map,reliable,reference_domains,
                                 {'file':str(Path(sequence_path).resolve()),'record_id':main.record_id,'alignment_method':a['method']},uncertain)
                assessed['domain_id']=f'{structural.domain_id}|FASTA:{main.domain_id}'
                if a['missing_reference'] or a['insertions'] or a['mismatches']:
                    assessed['risks'].append(risk('sequence_structure_disagreement','Supplied FASTA and structure-derived sequences differ; coordinate gaps are interpreted relative to the supplied FASTA, not silently repaired.',
                                                  certainty='reference_relative',evidence={'structure_sequence_basis':chain.sequence_basis,
                                                      'FASTA_positions_absent_from_structure_sequence':[comparison.start+i+1 for i in a['missing_reference']],
                                                      'structure_positions_absent_from_FASTA_sequence':[structural_compare.start+i+1 for i in a['insertions']],
                                                      'structure_residues_absent_from_FASTA_sequence':''.join(structural_compare.sequence[i] for i in a['insertions']),
                                                      'alignment_mismatches':a['mismatches']}))
                if a['insertions'] and chain.declared_sequence is not None:
                    relative,_=_reference_risks(main,[structural],paired_records[main.record_id])
                    for finding in relative:
                        finding['evidence']['reference_source']='structure_declared_sequence'
                        finding['evidence']['structure_file']=str(path.resolve())
                    assessed['risks'].extend(relative)
                domains.append(assessed)
            for d in paired_domains:
                if d.domain_id not in used_domains:
                    problems.append(risk('unmatched_fasta_domain',f'{d.domain_id}: not assigned to an input structure domain.',certainty='not_assessable'))
        else:
            for d in structure_domains:
                chain=source.chains[d.record_id]
                domains.append(_assess(d,input_records[d.record_id],chain,references=reference_domains))
    else:
        input_records,fasta_problems=read_audit_fasta(path,records)
        problems.extend(fasta_problems)
        for d in _number_records(input_records,numberer,problems):
            domains.append(_assess(d,input_records[d.record_id],references=reference_domains))
    all_risks=problems+[r for d in domains for r in d['risks']]
    supporting_files=[]
    for role,file in [('paired_sequence',sequence_path),('construct_reference',reference_path)]:
        if file:
            supporting_files.append({'role':role,'path':str(Path(file).resolve()),'sha256':sha256(Path(file).read_bytes()).hexdigest()})
    return {'schema_version':'1.0','status':'audited' if domains else 'not_assessable','input':str(path.resolve()),
            'input_sha256':sha256(path.read_bytes()).hexdigest(),'scheme':'imgt','scheme_source':IMGT_SOURCE,
            'sequence_input':str(Path(sequence_path).resolve()) if sequence_path else None,
            'reference_input':str(Path(reference_path).resolve()) if reference_path else None,
            'supporting_files':supporting_files,
            'model_index':model_index if source else None,
            'structure_warnings':source.warnings if source else [],
            'coordinate_provenance':source.model_provenance if source else None,
            'domain_count':len(domains),'domains':domains,'risks':problems,
            'risk_counts':dict(Counter(r['code'] for r in all_risks)),
            'limitations':['No residues or coordinates were added, deleted or altered.',
                           'IMGT numbering gaps are not themselves evidence of missing sequence.',
                           'An exact construct reference is needed to establish terminal completeness; a typical FR length does not prove completeness.',
                           'Reference differences do not establish accidental truncation or expression/functional risk.',
                           'Coordinates may be experimental or modeled; occupancy and B-factor alone do not establish their provenance.']}


def write_audit(report, out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    (out/'audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    frs=[];mapping=[];risks=[];fr_fasta=[];vd_fasta=[]
    for r in report['risks']:
        risks.append({'domain_id':'',**r})
    for d in report['domains']:
        vd_fasta.extend([f">{d['domain_id']}",d['sequence']])
        mapping.extend({'domain_id':d['domain_id'],**r} for r in d['rows'])
        risks.extend({'domain_id':d['domain_id'],**r} for r in d['risks'])
        for region,info in d['regions'].items():
            frs.append({'domain_id':d['domain_id'],'chain_type':d['chain_type'],'region':region,**info})
            fr_fasta.extend([f">{d['domain_id']}|{region}",info['sequence']])
    for name,data in [('FR_summary.tsv',frs),('residue_map.tsv',mapping),('risks.tsv',risks)]:
        keys=list(dict.fromkeys(k for r in data for k in r)) or ['status']
        with (out/name).open('w',encoding='utf-8-sig',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=keys,delimiter='\t');writer.writeheader()
            for row in data:
                writer.writerow({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in row.items()})
    (out/'variable_domains.fasta').write_text('\n'.join(vd_fasta)+'\n',encoding='utf-8')
    (out/'FR1_FR2_FR3_FR4.fasta').write_text('\n'.join(fr_fasta)+'\n',encoding='utf-8')
    lines=['# IMGT framework integrity audit','',f"Input: `{report['input']}`",'',f"Status: **{report['status']}**; domains: {report['domain_count']}.",'',
           'This report separates coordinate observations, differences relative to a supplied reference, and unconfirmed screening hints. It does not repair the input.','',
           '| Domain | FR | Length | Coordinate assessment | Missing CA | Sequence completeness |','|---|---|---:|---|---:|---|']
    for row in frs:
        uncertain=row['coordinate_mapping_uncertain_count']
        status='partially assessed' if row['coordinate_assessed'] and uncertain else 'assessed' if row['coordinate_assessed'] else 'not assessable'
        count=str(row['missing_positive_CA_count']) if row['missing_positive_CA_count'] is not None else 'unknown'
        if uncertain:
            count+=f' known; {uncertain} uncertain'
        lines.append(f"| {row['domain_id']} | {row['region']} | {row['length']} | {status} | {count} | {row['sequence_completeness']} |")
    if report['coordinate_provenance'] or report['structure_warnings']:
        lines.extend(['','## Coordinate provenance and input warnings','',
                      '```json',json.dumps({'provenance':report['coordinate_provenance'],'warnings':report['structure_warnings']},ensure_ascii=False,indent=2),'```'])
    lines.extend(['','## Findings',''])
    for r in risks:
        lines.append(f"- **{r['code']}** [{r['severity']}; {r['certainty']}] {r['domain_id']} {r.get('region') or ''}: {r['message']} Positions: {r['sequence_positions']}; IMGT: {r['imgt_labels']}. {r['recommendation']}")
    lines.extend(['','## Limits','']+['- '+x for x in report['limitations']]+['',f'[IMGT region definitions]({IMGT_SOURCE})',''])
    (out/'README_audit.md').write_text('\n'.join(lines),encoding='utf-8-sig')

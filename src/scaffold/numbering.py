"""Antibody domain identification and selectable numbering/CDR projections."""
from __future__ import annotations
from dataclasses import dataclass, asdict
from importlib.metadata import version
from pathlib import Path
import contextlib
import io
import re
from Bio import SeqIO
from .errors import ScaffoldError
from .schemes import resolve_policy,project_numbering

AA = set("ACDEFGHIKLMNPQRSTVWY")
REGIONS = [(1,26,'FR1'),(27,38,'CDR1'),(39,55,'FR2'),(56,65,'CDR2'),(66,104,'FR3'),(105,117,'CDR3'),(118,129,'FR4')]
LINKERS = [r'(?:GGGGS){1,}', r'EGKSSGSGSESKST', r'KESGSVSSEQLAQFRSLD', r'GSAGSAAGSGEF']

def region_for(position: int) -> str:
    for start, end, name in REGIONS:
        if start <= position <= end:
            return name
    raise ScaffoldError(f'Unsupported IMGT base position {position}')

@dataclass
class Domain:
    record_id: str
    domain_id: str
    chain_type: str
    start: int   # original sequence, zero-based inclusive
    end: int     # original sequence, zero-based exclusive
    sequence: str
    rows: list[dict]
    score: float | None

    def to_dict(self):
        return asdict(self)

def read_fasta(path: str | Path) -> dict[str, str]:
    records = {}
    for rec in SeqIO.parse(str(path), 'fasta'):
        if rec.id in records:
            raise ScaffoldError(f'Duplicate FASTA record ID: {rec.id}')
        seq = str(rec.seq).upper()
        if not seq or set(seq) - AA:
            raise ScaffoldError(f'{rec.id}: sequence must contain only the 20 standard amino acids; invalid={sorted(set(seq)-AA)}')
        records[rec.id] = seq
    if not records:
        raise ScaffoldError('No FASTA records')
    return records

class Numberer:
    def __init__(self, engine='anarcii', ncpu=2, scheme='imgt', cdr_definition=None):
        self.engine = engine
        self.ncpu = ncpu
        self.policy = resolve_policy(scheme,cdr_definition)
        self._model = None
        self._cache = {}
        self.version = version(engine)

    def _number_one(self, sequence):
        if sequence in self._cache:
            return self._cache[sequence]
        if self.engine == 'anarcii':
            from anarcii import Anarcii
            with contextlib.redirect_stdout(io.StringIO()):
                if self._model is None:
                    self._model = Anarcii(seq_type='antibody', mode='accuracy', cpu=True, ncpu=self.ncpu, verbose=False)
                results = self._model.number([('query', sequence)])
            found = []
            for item in results.values():
                if not item.get('error') and item.get('numbering') and item.get('chain_type') in ('H','K','L'):
                    found.append((item['numbering'], item['chain_type'], item.get('score'), item.get('query_start')))
        elif self.engine == 'anarci':
            from anarci import anarci
            numbered, details, _ = anarci([('query', sequence)], scheme='imgt', ncpu=self.ncpu)
            found = []
            if numbered and numbered[0]:
                for (positions,start,end), detail in zip(numbered[0],details[0]):
                    if detail['chain_type'] in ('H','K','L'):
                        found.append((positions, detail['chain_type'], detail.get('bitscore'), start))
        else:
            raise ScaffoldError(f'Unknown numbering engine: {self.engine}')
        self._cache[sequence] = found
        return found

    def domains(self, record_id, sequence) -> list[Domain]:
        """Locate domains, then scan unnumbered flanks to recover tandem VH/VL."""
        pending = [(0, len(sequence))]
        found = []
        while pending:
            left, right = pending.pop(0)
            if right-left < 80:
                continue
            fragment = sequence[left:right]
            numbered = self._number_one(fragment)
            located = []
            for positions, chain_type, score, hint in numbered:
                real = [(p,aa) for p,aa in positions if aa != '-']
                domain_seq = ''.join(aa for p,aa in real)
                if len(domain_seq) < 80:
                    continue
                starts = [m.start() for m in re.finditer(re.escape(domain_seq), fragment)]
                if hint is not None and int(hint) in starts:
                    start = int(hint)
                elif len(starts) == 1:
                    start = starts[0]
                else:
                    raise ScaffoldError(f'{record_id}: numbered domain cannot be uniquely mapped back to sequence')
                start += left
                end = start + len(domain_seq)
                rows = []
                projected = project_numbering(positions,chain_type,score,self.policy)
                seen = set()
                for i, ((pos, ins), aa) in enumerate(real):
                    ins = str(ins).strip()
                    identity = (int(pos), ins)
                    if identity in seen:
                        raise ScaffoldError(f'{record_id}: duplicate numbering label {identity}')
                    seen.add(identity)
                    rows.append({'record_id':record_id,'input_seq_index':start+i+1,'aa':aa,'imgt_position':int(pos),'imgt_insertion':ins,'imgt_label':f'{pos}{ins}',**projected[i]})
                # Require all framework regions and all CDRs, but permit terminal truncation.
                if {x['region'] for x in rows} != {r[2] for r in REGIONS}:
                    raise ScaffoldError(f'{record_id}: incomplete or unrecognized variable domain')
                found.append(Domain(record_id,'',chain_type,start,end,domain_seq,rows,float(score) if score is not None else None))
                located.append((start,end))
            if located:
                cursor = left
                for start,end in sorted(located):
                    if start < cursor:
                        raise ScaffoldError(f'{record_id}: overlapping domain predictions')
                    if start-cursor >= 80:
                        pending.append((cursor,start))
                    cursor = end
                if right-cursor >= 80:
                    pending.append((cursor,right))
        found.sort(key=lambda d:d.start)
        for i, domain in enumerate(found,1):
            domain.domain_id = f'{record_id}:domain{i}'
            for row in domain.rows:
                row['domain_id'] = domain.domain_id
                row['chain_type'] = domain.chain_type
        return found

def classify(records: dict[str,str], numberer: Numberer, requested='auto') -> dict:
    domains = [d for rid,seq in records.items() for d in numberer.domains(rid,seq)]
    heavy = [d for d in domains if d.chain_type == 'H']
    light = [d for d in domains if d.chain_type in ('K','L')]
    warnings = []
    linkers = []
    if len(heavy) == len(light) == 1 and len(domains) == 2:
        kind = 'scfv' if heavy[0].record_id == light[0].record_id else 'fab'
        selected = [heavy[0], light[0]]
        if kind == 'scfv':
            first, second = sorted(selected, key=lambda d:d.start)
            sequence = records[first.record_id]
            linker = sequence[first.end:second.start]
            if not 1 <= len(linker) <= 80:
                raise ScaffoldError(f'scFv inter-domain segment has {len(linker)} residues; expected 1..80. Inspect possible fusion.')
            pattern = next((p for p in LINKERS if re.search(p,linker)), None)
            linkers.append({'record_id':first.record_id,'start':first.end+1,'end':second.start,'sequence':linker,'matched_pattern':pattern,'orientation':f'{first.chain_type}-{second.chain_type}','validated_by':'one heavy and one light numbered domain flanking linker'})
    elif len(heavy)==1 and not light and len(domains)==1:
        kind = 'vhh'
        selected = heavy
        if requested == 'auto':
            warnings.append('Single heavy variable domain: VHH candidate; sequence-domain typing alone does not establish VHH origin. Use --type vhh when known.')
    else:
        summary = [(d.record_id,d.chain_type,d.start+1,d.end) for d in domains]
        raise ScaffoldError(f'Cannot select one Fab/scFv/VHH unambiguously: domains={summary}. Select FASTA records or structure chains explicitly.')
    if requested != 'auto' and requested != kind:
        raise ScaffoldError(f'Requested {requested}, identified {kind}')
    used = {d.record_id for d in selected}
    return {'molecule_type':kind,'domains':selected,'linkers':linkers,'warnings':warnings+list(numberer.policy.notes),'ignored_records':[rid for rid in records if rid not in used],'numbering_engine':numberer.engine,'numbering_version':numberer.version,'numbering_converter':'anarcii','numbering_converter_version':version('anarcii'),**numberer.policy.metadata()}

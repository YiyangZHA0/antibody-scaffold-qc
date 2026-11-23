"""Number labels and CDR definitions are separate, explicitly recorded choices.

Boundaries: Andrew C. R. Martin's CDR table (Chothia-numbered AbM/Contact),
Oxford ANARCI ImmunoPDB.py, and IMGT. See docs/SCHEMES.md for source links.
AHo is a numbering scheme; our AHo preset explicitly uses IMGT CDRs.
"""
from dataclasses import dataclass
from copy import deepcopy
from .errors import ScaffoldError

PRESETS=('imgt','kabat','chothia','abm','maccallum','aho')
DEFINITIONS=('imgt','kabat','chothia','abm','maccallum','chothia-consensus')
ALIASES={'contact':'maccallum','mac-callum':'maccallum'}
BASE_NUMBERING={'imgt':'imgt','kabat':'kabat','chothia':'chothia','abm':'chothia','maccallum':'chothia','chothia-consensus':'chothia'}
CDR_BOUNDS={
    'imgt':{'H':((27,38),(56,65),(105,117)),'L':((27,38),(56,65),(105,117))},
    'kabat':{'H':((31,35),(50,65),(95,102)),'L':((24,34),(50,56),(89,97))},
    'chothia':{'H':((26,32),(52,56),(95,102)),'L':((24,34),(50,56),(89,97))},
    'chothia-consensus':{'H':((26,32),(52,56),(96,101)),'L':((26,32),(50,52),(91,96))},
    'abm':{'H':((26,35),(50,58),(95,102)),'L':((24,34),(50,56),(89,97))},
    'maccallum':{'H':((30,35),(47,58),(93,101)),'L':((30,36),(46,55),(89,96))},
}


def normalize_name(value):
    value=value.strip().lower()
    return ALIASES.get(value,value)


@dataclass(frozen=True)
class SchemePolicy:
    preset: str
    numbering_scheme: str
    cdr_definition: str
    cdr_numbering_scheme: str
    notes: tuple[str,...]=()

    def metadata(self):
        ranges={}
        for chain,bounds in CDR_BOUNDS[self.cdr_definition].items():
            ranges[chain]={f'CDR{i}':f'{a}..{b}' for i,(a,b) in enumerate(bounds,1)}
        if self.cdr_definition=='kabat':
            ranges['H']['CDR1']='31..35B (35C and later are FR2)'
        return {'scheme_preset':self.preset,'numbering_scheme':self.numbering_scheme,'cdr_definition':self.cdr_definition,'cdr_numbering_scheme':self.cdr_numbering_scheme,'cdr_ranges':ranges,'scheme_notes':list(self.notes)}


def resolve_policy(scheme='imgt',cdr_definition=None):
    preset=normalize_name(scheme)
    if preset not in PRESETS:
        raise ScaffoldError(f'Unsupported scheme {scheme!r}; choose {", ".join(PRESETS)}')
    numbering='aho' if preset=='aho' else BASE_NUMBERING[preset]
    definition=normalize_name(cdr_definition) if cdr_definition else ('imgt' if preset=='aho' else preset)
    if definition not in DEFINITIONS:
        raise ScaffoldError(f'Unsupported CDR definition {definition!r}; AHo needs a separate CDR definition')
    notes=[]
    if preset in ('abm','maccallum'):
        notes.append(f'{preset} is a CDR-definition preset using Chothia number labels; not a standalone numbering algorithm.')
    if preset=='aho':
        notes.append(f'AHo number labels with explicit {definition} CDR definition; not an independent AHo CDR boundary convention.')
    if definition=='maccallum':
        notes.append('MacCallum/Contact uses published Chothia-numbered boundary intervals, not distances to an antigen in this particular structure.')
    if definition=='chothia':
        notes.append('Chothia CDRs use the historical wider ANARCI definition; --cdr-definition chothia-consensus selects the UCL 2021 consensus interpretation.')
    return SchemePolicy(preset,numbering,definition,BASE_NUMBERING[definition],tuple(notes))


def region_for_label(position,insertion,chain_type,definition):
    """Interpret labels only in the native numbering basis of the definition."""
    if chain_type not in ('H','K','L') or definition not in DEFINITIONS:
        raise ScaffoldError('Unknown chain type or CDR definition')
    if position<1 or (definition=='imgt' and position>129):
        raise ScaffoldError(f'Invalid {definition} position {position}')
    role='H' if chain_type=='H' else 'L'
    for i,(start,end) in enumerate(CDR_BOUNDS[definition][role],1):
        if position<start:
            return f'FR{i}'
        if position<=end:
            # Historical Kabat H1 ends at H35B, not every insertion at H35.
            if definition=='kabat' and role=='H' and position==35 and insertion.strip() not in ('','A','B'):
                return 'FR2'
            return f'CDR{i}'
    return 'FR4'


def project_numbering(imgt_positions,chain_type,score,policy):
    """Project the same real residues; fail on a conversion that loses any AA.

    CDR membership is computed in its own native scheme, then mapped by
    residue order. This avoids applying e.g. Chothia cutoffs to Kabat labels.
    """
    original=[(p,aa) for p,aa in imgt_positions if aa!='-']
    original_sequence=''.join(aa for _,aa in original)
    schemes={'imgt':original}
    required={policy.numbering_scheme,policy.cdr_numbering_scheme}-{'imgt'}
    if required:
        from anarcii.output_data_processing.schemes import convert_number_scheme
        payload={'query':{'numbering':imgt_positions,'chain_type':chain_type,'score':score,'query_start':0,'query_end':len(original)-1,'error':None,'scheme':'imgt'}}
        for scheme in sorted(required):
            try:
                item=convert_number_scheme(deepcopy(payload),scheme)['query']
                real=[(p,aa) for p,aa in item['numbering'] if aa!='-']
            except Exception as exc:
                raise ScaffoldError(f'Could not convert IMGT to {scheme}: {exc}') from exc
            if ''.join(aa for _,aa in real)!=original_sequence:
                raise ScaffoldError(f'{scheme}: numbering conversion changed the variable-domain sequence')
            labels=[(int(p[0]),str(p[1]).strip()) for p,_ in real]
            if len(set(labels))!=len(labels):
                raise ScaffoldError(f'{scheme}: numbering conversion returned duplicate labels')
            schemes[scheme]=real
    rows=[]
    for (number,_),(cdr,_) in zip(schemes[policy.numbering_scheme],schemes[policy.cdr_numbering_scheme]):
        pos,ins=int(number[0]),str(number[1]).strip()
        cpos,cins=int(cdr[0]),str(cdr[1]).strip()
        rows.append({'numbering_scheme':policy.numbering_scheme,'numbering_position':pos,'numbering_insertion':ins,'numbering_label':f'{pos}{ins}','cdr_definition':policy.cdr_definition,'cdr_numbering_scheme':policy.cdr_numbering_scheme,'cdr_position':cpos,'cdr_insertion':cins,'cdr_label':f'{cpos}{cins}','region':region_for_label(cpos,cins,chain_type,policy.cdr_definition)})
    return rows

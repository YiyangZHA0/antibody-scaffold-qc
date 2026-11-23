"""Acquire an actual matching structure; never synthesize coordinates from sequence."""
from __future__ import annotations
from hashlib import sha256
from itertools import product
from pathlib import Path
import re
import requests
from .errors import ScaffoldError
from .structure import read_structure

def download_pdb(pdb_id, cache):
    if not re.fullmatch(r'[0-9][A-Za-z0-9]{3}',pdb_id):
        raise ScaffoldError(f'Invalid PDB ID: {pdb_id}')
    cache = Path(cache)
    cache.mkdir(parents=True,exist_ok=True)
    path = cache / f'{pdb_id.lower()}.cif'
    url = f'https://files.rcsb.org/download/{pdb_id.upper()}.cif'
    if not path.exists():
        response = requests.get(url,timeout=60)
        response.raise_for_status()
        if not response.content.lstrip().startswith(b'data_'):
            raise ScaffoldError(f'Non-mmCIF response from {url}')
        path.write_bytes(response.content)
    return path, {'structure_url':url,'structure_sha256':sha256(path.read_bytes()).hexdigest(),'resolution_method':'PDB identifier followed by exact variable-domain matching'}

def search_rcsb(domains, limit=25):
    """Explicit opt-in sequence transmission to RCSB's public search API."""
    common = None
    for domain in domains:
        query = {'query':{'type':'terminal','service':'sequence','parameters':{'evalue_cutoff':0.1,'identity_cutoff':1.0,'sequence_type':'protein','value':domain.sequence}},'return_type':'polymer_entity','request_options':{'paginate':{'start':0,'rows':limit},'results_content_type':['experimental']}}
        response = requests.post('https://search.rcsb.org/rcsbsearch/v2/query',json=query,timeout=60)
        if response.status_code == 204:
            return []
        response.raise_for_status()
        ids = {r['identifier'].split('_')[0].upper() for r in response.json().get('result_set',[])}
        common = ids if common is None else common & ids
    return sorted(common or [])

def select_structure_domains(source, domains, chain_map=None):
    candidates = []
    for domain in domains:
        role = 'heavy' if domain.chain_type=='H' else 'light'
        explicit = (chain_map or {}).get(role)
        possible = [key for key,chain in source.chains.items() if chain.sequence.count(domain.sequence)==1 and (explicit is None or key==explicit)]
        if not possible:
            raise ScaffoldError(f'No exact structural sequence match for {role} {domain.domain_id}; a similar sequence is insufficient')
        candidates.append(possible)
    assignments = list(product(*candidates))
    # For Fab, two domains originated in different input chains and must retain that distinction.
    if len(domains)==2 and domains[0].record_id != domains[1].record_id:
        assignments = [a for a in assignments if a[0]!=a[1]]
    if len(assignments)!=1:
        raise ScaffoldError(f'Ambiguous structural chains: {assignments}. Specify --heavy-chain/--light-chain or --vhh-chain.')
    chosen = assignments[0]
    return [(domain,key,'H' if domain.chain_type=='H' else 'L') for domain,key in zip(domains,chosen)]

def resolve_fasta_structure(annotation, records, cache, structure_path=None, pdb_id=None, search=False, chain_map=None, model_index=1, namespace='auth',input_path=None,offline=False):
    selected=list((chain_map or {}).values()) or None
    if not structure_path and not pdb_id and not search and input_path:
        # Local companion discovery is deterministic and requires no network.
        for suffix in ('.cif','.mmcif','.pdb'):
            companion=Path(input_path).with_suffix(suffix)
            if companion.is_file():
                source=read_structure(companion,model_index=model_index,selected_chains=selected,namespace=namespace)
                selections=select_structure_domains(source,annotation['domains'],chain_map)
                return source,selections,{'resolution_method':'same-basename local companion, exact variable-domain sequence match'}
    if structure_path:
        source=read_structure(structure_path,model_index=model_index,selected_chains=selected,namespace=namespace)
        selections=select_structure_domains(source,annotation['domains'],chain_map)
        return source,selections,{'resolution_method':'explicit structure, exact variable-domain sequence match'}
    if offline:
        raise ScaffoldError('Offline FASTA conversion needs coordinates from --structure or a same-stem local PDB/CIF; no remote request was made.')
    inferred = {m.group(1).upper() for rid in records if (m:=re.match(r'^([0-9][A-Za-z0-9]{3})[_|:]',rid))}
    if pdb_id:
        candidate_ids=[pdb_id.upper()]
    elif len(inferred)==1:
        candidate_ids=sorted(inferred)
    elif search:
        candidate_ids=search_rcsb(annotation['domains'])
    else:
        raise ScaffoldError('FASTA annotation is complete, but structural scaffold needs coordinates. Supply --structure, --pdb-id, RCSB FASTA headers, or --search-rcsb (sends selected sequences to RCSB). No coordinates were invented.')
    errors=[]
    for entry in candidate_ids:
        try:
            path,provenance=download_pdb(entry,cache)
            source=read_structure(path,model_index=model_index,selected_chains=selected,namespace=namespace)
            selections=select_structure_domains(source,annotation['domains'],chain_map)
            provenance['candidate_pdb_ids']=candidate_ids
            return source,selections,provenance
        except (ScaffoldError,requests.RequestException) as exc:
            errors.append(f'{entry}: {exc}')
    raise ScaffoldError('No unambiguous exact matching experimental structure found. '+ '; '.join(errors))

"""Refresh independent CDR-boundary fixtures from primary sources.

This maintenance command uses the network. Tests only read its saved JSON.
It parses literal values/HTML cells; it never executes downloaded Python.
"""
from datetime import datetime, timezone
from hashlib import sha256
from html import unescape
from pathlib import Path
import ast
import json
import re
import requests

ROOT = Path(__file__).resolve().parents[1]
URLS = {
    'anarci': 'https://raw.githubusercontent.com/oxpig/ANARCI/master/Example_scripts_and_sequences/ImmunoPDB.py',
    'ucl': 'http://www.bioinf.org.uk/abs/info.html',
}


def main():
    texts, sources = {}, {}
    for key, url in URLS.items():
        r = requests.get(url, timeout=60)
        r.raise_for_status()
        texts[key] = r.text
        sources[key] = {'url': url, 'retrieved_utc': datetime.now(timezone.utc).isoformat(),
                        'sha256': sha256(r.content).hexdigest()}
    node = next(n for n in ast.parse(texts['anarci']).body if isinstance(n, ast.Assign)
                and any(isinstance(t, ast.Name) and t.id == 'regions_in_chothia' for t in n.targets))
    reference = ast.literal_eval(node.value)
    boundaries = {'chothia': {loop: values['chothia'] for loop, values in reference.items()},
                  'maccallum': {loop: values['contact'] for loop, values in reference.items()}}
    table = next(t for t in re.findall(r'<table\b.*?</table>', texts['ucl'], flags=re.S | re.I) if 'AbM' in t)
    ucl = {'abm': {}, 'chothia-consensus': {}, 'maccallum': {}}
    for row in re.findall(r'<tr\b.*?</tr>', table, flags=re.S | re.I):
        cells = [unescape(re.sub(r'<[^>]+>', '', c)).strip() for c in re.findall(r'<td\b[^>]*>(.*?)</td>', row, flags=re.S | re.I)]
        if not cells:
            continue
        if cells[0] == 'H1' and 'Chothia/Martin' not in cells[1]:
            continue  # Use the explicit Chothia/Martin row, not Kabat H35B.
        for scheme, col in [('abm', 2), ('chothia-consensus', 3), ('maccallum', 4)]:
            match = re.fullmatch(r'[HL](\d+)-[HL](\d+)', cells[col])
            if not match:
                raise ValueError(f'Unexpected source table cell: {cells[col]}')
            ucl[scheme][cells[0]] = [int(x) for x in match.groups()]
    if {k: list(v) for k, v in boundaries['maccallum'].items()} != ucl['maccallum']:
        raise ValueError('Independent Oxford/UCL Contact definitions disagree')
    boundaries.update(ucl)
    if any(len(v) != 6 for v in boundaries.values()):
        raise ValueError('Incomplete six-loop reference')
    path = ROOT / 'tests/reference/cdr_boundaries.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({'sources': sources, 'numbering_basis': 'chothia', 'boundaries': boundaries}, indent=2) + '\n', encoding='utf-8')
    print(path)


if __name__ == '__main__':
    main()

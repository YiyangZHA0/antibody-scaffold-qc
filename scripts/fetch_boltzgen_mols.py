"""Fetch and verify BoltzGen's official molecular dictionary (no design weights)."""
from pathlib import Path
from urllib.request import Request, urlopen
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]

def main():
    cache = ROOT / '.cache'
    cache.mkdir(exist_ok=True)
    url = 'https://huggingface.co/datasets/boltzgen/inference-data/resolve/main/mols.zip'
    with urlopen('https://huggingface.co/api/datasets/boltzgen/inference-data/tree/main', timeout=60) as response:
        info = next(x for x in json.load(response) if x['path'] == 'mols.zip')
    path = cache / 'mols.zip'
    if not path.exists():
        temporary = path.with_suffix('.zip.part')
        with urlopen(Request(url, headers={'User-Agent': 'scaffold/0.7'}), timeout=120) as response, temporary.open('wb') as out:
            total = 0
            while block := response.read(1024 * 1024):
                out.write(block)
                total += len(block)
                if total % (64 * 1024 * 1024) == 0:
                    print(f'Downloaded {total // (1024*1024)} MiB', flush=True)
        temporary.replace(path)
    with path.open('rb') as handle:
        digest = hashlib.file_digest(handle, 'sha256').hexdigest()
    expected = info.get('lfs', {}).get('oid')
    if path.stat().st_size != info['size'] or (expected and digest != expected):
        raise ValueError('Official mols.zip size or SHA256 mismatch')
    record = {'url': url, 'bytes': path.stat().st_size, 'sha256': digest, 'lfs_sha256': expected}
    (cache / 'mols.source.json').write_text(json.dumps(record, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(record), flush=True)

if __name__ == '__main__':
    main()

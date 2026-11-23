"""Download unmodified public RCSB PDB, mmCIF and FASTA regression inputs."""
from __future__ import annotations
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import argparse
import time
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / 'tests/data'
CASES = {
    "5Y9K": "belimumab Fab",
    "7EOW": "caplacizumab VHH in vWF complex",
    "6CR1": "adalimumab EFab; test VH/VL only",
    "1X9Q": "4m5.3 anti-fluorescein scFv; VL-linker-VH",
    "7XL0": "vobarilizumab VHH; multiple copies",
    "6WGB": "dupilumab Fab; incomplete structure negative control",
    "5J13": "tezepelumab Fab with antigen",
}


def fetch(job):
    pdb_id, suffix, url = job
    path = DATA_DIR / f"{pdb_id.lower()}.{suffix}"
    path.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(3):
        try:
            with urlopen(Request(url, headers={"User-Agent": "scaffold/0.7 research-validation"}), timeout=60) as response:
                data = response.read()
            if suffix == "fasta" and not data.lstrip().startswith(b">"):
                raise ValueError("RCSB response is not FASTA")
            if suffix == "pdb" and b"ATOM  " not in data:
                raise ValueError("RCSB response has no PDB atoms")
            if suffix == "cif" and not data.lstrip().startswith(b"data_"):
                raise ValueError("RCSB response is not mmCIF")
            path.write_bytes(data)
            return {"pdb_id": pdb_id, "description": CASES[pdb_id], "file": path.relative_to(ROOT).as_posix(), "url": url, "sha256": sha256(data).hexdigest(), "bytes": len(data), "retrieved_utc": datetime.now(timezone.utc).isoformat()}
        except Exception:
            if attempt == 2:
                raise
            time.sleep(2 ** attempt)

def main():
    global DATA_DIR
    p=argparse.ArgumentParser()
    p.add_argument('--ids',nargs='+',choices=list(CASES),default=list(CASES))
    p.add_argument('--output-dir',type=Path,default=DATA_DIR)
    args=p.parse_args()
    DATA_DIR=args.output_dir.resolve()
    DATA_DIR.relative_to(ROOT)
    jobs = []
    for pdb_id in args.ids:
        for suffix in ("pdb", "cif"):
            jobs.append((pdb_id, suffix, f"https://files.rcsb.org/download/{pdb_id}.{suffix}"))
        jobs.append((pdb_id, "fasta", f"https://www.rcsb.org/fasta/entry/{pdb_id}/display"))
    with ThreadPoolExecutor(max_workers=3) as pool:
        records = list(pool.map(fetch, jobs))
    manifest=DATA_DIR / 'sources.json'
    old=json.loads(manifest.read_text()) if manifest.exists() else []
    records=[r for r in old if r['pdb_id'] not in args.ids]+records
    manifest.write_text(json.dumps(records, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(records, indent=2))

if __name__ == "__main__":
    main()

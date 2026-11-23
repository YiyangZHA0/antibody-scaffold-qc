"""Build the bundled endpoint panel from immutable public downloads (no query sequences).

Run with the project's Python. A pre-downloaded file can be supplied to each
option for a fully offline rebuild; both inputs have pinned SHA256 digests.
"""
from __future__ import annotations

import argparse
import ast
from collections import Counter, defaultdict
import gzip
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
import zipfile

from Bio.Seq import Seq
import requests


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = "https://github.com/repseqio/library-imgt/releases/download/v8/imgt.202312-3.sv8.json.gz"
ARCHIVE_HASH = "694b6beb7bb7d3a0d78f13e7b674c67ab1c41b043b2df722938507baa9c96b8b"
WHEEL = ("https://files.pythonhosted.org/packages/1e/50/1345543216a281cab84030a73ee4918b87e7249679926392d8fa44cb9420/"
         "anarci-2026.2.13.2-py3-none-any.whl")
WHEEL_HASH = "8f88d49857abc2c5427c6cc903ebdac7b35e6636635eaa9d60bd8f0d4e7d14b7"
SPECIES = {9606: ("Homo sapiens", "human"), 30538: ("Vicugna pacos", "alpaca")}
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


def pinned(path, url, expected):
    if path:
        data = path.read_bytes()
    else:
        response = requests.get(url, timeout=45)
        response.raise_for_status()
        data = response.content
    if sha256(data).hexdigest() != expected:
        raise ValueError(f"SHA256 mismatch: {path or url}")
    return data


def read_gapped_v(wheel_bytes):
    # Do not install or execute downloaded Python code.
    with zipfile.ZipFile(BytesIO(wheel_bytes)) as wheel:
        tree = ast.parse(wheel.read("anarci/germlines.py"))
        license_text = wheel.read("anarci-2026.2.13.2.dist-info/licenses/LICENCE").decode("utf-8")
    for node in tree.body:
        if isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name) and node.targets[0].id == "all_germlines":
            return ast.literal_eval(node.value)["V"], license_text
    raise ValueError("V alignment table missing")


def endpoint(gene, nt, role, species_key, grouped_v):
    anchor = gene["anchorPoints"]
    if gene["geneType"] == "V":
        # Only use original gapped coordinates independently confirmed by the
        # functional genomic sequence from the mature start to the CDR3 start.
        start, end = anchor.get("FR1Begin"), anchor.get("CDR3Begin")
        if start != 0 or end is None or end <= 0 or end % 3:
            return None
        # RepSeq's CDR3Begin starts at Cys104; include its codon for IMGT FR3.
        germline = str(Seq(nt[start:end + 3]).translate())
        aligned = grouped_v.get(role, {}).get(species_key, {}).get(gene["name"], "")[:104].replace("-", ".")
        if len(aligned) != 104 or aligned.replace(".", "") != germline or set(germline) - AMINO_ACIDS:
            return None
        return aligned
    # RepSeq's FR4Begin follows the W/F codon (it includes that codon
    # in its CDR3). IMGT FR4 includes W/F118. Do not translate the dangling
    # J splice nucleotide(s) by guessing bases from the constant domain.
    begin = anchor.get("FR4Begin")
    if begin is None or begin < 3 or anchor.get("FR4End") != len(nt):
        return None
    begin -= 3
    codon_start = int(gene["meta"]["comments"].split("|")[7])
    if begin % 3 != codon_start - 1:
        return None
    triplets = (len(nt) - begin) // 3
    aa = str(Seq(nt[begin:begin + 3 * triplets]).translate())
    if not aa.startswith(("W", "F")) or not 10 <= len(aa) <= 12 or set(aa) - AMINO_ACIDS:
        return None
    return aa


def build(archive, wheel_bytes):
    grouped_v, license_text = read_gapped_v(wheel_bytes)
    groups = defaultdict(list)
    source_counts = Counter()
    for item in json.loads(gzip.decompress(archive)):
        if item["taxonId"] not in SPECIES:
            continue
        species, species_key = SPECIES[item["taxonId"]]
        fragments = {f["uri"]: f for f in item["sequenceFragments"]}
        for gene in item["genes"]:
            segment = gene["geneType"]
            chains = gene["chains"]
            if segment not in {"V", "J"} or len(chains) != 1 or chains[0] not in {"IGH", "IGK", "IGL"}:
                continue
            role = chains[0][-1]
            if species_key == "alpaca" and role != "H":
                continue
            header = gene.get("meta", {}).get("comments", "").split("|")
            if len(header) < 14 or not gene["isFunctional"] or header[3] != "F" or "partial" in "|".join(header).lower():
                continue
            key = (species, role, segment)
            source_counts[key] += 1
            fragment = fragments.get(gene["baseSequence"])
            if fragment is None or fragment["range"]["from"] != 0:
                continue
            aligned = endpoint(gene, fragment["sequence"], role, species_key, grouped_v)
            if aligned:
                groups[(*key, aligned)].append((gene["name"], header[0]))
    entries, coverage = [], []
    for (species, role, segment, aligned), members in sorted(groups.items()):
        ordered = sorted(set(members))
        entries.append(dict(segment=segment, role=role, species=species, aligned=aligned,
                            genes=[g for g, _ in ordered], accessions=[a for _, a in ordered]))
    for (species, role, segment), total in sorted(source_counts.items()):
        entries_here = [e for e in entries if (e["species"], e["role"], e["segment"]) == (species, role, segment)]
        coverage.append(dict(species=species, role=role, segment=segment, source_alleles=total,
                             retained_alleles=sum(len(e["genes"]) for e in entries_here), unique_proteins=len(entries_here)))
    required = {(sp, r, g) for sp, r in [("Homo sapiens",c) for c in "HKL"]+[("Vicugna pacos","H")] for g in "VJ"}
    if required != {(e["species"], e["role"], e["segment"]) for e in entries}:
        raise ValueError("Missing required species, chain or segment")
    return dict(version="IMGT-202312-3_oxford-2026.2.13.2_v1",
                attribution="Protein features adapted from IMGT®, the international ImMunoGeneTics information system®, https://www.imgt.org, Montpellier, France (founder: Marie-Paule Lefranc). DOI: 10.1093/nar/gkn838.",
                sources=[
                    dict(url=ARCHIVE, sha256=ARCHIVE_HASH, release="IMGT 202312-3; RepSeq.IO v8",
                         archive_notes="https://github.com/repseqio/library-imgt",
                         license="CC BY 4.0 for IMGT data under the publisher's current terms; see the release notes for historical mirror terms.",
                         license_url="https://creativecommons.org/licenses/by/4.0/",
                         data_terms="https://www.imgt.org/about/termsofuse.php"),
                    dict(url=WHEEL, sha256=WHEEL_HASH, release="anarci 2026.2.13.2; gapped V coordinate corroboration only",
                         upstream="https://github.com/oxpig/ANARCI", license="BSD-3-Clause (ANARCI software); IMGT-derived sequence data",
                         license_notice=license_text)],
                scope="Functional complete source records. V includes only alleles with agreeing published gapped V alignment and translated FR1-through-FR3. J includes complete codons from W/F118; partial splice codons excluded. Coverage is a dated subset, not species identification.",
                coverage=coverage, entries=entries)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--imgt-archive", type=Path)
    parser.add_argument("--anarci-wheel", type=Path)
    parser.add_argument("--out", type=Path, default=ROOT/"src/scaffold/static/germline_endpoints.json")
    args = parser.parse_args()
    result = build(pinned(args.imgt_archive, ARCHIVE, ARCHIVE_HASH), pinned(args.anarci_wheel, WHEEL, WHEEL_HASH))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(result, ensure_ascii=False, separators=(",", ":")) + "\n"
    args.out.write_bytes(payload.encode("utf-8"))
    for row in result["coverage"]:
        print(row)


if __name__ == "__main__":
    main()

"""Post-extraction evidence, using the public 7XL0 FASTA and production writers."""
import json
from pathlib import Path
import re

import pytest

from scaffold.numbering import Numberer, classify, read_fasta
from scaffold.output import write_annotation
from scaffold.postcheck_sequence import postcheck_sequences
from scaffold.schemes import resolve_policy
from scaffold.sequence import write_sequence_scaffold

pytestmark = pytest.mark.integration
DATA = Path(__file__).parent / "data" / "7xl0.fasta"


@pytest.fixture(scope="module")
def converted(tmp_path_factory):
    records = read_fasta(DATA)
    numberer = Numberer()
    root = tmp_path_factory.mktemp("postcheck_numbered_7xl0")
    for preset in ("imgt", "maccallum", "aho"):
        numberer.policy = resolve_policy(preset)
        annotation = classify(records, numberer, "vhh")
        path = root / preset
        metadata = write_annotation(path, annotation, records, DATA)
        write_sequence_scaffold(path, metadata, annotation, records, ["CDR1", "CDR2", "CDR3"], ["heavy"])
    return root


def read_case(path):
    metadata = json.loads((path / "metadata.json").read_text(encoding="utf-8"))
    rows = [json.loads(line) for line in (path / "residue_map.jsonl").read_text(encoding="utf-8").splitlines()]
    return metadata, rows, path / "scaffold.yaml"


def cut_fr4_header(path, destination, remove_indices):
    lines = path.read_text(encoding="utf-8").splitlines()
    index = next(i for i, line in enumerate(lines) if line.startswith("# H-FR4:"))
    line = lines[index].replace("sequence=WGQGTLVTVSS;", "sequence=WGQGTLVTVS;")
    assert line != lines[index]
    if remove_indices:
        line = line.replace(",127,128;", ",127;")
        line = re.sub(r"(res_index=\d+\.\.)(\d+);", lambda m: f"{m[1]}{int(m[2])-1};", line)
    lines[index] = line
    destination.write_text("\n".join(lines), encoding="utf-8")


def test_tampered_fr4_in_header(converted, tmp_path):
    metadata, rows, yaml_path = read_case(converted / "imgt")
    altered = tmp_path / "scaffold.yaml"
    cut_fr4_header(yaml_path, altered, remove_indices=False)
    result = postcheck_sequences(metadata, rows, altered)
    problems = [f for f in result["findings"] if f["status"] == "observed"]
    assert [(f["chain"], f["region"], f["code"]) for f in problems] == [("H", "FR4", "framework_header_mismatch")]
    difference = problems[0]["evidence"]["differences"]["sequence"]
    assert difference == {"header": "WGQGTLVTVS", "map": "WGQGTLVTVSS"}


def test_fr4_row_and_header_both_truncated(converted, tmp_path):
    metadata, rows, yaml_path = read_case(converted / "imgt")
    last = rows.pop()
    assert (last["region"], last["aa"], last["imgt_label"]) == ("FR4", "S", "128")
    altered = tmp_path / "scaffold.yaml"
    cut_fr4_header(yaml_path, altered, remove_indices=True)
    result = postcheck_sequences(metadata, rows, altered)
    problems = [f for f in result["findings"] if f["status"] == "observed"]
    assert [f["code"] for f in problems] == ["framework_annotation_mismatch"]
    assert problems[0]["evidence"]["annotated_sequence"] == "WGQGTLVTVSS"
    assert problems[0]["evidence"]["mapped_sequence"] == "WGQGTLVTVS"


@pytest.mark.parametrize("preset,definition", [("imgt", "imgt"), ("maccallum", "maccallum"), ("aho", "imgt")])
def test_real_reference_uses_selected_definition(converted, preset, definition):
    metadata, rows, path = read_case(converted / preset)
    report = postcheck_sequences(metadata, rows, path, read_fasta(DATA))
    regions = {fr["region"]: fr for fr in report["frameworks"]}
    assert set(regions) == {"FR1", "FR2", "FR3", "FR4"}
    assert all(fr["cdr_definition"] == definition for fr in regions.values())
    # The Contact definition places extra residues in FR1; it should not use the IMGT length rule.
    assert (regions["FR1"]["length"] > 26) == (definition == "maccallum")
    assert all(f["status"] != "observed" and f["code"] != "atypical_imgt_length" for f in report["findings"])


def test_original_record_disagrees_with_numbered_residue(converted):
    metadata, rows, path = read_case(converted / "imgt")
    records = read_fasta(DATA)
    record_id = metadata["chains"][0]["record_id"]
    seq = records[record_id]
    index = rows[-1]["input_seq_index"] - 1
    assert seq[index] == "S"
    records[record_id] = seq[:index] + "A" + seq[index + 1:]
    report = postcheck_sequences(metadata, rows, path, records)
    problems = [f for f in report["findings"] if f["code"] == "source_sequence_mismatch"]
    assert len(problems) == 1
    assert problems[0]["evidence"] == {"basis": "original_file", "field": "input_seq_index",
                                      "items": [{"position": index + 1, "mapped_aa": "S", "source_aa": "A"}]}

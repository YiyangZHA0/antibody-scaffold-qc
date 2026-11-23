"""Check the FRs recorded by a completed conversion, in its original CDR definition."""
from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path
import re

from .audit import AA, TYPICAL, region_for
from .errors import ScaffoldError
from .output import ranges
from .schemes import BASE_NUMBERING, region_for_label


_HEADER = re.compile(r"^#\s+(\S+)-(FR[1-4]):\s*(.*)$")
_ANCHORS = {"23": "C", "41": "W", "89": "ILVFM", "104": "C", "118": "FW",
            "119": "G", "121": "G"}


def _headers(path):
    result = defaultdict(list)
    for line in Path(path).read_text(encoding="utf-8-sig").splitlines():
        match = _HEADER.match(line)
        if match:
            fields = {}
            for item in match[3].split(";"):
                key, sep, value = item.strip().partition("=")
                if sep:
                    fields[key] = value
            result[match[1], match[2]].append(fields)
    return result


def _integer(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _label(row, prefix):
    number = _integer(row.get(f"{prefix}_position"))
    if number is None or number <= 0:
        return None
    return f"{number}{row.get(f'{prefix}_insertion') or ''}"


def _included(row):
    return row.get("included") in (True, "True", "true", 1, "1")


def _sequence(rows):
    return "".join(str(row.get("aa") or "") for row in rows)


def _index_errors(rows, field, sequence):
    """Check against the corresponding full record, not its numbered substring."""
    problems = []
    for row in rows:
        position = _integer(row.get(field))
        letter = sequence[position - 1] if position and 1 <= position <= len(sequence) else None
        if letter != row.get("aa") or letter is None:
            problems.append({"position": position, "mapped_aa": row.get("aa"), "source_aa": letter})
    return problems


def _check_labels(rows, chain_type, definition, scheme, record):
    invalid = []
    misplaced = []
    seen = defaultdict(set)
    for row in rows:
        position = _integer(row.get("input_seq_index"))
        for prefix in ("numbering", "imgt", "cdr"):
            label = row.get(f"{prefix}_label")
            if not label or label != _label(row, prefix) or label in seen[prefix]:
                invalid.append({"input_position": position, "namespace": prefix, "label": label})
            seen[prefix].add(label)
        if (row.get("numbering_scheme") != scheme or row.get("cdr_definition") != definition
                or row.get("cdr_numbering_scheme") != BASE_NUMBERING.get(definition)):
            invalid.append({"input_position": position, "row_definition": row.get("cdr_definition"),
                            "row_scheme": row.get("numbering_scheme"),
                            "row_boundary_basis": row.get("cdr_numbering_scheme")})
        try:
            number = int(row["cdr_position"])
            expected = (region_for(number) if definition == "imgt"
                        else region_for_label(number, row.get("cdr_insertion") or "", chain_type, definition))
        except (KeyError, TypeError, ValueError, ScaffoldError):
            expected = "unassigned"
        if expected != row.get("region"):
            misplaced.append({"input_position": position, "label": row.get("cdr_label"),
                              "mapped_region": row.get("region"), "expected_region": expected})
    if invalid:
        record(None, "numbering_mapping_inconsistent", "observed",
               "A label, numbering basis or duplicate conflicts with the recorded convention.", items=invalid)
    if misplaced:
        record(None, "region_definition_mismatch", "observed",
               "The mapped region differs from the conversion's declared boundary rule.", items=misplaced)


def _terminal_checks(domain, source_sequence, record):
    rows = domain.get("rows") or []
    if not rows or not source_sequence:
        return
    first, last = rows[0], rows[-1]
    start, end = _integer(first.get("input_seq_index")), _integer(last.get("input_seq_index"))
    first_imgt = _integer(first.get("imgt_position"))
    last_imgt = _integer(last.get("imgt_position"))
    if start == 1 and first_imgt is not None and first_imgt > 1:
        record("FR1", "possible_n_terminal_truncation", "suspected",
               "The supplied record begins after IMGT position 1.", first_imgt=first_imgt)
    if end == len(source_sequence) and last_imgt is not None and last_imgt < 118:
        record("FR4", "possible_c_terminal_truncation", "suspected",
               "The sequence ends before the usual start of FR4 in IMGT numbering.", last_imgt=last_imgt)
    for region, flank in (("FR1", source_sequence[:start - 1] if start else ""),
                          ("FR4", source_sequence[end:] if end else "")):
        if 1 <= len(flank) <= 2:
            record(region, "short_unnumbered_terminal_flank", "suspected",
                   "One or two residues in the supplied record flank the numbered domain. "
                   "Review the domain boundary and construct annotations.", flank=flank)


def _check_framework(chain, region, local, expected, headers, metadata, record):
    sequence = _sequence(local)
    included = [row for row in local if _included(row)]
    numbers = [_integer(row.get("boltzgen_res_index")) for row in included]
    indices_valid = all(number is not None and number > 0 for number in numbers)
    scheme = metadata["numbering_scheme"]
    header_sequence = headers[0].get("sequence") if len(headers) == 1 else None
    if len(headers) != 1:
        record(region, "framework_header_count", "observed",
               "Expected exactly one exported FR header in scaffold.yaml.", count=len(headers))
    else:
        expected_fields = {
            "sequence": sequence,
            scheme: ",".join(str(row.get("numbering_label") or "") for row in local),
            "IMGT_reference": ",".join(str(row.get("imgt_label") or "") for row in local),
            "res_index": ranges(numbers) if indices_valid else None,
        }
        different = {key: {"header": headers[0].get(key), "map": value}
                     for key, value in expected_fields.items() if headers[0].get(key) != value}
        if different or not indices_valid:
            record(region, "framework_header_mismatch", "observed",
                   "Exported FR sequence, labels or included indices disagree with the residue map.",
                   differences=different, included_indices_valid=indices_valid)
    if expected is not None:
        fields = ("input_seq_index", "aa", "imgt_label", "numbering_label", "cdr_label")
        def identity(row):
            return tuple(str(row.get(key, "")) for key in fields)
        if list(map(identity, local)) != list(map(identity, expected)):
            record(region, "framework_annotation_mismatch", "observed",
                   "The residue map and original numbering annotation differ in residue choice, order or labels.",
                   annotated_sequence=_sequence(expected), mapped_sequence=sequence,
                   annotated_positions=[r.get("input_seq_index") for r in expected],
                   mapped_positions=[r.get("input_seq_index") for r in local])
    if not local:
        record(region, "framework_empty", "observed", "The exported map contains no residues for this FR.")
    nonstandard = [row.get("input_seq_index") for row in local if row.get("aa") not in AA]
    if nonstandard:
        record(region, "nonstandard_framework_residue", "observed", "The FR contains unknown or nonstandard letters.",
               input_positions=nonstandard)
    if local and metadata.get("cdr_definition") == "imgt":
        minimum, maximum = TYPICAL[region]
        if not minimum <= len(sequence) <= maximum:
            record(region, "atypical_imgt_length", "suspected",
                   "FR length is outside the typical IMGT range.", length=len(sequence), typical=[minimum, maximum])
    for row in local:
        label = row.get("imgt_label")
        allowed = _ANCHORS.get(label)
        if allowed and _label(row, "imgt") == label and row.get("aa") in AA and row["aa"] not in allowed:
            record(region, "noncanonical_anchor", "suspected",
                   "A numbered conserved position carries an atypical residue. Engineering or natural variation is possible.",
                   imgt_label=label, aa=row["aa"], typical=list(allowed), input_position=row.get("input_seq_index"))
    return {"chain": chain, "region": region, "sequence": sequence, "length": len(sequence),
            "exported_header_sequence": header_sequence, "included_sequence": _sequence(included),
            "included_length": len(included), "numbering_scheme": scheme, "cdr_definition": metadata["cdr_definition"]}


def postcheck_sequences(metadata, rows, yaml_path, source_records=None):
    """Return 'frameworks' and 'findings'; source_records are original input records by ID.

    Structural template residues (source_seq_index) and the numbered input record
    (input_seq_index) can have different offsets for structure-backed FASTA runs.
    """
    headers = _headers(yaml_path)
    findings, frameworks = [], []
    chains = metadata["chains"]
    inputs = metadata.get("full_input_sequences") or {}
    domains = {d["domain_id"]: d for d in metadata.get("domains") or []}
    grouped = defaultdict(list)
    for row in rows:
        grouped[row.get("output_chain_id")].append(row)

    for summary in chains:
        chain = summary["output_chain_id"]
        mapped = grouped[chain]
        def record(region, code, status, message, **evidence):
            findings.append({"chain": chain, "region": region, "code": code, "status": status,
                             "message": message, "evidence": evidence})
        domain = domains.get(summary.get("domain_id"), {})
        annotation_rows = domain.get("rows")
        if annotation_rows is None:
            record(None, "original_numbering_unavailable", "not_assessable",
                   "The conversion contains no original domain rows for comparison.")
        positions = [_integer(row.get("input_seq_index")) for row in mapped]
        if any(p is None or p <= 0 for p in positions) or any(
                a is not None and b is not None and b != a + 1 for a, b in zip(positions, positions[1:])):
            record(None, "nonconsecutive_source_positions", "observed",
                   "The exported domain rows skip, repeat or reverse input sequence positions.", positions=positions)
        _check_labels(mapped, summary.get("chain_type"), metadata["cdr_definition"], metadata["numbering_scheme"], record)
        record_id = summary.get("record_id")
        input_sequence = inputs.get(record_id)
        for field, reference, basis in (
                ("input_seq_index", input_sequence, "saved_input_record"),
                ("source_seq_index", summary.get("full_source_sequence"), "saved_source_chain"),
                ("input_seq_index", (source_records or {}).get(record_id), "original_file")):
            if reference is None:
                if basis == "saved_input_record" or (basis == "original_file" and source_records is not None):
                    record(None, "source_record_unavailable", "not_assessable",
                           "Input sequence is unavailable for this chain.", record_id=record_id, basis=basis)
                continue
            errors = _index_errors(mapped, field, reference)
            if errors:
                record(None, "source_sequence_mismatch", "observed",
                       "Mapped amino acids disagree with the source record at their stated indices.",
                       basis=basis, field=field, items=errors)
        _terminal_checks(domain, (source_records or {}).get(record_id, input_sequence), record)
        for region in TYPICAL:
            local = [row for row in mapped if row.get("region") == region]
            expected = (None if annotation_rows is None else
                        [row for row in annotation_rows if row.get("region") == region])
            frameworks.append(_check_framework(chain, region, local, expected, headers.get((chain, region), []),
                                               metadata, record))
    chain_ids = [c["output_chain_id"] for c in chains]
    unknown = (set(grouped) | {chain for chain, _ in headers}) - set(chain_ids)
    duplicates = [chain for chain, count in Counter(chain_ids).items() if count > 1]
    if unknown or duplicates:
        findings.append({"chain": None, "region": None, "code": "unexpected_output_chains", "status": "observed",
                         "message": "Chain IDs disagree between conversion artifacts.",
                         "evidence": {"unknown": list(map(str, unknown)), "duplicate_metadata_ids": duplicates}})
    return {"frameworks": frameworks, "findings": findings}

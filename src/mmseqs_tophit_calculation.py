#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import io
import sys
from pathlib import Path
from typing import Dict, List, Tuple

# Column indices for MMseqs m8 format (similar to BLAST output fmt 6).
# These positions are locked: phidra_run.py always emits these 12 columns
# first via --format-output, and any user-requested extra columns (e.g.
# qcov,tcov or qaln,taln) are appended after index 11, never inserted
QID, SID, PIDENT, LENGTH, MISMATCH, GAPOPEN, QSTART, QEND, SSTART, SEND, EVALUE, BITS = range(12)

BASE_HEADERS: List[str] = [
    "Query_ID",
    "Target_ID",
    "Sequence_Identity",
    "Alignment_Length",
    "Mismatches",
    "Gap_Openings",
    "Query_Start",
    "Query_End",
    "Target_Start",
    "Target_End",
    "Evalue",
    "Bit_Score",
]

Row = List[str]
Rec = Dict[str, str]


def _open_maybe_gz(path: Path) -> io.TextIOBase:
    """Open text file, accepting .gz/.bgz. Why: large outputs are often compressed."""
    p = str(path)
    if p.endswith((".gz", ".bgz")):
        return io.TextIOWrapper(gzip.open(p, mode="rb"), encoding="utf-8", newline="")
    return open(p, mode="r", encoding="utf-8", newline="")


def _to_rec(row: Row, extra_columns: List[str]) -> Rec:
    rec = {
        "Query_ID": row[QID],
        "Target_ID": row[SID],
        "Sequence_Identity": row[PIDENT],
        "Alignment_Length": row[LENGTH],
        "Mismatches": row[MISMATCH],
        "Gap_Openings": row[GAPOPEN],
        "Query_Start": row[QSTART],
        "Query_End": row[QEND],
        "Target_Start": row[SSTART],
        "Target_End": row[SEND],
        "Evalue": row[EVALUE],
        "Bit_Score": row[BITS],
    }

    for i, col in enumerate(extra_columns):
        rec[col] = row[12 + i]
    return rec


def _parse_numeric(row: Row) -> Tuple[float, float, float]:
    """Return (evalue, bitscore, pident) or raise ValueError. Why: centralize coercion."""
    return float(row[EVALUE]), float(row[BITS]), float(row[PIDENT])


# Tie breakers below use (evalue, bitscore, identity) to resolve exact ties
# on the primary sort key, so two hits with identical evalue or bitscore
# still get a deterministic winner instead of depending on row order.
def _better_by_evalue(new: Tuple[float, float, float], best: Tuple[float, float, float]) -> bool:
    # Lower evalue - higher bitscore - higher identity
    if new[0] != best[0]:
        return new[0] < best[0]
    if new[1] != best[1]:
        return new[1] > best[1]
    return new[2] > best[2]


def _better_by_bitscore(new: Tuple[float, float, float], best: Tuple[float, float, float]) -> bool:
    # Higher bitscore - lower evalue - higher identity
    if new[1] != best[1]:
        return new[1] > best[1]
    if new[0] != best[0]:
        return new[0] < best[0]
    return new[2] > best[2]


def compute_top_hits(input_path: Path, extra_columns: List[str]) -> Tuple[List[Rec], List[Rec]]:
    top_e: Dict[str, Tuple[Tuple[float, float, float], Rec]] = {}
    top_b: Dict[str, Tuple[Tuple[float, float, float], Rec]] = {}

    total = 0
    skipped = 0
    min_len = 12 + len(extra_columns)

    with _open_maybe_gz(input_path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            total += 1
            if not row or row[0].startswith("#"):
                continue
            if len(row) < min_len:
                skipped += 1
                continue
            try:
                evalue, bits, pident = _parse_numeric(row)
            except ValueError:
                skipped += 1
                continue

            rec = _to_rec(row, extra_columns)
            qid = row[QID]
            triple = (evalue, bits, pident)

            if qid not in top_e or _better_by_evalue(triple, top_e[qid][0]):
                top_e[qid] = (triple, rec)
            if qid not in top_b or _better_by_bitscore(triple, top_b[qid][0]):
                top_b[qid] = (triple, rec)

    # Deterministic order for outputs
    evalue_records = [pair[1] for q, pair in sorted(top_e.items(), key=lambda kv: kv[0])]
    bits_records = [pair[1] for q, pair in sorted(top_b.items(), key=lambda kv: kv[0])]

    # Minimal diagnostics to stderr (why: help users catch format issues quickly)
    print(
        f"Processed {total} lines; kept {len(top_e)} queries; skipped {skipped} malformed/invalid lines.",
        file=sys.stderr,
    )

    return evalue_records, bits_records


def write_tsv(path: Path, rows: List[Rec], headers: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=headers, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description="Compute per-query top hits by evalue and bitscore from MMseqs m8.")
    ap.add_argument("input_file", help="MMseqs m8 file (optionally .gz)")
    ap.add_argument("output_evalue_file", help="Output TSV for best-by-evalue")
    ap.add_argument("output_bitscore_file", help="Output TSV for best-by-bitscore")
    ap.add_argument(
        "--extra_columns", default="",
        help=(
            "Comma-separated extra mmseqs column names, in the same order passed "
            "to mmseqs --format-output (e.g. qcov,tcov). Must match phidra_run.py's "
            "--extra_columns for the same run."
        ),
    )
    args = ap.parse_args()

    inp = Path(args.input_file)
    if not inp.exists():
        sys.exit(f"Input file not found: {inp}")

    extra_columns = [c.strip() for c in args.extra_columns.split(",") if c.strip()]
    headers = BASE_HEADERS + extra_columns

    e_rows, b_rows = compute_top_hits(inp, extra_columns)
    write_tsv(Path(args.output_evalue_file), e_rows, headers)
    write_tsv(Path(args.output_bitscore_file), b_rows, headers)


if __name__ == "__main__":
    main()

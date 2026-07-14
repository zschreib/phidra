#!/usr/bin/env python3

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path
from shlex import quote as _q
from typing import Dict

from tqdm import tqdm

BIN = Path(__file__).parent / "src"
LOG_DIR: Path | None = None

# Locked base columns for mmseqs --format-output. Order and identity must
# never change: mmseqs_tophit_calculation.py reads these by fixed position.
# Any user-requested --extra_columns are appended after
MMSEQS_BASE_COLS = [
    "query", "target", "pident", "alnlen", "mismatch", "gapopen",
    "qstart", "qend", "tstart", "tend", "evalue", "bits",
]

# Extra columns that only get populated with a sufficiently high alignment
# mode (they depend on the backtrace, which lower modes skip for speed).
ALIGNMENT_DEPENDENT_COLS = {"qaln", "taln", "cigar"}

# Only works with ORF input at the moment
MMSEQS_SEARCH_TYPE = 1


def build_format_output(extra_columns: list[str]) -> str:
    """Build a --format-output string: locked base columns + user extras appended.

    Rejects re-listing a base column or duplicate extras
    """
    dupes_with_base = set(extra_columns) & set(MMSEQS_BASE_COLS)
    if dupes_with_base:
        sys.exit(
            f"--extra_columns cannot re-list base columns already included by "
            f"default: {sorted(dupes_with_base)}"
        )
    if len(extra_columns) != len(set(extra_columns)):
        sys.exit("--extra_columns contains duplicate entries.")
    return ",".join(MMSEQS_BASE_COLS + extra_columns)


def mmseqs_extra_flags(args: argparse.Namespace) -> str:
    """Build mmseqs easy-search flags: user-configurable ones plus the locked search-type.

    """
    parts = [f"--search-type {MMSEQS_SEARCH_TYPE}"]
    if args.alignment_mode is not None:
        parts.append(f"--alignment-mode {args.alignment_mode}")
    if args.sensitivity is not None:
        parts.append(f"-s {args.sensitivity}")
    return " ".join(parts)


def q(x: object) -> str:
    return _q(str(x))


def set_log_dir(d: Path) -> None:
    """Set global log directory and ensure it exists."""
    global LOG_DIR
    LOG_DIR = d
    d.mkdir(parents=True, exist_ok=True)


def log(msg: str) -> None:
    """Write to console and persistent log. Fallback to stdout if LOG_DIR unset."""
    tqdm.write(msg)
    if LOG_DIR is None:
        return
    log_file = LOG_DIR / "job_logger.txt"
    log_file.parent.mkdir(parents=True, exist_ok=True)
    with log_file.open("a", encoding="utf-8") as f:
        f.write(msg + "\n")


def log_job_details(args: argparse.Namespace) -> None:
    """Emit a single, clearly delimited header with all CLI inputs."""
    # User need a reproducible record of every invocation parameter at the very top
    sep = "=" * 28
    lines = [
        f"{sep} JOB Details {sep}",
        f"Running {args.input_fasta} with the following parameters:",
        f"--input_fasta      {args.input_fasta}",
        f"--subject_db       {args.subject_db}",
        f"--pfam_hmm_db      {args.pfam_hmm_db}",
        f"--ida_file         {args.ida_file}",
        f"--function         {args.function}",
        f"--output_dir       {args.output_dir}",
        f"--threads          {args.threads}",
        f"--evalue           {args.evalue}",
        f"--extra_columns    {args.extra_columns or '(none)'}",
        f"--alignment_mode   {args.alignment_mode if args.alignment_mode is not None else '(mmseqs default)'}",
        f"--search_type      {MMSEQS_SEARCH_TYPE} (fixed: protein input, not user-configurable)",
        f"--sensitivity      {args.sensitivity if args.sensitivity is not None else '(mmseqs default)'}",
        f"{sep*2}",
    ]
    for line in lines:
        log(line)


def run(cmd: str) -> str:
    """Run a shell command; capture stdout; log stderr on failure."""
    try:
        return subprocess.run(
            cmd, shell=True, check=True, capture_output=True, text=True
        ).stdout
    except subprocess.CalledProcessError as e:
        log(f"ERROR: {e.stderr.strip()}")
        sys.exit(1)


def safe_rmtree(p: Path) -> None:
    """Best-effort recursive delete without shell."""
    try:
        shutil.rmtree(p, ignore_errors=True)
    except Exception as e:
        log(f"WARNING: failed to remove {p}: {e}")


def make_dirs(base: Path, func: str) -> Dict[str, Path]:
    """Create and return all working directories."""
    p = base / func
    dirs = {
        "db": p / "databases" / "subject",
        "rdb": p / "databases" / "recursive",
        "mm_i": p / "mmseqs" / "initial",
        "mm_r": p / "mmseqs" / "recursive",
        "pf_i": p / "pfam" / "initial",
        "pf_r": p / "pfam" / "recursive",
        "final": p / "final_results",
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)
    dirs["val"] = dirs["final"] / "validated_ida_pfams"
    dirs["unv"] = dirs["final"] / "unvalidated_ida_pfams"
    dirs["val"].mkdir(exist_ok=True)
    dirs["unv"].mkdir(exist_ok=True)
    return dirs


def fasta_count(path: Path) -> int:
    """Count FASTA records by header lines."""
    if not path.exists():
        return 0
    return sum(1 for l in path.open(encoding="utf-8", errors="ignore") if l.startswith(">"))


def tsv_count(path: Path) -> int:
    """Count data rows in TSV (skip header if present)."""
    if not path.exists():
        return 0
    return max(0, sum(1 for _ in path.open(encoding="utf-8", errors="ignore")) - 1)


def preflight(args: argparse.Namespace) -> None:
    # Files
    for fp, desc in [
        (Path(args.input_fasta), "input FASTA"),
        (Path(args.subject_db), "subject FASTA for mmseqs createdb"),
        (Path(args.pfam_hmm_db), "Pfam HMM DB"),
        (Path(args.ida_file), "IDA file"),
    ]:
        if not fp.exists() or fp.stat().st_size == 0:
            sys.exit(f"{desc} not found or empty: {fp}")

    # Binaries
    if shutil.which("mmseqs") is None:
        sys.exit("Required executable not found in PATH: mmseqs")

    # Helper scripts
    for script in [
        BIN / "mmseqs_tophit_calculation.py",
        BIN / "generate_tophit_fasta.py",
        BIN / "pfam_scan.py",
        BIN / "validate_pfam.py",
        BIN / "remove_duplicates.py",
    ]:
        if not script.exists():
            sys.exit(f"Required helper script missing: {script}")

    # Threads
    if args.threads < 1:
        sys.exit("--threads must be >= 1")

    args.extra_columns_list = [c.strip() for c in args.extra_columns.split(",") if c.strip()]

    needs_aln = ALIGNMENT_DEPENDENT_COLS & set(args.extra_columns_list)
    if needs_aln and (args.alignment_mode is None or args.alignment_mode < 3):
        sys.exit(
            f"--extra_columns includes {sorted(needs_aln)}, which require "
            f"--alignment_mode 3 or higher (got {args.alignment_mode})."
        )

    build_format_output(args.extra_columns_list)


def process(args: argparse.Namespace) -> None:
    base = Path(args.output_dir)
    set_log_dir(base / args.function)
    dirs = make_dirs(base, args.function)

    py = q(sys.executable)

    fmt = build_format_output(args.extra_columns_list)
    extra_flags = mmseqs_extra_flags(args)
    extra_flags_part = f"{extra_flags} " if extra_flags else ""
    tophit_extra_arg = (
        f" --extra_columns {q(args.extra_columns)}" if args.extra_columns_list else ""
    )

    # 1) Initial workflow
    steps = [
        (
            "Initial DB build",
            f"mmseqs createdb {q(args.subject_db)} {q(dirs['db'] / 'db')}",
        ),
        (
            "Index",
            f"mmseqs createindex {q(dirs['db'] / 'db')} {q(dirs['db'] / 'tmp')}",
        ),
        (
            "Search",
            f"mmseqs easy-search -e {q(args.evalue)} --threads {args.threads} "
            f"{extra_flags_part}--format-output {q(fmt)} "
            f"{q(args.input_fasta)} {q(dirs['db'] / 'db')} {q(dirs['mm_i'] / 'res.m8')} {q(dirs['mm_i'] / 'tmp')}",
        ),
        (
            "TopHits",
            f"{py} {q(BIN / 'mmseqs_tophit_calculation.py')} "
            f"{q(dirs['mm_i'] / 'res.m8')} {q(dirs['mm_i'] / 'hits.tsv')} {q(dirs['mm_i'] / 'bits.tsv')}{tophit_extra_arg}",
        ),
        (
            "MakeFASTA",
            f"{py} {q(BIN / 'generate_tophit_fasta.py')} "
            f"{q(dirs['mm_i'] / 'hits.tsv')} {q(args.input_fasta)} {q(dirs['mm_i'] / 'hits.fa')}",
        ),
        (
            "PfamScan",
            f"{py} {q(BIN / 'pfam_scan.py')} -out {q(dirs['pf_i'] / 'pfam_coverage_report.tsv')} "
            f"-outfmt tsv -cpu {args.threads} {q(dirs['mm_i'] / 'hits.fa')} {q(args.pfam_hmm_db)}",
        ),
        (
            "Validate",
            f"{py} {q(BIN / 'validate_pfam.py')} "
            f"{q(dirs['pf_i'] / 'pfam_coverage_report.tsv')} {q(args.input_fasta)} {q(args.ida_file)} {q(dirs['pf_i'])}",
        ),
        (
            "FilterOut",
            f"{py} {q(BIN / 'remove_duplicates.py')} "
            f"{q(args.input_fasta)} "
            f"{q(dirs['mm_i'] / 'hits.fa')} "
            f"{q(dirs['rdb'] / 'filtered.fa')}",
        ),
    ]
    for name, cmd in steps:
        log(f"[{name}] {cmd}")
        run(cmd)

    # 2) Recursive workflow
    rem = dirs["rdb"] / "filtered.fa"
    if rem.exists() and rem.stat().st_size > 0:
        log("[Recursive] Starting recursive search")
        recurse_out = dirs["mm_r"] / "res.m8"
        recurse_tmp = dirs["mm_r"] / "tmp"
        cmd = (
            f"mmseqs easy-search -e {q(args.evalue)} --threads {args.threads} "
            f"{extra_flags_part}--format-output {q(fmt)} "
            f"{q(rem)} "
            f"{q(dirs['pf_i'] / 'validated_ida_report' / 'full_proteins.fa')} "
            f"{q(recurse_out)} {q(recurse_tmp)}"
        )
        run(cmd)

        if recurse_out.exists() and recurse_out.stat().st_size > 0:
            rec_steps = [
                (
                    "TopHits_R",
                    f"{py} {q(BIN / 'mmseqs_tophit_calculation.py')} "
                    f"{q(recurse_out)} {q(dirs['mm_r'] / 'hits.tsv')} {q(dirs['mm_r'] / 'bits.tsv')}{tophit_extra_arg}",
                ),
                (
                    "MakeFASTA_R",
                    f"{py} {q(BIN / 'generate_tophit_fasta.py')} "
                    f"{q(dirs['mm_r'] / 'hits.tsv')} {q(args.input_fasta)} {q(dirs['mm_r'] / 'hits.fa')}",
                ),
                (
                    "PfamScan_R",
                    f"{py} {q(BIN / 'pfam_scan.py')} -out {q(dirs['pf_r'] / 'pfam_coverage_report.tsv')} "
                    f"-outfmt tsv -cpu {args.threads} {q(dirs['mm_r'] / 'hits.fa')} {q(args.pfam_hmm_db)}",
                ),
                (
                    "Validate_R",
                    f"{py} {q(BIN / 'validate_pfam.py')} "
                    f"{q(dirs['pf_r'] / 'pfam_coverage_report.tsv')} {q(args.input_fasta)} {q(args.ida_file)} {q(dirs['pf_r'])}",
                ),
            ]
            for name, cmd in rec_steps:
                log(f"[{name}] {cmd}")
                run(cmd)
        else:
            log("Recursive search produced no hits; skipping recursive steps.")
            with open(recurse_out, "w", encoding="utf-8") as f:
                f.write("No significant hits.\n")

    # 3) Merge validated/unvalidated reports & sequences
    for status in ("validated", "unvalidated"):
        src_i = dirs["pf_i"] / f"{status}_ida_report"
        src_r = dirs["pf_r"] / f"{status}_ida_report"
        tgt = dirs["final"] / f"{status}_ida_pfams"

        # a) merged validation report
        rep_i = src_i / f"pfam_{status}_report.tsv"
        rep_r = src_r / f"pfam_{status}_report.tsv"
        out_rep = tgt / f"pfam_{status}_merged_report.tsv"
        with open(out_rep, "wb") as w:
            if rep_i.exists():
                w.write(rep_i.read_bytes())
            if rep_r.exists():
                lines = rep_r.read_text(encoding="utf-8").splitlines(True)[1:]
                w.writelines(l.encode() for l in lines)

        # b) merged domains.fa
        out_dom = tgt / "domains.fa"
        with open(out_dom, "wb") as w:
            for src in (src_i / "domains.fa", src_r / "domains.fa"):
                if src.exists():
                    w.write(src.read_bytes())

        # c) merged full_proteins.fa
        out_full = tgt / "full_proteins.fa"
        with open(out_full, "wb") as w:
            for src in (src_i / "full_proteins.fa", src_r / "full_proteins.fa"):
                if src.exists():
                    w.write(src.read_bytes())

    # 4) Merge pfam_coverage_report.tsv into final_results/
    cov_i = dirs["pf_i"] / "pfam_coverage_report.tsv"
    cov_r = dirs["pf_r"] / "pfam_coverage_report.tsv"
    final_cov = dirs["final"] / "pfam_coverage_report.tsv"
    with open(final_cov, "wb") as w:
        if cov_i.exists():
            w.write(cov_i.read_bytes())
        if cov_r.exists():
            lines = cov_r.read_text(encoding="utf-8").splitlines(True)[1:]
            w.writelines(l.encode() for l in lines)

    # 5) Build summary table
    initial_input = fasta_count(Path(args.input_fasta))
    initial_homology_hits = fasta_count(dirs["mm_i"] / "hits.fa")
    recursive_homology_hits = fasta_count(dirs["mm_r"] / "hits.fa")
    initial_validated = tsv_count(dirs["pf_i"] / "validated_ida_report" / "pfam_validated_report.tsv")
    initial_unvalidated = tsv_count(dirs["pf_i"] / "unvalidated_ida_report" / "pfam_unvalidated_report.tsv")
    recursive_validated = tsv_count(dirs["pf_r"] / "validated_ida_report" / "pfam_validated_report.tsv")
    recursive_unvalidated = tsv_count(dirs["pf_r"] / "unvalidated_ida_report" / "pfam_unvalidated_report.tsv")

    summary = {
        "initial_input": initial_input,
        "initial_homology_hits": initial_homology_hits,
        "recursive_homology_hits": recursive_homology_hits,
        "initial_validated": initial_validated,
        "initial_unvalidated": initial_unvalidated,
        "recursive_validated": recursive_validated,
        "recursive_unvalidated": recursive_unvalidated,
        "sequences_tossed": initial_input - (
            initial_validated
            + initial_unvalidated
            + recursive_validated
            + recursive_unvalidated
        ),
    }

    summary_file = dirs["final"] / "summary.tsv"
    with summary_file.open("w", encoding="utf-8") as f:
        f.write("metric\tcount\n")
        for m, c in summary.items():
            f.write(f"{m}\t{c}\n")

    log("Summary:")
    for m, c in summary.items():
        log(f"  {m.replace('_', ' ')}: {c}")

    # Cleanup (including databases directory)
    for d in (
        dirs["mm_i"] / "tmp",
        dirs["mm_r"] / "tmp",
        dirs["db"],
        dirs["rdb"],
        dirs["db"].parent,  # removes the "databases" directory
    ):
        safe_rmtree(d)

    log("Job complete!")


def main() -> None:
    p = argparse.ArgumentParser(description="Identifies homologous proteins and associated Pfam domains from input protein sequences, while comparing against InterPro Domain Architectures to analyze domain-level similarities and functional relationships.", 
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    add_help=False)
    
    help_args = p.add_argument_group("Help")
    help_args.add_argument(
        "-h", "--help", action="help", help="Show this help message and exit"
    )
    help_args.add_argument(
        "-v", "--version", action="version", version="PHIDRA v2.1.0",
        help="Show program version and exit"
    )
    
    required = p.add_argument_group("Required arguments")
    required.add_argument("-i", "--input_fasta", required=True, help="Query FASTA for mmseqs search")
    required.add_argument("-db", "--subject_db", required=True, help="Subject FASTA for mmseqs createdb")
    required.add_argument("-pfam", "--pfam_hmm_db", required=True, help="Pfam HMM format database path")
    required.add_argument("-ida", "--ida_file", required=True, help="IDA TSV file")
    required.add_argument("-f", "--function", required=True, help="User label for this run")
    required.add_argument("-o", "--output_dir", required=True, help="Base output directory")
    
    optional = p.add_argument_group("Optional arguments")
    optional.add_argument("-t", "--threads", type=int, default=1, help="Threads for tools supporting -cpu/--threads")
    optional.add_argument(
        "-e",
        "--evalue",
        default="1E-3",
        help="E-value threshold for mmseqs easy-search (e.g., 1E-3, 1e-5)",
    )
    optional.add_argument(
        "-x", "--extra_columns", default="", metavar="COL1,COL2,...",
        help=(
            "Comma-separated extra mmseqs --format-output columns appended after "
            "the fixed base columns (e.g. qcov,tcov or qaln,taln). See "
            "https://github.com/soedinglab/MMseqs2/wiki#custom-alignment-format-with-convertalis "
            "for valid column names."
        ),
    )
    optional.add_argument(
        "-a", "--alignment_mode", type=int, choices=[0, 1, 2, 3, 4], default=None,
        help=(
            "mmseqs --alignment-mode (default: mmseqs' own default for easy-search). "
            "Use 3 or higher for real (not estimated) sequence identity, or if "
            "--extra_columns requests alignment-dependent fields like qaln/taln/cigar."
        ),
    )
    optional.add_argument(
        "-s", "--sensitivity", type=float, default=None,
        help="mmseqs -s sensitivity (default: mmseqs' own default)",
    )
    args = p.parse_args()

    out = Path(args.output_dir) / args.function
    if out.exists():
        sys.exit(f"{out} exists.")

    set_log_dir(out)
    log_job_details(args)  # header printed before any processing
    preflight(args)
    process(args)


if __name__ == "__main__":
    main()

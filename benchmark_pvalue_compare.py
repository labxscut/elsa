#!/usr/bin/env python3
"""Benchmark helper for comparing two p-value modes.

This script is intentionally small and generic:
- it can time two shell commands that produce result files;
- it can compare the resulting p-values row-by-row;
- it reports absolute and relative precision gaps.

The default row key is auto-detected from the common identifier columns in the
result files, preferring X/Y/Z and falling back to X/Y when Z is absent.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import statistics
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


KEY_CANDIDATES = ("X", "Y", "Z")
P_CANDIDATES = ("P", "p", "P-value", "P_value", "pvalue")


@dataclass
class RunResult:
    seconds: float
    returncode: int


def run_shell_command(command: str) -> RunResult:
    start = time.perf_counter()
    completed = subprocess.run(command, shell=True, check=False)
    elapsed = time.perf_counter() - start
    if completed.returncode != 0:
        raise RuntimeError(f"Command failed with exit code {completed.returncode}: {command}")
    return RunResult(seconds=elapsed, returncode=completed.returncode)


def split_row(line: str) -> List[str]:
    return [part for part in re.split(r"\s+", line.strip()) if part]


def load_table(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    with path.open("r", encoding="utf-8") as handle:
        lines = [line for line in handle if line.strip()]

    if not lines:
        raise ValueError(f"Empty result file: {path}")

    headers = split_row(lines[0])
    rows: List[Dict[str, str]] = []
    for raw_line in lines[1:]:
        values = split_row(raw_line)
        if len(values) < len(headers):
            raise ValueError(
                f"Row has fewer fields than header in {path}: {raw_line.rstrip()}"
            )
        row = {headers[i]: values[i] for i in range(len(headers))}
        rows.append(row)
    return headers, rows


def detect_key_columns(headers_a: Sequence[str], headers_b: Sequence[str]) -> List[str]:
    common = [name for name in KEY_CANDIDATES if name in headers_a and name in headers_b]
    if not common:
        raise ValueError(
            "Could not detect shared key columns. Expected at least X and Y in both files."
        )
    return common


def detect_p_column(headers: Sequence[str]) -> str:
    for candidate in P_CANDIDATES:
        if candidate in headers:
            return candidate
    raise ValueError(f"Could not find a p-value column in headers: {headers}")


def row_key(row: Dict[str, str], key_columns: Sequence[str]) -> Tuple[str, ...]:
    return tuple(row[column] for column in key_columns)


def to_float(value: str) -> float:
    try:
        return float(value)
    except ValueError as exc:
        raise ValueError(f"Cannot parse float value: {value!r}") from exc


def compare_tables(
    theo_path: Path,
    perm_path: Path,
    tolerance: float,
) -> Dict[str, object]:
    theo_headers, theo_rows = load_table(theo_path)
    perm_headers, perm_rows = load_table(perm_path)

    key_columns = detect_key_columns(theo_headers, perm_headers)
    theo_p_column = detect_p_column(theo_headers)
    perm_p_column = detect_p_column(perm_headers)

    theo_map = {row_key(row, key_columns): row for row in theo_rows}
    perm_map = {row_key(row, key_columns): row for row in perm_rows}

    common_keys = sorted(set(theo_map) & set(perm_map))
    missing_in_perm = sorted(set(theo_map) - set(perm_map))
    missing_in_theo = sorted(set(perm_map) - set(theo_map))

    diffs: List[Dict[str, object]] = []
    abs_diffs: List[float] = []
    rel_diffs: List[float] = []
    within_tolerance = 0

    for key in common_keys:
        theo_p = to_float(theo_map[key][theo_p_column])
        perm_p = to_float(perm_map[key][perm_p_column])
        abs_diff = abs(theo_p - perm_p)
        rel_diff = abs_diff / max(abs(theo_p), abs(perm_p), 1e-12)
        abs_diffs.append(abs_diff)
        rel_diffs.append(rel_diff)
        if abs_diff <= tolerance:
            within_tolerance += 1
        diffs.append(
            {
                "key": key,
                "theo_p": theo_p,
                "perm_p": perm_p,
                "abs_diff": abs_diff,
                "rel_diff": rel_diff,
            }
        )

    def safe_mean(values: Sequence[float]) -> float:
        return statistics.fmean(values) if values else float("nan")

    summary: Dict[str, object] = {
        "key_columns": key_columns,
        "theo_rows": len(theo_rows),
        "perm_rows": len(perm_rows),
        "matched_rows": len(common_keys),
        "missing_in_perm": len(missing_in_perm),
        "missing_in_theo": len(missing_in_theo),
        "mean_abs_diff": safe_mean(abs_diffs),
        "median_abs_diff": statistics.median(abs_diffs) if abs_diffs else float("nan"),
        "max_abs_diff": max(abs_diffs) if abs_diffs else float("nan"),
        "mean_rel_diff": safe_mean(rel_diffs),
        "within_tolerance": within_tolerance,
        "within_tolerance_rate": within_tolerance / len(common_keys) if common_keys else float("nan"),
        "tolerance": tolerance,
        "top_diffs": sorted(diffs, key=lambda item: item["abs_diff"], reverse=True)[:10],
    }
    return summary


def write_diff_csv(path: Path, theo_path: Path, perm_path: Path) -> None:
    theo_headers, theo_rows = load_table(theo_path)
    perm_headers, perm_rows = load_table(perm_path)
    key_columns = detect_key_columns(theo_headers, perm_headers)
    theo_p_column = detect_p_column(theo_headers)
    perm_p_column = detect_p_column(perm_headers)

    theo_map = {row_key(row, key_columns): row for row in theo_rows}
    perm_map = {row_key(row, key_columns): row for row in perm_rows}
    common_keys = sorted(set(theo_map) & set(perm_map))

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["key", "theo_p", "perm_p", "abs_diff", "rel_diff"])
        for key in common_keys:
            theo_p = to_float(theo_map[key][theo_p_column])
            perm_p = to_float(perm_map[key][perm_p_column])
            abs_diff = abs(theo_p - perm_p)
            rel_diff = abs_diff / max(abs(theo_p), abs(perm_p), 1e-12)
            writer.writerow(
                [
                    "|".join(key),
                    f"{theo_p:.10g}",
                    f"{perm_p:.10g}",
                    f"{abs_diff:.10g}",
                    f"{rel_diff:.10g}",
                ]
            )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Time and compare theo vs perm p-value runs."
    )
    parser.add_argument(
        "--theo-cmd",
        default=None,
        help="Shell command that generates the theo result file.",
    )
    parser.add_argument(
        "--perm-cmd",
        default=None,
        help="Shell command that generates the perm result file.",
    )
    parser.add_argument(
        "--theo-result",
        required=True,
        help="Path to the theo result file.",
    )
    parser.add_argument(
        "--perm-result",
        required=True,
        help="Path to the perm result file.",
    )
    parser.add_argument(
        "--summary-out",
        default="benchmark_pvalue_compare.summary.json",
        help="Path to the JSON summary output.",
    )
    parser.add_argument(
        "--diff-out",
        default="benchmark_pvalue_compare.diff.csv",
        help="Path to the detailed diff CSV output.",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=1e-3,
        help="Absolute error tolerance used for the pass rate.",
    )
    args = parser.parse_args()

    timings: Dict[str, float] = {}

    if args.theo_cmd:
        timings["theo_seconds"] = run_shell_command(args.theo_cmd).seconds
    if args.perm_cmd:
        timings["perm_seconds"] = run_shell_command(args.perm_cmd).seconds

    theo_path = Path(args.theo_result)
    perm_path = Path(args.perm_result)
    summary = compare_tables(theo_path, perm_path, args.tolerance)
    summary.update(timings)

    summary_path = Path(args.summary_out)
    diff_path = Path(args.diff_out)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_diff_csv(diff_path, theo_path, perm_path)

    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"Summary written to: {summary_path}")
    print(f"Detailed diffs written to: {diff_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
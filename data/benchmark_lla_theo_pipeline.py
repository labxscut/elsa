#!/usr/bin/env python3
"""CPU-only simulation pipeline for comparing LLA theo vs perm p-values.

This script generates reproducible synthetic triplets, runs lla_compute.py in
the two p-value modes, compares the outputs, and aggregates results across
multiple sequence lengths and replicates.

It is designed to be run inside the same Docker environment used for the main
LLA analysis, but it forces ELSA_COMPCORE_BACKEND=cpu for fairness and to keep
timings comparable across runs.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from gen_unified_triplets import generate_unified_triplet, write_triplet_to_file
from benchmark_pvalue_compare import compare_tables, detect_p_column, load_table


def resolve_window(n: int, window_fraction: float) -> Tuple[Optional[int], Optional[int]]:
    if window_fraction >= 1.0:
        return None, None

    window_length = max(10, int(round(n * window_fraction)))
    window_length = min(window_length, n)
    start = max(0, (n - window_length) // 2)
    end = start + window_length
    return start, end


def run_lla_compute(
    input_file: Path,
    output_file: Path,
    n: int,
    pvalue_method: str,
    delay_limit: int,
    precision: int,
    keep_trace: bool,
) -> float:
    cmd = [
        sys.executable,
        "-m",
        "lla.lla_compute",
        str(input_file),
        str(output_file),
        "-r",
        "1",
        "-s",
        str(n),
        "-d",
        str(delay_limit),
        "-p",
        pvalue_method,
        "-x",
        str(precision),
        "-n",
        "pnz",
        "-f",
        "linear",
    ]
    if keep_trace:
        cmd.append("--keep-trace")

    env = os.environ.copy()
    env["ELSA_COMPCORE_BACKEND"] = "cpu"
    env["PYTHONPATH"] = os.pathsep.join(
        [str(REPO_ROOT), env.get("PYTHONPATH", "")]
        if env.get("PYTHONPATH")
        else [str(REPO_ROOT)]
    )

    start = time.perf_counter()
    completed = subprocess.run(
        cmd,
        cwd=str(REPO_ROOT),
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    elapsed = time.perf_counter() - start
    if completed.returncode != 0:
        raise RuntimeError(
            f"lla_compute.py failed for mode={pvalue_method}, n={n}:\n"
            f"STDOUT:\n{completed.stdout}\nSTDERR:\n{completed.stderr}"
        )
    return elapsed


def extract_first_pvalue(result_path: Path) -> float:
    headers, rows = load_table(result_path)
    if not rows:
        raise ValueError(f"No result rows found in {result_path}")
    p_column = detect_p_column(headers)
    return float(rows[0][p_column])


def safe_mean(values: Sequence[float]) -> float:
    return statistics.fmean(values) if values else float("nan")


def safe_stdev(values: Sequence[float]) -> float:
    return statistics.stdev(values) if len(values) > 1 else 0.0


def write_csv(path: Path, rows: Sequence[Dict[str, object]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return

    fieldnames = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def run_case(
    *,
    n: int,
    replicate_index: int,
    base_seed: int,
    output_dir: Path,
    delay_limit: int,
    precision: int,
    effect_size: float,
    noise_sd: float,
    sigma_y: float,
    delay_yz: int,
    window_fraction: float,
    z_encoding: str,
    keep_trace: bool,
    control: bool,
) -> Dict[str, object]:
    seed = base_seed + n * 10000 + replicate_index
    window_start, window_end = resolve_window(n, window_fraction)

    case_id = f"n{n}_r{replicate_index:03d}_seed{seed}"
    input_dir = output_dir / "inputs"
    result_dir = output_dir / "results"
    report_dir = output_dir / "reports"
    input_dir.mkdir(parents=True, exist_ok=True)
    result_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)

    input_file = input_dir / f"{case_id}.tsv"
    theo_result = result_dir / f"{case_id}_theo.txt"
    perm_result = result_dir / f"{case_id}_perm.txt"
    diff_csv = report_dir / f"{case_id}.diff.csv"
    summary_json = report_dir / f"{case_id}.json"

    X, Y, Z, metadata = generate_unified_triplet(
        n=n,
        method="additive",
        effect_size=effect_size,
        noise_sd=noise_sd,
        sigma_y=sigma_y,
        z_encoding=z_encoding,
        window_start=window_start,
        window_end=window_end,
        delay_yz=delay_yz,
        is_control=control,
        seed=seed,
    )

    write_triplet_to_file(str(input_file), X, Y, Z, labels=("S1", "S2", "S3"))

    theo_seconds = run_lla_compute(
        input_file=input_file,
        output_file=theo_result,
        n=n,
        pvalue_method="theo",
        delay_limit=delay_limit,
        precision=precision,
        keep_trace=keep_trace,
    )
    perm_seconds = run_lla_compute(
        input_file=input_file,
        output_file=perm_result,
        n=n,
        pvalue_method="perm",
        delay_limit=delay_limit,
        precision=precision,
        keep_trace=keep_trace,
    )

    compare_summary = compare_tables(theo_result, perm_result, tolerance=1e-3)

    theo_p = extract_first_pvalue(theo_result)
    perm_p = extract_first_pvalue(perm_result)
    abs_diff = abs(theo_p - perm_p)
    rel_diff = abs_diff / max(abs(theo_p), abs(perm_p), 1e-12)

    case_summary: Dict[str, object] = {
        "case_id": case_id,
        "n": n,
        "replicate_index": replicate_index,
        "seed": seed,
        "delay_limit": delay_limit,
        "precision": precision,
        "effect_size": effect_size,
        "noise_sd": noise_sd,
        "sigma_y": sigma_y,
        "delay_yz": delay_yz,
        "window_fraction": window_fraction,
        "window_start": metadata["window_start"],
        "window_end": metadata["window_end"],
        "control": control,
        "theo_p": theo_p,
        "perm_p": perm_p,
        "abs_diff": abs_diff,
        "rel_diff": rel_diff,
        "theo_seconds": theo_seconds,
        "perm_seconds": perm_seconds,
        "compare_mean_abs_diff": compare_summary["mean_abs_diff"],
        "compare_max_abs_diff": compare_summary["max_abs_diff"],
        "compare_within_tolerance_rate": compare_summary["within_tolerance_rate"],
        "compare_matched_rows": compare_summary["matched_rows"],
    }

    summary_json.write_text(json.dumps(case_summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_csv(diff_csv, [
        {
            "n": n,
            "replicate_index": replicate_index,
            "seed": seed,
            "theo_p": theo_p,
            "perm_p": perm_p,
            "abs_diff": abs_diff,
            "rel_diff": rel_diff,
            "theo_seconds": theo_seconds,
            "perm_seconds": perm_seconds,
            "compare_mean_abs_diff": compare_summary["mean_abs_diff"],
            "compare_max_abs_diff": compare_summary["max_abs_diff"],
            "compare_within_tolerance_rate": compare_summary["within_tolerance_rate"],
        }
    ])

    return case_summary


def aggregate_by_n(records: Sequence[Dict[str, object]]) -> List[Dict[str, object]]:
    grouped: Dict[int, List[Dict[str, object]]] = {}
    for record in records:
        grouped.setdefault(int(record["n"]), []).append(record)

    summary_rows: List[Dict[str, object]] = []
    for n in sorted(grouped):
        subset = grouped[n]
        theo_ps = [float(row["theo_p"]) for row in subset]
        perm_ps = [float(row["perm_p"]) for row in subset]
        abs_diffs = [float(row["abs_diff"]) for row in subset]
        rel_diffs = [float(row["rel_diff"]) for row in subset]
        theo_seconds = [float(row["theo_seconds"]) for row in subset]
        perm_seconds = [float(row["perm_seconds"]) for row in subset]

        summary_rows.append(
            {
                "n": n,
                "replicates": len(subset),
                "mean_theo_p": safe_mean(theo_ps),
                "std_theo_p": safe_stdev(theo_ps),
                "mean_perm_p": safe_mean(perm_ps),
                "std_perm_p": safe_stdev(perm_ps),
                "mean_abs_diff": safe_mean(abs_diffs),
                "std_abs_diff": safe_stdev(abs_diffs),
                "mean_rel_diff": safe_mean(rel_diffs),
                "std_rel_diff": safe_stdev(rel_diffs),
                "mean_theo_seconds": safe_mean(theo_seconds),
                "std_theo_seconds": safe_stdev(theo_seconds),
                "mean_perm_seconds": safe_mean(perm_seconds),
                "std_perm_seconds": safe_stdev(perm_seconds),
                "theo_detect_rate": sum(p < 0.05 for p in theo_ps) / len(theo_ps),
                "perm_detect_rate": sum(p < 0.05 for p in perm_ps) / len(perm_ps),
            }
        )

    return summary_rows


def maybe_plot(summary_rows: Sequence[Dict[str, object]], output_dir: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return

    if not summary_rows:
        return

    ns = [int(row["n"]) for row in summary_rows]
    mean_abs_diff = [float(row["mean_abs_diff"]) for row in summary_rows]
    mean_theo_seconds = [float(row["mean_theo_seconds"]) for row in summary_rows]
    mean_perm_seconds = [float(row["mean_perm_seconds"]) for row in summary_rows]
    mean_theo_p = [float(row["mean_theo_p"]) for row in summary_rows]
    mean_perm_p = [float(row["mean_perm_p"]) for row in summary_rows]

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)

    axes[0, 0].plot(ns, mean_theo_p, marker="o", label="theo")
    axes[0, 0].plot(ns, mean_perm_p, marker="s", label="perm")
    axes[0, 0].set_title("Mean P-value vs n")
    axes[0, 0].set_xlabel("n")
    axes[0, 0].set_ylabel("P-value")
    axes[0, 0].legend()

    axes[0, 1].plot(ns, mean_abs_diff, marker="o", color="tab:red")
    axes[0, 1].set_title("Mean Absolute Difference vs n")
    axes[0, 1].set_xlabel("n")
    axes[0, 1].set_ylabel("|theo - perm|")

    axes[1, 0].plot(ns, mean_theo_seconds, marker="o", label="theo")
    axes[1, 0].plot(ns, mean_perm_seconds, marker="s", label="perm")
    axes[1, 0].set_title("Mean Runtime vs n")
    axes[1, 0].set_xlabel("n")
    axes[1, 0].set_ylabel("seconds")
    axes[1, 0].legend()

    axes[1, 1].plot(ns, [float(row["theo_detect_rate"]) for row in summary_rows], marker="o", label="theo")
    axes[1, 1].plot(ns, [float(row["perm_detect_rate"]) for row in summary_rows], marker="s", label="perm")
    axes[1, 1].set_title("Detection Rate (p < 0.05)")
    axes[1, 1].set_xlabel("n")
    axes[1, 1].set_ylabel("rate")
    axes[1, 1].legend()

    fig.savefig(output_dir / "benchmark_lla_theo_pipeline.png", dpi=200)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate and compare CPU-only simulated LLA theo vs perm runs."
    )
    parser.add_argument("--output-dir", default="data/lla_theo_pipeline", help="Output directory.")
    parser.add_argument("--n-values", nargs="+", type=int, default=[20, 40, 60, 80, 100], help="Sequence lengths to benchmark.")
    parser.add_argument("--replicates", type=int, default=20, help="Replicates per n.")
    parser.add_argument("--base-seed", type=int, default=12345, help="Base random seed.")
    parser.add_argument("--delay-limit", type=int, default=0, help="LLA delay limit. Use 0 for the current benchmark.")
    parser.add_argument("--precision", type=int, default=1000, help="Permutation precision for perm mode and theo lookup granularity.")
    parser.add_argument("--effect-size", type=float, default=1.0, help="Signal amplitude for the additive simulator.")
    parser.add_argument("--noise-sd", type=float, default=0.1, help="Noise SD for the additive simulator.")
    parser.add_argument("--sigma-y", type=float, default=1.0, help="Baseline SD of Y.")
    parser.add_argument("--delay-yz", type=int, default=0, help="Y-Z delay in the simulator.")
    parser.add_argument("--window-fraction", type=float, default=0.8, help="Regulation window fraction of n. Use <1.0 for local signal; 1.0 makes Z constant and is not informative under pnz normalization.")
    parser.add_argument("--z-encoding", choices=["01", "-11"], default="01", help="Encoding for Z.")
    parser.add_argument("--keep-trace", action="store_true", help="Keep trace output from lla_compute.py.")
    parser.add_argument("--control", action="store_true", help="Generate control triplets instead of associated ones.")
    parser.add_argument("--save-intermediates", action="store_true", help="Keep generated input and result files.")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    all_records: List[Dict[str, object]] = []
    for n in args.n_values:
        for rep in range(args.replicates):
            record = run_case(
                n=n,
                replicate_index=rep,
                base_seed=args.base_seed,
                output_dir=output_dir,
                delay_limit=args.delay_limit,
                precision=args.precision,
                effect_size=args.effect_size,
                noise_sd=args.noise_sd,
                sigma_y=args.sigma_y,
                delay_yz=args.delay_yz,
                window_fraction=args.window_fraction,
                z_encoding=args.z_encoding,
                keep_trace=args.keep_trace,
                control=args.control,
            )
            all_records.append(record)

            if not args.save_intermediates:
                case_id = record["case_id"]
                for subdir in ("inputs", "results", "reports"):
                    # The per-case files are useful for auditing; delete them only when requested.
                    for suffix in (".tsv", "_theo.txt", "_perm.txt", ".diff.csv", ".json"):
                        candidate = output_dir / subdir / f"{case_id}{suffix}"
                        if candidate.exists():
                            candidate.unlink()

    summary_rows = aggregate_by_n(all_records)
    summary_csv = output_dir / "benchmark_lla_theo_pipeline.summary.csv"
    summary_json = output_dir / "benchmark_lla_theo_pipeline.summary.json"
    summary_csv_rows: List[Dict[str, object]] = []
    for row in summary_rows:
        summary_csv_rows.append({key: row[key] for key in row.keys()})

    write_csv(summary_csv, summary_csv_rows)
    summary_json.write_text(json.dumps(summary_rows, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    maybe_plot(summary_rows, output_dir)

    print(json.dumps(summary_rows, indent=2, sort_keys=True))
    print(f"Summary CSV: {summary_csv}")
    print(f"Summary JSON: {summary_json}")
    if (output_dir / "benchmark_lla_theo_pipeline.png").exists():
        print(f"Plot: {output_dir / 'benchmark_lla_theo_pipeline.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
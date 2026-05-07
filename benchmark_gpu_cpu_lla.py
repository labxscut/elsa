#!/usr/bin/env python3
"""Benchmark GPU vs CPU runtime for LLA on synthetic mix-mode workloads.

The script generates one full input matrix per benchmark point, with a planted
X/Y/Z triplet in the first three rows and additional independent noise rows to
scale the workload size.

It supports the two benchmark families you asked for:

- fixed sequence length n=100, varying factor count m;
- fixed factor count m=2000, varying sequence length n.

For each generated input, it times the same workload under CPU and GPU backends
and runs both `theo` and `perm` p-value modes on the same synthetic data.
"""

from __future__ import annotations

import argparse
import csv
import os
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np


REPO_ROOT = Path(__file__).resolve().parent


def parse_int_list(raw_value: str | None) -> List[int]:
    if not raw_value:
        return []
    values: List[int] = []
    for item in raw_value.split(","):
        item = item.strip()
        if item:
            values.append(int(item))
    return values


def unique_sorted(values: Iterable[int]) -> List[int]:
    return sorted({int(value) for value in values})


def build_default_points(limit: int) -> List[int]:
    if limit <= 3:
        return [limit]
    points = [max(3, limit // 4), max(3, limit // 2), limit]
    return unique_sorted(point for point in points if 1 <= point <= limit)


def build_seed(base_seed: int, sweep_index: int, value: int, repeat_index: int, method_index: int) -> int:
    return (
        int(base_seed)
        + int(sweep_index) * 10_000_000
        + int(value) * 10_000
        + int(repeat_index) * 100
        + int(method_index)
    )


def generate_mix_matrix(
    *,
    m: int,
    n: int,
    alpha: float,
    window_fraction: float,
    seed: int,
) -> Tuple[List[str], np.ndarray, Dict[str, int]]:
    if m < 3:
        raise ValueError("m must be at least 3 because the planted triplet uses three rows")
    if n < 1:
        raise ValueError("n must be positive")

    rng = np.random.default_rng(seed)

    window_len = max(1, int(np.floor(n * window_fraction)))
    window_len = min(window_len, n)
    window_start = int(np.floor((n - window_len) / 2))
    window_end = window_start + window_len

    y = rng.normal(0.0, 1.0, size=n)
    x = rng.normal(0.0, 1.0, size=n)
    inside_eps = rng.normal(0.0, 0.01, size=window_len)
    x[window_start:window_end] = alpha * y[window_start:window_end] + (1.0 - alpha) * inside_eps

    z = np.zeros(n, dtype=float)
    z[window_start:window_end] = 1.0

    matrix = rng.normal(0.0, 1.0, size=(m, n))
    matrix[0] = x
    matrix[1] = y
    matrix[2] = z

    labels = [f"S{i + 1}" for i in range(m)]
    metadata = {
        "window_start": window_start,
        "window_end": window_end,
        "window_length": window_len,
        "seed": seed,
        "m": m,
        "n": n,
    }
    return labels, matrix, metadata


def write_matrix_to_file(path: Path, labels: Sequence[str], values: np.ndarray) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["#"] + [f"T{i + 1}" for i in range(values.shape[1])])
        for label, row in zip(labels, values):
            writer.writerow([label] + [f"{float(item):.6f}" for item in row])


def run_lla_compute(
    input_file: Path,
    output_file: Path,
    *,
    backend: str,
    pvalue_method: str,
    n: int,
    delay_limit: int,
    precision: int,
    keep_trace: bool,
    fill_method: str,
    norm_method: str,
) -> float:
    print(
        f"[run start] backend={backend} pvalue={pvalue_method} n={n} input={input_file.name}",
        file=sys.stderr,
        flush=True,
    )
    cmd = [
        sys.executable,
        "-u",
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
        norm_method,
        "-f",
        fill_method,
    ]
    if keep_trace:
        cmd.append("--keep-trace")

    env = os.environ.copy()
    env["ELSA_COMPCORE_BACKEND"] = backend
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
    if completed.stdout:
        sys.stdout.write(completed.stdout)
        if not completed.stdout.endswith("\n"):
            sys.stdout.write("\n")
        sys.stdout.flush()
    if completed.stderr:
        sys.stderr.write(completed.stderr)
        if not completed.stderr.endswith("\n"):
            sys.stderr.write("\n")
        sys.stderr.flush()
    if completed.returncode != 0:
        raise RuntimeError(
            f"lla_compute failed for backend={backend}, pvalue_method={pvalue_method}, n={n}:\n"
            f"STDOUT:\n{completed.stdout}\nSTDERR:\n{completed.stderr}"
        )
    print(
        f"[run done] backend={backend} pvalue={pvalue_method} seconds={elapsed:.3f}",
        file=sys.stderr,
        flush=True,
    )
    return elapsed


def write_csv(path: Path, rows: Sequence[Dict[str, object]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return

    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def aggregate_runs(rows: Sequence[Dict[str, object]]) -> List[Dict[str, object]]:
    collapsed: Dict[Tuple[str, int, str, str], List[float]] = {}
    for row in rows:
        key = (
            str(row["sweep_mode"]),
            int(row["value"]),
            str(row["backend"]),
            str(row["pvalue_method"]),
        )
        collapsed.setdefault(key, []).append(float(row["seconds"]))

    summary_rows: List[Dict[str, object]] = []
    cpu_means: Dict[Tuple[str, int, str], float] = {}

    for (sweep_mode, value, backend, pmethod), seconds_list in sorted(collapsed.items()):
        mean_seconds = statistics.fmean(seconds_list)
        std_seconds = statistics.stdev(seconds_list) if len(seconds_list) > 1 else 0.0
        summary_rows.append(
            {
                "sweep_mode": sweep_mode,
                "value": value,
                "backend": backend,
                "pvalue_method": pmethod,
                "runs": len(seconds_list),
                "mean_seconds": mean_seconds,
                "std_seconds": std_seconds,
                "min_seconds": min(seconds_list),
                "max_seconds": max(seconds_list),
            }
        )
        if backend == "cpu":
            cpu_means[(sweep_mode, value, pmethod)] = mean_seconds

    for row in summary_rows:
        if row["backend"] == "gpu":
            cpu_mean = cpu_means.get((str(row["sweep_mode"]), int(row["value"]), str(row["pvalue_method"])))
            row["speedup_vs_cpu"] = cpu_mean / float(row["mean_seconds"]) if cpu_mean else ""
        else:
            row["speedup_vs_cpu"] = ""

    return summary_rows


def plot_summary(summary_rows: Sequence[Dict[str, object]], output_path: Path) -> None:
    if not summary_rows:
        return

    sweep_modes = sorted({str(row["sweep_mode"]) for row in summary_rows})
    pvalue_methods = sorted({str(row["pvalue_method"]) for row in summary_rows})

    fig, axes = plt.subplots(
        len(sweep_modes),
        len(pvalue_methods),
        figsize=(7 * len(pvalue_methods), 4.8 * len(sweep_modes)),
        squeeze=False,
    )
    fig.suptitle("LLA CPU vs GPU runtime on synthetic mix workloads", fontsize=15, fontweight="bold")

    x_labels = {
        "rows": "m (factor count)",
        "length": "n (sequence length)",
    }

    for i, sweep_mode in enumerate(sweep_modes):
        for j, pmethod in enumerate(pvalue_methods):
            axis = axes[i][j]
            panel_rows = [
                row
                for row in summary_rows
                if str(row["sweep_mode"]) == sweep_mode and str(row["pvalue_method"]) == pmethod
            ]
            values = sorted({int(row["value"]) for row in panel_rows})
            cpu_mean: List[float] = []
            cpu_std: List[float] = []
            gpu_mean: List[float] = []
            gpu_std: List[float] = []

            for value in values:
                cpu_row = next(
                    (row for row in panel_rows if int(row["value"]) == value and row["backend"] == "cpu"),
                    None,
                )
                gpu_row = next(
                    (row for row in panel_rows if int(row["value"]) == value and row["backend"] == "gpu"),
                    None,
                )
                cpu_mean.append(float(cpu_row["mean_seconds"]) if cpu_row else np.nan)
                cpu_std.append(float(cpu_row["std_seconds"]) if cpu_row else np.nan)
                gpu_mean.append(float(gpu_row["mean_seconds"]) if gpu_row else np.nan)
                gpu_std.append(float(gpu_row["std_seconds"]) if gpu_row else np.nan)

            axis.errorbar(values, cpu_mean, yerr=cpu_std, marker="o", linewidth=2, capsize=4, label="CPU")
            axis.errorbar(values, gpu_mean, yerr=gpu_std, marker="s", linewidth=2, capsize=4, label="GPU")
            axis.set_xlabel(x_labels.get(sweep_mode, "workload size"))
            axis.set_ylabel("Wall time (seconds)")
            axis.set_title(f"{sweep_mode} sweep, pvalue={pmethod}")
            axis.grid(alpha=0.3)
            axis.legend()

    plt.tight_layout(rect=[0, 0.02, 1, 0.95])
    fig.savefig(output_path, dpi=300, bbox_inches="tight")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare GPU and CPU runtime for LLA using synthetic mix-mode workloads."
    )
    parser.add_argument(
        "--row-points",
        type=str,
        default="100,200,500,1000,2000",
        help="Comma-separated m values for the fixed-n sweep (default: 100,200,500,1000,2000).",
    )
    parser.add_argument(
        "--length-points",
        type=str,
        default="100,200,500,1000",
        help="Comma-separated n values for the fixed-m sweep (default: 100,200,500,1000).",
    )
    parser.add_argument(
        "--fixed-length",
        type=int,
        default=100,
        help="Fixed n used while sweeping m (default: 100).",
    )
    parser.add_argument(
        "--fixed-rows",
        type=int,
        default=2000,
        help="Fixed m used while sweeping n (default: 2000).",
    )
    parser.add_argument(
        "--repeats",
        type=int,
        default=2,
        help="Number of generated inputs per benchmark point (default: 2).",
    )
    parser.add_argument(
        "--pvalue-methods",
        type=str,
        default="theo,perm",
        help="Comma-separated p-value methods to benchmark (default: theo,perm).",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.8,
        help="Association strength inside the window (default: 0.6).",
    )
    parser.add_argument(
        "--window-fraction",
        type=float,
        default=0.6,
        help="Fraction of n used for the centered local window (default: 0.8).",
    )
    parser.add_argument(
        "--base-seed",
        type=int,
        default=12345,
        help="Base seed for synthetic data generation (default: 12345).",
    )
    parser.add_argument(
        "--delay-limit",
        type=int,
        default=0,
        help="Delay search limit passed to lla_compute (default: 0).",
    )
    parser.add_argument(
        "--precision",
        type=int,
        default=1000,
        help="P-value precision passed to lla_compute (default: 1000).",
    )
    parser.add_argument(
        "--fill-method",
        type=str,
        default="linear",
        choices=["none", "zero", "linear", "quadratic", "cubic", "slinear", "nearest"],
        help="Missing-value fill method passed to lla_compute (default: linear).",
    )
    parser.add_argument(
        "--norm-method",
        type=str,
        default="pnz",
        choices=["percentile", "pnz", "none"],
        help="Normalization method passed to lla_compute (default: pnz).",
    )
    parser.add_argument(
        "--backends",
        type=str,
        default="cpu,gpu",
        help="Comma-separated backends to run (default: cpu,gpu).",
    )
    keep_trace_group = parser.add_mutually_exclusive_group()
    keep_trace_group.add_argument(
        "--keep-trace",
        dest="keep_trace",
        action="store_true",
        default=True,
        help="Forward --keep-trace to lla_compute (default: enabled).",
    )
    keep_trace_group.add_argument(
        "--no-keep-trace",
        dest="keep_trace",
        action="store_false",
        help="Disable --keep-trace for lla_compute.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("gpu_cpu_benchmark_results"),
        help="Output directory for generated inputs, CSVs, and plots (default: gpu_cpu_benchmark_results).",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Save a PNG plot in addition to the CSV files.",
    )

    args = parser.parse_args()

    row_points = unique_sorted(parse_int_list(args.row_points))
    length_points = unique_sorted(parse_int_list(args.length_points))
    if not row_points:
        row_points = build_default_points(args.fixed_rows)
    if not length_points:
        length_points = build_default_points(args.fixed_length)

    backends = [backend.strip().lower() for backend in args.backends.split(",") if backend.strip()]
    if not backends:
        raise SystemExit("At least one backend must be provided")
    invalid_backends = [backend for backend in backends if backend not in {"cpu", "gpu"}]
    if invalid_backends:
        raise SystemExit(f"Unsupported backend values: {', '.join(invalid_backends)}")

    pvalue_methods = [method.strip().lower() for method in args.pvalue_methods.split(",") if method.strip()]
    if not pvalue_methods:
        raise SystemExit("At least one p-value method must be provided")
    invalid_pmethods = [method for method in pvalue_methods if method not in {"theo", "perm", "mix"}]
    if invalid_pmethods:
        raise SystemExit(f"Unsupported p-value modes: {', '.join(invalid_pmethods)}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    inputs_dir = args.output_dir / "inputs"
    inputs_dir.mkdir(parents=True, exist_ok=True)
    raw_csv = args.output_dir / "gpu_cpu_lla_raw.csv"
    summary_csv = args.output_dir / "gpu_cpu_lla_summary.csv"
    plot_png = args.output_dir / "gpu_cpu_lla_summary.png"

    print(f"Fixed length for row sweep: n={args.fixed_length}", file=sys.stderr)
    print(f"Fixed rows for length sweep: m={args.fixed_rows}", file=sys.stderr)
    print(f"Row sweep points (m): {row_points}", file=sys.stderr)
    print(f"Length sweep points (n): {length_points}", file=sys.stderr)
    print(f"Backends: {backends}", file=sys.stderr)
    print(f"P-value methods: {pvalue_methods}", file=sys.stderr)

    raw_rows: List[Dict[str, object]] = []
    start_time = time.perf_counter()

    sweep_specs = [
        ("rows", row_points, args.fixed_length),
        ("length", length_points, args.fixed_rows),
    ]

    for sweep_index, (sweep_mode, points, fixed_value) in enumerate(sweep_specs):
        for value in points:
            m = value if sweep_mode == "rows" else args.fixed_rows
            n = value if sweep_mode == "length" else args.fixed_length
            if m < 3:
                continue

            for repeat_index in range(args.repeats):
                seed = build_seed(args.base_seed, sweep_index, value, repeat_index, 0)
                labels, matrix, metadata = generate_mix_matrix(
                    m=m,
                    n=n,
                    alpha=args.alpha,
                    window_fraction=args.window_fraction,
                    seed=seed,
                )

                input_file = inputs_dir / f"{sweep_mode}_v{value}_rep{repeat_index}_seed{seed}.tsv"
                write_matrix_to_file(input_file, labels, matrix)

                for pmethod in pvalue_methods:
                    for backend in backends:
                        output_file = args.output_dir / f"{sweep_mode}_v{value}_rep{repeat_index}_{backend}_{pmethod}.txt"
                        elapsed = run_lla_compute(
                            input_file,
                            output_file,
                            backend=backend,
                            pvalue_method=pmethod,
                            n=n,
                            delay_limit=args.delay_limit,
                            precision=args.precision,
                            keep_trace=args.keep_trace,
                            fill_method=args.fill_method,
                            norm_method=args.norm_method,
                        )
                        raw_rows.append(
                            {
                                "sweep_mode": sweep_mode,
                                "value": value,
                                "fixed_value": fixed_value,
                                "backend": backend,
                                "pvalue_method": pmethod,
                                "repeat": repeat_index,
                                "seconds": elapsed,
                                "seed": seed,
                                "alpha": args.alpha,
                                "window_fraction": args.window_fraction,
                                "window_start": metadata["window_start"],
                                "window_end": metadata["window_end"],
                                "window_length": metadata["window_length"],
                                "m": m,
                                "n": n,
                                "delay_limit": args.delay_limit,
                                "precision": args.precision,
                                "keep_trace": args.keep_trace,
                                "fill_method": args.fill_method,
                                "norm_method": args.norm_method,
                                "input_file": str(input_file),
                                "output_file": str(output_file),
                            }
                        )

    summary_rows = aggregate_runs(raw_rows)
    elapsed = time.perf_counter() - start_time

    write_csv(raw_csv, raw_rows)
    write_csv(summary_csv, summary_rows)
    if args.plot:
        plot_summary(summary_rows, plot_png)

    print(f"Saved raw timings to: {raw_csv}")
    print(f"Saved summary to: {summary_csv}")
    if args.plot:
        print(f"Saved plot to: {plot_png}")
    print(f"Total elapsed time: {elapsed:.2f} seconds")

    for row in summary_rows:
        if row["backend"] == "gpu" and row["speedup_vs_cpu"] != "":
            print(
                f"[{row['sweep_mode']}={row['value']}, {row['pvalue_method']}] "
                f"CPU/GPU speedup={row['speedup_vs_cpu']:.2f}x"
            )


if __name__ == "__main__":
    main()

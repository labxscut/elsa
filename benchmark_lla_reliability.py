#!/usr/bin/env python
"""
benchmark_lla_reliability.py

Quantify how accurately and consistently LLA finds the true local association
interval when delay = 0.

Metrics:
- Jaccard (overlap ratio): J = L_int / L_union
- Coverage: fraction of true region captured
- Boundary errors: MAE for start and end positions
- Power: detection rate at alpha = 0.05

Usage:
    python benchmark_lla_reliability.py --config-csv configs.csv --output-csv results.csv
    
Or use default configuration:
    python benchmark_lla_reliability.py
"""

import csv
import os
import sys
import tempfile
import subprocess
import argparse
from pathlib import Path
from statistics import mean, stdev
from math import sqrt
import time


def jaccard_and_coverage(s_true, e_true, s_pred, e_pred):
    """
    Compute Jaccard and Coverage metrics.
    
    Args:
        s_true, e_true: True interval (1-based, inclusive)
        s_pred, e_pred: Predicted interval (1-based, inclusive)
    
    Returns:
        J, Cov, L_int, L_union
    """
    # Handle empty prediction
    if s_pred is None or e_pred is None or e_pred < s_pred:
        return 0.0, 0.0, 0, (e_true - s_true + 1)
    
    # Compute intersection length
    L_int = max(0, min(e_true, e_pred) - max(s_true, s_pred) + 1)
    
    # Compute union length
    L_true = e_true - s_true + 1
    L_pred = e_pred - s_pred + 1
    L_union = L_true + L_pred - L_int
    
    # Compute metrics
    J = 0.0 if L_union == 0 else L_int / L_union
    Cov = 0.0 if L_true == 0 else L_int / L_true
    
    return J, Cov, L_int, L_union


def parse_result_tsv(path):
    """
    Parse LLA result file and extract metrics.
    
    Returns:
        dict with keys: Start_X, Start_Y, Start_Z, End_X, End_Y, End_Z,
                        Delay_Y-X, Delay_Z-Y, P, LA, or None if parsing fails
    """
    try:
        with open(path, 'r', newline='', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter='\t')
            header = next(reader)
            row = next(reader)  # Only one triple S1, S2, S3
        
        col = {name: i for i, name in enumerate(header)}
        
        def to_int(v):
            if v in ("", "nan", "NaN"):
                return None
            try:
                return int(v)
            except ValueError:
                return None
        
        def to_float(v):
            if v in ("", "nan", "NaN"):
                return float("nan")
            try:
                return float(v)
            except ValueError:
                return float("nan")
        
        res = {
            "Start_X": to_int(row[col["Start_X"]]),
            "Start_Y": to_int(row[col["Start_Y"]]),
            "Start_Z": to_int(row[col["Start_Z"]]),
            "End_X": to_int(row[col["End_X"]]),
            "End_Y": to_int(row[col["End_Y"]]),
            "End_Z": to_int(row[col["End_Z"]]),
            "Delay_Y-X": to_int(row[col["Delay_Y-X"]]) if "Delay_Y-X" in col else 0,
            "Delay_Z-Y": to_int(row[col["Delay_Z-Y"]]) if "Delay_Z-Y" in col else 0,
            "P": to_float(row[col["P"]]) if "P" in col else float("nan"),
            "LA": to_float(row[col["LA"]]) if "LA" in col else float("nan"),
        }
        return res
    except Exception as e:
        print(f"Error parsing result file {path}: {e}", file=sys.stderr)
        return None


def consensus_interval_no_delay(r):
    """
    Compute consensus predicted interval from Start_X/Y/Z and End_X/Y/Z.
    
    For no-delay case, the consensus is the intersection of all three windows.
    
    Returns:
        s_pred, e_pred (1-based, inclusive), or (None, None) if invalid
    """
    if any(v is None for v in [r["Start_X"], r["Start_Y"], r["Start_Z"],
                                 r["End_X"], r["End_Y"], r["End_Z"]]):
        return None, None
    
    s_pred = max(r["Start_X"], r["Start_Y"], r["Start_Z"])
    e_pred = min(r["End_X"], r["End_Y"], r["End_Z"])
    
    # If no overlap, return None
    if e_pred < s_pred:
        return None, None
    
    return s_pred, e_pred


def run_one(n, corr_start0, corr_end0, seed, 
            python_exe=sys.executable, alpha=0.05, 
            delayLimit=3, precision=1000):
    """
    Run one simulation + LLA analysis iteration.
    
    Args:
        n: Number of time points
        corr_start0: Correlation start index (0-based, inclusive)
        corr_end0: Correlation end index (0-based, inclusive)
        seed: Random seed
        python_exe: Python executable path
        alpha: Significance level for power calculation
        delayLimit: Maximum delay for LLA
        precision: Number of permutations for p-value
    
    Returns:
        dict with metrics, or None if failed
    """
    # Convert generator's 0-based indices to 1-based for comparison
    s_true = corr_start0 + 1
    e_true = corr_end0 + 1
    
    # Get script directory (where benchmark_lla_reliability.py is located)
    script_dir = Path(__file__).parent.absolute()
    
    with tempfile.TemporaryDirectory() as td:
        sim_path = os.path.join(td, "sim.txt")
        out_path = os.path.join(td, "res.txt")
        
        # 1) Generate data (no delay)
        gen_script = script_dir / "localsim_with_delay.py"
        gen_cmd = [
            python_exe, str(gen_script),
            "--out", sim_path,
            "--n", str(n),
            "--corr-start", str(corr_start0),
            "--corr-end", str(corr_end0),
            "--seed", str(seed),
            "--delay_xy", "0",
            "--delay_yz", "0"
        ]
        
        try:
            subprocess.run(gen_cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            print(f"Generation failed for seed={seed}: {e.stderr}", file=sys.stderr)
            return None
        
        # 2) Run LLA
        lla_script = script_dir / "lla" / "lla_compute.py"
        lla_cmd = [
            python_exe, str(lla_script),
            sim_path, out_path,
            "-s", str(n), "-r", "1",
            "-d", str(delayLimit),
            "-p", "perm", "-x", str(precision),
            "-n", "pnz", "-f", "linear",
            "--keep-trace"  # Required to get Start/End positions
        ]
        
        try:
            subprocess.run(lla_cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            print(f"LLA failed for seed={seed}: {e.stderr}", file=sys.stderr)
            return None
        
        # 3) Parse result
        r = parse_result_tsv(out_path)
        if r is None:
            return None
        
        s_pred, e_pred = consensus_interval_no_delay(r)
        
        # 4) Compute metrics
        J, Cov, L_int, L_union = jaccard_and_coverage(s_true, e_true, s_pred, e_pred)
        
        mae_start = abs(s_pred - s_true) if s_pred is not None else None
        mae_end = abs(e_pred - e_true) if e_pred is not None else None
        
        # Power: significant at alpha?
        power_hit = (r["P"] <= alpha) if r["P"] == r["P"] else False  # handle NaN
        
        return {
            "n": n,
            "seed": seed,
            "corr_start0": corr_start0,
            "corr_end0": corr_end0,
            "s_true": s_true,
            "e_true": e_true,
            "s_pred": s_pred,
            "e_pred": e_pred,
            "J": J,
            "Coverage": Cov,
            "L_int": L_int,
            "L_union": L_union,
            "MAE_start": mae_start,
            "MAE_end": mae_end,
            "P": r["P"],
            "LA": r["LA"],
            "power_hit": power_hit
        }


def summarize(rows, alpha=0.05):
    """
    Compute summary statistics across multiple runs.
    
    Returns:
        dict with mean/SE for each metric
    """
    if not rows:
        return {}
    
    R = len(rows)
    
    def safe_mean(key):
        vals = [r[key] for r in rows if r[key] is not None and r[key] == r[key]]  # not None and not NaN
        return mean(vals) if vals else float("nan")
    
    def safe_stdev(key):
        vals = [r[key] for r in rows if r[key] is not None and r[key] == r[key]]
        return stdev(vals) if len(vals) > 1 else 0.0
    
    # Power (proportion)
    power_vals = [1.0 if r["power_hit"] else 0.0 for r in rows]
    p = mean(power_vals)
    se_power = sqrt(p * (1 - p) / R) if R > 0 else float("nan")
    
    summary = {
        "n": rows[0]["n"],
        "corr_len": rows[0]["e_true"] - rows[0]["s_true"] + 1,
        "R": R,
        "J_mean": safe_mean("J"),
        "J_std": safe_stdev("J"),
        "Coverage_mean": safe_mean("Coverage"),
        "Coverage_std": safe_stdev("Coverage"),
        "MAE_start_mean": safe_mean("MAE_start"),
        "MAE_start_std": safe_stdev("MAE_start"),
        "MAE_end_mean": safe_mean("MAE_end"),
        "MAE_end_std": safe_stdev("MAE_end"),
        "Power": p,
        "Power_SE": se_power,
        "LA_mean": safe_mean("LA"),
        "LA_std": safe_stdev("LA")
    }
    
    return summary


def run_benchmark(configs, python_exe=sys.executable, alpha=0.05, 
                  delayLimit=0, precision=1000, per_run_csv=None):
    """
    Run benchmark across multiple configurations.
    
    Args:
        configs: List of dicts with keys: n, corr_len, R
        python_exe: Python executable
        alpha: Significance level
        delayLimit: Max delay for LLA
        precision: Permutations
        per_run_csv: Optional path to save per-run results
    
    Returns:
        List of summary dicts
    """
    summaries = []
    all_runs = []
    
    for cfg in configs:
        n = cfg["n"]
        corr_len = cfg["corr_len"]
        R = cfg["R"]
        
        # Validate that corr_len fits within n
        if corr_len >= n:
            print(f"WARNING: corr_len={corr_len} >= n={n}, skipping this config", file=sys.stderr)
            continue
        
        # Center the correlation window
        corr_start0 = n // 2 - corr_len // 2
        corr_end0 = corr_start0 + corr_len - 1
        
        # Ensure it fits within valid range [0, n)
        if corr_end0 >= n:
            corr_end0 = n - 1
            corr_start0 = corr_end0 - corr_len + 1
        
        print(f"\n=== Running n={n}, corr_len={corr_len}, R={R} ===", file=sys.stderr)
        print(f"    True region (0-based): [{corr_start0}, {corr_end0}]", file=sys.stderr)
        print(f"    True region (1-based): [{corr_start0 + 1}, {corr_end0 + 1}]", file=sys.stderr)
        
        rows = []
        for seed in range(1, R + 1):
            if seed % 10 == 0:
                print(f"    Progress: {seed}/{R}", file=sys.stderr)
            
            result = run_one(n, corr_start0, corr_end0, seed, 
                            python_exe=python_exe, alpha=alpha,
                            delayLimit=delayLimit, precision=precision)
            
            if result is not None:
                rows.append(result)
                all_runs.append(result)
        
        if rows:
            summary = summarize(rows, alpha=alpha)
            summaries.append(summary)
            
            print(f"\n    Results for n={n}:", file=sys.stderr)
            print(f"      Jaccard: {summary['J_mean']:.4f} ± {summary['J_std']:.4f}", file=sys.stderr)
            print(f"      Coverage: {summary['Coverage_mean']:.4f} ± {summary['Coverage_std']:.4f}", file=sys.stderr)
            print(f"      MAE_start: {summary['MAE_start_mean']:.2f} ± {summary['MAE_start_std']:.2f}", file=sys.stderr)
            print(f"      MAE_end: {summary['MAE_end_mean']:.2f} ± {summary['MAE_end_std']:.2f}", file=sys.stderr)
            print(f"      Power@{alpha}: {summary['Power']:.4f} ± {summary['Power_SE']:.4f}", file=sys.stderr)
        else:
            print(f"    WARNING: No successful runs for n={n}", file=sys.stderr)
    
    # Optionally save per-run data
    if per_run_csv and all_runs:
        with open(per_run_csv, 'w', newline='', encoding='utf-8') as f:
            fieldnames = ["n", "seed", "corr_start0", "corr_end0", 
                         "s_true", "e_true", "s_pred", "e_pred",
                         "J", "Coverage", "L_int", "L_union",
                         "MAE_start", "MAE_end", "P", "LA", "power_hit"]
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_runs)
        print(f"\nPer-run results saved to: {per_run_csv}", file=sys.stderr)
    
    return summaries


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark LLA reliability for no-delay case",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
    python benchmark_lla_reliability.py --output-csv summary.csv --per-run-csv runs.csv
    
This will run the default configuration:
    n in {20, 40, 80, 100}
    corr_len = 8
    R = 50 runs per configuration
        """
    )
    
    parser.add_argument("--output-csv", dest="output_csv", 
                       default="lla_reliability_summary.csv",
                       help="Output CSV for summary statistics (default: lla_reliability_summary.csv)")
    parser.add_argument("--per-run-csv", dest="per_run_csv",
                       default=None,
                       help="Optional output CSV for per-run details")
    parser.add_argument("--alpha", type=float, default=0.05,
                       help="Significance level for power (default: 0.05)")
    parser.add_argument("--delay-limit", dest="delayLimit", type=int, default=0,
                       help="Maximum delay for LLA (default: 0)")
    parser.add_argument("--precision", type=int, default=1000,
                       help="Permutations for p-value (default: 1000)")
    parser.add_argument("--python", dest="python_exe", default=sys.executable,
                       help=f"Python executable (default: {sys.executable})")
    
    args = parser.parse_args()
    
    # Default configuration grid - varying both n and corr_len
    # Test how performance changes with:
    # 1. Absolute signal length (corr_len)
    # 2. Relative signal proportion (corr_len/n)
    configs = [
        # Short series (n=20)
        # {"n": 20, "corr_len": 6, "R": 50},   # 30% of series
        # {"n": 20, "corr_len": 8, "R": 50},   # 40% of series
        # {"n": 20, "corr_len": 10, "R": 50},  # 50% of series
        
        # Medium series (n=40)
        {"n": 40, "corr_len": 6, "R": 50},   # 15% of series
        {"n": 40, "corr_len": 8, "R": 50},   # 20% of series
        {"n": 40, "corr_len": 12, "R": 50},  # 30% of series
        {"n": 40, "corr_len": 16, "R": 50},  # 40% of series
        
        # # Long series (n=80)
        # {"n": 80, "corr_len": 8, "R": 50},   # 10% of series
        # {"n": 80, "corr_len": 16, "R": 50},  # 20% of series
        # {"n": 80, "corr_len": 24, "R": 50},  # 30% of series
        
        # # Very long series (n=100)
        # {"n": 100, "corr_len": 10, "R": 50}, # 10% of series
        # {"n": 100, "corr_len": 20, "R": 50}, # 20% of series
        # {"n": 100, "corr_len": 30, "R": 50}, # 30% of series
    ]
    
    print("=" * 70, file=sys.stderr)
    print("LLA Reliability Benchmark (No-Delay Case)", file=sys.stderr)
    print("=" * 70, file=sys.stderr)
    print(f"Alpha: {args.alpha}", file=sys.stderr)
    print(f"Delay limit: {args.delayLimit}", file=sys.stderr)
    print(f"Precision: {args.precision}", file=sys.stderr)
    print(f"Python: {args.python_exe}", file=sys.stderr)
    print(f"Configurations: {len(configs)}", file=sys.stderr)
    print(f"\nConfiguration matrix:", file=sys.stderr)
    print(f"  n values: {sorted(set(c['n'] for c in configs))}", file=sys.stderr)
    print(f"  corr_len values: {sorted(set(c['corr_len'] for c in configs))}", file=sys.stderr)
    
    proportions = [f"{c['corr_len']}/{c['n']}" for c in configs]
    print(f"  Signal proportions tested: {sorted(set(proportions))}", file=sys.stderr)
    
    start_time = time.time()
    
    summaries = run_benchmark(
        configs,
        python_exe=args.python_exe,
        alpha=args.alpha,
        delayLimit=args.delayLimit,
        precision=args.precision,
        per_run_csv=args.per_run_csv
    )
    
    # Save summaries
    if summaries:
        with open(args.output_csv, 'w', newline='', encoding='utf-8') as f:
            fieldnames = ["n", "corr_len", "R", 
                         "J_mean", "J_std",
                         "Coverage_mean", "Coverage_std",
                         "MAE_start_mean", "MAE_start_std",
                         "MAE_end_mean", "MAE_end_std",
                         "Power", "Power_SE",
                         "LA_mean", "LA_std"]
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(summaries)
        
        print(f"\n{'=' * 70}", file=sys.stderr)
        print(f"Summary results saved to: {args.output_csv}", file=sys.stderr)
        print(f"Elapsed time: {time.time() - start_time:.1f} seconds", file=sys.stderr)
        print(f"{'=' * 70}", file=sys.stderr)
    else:
        print("\nNo results to save!", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

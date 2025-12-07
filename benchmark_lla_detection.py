#!/usr/bin/env python
"""
benchmark_lla_detection.py

Test the SIMPLEST question: Can LLA distinguish associated vs unassociated triplets?

This is a sanity check that:
1. LLA detects global association (Power at alpha=0.05)
2. LLA correctly rejects pure noise (FPR ≈ alpha)

Usage:
    python benchmark_lla_detection.py --output-csv detection_results.csv
"""

import argparse
import csv
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from statistics import mean
from math import sqrt


def parse_result_tsv(path):
    """Parse LLA result file and extract P and LA score."""
    try:
        with open(path, 'r', newline='', encoding='utf-8') as f:
            lines = f.readlines()
        
        if len(lines) < 2:
            return None
        
        # Split by whitespace (handles multiple spaces)
        header = lines[0].strip().split()
        row = lines[1].strip().split()
        
        # Find columns
        col = {name: i for i, name in enumerate(header)}
        
        def to_float(v):
            if v in ("", "nan", "NaN"):
                return float("nan")
            try:
                return float(v)
            except ValueError:
                return float("nan")
        
        p_val = to_float(row[col["P"]]) if "P" in col else float("nan")
        la_score = to_float(row[col["LA"]]) if "LA" in col else float("nan")
        
        return {"P": p_val, "LA": la_score}
    except Exception as e:
        print(f"Error parsing result file {path}: {e}", file=sys.stderr)
        return None


def run_one(mode, n, seed, python_exe=sys.executable, 
            precision=1000, alpha=0.05):
    """
    Run one LLA detection test.
    
    Args:
        mode: "global_assoc" or "null"
        n: Number of time points
        seed: Random seed
        python_exe: Python executable path
        precision: Number of permutations
        alpha: Significance level
    
    Returns:
        dict with: mode, seed, P, LA, detected
    """
    script_dir = Path(__file__).parent.absolute()
    
    with tempfile.TemporaryDirectory() as td:
        sim_path = os.path.join(td, "sim.txt")
        out_path = os.path.join(td, "res.txt")
        
        # 1) Generate data
        gen_script = script_dir / "gen_global_triplets.py"
        gen_cmd = [
            python_exe, str(gen_script),
            "--out", sim_path,
            "--mode", mode,
            "--n", str(n),
            "--seed", str(seed)
        ]
        
        try:
            subprocess.run(gen_cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            print(f"Generation failed for {mode}, seed={seed}: {e.stderr}", 
                  file=sys.stderr)
            return None
        
        # 2) Run LLA with delayLimit=0 (no delay)
        lla_script = script_dir / "lla" / "lla_compute.py"
        lla_cmd = [
            python_exe, str(lla_script),
            sim_path, out_path,
            "-s", str(n), "-r", "1",
            "-d", "0",  # delayLimit=0 for simplest case
            "-p", "perm", "-x", str(precision),
            "-n", "pnz", "-f", "linear",
            "--keep-trace"
        ]
        
        try:
            subprocess.run(lla_cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            print(f"LLA failed for {mode}, seed={seed}: {e.stderr}", 
                  file=sys.stderr)
            return None
        
        # 3) Parse result
        r = parse_result_tsv(out_path)
        if r is None:
            return None
        
        # 4) Determine if detected
        detected = (r["P"] <= alpha) if r["P"] == r["P"] else False  # handle NaN
        
        return {
            "mode": mode,
            "n": n,
            "seed": seed,
            "P": r["P"],
            "LA": r["LA"],
            "detected": detected
        }


def run_benchmark(n, R_assoc, R_null, python_exe=sys.executable,
                  precision=1000, alpha=0.05):
    """
    Run full detection benchmark.
    
    Args:
        n: Number of time points
        R_assoc: Number of associated triplets to test
        R_null: Number of null triplets to test
        python_exe: Python executable
        precision: Permutations
        alpha: Significance level
    
    Returns:
        (all_runs, summary)
    """
    all_runs = []
    
    # Test associated triplets
    print(f"\nTesting {R_assoc} associated triplets...", file=sys.stderr)
    for seed in range(1, R_assoc + 1):
        if seed % 10 == 0:
            print(f"  Associated: {seed}/{R_assoc}", file=sys.stderr)
        
        result = run_one("global_assoc", n, seed, python_exe, precision, alpha)
        if result:
            all_runs.append(result)
    
    # Test null triplets
    print(f"\nTesting {R_null} null triplets...", file=sys.stderr)
    for seed in range(1, R_null + 1):
        if seed % 10 == 0:
            print(f"  Null: {seed}/{R_null}", file=sys.stderr)
        
        result = run_one("null", n, seed, python_exe, precision, alpha)
        if result:
            all_runs.append(result)
    
    # Compute summary
    assoc_runs = [r for r in all_runs if r["mode"] == "global_assoc"]
    null_runs = [r for r in all_runs if r["mode"] == "null"]
    
    # Power (for associated)
    power_vals = [1.0 if r["detected"] else 0.0 for r in assoc_runs]
    power = mean(power_vals) if power_vals else float("nan")
    power_se = sqrt(power * (1 - power) / len(power_vals)) if power_vals else float("nan")
    
    # FPR (for null)
    fpr_vals = [1.0 if r["detected"] else 0.0 for r in null_runs]
    fpr = mean(fpr_vals) if fpr_vals else float("nan")
    fpr_se = sqrt(fpr * (1 - fpr) / len(fpr_vals)) if fpr_vals else float("nan")
    
    # Mean LA scores
    la_assoc = [r["LA"] for r in assoc_runs if r["LA"] == r["LA"]]
    la_null = [r["LA"] for r in null_runs if r["LA"] == r["LA"]]
    
    summary = {
        "n": n,
        "R_assoc": len(assoc_runs),
        "R_null": len(null_runs),
        "alpha": alpha,
        "precision": precision,
        "Power": power,
        "Power_SE": power_se,
        "FPR": fpr,
        "FPR_SE": fpr_se,
        "LA_mean_assoc": mean(la_assoc) if la_assoc else float("nan"),
        "LA_mean_null": mean(la_null) if la_null else float("nan")
    }
    
    return all_runs, summary


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark LLA detection: can it distinguish associated vs null?",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
    python benchmark_lla_detection.py --output-dir results/
    
This tests the SIMPLEST question:
    - Can LLA detect global association? (Power)
    - Does LLA control false positives on pure noise? (FPR ≈ alpha)
    
Will run multiple n configurations and save results separately.
        """
    )
    
    parser.add_argument("--output-dir", dest="output_dir",
                       default=".",
                       help="Output directory for results (default: current dir)")
    parser.add_argument("--save-per-run", dest="save_per_run", action='store_true',
                       help="Save per-run details for each n")
    parser.add_argument("--alpha", type=float, default=0.05,
                       help="Significance level (default: 0.05)")
    parser.add_argument("--precision", type=int, default=1000,
                       help="Permutations for p-value (default: 1000)")
    
    args = parser.parse_args()
    
    # Configuration matrix
    configs = [
        {"n": 40,  "R_assoc": 50,  "R_null": 50},
        {"n": 80,  "R_assoc": 50,  "R_null": 50},
        {"n": 100, "R_assoc": 50,  "R_null": 50}
    ]
    
    # Use current Python executable
    python_exe = sys.executable
    
    # Create output directory if needed
    os.makedirs(args.output_dir, exist_ok=True)
    
    print("=" * 70, file=sys.stderr)
    print("LLA Detection Benchmark (Batch Mode)", file=sys.stderr)
    print("=" * 70, file=sys.stderr)
    print(f"Configurations: {len(configs)}", file=sys.stderr)
    print(f"Alpha: {args.alpha}", file=sys.stderr)
    print(f"Precision: {args.precision}", file=sys.stderr)
    print(f"Output dir: {args.output_dir}", file=sys.stderr)
    print(f"Python: {python_exe}", file=sys.stderr)
    print("=" * 70, file=sys.stderr)
    
    total_start = time.time()
    all_summaries = []
    
    for i, cfg in enumerate(configs, 1):
        n = cfg["n"]
        R_assoc = cfg["R_assoc"]
        R_null = cfg["R_null"]
        
        print(f"\n{'='*70}", file=sys.stderr)
        print(f"Config {i}/{len(configs)}: n={n}, R_assoc={R_assoc}, R_null={R_null}", file=sys.stderr)
        print(f"{'='*70}", file=sys.stderr)
        
        start_time = time.time()
        
        all_runs, summary = run_benchmark(
            n=n,
            R_assoc=R_assoc,
            R_null=R_null,
            python_exe=python_exe,
            precision=args.precision,
            alpha=args.alpha
        )
        
        # Save per-run results if requested
        if args.save_per_run and all_runs:
            per_run_path = os.path.join(args.output_dir, f"global_lla_n{n}.csv")
            with open(per_run_path, 'w', newline='', encoding='utf-8') as f:
                fieldnames = ["mode", "n", "seed", "P", "LA", "detected"]
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(all_runs)
            print(f"\nPer-run results saved to: {per_run_path}", file=sys.stderr)
        
        # Save individual summary
        summary_path = os.path.join(args.output_dir, f"global_lla_n{n}_sum.csv")
        with open(summary_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=summary.keys())
            writer.writeheader()
            writer.writerow(summary)
        
        all_summaries.append(summary)
        
        # Print results
        print("\n" + "-" * 70, file=sys.stderr)
        print(f"RESULTS for n={n}", file=sys.stderr)
        print("-" * 70, file=sys.stderr)
        print(f"Power (associated detected): {summary['Power']:.3f} ± {summary['Power_SE']:.3f}", 
              file=sys.stderr)
        print(f"FPR (null detected):         {summary['FPR']:.3f} ± {summary['FPR_SE']:.3f}", 
              file=sys.stderr)
        print(f"Expected FPR (alpha):        {args.alpha:.3f}", file=sys.stderr)
        print(f"\nMean LA (associated):        {summary['LA_mean_assoc']:.4f}", file=sys.stderr)
        print(f"Mean LA (null):              {summary['LA_mean_null']:.4f}", file=sys.stderr)
        print(f"\nSummary saved to: {summary_path}", file=sys.stderr)
        print(f"Time for n={n}: {time.time() - start_time:.1f} seconds", file=sys.stderr)
    
    # Save combined summary
    combined_path = os.path.join(args.output_dir, "global_lla_all_summary.csv")
    with open(combined_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=all_summaries[0].keys())
        writer.writeheader()
        writer.writerows(all_summaries)
    
    # Final summary
    print("\n" + "=" * 70, file=sys.stderr)
    print("ALL CONFIGURATIONS COMPLETE", file=sys.stderr)
    print("=" * 70, file=sys.stderr)
    print(f"Combined summary saved to: {combined_path}", file=sys.stderr)
    print(f"Total elapsed time: {time.time() - total_start:.1f} seconds", file=sys.stderr)
    print("=" * 70, file=sys.stderr)


if __name__ == "__main__":
    main()

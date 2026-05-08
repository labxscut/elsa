#!/usr/bin/env python
"""
Comprehensive benchmark for LLA detection, localization, and delay estimation.

Tests:
1. Detection Power: Can LLA distinguish associated vs null triplets?
2. Localization Accuracy: Can LLA identify the correct regulation window?
3. Delay Estimation: Can LLA recover the true Y-Z delay?
4. Robustness: How does performance vary with α, n, window size, delay?
5. P-value method comparison: Compare theo vs perm on the same synthetic triplets.
"""

import numpy as np
import pandas as pd
import subprocess
import sys
import os
import tempfile
import argparse
import time
from pathlib import Path
from typing import List, Dict, Tuple
import json

# Import the unified generator
try:
    from gen_unified_triplets import generate_unified_triplet, write_triplet_to_file
except ImportError:
    print("Error: Cannot import gen_unified_triplets.py", file=sys.stderr)
    print("Make sure gen_unified_triplets.py is in the same directory", file=sys.stderr)
    sys.exit(1)


def build_seed(
    base_seed: int,
    n: int,
    effect_index: int,
    delay_index: int,
    window_index: int,
    replicate_index: int,
    is_control: bool,
) -> int:
    """Build a deterministic seed for a benchmark case."""

    seed = (
        int(base_seed)
        + int(n) * 1_000_000
        + int(effect_index) * 10_000
        + int(delay_index) * 1_000
        + int(window_index) * 100
        + int(replicate_index)
    )
    if is_control:
        seed += 50_000_000
    return seed


def run_lla_analysis(
    input_file: str,
    n: int,
    pvalue_method: str = 'perm',
    delay_limit: int = 0,
    precision: int = 1000,
    keep_trace: bool = True
) -> pd.DataFrame:
    """Run lla_compute.py and parse results."""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tmp_out:
        output_file = tmp_out.name
    
    try:
        cmd = [
            sys.executable, '-m', 'lla.lla_compute',
            input_file, output_file,
            '-s', str(n),
            '-r', '1',
            '-d', str(delay_limit),
            '-p', pvalue_method,
            '-x', str(precision),
            '-n', 'pnz',
            '-f', 'linear'
        ]
        
        if keep_trace:
            cmd.append('--keep-trace')

        env = os.environ.copy()
        repo_root = str(Path(__file__).resolve().parent)
        env['PYTHONPATH'] = os.pathsep.join(
            [repo_root, env.get('PYTHONPATH', '')] if env.get('PYTHONPATH') else [repo_root]
        )
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,
            cwd=repo_root,
            env=env,
        )
        
        if result.returncode != 0:
            print(f"LLA command failed: {result.stderr}", file=sys.stderr)
            return None
        
        # Read with header=0 to use first line as column names
        df = pd.read_csv(output_file, sep=r'\s+', engine='python', header=0)
        return df
        
    except subprocess.TimeoutExpired:
        print("LLA analysis timeout", file=sys.stderr)
        return None
    except Exception as e:
        print(f"Error running LLA: {e}", file=sys.stderr)
        return None
    finally:
        if os.path.exists(output_file):
            os.remove(output_file)


def evaluate_detection(
    df: pd.DataFrame,
    alpha_threshold: float = 0.05,
    is_associated: bool = True
) -> Dict:
    """Evaluate detection performance."""
    
    if df is None or len(df) == 0:
        return {
            'detected': False,
            'p_value': 1.0,
            'lla_score': 0.0
        }
    
    # Filter for the S1-S2-S3 triplet (ignore header row if it exists)
    # The actual triplet should have X/Y/Z columns containing 'S1', 'S2', 'S3'
    triplet_row = df[(df['X'] == 'S1') & (df['Y'] == 'S2') & (df['Z'] == 'S3')]
    
    if len(triplet_row) == 0:
        # No valid triplet found
        return {
            'detected': False,
            'p_value': 1.0,
            'lla_score': 0.0
        }
    
    # Take the first (should be only) row
    best_row = triplet_row.iloc[0]
    
    # Handle both LLA and LA column names (version compatibility)
    la_col = 'LLA' if 'LLA' in df.columns else 'LA'
    la_score = float(best_row[la_col]) if la_col in df.columns else 0.0
    p_value = float(best_row['P']) if 'P' in df.columns else 1.0
    
    detected = p_value < alpha_threshold
    
    return {
        'detected': detected,
        'p_value': p_value,
        'lla_score': la_score
    }


def evaluate_localization(
    df: pd.DataFrame,
    true_window_start: int,
    true_window_end: int,
    alpha_threshold: float = 0.05
) -> Dict:
    """Evaluate localization accuracy."""
    
    if df is None or len(df) == 0:
        return {
            'localization_correct': False,
            'start_error': None,
            'end_error': None,
            'window_overlap': 0.0
        }
    
    # Filter for the S1-S2-S3 triplet
    triplet_row = df[(df['X'] == 'S1') & (df['Y'] == 'S2') & (df['Z'] == 'S3')]
    
    if len(triplet_row) == 0:
        return {
            'localization_correct': False,
            'start_error': None,
            'end_error': None,
            'window_overlap': 0.0
        }
    
    best = triplet_row.iloc[0]
    p_value = float(best['P']) if 'P' in df.columns else 1.0
    
    # Check if detected
    if p_value >= alpha_threshold:
        return {
            'localization_correct': False,
            'start_error': None,
            'end_error': None,
            'window_overlap': 0.0
        }
    
    # Extract detected window (use Start_Z and End_Z as they define regulation region)
    detected_start = int(best['Start_Z']) if 'Start_Z' in df.columns else -1
    detected_end = int(best['End_Z']) if 'End_Z' in df.columns else -1
    
    # Handle -1 (no trace)
    if detected_start == -1 or detected_end == -1:
        return {
            'localization_correct': False,
            'start_error': None,
            'end_error': None,
            'window_overlap': 0.0
        }
    
    # Calculate errors (1-based indices in trace)
    start_error = detected_start - (true_window_start + 1)
    end_error = detected_end - true_window_end
    
    # Calculate overlap (IoU - Intersection over Union)
    true_set = set(range(true_window_start, true_window_end))
    detected_set = set(range(detected_start - 1, detected_end))  # Convert to 0-based
    
    intersection = len(true_set & detected_set)
    union = len(true_set | detected_set)
    overlap = intersection / union if union > 0 else 0.0
    
    return {
        'localization_correct': (start_error == 0 and end_error == 0),
        'start_error': start_error,
        'end_error': end_error,
        'window_overlap': overlap
    }


def evaluate_delay(
    df: pd.DataFrame,
    true_delay: int,
    alpha_threshold: float = 0.05
) -> Dict:
    """Evaluate delay estimation accuracy."""
    
    if df is None or len(df) == 0:
        return {
            'delay_correct': False,
            'detected_delay': None,
            'delay_error': None
        }
    
    # Filter for the S1-S2-S3 triplet
    triplet_row = df[(df['X'] == 'S1') & (df['Y'] == 'S2') & (df['Z'] == 'S3')]
    
    if len(triplet_row) == 0:
        return {
            'delay_correct': False,
            'detected_delay': None,
            'delay_error': None
        }
    
    best = triplet_row.iloc[0]
    p_value = float(best['P']) if 'P' in df.columns else 1.0
    
    # Check if detected
    if p_value >= alpha_threshold:
        return {
            'delay_correct': False,
            'detected_delay': None,
            'delay_error': None
        }
    
    detected_delay = int(best['Delay']) if 'Delay' in df.columns else 0
    delay_error = detected_delay - true_delay
    
    return {
        'delay_correct': (delay_error == 0),
        'detected_delay': detected_delay,
        'delay_error': delay_error
    }


def run_single_experiment(
    n: int,
    effect_size: float,
    window_start: int,
    window_end: int,
    delay_yz: int,
    is_control: bool,
    seed: int,
    delay_limit: int,
    precision: int,
    method: str = 'additive',
    noise_sd: float = 0.1,
    sigma_y: float = 1.0,
    z_encoding: str = '01',
    pvalue_method: str = 'perm'
) -> Tuple[Dict, float]:
    """Run a single experiment and evaluate all metrics.
    
    Returns:
        Tuple of (result_dict, computation_time_seconds)
    """
    
    # Generate triplet
    X, Y, Z, metadata = generate_unified_triplet(
        n=n,
        alpha=1.0,  # Not used for additive method
        method=method,
        effect_size=effect_size,
        noise_sd=noise_sd,
        sigma_y=sigma_y,
        z_encoding=z_encoding,
        window_start=window_start,
        window_end=window_end,
        delay_yz=delay_yz,
        is_control=is_control,
        seed=seed
    )
    
    # Write to temp file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tmp_in:
        input_file = tmp_in.name
    
    try:
        write_triplet_to_file(input_file, X, Y, Z)
        
        # Run LLA analysis and time it
        start_time = time.time()
        df = run_lla_analysis(
            input_file=input_file,
            n=n,
            pvalue_method=pvalue_method,
            delay_limit=delay_limit,
            precision=precision,
            keep_trace=True
        )
        computation_time = time.time() - start_time
        
        # Evaluate detection
        detection_metrics = evaluate_detection(df, is_associated=(not is_control))
        
        # Evaluate localization (only for experimental group)
        # Use window values from metadata (which has actual values, not None)
        if not is_control and df is not None:
            localization_metrics = evaluate_localization(
                df,
                true_window_start=metadata['window_start'],
                true_window_end=metadata['window_end']
            )
        else:
            localization_metrics = {
                'localization_correct': None,
                'start_error': None,
                'end_error': None,
                'window_overlap': None
            }
        
        # Evaluate delay (only if delay_yz > 0 and experimental)
        if not is_control and delay_yz != 0 and df is not None:
            delay_metrics = evaluate_delay(df, true_delay=delay_yz)
        else:
            delay_metrics = {
                'delay_correct': None,
                'detected_delay': None,
                'delay_error': None
            }
        
        # Combine all metrics
        result = {
            **metadata,
            'pvalue_method': pvalue_method,
            **detection_metrics,
            **localization_metrics,
            **delay_metrics
        }
        
        return result, computation_time
        
    finally:
        if os.path.exists(input_file):
            os.remove(input_file)


def run_benchmark_suite(
    n_values: List[int],
    effect_size_values: List[float],
    delay_values: List[int],
    window_fractions: List[float],  # Window size as fraction of n
    n_replicates: int,
    delay_limit: int,
    precision: int,
    output_dir: str,
    method: str = 'additive',
    noise_sd: float = 0.1,
    sigma_y: float = 1.0,
    z_encoding: str = '01',
    pvalue_methods: List[str] = None,
    base_seed: int = 12345,
) -> pd.DataFrame:
    """Run comprehensive benchmark suite."""

    if pvalue_methods is None:
        pvalue_methods = ['perm']
    pvalue_methods = list(dict.fromkeys(pvalue_methods))
    
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    results = []
    total_experiments = 0
    
    # Count total experiments
    for n in n_values:
        for effect_size in effect_size_values:
            for delay in delay_values:
                for frac in window_fractions:
                    # Experimental group + control group, for every p-value method.
                    total_experiments += n_replicates * 2 * len(pvalue_methods)
    
    print(f"Running {total_experiments} experiments...", file=sys.stderr)
    completed = 0
    
    # Track batch timing per condition
    batch_timings = []
    
    for n_index, n in enumerate(n_values):
        for effect_index, effect_size in enumerate(effect_size_values):
            for delay_index, delay_yz in enumerate(delay_values):
                for window_index, window_frac in enumerate(window_fractions):
                    # Define window
                    # Note: argparse may parse "1" as int, so check both 1.0 and 1
                    if window_frac == 1.0 or window_frac == 1:
                        # Global case: pass None to use entire sequence
                        window_start = None
                        window_end = None
                        is_global = True
                    else:
                        # Local case: centered window
                        window_length = max(10, int(n * window_frac))
                        window_start = (n - window_length) // 2
                        window_end = window_start + window_length
                        is_global = False
                    for pvalue_method in pvalue_methods:
                        # Track times for this batch and method
                        batch_times = []
                        batch_start = time.time()

                        # Run experimental replicates
                        for rep in range(n_replicates):
                            seed = build_seed(
                                base_seed=base_seed,
                                n=n,
                                effect_index=effect_index,
                                delay_index=delay_index,
                                window_index=window_index,
                                replicate_index=rep,
                                is_control=False,
                            )

                            result, comp_time = run_single_experiment(
                                n=n,
                                effect_size=effect_size,
                                window_start=window_start,
                                window_end=window_end,
                                delay_yz=delay_yz,
                                is_control=False,
                                seed=seed,
                                delay_limit=delay_limit,
                                precision=precision,
                                method=method,
                                noise_sd=noise_sd,
                                sigma_y=sigma_y,
                                z_encoding=z_encoding,
                                pvalue_method=pvalue_method,
                            )
                            result['replicate'] = rep
                            result['computation_time'] = comp_time
                            results.append(result)
                            batch_times.append(comp_time)

                            completed += 1
                            if completed % 10 == 0:
                                print(f"Progress: {completed}/{total_experiments}", file=sys.stderr)

                        # Run control replicates (for both global and local cases)
                        # Control: same Z structure, but X and Y are independent noise
                        for rep in range(n_replicates):
                            seed = build_seed(
                                base_seed=base_seed,
                                n=n,
                                effect_index=effect_index,
                                delay_index=delay_index,
                                window_index=window_index,
                                replicate_index=rep,
                                is_control=True,
                            )

                            result, comp_time = run_single_experiment(
                                n=n,
                                effect_size=effect_size,
                                window_start=window_start,
                                window_end=window_end,
                                delay_yz=delay_yz,
                                is_control=True,
                                seed=seed,
                                delay_limit=delay_limit,
                                precision=precision,
                                method=method,
                                noise_sd=noise_sd,
                                sigma_y=sigma_y,
                                z_encoding=z_encoding,
                                pvalue_method=pvalue_method,
                            )
                            result['replicate'] = rep
                            result['computation_time'] = comp_time
                            results.append(result)
                            batch_times.append(comp_time)

                            completed += 1
                            if completed % 10 == 0:
                                print(f"Progress: {completed}/{total_experiments}", file=sys.stderr)

                        # Record batch timing statistics per method.
                        batch_total = time.time() - batch_start
                        batch_timings.append({
                            'n': n,
                            'effect_size': effect_size,
                            'delay_yz': delay_yz,
                            'is_global': is_global,
                            'window_fraction': window_frac,
                            'pvalue_method': pvalue_method,
                            'batch_time_seconds': batch_total,
                            'n_runs': n_replicates * 2,
                        })

                        print(f"Batch (n={n}, E={effect_size}, global={is_global}, method={pvalue_method}): "
                              f"{batch_total:.1f}s total, "
                              f"{np.median(batch_times):.2f}s median/run", 
                              file=sys.stderr)
    
    # Convert to DataFrame
    df = pd.DataFrame(results)
    
    # Add batch timing information to results
    # Create a mapping from condition to batch time
    batch_timing_map = {}
    for bt in batch_timings:
        key = (bt['n'], bt['effect_size'], bt['delay_yz'], bt['is_global'], bt['pvalue_method'])
        batch_timing_map[key] = bt['batch_time_seconds']
    
    # Add batch_time column to each row
    df['batch_time_seconds'] = df.apply(
        lambda row: batch_timing_map.get(
            (row['n'], row['effect_size'], row['delay_yz'], row['is_global'], row['pvalue_method']), 
            None
        ), 
        axis=1
    )
    
    # Save raw results
    raw_output = os.path.join(output_dir, 'benchmark_raw_results.csv')
    df.to_csv(raw_output, index=False)
    print(f"Raw results saved to: {raw_output}", file=sys.stderr)
    
    return df, batch_timings


def summarize_results(df: pd.DataFrame, batch_timings: List[Dict], output_dir: str):
    """Generate summary statistics."""

    if df is None or df.empty:
        print("Error: DataFrame is empty. No results to summarize.", file=sys.stderr)
        return

    required_columns = {'n', 'effect_size', 'delay_yz', 'is_global', 'is_control', 'pvalue_method'}
    if not required_columns.issubset(df.columns):
        missing = required_columns - set(df.columns)
        print(f"Error: Missing required columns in DataFrame: {missing}", file=sys.stderr)
        return

    # Detection power by condition
    detection_summary = df.groupby(['n', 'effect_size', 'delay_yz', 'is_global', 'pvalue_method', 'is_control']).agg({
        'detected': ['mean', 'std', 'count'],
        'lla_score': ['mean', 'std']
    }).round(4)
    
    # Flatten column names for easier manipulation
    detection_summary.columns = ['_'.join(col).strip() for col in detection_summary.columns.values]
    detection_summary = detection_summary.reset_index()
    
    # FDR calculation for simulation data
    # Calculate FDR for each condition
    # FDR = FP / (TP + FP) where:
    # - TP (True Positives) = detections in experimental group (is_control=False)
    # - FP (False Positives) = detections in control group (is_control=True)
    
    # Separate experimental and control for Recall calculation
    exp_detections = detection_summary[detection_summary['is_control'] == False].copy()
    ctrl_detections = detection_summary[detection_summary['is_control'] == True].copy()
    
    # Prepare experimental detections with Recall (TPR)
    exp_detections = exp_detections.rename(columns={
        'detected_mean': 'Recall',
        'detected_std': 'Recall_std',
        'detected_count': 'n_exp',
        'lla_score_mean': 'lla_score_mean',
        'lla_score_std': 'lla_score_std'
    })
    
    ctrl_detections = ctrl_detections.rename(columns={
        'detected_mean': 'FPR',
        'detected_std': 'FPR_std',
        'detected_count': 'n_ctrl'
    })
    
    # Merge on condition keys
    merge_keys = ['n', 'effect_size', 'delay_yz', 'is_global', 'pvalue_method']
    detection_summary_final = exp_detections[merge_keys + ['Recall', 'Recall_std', 'n_exp', 'lla_score_mean', 'lla_score_std']].merge(
        ctrl_detections[merge_keys + ['FPR', 'FPR_std', 'n_ctrl']],
        on=merge_keys,
        how='left'
    )
    
    # FPR is already calculated as the detection rate in control group
    # FPR = (# false positives) / (# true negatives + # false positives)
    # In our case: FPR = detected_mean in control group = FPR column
    
    # Note: FDR would be calculated as:
    # FDR = FP / (TP + FP) where:
    # FP = FPR * n_ctrl (number of false positives)
    # TP = Recall * n_exp (number of true positives)
    # But we use FPR here as requested by user
    
    # Add batch timing information
    batch_timing_map = {}
    for bt in batch_timings:
        key = (bt['n'], bt['effect_size'], bt['delay_yz'], bt['is_global'], bt['pvalue_method'])
        batch_timing_map[key] = bt['batch_time_seconds']
    
    detection_summary_final['batch_time_seconds'] = detection_summary_final.apply(
        lambda row: batch_timing_map.get(
            (row['n'], row['effect_size'], row['delay_yz'], row['is_global'], row['pvalue_method']), 
            None
        ), 
        axis=1
    )
    
    # Round metrics
    detection_summary_final['FPR'] = detection_summary_final['FPR'].round(4)
    detection_summary_final['Recall'] = detection_summary_final['Recall'].round(4)
    detection_summary_final['batch_time_seconds'] = detection_summary_final['batch_time_seconds'].round(2)

    detection_output = os.path.join(output_dir, 'detection_summary.csv')
    detection_summary_final.to_csv(detection_output, index=False)
    print(f"Detection summary saved to: {detection_output}", file=sys.stderr)

    # Localization accuracy (experimental only)
    exp_only = df[~df['is_control']]
    if len(exp_only) > 0:
        localization_summary = exp_only.groupby(['n', 'effect_size', 'delay_yz', 'is_global', 'pvalue_method']).agg({
            'window_overlap': ['mean', 'std'],
            'start_error': ['mean', 'std'],
            'end_error': ['mean', 'std']
        }).round(4)

        localization_output = os.path.join(output_dir, 'localization_summary.csv')
        localization_summary.to_csv(localization_output)
        print(f"Localization summary saved to: {localization_output}", file=sys.stderr)

    # Delay estimation accuracy (delayed cases only)
    delayed_exp = exp_only[exp_only['delay_yz'] != 0]
    if len(delayed_exp) > 0:
        delay_summary = delayed_exp.groupby(['n', 'effect_size', 'delay_yz', 'pvalue_method']).agg({
            'delay_correct': ['mean', 'std'],
            'delay_error': ['mean', 'std']
        }).round(4)

        delay_output = os.path.join(output_dir, 'delay_summary.csv')
        delay_summary.to_csv(delay_output)
        print(f"Delay summary saved to: {delay_output}", file=sys.stderr)

    # theo-vs-perm comparison summary (only when both methods are present)
    available_methods = set(str(v) for v in df['pvalue_method'].dropna().unique())
    if {'theo', 'perm'}.issubset(available_methods):
        compare_keys = ['n', 'effect_size', 'delay_yz', 'is_global', 'is_control', 'replicate']
        compare_df = df[df['pvalue_method'].isin(['theo', 'perm'])].copy()
        compare_pivot = compare_df.pivot_table(
            index=compare_keys,
            columns='pvalue_method',
            values=['p_value', 'detected', 'computation_time'],
            aggfunc='first'
        )
        compare_pivot.columns = [f'{left}_{right}' for left, right in compare_pivot.columns.to_flat_index()]
        compare_pivot = compare_pivot.reset_index()

        if {'p_value_theo', 'p_value_perm'}.issubset(compare_pivot.columns):
            compare_pivot['abs_p_diff'] = (compare_pivot['p_value_theo'] - compare_pivot['p_value_perm']).abs()
            compare_pivot['rel_p_diff'] = compare_pivot['abs_p_diff'] / compare_pivot[['p_value_theo', 'p_value_perm']].abs().max(axis=1).clip(lower=1e-12)
        if {'detected_theo', 'detected_perm'}.issubset(compare_pivot.columns):
            compare_pivot['detected_agree'] = (compare_pivot['detected_theo'] == compare_pivot['detected_perm']).astype(float)
        if {'computation_time_theo', 'computation_time_perm'}.issubset(compare_pivot.columns):
            compare_pivot['time_ratio_perm_over_theo'] = compare_pivot['computation_time_perm'] / compare_pivot['computation_time_theo'].clip(lower=1e-12)

        compare_summary = compare_pivot.groupby(['n', 'effect_size', 'delay_yz', 'is_global', 'is_control']).agg({
            'abs_p_diff': ['mean', 'std'],
            'rel_p_diff': ['mean', 'std'],
            'detected_agree': ['mean', 'std'],
            'time_ratio_perm_over_theo': ['mean', 'std'],
        }).round(4)
        compare_output = os.path.join(output_dir, 'pvalue_method_comparison_summary.csv')
        compare_summary.to_csv(compare_output)
        print(f"P-value comparison summary saved to: {compare_output}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description='Comprehensive LLA benchmark')
    
    parser.add_argument('--n_values', nargs='*', type=int, default=[100],
                       help='Sequence lengths to test (default: 100)')
    parser.add_argument('--effect_size_values', nargs='*', type=float, default=[0.8, 1.0],
                       help='Effect sizes to test (default: 0.8 1.0)')
    parser.add_argument('--effect_sweep', action='store_true',
                       help='Run effect_size sweep from 0.8 to 1.0 with step 0.02 (overrides --effect_size_values)')
    parser.add_argument('--effect_start', type=float, default=0.8,
                       help='Effect sweep start value (default: 0.8, used with --effect_sweep)')
    parser.add_argument('--effect_end', type=float, default=1.0,
                       help='Effect sweep end value (default: 1.0, used with --effect_sweep)')
    parser.add_argument('--effect_step', type=float, default=0.02,
                       help='Effect sweep step size (default: 0.02, used with --effect_sweep)')
    parser.add_argument('--method', type=str, default='additive', choices=['mixing', 'additive'],
                       help='Generation method (default: additive)')
    parser.add_argument('--noise_sd', type=float, default=0.1,
                       help='Baseline noise SD (default: 0.1)')
    parser.add_argument('--sigma_y', type=float, default=1.0,
                       help='Y baseline SD (default: 1.0)')
    parser.add_argument('--z_encoding', type=str, default='01', choices=['01', '-11'],
                       help='Z regulator encoding: 01 for 0/1 binary, -11 for -1/1 binary (default: 01)')
    parser.add_argument('--delay_values', nargs='*', type=int, default=[0],
                       help='Y-Z delays to test (default: 0)')
    parser.add_argument('--window_fractions', nargs='*', type=float, default=[0.8],
                       help='Window sizes as fraction of n (default: 0.8; 1.0=global)')
    parser.add_argument('--n_replicates', type=int, default=50,
                       help='Number of replicates per condition (default: 50)')
    parser.add_argument('--pvalue_methods', nargs='+', default=['perm'], choices=['perm', 'theo', 'mix'],
                       help='P-value methods to run (default: perm; use perm theo to compare)')
    parser.add_argument('--base_seed', type=int, default=12345,
                       help='Base seed for deterministic synthetic data generation (default: 12345)')
    parser.add_argument('--delay_limit', type=int, default=0,
                       help='Maximum delay for LLA search (default: 0)')
    parser.add_argument('--precision', type=int, default=1000,
                       help='Permutation precision (default: 1000)')
    parser.add_argument('--output_dir', type=str, default='additive_results',
                       help='Output directory (default: additive_results)')
    
    args = parser.parse_args()
    
    # Handle effect_size sweep mode
    if args.effect_sweep:
        effect_size_values = np.arange(args.effect_start, args.effect_end + args.effect_step/2, args.effect_step)
        effect_size_values = np.round(effect_size_values, 3).tolist()  # Avoid floating point errors
        print(f"Effect size sweep mode: testing {len(effect_size_values)} values from {args.effect_start} to {args.effect_end}", 
              file=sys.stderr)
        print(f"Effect size values: {effect_size_values}", file=sys.stderr)
    else:
        effect_size_values = args.effect_size_values
    
    # Run benchmark
    df, batch_timings = run_benchmark_suite(
        n_values=args.n_values,
        effect_size_values=effect_size_values,
        delay_values=args.delay_values,
        window_fractions=args.window_fractions,
        n_replicates=args.n_replicates,
        delay_limit=args.delay_limit,
        precision=args.precision,
        output_dir=args.output_dir,
        method=args.method,
        noise_sd=args.noise_sd,
        sigma_y=args.sigma_y,
        z_encoding=args.z_encoding,
        pvalue_methods=args.pvalue_methods,
        base_seed=args.base_seed,
    )
    
    # Generate summaries
    summarize_results(df, batch_timings, args.output_dir)
    
    print("Benchmark complete!", file=sys.stderr)


if __name__ == '__main__':
    main()

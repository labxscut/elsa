#!/usr/bin/env python
"""Visualize combined effect size sweep results from multiple directories.

This script loads detection_summary.csv and localization_summary.csv from multiple
result directories, combines them, and generates a unified visualization.

Usage:
    python visualize_combined_results.py --dirs additive_results additive_results2 --save combined_plot.png
    
Or specify all directories with same prefix:
    python visualize_combined_results.py --prefix additive_results --save combined_plot.png
"""

import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from glob import glob


def load_detection_summary(directory: str) -> pd.DataFrame:
    """Load detection summary from a single directory."""
    det_path = os.path.join(directory, 'detection_summary.csv')
    if not os.path.exists(det_path):
        raise FileNotFoundError(f"Missing detection_summary.csv in {directory}")
    return pd.read_csv(det_path)


def load_localization_summary(directory: str) -> pd.DataFrame:
    """Load localization summary from a single directory."""
    loc_path = os.path.join(directory, 'localization_summary.csv')
    if not os.path.exists(loc_path):
        raise FileNotFoundError(f"Missing localization_summary.csv in {directory}")
    
    # Read the raw file to parse the 3-row header
    with open(loc_path, 'r') as f:
        lines = [f.readline().strip() for _ in range(3)]
    
    # Parse the three header rows
    row0 = lines[0].split(',')  # Metric names
    row1 = lines[1].split(',')  # Statistics (mean/std)
    row2 = lines[2].split(',')  # Index names
    
    # Build proper column names
    columns = []
    for i, (r0, r1, r2) in enumerate(zip(row0, row1, row2)):
        if r2:  # Has index name - this is an index column
            columns.append(r2)
        elif r0 and r1:  # Has metric and stat - this is a data column
            columns.append(f"{r0}_{r1}")
        else:
            columns.append(f"col_{i}")
    
    # Read data with skiprows and custom column names
    loc = pd.read_csv(loc_path, skiprows=3, names=columns)
    return loc


def combine_results(directories: list) -> tuple:
    """Combine results from multiple directories (or load single directory).
    
    Returns:
        (combined_detection_df, combined_localization_df)
    """
    all_detection = []
    all_localization = []
    
    for directory in directories:
        print(f"Loading data from: {directory}")
        try:
            det = load_detection_summary(directory)
            loc = load_localization_summary(directory)
            
            # Add source directory for tracking
            det['source_dir'] = directory
            loc['source_dir'] = directory
            
            all_detection.append(det)
            all_localization.append(loc)
            
            print(f"  Loaded {len(det)} detection rows, {len(loc)} localization rows")
            
        except FileNotFoundError as e:
            print(f"  Warning: {e}")
            continue
    
    if not all_detection:
        raise ValueError("No valid data found in any directory")
    
    # Combine all dataframes
    combined_det = pd.concat(all_detection, ignore_index=True)
    combined_loc = pd.concat(all_localization, ignore_index=True)
    
    # Remove duplicates (keep first occurrence) based on key columns
    key_cols = ['n', 'effect_size', 'delay_yz', 'is_global']
    combined_det = combined_det.drop_duplicates(subset=key_cols, keep='first')
    combined_loc = combined_loc.drop_duplicates(subset=key_cols, keep='first')
    
    # Sort by effect_size
    combined_det = combined_det.sort_values('effect_size').reset_index(drop=True)
    combined_loc = combined_loc.sort_values('effect_size').reset_index(drop=True)
    
    if len(directories) == 1:
        print(f"\nLoaded results from single directory:")
    else:
        print(f"\nCombined results from {len(directories)} directories:")
    print(f"  Detection: {len(combined_det)} rows")
    print(f"  Localization: {len(combined_loc)} rows")
    print(f"  Effect sizes: {sorted(combined_det['effect_size'].unique())}")
    
    return combined_det, combined_loc


def compute_derived(merged: pd.DataFrame, window_fraction: float = 0.8) -> pd.DataFrame:
    """Compute derived metrics."""
    m = merged.copy()
    true_len = window_fraction * m['n']
    m['detected_len_est'] = true_len + (-m['start_error_mean'] + m['end_error_mean'])
    m['inflation_factor'] = m['detected_len_est'] / true_len
    return m


def plot_combined_effect_sweep(det: pd.DataFrame, loc: pd.DataFrame, merged: pd.DataFrame, 
                                derived: pd.DataFrame, window_fraction: float = 0.8, 
                                n: int = 100, save_path: str = None):
    """Plot combined effect size sweep results."""
    effect_sizes = det['effect_size']
    fig, axes = plt.subplots(3, 2, figsize=(14, 12))
    
    # Get number of data points for title
    n_points = len(effect_sizes)
    effect_range = f"{effect_sizes.min():.2f}-{effect_sizes.max():.2f}"
    
    fig.suptitle(f'Combined Effect Size Sweep (n={n_points} points, E={effect_range})\n'
                 f'Additive Method: n={n}, window_fraction={window_fraction}, noise_sd=0.1',
                 fontsize=12)

    # 1. LLA score
    ax = axes[0,0]
    ax.errorbar(effect_sizes, det['lla_score_mean'], yerr=det['lla_score_std'], 
                fmt='o-', capsize=4, linewidth=2, markersize=6, color='tab:blue', 
                label='Mean ± 1 SD (n=50)')
    ax.set_xlabel('Effect Size')
    ax.set_ylabel('LLA Score')
    ax.set_title('Association Strength vs Effect Size')
    ax.grid(alpha=0.3)
    ax.legend(loc='best', fontsize=9)
    if det['lla_score_mean'].min() < 0:
        ax.axhline(0, color='gray', linestyle=':', linewidth=1, alpha=0.5)

    # 2. Window overlap
    ax = axes[0,1]
    ax.errorbar(effect_sizes, loc['window_overlap_mean'], yerr=loc['window_overlap_std'], 
                fmt='s--', capsize=4, color='tab:green', label='Mean ± 1 SD (n=50)')
    ax.axhline(0.9, color='darkblue', linestyle=':', linewidth=1, alpha=0.6, label='Excellent (0.9)')
    ax.set_xlabel('Effect Size')
    ax.set_ylabel('Window Overlap (IoU)')
    ax.set_ylim(0, 1.0)
    ax.set_title('Localization Overlap vs Effect Size')
    ax.grid(alpha=0.3)
    ax.legend(loc='best', fontsize=8)

    # 3. Boundary errors
    ax = axes[1,0]
    ax.plot(effect_sizes, loc['start_error_mean'], 'o-', label='Start Error', 
            color='tab:orange', linewidth=2)
    ax.plot(effect_sizes, loc['end_error_mean'], 'o-', label='End Error', 
            color='tab:purple', linewidth=2)
    ax.fill_between(effect_sizes, 
                     loc['start_error_mean']-loc['start_error_std'], 
                     loc['start_error_mean']+loc['start_error_std'], 
                     color='tab:orange', alpha=0.2)
    ax.fill_between(effect_sizes, 
                     loc['end_error_mean']-loc['end_error_std'], 
                     loc['end_error_mean']+loc['end_error_std'], 
                     color='tab:purple', alpha=0.2)
    ax.axhline(0, color='black', linewidth=1, linestyle='--', alpha=0.7, label='Perfect (0)')
    ax.set_xlabel('Effect Size')
    ax.set_ylabel('Boundary Error (indices)')
    ax.set_title('Boundary Localization Accuracy vs Effect Size')
    ax.grid(alpha=0.3)
    ax.legend(loc='upper right', fontsize=8, ncol=2)

    # 4. FPR (False Positive Rate)
    ax = axes[1,1]
    if 'FPR' in det.columns:
        fpr_values = det['FPR']
        y_label = 'False Positive Rate (FPR)'
        title = 'False Positive Control vs Effect Size'
    else:
        fpr_values = pd.Series([0]*len(effect_sizes))
        y_label = 'False Positive Rate (FPR, unavailable)'
        title = 'FPR vs Effect Size (data unavailable)'
    
    ax.bar(effect_sizes, fpr_values, width=0.012, color='tab:red', alpha=0.7, label='Observed FPR')
    ax.axhline(0.05, color='darkred', linestyle='--', linewidth=1.5, alpha=0.7, label='Nominal α=0.05')
    
    ax.set_xlabel('Effect Size')
    ax.set_ylabel(y_label)
    y_max = max(0.15, max(fpr_values)*1.5 + 0.02) if max(fpr_values) > 0 else 0.15
    ax.set_ylim(0, y_max)
    ax.set_title(title)
    ax.grid(axis='y', alpha=0.3)
    ax.legend(loc='best', fontsize=9)

    # 5. Batch time
    ax = axes[2,0]
    ax.plot(effect_sizes, det['batch_time_seconds'], 'd-', color='tab:blue', 
            linewidth=2, markersize=7)
    ax.set_xlabel('Effect Size')
    ax.set_ylabel('Batch Time (s)')
    ax.set_title('Computational Efficiency vs Effect Size\n(Time for 100 experiments)')
    ax.grid(alpha=0.3)
    ax.set_ylim(bottom=0)

    # 6. Inflation factor
    ax = axes[2,1]
    # 注释掉绿色阴影
    # ax.axhspan(1.0, 1.1, color='lightgreen', alpha=0.2, label='Excellent (<10% inflation)')
    ax.plot(effect_sizes, derived['inflation_factor'], 'o-', color='tab:brown', 
            linewidth=2, markersize=6, label='Detected / True Length')
    # 注释掉Perfect虚线
    # ax.axhline(1.0, color='darkgreen', linestyle='--', linewidth=1.5, alpha=0.8, label='Perfect (1.0)')
    ax.set_xlabel('Effect Size')
    ax.set_ylabel('Detected Length / True Length')
    ax.set_ylim(1.0, 2.0)
    ax.set_title('Window Size Accuracy vs Effect Size')
    ax.grid(alpha=0.3)
    # 仅保留折线的图例
    ax.legend(loc='best', fontsize=8)
    plt.tight_layout(rect=[0,0,1,0.97])
    if save_path:
        fig.savefig(save_path, dpi=200)
        print(f"\nSaved combined figure to: {save_path}")
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser(
        description='Visualize effect size sweep results from one or multiple directories',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Single directory
    python visualize_combined_results.py --dirs additive_results --save plot.png
    
    # Combine two specific directories
    python visualize_combined_results.py --dirs additive_results additive_results2 --save combined.png
    
    # Auto-discover all directories with prefix
    python visualize_combined_results.py --prefix additive_results --save combined.png
    
    # Display without saving
    python visualize_combined_results.py --dirs additive_results
        """
    )
    
    parser.add_argument('--dirs', nargs='+', type=str, default=None,
                       help='List of result directories to combine')
    parser.add_argument('--prefix', type=str, default=None,
                       help='Auto-discover directories starting with this prefix')
    parser.add_argument('--save', type=str, default=None, 
                       help='Path to save output figure')
    parser.add_argument('--window_fraction', type=float, default=0.8,
                       help='Window fraction used in benchmark (default: 0.8)')
    parser.add_argument('--n', type=int, default=100,
                       help='Sequence length used in benchmark (default: 100)')
    
    args = parser.parse_args()
    
    # Determine directories to process
    directories = []
    if args.dirs:
        directories = args.dirs
    elif args.prefix:
        # Find all directories starting with prefix
        pattern = f"{args.prefix}*"
        found = [d for d in glob(pattern) if os.path.isdir(d)]
        if not found:
            print(f"Error: No directories found matching pattern '{pattern}'")
            return
        directories = sorted(found)
        print(f"Auto-discovered {len(directories)} directories with prefix '{args.prefix}':")
        for d in directories:
            print(f"  - {d}")
    else:
        print("Error: Must specify either --dirs or --prefix")
        parser.print_help()
        return
    
    # Validate directories exist
    for d in directories:
        if not os.path.isdir(d):
            print(f"Warning: Directory does not exist: {d}")
    
    directories = [d for d in directories if os.path.isdir(d)]
    if not directories:
        print("Error: No valid directories found")
        return
    
    print(f"\nCombining results from {len(directories)} directories...")
    
    # Load and combine data
    combined_det, combined_loc = combine_results(directories)
    
    # Merge detection and localization
    merged = pd.merge(combined_det, combined_loc, 
                     on=['n', 'effect_size', 'delay_yz', 'is_global'], 
                     suffixes=('_det', '_loc'))
    
    # Compute derived metrics
    derived = compute_derived(merged, window_fraction=args.window_fraction)
    
    # Plot
    plot_combined_effect_sweep(combined_det, combined_loc, merged, derived,
                               window_fraction=args.window_fraction,
                               n=args.n,
                               save_path=args.save)


if __name__ == '__main__':
    main()

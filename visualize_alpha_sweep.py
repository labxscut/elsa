#!/usr/bin/env python
"""Visualize alpha sweep benchmark results.

Generates a multi-panel figure summarizing:
1. LLA score vs alpha (mean ± std)
2. Window overlap (IoU) vs alpha (mean ± std)
3. Start and end boundary errors vs alpha (mean ± std shaded)
4. False Positive Rate vs alpha (bar + count annotation)
5. Batch computation time vs alpha
6. Aggregated error inflation factor vs alpha

Usage (from repo root):
    python visualize_alpha_sweep.py --dir alpha_win8 --save alpha_win8/alpha_sweep_summary.png

Requires: pandas, matplotlib (present in requirements.txt)
"""

import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def load_data(directory: str):
    det_path = os.path.join(directory, 'detection_summary.csv')
    loc_path = os.path.join(directory, 'localization_summary.csv')
    if not os.path.exists(det_path) or not os.path.exists(loc_path):
        raise FileNotFoundError(f"Missing required CSV files in {directory}")
    
    # Load detection summary (simple header)
    det = pd.read_csv(det_path)
    
    # Load localization summary
    # This CSV has a 3-row header from pandas groupby.to_csv()
    # Row 0: metric names (window_overlap, window_overlap, start_error, ...)
    # Row 1: statistic names (mean, std, mean, std, ...)
    # Row 2: index column names (n, alpha, delay_yz, is_global, ...)
    
    # Read the raw file to parse manually
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
    
    # Now read the data with skiprows and custom column names
    loc = pd.read_csv(loc_path, skiprows=3, names=columns)
    
    print(f"Loaded localization summary with columns: {list(loc.columns)}")
    
    # Ensure sorted by alpha
    det = det.sort_values('alpha').reset_index(drop=True)
    loc = loc.sort_values('alpha').reset_index(drop=True)
    
    # Merge on shared keys
    merged = pd.merge(det, loc, on=['n','alpha','delay_yz','is_global'], suffixes=('_det','_loc'))
    return det, loc, merged


def compute_derived(merged: pd.DataFrame) -> pd.DataFrame:
    m = merged.copy()
    # True window length from generation: window_fraction * n. We infer from overlap stats (given) or from run parameters.
    # Benchmark used window_fraction=0.8 with n=100 -> true length=80.
    true_len = 0.8 * m['n']
    # Detected length approximation: true length / overlap * intersection length.
    # We do not have intersection length directly; overshoot approximation via start/end errors: detected_len = true_len + |start_error_mean| + |end_error_mean|.
    m['detected_len_est'] = true_len + (-m['start_error_mean'] + m['end_error_mean'])  # start_error negative; end_error positive.
    m['inflation_factor'] = m['detected_len_est'] / true_len
    return m


def plot_alpha_sweep(det: pd.DataFrame, loc: pd.DataFrame, merged: pd.DataFrame, derived: pd.DataFrame, save_path: str=None):
    alphas = det['alpha']
    fig, axes = plt.subplots(3, 2, figsize=(14, 12))
    fig.suptitle('Alpha Sweep Benchmark (n=100, window_fraction=0.8, delay=0)\nFPR-based Analysis')

    # 1. LLA score
    ax = axes[0,0]
    ax.errorbar(alphas, det['lla_score_mean'], yerr=det['lla_score_std'], fmt='o-', capsize=4, 
                linewidth=2, markersize=6, color='tab:blue', label='Mean ± 1 SD (n=50)')
    ax.set_xlabel('Alpha')
    ax.set_ylabel('LLA Score')
    ax.set_title('Association Strength vs Alpha')
    ax.grid(alpha=0.3)
    ax.legend(loc='best', fontsize=9)
    # Add zero reference line if scores can be negative
    if det['lla_score_mean'].min() < 0:
        ax.axhline(0, color='gray', linestyle=':', linewidth=1, alpha=0.5)

    # 2. Window overlap
    ax = axes[0,1]
    ax.errorbar(alphas, loc['window_overlap_mean'], yerr=loc['window_overlap_std'], fmt='s--', capsize=4, color='tab:green', label='Mean ± 1 SD (n=50)')
    # Add reference lines for context
    # ax.axhline(0.7, color='orange', linestyle=':', linewidth=1, alpha=0.6, label='Good (0.7)')
    ax.axhline(0.9, color='darkblue', linestyle=':', linewidth=1, alpha=0.6, label='Excellent (0.9)')
    ax.set_xlabel('Alpha')
    ax.set_ylabel('Window Overlap (IoU)')
    ax.set_ylim(0.5, 1.0)  # Changed from 0.8 to 0.5 to show full quality spectrum
    ax.set_title('Localization Overlap vs Alpha')
    ax.grid(alpha=0.3)
    ax.legend(loc='best', fontsize=8)

    # 3. Boundary errors
    ax = axes[1,0]
    # Add acceptable error zone (±2 indices is quite good)
    ax.plot(alphas, loc['start_error_mean'], 'o-', label='Start Error', color='tab:orange', linewidth=2)
    ax.plot(alphas, loc['end_error_mean'], 'o-', label='End Error', color='tab:purple', linewidth=2)
    ax.fill_between(alphas, loc['start_error_mean']-loc['start_error_std'], loc['start_error_mean']+loc['start_error_std'], color='tab:orange', alpha=0.2)
    ax.fill_between(alphas, loc['end_error_mean']-loc['end_error_std'], loc['end_error_mean']+loc['end_error_std'], color='tab:purple', alpha=0.2)
    ax.axhline(0, color='black', linewidth=1, linestyle='--', alpha=0.7, label='Perfect (0)')
    ax.set_xlabel('Alpha')
    ax.set_ylabel('Boundary Error (indices)')
    ax.set_title('Boundary Localization Accuracy vs Alpha\n(Errors within ±2 indices = Excellent)')
    ax.grid(alpha=0.3)
    ax.legend(loc='upper right', fontsize=8, ncol=2)

    # 4. FPR (False Positive Rate)
    ax = axes[1,1]
    # Use FPR column from detection_summary.csv
    if 'FPR' in det.columns:
        fpr_values = det['FPR']
        y_label = 'False Positive Rate (FPR)'
        title = 'False Positive Control vs Alpha'
    else:
        # Fallback if FPR not available
        fpr_values = pd.Series([0]*len(alphas))
        y_label = 'False Positive Rate (FPR, unavailable)'
        title = 'FPR vs Alpha (data unavailable)'
    
    ax.bar(alphas, fpr_values, width=0.012, color='tab:red', alpha=0.7, label='Observed FPR')
    # Add reference line for nominal α=0.05 level
    ax.axhline(0.05, color='darkred', linestyle='--', linewidth=1.5, alpha=0.7, label='Nominal α=0.05')
    
    # Annotate with actual counts if possible
    if 'n_ctrl' in det.columns:
        for a, fpr, nctrl in zip(alphas, fpr_values, det['n_ctrl']):
            fp_count = int(round(fpr * nctrl))
            ax.text(a, fpr + 0.005, f'{fp_count}', ha='center', va='bottom', fontsize=7, fontweight='bold')
    
    ax.set_xlabel('Alpha')
    ax.set_ylabel(y_label)
    # Set reasonable y-axis: max 0.15 or 1.5x max FPR, whichever is larger
    y_max = max(0.15, max(fpr_values)*1.5 + 0.02) if max(fpr_values) > 0 else 0.15
    ax.set_ylim(0, y_max)
    ax.set_title(title)
    ax.grid(axis='y', alpha=0.3)
    ax.legend(loc='best', fontsize=9)

    # 5. Batch time
    ax = axes[2,0]
    ax.plot(alphas, det['batch_time_seconds'], 'd-', color='tab:blue', linewidth=2, markersize=7)
    ax.set_xlabel('Alpha')
    ax.set_ylabel('Batch Time (s)')
    ax.set_title('Computational Efficiency vs Alpha\n(Time for 100 experiments)')
    ax.grid(alpha=0.3)
    # Start y-axis from 0 for better visual scaling
    ax.set_ylim(bottom=0)

    # 6. Inflation factor
    ax = axes[2,1]
    # Add acceptable zone (inflation < 1.1 is excellent)
    ax.axhspan(1.0, 1.1, color='lightgreen', alpha=0.2, label='Excellent (<10% inflation)')
    ax.plot(alphas, derived['inflation_factor'], 'o-', color='tab:brown', linewidth=2, markersize=6, label='Detected / True Length')
    ax.axhline(1.0, color='darkgreen', linestyle='--', linewidth=1.5, alpha=0.8, label='Perfect (1.0)')
    ax.set_xlabel('Alpha')
    ax.set_ylabel('Detected Length / True Length')
    ax.set_ylim(1.0, 2.0)  # Fixed range to avoid exaggerating small inflations
    ax.set_title('Window Size Accuracy vs Alpha\n(Values near 1.0 = Accurate localization)')
    ax.grid(alpha=0.3)
    ax.legend(loc='best', fontsize=8)

    plt.tight_layout(rect=[0,0,1,0.97])
    if save_path:
        fig.savefig(save_path, dpi=200)
        print(f"Saved figure to {save_path}")
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser(description='Visualize alpha sweep benchmark results')
    parser.add_argument('--dir', type=str, default='alpha_win8', help='Directory containing CSV summaries')
    parser.add_argument('--save', type=str, default=None, help='Optional path to save figure')
    args = parser.parse_args()

    det, loc, merged = load_data(args.dir)
    derived = compute_derived(merged)
    plot_alpha_sweep(det, loc, merged, derived, save_path=args.save)


if __name__ == '__main__':
    main()

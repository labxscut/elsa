#!/usr/bin/env python
"""
Visualize LLA Benchmark Results

This script reads CSV results from benchmark_lla_unified.py and generates comprehensive plots:
1. Recall and FPR vs alpha (detection performance)
2. Localization IoU vs alpha
3. Start error and End error vs alpha (with scatter plot to check correlation)
4. Combined multi-panel figures for publication
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import sys
import os
from pathlib import Path
from typing import Optional


def load_results(results_dir: str) -> tuple:
    """Load benchmark results from CSV files.
    
    Args:
        results_dir: Directory containing benchmark results
        
    Returns:
        Tuple of (raw_df, detection_summary_df, localization_summary_df)
    """
    raw_path = os.path.join(results_dir, 'benchmark_raw_results.csv')
    detection_path = os.path.join(results_dir, 'detection_summary.csv')
    localization_path = os.path.join(results_dir, 'localization_summary.csv')
    
    if not os.path.exists(raw_path):
        raise FileNotFoundError(f"Raw results not found: {raw_path}")
    
    raw_df = pd.read_csv(raw_path)
    
    detection_df = None
    if os.path.exists(detection_path):
        detection_df = pd.read_csv(detection_path)
    
    localization_df = None
    if os.path.exists(localization_path):
        localization_df = pd.read_csv(localization_path)
    
    return raw_df, detection_df, localization_df


def plot_detection_performance(detection_df: pd.DataFrame, output_dir: str, 
                               filter_params: Optional[dict] = None):
    """Plot Recall and FPR vs alpha.
    
    Args:
        detection_df: Detection summary DataFrame
        output_dir: Output directory for plots
        filter_params: Optional dict to filter data (e.g., {'n': 40, 'is_global': True})
    """
    if detection_df is None:
        print("No detection summary data available", file=sys.stderr)
        return
    
    # Apply filters
    df = detection_df.copy()
    if filter_params:
        for key, value in filter_params.items():
            if key in df.columns:
                df = df[df[key] == value]
    
    if len(df) == 0:
        print(f"No data after filtering with {filter_params}", file=sys.stderr)
        return
    
    # Sort by alpha
    df = df.sort_values('alpha')
    
    # Create two-panel plot
    fig, axes = plt.subplots(2, 1, figsize=(10, 10))
    
    alpha_values = df['alpha'].values
    
    # Plot 1: Recall vs Alpha
    ax1 = axes[0]
    ax1.plot(alpha_values, df['Recall'], 'b-o', linewidth=2, markersize=6, label='Recall (TPR)')
    if 'Recall_std' in df.columns:
        ax1.fill_between(alpha_values,
                          df['Recall'] - df['Recall_std'],
                          df['Recall'] + df['Recall_std'],
                          alpha=0.2, color='blue')
    ax1.set_xlabel('Association Strength (α)', fontsize=12)
    ax1.set_ylabel('Recall (True Positive Rate)', fontsize=12)
    ax1.set_title('LLA Detection Recall vs Association Strength', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(-0.05, 1.05)
    ax1.axhline(y=0.8, color='gray', linestyle='--', alpha=0.5, label='80% threshold')
    ax1.legend()
    
    # Plot 2: FPR vs Alpha
    ax2 = axes[1]
    ax2.plot(alpha_values, df['FPR'], 'r-s', linewidth=2, markersize=6, label='FPR')
    if 'FPR_std' in df.columns:
        ax2.fill_between(alpha_values,
                          df['FPR'] - df['FPR_std'],
                          df['FPR'] + df['FPR_std'],
                          alpha=0.2, color='red')
    ax2.set_xlabel('Association Strength (α)', fontsize=12)
    ax2.set_ylabel('False Positive Rate', fontsize=12)
    ax2.set_title('LLA False Positive Rate vs Association Strength', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(-0.05, 0.25)
    ax2.axhline(y=0.05, color='gray', linestyle='--', alpha=0.5, label='5% threshold')
    ax2.legend()
    
    plt.tight_layout()
    
    # Save
    suffix = '_'.join([f"{k}{v}" for k, v in (filter_params or {}).items()])
    filename = f'detection_vs_alpha_{suffix}.png' if suffix else 'detection_vs_alpha.png'
    output_path = os.path.join(output_dir, filename)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Detection plot saved: {output_path}", file=sys.stderr)
    
    # Combined plot
    fig2, ax = plt.subplots(1, 1, figsize=(10, 6))
    ax.plot(alpha_values, df['Recall'], 'b-o', linewidth=2, markersize=6, label='Recall (TPR)')
    ax.plot(alpha_values, df['FPR'], 'r-s', linewidth=2, markersize=6, label='FPR')
    
    if 'Recall_std' in df.columns:
        ax.fill_between(alpha_values,
                         df['Recall'] - df['Recall_std'],
                         df['Recall'] + df['Recall_std'],
                         alpha=0.2, color='blue')
    if 'FPR_std' in df.columns:
        ax.fill_between(alpha_values,
                         df['FPR'] - df['FPR_std'],
                         df['FPR'] + df['FPR_std'],
                         alpha=0.2, color='red')
    
    ax.set_xlabel('Association Strength (α)', fontsize=12)
    ax.set_ylabel('Rate', fontsize=12)
    ax.set_title('LLA Detection Performance vs Association Strength', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-0.05, 1.05)
    ax.axhline(y=0.05, color='gray', linestyle='--', alpha=0.5, linewidth=1)
    ax.legend(fontsize=11, loc='right')
    
    plt.tight_layout()
    
    filename_combined = f'detection_combined_{suffix}.png' if suffix else 'detection_combined.png'
    output_path_combined = os.path.join(output_dir, filename_combined)
    plt.savefig(output_path_combined, dpi=300, bbox_inches='tight')
    print(f"Combined detection plot saved: {output_path_combined}", file=sys.stderr)
    
    plt.close('all')


def plot_localization_performance(localization_df: pd.DataFrame, output_dir: str,
                                  filter_params: Optional[dict] = None):
    """Plot localization metrics vs alpha.
    
    Args:
        localization_df: Localization summary DataFrame
        output_dir: Output directory for plots
        filter_params: Optional dict to filter data
    """
    if localization_df is None:
        print("No localization summary data available", file=sys.stderr)
        return
    
    # Reset index to access multi-level columns
    df = localization_df.reset_index()
    
    # Apply filters
    if filter_params:
        for key, value in filter_params.items():
            if key in df.columns:
                df = df[df[key] == value]
    
    if len(df) == 0:
        print(f"No localization data after filtering with {filter_params}", file=sys.stderr)
        return
    
    # Sort by alpha
    df = df.sort_values('alpha')
    alpha_values = df['alpha'].values
    
    # Create three-panel plot
    fig, axes = plt.subplots(3, 1, figsize=(10, 14))
    
    # Plot 1: IoU (window overlap) vs Alpha
    ax1 = axes[0]
    if ('window_overlap', 'mean') in localization_df.columns:
        iou_mean = df[('window_overlap', 'mean')]
        iou_std = df[('window_overlap', 'std')]
    else:
        # Try flattened column names
        iou_mean = df['window_overlap_mean'] if 'window_overlap_mean' in df.columns else None
        iou_std = df['window_overlap_std'] if 'window_overlap_std' in df.columns else None
    
    if iou_mean is not None:
        ax1.plot(alpha_values, iou_mean, 'g-o', linewidth=2, markersize=6, label='IoU')
        if iou_std is not None:
            ax1.fill_between(alpha_values,
                              iou_mean - iou_std,
                              iou_mean + iou_std,
                              alpha=0.2, color='green')
        ax1.set_ylabel('Intersection over Union', fontsize=12)
        ax1.set_title('Localization Accuracy (IoU) vs Association Strength', fontsize=14, fontweight='bold')
        ax1.grid(True, alpha=0.3)
        ax1.set_ylim(-0.05, 1.05)
        ax1.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5, label='50% threshold')
        ax1.legend()
    
    # Plot 2: Start Error vs Alpha
    ax2 = axes[1]
    if ('start_error', 'mean') in localization_df.columns:
        start_mean = df[('start_error', 'mean')]
        start_std = df[('start_error', 'std')]
    else:
        start_mean = df['start_error_mean'] if 'start_error_mean' in df.columns else None
        start_std = df['start_error_std'] if 'start_error_std' in df.columns else None
    
    if start_mean is not None:
        ax2.plot(alpha_values, start_mean, 'b-s', linewidth=2, markersize=6, label='Start Error')
        if start_std is not None:
            ax2.fill_between(alpha_values,
                              start_mean - start_std,
                              start_mean + start_std,
                              alpha=0.2, color='blue')
        ax2.set_ylabel('Start Position Error (time steps)', fontsize=12)
        ax2.set_title('Window Start Detection Error vs Association Strength', fontsize=14, fontweight='bold')
        ax2.grid(True, alpha=0.3)
        ax2.axhline(y=0, color='gray', linestyle='--', alpha=0.5, label='Perfect')
        ax2.legend()
    
    # Plot 3: End Error vs Alpha
    ax3 = axes[2]
    if ('end_error', 'mean') in localization_df.columns:
        end_mean = df[('end_error', 'mean')]
        end_std = df[('end_error', 'std')]
    else:
        end_mean = df['end_error_mean'] if 'end_error_mean' in df.columns else None
        end_std = df['end_error_std'] if 'end_error_std' in df.columns else None
    
    if end_mean is not None:
        ax3.plot(alpha_values, end_mean, 'r-^', linewidth=2, markersize=6, label='End Error')
        if end_std is not None:
            ax3.fill_between(alpha_values,
                              end_mean - end_std,
                              end_mean + end_std,
                              alpha=0.2, color='red')
        ax3.set_xlabel('Association Strength (α)', fontsize=12)
        ax3.set_ylabel('End Position Error (time steps)', fontsize=12)
        ax3.set_title('Window End Detection Error vs Association Strength', fontsize=14, fontweight='bold')
        ax3.grid(True, alpha=0.3)
        ax3.axhline(y=0, color='gray', linestyle='--', alpha=0.5, label='Perfect')
        ax3.legend()
    
    plt.tight_layout()
    
    # Save
    suffix = '_'.join([f"{k}{v}" for k, v in (filter_params or {}).items()])
    filename = f'localization_vs_alpha_{suffix}.png' if suffix else 'localization_vs_alpha.png'
    output_path = os.path.join(output_dir, filename)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Localization plot saved: {output_path}", file=sys.stderr)
    
    plt.close('all')


def plot_error_correlation(raw_df: pd.DataFrame, output_dir: str,
                           filter_params: Optional[dict] = None):
    """Plot scatter plot of start_error vs end_error to check correlation.
    
    Args:
        raw_df: Raw results DataFrame
        output_dir: Output directory for plots
        filter_params: Optional dict to filter data
    """
    # Filter experimental group only
    df = raw_df[~raw_df['is_control']].copy()
    
    # Apply additional filters
    if filter_params:
        for key, value in filter_params.items():
            if key in df.columns:
                df = df[df[key] == value]
    
    if len(df) == 0:
        print(f"No data for error correlation after filtering", file=sys.stderr)
        return
    
    # Remove NaN values
    df = df.dropna(subset=['start_error', 'end_error'])
    
    if len(df) == 0:
        print("No valid error data for correlation plot", file=sys.stderr)
        return
    
    # Group by alpha for color coding
    alphas = sorted(df['alpha'].unique())
    
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    
    # Create colormap
    cmap = plt.cm.viridis
    colors = [cmap(i / len(alphas)) for i in range(len(alphas))]
    
    for i, alpha in enumerate(alphas):
        alpha_data = df[df['alpha'] == alpha]
        ax.scatter(alpha_data['start_error'], alpha_data['end_error'],
                  c=[colors[i]], label=f'α={alpha:.2f}', alpha=0.6, s=50)
    
    # Add diagonal line (perfect correlation)
    all_errors = pd.concat([df['start_error'], df['end_error']])
    min_err, max_err = all_errors.min(), all_errors.max()
    ax.plot([min_err, max_err], [min_err, max_err], 'k--', alpha=0.3, label='y=x')
    
    # Calculate and display correlation
    corr = df['start_error'].corr(df['end_error'])
    ax.text(0.05, 0.95, f'Correlation: {corr:.3f}', 
            transform=ax.transAxes, fontsize=12,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    ax.set_xlabel('Start Error (time steps)', fontsize=12)
    ax.set_ylabel('End Error (time steps)', fontsize=12)
    ax.set_title('Correlation between Start and End Detection Errors', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
    ax.axhline(y=0, color='gray', linestyle='-', alpha=0.3, linewidth=0.5)
    ax.axvline(x=0, color='gray', linestyle='-', alpha=0.3, linewidth=0.5)
    
    plt.tight_layout()
    
    # Save
    suffix = '_'.join([f"{k}{v}" for k, v in (filter_params or {}).items()])
    filename = f'error_correlation_{suffix}.png' if suffix else 'error_correlation.png'
    output_path = os.path.join(output_dir, filename)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Error correlation plot saved: {output_path}", file=sys.stderr)
    
    plt.close('all')


def create_comprehensive_figure(raw_df: pd.DataFrame, detection_df: pd.DataFrame,
                                localization_df: pd.DataFrame, output_dir: str,
                                filter_params: Optional[dict] = None):
    """Create a comprehensive multi-panel figure combining all metrics.
    
    Args:
        raw_df: Raw results DataFrame
        detection_df: Detection summary DataFrame
        localization_df: Localization summary DataFrame
        output_dir: Output directory
        filter_params: Optional dict to filter data
    """
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.3)
    
    # Filter detection data
    det_df = detection_df.copy() if detection_df is not None else None
    if det_df is not None and filter_params:
        for key, value in filter_params.items():
            if key in det_df.columns:
                det_df = det_df[det_df[key] == value]
        det_df = det_df.sort_values('alpha')
    
    # Filter localization data
    loc_df = localization_df.reset_index() if localization_df is not None else None
    if loc_df is not None and filter_params:
        for key, value in filter_params.items():
            if key in loc_df.columns:
                loc_df = loc_df[loc_df[key] == value]
        loc_df = loc_df.sort_values('alpha')
    
    # Panel 1: Recall vs Alpha
    ax1 = fig.add_subplot(gs[0, 0])
    if det_df is not None and len(det_df) > 0:
        alpha_vals = det_df['alpha'].values
        ax1.plot(alpha_vals, det_df['Recall'], 'b-o', linewidth=2, markersize=5)
        if 'Recall_std' in det_df.columns:
            ax1.fill_between(alpha_vals, det_df['Recall'] - det_df['Recall_std'],
                            det_df['Recall'] + det_df['Recall_std'], alpha=0.2, color='blue')
        ax1.set_ylabel('Recall (TPR)', fontsize=11)
        ax1.set_title('Detection Recall', fontsize=12, fontweight='bold')
        ax1.grid(True, alpha=0.3)
        ax1.set_ylim(-0.05, 1.05)
    
    # Panel 2: FPR vs Alpha
    ax2 = fig.add_subplot(gs[0, 1])
    if det_df is not None and len(det_df) > 0:
        ax2.plot(alpha_vals, det_df['FPR'], 'r-s', linewidth=2, markersize=5)
        if 'FPR_std' in det_df.columns:
            ax2.fill_between(alpha_vals, det_df['FPR'] - det_df['FPR_std'],
                            det_df['FPR'] + det_df['FPR_std'], alpha=0.2, color='red')
        ax2.set_ylabel('False Positive Rate', fontsize=11)
        ax2.set_title('False Positive Rate', fontsize=12, fontweight='bold')
        ax2.grid(True, alpha=0.3)
        ax2.axhline(y=0.05, color='gray', linestyle='--', alpha=0.5)
    
    # Panel 3: IoU vs Alpha
    ax3 = fig.add_subplot(gs[1, 0])
    if loc_df is not None and len(loc_df) > 0:
        alpha_vals_loc = loc_df['alpha'].values
        iou_mean = loc_df.get(('window_overlap', 'mean'), loc_df.get('window_overlap_mean'))
        iou_std = loc_df.get(('window_overlap', 'std'), loc_df.get('window_overlap_std'))
        if iou_mean is not None:
            ax3.plot(alpha_vals_loc, iou_mean, 'g-o', linewidth=2, markersize=5)
            if iou_std is not None:
                ax3.fill_between(alpha_vals_loc, iou_mean - iou_std,
                                iou_mean + iou_std, alpha=0.2, color='green')
        ax3.set_ylabel('IoU (Overlap)', fontsize=11)
        ax3.set_title('Localization Accuracy', fontsize=12, fontweight='bold')
        ax3.grid(True, alpha=0.3)
        ax3.set_ylim(-0.05, 1.05)
    
    # Panel 4: Start & End Errors vs Alpha
    ax4 = fig.add_subplot(gs[1, 1])
    if loc_df is not None and len(loc_df) > 0:
        start_mean = loc_df.get(('start_error', 'mean'), loc_df.get('start_error_mean'))
        end_mean = loc_df.get(('end_error', 'mean'), loc_df.get('end_error_mean'))
        if start_mean is not None:
            ax4.plot(alpha_vals_loc, start_mean, 'b-s', linewidth=2, markersize=5, label='Start Error')
        if end_mean is not None:
            ax4.plot(alpha_vals_loc, end_mean, 'r-^', linewidth=2, markersize=5, label='End Error')
        ax4.set_ylabel('Position Error (steps)', fontsize=11)
        ax4.set_title('Window Boundary Errors', fontsize=12, fontweight='bold')
        ax4.grid(True, alpha=0.3)
        ax4.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
        ax4.legend(fontsize=10)
    
    # Panel 5: Error Correlation Scatter
    ax5 = fig.add_subplot(gs[2, :])
    exp_df = raw_df[~raw_df['is_control']].copy()
    if filter_params:
        for key, value in filter_params.items():
            if key in exp_df.columns:
                exp_df = exp_df[exp_df[key] == value]
    exp_df = exp_df.dropna(subset=['start_error', 'end_error'])
    
    if len(exp_df) > 0:
        alphas = sorted(exp_df['alpha'].unique())
        cmap = plt.cm.viridis
        colors = [cmap(i / len(alphas)) for i in range(len(alphas))]
        
        for i, alpha in enumerate(alphas):
            alpha_data = exp_df[exp_df['alpha'] == alpha]
            ax5.scatter(alpha_data['start_error'], alpha_data['end_error'],
                       c=[colors[i]], label=f'α={alpha:.2f}', alpha=0.6, s=30)
        
        all_errors = pd.concat([exp_df['start_error'], exp_df['end_error']])
        if len(all_errors) > 0:
            min_err, max_err = all_errors.min(), all_errors.max()
            ax5.plot([min_err, max_err], [min_err, max_err], 'k--', alpha=0.3, linewidth=1)
        
        corr = exp_df['start_error'].corr(exp_df['end_error'])
        ax5.text(0.02, 0.98, f'Corr: {corr:.3f}', transform=ax5.transAxes,
                fontsize=11, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        ax5.set_xlabel('Start Error (time steps)', fontsize=11)
        ax5.set_ylabel('End Error (time steps)', fontsize=11)
        ax5.set_title('Start vs End Error Correlation', fontsize=12, fontweight='bold')
        ax5.grid(True, alpha=0.3)
        ax5.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8, ncol=2)
    
    # Add overall title
    title_parts = []
    if filter_params:
        for k, v in filter_params.items():
            title_parts.append(f"{k}={v}")
    title = f"LLA Benchmark Results" + (f" ({', '.join(title_parts)})" if title_parts else "")
    fig.suptitle(title, fontsize=16, fontweight='bold', y=0.995)
    
    # Save
    suffix = '_'.join([f"{k}{v}" for k, v in (filter_params or {}).items()])
    filename = f'comprehensive_results_{suffix}.png' if suffix else 'comprehensive_results.png'
    output_path = os.path.join(output_dir, filename)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Comprehensive figure saved: {output_path}", file=sys.stderr)
    
    plt.close('all')


def main():
    parser = argparse.ArgumentParser(description='Visualize LLA Benchmark Results')
    
    parser.add_argument('results_dir', type=str,
                       help='Directory containing benchmark results CSV files')
    parser.add_argument('--output_dir', type=str, default=None,
                       help='Output directory for plots (default: same as results_dir)')
    parser.add_argument('--n', type=int, default=None,
                       help='Filter by sequence length n')
    parser.add_argument('--delay', type=int, default=None,
                       help='Filter by delay value')
    parser.add_argument('--global_only', action='store_true',
                       help='Show only global window results')
    parser.add_argument('--local_only', action='store_true',
                       help='Show only local window results')
    parser.add_argument('--no_comprehensive', action='store_true',
                       help='Skip comprehensive multi-panel figure')
    
    args = parser.parse_args()
    
    # Set output directory
    output_dir = args.output_dir if args.output_dir else args.results_dir
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    print(f"Loading results from: {args.results_dir}", file=sys.stderr)
    
    # Load data
    try:
        raw_df, detection_df, localization_df = load_results(args.results_dir)
        print(f"Loaded {len(raw_df)} raw results", file=sys.stderr)
    except Exception as e:
        print(f"Error loading results: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Build filter parameters
    filter_params = {}
    if args.n is not None:
        filter_params['n'] = args.n
    if args.delay is not None:
        filter_params['delay_yz'] = args.delay
    if args.global_only:
        filter_params['is_global'] = True
    elif args.local_only:
        filter_params['is_global'] = False
    
    print(f"Filters: {filter_params if filter_params else 'None'}", file=sys.stderr)
    
    # Generate plots
    print("\nGenerating plots...", file=sys.stderr)
    
    if detection_df is not None:
        plot_detection_performance(detection_df, output_dir, filter_params)
    
    if localization_df is not None:
        plot_localization_performance(localization_df, output_dir, filter_params)
    
    plot_error_correlation(raw_df, output_dir, filter_params)
    
    if not args.no_comprehensive:
        create_comprehensive_figure(raw_df, detection_df, localization_df,
                                    output_dir, filter_params)
    
    print("\n" + "="*60, file=sys.stderr)
    print("Visualization complete!", file=sys.stderr)
    print(f"Plots saved to: {output_dir}", file=sys.stderr)
    print("="*60, file=sys.stderr)


if __name__ == '__main__':
    main()

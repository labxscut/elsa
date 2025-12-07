#!/usr/bin/env python
"""
visualize_benchmark_results.py

Create plots from benchmark results to visualize LLA reliability metrics.

Usage:
    python visualize_benchmark_results.py lla_reliability_summary.csv
"""

import sys
import csv
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path


def load_summary(csv_path):
    """Load summary CSV into list of dicts."""
    with open(csv_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        rows = []
        for row in reader:
            # Convert numeric fields
            for key in ['n', 'corr_len', 'R']:
                row[key] = int(row[key])
            for key in ['J_mean', 'J_std', 'Coverage_mean', 'Coverage_std',
                       'MAE_start_mean', 'MAE_start_std', 'MAE_end_mean', 'MAE_end_std',
                       'Power', 'Power_SE', 'LA_mean', 'LA_std']:
                row[key] = float(row[key])
            rows.append(row)
    return rows


def plot_metrics(data, output_prefix="lla_reliability"):
    """Create plots for benchmark metrics."""
    
    # Extract data
    n_vals = [r['n'] for r in data]
    
    j_mean = [r['J_mean'] for r in data]
    j_std = [r['J_std'] for r in data]
    
    cov_mean = [r['Coverage_mean'] for r in data]
    cov_std = [r['Coverage_std'] for r in data]
    
    mae_start_mean = [r['MAE_start_mean'] for r in data]
    mae_start_std = [r['MAE_start_std'] for r in data]
    
    mae_end_mean = [r['MAE_end_mean'] for r in data]
    mae_end_std = [r['MAE_end_std'] for r in data]
    
    power = [r['Power'] for r in data]
    power_se = [r['Power_SE'] for r in data]
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('LLA Reliability Metrics (No-Delay Case)', fontsize=16, fontweight='bold')
    
    # Plot 1: Jaccard and Coverage
    ax = axes[0, 0]
    ax.errorbar(n_vals, j_mean, yerr=j_std, marker='o', label='Jaccard', capsize=5, linewidth=2)
    ax.errorbar(n_vals, cov_mean, yerr=cov_std, marker='s', label='Coverage', capsize=5, linewidth=2)
    ax.axhline(y=0.7, color='gray', linestyle='--', alpha=0.5, label='Good threshold (0.7)')
    ax.set_xlabel('Time points (n)', fontsize=12)
    ax.set_ylabel('Metric value', fontsize=12)
    ax.set_title('Spatial Accuracy', fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)
    ax.set_ylim([0, 1.05])
    
    # Plot 2: Boundary Errors
    ax = axes[0, 1]
    ax.errorbar(n_vals, mae_start_mean, yerr=mae_start_std, marker='o', label='MAE Start', capsize=5, linewidth=2)
    ax.errorbar(n_vals, mae_end_mean, yerr=mae_end_std, marker='s', label='MAE End', capsize=5, linewidth=2)
    ax.axhline(y=2.0, color='gray', linestyle='--', alpha=0.5, label='Good threshold (≤2)')
    ax.set_xlabel('Time points (n)', fontsize=12)
    ax.set_ylabel('Mean Absolute Error (time points)', fontsize=12)
    ax.set_title('Boundary Precision', fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)
    
    # Plot 3: Statistical Power
    ax = axes[1, 0]
    ax.errorbar(n_vals, power, yerr=power_se, marker='o', color='darkgreen', capsize=5, linewidth=2, markersize=8)
    ax.axhline(y=0.8, color='gray', linestyle='--', alpha=0.5, label='Good threshold (0.8)')
    ax.axhline(y=0.05, color='red', linestyle=':', alpha=0.5, label='Type I error (α=0.05)')
    ax.set_xlabel('Time points (n)', fontsize=12)
    ax.set_ylabel('Power (detection rate)', fontsize=12)
    ax.set_title('Statistical Power @ α=0.05', fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)
    ax.set_ylim([0, 1.05])
    
    # Plot 4: Summary table
    ax = axes[1, 1]
    ax.axis('off')
    
    # Create table data
    table_data = [['n', 'Jaccard', 'Coverage', 'MAE', 'Power']]
    for i, row in enumerate(data):
        table_data.append([
            f"{row['n']}",
            f"{row['J_mean']:.3f}±{row['J_std']:.3f}",
            f"{row['Coverage_mean']:.3f}±{row['Coverage_std']:.3f}",
            f"{row['MAE_start_mean']:.1f}/{row['MAE_end_mean']:.1f}",
            f"{row['Power']:.3f}±{row['Power_SE']:.3f}"
        ])
    
    table = ax.table(cellText=table_data, cellLoc='center', loc='center',
                    colWidths=[0.12, 0.22, 0.22, 0.22, 0.22])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2)
    
    # Style header row
    for i in range(len(table_data[0])):
        table[(0, i)].set_facecolor('#4CAF50')
        table[(0, i)].set_text_props(weight='bold', color='white')
    
    # Style data rows
    for i in range(1, len(table_data)):
        for j in range(len(table_data[0])):
            table[(i, j)].set_facecolor('#f0f0f0' if i % 2 == 0 else 'white')
    
    ax.set_title('Summary Statistics', fontsize=14, fontweight='bold', pad=20)
    
    # Adjust layout and save
    plt.tight_layout(rect=[0, 0.03, 1, 0.97])
    
    output_path = f"{output_prefix}_plots.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved plot to: {output_path}")
    
    plt.show()


def main():
    if len(sys.argv) < 2:
        print("Usage: python visualize_benchmark_results.py <summary_csv>", file=sys.stderr)
        print("\nExample:", file=sys.stderr)
        print("  python visualize_benchmark_results.py lla_reliability_summary.csv", file=sys.stderr)
        sys.exit(1)
    
    csv_path = sys.argv[1]
    
    if not Path(csv_path).exists():
        print(f"Error: File not found: {csv_path}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Loading data from: {csv_path}")
    data = load_summary(csv_path)
    
    if not data:
        print("Error: No data loaded", file=sys.stderr)
        sys.exit(1)
    
    print(f"Loaded {len(data)} configurations")
    print("Creating plots...")
    
    output_prefix = Path(csv_path).stem
    plot_metrics(data, output_prefix=output_prefix)
    
    print("Done!")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""
visualize_benchmark_advanced.py

Create advanced visualizations for benchmark results with varying corr_len.
Includes heat maps, proportion analysis, and multi-panel comparisons.

Usage:
    python visualize_benchmark_advanced.py lla_reliability_summary.csv
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
            # Add derived fields
            row['proportion'] = row['corr_len'] / row['n']
            rows.append(row)
    return rows


def create_heatmap(data, metric_key, metric_name, output_prefix):
    """Create heat map of metric vs (n, corr_len)."""
    # Get unique values
    n_vals = sorted(set(r['n'] for r in data))
    corr_lens = sorted(set(r['corr_len'] for r in data))
    
    # Create matrix
    matrix = np.full((len(n_vals), len(corr_lens)), np.nan)
    
    for i, n in enumerate(n_vals):
        for j, cl in enumerate(corr_lens):
            matches = [r for r in data if r['n'] == n and r['corr_len'] == cl]
            if matches:
                matrix[i, j] = matches[0][metric_key]
    
    # Plot
    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(matrix, aspect='auto', cmap='RdYlGn', vmin=0, vmax=1 if 'Power' in metric_key or 'J' in metric_key or 'Coverage' in metric_key else None)
    
    ax.set_xticks(range(len(corr_lens)))
    ax.set_yticks(range(len(n_vals)))
    ax.set_xticklabels(corr_lens)
    ax.set_yticklabels(n_vals)
    
    ax.set_xlabel('Correlation Length (corr_len)', fontsize=12)
    ax.set_ylabel('Time Points (n)', fontsize=12)
    ax.set_title(f'{metric_name} vs (n, corr_len)', fontsize=14, fontweight='bold')
    
    # Add text annotations
    for i in range(len(n_vals)):
        for j in range(len(corr_lens)):
            if not np.isnan(matrix[i, j]):
                text = ax.text(j, i, f'{matrix[i, j]:.2f}',
                             ha="center", va="center", color="black", fontsize=9)
    
    plt.colorbar(im, ax=ax, label=metric_name)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_heatmap_{metric_key}.png", dpi=300, bbox_inches='tight')
    print(f"Saved heat map: {output_prefix}_heatmap_{metric_key}.png")
    plt.close()


def plot_proportion_analysis(data, output_prefix):
    """Plot metrics vs signal proportion (corr_len/n)."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Performance vs Signal Proportion (corr_len/n)', fontsize=16, fontweight='bold')
    
    # Group by n for color coding
    n_vals = sorted(set(r['n'] for r in data))
    colors = plt.cm.viridis(np.linspace(0, 1, len(n_vals)))
    
    # Plot 1: Jaccard vs proportion
    ax = axes[0, 0]
    for n, color in zip(n_vals, colors):
        subset = [r for r in data if r['n'] == n]
        props = [r['proportion'] for r in subset]
        jaccards = [r['J_mean'] for r in subset]
        ax.scatter(props, jaccards, label=f'n={n}', color=color, s=80, alpha=0.7)
        ax.plot(props, jaccards, color=color, alpha=0.3, linestyle='--')
    ax.set_xlabel('Signal Proportion (corr_len/n)', fontsize=11)
    ax.set_ylabel('Jaccard', fontsize=11)
    ax.set_title('Spatial Accuracy vs Proportion', fontsize=12, fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)
    
    # Plot 2: Coverage vs proportion
    ax = axes[0, 1]
    for n, color in zip(n_vals, colors):
        subset = [r for r in data if r['n'] == n]
        props = [r['proportion'] for r in subset]
        covs = [r['Coverage_mean'] for r in subset]
        ax.scatter(props, covs, label=f'n={n}', color=color, s=80, alpha=0.7)
        ax.plot(props, covs, color=color, alpha=0.3, linestyle='--')
    ax.set_xlabel('Signal Proportion (corr_len/n)', fontsize=11)
    ax.set_ylabel('Coverage', fontsize=11)
    ax.set_title('Coverage vs Proportion', fontsize=12, fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)
    
    # Plot 3: MAE vs proportion
    ax = axes[1, 0]
    for n, color in zip(n_vals, colors):
        subset = [r for r in data if r['n'] == n]
        props = [r['proportion'] for r in subset]
        maes = [(r['MAE_start_mean'] + r['MAE_end_mean'])/2 for r in subset]
        ax.scatter(props, maes, label=f'n={n}', color=color, s=80, alpha=0.7)
        ax.plot(props, maes, color=color, alpha=0.3, linestyle='--')
    ax.set_xlabel('Signal Proportion (corr_len/n)', fontsize=11)
    ax.set_ylabel('Average MAE', fontsize=11)
    ax.set_title('Boundary Precision vs Proportion', fontsize=12, fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)
    
    # Plot 4: Power vs proportion
    # ax = axes[1, 1]
    # for n, color in zip(n_vals, colors):
    #     subset = [r for r in data if r['n'] == n]
    #     props = [r['proportion'] for r in subset]
    #     powers = [r['Power'] for r in subset]
    #     ax.scatter(props, powers, label=f'n={n}', color=color, s=80, alpha=0.7)
    #     ax.plot(props, powers, color=color, alpha=0.3, linestyle='--')
    # ax.axhline(y=0.8, color='red', linestyle='--', alpha=0.5, label='Target (0.8)')
    # ax.set_xlabel('Signal Proportion (corr_len/n)', fontsize=11)
    # ax.set_ylabel('Power @ α=0.05', fontsize=11)
    # ax.set_title('Detection Power vs Proportion', fontsize=12, fontweight='bold')
    # ax.legend()
    # ax.grid(alpha=0.3)
    
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.savefig(f"{output_prefix}_proportion_analysis.png", dpi=300, bbox_inches='tight')
    print(f"Saved proportion analysis: {output_prefix}_proportion_analysis.png")
    plt.close()


def plot_absolute_vs_relative(data, output_prefix):
    """Compare absolute signal length vs relative proportion effects."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle('Absolute vs Relative Signal Effects', fontsize=16, fontweight='bold')
    
    # Plot 1: Power vs absolute corr_len (grouped by n)
    ax = axes[0]
    n_vals = sorted(set(r['n'] for r in data))
    colors = plt.cm.viridis(np.linspace(0, 1, len(n_vals)))
    
    for n, color in zip(n_vals, colors):
        subset = [r for r in data if r['n'] == n]
        corr_lens = [r['corr_len'] for r in subset]
        powers = [r['Power'] for r in subset]
        ax.plot(corr_lens, powers, marker='o', label=f'n={n}', color=color, linewidth=2, markersize=8)
    
    ax.axhline(y=0.8, color='red', linestyle='--', alpha=0.5, label='Target (0.8)')
    ax.set_xlabel('Absolute Signal Length (corr_len)', fontsize=12)
    ax.set_ylabel('Power @ α=0.05', fontsize=12)
    ax.set_title('Effect of Absolute Signal Duration', fontsize=13, fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)
    
    # Plot 2: Power vs proportion (showing same data differently)
    ax = axes[1]
    for n, color in zip(n_vals, colors):
        subset = [r for r in data if r['n'] == n]
        props = [r['proportion'] for r in subset]
        powers = [r['Power'] for r in subset]
        ax.plot(props, powers, marker='s', label=f'n={n}', color=color, linewidth=2, markersize=8)
    
    ax.axhline(y=0.8, color='red', linestyle='--', alpha=0.5, label='Target (0.8)')
    ax.set_xlabel('Relative Signal Proportion (corr_len/n)', fontsize=12)
    ax.set_ylabel('Power @ α=0.05', fontsize=12)
    ax.set_title('Effect of Relative Signal Proportion', fontsize=13, fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)
    
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(f"{output_prefix}_absolute_vs_relative.png", dpi=300, bbox_inches='tight')
    print(f"Saved absolute vs relative: {output_prefix}_absolute_vs_relative.png")
    plt.close()


def main():
    if len(sys.argv) < 2:
        print("Usage: python visualize_benchmark_advanced.py <summary_csv>", file=sys.stderr)
        print("\nExample:", file=sys.stderr)
        print("  python visualize_benchmark_advanced.py lla_reliability_summary.csv", file=sys.stderr)
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
    print(f"  n values: {sorted(set(r['n'] for r in data))}")
    print(f"  corr_len values: {sorted(set(r['corr_len'] for r in data))}")
    
    output_prefix = Path(csv_path).stem
    
    print("\nCreating visualizations...")
    
    # Heat maps for key metrics
    print("\n1. Creating heat maps...")
    create_heatmap(data, 'J_mean', 'Jaccard', output_prefix)
    create_heatmap(data, 'Coverage_mean', 'Coverage', output_prefix)
    create_heatmap(data, 'Power', 'Power @ α=0.05', output_prefix)
    
    # Proportion analysis
    print("\n2. Creating proportion analysis...")
    plot_proportion_analysis(data, output_prefix)
    
    # Absolute vs relative comparison
    print("\n3. Creating absolute vs relative comparison...")
    plot_absolute_vs_relative(data, output_prefix)
    
    print("\n✓ All visualizations created successfully!")
    print(f"\nOutput files:")
    print(f"  - {output_prefix}_heatmap_*.png (3 heat maps)")
    print(f"  - {output_prefix}_proportion_analysis.png")
    print(f"  - {output_prefix}_absolute_vs_relative.png")


if __name__ == "__main__":
    main()

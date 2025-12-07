#!/usr/bin/env python
"""Plot Recall vs Effect Size from effect size sweep results.

This script loads detection_summary.csv from one or multiple result directories
and plots only the Recall curve.

Usage:
    python plot_recall_only.py --dirs additive_results additive_results2 --save recall_plot.png
    python plot_recall_only.py --prefix additive_results --save recall_plot.png
"""

import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob


def load_detection_summary(directory: str) -> pd.DataFrame:
    """Load detection summary from a single directory."""
    det_path = os.path.join(directory, 'detection_summary.csv')
    if not os.path.exists(det_path):
        raise FileNotFoundError(f"Missing detection_summary.csv in {directory}")
    return pd.read_csv(det_path)


def combine_detection_results(directories: list) -> pd.DataFrame:
    """Combine detection summaries from multiple directories."""
    all_detection = []
    for directory in directories:
        print(f"Loading data from: {directory}")
        try:
            det = load_detection_summary(directory)
            det['source_dir'] = directory
            all_detection.append(det)
            print(f"  Loaded {len(det)} detection rows")
        except FileNotFoundError as e:
            print(f"  Warning: {e}")
            continue

    if not all_detection:
        raise ValueError("No valid detection data found in any directory")

    combined_det = pd.concat(all_detection, ignore_index=True)

    # Deduplicate and sort
    key_cols = ['n', 'effect_size', 'delay_yz', 'is_global']
    combined_det = combined_det.drop_duplicates(subset=key_cols, keep='first')
    combined_det = combined_det.sort_values('effect_size').reset_index(drop=True)

    print(f"\nCombined detection results from {len(directories)} directories:")
    print(f"  Total rows: {len(combined_det)}")
    print(f"  Effect sizes: {sorted(combined_det['effect_size'].unique())}")
    return combined_det


def plot_recall(det: pd.DataFrame, n: int = 100, window_fraction: float = 0.8, save_path: str = None):
    """Plot Recall vs Effect Size."""
    effect_sizes = det['effect_size']

    fig, ax = plt.subplots(1, 1, figsize=(7, 5))

    n_points = len(effect_sizes)
    effect_range = f"{effect_sizes.min():.2f}-{effect_sizes.max():.2f}"
    fig.suptitle(f'Recall vs Effect Size (n={n_points} points, E={effect_range})\n'
                 f'Additive Method: n={n}, window_fraction={window_fraction}, noise_sd=0.1',
                 fontsize=12)

    if 'Recall' not in det.columns:
        print("Warning: 'Recall' column not found in detection_summary.csv")
        ax.text(0.5, 0.5, 'Recall\nunavailable', transform=ax.transAxes,
                ha='center', va='center', fontsize=14, color='gray')
        ax.set_xlabel('Effect Size')
        ax.set_ylabel('Recall')
        ax.set_ylim(0, 1)
    else:
        recall_mean = det['Recall']
        has_std = 'Recall_std' in det.columns

        if has_std:
            recall_std = det['Recall_std']
            ax.errorbar(effect_sizes, recall_mean, yerr=recall_std,
                        fmt='o-', capsize=4, linewidth=2, markersize=6,
                        color='tab:olive', label='Mean ± 1 SD (n=50)')
        else:
            ax.plot(effect_sizes, recall_mean, 'o-',
                    linewidth=2, markersize=6, color='tab:olive',
                    label='Recall')

        # Smart y-lim to handle all-ones case
        ymax = 1.05 if recall_mean.max() >= 1.0 else 1.0
        ax.set_ylim(0, ymax)

        ax.set_ylabel('Recall')
        ax.set_xlabel('Effect Size')
        ax.set_title('Recall vs Effect Size')
        ax.grid(alpha=0.3)
        ax.legend(loc='lower right', fontsize=9)

    plt.tight_layout(rect=[0, 0, 1, 0.92])
    if save_path:
        fig.savefig(save_path, dpi=200, bbox_inches='tight')
        print(f"\nSaved Recall plot to: {save_path}")
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser(
        description='Plot Recall vs Effect Size from sweep results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python plot_recall_only.py --dirs additive_results --save recall.png
    python plot_recall_only.py --prefix additive_results --save recall.png
        """
    )
    parser.add_argument('--dirs', nargs='+', type=str, default=None,
                        help='List of result directories')
    parser.add_argument('--prefix', type=str, default=None,
                        help='Auto-discover directories with this prefix')
    parser.add_argument('--save', type=str, default=None,
                        help='Path to save the figure')
    parser.add_argument('--window_fraction', type=float, default=0.8,
                        help='Window fraction used in benchmark (default: 0.8)')
    parser.add_argument('--n', type=int, default=100,
                        help='Sequence length used in benchmark (default: 100)')

    args = parser.parse_args()

    directories = []
    if args.dirs:
        directories = args.dirs
    elif args.prefix:
        pattern = f"{args.prefix}*"
        found = [d for d in glob(pattern) if os.path.isdir(d)]
        if not found:
            print(f"Error: No directories found matching '{pattern}'")
            return
        directories = sorted(found)
        print(f"Auto-discovered {len(directories)} directories:")
        for d in directories:
            print(f"  - {d}")
    else:
        print("Error: Must specify --dirs or --prefix")
        parser.print_help()
        return

    directories = [d for d in directories if os.path.isdir(d)]
    if not directories:
        print("Error: No valid directories found")
        return

    print(f"\nProcessing {len(directories)} directories for Recall plot...")
    combined_det = combine_detection_results(directories)
    plot_recall(combined_det, n=args.n, window_fraction=args.window_fraction, save_path=args.save)


if __name__ == '__main__':
    main()
#!/usr/bin/env python
"""
Quick helper to process benchmark results and generate FDR-based visualizations.

Usage:
    python process_benchmark_results.py --dir alpha_win8

This will:
1. Recalculate detection_summary.csv with FDR from raw results
2. Generate visualization with FDR
3. Create output in {dir}_fdr/
"""

import argparse
import os
import sys
import subprocess


def main():
    parser = argparse.ArgumentParser(description='Process benchmark results with FDR')
    parser.add_argument('--dir', type=str, required=True,
                       help='Directory containing benchmark_raw_results.csv')
    parser.add_argument('--skip-recalc', action='store_true',
                       help='Skip FDR recalculation (use if already done)')
    args = parser.parse_args()
    
    input_dir = args.dir
    output_dir = f"{input_dir}_fdr"
    raw_file = os.path.join(input_dir, 'benchmark_raw_results.csv')
    loc_file = os.path.join(input_dir, 'localization_summary.csv')
    
    # Validate input
    if not os.path.exists(raw_file):
        print(f"❌ Error: {raw_file} not found")
        return 1
    
    if not os.path.exists(loc_file):
        print(f"❌ Error: {loc_file} not found")
        return 1
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Step 1: Recalculate FDR
    if not args.skip_recalc:
        print(f"📊 Step 1: Recalculating FDR from {raw_file}...")
        result = subprocess.run([
            sys.executable, 'recalculate_fdr.py',
            '--raw', raw_file,
            '--output-dir', output_dir
        ])
        if result.returncode != 0:
            print("❌ FDR recalculation failed")
            return 1
        print("✅ FDR recalculation complete\n")
    else:
        print("⏭️  Skipping FDR recalculation\n")
    
    # Step 2: Copy localization summary
    print(f"📋 Step 2: Copying localization summary...")
    import shutil
    shutil.copy(loc_file, os.path.join(output_dir, 'localization_summary.csv'))
    print("✅ Localization summary copied\n")
    
    # Step 3: Generate visualization
    print(f"📈 Step 3: Generating FDR-based visualization...")
    viz_output = os.path.join(output_dir, 'alpha_sweep_FDR_summary.png')
    result = subprocess.run([
        sys.executable, 'visualize_alpha_sweep.py',
        '--dir', output_dir,
        '--save', viz_output
    ])
    if result.returncode != 0:
        print("❌ Visualization failed")
        return 1
    print(f"✅ Visualization saved to {viz_output}\n")
    
    # Summary
    print("=" * 60)
    print("✨ Processing complete!")
    print("=" * 60)
    print(f"Results directory: {output_dir}/")
    print(f"  - detection_summary.csv (with FDR)")
    print(f"  - localization_summary.csv")
    print(f"  - alpha_sweep_FDR_summary.png")
    print()
    print("To view the figure:")
    print(f"  Start-Process {viz_output}")
    print()
    
    return 0


if __name__ == '__main__':
    exit(main())

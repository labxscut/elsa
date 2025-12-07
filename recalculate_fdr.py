#!/usr/bin/env python
"""
Recalculate FDR from existing benchmark_raw_results.csv
This allows updating detection_summary.csv with FDR instead of FPR
without re-running the entire benchmark.
"""

import pandas as pd
import argparse
import os


def recalculate_detection_summary(raw_results_path: str, output_dir: str):
    """Recalculate detection summary with FDR from raw results."""
    
    # Load raw results
    df = pd.read_csv(raw_results_path)
    
    print(f"Loaded {len(df)} raw results from {raw_results_path}")
    
    # Detection power by condition
    detection_summary = df.groupby(['n', 'alpha', 'delay_yz', 'is_global', 'is_control']).agg({
        'detected': ['mean', 'std', 'count'],
        'lla_score': ['mean', 'std']
    }).round(4)
    
    # Flatten column names
    detection_summary.columns = ['_'.join(col).strip() for col in detection_summary.columns.values]
    detection_summary = detection_summary.reset_index()
    
    # Separate experimental and control
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
        'detected_mean': 'FP_rate',
        'detected_std': 'FP_rate_std',
        'detected_count': 'n_ctrl'
    })
    
    # Merge on condition keys
    merge_keys = ['n', 'alpha', 'delay_yz', 'is_global']
    detection_summary_final = exp_detections[merge_keys + ['Recall', 'Recall_std', 'n_exp', 'lla_score_mean', 'lla_score_std']].merge(
        ctrl_detections[merge_keys + ['FP_rate', 'FP_rate_std', 'n_ctrl']],
        on=merge_keys,
        how='left'
    )
    
    # Calculate FDR = FP / (TP + FP)
    detection_summary_final['FDR'] = detection_summary_final.apply(
        lambda row: (row['FP_rate'] * row['n_ctrl']) / 
                    (row['Recall'] * row['n_exp'] + row['FP_rate'] * row['n_ctrl'])
                    if (row['Recall'] * row['n_exp'] + row['FP_rate'] * row['n_ctrl']) > 0 
                    else 0.0,
        axis=1
    )
    
    # Add batch timing if available
    if 'batch_time_seconds' in df.columns:
        # Get unique batch times per condition
        batch_times = df.groupby(merge_keys)['batch_time_seconds'].first()
        detection_summary_final = detection_summary_final.merge(
            batch_times.reset_index(),
            on=merge_keys,
            how='left'
        )
    
    # Round metrics
    detection_summary_final['FDR'] = detection_summary_final['FDR'].round(4)
    detection_summary_final['Recall'] = detection_summary_final['Recall'].round(4)
    if 'batch_time_seconds' in detection_summary_final.columns:
        detection_summary_final['batch_time_seconds'] = detection_summary_final['batch_time_seconds'].round(2)
    
    # Save updated detection summary
    detection_output = os.path.join(output_dir, 'detection_summary.csv')
    detection_summary_final.to_csv(detection_output, index=False)
    print(f"✓ Updated detection summary saved to: {detection_output}")
    print(f"  Columns: {list(detection_summary_final.columns)}")
    print(f"\nSample FDR values:")
    print(detection_summary_final[['alpha', 'Recall', 'FP_rate', 'FDR']].head(5))
    
    return detection_summary_final


def main():
    parser = argparse.ArgumentParser(description='Recalculate FDR from raw benchmark results')
    parser.add_argument('--raw', type=str, required=True,
                       help='Path to benchmark_raw_results.csv')
    parser.add_argument('--output-dir', type=str, required=True,
                       help='Output directory for updated detection_summary.csv')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.raw):
        print(f"Error: Raw results file not found: {args.raw}")
        return 1
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    recalculate_detection_summary(args.raw, args.output_dir)
    print("\n✓ FDR recalculation complete!")
    
    return 0


if __name__ == '__main__':
    exit(main())

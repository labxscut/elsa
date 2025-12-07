# LLA Benchmark Workflow Guide

## Quick Start: Alpha Sweep Analysis

### Step 1: Run Benchmark with Alpha Sweep

```powershell
# Simple alpha sweep (α from 0.8 to 1.0, step 0.02)
python benchmark_lla_unified.py --alpha_sweep --n_values 40 --window_fractions 1.0 --output_dir alpha_sweep_n40

# Custom alpha range
python benchmark_lla_unified.py --alpha_sweep --alpha_start 0.7 --alpha_end 1.0 --alpha_step 0.01 --n_values 40 --window_fractions 1.0 --output_dir alpha_sweep_custom

# With more replicates for smoother curves
python benchmark_lla_unified.py --alpha_sweep --n_replicates 100 --n_values 40 --window_fractions 1.0 --output_dir alpha_sweep_n40_r100
```

### Step 2: Visualize Results

```powershell
# Generate all plots from the results
python visualize_lla_results.py alpha_sweep_n40

# Filter to specific parameters
python visualize_lla_results.py alpha_sweep_n40 --n 40 --global_only

# Save plots to different directory
python visualize_lla_results.py alpha_sweep_n40 --output_dir my_plots
```

## Generated Plots

The visualization script creates:

1. **`detection_vs_alpha.png`** - Two-panel figure:
   - Recall (TPR) vs α
   - FPR vs α

2. **`detection_combined.png`** - Single plot with both Recall and FPR

3. **`localization_vs_alpha.png`** - Three-panel figure:
   - IoU (window overlap) vs α
   - Start error vs α
   - End error vs α

4. **`error_correlation.png`** - Scatter plot:
   - Start error vs End error (colored by α)
   - Shows correlation coefficient
   - Helps identify if errors are systematic

5. **`comprehensive_results.png`** - All-in-one publication figure:
   - 5-panel layout with all key metrics
   - Perfect for presentations/papers

## Example Workflows

### Workflow 1: Compare Different Sequence Lengths

```powershell
# Run for n=20
python benchmark_lla_unified.py --alpha_sweep --n_values 20 --window_fractions 1.0 --output_dir results_n20

# Run for n=40
python benchmark_lla_unified.py --alpha_sweep --n_values 40 --window_fractions 1.0 --output_dir results_n40

# Run for n=60
python benchmark_lla_unified.py --alpha_sweep --n_values 60 --window_fractions 1.0 --output_dir results_n60

# Visualize each
python visualize_lla_results.py results_n20
python visualize_lla_results.py results_n40
python visualize_lla_results.py results_n60
```

### Workflow 2: Compare Global vs Local Windows

```powershell
# Global window
python benchmark_lla_unified.py --alpha_sweep --n_values 40 --window_fractions 1.0 --output_dir results_global

# Local window (50% of sequence)
python benchmark_lla_unified.py --alpha_sweep --n_values 40 --window_fractions 0.5 --output_dir results_local

# Visualize with filters
python visualize_lla_results.py results_global --global_only
python visualize_lla_results.py results_local --local_only
```

### Workflow 3: Comprehensive Multi-Parameter Study

```powershell
# Run full benchmark (default: multiple n, α, delays, windows)
python benchmark_lla_unified.py --output_dir full_benchmark

# Visualize specific slices
python visualize_lla_results.py full_benchmark --n 40 --delay 0 --global_only
python visualize_lla_results.py full_benchmark --n 60 --delay 2 --local_only
```

## Understanding the Plots

### Detection Performance
- **Recall (TPR)**: Should increase with α (stronger association = easier detection)
- **FPR**: Should stay low (~0.05) regardless of α (controlled false positives)

### Localization Performance
- **IoU**: Should increase with α (better window detection)
- **Start/End Errors**: Should decrease (toward 0) with stronger α
  - Positive error = detected later than true
  - Negative error = detected earlier than true

### Error Correlation
- **High correlation**: Errors are systematic (consistent bias)
- **Low correlation**: Errors are random
- **Points near y=x line**: Start and end errors are similar
- **Color gradient**: Shows how correlation changes with α

## Tips

1. **Start simple**: Run with default parameters first
2. **Use alpha sweep**: Better for understanding detection threshold behavior
3. **Check error correlation**: Reveals if algorithm has systematic biases
4. **Multiple replicates**: Use 50-100 for smooth curves
5. **Filter when visualizing**: Focus on specific conditions of interest

## Files Organization

```
your_output_dir/
├── benchmark_raw_results.csv       # All individual experiment results
├── detection_summary.csv            # Recall, FPR by condition
├── localization_summary.csv         # IoU, errors by condition
├── detection_vs_alpha.png          # Detection plots
├── localization_vs_alpha.png       # Localization plots
├── error_correlation.png           # Error scatter plot
└── comprehensive_results.png       # All-in-one figure
```

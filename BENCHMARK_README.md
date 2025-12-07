# LLA Reliability Benchmark

This document explains how to use the `benchmark_lla_reliability.py` script to quantify the accuracy and consistency of Local Liquid Association (LLA) analysis for the **no-delay case** (delay_xy=0, delay_yz=0).

## Quick Start

### Minimal Test (3 runs, fast)
```powershell
python test_benchmark_mini.py
```

This runs 3 iterations at n=20 with reduced precision (100 permutations) to verify everything works.

### Full Benchmark (default: 50 runs × 4 configurations)
```powershell
python benchmark_lla_reliability.py
```

This runs:
- n ∈ {20, 40, 80, 100}
- corr_len = 8 (fixed)
- R = 50 runs per configuration
- Total: 200 simulations

Outputs:
- `lla_reliability_summary.csv` - Summary statistics per configuration
- Per-run details (optional, use `--per-run-csv`)

### Custom Configuration
```powershell
python benchmark_lla_reliability.py `
    --output-csv my_summary.csv `
    --per-run-csv my_runs.csv `
    --alpha 0.01 `
    --precision 500
```

## Metrics Explained

### Index Conventions

**Important**: The generator (`localsim_with_delay.py`) uses **0-based** array indices for `corr_start` and `corr_end`, but LLA output reports **1-based** indices (matching T1, T2, ... header convention).

The benchmark script automatically converts:
- Generator: `corr_start0=10, corr_end0=17` (0-based, inclusive)
- Truth for comparison: `s_true=11, e_true=18` (1-based, inclusive)
- LLA output: `Start_Z=11, End_Z=18` (1-based, inclusive) ✓ matches

### Predicted Interval (Consensus)

For the **no-delay case**, we define the predicted interval as the **intersection** of the three reported windows:

```
s_pred = max(Start_X, Start_Y, Start_Z)
e_pred = min(End_X, End_Y, End_Z)
```

If `e_pred < s_pred`, the predicted interval is **empty** (no consensus).

### Overlap Metrics

1. **Intersection length**:
   ```
   L_int = max(0, min(e_true, e_pred) - max(s_true, s_pred) + 1)
   ```

2. **Union length**:
   ```
   L_union = (e_true - s_true + 1) + (e_pred - s_pred + 1) - L_int
   ```

3. **Jaccard (overlap ratio)**:
   ```
   J = L_int / L_union
   ```
   - Range: [0, 1]
   - J=1: Perfect overlap
   - J=0: No overlap

4. **Coverage (recall)**:
   ```
   Coverage = L_int / (e_true - s_true + 1)
   ```
   - Fraction of true region captured by prediction
   - Range: [0, 1]

### Boundary Errors

1. **Start error**:
   ```
   e_start = s_pred - s_true
   ```
   - Positive: predicted start is late (misses early points)
   - Negative: predicted start is early (false positives)

2. **End error**:
   ```
   e_end = e_pred - e_true
   ```
   - Positive: predicted end is late (false positives)
   - Negative: predicted end is early (misses late points)

3. **MAE (Mean Absolute Error)**:
   - `MAE_start = mean(|e_start|)` across all runs
   - `MAE_end = mean(|e_end|)` across all runs

### Power (Detection Rate)

```
Power@α = fraction of runs where P ≤ α
```

- Measures how often LLA detects the association at significance level α (default: 0.05)
- Power = 1.0: Always detects (100% sensitivity)
- Power = 0.0: Never detects

**Standard Error**:
```
SE_Power ≈ sqrt(Power × (1 - Power) / R)
```

For R=400 and Power≈0.5 (worst case), SE ≈ 0.025 (±0.05 with ~95% CI).

## Output Files

### Summary CSV (`lla_reliability_summary.csv`)

Columns:
- `n`: Number of time points
- `corr_len`: Length of true correlation window
- `R`: Number of runs
- `J_mean`, `J_std`: Jaccard mean ± std
- `Coverage_mean`, `Coverage_std`: Coverage mean ± std
- `MAE_start_mean`, `MAE_start_std`: Start boundary error mean ± std
- `MAE_end_mean`, `MAE_end_std`: End boundary error mean ± std
- `Power`: Detection rate at α
- `Power_SE`: Standard error of power
- `LA_mean`, `LA_std`: LA score mean ± std

### Per-Run CSV (optional)

Columns:
- `n`, `seed`: Configuration and random seed
- `corr_start0`, `corr_end0`: True region (0-based)
- `s_true`, `e_true`: True region (1-based)
- `s_pred`, `e_pred`: Predicted region (1-based)
- `J`, `Coverage`: Overlap metrics
- `L_int`, `L_union`: Lengths
- `MAE_start`, `MAE_end`: Boundary errors
- `P`: P-value from LLA
- `LA`: LA score
- `power_hit`: Boolean (True if P ≤ α)

## Command-Line Options

```
python benchmark_lla_reliability.py --help
```

Options:
- `--output-csv`: Summary output path (default: `lla_reliability_summary.csv`)
- `--per-run-csv`: Per-run details output path (optional)
- `--alpha`: Significance level for power (default: 0.05)
- `--delay-limit`: Maximum delay for LLA `-d` parameter (default: 3)
- `--precision`: Permutations for p-value `-x` parameter (default: 1000)
- `--python`: Python executable (default: current interpreter)

## Interpretation Guidelines

### Good Performance

- **Jaccard ≥ 0.7**: Strong overlap
- **Coverage ≥ 0.8**: Captures most of true region
- **MAE_start, MAE_end ≤ 2**: Boundaries within ±2 time points
- **Power ≥ 0.8**: Detects 80%+ of cases

### Factors Affecting Performance

1. **Series length (n)**: Longer series → more stable estimation
2. **Correlation window length**: Longer windows → easier to detect
3. **Effect size (LA magnitude)**: Stronger association → higher power
4. **Noise level**: Lower SNR → reduced sensitivity

### Expected Trends

- As **n increases**: J, Coverage, and Power should increase; MAE should decrease
- As **corr_len increases**: Easier detection, higher power
- As **precision increases**: More accurate p-values, but slower runtime

## Customizing Configurations

To test different configurations, edit the `configs` list in `benchmark_lla_reliability.py`:

```python
configs = [
    {"n": 20, "corr_len": 6, "R": 100},
    {"n": 40, "corr_len": 8, "R": 100},
    {"n": 80, "corr_len": 10, "R": 100},
]
```

Or create a separate driver script:

```python
from benchmark_lla_reliability import run_benchmark

configs = [
    {"n": 30, "corr_len": 10, "R": 200},
]

summaries = run_benchmark(configs, per_run_csv="custom_runs.csv")
```

## Runtime Estimates

Per-run time ≈ 5-10 seconds (depends on n and precision):
- Mini-test (3 runs): ~20 seconds
- Default (200 runs): ~20-30 minutes
- Full (400 runs × 4 configs): ~2-3 hours

Tips for faster testing:
- Reduce `--precision` (e.g., 100 instead of 1000) for development
- Use smaller R (e.g., 10-20) for quick checks
- Run on subset of n values

## Extending to Delay Cases

For future work with delays (delay_xy ≠ 0 or delay_yz ≠ 0):

1. **Alignment to Z axis**:
   ```python
   # Map X and Y windows back to Z timeline
   d_yx_hat = result["Delay_Y-X"]
   d_zy_hat = result["Delay_Z-Y"]
   
   X_on_Z_start = result["Start_X"] - (d_zy_hat + d_yx_hat)
   X_on_Z_end = result["End_X"] - (d_zy_hat + d_yx_hat)
   
   Y_on_Z_start = result["Start_Y"] - d_zy_hat
   Y_on_Z_end = result["End_Y"] - d_zy_hat
   
   # Consensus on Z axis
   s_pred = max(X_on_Z_start, Y_on_Z_start, result["Start_Z"])
   e_pred = min(X_on_Z_end, Y_on_Z_end, result["End_Z"])
   ```

2. **Delay accuracy metrics**:
   ```python
   MAE_delay_yx = mean(|d_yx_hat - d_yx_true|)
   MAE_delay_zy = mean(|d_zy_hat - d_zy_true|)
   ```

3. **Exact-match rate**: Fraction of runs where predicted delays exactly match true delays

## Troubleshooting

### Import Errors
Ensure you run from the repository root:
```powershell
cd d:\wsl\elsa
python benchmark_lla_reliability.py
```

### Subprocess Failures
Check that:
- `localsim_with_delay.py` is in the current directory
- `lla/lla_compute.py` exists and is executable
- Python can import `numpy` and other dependencies

### All Runs Fail
- Verify `lla_compute.py` works standalone
- Check error output in stderr
- Try with `--precision 10` for ultra-fast debugging

### Unexpected Metrics
- Check index conventions in your generator
- Verify LLA output format matches expected columns
- Use `--per-run-csv` to inspect individual runs

## Citation and References

When reporting these metrics, describe them as:

> "We quantified LLA reliability using Jaccard overlap ratio (intersection-over-union of predicted and true intervals), coverage (fraction of true region captured), mean absolute boundary errors, and statistical power (detection rate at α=0.05) across R repeated simulations with varying time series lengths."

This provides a comprehensive, interpretable, and statistically principled assessment of your LLA method's performance.

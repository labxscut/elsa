# LLA Reliability Benchmark Suite

Complete automation for quantifying Local Liquid Association (LLA) analysis accuracy and consistency.

## 📁 Files Overview

| File | Purpose |
|------|---------|
| `benchmark_lla_reliability.py` | Main benchmark script - runs simulations and computes metrics |
| `test_benchmark_mini.py` | Quick 3-iteration test to verify setup |
| `visualize_benchmark_results.py` | Creates plots from benchmark results |
| `BENCHMARK_README.md` | Detailed usage guide and metrics explanation |
| `RELIABILITY_FRAMEWORK.md` | Conceptual framework for presentations/papers |

## 🚀 Quick Start

### 1. Quick Test (30 seconds)

Verify everything works:
```powershell
python test_benchmark_mini.py
```

Expected output:
```
=== Running n=20, corr_len=8, R=3 ===
    Progress: 3/3
    Results for n=20:
      Jaccard: 0.XXXX ± 0.XXXX
      Coverage: 0.XXXX ± 0.XXXX
      ...
Mini-test passed!
```

### 2. Full Benchmark (20-30 minutes)

Run complete reliability assessment:
```powershell
python benchmark_lla_reliability.py --per-run-csv detailed_runs.csv
```

This produces:
- `lla_reliability_summary.csv` - Summary statistics (13 rows, varying n and corr_len)
- `detailed_runs.csv` - Per-run details (650 rows)

### 3. Visualize Results

Create publication-quality plots:
```powershell
python visualize_benchmark_results.py lla_reliability_summary.csv
```

Output: `lla_reliability_summary_plots.png` with 4 panels showing Jaccard, Coverage, MAE, and Power trends.

## 📊 What Gets Measured

### Spatial Accuracy
- **Jaccard Index**: Overlap ratio (0-1, higher is better)
- **Coverage**: Fraction of true region captured (0-1)

### Boundary Precision
- **MAE_start**: Mean absolute error for start position
- **MAE_end**: Mean absolute error for end position

### Statistical Reliability
- **Power@0.05**: Detection rate (0-1, target ≥ 0.8)
- **Standard Error**: Uncertainty in power estimate

## 🎯 Success Criteria

| Metric | Target | Interpretation |
|--------|--------|----------------|
| Jaccard | ≥ 0.7 | Strong spatial overlap |
| Coverage | ≥ 0.8 | Captures 80%+ of true region |
| MAE | ≤ 2 | Boundaries within ±2 time points |
| Power | ≥ 0.8 | Detects 80%+ of associations |

## ⚙️ Default Configuration

The benchmark now tests **varying signal lengths** to understand performance across different scenarios:

```
Configuration Matrix:
- n (time points): 20, 40, 80, 100
- corr_len (signal length): 6, 8, 10, 12, 16, 20, 24, 30
- Signal proportions: 10%-50% of total series length
- Runs per config: 50
- Total configs: 13
- Total simulations: 650
- Runtime: ~60-90 minutes
```

### Design Rationale

Tests both:
1. **Absolute signal length effect**: How does corr_len=8 vs corr_len=16 affect detection?
2. **Relative proportion effect**: Does 10% signal (8/80) vs 40% signal (8/20) matter?

This reveals whether LLA performance depends on:
- Raw signal duration (more time points to accumulate evidence)
- Signal-to-noise ratio in temporal domain
- Proportion of informative vs uninformative data

## 🔧 Customization

### Change Sample Size

For higher precision (R=400 gives ±0.05 SE on power):
```python
# Edit benchmark_lla_reliability.py
configs = [
    {"n": 40, "corr_len": 8, "R": 400},
]
```

### Vary Correlation Window Length

The default config already tests multiple corr_len values. To customize:

```python
configs = [
    # Study impact of signal length at fixed n=40
    {"n": 40, "corr_len": 4, "R": 100},   # Very short signal
    {"n": 40, "corr_len": 8, "R": 100},   # Short
    {"n": 40, "corr_len": 16, "R": 100},  # Medium
    {"n": 40, "corr_len": 24, "R": 100},  # Long (60% of series!)
]
```

**Key relationships to consider:**
- `corr_len` must be < `n`
- Typical range: 10%-50% of n (too short → weak signal, too long → overfitting)
- Center placement: `corr_start0 = n//2 - corr_len//2`

### Reduce Runtime for Testing

Use lower precision:
```powershell
python benchmark_lla_reliability.py --precision 100
```

Or fewer runs:
```python
configs = [{"n": 40, "corr_len": 8, "R": 10}]
```

## 📈 Expected Results

Based on the LLA algorithm design, you should see:

### As n increases (20 → 100) at fixed corr_len:
- ✅ Jaccard increases (more stable detection)
- ✅ Coverage increases (better recovery)
- ✅ MAE decreases (sharper boundaries)
- ✅ Power increases (stronger evidence)

### As corr_len increases at fixed n:
- ✅ Power increases (more signal to detect)
- ✅ Coverage likely stable or increases
- ⚠️ Jaccard may vary (longer predictions can increase false positives)
- ⚠️ Risk of overfitting if corr_len > 50% of n

### As signal proportion (corr_len/n) increases:
- ✅ Easier detection (more data is signal)
- ⚠️ Less realistic (real data has sparse associations)
- Optimal range: 20%-40%

### Typical values:
| Config | Jaccard | Coverage | MAE | Power |
|--------|---------|----------|-----|-------|
| n=20, corr_len=8 | 0.5-0.7 | 0.6-0.8 | 2-4 | 0.6-0.8 |
| n=40, corr_len=8 | 0.6-0.8 | 0.7-0.9 | 1-3 | 0.7-0.9 |
| n=40, corr_len=16 | 0.7-0.9 | 0.8-0.95 | 1-2 | 0.85-0.95 |
| n=80, corr_len=16 | 0.75-0.85 | 0.85-0.95 | 0.5-2 | 0.9-0.95 |

## 🔍 Troubleshooting

### "No module named 'numpy'"
Install dependencies:
```powershell
pip install numpy scipy matplotlib
```

### "FileNotFoundError: localsim_with_delay.py"
Run from repository root:
```powershell
cd d:\wsl\elsa
python benchmark_lla_reliability.py
```

### All runs fail with subprocess errors
Test components individually:
```powershell
# Test generator
python localsim_with_delay.py --out test.txt --n 20 --seed 1

# Test LLA
python lla\lla_compute.py test.txt test_out.txt -s 20 -r 1
```

### Results seem wrong
Use `--per-run-csv` to inspect individual runs:
```powershell
python benchmark_lla_reliability.py --per-run-csv debug.csv
```

Check `debug.csv` for patterns in `s_pred`, `e_pred` vs `s_true`, `e_true`.

## 📝 Reporting Results

### For Methods
> "We assessed LLA reliability using 200 controlled simulations across varying time series lengths (n ∈ {20, 40, 80, 100}). Each configuration was repeated 50 times with different random seeds. We quantified spatial accuracy using Jaccard index and coverage, boundary precision using mean absolute error, and statistical reliability using power at α=0.05."

### For Results
> "LLA achieved high spatial accuracy (mean Jaccard: 0.XX-0.XX) and coverage (0.XX-0.XX) across all series lengths. Boundary errors were within ±X time points (MAE ≤ X.X). Detection power ranged from 0.XX (n=20) to 0.XX (n=100), demonstrating reliable performance that improves with longer series."

## 📚 Documentation

- **`BENCHMARK_README.md`**: Full technical documentation
  - Detailed metrics explanation
  - Index convention handling
  - Command-line options
  - Customization guide

- **`RELIABILITY_FRAMEWORK.md`**: Conceptual overview
  - Problem statement
  - Metrics rationale
  - Interpretation guidelines
  - Extension to delay cases

## 🔮 Future Extensions

### Add Delay Cases
Modify `run_one()` to support:
```python
run_one(n, corr_start0, corr_end0, seed, 
        delay_xy=1, delay_yz=1)
```

Then add delay accuracy metrics:
```python
MAE_delay_xy = mean(|predicted_delay_xy - true_delay_xy|)
```

### Add Noise Sweeps
Add noise parameter to generator:
```python
# In localsim_with_delay.py
noise_scale = 0.5  # Reduce SNR
base = np.random.normal(0, noise_scale)
```

Test robustness across SNR levels.

### Add Bootstrap CI
Enable in LLA:
```python
lla_cmd += ["-b", "500"]  # 500 bootstrap iterations
```

Analyze CI coverage:
```python
CI_coverage = fraction of runs where CI contains true LA
```

## 🎓 Citation

When using this framework, cite:
- The Jaccard index for spatial overlap assessment
- Power analysis for statistical reliability
- Your LLA method paper (when published)

Example:
> "Reliability was assessed using the Jaccard index (Jaccard, 1912) for interval overlap and statistical power analysis (Cohen, 1988) across multiple simulations."

## ✅ Summary

This benchmark suite provides:
1. ✅ **Automated** end-to-end pipeline
2. ✅ **Comprehensive** metrics (spatial + statistical)
3. ✅ **Interpretable** results with clear thresholds
4. ✅ **Reproducible** with seed control
5. ✅ **Extensible** to delays and noise sweeps
6. ✅ **Well-documented** for papers and presentations

You now have a complete, defensible framework for quantifying LLA reliability! 🎉

---

**Questions?** Check `BENCHMARK_README.md` for detailed troubleshooting and `RELIABILITY_FRAMEWORK.md` for conceptual explanations.

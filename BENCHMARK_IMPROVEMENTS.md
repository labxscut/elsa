# Benchmark Script Improvements

## Changes Made

### 1. ✅ Variable Signal Length (`corr_len`)

**Problem**: Original config had fixed `corr_len=8`, preventing investigation of signal length effects.

**Solution**: Extended configuration matrix to test multiple `corr_len` values:

```python
configs = [
    # Short series (n=20)
    {"n": 20, "corr_len": 6, "R": 50},   # 30% of series
    {"n": 20, "corr_len": 8, "R": 50},   # 40% of series
    {"n": 20, "corr_len": 10, "R": 50},  # 50% of series
    
    # Medium series (n=40)
    {"n": 40, "corr_len": 6, "R": 50},   # 15% of series
    {"n": 40, "corr_len": 8, "R": 50},   # 20% of series
    {"n": 40, "corr_len": 12, "R": 50},  # 30% of series
    {"n": 40, "corr_len": 16, "R": 50},  # 40% of series
    
    # Long series (n=80)
    {"n": 80, "corr_len": 8, "R": 50},   # 10% of series
    {"n": 80, "corr_len": 16, "R": 50},  # 20% of series
    {"n": 80, "corr_len": 24, "R": 50},  # 30% of series
    
    # Very long series (n=100)
    {"n": 100, "corr_len": 10, "R": 50}, # 10% of series
    {"n": 100, "corr_len": 20, "R": 50}, # 20% of series
    {"n": 100, "corr_len": 30, "R": 50}, # 30% of series
]
```

**Benefits**:
- Tests **13 configurations** instead of 4
- Investigates both **absolute** signal length (6 vs 30 time points)
- Investigates **relative** signal proportion (10% vs 50% of series)
- Reveals whether performance depends on raw duration or proportion

**Validation Added**:
```python
# Validate that corr_len fits within n
if corr_len >= n:
    print(f"WARNING: corr_len={corr_len} >= n={n}, skipping", file=sys.stderr)
    continue

# Ensure window fits within valid range [0, n)
if corr_end0 >= n:
    corr_end0 = n - 1
    corr_start0 = corr_end0 - corr_len + 1
```

### 2. ✅ Fixed File Path Resolution

**Problem**: Original code used relative paths that assumed:
- `benchmark_lla_reliability.py` in root
- `localsim_with_delay.py` in root
- `lla/lla_compute.py` in subdirectory

But paths were hardcoded as strings, causing failures if script run from different directories.

**Solution**: Use `Path(__file__).parent` to resolve paths dynamically:

```python
# Before (fragile):
gen_cmd = [python_exe, "localsim_with_delay.py", ...]
lla_cmd = [python_exe, os.path.join("lla", "lla_compute.py"), ...]

# After (robust):
script_dir = Path(__file__).parent.absolute()
gen_script = script_dir / "localsim_with_delay.py"
lla_script = script_dir / "lla" / "lla_compute.py"

gen_cmd = [python_exe, str(gen_script), ...]
lla_cmd = [python_exe, str(lla_script), ...]
```

**Benefits**:
- Works regardless of current working directory
- Portable across systems (handles Windows/Unix path separators)
- Explicit about file locations
- Fails fast with clear error if scripts missing

### 3. ✅ Enhanced Progress Reporting

**Problem**: User couldn't see configuration matrix details at startup.

**Solution**: Added comprehensive startup summary:

```
======================================================================
LLA Reliability Benchmark (No-Delay Case)
======================================================================
Alpha: 0.05
Delay limit: 3
Precision: 1000
Python: C:\...\python.exe
Configurations: 13

Configuration matrix:
  n values: [20, 40, 80, 100]
  corr_len values: [6, 8, 10, 12, 16, 20, 24, 30]
  Signal proportions tested: ['6/20', '8/20', '8/40', '10/20', ...]
```

## Research Questions Enabled

With these changes, you can now investigate:

### Q1: Does Absolute Signal Length Matter?

Compare same `corr_len` across different `n`:
- `n=40, corr_len=8` vs `n=80, corr_len=8`
- If performance differs → series length provides context beyond signal
- If similar → detection is primarily about signal duration

### Q2: Does Relative Signal Proportion Matter?

Compare different proportions:
- `n=20, corr_len=8` (40%) vs `n=40, corr_len=8` (20%)
- Higher proportion → easier detection?
- Or does absolute length dominate?

### Q3: What's the Optimal Signal Length Range?

By testing corr_len from 6 to 30:
- Find minimum viable signal length for reliable detection
- Identify point of diminishing returns (more signal ≠ better if overfitting)
- Establish guidelines: "For n=X, use corr_len ≥ Y for Power ≥ 0.8"

### Q4: Is There a Universal "Sweet Spot"?

Look for consistent pattern:
- Does 20-30% proportion work best across all n?
- Or does it depend on absolute values?

## Expected Findings

### Hypothesis 1: Longer Absolute Signal → Better Performance
At fixed proportion (e.g., 20%):
- `n=40, corr_len=8` < `n=80, corr_len=16` < `n=100, corr_len=20`
- More time points to accumulate evidence

### Hypothesis 2: Higher Proportion → Better Detection (up to limit)
At fixed n=40:
- `corr_len=6` (15%) < `corr_len=8` (20%) < `corr_len=16` (40%)
- But may plateau or degrade beyond 50% (overfitting)

### Hypothesis 3: Interaction Effects
Best performance when **both** n and corr_len are large:
- `n=100, corr_len=30` should dominate all other configs
- Worst: `n=20, corr_len=6` (short series, short signal)

### Hypothesis 4: Boundary Precision Improves with Context
MAE should decrease as `n - corr_len` increases:
- More uninformative background → clearer transition points
- `n=80, corr_len=8` should have sharper boundaries than `n=20, corr_len=8`

## Analysis Recommendations

### 1. Create Heat Maps

Plot metrics as function of (n, corr_len):
```python
import matplotlib.pyplot as plt
import numpy as np

# After loading summary CSV
n_vals = sorted(set(r['n'] for r in data))
corr_lens = sorted(set(r['corr_len'] for r in data))

# Create matrix
J_matrix = np.zeros((len(n_vals), len(corr_lens)))
for i, n in enumerate(n_vals):
    for j, cl in enumerate(corr_lens):
        row = [r for r in data if r['n']==n and r['corr_len']==cl]
        J_matrix[i,j] = row[0]['J_mean'] if row else np.nan

plt.imshow(J_matrix, aspect='auto', cmap='viridis')
plt.xlabel('corr_len')
plt.ylabel('n')
plt.colorbar(label='Jaccard')
plt.title('Jaccard vs (n, corr_len)')
```

### 2. Proportion Analysis

Add derived column:
```python
for row in data:
    row['proportion'] = row['corr_len'] / row['n']
```

Plot metrics vs proportion to see if there's a universal optimal range.

### 3. Regression Models

Fit model to predict performance:
```python
from sklearn.linear_model import LinearRegression

X = [[r['n'], r['corr_len'], r['n']*r['corr_len']] for r in data]
y = [r['Power'] for r in data]

model = LinearRegression().fit(X, y)
# Coefficients reveal relative importance
```

### 4. Practical Guidelines

Based on results, produce table:
```
| n   | Min corr_len for Power≥0.8 | Recommended corr_len |
|-----|----------------------------|---------------------|
| 20  | 8                          | 8-10                |
| 40  | 6                          | 8-12                |
| 80  | 6                          | 12-20               |
| 100 | 8                          | 15-25               |
```

## Documentation Updates

All docs updated to reflect new design:
- ✅ `BENCHMARK_QUICKSTART.md` - Updated config matrix and expectations
- ✅ `benchmark_lla_reliability.py` - Inline comments explain validation
- ✅ `test_benchmark_mini.py` - Shows proper config structure
- ✅ This file (`BENCHMARK_IMPROVEMENTS.md`) - Explains rationale

## Runtime Impact

- **Before**: 4 configs × 50 runs = 200 simulations (~20-30 min)
- **After**: 13 configs × 50 runs = 650 simulations (~60-90 min)

**Trade-off**: 3× longer runtime for much richer insights.

**Mitigation**: 
- Use `--precision 100` for quick exploratory runs
- Reduce R to 20 for testing
- Run overnight for final publication results with R=400

## Backward Compatibility

Old configs still work:
```python
configs = [{"n": 40, "corr_len": 8, "R": 50}]
```

Just now supports more flexible exploration.

## Future Extensions

With this infrastructure, easy to add:
1. **Noise sweeps**: Add `noise_scale` to generator
2. **Delay cases**: Modify `run_one()` to accept delay_xy, delay_yz
3. **Effect size**: Vary strength of X-Y-Z association
4. **Different placements**: Test edge vs center windows
5. **Asymmetric windows**: corr_len varies for X, Y, Z

All enabled by the flexible config dict structure!

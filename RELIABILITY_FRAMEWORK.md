# LLA Reliability Measurement Framework

## Overview

This document provides a concise explanation of how we measure the reliability of Local Liquid Association (LLA) analysis, suitable for presentations or papers.

---

## Problem Statement

**Goal**: Quantify how accurately and consistently LLA detects the true local association interval in time series data.

**Challenge**: Unlike global methods that report a single score, LLA reports:
- A time interval where the association is strongest
- Directional delays between variables
- Statistical significance (p-value)

We need metrics that assess **spatial accuracy** (interval overlap), **boundary precision**, and **statistical power** (detection reliability).

---

## Ground Truth Definition

### Simulation Design

We generate controlled synthetic data where the "truth" is known:

1. **Three time series**: X, Y, Z (length n)
2. **True correlation region**: Time points [s_true, e_true] (1-based, inclusive)
3. **Within correlation region**: X, Y, Z exhibit liquid association
   - X and Y are modulated by Z
   - Constraint: X·Y·Z has consistent sign (positive or negative)
4. **Outside correlation region**: Independent Gaussian noise
5. **No delay initially**: delay_xy = 0, delay_yz = 0

**Index Convention**: 
- Generator uses 0-based Python indices internally
- Output and comparisons use 1-based indices (matching T1, T2, ... convention)

---

## Predicted Interval Extraction

From LLA output, we extract:
- `Start_X`, `Start_Y`, `Start_Z`: Start positions for each variable
- `End_X`, `End_Y`, `End_Z`: End positions for each variable

**Consensus predicted interval** (no-delay case):
```
s_pred = max(Start_X, Start_Y, Start_Z)
e_pred = min(End_X, End_Y, End_Z)
```

Rationale: The intersection represents the region where **all three variables** show strong local association.

If `e_pred < s_pred`, the prediction is considered **empty** (no consensus).

---

## Primary Metrics

### 1. Jaccard Index (Overlap Ratio)

**Definition**:
$$
J = \frac{L_{\text{int}}}{L_{\text{union}}}
$$

Where:
- $L_{\text{int}} = \max(0, \min(e_{\text{true}}, e_{\text{pred}}) - \max(s_{\text{true}}, s_{\text{pred}}) + 1)$
- $L_{\text{union}} = L_{\text{true}} + L_{\text{pred}} - L_{\text{int}}$

**Range**: [0, 1]
- J = 1: Perfect overlap (predicted = true)
- J = 0: No overlap (complete miss)

**Interpretation**: Balanced measure of spatial accuracy (penalizes both false positives and false negatives).

### 2. Coverage (Recall)

**Definition**:
$$
\text{Coverage} = \frac{L_{\text{int}}}{e_{\text{true}} - s_{\text{true}} + 1}
$$

**Range**: [0, 1]
- Coverage = 1: Captures entire true region
- Coverage = 0: Misses entire true region

**Interpretation**: Sensitivity to the true association region (fraction recovered).

### 3. Boundary Errors

**Start error**:
$$
e_{\text{start}} = s_{\text{pred}} - s_{\text{true}}
$$

**End error**:
$$
e_{\text{end}} = e_{\text{pred}} - e_{\text{true}}
$$

**Aggregated metric** (across R runs):
$$
\text{MAE}_{\text{start}} = \frac{1}{R}\sum_{i=1}^{R}|e_{\text{start}}^{(i)}|
$$
$$
\text{MAE}_{\text{end}} = \frac{1}{R}\sum_{i=1}^{R}|e_{\text{end}}^{(i)}|
$$

**Interpretation**: Average boundary precision (in time points). MAE ≤ 2 indicates boundaries are accurate within ±2 time points.

### 4. Statistical Power

**Definition**:
$$
\text{Power}(\alpha) = P(p\text{-value} \leq \alpha \mid \text{true association exists})
$$

**Empirical estimate**:
$$
\hat{\text{Power}} = \frac{\#\{\text{runs with } p \leq \alpha\}}{R}
$$

**Standard error** (binomial):
$$
SE = \sqrt{\frac{\hat{\text{Power}}(1 - \hat{\text{Power}})}{R}}
$$

**Interpretation**: Detection reliability. Power = 0.8 means 80% of true associations are detected at significance level α.

**Sample size for precision**: For ±0.05 margin with 95% confidence and worst-case Power ≈ 0.5, need R ≥ 384 (we use R = 400+).

---

## Experimental Design

### Configuration Grid

| Parameter | Values | Notes |
|-----------|--------|-------|
| n (time points) | 20, 40, 80, 100 | Vary series length |
| corr_len | 8 (fixed) | True correlation window length |
| corr_position | Centered: $\lfloor n/2 \rfloor - \text{corr\_len}/2$ | Consistent placement |
| R (repetitions) | 50-400 | Higher R → better SE |
| α (significance) | 0.05 | Standard threshold |
| delay | 0 (initially) | Extend later |

### Procedure

For each configuration (n, corr_len):
1. **Loop** over R random seeds
2. **Generate** synthetic data with `localsim_with_delay.py`
3. **Analyze** with `lla_compute.py` (1000 permutations, delay_limit=3)
4. **Parse** output: extract predicted intervals, p-value, LA score
5. **Compute** metrics: J, Coverage, MAE_start, MAE_end, Power
6. **Aggregate** across runs: mean ± std for each metric

---

## Interpretation Guidelines

### Performance Thresholds

| Metric | Good | Acceptable | Poor |
|--------|------|------------|------|
| Jaccard | ≥ 0.7 | 0.5-0.7 | < 0.5 |
| Coverage | ≥ 0.8 | 0.6-0.8 | < 0.6 |
| MAE (boundary) | ≤ 2 | 2-5 | > 5 |
| Power@0.05 | ≥ 0.8 | 0.6-0.8 | < 0.6 |

### Expected Trends

1. **As n increases** (longer series):
   - Jaccard ↑ (more stable interval detection)
   - Coverage ↑ (better recovery of true region)
   - MAE ↓ (more precise boundaries)
   - Power ↑ (stronger statistical evidence)

2. **As corr_len increases** (longer true window):
   - Power ↑ (easier to detect)
   - Jaccard may vary (depends on prediction strategy)

3. **As noise increases**:
   - All metrics degrade (lower signal-to-noise ratio)

---

## Reporting Template

### For Methods Section

> "We assessed LLA reliability using controlled simulations with known ground truth. For each configuration (time series length n and correlation window length corr_len), we generated R=400 synthetic datasets with randomized seeds. We quantified accuracy using: (1) Jaccard index (intersection-over-union of predicted and true intervals), (2) coverage (fraction of true region captured), (3) mean absolute error for start/end boundaries, and (4) statistical power (detection rate at α=0.05). Standard errors for power were computed using binomial variance."

### For Results Section

> "Across series lengths n ∈ {20, 40, 80, 100}, LLA achieved mean Jaccard scores of 0.XX-0.XX (indicating strong spatial overlap), coverage of 0.XX-0.XX (capturing XX-XX% of true associations), and boundary errors within ±X.X time points. Statistical power ranged from 0.XX to 0.XX, demonstrating reliable detection at α=0.05. Performance improved with longer series (n=100: J=0.XX, Power=0.XX), consistent with increased statistical evidence in longer observations."

---

## Extension to Delay Cases

For delay ≠ 0:

1. **Align predicted windows to Z timeline**:
   - Subtract predicted delays from X and Y windows
   - Compute consensus on Z axis

2. **Additional metrics**:
   - Delay accuracy: MAE for delay_xy and delay_zy
   - Exact-match rate: fraction with correct delays

3. **Challenge**: Delay estimation adds complexity (error propagates to interval alignment)

---

## Advantages of This Framework

1. **Comprehensive**: Covers spatial accuracy, boundary precision, and statistical reliability
2. **Interpretable**: Each metric has clear operational meaning
3. **Statistically principled**: SE for power, proper handling of uncertainty
4. **Extendable**: Can add noise sweeps, delay cases, and effect size variations
5. **Reproducible**: Automated pipeline from generation → analysis → metrics → aggregation

---

## References for Metrics

- **Jaccard Index**: Standard in information retrieval and image segmentation
- **Coverage (Recall)**: Common in classification evaluation
- **MAE**: Standard regression metric for boundary precision
- **Power**: Fundamental statistical concept (Neyman-Pearson framework)

This multi-faceted approach provides a thorough, defensible assessment of LLA performance that addresses both spatial accuracy and statistical reliability.

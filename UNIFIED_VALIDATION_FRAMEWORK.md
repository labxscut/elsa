# Unified LLA Validation Framework

## Overview

This framework provides a comprehensive system for validating LLA's capabilities across three key dimensions:

1. **Detection Power**: Distinguish true associations from noise
2. **Localization Accuracy**: Identify correct regulation windows
3. **Delay Estimation**: Recover true time delays

## Core Components

### 1. Data Generation (`gen_unified_triplets.py`)

Generates (X, Y, Z) triplets with controlled properties:

#### Parameters

- **n**: Sequence length (20, 40, 60, 80, 100)
- **α**: Association strength in [0, 1]
- **Window I**: Regulation region
  - Global: I = [0, n) (entire sequence)
  - Local: I ⊂ [0, n) with |I| ≥ min(10, n/2)
- **delay_yz**: Y-Z time delay (X-Y always synchronous)
- **is_control**: Generate null data (pure noise)

#### Data Model

**Experimental Group** (with association):

```
For each position i:
  If i ∈ I (inside regulation window):
    y_idx = i + delay_yz
    If y_idx ∈ I and valid:
      X[i] = α * Y[y_idx] + (1-α) * ε[i]
    Else:
      X[i] = ε'[i]  (independent noise)
  Else (outside window):
    X[i] = ε'[i]  (independent noise)

Where:
  Y[i] ~ N(0,1)       (base variable)
  ε[i] ~ N(0,1)       (association noise)
  ε'[i] ~ N(0,1)      (independent noise)
  Z[i] = 1 if i ∈ I, else 0  (regulation indicator: 1 inside window, 0 outside)
  
Note: 
  - Global case: Z = [1, 1, ..., 1] (all ones, since I = entire sequence)
  - Local case: Z = [0, ..., 1, 1, ..., 0] (ones only in window I)
```

**Control Group** (null):

```
X[i] = ε'[i]  for all i  (independent of Y)
Z[i] = 1 if i ∈ I, else 0  (same structure as experimental)
```

#### Usage Examples

```bash
# Global association, α=0.8, no delay
python gen_unified_triplets.py --n 40 --alpha 0.8 --out global_a08.txt

# Local association, window [10, 30), α=1.0
python gen_unified_triplets.py --n 40 --alpha 1.0 \
  --window-start 10 --window-end 30 --out local_w1030.txt

# Local with Y-Z delay=3
python gen_unified_triplets.py --n 60 --alpha 0.6 \
  --window-start 15 --window-end 45 --delay-yz 3 --out local_delay3.txt

# Control group (pure noise)
python gen_unified_triplets.py --n 40 --control --out control.txt
```

### 2. Benchmark Framework (`benchmark_lla_unified.py`)

Systematically tests LLA across parameter space and evaluates:

#### Metrics

**Detection Metrics**:

- **Power**: P(detected | experimental group)
- **FPR**: P(detected | control group)
- **|LLA score|**: Effect size

**Localization Metrics**:

- **Window overlap**: IoU between true and detected windows
- **Start error**: Detected_start - True_start
- **End error**: Detected_end - True_end
- **Exact match rate**: P(perfect localization)

**Delay Metrics** (when delay_yz ≠ 0):

- **Delay accuracy**: P(detected_delay == true_delay)
- **Delay error**: Mean absolute error
- **Delay distribution**: Histogram of errors

#### Usage

```bash
# Quick test (small parameter space)
python benchmark_lla_unified.py \
  --n-values 20 40 \
  --alpha-values 0.8 1.0 \
  --delay-values 0 2 \
  --window-fractions 1.0 0.5 \
  --n-replicates 30 \
  --output-dir results_quick

# Full benchmark (comprehensive)
python benchmark_lla_unified.py \
  --n-values 20 40 60 80 100 \
  --alpha-values 0.2 0.4 0.6 0.8 1.0 \
  --delay-values 0 1 2 3 \
  --window-fractions 1.0 0.75 0.5 \
  --n-replicates 50 \
  --delay-limit 4 \
  --precision 1000 \
  --output-dir results_full
```

#### Output Files

1. **benchmark_raw_results.csv**: All individual experiments

   - Columns: n, alpha, delay_yz, is_global, is_control, detected, p_value, lla_score, window_overlap, delay_error, etc.
2. **detection_summary.csv**: Detection power by condition

   - Grouped by: n, α, delay, window type, control/experimental
   - Metrics: mean/std of detection rate, |LLA score|
3. **localization_summary.csv**: Localization accuracy

   - Grouped by: n, α, delay, window type
   - Metrics: mean/std of window overlap, start/end errors
4. **delay_summary.csv**: Delay estimation accuracy

   - Grouped by: n, α, delay
   - Metrics: mean/std of delay accuracy, delay errors

## Validation Strategy

### Phase 1: Detection Validation

**Goal**: Verify LLA can distinguish signal from noise

**Test Cases**:

- Global, α=1.0, delay=0 (strongest signal)
- Global, α=0.4, delay=0 (weak signal)
- Control groups (pure noise)

**Expected Results**:

- Power ≈ 100% for α=1.0
- Power > 80% for α≥0.6
- FPR ≈ 5% (α=0.05)

### Phase 2: Localization Validation

**Goal**: Verify LLA identifies correct regulation windows

**Test Cases**:

- Local, centered windows of varying sizes
- α ∈ {0.6, 0.8, 1.0}
- delay=0 initially

**Expected Results**:

- Window overlap > 80% for α≥0.8
- Start/end errors within ±2 positions
- Better localization with larger α

### Phase 3: Delay Validation

**Goal**: Verify LLA recovers true time delays

**Test Cases**:

- Global/Local with delay ∈ {1, 2, 3}
- α ∈ {0.6, 0.8, 1.0}
- Various n values

**Expected Results**:

- Delay accuracy > 90% for α≥0.8
- Mean delay error ≈ 0
- Accuracy decreases with smaller α

### Phase 4: Robustness Testing

**Goal**: Characterize performance across parameter space

**Test Cases**:

- All combinations of n, α, delay, window size
- Identify failure modes
- Determine minimum detectable α

**Expected Results**:

- Monotonic relationship: Power ↑ as α ↑, n ↑
- Localization improves with longer windows
- Delay estimation robust when Power is high

## Comparison with Old Scripts

### `localsim.py` (Old)

- **Scope**: Fixed local window [10, 20], n=20, binary Z ∈ {-1, +1}
- **Limitation**: Single scenario, no parameter sweep
- **Use Case**: Quick sanity check

### `localsim_with_delay.py` (Old)

- **Scope**: Adds delay_xy + delay_yz
- **Limitation**: Still supports X-Y delay (now removed in C++ code)
- **Use Case**: Testing delay mechanism

### `gen_global_triplets.py` & `benchmark_lla_detection.py`(Old)
- only global case with $X=\pm Y$.

### `gen_unified_triplets.py` & `benchmark_lla_unified.py`(New)

- **Scope**: Unified framework for global/local, variable α, delay_yz only
- **Advantages**:
  - Consistent with simplified delay model (X-Y synchronous)
  - Supports full parameter space
  - Clear experimental/control distinction
  - Validates window constraints
- **Use Case**: Production validation experiments

## Interpretation Guide

### Detection Results

| Power  | FPR  | Interpretation                      |
| ------ | ---- | ----------------------------------- |
| >95%   | ≈5% | Excellent - strong signal detection |
| 80-95% | ≈5% | Good - reliable detection           |
| 60-80% | ≈5% | Moderate - signal exists but weak   |
| <60%   | ≈5% | Poor - insufficient power           |
| Any    | >10% | Problem - inflated false positives  |

### Localization Results

| Window Overlap (IoU) | Interpretation                |
| -------------------- | ----------------------------- |
| >0.9                 | Excellent localization        |
| 0.7-0.9              | Good - slight boundary errors |
| 0.5-0.7              | Moderate - significant errors |
| <0.5                 | Poor - wrong window detected  |

### Delay Results

| Delay Accuracy | Interpretation                        |
| -------------- | ------------------------------------- |
| >95%           | Excellent - precise delay recovery    |
| 85-95%         | Good - mostly correct with few errors |
| 70-85%         | Moderate - systematic bias possible   |
| <70%           | Poor - unreliable delay estimates     |

## Next Steps

1. **Run Quick Validation**:

   ```bash
   python benchmark_lla_unified.py --n-values 40 --alpha-values 1.0 \
     --delay-values 0 --window-fractions 1.0 --n-replicates 20
   ```

   - Verify Power ≈ 100%, FPR ≈ 5%
2. **Test Localization**:

   ```bash
   python benchmark_lla_unified.py --n-values 40 --alpha-values 0.8 1.0 \
     --delay-values 0 --window-fractions 1.0 0.5 --n-replicates 30
   ```

   - Check window overlap metrics
3. **Test Delay Estimation**:

   ```bash
   python benchmark_lla_unified.py --n-values 60 --alpha-values 0.8 1.0 \
     --delay-values 0 2 3 --window-fractions 1.0 --n-replicates 30 \
     --delay-limit 4
   ```

   - Verify delay accuracy
4. **Full Benchmark**:

   - Run complete parameter sweep
   - Generate visualizations
   - Publish validation report

## Troubleshooting

**Low Power (<80%) for α=1.0**:

- Check normalization (should use 'pnz')
- Verify window constraints
- Increase precision (permutations)

**High FPR (>10%)**:

- Check p-value calculation
- Verify control group generation
- Review permutation implementation

**Poor Localization**:

- Ensure keep_trace=True
- Check Start_Z/End_Z parsing
- Verify window is within valid range

**Incorrect Delay**:

- Confirm delay_limit ≥ true delay
- Check that delay is within window
- Verify trace backtracking logic

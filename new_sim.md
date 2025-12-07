## Goals
Produce X, Y, Z triplets for detection/localization/delay benchmarking of LLA.
Keep baseline noise and injected signal cleanly separable.
Make SNR explicit and easy to tune.
Support continuous and discrete (Poisson, NegBin) extensions later.
Allow reproducible experiments (seeded RNG) and both experimental and control groups.
Core concepts and notation
n: sequence length (timepoints).
Z: regulation indicator (binary mask) — Z[i]=1 when position i is in the regulated window I, else 0.
Y_base: baseline Y sequence (noise, uncorrelated except where we optionally add signal).
X_base: baseline X sequence (noise).
effect_size (E): amplitude of the signal added to X in regulated window, relative to Y_base.
noise_sd: standard deviation of baseline noise for X and Y (choose noise_sd << effect amplitude).
delay_yz: integer delay between Y and Z (Y index offset relative to Z).
is_control: boolean — if True then do not add signal (control / null dataset).

default：n=100, delay=0

## Method A — Preferred (Baseline noise + additive signal)
Rationale: Clean separation of noise and signal; easier to tune SNR.

### Steps:
Generate Y_base from a target distribution:
Continuous: Y_base ~ Normal(0, sigma_Y) (e.g., sigma_Y = 1)
Discrete (future): Y_base ~ Poisson(lambda) or NegBinom
Generate X_base ~ Normal(0, noise_sd) with noise_sd << effect amplitude (e.g., noise_sd = 0.1)
Define Z as a binary mask for window [ws, we).
For i in 0..n-1:
Let y_idx = i + delay_yz, and note X is always align with Y(i.e. if delay exits, then X and Y aligns but delay Z). Now consider delay_yz = 0 first.
If is_control: X[i] = X_base[i] (no signal)
Else if Z[i] == 1 and y_idx in window (or within bounds):
X[i] = X_base[i] + E * Y_base[i]
Else:
X[i] = X_base[i]
Y[i] = Y_base[i] (or optionally add a small observation noise)
- Actually, in this case we have **X[i] = X_base[i] + E * Z[i]*Y_base[i], globally**.
vary E (0.5, 0,8, 1.0, 2.0) to sweep difficulty while keeping noise small.

## For future experiments with counts: (no need to achieve in the codes now)
Generate baseline counts, e.g., Y_base ~ Poisson(lambda_y)
X_base ~ Poisson(lambda_x)
Add signal by increasing counts inside window, e.g. X[i] = X_base[i] + Poisson(lambda_signal * Z[i] * g(Y_base[y_idx]))
Or use multiplicative factor: X[i] = Poisson( lambda_x + effect_size * Y_base[y_idx] )
Negative binomial extension: same idea but sample NB instead of Poisson.

## Control group (null)
Keep same Z mask and same marginal distributions for X and Y, but do not add signal:
X_control = X_base
Y_control = Y_base
Parameters to expose in generator API
n (int)
window_start, window_end (int or None for global)
delay_yz (int)
is_control (bool)
seed (int)
method (str): one of {"additive", "mixing", "count_poisson", "count_nb"}
effect_size / E (float)
noise_sd (float)
sigma_Y (float) — std of Y_base (defaults: sigma_Y=1)
count_args (dict) — for Poisson/NB parameters
Recommended default values for initial runs
n = 100
window_fraction = 0.8 → window_length = 80 (so window_start = 10, window_end = 90)
delay_yz = 0
is_control = 50 real associated and 50 noise samples. (maybe it should always produced same amount of real and noises so this para is not used anymore)
seed: deterministic for reproducibility
method = "additive"
sigma_Y = 1.0
effect_size = 0.8，0.82，...,0.9,0.92,...,1 
noise_sd = 0.1 (order-of-magnitude smaller than signal)
n_replicates = 50 for power curves
precision (permutations) = 1000 for p-values

## Tests to run after implementing
- make the result format same as `benchmark_lla_unified.py` generalized so can still used `visualize_alpha_sweep.py` to draw graph under new simulation methods.
- fpr,recall, time, iou, score, start/end error, inflation should be recorded.
- Adjust visualize_alpha_sweep.py defaults if needed for new parameter scales.

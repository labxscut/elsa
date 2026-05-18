import os
from pathlib import Path
#!/usr/bin/env python
"""
Unified triplet generator for LLA validation experiments.

Supports both global and local association patterns with:
- Variable association strength α ∈ [0, 1]
- Global regulation (window covers entire sequence)
- Local regulation (window is subsequence)
- Optional Y-Z delay (X-Y always synchronous)
- Control group (pure noise with same Z structure)
"""

import numpy as np
import argparse
import sys
import json
from typing import Optional, Tuple


def generate_unified_triplet(
    n: int,
    alpha: float = 1.0,
    window_start: Optional[int] = None,
    window_end: Optional[int] = None,
    delay_yz: int = 0,
    is_control: bool = False,
    seed: Optional[int] = None,
    method: str = 'mixing',
    effect_size: float = 1.0,
    noise_sd: float = 0.1,
    sigma_y: float = 1.0,
    z_encoding: str = '01'
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, dict]:
    """
    Generate (X, Y, Z) triplet with unified framework.
    
    Parameters
    ----------
    n : int
        Sequence length
    alpha : float
        Association strength in [0, 1]. Controls X = α*Y + (1-α)*ε (method='mixing' only)
    window_start : int, optional
        Start of regulation window (inclusive). If None, uses global window [0, n)
    window_end : int, optional
        End of regulation window (exclusive). If None, uses global window [0, n)
    delay_yz : int
        Y-Z delay (Y lags Z by this amount). X-Y are always synchronous.
    is_control : bool
        If True, generate control group (pure noise, no X-Y association)
    seed : int, optional
        Random seed for reproducibility
    method : str
        Generation method: 'mixing' (old alpha-based) or 'additive' (new signal+noise)
    effect_size : float
        Signal amplitude for additive method (default: 1.0)
    noise_sd : float
        Baseline noise standard deviation for additive method (default: 0.1)
    sigma_y : float
        Standard deviation of Y baseline (default: 1.0)
    z_encoding : str
        Encoding for Z regulator: '01' for binary 0/1, '-11' for binary -1/1 (default: '01')
        
    Returns
    -------
    x, y, z : np.ndarray
        Generated time series
        Note: Z is the regulation indicator:
              For z_encoding='01': Z[i] = 1 if i is in regulation window I, else 0
              For z_encoding='-11': Z[i] = 1 if i is in regulation window I, else -1
              Global case: Z = [1, 1, ..., 1] (or all 1 for -11)
              Local case: Z = [0/1, ..., 1, 1, ..., 0/1] or [-1/1, ..., 1, 1, ..., -1/1]
    metadata : dict
        Information about the generation parameters
    """
    if seed is not None:
        np.random.seed(seed)
    
    # Validate z_encoding
    if z_encoding not in ['01', '-11']:
        raise ValueError(f"z_encoding must be '01' or '-11', got {z_encoding}")
    
    # Define regulation window I
    if window_start is None or window_end is None:
        # Global case: entire sequence
        window_start = 0
        window_end = n
        is_global = True
    else:
        is_global = False
        # Validate window
        if window_start < 0 or window_end > n or window_start >= window_end:
            raise ValueError(f"Invalid window [{window_start}, {window_end}) for n={n}")
        
        window_length = window_end - window_start
        min_length = min(10, n // 2)
        if window_length < min_length:
            raise ValueError(f"Window length {window_length} < required minimum {min_length}")
    
    # Validate parameters based on method
    if method == 'mixing':
        if not (0 <= alpha <= 1):
            raise ValueError(f"alpha must be in [0, 1], got {alpha}")
    elif method == 'additive':
        if effect_size < 0:
            raise ValueError(f"effect_size must be >= 0, got {effect_size}")
        if noise_sd < 0:
            raise ValueError(f"noise_sd must be >= 0, got {noise_sd}")
    else:
        raise ValueError(f"method must be 'mixing' or 'additive', got {method}")
    
    # Initialize Z (regulation indicator) with proper encoding
    if z_encoding == '01':
        Z = np.zeros(n)
        Z[window_start:window_end] = 1.0
    else:  # z_encoding == '-11'
        Z = np.ones(n) * -1.0
        Z[window_start:window_end] = 1.0
    
    if method == 'additive':
        # New additive method: baseline noise + signal
        # Generate Y baseline
        Y_base = np.random.normal(0, sigma_y, size=n)
        # Generate X baseline noise
        X_base = np.random.normal(0, noise_sd, size=n)
        
        # Start with baseline
        Y = Y_base.copy()
        X = X_base.copy()
        
        if not is_control:
            # Add signal to X inside window: X[i] += effect_size * Y[i] * Z[i]
            # Handle delay: y_idx = i + delay_yz
            for i in range(n):
                if Z[i] == 1.0:  # Inside window
                    y_idx = i + delay_yz
                    # Check bounds
                    if 0 <= y_idx < n:
                        X[i] += effect_size * Y_base[y_idx]
    else:
        # Old mixing method (backward compatible)
        # Initialize with noise
        Y = np.random.normal(0, 1, size=n)      # Base variable
        epsilon = np.random.normal(0, 1, size=n)    # Noise for association
        epsilon_prime = np.random.normal(0, 1, size=n)  # Independent noise
        
        # Initialize X
        X = np.zeros(n)
        
        if is_control:
            # Control group: X is independent noise everywhere
            X = epsilon_prime.copy()
        else:
            # Experimental group: X depends on Y within window (with delay)
            for i in range(n):
                if window_start <= i < window_end:  # i in window I (Z_i = 1)
                    # Calculate Y index considering delay
                    y_idx = i + delay_yz
                    
                    # Check if y_idx is valid and also in regulation window
                    if 0 <= y_idx < n and window_start <= y_idx < window_end:
                        # X_i = α * Y_{i+delay} + (1-α) * ε_i
                        X[i] = alpha * Y[y_idx] + (1 - alpha) * epsilon[i]
                    else:
                        # Y index out of bounds or outside window: independent noise
                        X[i] = epsilon_prime[i]
                else:
                    # Outside window (Z_i = 0): independent noise
                    X[i] = epsilon_prime[i]
    
    # Prepare metadata
    metadata = {
        'n': n,
        'alpha': alpha,
        'method': method,
        'effect_size': effect_size,
        'noise_sd': noise_sd,
        'sigma_y': sigma_y,
        'z_encoding': z_encoding,
        'window_start': window_start,
        'window_end': window_end,
        'window_length': window_end - window_start,
        'is_global': is_global,
        'delay_yz': delay_yz,
        'is_control': is_control,
        'seed': seed
    }
    
    return X, Y, Z, metadata


def write_triplet_to_file(
    filepath: str,
    X: np.ndarray,
    Y: np.ndarray,
    Z: np.ndarray,
    labels: Tuple[str, str, str] = ('S1', 'S2', 'S3')
) -> None:
    """Write triplet to tab-delimited file compatible with lla_compute.py"""
    n = len(X)
    
    with open(filepath, 'w', encoding='utf-8') as f:
        # Header: # T1 T2 ... Tn
        header = ['#'] + [f'T{i+1}' for i in range(n)]
        f.write('\t'.join(header) + '\n')
        
        # Row for X
        row_x = [labels[0]] + [f'{v:.6f}' for v in X]
        f.write('\t'.join(row_x) + '\n')
        
        # Row for Y
        row_y = [labels[1]] + [f'{v:.6f}' for v in Y]
        f.write('\t'.join(row_y) + '\n')
        
        # Row for Z
        row_z = [labels[2]] + [f'{v:.6f}' for v in Z]
        f.write('\t'.join(row_z) + '\n')


def generate_and_save(
    out_dir: str,
    n: int,
    effect_size: float = 1.0,
    window_start: Optional[int] = None,
    window_end: Optional[int] = None,
    delay_yz: int = 0,
    is_control: bool = False,
    seed: Optional[int] = None,
    method: str = 'additive',
    noise_sd: float = 0.1,
    sigma_y: float = 1.0,
    z_encoding: str = '01'
) -> str:
    """Generate a triplet and save to out_dir with metadata JSON. Returns filepath."""
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    if seed is None:
        seed = np.random.randint(0, 100000000)
    X, Y, Z, metadata = generate_unified_triplet(
        n=n,
        alpha=1.0,
        method=method,
        effect_size=effect_size,
        noise_sd=noise_sd,
        sigma_y=sigma_y,
        z_encoding=z_encoding,
        window_start=window_start,
        window_end=window_end,
        delay_yz=delay_yz,
        is_control=is_control,
        seed=seed
    )
    # Compose filename
    kind = 'ctrl' if is_control else 'exp'
    window_tag = 'global' if metadata['is_global'] else f"w{metadata['window_start']}_{metadata['window_end']}"
    name = f"n{n}_E{effect_size}_d{delay_yz}_{window_tag}_{kind}_s{seed}.tsv"
    out_path = os.path.join(out_dir, name)
    write_triplet_to_file(out_path, X, Y, Z)
    # Save metadata
    meta_path = out_path + '.meta.json'
    with open(meta_path, 'w') as mf:
        json.dump(metadata, mf)
    return out_path


def main():
    parser = argparse.ArgumentParser(
        description='Generate unified triplets for LLA validation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Global association, α=0.8, no delay
  python gen_unified_triplets.py --n 40 --alpha 0.8 --out global_a08.txt
  
  # Local association, α=1.0, window [10, 30), no delay
  python gen_unified_triplets.py --n 40 --alpha 1.0 --window-start 10 --window-end 30 --out local_w1030.txt
  
  # Local with Y-Z delay=3
  python gen_unified_triplets.py --n 60 --alpha 0.6 --window-start 15 --window-end 45 --delay-yz 3 --out local_delay3.txt
  
  # Control group (pure noise)
  python gen_unified_triplets.py --n 40 --control --out control.txt
        """
    )
    
    parser.add_argument('--n', type=int, required=True,
                       help='Sequence length (e.g., 20, 40, 60, 80, 100)')
    parser.add_argument('--alpha', type=float, default=1.0,
                       help='Association strength in [0, 1] (default: 1.0, for method=mixing)')
    parser.add_argument('--method', type=str, default='additive', choices=['mixing', 'additive'],
                       help='Generation method (default: additive)')
    parser.add_argument('--effect-size', type=float, default=1.0,
                       help='Signal amplitude for additive method (default: 1.0)')
    parser.add_argument('--noise-sd', type=float, default=0.1,
                       help='Baseline noise SD for additive method (default: 0.1)')
    parser.add_argument('--sigma-y', type=float, default=1.0,
                       help='Y baseline standard deviation (default: 1.0)')
    parser.add_argument('--z-encoding', type=str, default='01', choices=['01', '-11'],
                       help='Z regulator encoding: 01 for 0/1 binary, -11 for -1/1 binary (default: 01)')
    parser.add_argument('--window-start', type=int, default=None,
                       help='Start of regulation window (inclusive). If not set, uses global window.')
    parser.add_argument('--window-end', type=int, default=None,
                       help='End of regulation window (exclusive). If not set, uses global window.')
    parser.add_argument('--delay-yz', type=int, default=0,
                       help='Y-Z delay: Y lags Z by this amount (default: 0)')
    parser.add_argument('--control', action='store_true',
                       help='Generate control group (pure noise, no association)')
    parser.add_argument('--seed', type=int, default=None,
                       help='Random seed (default: random)')
    parser.add_argument('--out', type=str, required=True,
                       help='Output file path')
    
    args = parser.parse_args()
    
    # Set default seed if not provided
    if args.seed is None:
        args.seed = np.random.randint(0, 1000000)
    
    try:
        # Generate triplet
        X, Y, Z, metadata = generate_unified_triplet(
            n=args.n,
            alpha=args.alpha,
            method=args.method,
            effect_size=args.effect_size,
            noise_sd=args.noise_sd,
            sigma_y=args.sigma_y,
            z_encoding=args.z_encoding,
            window_start=args.window_start,
            window_end=args.window_end,
            delay_yz=args.delay_yz,
            is_control=args.control,
            seed=args.seed
        )
        
        # Write to file
        write_triplet_to_file(args.out, X, Y, Z)
        
        # Print metadata
        print(f"Generated: {args.out}", file=sys.stderr)
        print(f"  Type: {'Control (null)' if metadata['is_control'] else 'Experimental'}", file=sys.stderr)
        print(f"  Mode: {'Global' if metadata['is_global'] else 'Local'}", file=sys.stderr)
        print(f"  Method: {metadata['method']}", file=sys.stderr)
        print(f"  Z encoding: {metadata['z_encoding']}", file=sys.stderr)
        print(f"  n = {metadata['n']}", file=sys.stderr)
        if not metadata['is_control']:
            if metadata['method'] == 'additive':
                print(f"  Effect size: {metadata['effect_size']}", file=sys.stderr)
                print(f"  Noise SD: {metadata['noise_sd']}", file=sys.stderr)
                print(f"  Sigma Y: {metadata['sigma_y']}", file=sys.stderr)
            else:
                print(f"  α = {metadata['alpha']}", file=sys.stderr)
        print(f"  Window: [{metadata['window_start']}, {metadata['window_end']}) (length={metadata['window_length']})", file=sys.stderr)
        print(f"  Delay Y-Z: {metadata['delay_yz']}", file=sys.stderr)
        print(f"  Seed: {metadata['seed']}", file=sys.stderr)
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()

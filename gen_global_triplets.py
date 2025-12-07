#!/usr/bin/env python
"""
LLA能否区分关联 vs 非关联三元组？

Generate the simplest possible triplets for LLA detection testing:
- global_assoc: entire series has X-Y correlation modulated by Z
- null: pure independent Gaussian noise

This is the SIMPLEST scenario to test if LLA can detect association at all.
"""

import numpy as np
import argparse
import sys
from typing import Tuple


def generate_global_assoc(n: int, seed: int = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Generate globally associated triplet: Z modulates X-Y across entire series.
    
    For each time point i:
        base = N(0,1)
        Y[i] = base
        X[i] = base * Z[i]
    where Z[i] ∈ {-1, +1} randomly
    
    This gives strong global LA signal across the whole series.
    """
    if seed is not None:
        np.random.seed(seed)
    
    z = np.random.choice([-1, 1], size=n)
    base = np.random.normal(0, 1, size=n)
    y = base.copy()
    x = base * z
    
    return x, y, z


def generate_null(n: int, seed: int = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Generate null triplet: pure independent Gaussian noise.
    
    X, Y, Z are all independent N(0,1).
    """
    if seed is not None:
        np.random.seed(seed)
    
    x = np.random.normal(0, 1, size=n)
    y = np.random.normal(0, 1, size=n)
    z = np.random.choice([-1, 1], size=n)  # Keep Z binary for consistency
    
    return x, y, z


def write_triplet_file(x: np.ndarray, y: np.ndarray, z: np.ndarray, 
                       output_path: str) -> None:
    """Write triplet to tab-delimited file for lla_compute."""
    n = len(x)
    
    with open(output_path, "w", encoding="utf-8") as f:
        # Header: # T1 T2 ... Tn (tab-separated)
        header = ["#"] + [f"T{i+1}" for i in range(n)]
        f.write("\t".join(header) + "\n")
        
        # S1: X (tab-separated)
        row_x = ["S1"] + [f"{v:.6f}" for v in x]
        f.write("\t".join(row_x) + "\n")
        
        # S2: Y (tab-separated)
        row_y = ["S2"] + [f"{v:.6f}" for v in y]
        f.write("\t".join(row_y) + "\n")
        
        # S3: Z (tab-separated)
        row_z = ["S3"] + [f"{float(v):.6f}" for v in z]
        f.write("\t".join(row_z) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Generate simplest triplets for LLA detection testing"
    )
    parser.add_argument("--out", dest="out_path", default="triplet.txt",
                       help="Output file path (default: triplet.txt)")
    parser.add_argument("--mode", dest="mode", required=True,
                       choices=["global_assoc", "null"],
                       help="Mode: global_assoc (associated) or null (no association)")
    parser.add_argument("--n", dest="n", type=int, default=40,
                       help="Number of time points (default: 40)")
    parser.add_argument("--seed", dest="seed", type=int, default=None,
                       help="Random seed (default: None)")
    
    args = parser.parse_args()
    
    # Generate data
    if args.mode == "global_assoc":
        x, y, z = generate_global_assoc(args.n, args.seed)
        print(f"Generated: {args.out_path} (global_assoc, n={args.n}, seed={args.seed})")
    else:
        x, y, z = generate_null(args.n, args.seed)
        print(f"Generated: {args.out_path} (null, n={args.n}, seed={args.seed})")
    
    # Write to file
    try:
        write_triplet_file(x, y, z, args.out_path)
    except OSError as e:
        print(f"Write failed: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

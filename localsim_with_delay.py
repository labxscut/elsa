import numpy as np
import argparse
import sys
from typing import Optional, Tuple

def generate_delayed_series(
    n: int,
    corr_start: int,
    corr_end: int,
    delay_xy: int = 0,
    delay_yz: int = 0,
    seed: Optional[int] = None
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    生成 (X, Y, Z) 三元组，满足：
      - Z 是固定的 ±1 序列
      - 仅当 k ∈ [corr_start, corr_end] 时，Z[k] 作为调控源
      - 它调控：
          Y[k + delay_yz] = base
          X[k + delay_yz + delay_xy] = base * Z[k]
      - 但要求：k ∈ [corr_start, corr_end] 
                AND y_idx ∈ [corr_start, corr_end]
                AND x_idx ∈ [corr_start, corr_end]
      - 所有未被调控的位置为独立高斯噪声
    """
    if seed is not None:
        np.random.seed(seed)

    # 初始化
    z = np.random.choice([-1, 1], size=n)
    x = np.random.normal(0, 1, size=n)
    y = np.random.normal(0, 1, size=n)

    total_delay = delay_yz + delay_xy

    for k in range(corr_start, corr_end):
        y_idx = k + delay_yz
        x_idx = k + total_delay

        # 检查所有索引是否在有效范围内
        if (0 <= y_idx < n and 
            0 <= x_idx < n and
            corr_start <= y_idx <= corr_end and
            corr_start <= x_idx <= corr_end):
            
            base = np.random.normal(0, 1)
            y[y_idx] = base
            x[x_idx] = base * z[k]

    return x, y, z


def main():
    parser = argparse.ArgumentParser(description="Generate local-sim data with delays for LLA computation")
    parser.add_argument("--out", dest="out_path", default="lla_input_delay.txt",
                        help="Output file path (tab-delimited, for lla_compute)")
    parser.add_argument("--n", dest="n", type=int, default=40,
                        help="Number of time points (default: 40)")
    parser.add_argument("--delay_xy", dest="delay_xy", type=int, default=0,
                        help="Delay: X lags Y by this many steps (default: 0)")
    parser.add_argument("--delay_yz", dest="delay_yz", type=int, default=0,
                        help="Delay: Y lags Z by this many steps (default: 0)")
    parser.add_argument("--corr-start", dest="corr_start", type=int, default=5,
                        help="Start index of correlated region (inclusive, default: 5)")
    parser.add_argument("--corr-end", dest="corr_end", type=int, default=15,
                        help="End index of correlated region (inclusive, default: 15)")
    parser.add_argument("--seed", dest="seed", type=int, default=33,
                        help="Random seed (default: 33)")
    args = parser.parse_args()

    # Validate corr range
    if args.corr_start > args.corr_end or args.corr_start < 0 or args.corr_end >= args.n:
        print("Error: Invalid corr_start/corr_end range.", file=sys.stderr)
        sys.exit(1)

    # Generate data
    x, y, z = generate_delayed_series(
        n=args.n,
        corr_start=args.corr_start,
        corr_end=args.corr_end,
        delay_xy=args.delay_xy,
        delay_yz=args.delay_yz,
        seed=args.seed
    )

    # Write to file (once!)
    try:
        with open(args.out_path, "w", encoding="utf-8") as f:
            # Header: # T1 T2 ... Tn
            header = ["#"] + [f"T{i+1}" for i in range(args.n)]
            f.write("\t".join(header) + "\n")

            # S1: X
            row_x = ["S1"] + [f"{v:.6f}" for v in x]
            f.write("\t".join(row_x) + "\n")

            # S2: Y
            row_y = ["S2"] + [f"{v:.6f}" for v in y]
            f.write("\t".join(row_y) + "\n")

            # S3: Z (as float -1.000000 / 1.000000)
            row_z = ["S3"] + [f"{float(v):.6f}" for v in z]
            f.write("\t".join(row_z) + "\n")
    except OSError as e:
        print(f"Write failed: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Generated: {args.out_path}")
    print(f"Delays: Y lags Z by {args.delay_yz}, X lags Y by {args.delay_xy} "
          f"(so X lags Z by {args.delay_yz + args.delay_xy})")
    print(f"Correlated region (for Z): [{args.corr_start}, {args.corr_end}]")
    print(f"Effective Y modulation range: [{args.corr_start + args.delay_yz}, {args.corr_end + args.delay_yz}] ∩ [{args.corr_start}, {args.corr_end}]")
    print(f"Effective X modulation range: [{args.corr_start + args.delay_yz + args.delay_xy}, {args.corr_end + args.delay_yz + args.delay_xy}] ∩ [{args.corr_start}, {args.corr_end}]")


if __name__ == "__main__":
    main()
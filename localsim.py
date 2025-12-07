import numpy as np
import argparse
import sys


def main():
    # CLI 参数（不改变模拟逻辑，仅控制规模与输出位置）
    parser = argparse.ArgumentParser(description="Generate local-sim data for lla_compute")
    parser.add_argument("--out", dest="out_path", default="lla_input.txt",
                        help="输出文件路径（制表符分隔，lla_compute可直接读取）")
    parser.add_argument("--n", dest="n", type=int, default=20,
                        help="时间点数量（默认 20）")
    parser.add_argument("--seed", dest="seed", type=int, default=33,
                        help="随机种子（默认 33）")
    args = parser.parse_args()

    # Set random seed for reproducibility
    np.random.seed(args.seed)

    # Generate z as {-1, +1}
    n = args.n
    z = np.random.choice([-1, 1], size=n)

    # generate x and y depending on z （no noise, local)
    x, y = np.zeros(n), np.zeros(n)  # 初始化x和y数组

    corr_start, corr_end = 10, 20  # 在以外的区间非随机（保持逻辑不变）

    for i in range(n):
        if corr_start <= i <= corr_end:
            # 区间内：使用版本1的正负相关逻辑（无噪声）
            if z[i] == 1:
                base = np.random.normal(0, 1)
                x[i] = base
                y[i] = base  # 正相关
            else:
                base = np.random.normal(0, 1)
                x[i] = base
                y[i] = -base  # 负相关
        else:
            # 区间外：x和y完全随机，无相关性
            x[i] = np.random.normal(0, 1)  # 随机值
            y[i] = np.random.normal(0, 1)  # 随机值（与x无关）

    # 输出为 lla_compute 需要的制表符文件格式：
    # 第一行：# \t T1 ... Tn
    # 后续三行：S1, S2, S3 对应 x, y, z
    try:
        with open(args.out_path, "w", encoding="utf-8") as f:
            # Header
            header = ["#"] + [f"T{i+1}" for i in range(n)]
            f.write("\t".join(header) + "\n")

            # Row S1 (x)
            row_x = ["S1"] + [f"{v:.6f}" for v in x]
            f.write("\t".join(row_x) + "\n")

            # Row S2 (y)
            row_y = ["S2"] + [f"{v:.6f}" for v in y]
            f.write("\t".join(row_y) + "\n")

            # Row S3 (z) as -1/1 with fixed decimals
            row_z = ["S3"] + [f"{float(v):.6f}" for v in z]
            f.write("\t".join(row_z) + "\n")
    except OSError as e:
        print(f"写入失败: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"已生成: {args.out_path}")


if __name__ == "__main__":
    main()
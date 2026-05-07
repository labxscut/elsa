#!/usr/bin/env python3

"""GPU smoke test for the LLA backend."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

import numpy as np

from lsa import compcore


def load_first_triplet_series(tsv_path: Path) -> tuple[list[str], np.ndarray]:
    with tsv_path.open("r", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        labels = []
        numeric_rows = []
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            try:
                numeric_rows.append([
                    float(value) if value not in {"NA", "na", ""} else np.nan
                    for value in row[1:]
                ])
            except ValueError:
                continue
            labels.append(row[0])

    values = np.array(numeric_rows, dtype=float)
    return labels, values


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", type=Path)
    parser.add_argument("--rep-num", type=int, default=1)
    parser.add_argument("--spot-num", type=int, default=None)
    parser.add_argument("--delay", type=int, default=0)
    args = parser.parse_args()

    labels, raw = load_first_triplet_series(args.input_file)
    if args.spot_num is None:
        spot_num = raw.shape[1] // args.rep_num
    else:
        spot_num = args.spot_num

    clean = np.ma.zeros((raw.shape[0], args.rep_num, spot_num), dtype=float)
    for i in range(raw.shape[0]):
        for j in range(args.rep_num):
            series = raw[i][j::args.rep_num][:spot_num]
            clean[i, j] = np.ma.array(series, mask=np.isnan(series))

    if raw.shape[0] < 3:
        raise SystemExit("Need at least 3 rows for LLA smoke test")

    x = clean[0]
    y = clean[1]
    z = clean[2]

    lla_data = compcore.LLA_Data(args.delay, x.tolist(), y.tolist(), z.tolist())
    lla_result = compcore.DP_lla(lla_data, True)

    print(f"backend={getattr(compcore, 'BACKEND_NAME', 'unknown')}")
    print(f"using_gpu={getattr(compcore, 'USING_GPU', False)}")
    print(f"score={lla_result.score}")
    print(f"trace={lla_result.trace}")


if __name__ == "__main__":
    main()

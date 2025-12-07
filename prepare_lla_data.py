#!/usr/bin/env python3
"""
Generate simulated data and transform it into LLA-compatible input format.
"""

import numpy as np
import os

def generate_simulation_data(n=20, corr_start=3, corr_end=10):
    """Generate simulated X, Y, Z sequences."""
    np.random.seed(33)

    # Z: ±1 alternating random sequence
    z = np.random.randint(0, 2, size=n) * 2 - 1

    # X and Y
    x, y = np.zeros(n), np.zeros(n)

    for i in range(n):
        if corr_start <= i <= corr_end:
            # Local correlation region
            base = np.random.normal(0, 1)
            if z[i] == 1:
                x[i] = base
                y[i] = base  # positive correlation
            else:
                x[i] = base
                y[i] = -base  # negative correlation
        else:
            # Outside region: uncorrelated random noise
            x[i] = np.random.normal(0, 1)
            y[i] = np.random.normal(0, 1)

    return x, y, z, n

def create_lla_input_file(x, y, z, n, filename="lla_input.txt"):
    """Create tab-delimited LLA input file."""
    time_labels = [f"T{i+1}" for i in range(n)]
    header = "#Factor\t" + "\t".join(time_labels)

    lines = []
    for label, arr in zip(["X", "Y", "Z"], [x, y, z]):
        arr_str = "\t".join(f"{val:.6f}" for val in arr)
        lines.append(f"{label}\t{arr_str}")

    with open(filename, "w") as f:
        f.write(header + "\n" + "\n".join(lines) + "\n")

    print(f"✅ LLA input file created: {filename}")
    print(f"   Shape: 3 factors × {n} time points")
    return filename

def create_lla_script(input_file, output_file="lla_results.txt"):
    """Generate a bash script to run LLA analysis."""
    script_content = f"""#!/bin/bash
# Auto-generated LLA run script
echo "Running LLA analysis..."

python3 lla/lla_compute.py {input_file} {output_file} \\
    -d 3 -p perm -x 1000 -r 1 -s 20 -m 0.5 -t simple -f linear -n pnz

echo "LLA analysis done. Results: {output_file}"
"""
    with open("run_lla_analysis.sh", "w") as f:
        f.write(script_content)
    os.chmod("run_lla_analysis.sh", 0o755)
    print("✅ Script created: run_lla_analysis.sh")

def main():
    print("=== Preparing LLA data ===")

    # 1. Generate data
    x, y, z, n = generate_simulation_data()

    # 2. Preview data
    print(f"\nPreview:")
    print(f"X[:5]: {x[:5]}")
    print(f"Y[:5]: {y[:5]}")
    print(f"Z[:5]: {z[:5]}")

    # 3. Create input file
    input_file = create_lla_input_file(x, y, z, n)

    # 4. Create LLA script
    create_lla_script(input_file)

    # 5. Show preview of file content
    print("\n=== Preview of generated file ===")
    with open(input_file) as f:
        for i, line in enumerate(f.readlines()[:4]):
            print(f"{i+1}: {line.strip()}")

if __name__ == "__main__":
    main()


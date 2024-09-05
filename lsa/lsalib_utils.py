import subprocess
import os

def safeCmd(cmd):
    try:
        return subprocess.check_output(cmd, shell=True, universal_newlines=True).strip()
    except subprocess.CalledProcessError as e:
        return f"Error: {e}"

# Add any other utility functions here that were in the original lsalib.py
# but not included in other new modules
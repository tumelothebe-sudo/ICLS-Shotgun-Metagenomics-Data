#!/usr/bin/env python3
"""Reproduce the two main Data in Brief non-normalized top-12 bar plots."""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = [
    ROOT / "scripts" / "species.py",
    ROOT / "scripts" / "kegg_l3.py",
]

def main():
    for script in SCRIPTS:
        print(f"Running {script.name} ...", flush=True)
        subprocess.run([sys.executable, str(script)], cwd=ROOT, check=True)
    print("All figures saved in:", ROOT / "output" / "figures")

if __name__ == "__main__":
    main()

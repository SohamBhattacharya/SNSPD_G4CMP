"""
load_phonon_steps.py
────────────────────
Diagnostic script — run this first to inspect your phonon_steps.txt
and understand its format before running the animation scripts.

Usage:
    python load_phonon_steps.py --input phonon_steps.txt
"""
import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', default='phonon_steps.txt')
    args = parser.parse_args()

    path = Path(args.input)
    if not path.exists():
        print(f"ERROR: {args.input} does not exist")
        return

    print(f"File: {args.input}")
    print(f"Size: {path.stat().st_size} bytes")
    print()

    lines = path.read_text(errors='replace').splitlines()
    print(f"Total lines: {len(lines)}")
    print()

    print("=== First 20 lines ===")
    for i, line in enumerate(lines[:20]):
        print(f"  [{i:3d}] {repr(line)}")
    print()

    # Count non-empty, non-comment lines
    data_lines = [l for l in lines if l.strip() and not l.strip().startswith('#')]
    print(f"Non-header lines: {len(data_lines)}")
    if data_lines:
        parts = data_lines[0].split()
        print(f"First data line has {len(parts)} columns: {parts}")

if __name__ == '__main__':
    main()

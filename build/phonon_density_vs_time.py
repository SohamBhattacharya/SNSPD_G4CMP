"""
phonon_density_vs_time.py
─────────────────────────
Reads RISQTutorial_hits.txt produced by G4CMP and plots
phonon hit density vs time for a user-defined spatial region of the slab.

Usage:
    python phonon_density_vs_time.py \
        --hits /path/to/RISQTutorial_hits.txt \
        --xmin -4 --xmax 4 \
        --ymin -4 --ymax 4 \
        --tmax 200 --bins 100

Units in the hits file are mm and ns (Geant4 default output).
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pathlib import Path


# ── Parse hits file ──────────────────────────────────────────────────────────

def load_hits(filepath):
    """
    RISQTutorial_hits.txt columns (space separated):
        EventID  TrackID  ParticleName  Edep(eV)  X(mm)  Y(mm)  Z(mm)  T(ns)

    Returns a dict of arrays keyed by column name.
    Skips header lines starting with '#' or non-numeric first token.
    """
    rows = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 8:
                continue
            try:
                event  = int(parts[0])
                track  = int(parts[1])
                name   = parts[2]
                edep   = float(parts[3])
                x      = float(parts[4])
                y      = float(parts[5])
                z      = float(parts[6])
                t      = float(parts[7])
                rows.append((event, track, name, edep, x, y, z, t))
            except ValueError:
                continue   # header or malformed line

    if not rows:
        raise RuntimeError(f"No valid hits found in {filepath}. "
                           "Check column order or run beamOn first.")

    rows = np.array(rows, dtype=object)
    return {
        'event':  rows[:, 0].astype(int),
        'track':  rows[:, 1].astype(int),
        'name':   rows[:, 2],
        'edep':   rows[:, 3].astype(float),
        'x':      rows[:, 4].astype(float),
        'y':      rows[:, 5].astype(float),
        'z':      rows[:, 6].astype(float),
        't':      rows[:, 7].astype(float),
    }


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Phonon hit density vs time")
    parser.add_argument('--hits',  default='RISQTutorial_hits.txt',
                        help='Path to hits file')
    # Spatial region (mm). Defaults = full 8x8 mm slab face.
    parser.add_argument('--xmin', type=float, default=-4.0)
    parser.add_argument('--xmax', type=float, default= 4.0)
    parser.add_argument('--ymin', type=float, default=-4.0)
    parser.add_argument('--ymax', type=float, default= 4.0)
    # Time axis
    parser.add_argument('--tmin', type=float, default=0.0,
                        help='Start time (ns)')
    parser.add_argument('--tmax', type=float, default=200.0,
                        help='End time (ns)')
    parser.add_argument('--bins', type=int,   default=100,
                        help='Number of time bins')
    # Particle filter
    parser.add_argument('--particle', default='all',
                        help='Particle name filter: all | phononL | phononTF | phononTS')
    parser.add_argument('--output', default=None,
                        help='Save figure to file instead of showing it')
    args = parser.parse_args()

    # ── Load ──────────────────────────────────────────────────────────────────
    hits = load_hits(args.hits)
    print(f"Loaded {len(hits['t'])} hits from {args.hits}")
    print(f"  Time range in file: {hits['t'].min():.2f} – {hits['t'].max():.2f} ns")
    print(f"  Particles: {set(hits['name'])}")

    # ── Spatial filter ────────────────────────────────────────────────────────
    mask = (
        (hits['x'] >= args.xmin) & (hits['x'] <= args.xmax) &
        (hits['y'] >= args.ymin) & (hits['y'] <= args.ymax)
    )

    # ── Particle filter ───────────────────────────────────────────────────────
    if args.particle != 'all':
        mask &= (hits['name'] == args.particle)

    t_sel = hits['t'][mask]
    print(f"  Hits in region after filter: {mask.sum()}")

    if mask.sum() == 0:
        print("No hits match the spatial/particle filter. "
              "Check --xmin/xmax/ymin/ymax and --particle.")
        return

    # ── Histogram in time ─────────────────────────────────────────────────────
    bin_edges = np.linspace(args.tmin, args.tmax, args.bins + 1)
    counts, _ = np.histogram(t_sel, bins=bin_edges)
    bin_centres = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_width   = bin_edges[1] - bin_edges[0]   # ns

    # Region area in mm²
    area_mm2 = (args.xmax - args.xmin) * (args.ymax - args.ymin)

    # Density = hits / (mm² · ns)
    density = counts / (area_mm2 * bin_width)

    # ── Plot ──────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 1, figsize=(9, 7),
                              gridspec_kw={'height_ratios': [3, 1]},
                              sharex=True)
    fig.suptitle('Phonon hit density vs time\n'
                 f'Region: X=[{args.xmin},{args.xmax}] mm, '
                 f'Y=[{args.ymin},{args.ymax}] mm  |  '
                 f'Particle: {args.particle}',
                 fontsize=11)

    # Top: density
    ax = axes[0]
    ax.fill_between(bin_centres, density, alpha=0.3, color='royalblue')
    ax.step(bin_centres, density, where='mid', color='royalblue', linewidth=1.2)
    ax.set_ylabel('Hit density  [hits mm⁻² ns⁻¹]')
    ax.set_xlim(args.tmin, args.tmax)
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    ax.grid(True, alpha=0.3)

    # Bottom: cumulative hits
    ax2 = axes[1]
    cumulative = np.cumsum(counts)
    ax2.plot(bin_centres, cumulative, color='darkorange', linewidth=1.2)
    ax2.set_ylabel('Cumulative hits')
    ax2.set_xlabel('Time  [ns]')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    if args.output:
        plt.savefig(args.output, dpi=150)
        print(f"Saved to {args.output}")
    else:
        plt.show()


if __name__ == '__main__':
    main()

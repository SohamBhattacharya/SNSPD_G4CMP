"""
phonon_radial_density.py
────────────────────────
X-axis: radial distance r from a reference point (mm)
Y-axis: time (ns)
Colour: phonon step density (steps / mm² / ns)

Usage:
    python phonon_radial_density.py --input phonon_steps.txt
    python phonon_radial_density.py --input phonon_steps.txt \\
        --cx 0 --cy 0 --rmax 7 --nr 80 --tmax 200 --nt 200 --log
"""

import argparse, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pathlib import Path

# ── shared loader (17-col format) ─────────────────────────────────────────────
PTYPE_MAP  = {'phononL': 1, 'phononTF': 2, 'phononTS': 3}
PTYPE_NAME = {1: 'phononL', 2: 'phononTF', 3: 'phononTS'}
PTYPE_CMAP = {'phononL': 'Blues', 'phononTF': 'Greens', 'phononTS': 'Reds'}

def load_steps(filepath):
    rows = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            p = line.split()
            if len(p) < 17:
                continue
            try:
                rows.append({'run': int(p[0]), 'event': int(p[1]),
                             'trackID': int(p[2]), 'particle': p[3],
                             'ptype': PTYPE_MAP.get(p[3], 0),
                             'preX': float(p[4]),  'preY': float(p[5]),
                             'preZ': float(p[6]),
                             'preKE': float(p[8]),
                             'postX': float(p[9]), 'postY': float(p[10]),
                             'postZ': float(p[11]),
                             'postKE': float(p[13]),
                             'preT': float(p[14]), 'postT': float(p[15]),
                             'process': p[16]})
            except (ValueError, IndexError):
                continue
    if not rows:
        raise RuntimeError("No valid rows found.")
    df = pd.DataFrame(rows)
    return df[df['ptype'] > 0]

# ── synthetic ─────────────────────────────────────────────────────────────────
def make_synthetic(seed=42):
    rng = np.random.default_rng(seed)
    rows = []
    for tid in range(1, 81):
        ptype = rng.choice([1,2,3])
        x, y = rng.uniform(-0.5,0.5,2)
        z    = rng.uniform(-0.25,0.25)
        v    = rng.normal(size=3); v /= np.linalg.norm(v)
        t    = 0.0
        for _ in range(400):
            dl = rng.exponential(0.5); dt = dl/5.0
            nx,ny,nz = x+v[0]*dl, y+v[1]*dl, z+v[2]*dl
            for i,(lo,hi) in enumerate([(-4,4),(-4,4),(-0.2625,0.2625)]):
                coord = [nx,ny,nz][i]
                if coord > hi: coord = 2*hi-coord; v[i]=-v[i]
                if coord < lo: coord = 2*lo-coord; v[i]=-v[i]
                if i==0: nx=coord
                elif i==1: ny=coord
                else: nz=coord
            rows.append({'trackID':tid,'ptype':ptype,
                         'preX':x,'preY':y,'preT':t,
                         'postX':nx,'postY':ny,'postT':t+dt})
            x,y,z,t = nx,ny,nz,t+dt
    return pd.DataFrame(rows)

# ── radial density ────────────────────────────────────────────────────────────
def radial_density(df, cx, cy, nr, rmax, nt, tmin, tmax):
    xm = 0.5*(df['preX'].values  + df['postX'].values)
    ym = 0.5*(df['preY'].values  + df['postY'].values)
    tm = 0.5*(df['preT'].values  + df['postT'].values)
    r  = np.sqrt((xm-cx)**2 + (ym-cy)**2)

    r_edges = np.linspace(0,    rmax, nr+1)
    t_edges = np.linspace(tmin, tmax, nt+1)
    dt      = t_edges[1]-t_edges[0]

    h, _, _ = np.histogram2d(r, tm, bins=[r_edges, t_edges],
                              range=[[0,rmax],[tmin,tmax]])

    # normalise by annular area × time bin
    r_lo, r_hi  = r_edges[:-1], r_edges[1:]
    annulus     = np.pi*(r_hi**2 - r_lo**2)
    density     = h / (annulus[:,np.newaxis] * dt)
    return density, r_edges, t_edges

# ── plot ──────────────────────────────────────────────────────────────────────
def plot(df, cx, cy, nr, rmax, nt, tmin, tmax, log_scale, output):
    modes = {PTYPE_NAME[p]: df[df['ptype']==p]
             for p in sorted(PTYPE_NAME) if (df['ptype']==p).any()}
    n_panels = 1 + len(modes)

    fig, axes = plt.subplots(1, n_panels,
                              figsize=(5.5*n_panels, 5.5),
                              constrained_layout=True)
    if n_panels == 1:
        axes = [axes]

    fig.suptitle(
        f'Phonon radial density  |  source ({cx:.1f}, {cy:.1f}) mm',
        fontsize=12, fontweight='bold')

    def make_panel(ax, density, r_edges, t_edges, title, cmap,
                   draw_speeds=False):
        vmin = density[density>0].min() if density.max()>0 else 1e-9
        vmax = density.max() if density.max()>0 else 1.0
        norm = mcolors.LogNorm(vmin=vmin, vmax=vmax) if log_scale \
               else mcolors.Normalize(vmin=0, vmax=vmax)

        pcm = ax.pcolormesh(r_edges, t_edges, density.T,
                            cmap=cmap, norm=norm, shading='flat')
        fig.colorbar(pcm, ax=ax, label='Steps mm⁻² ns⁻¹', pad=0.02)
        ax.set_xlabel('Radial distance  r  [mm]', fontsize=10)
        ax.set_ylabel('Time  [ns]', fontsize=10)
        ax.set_title(title, fontsize=10)
        ax.set_xlim(0, rmax)
        ax.set_ylim(tmin, tmax)

        if draw_speeds:
            speeds = {'$v_L$=8.97': 8.97, '$v_{TF}$=5.77': 5.77,
                      '$v_{TS}$=3.35': 3.35}
            cols   = ['#aec7e8','#98df8a','#ff9896']
            t_arr  = np.array([tmin, tmax])
            for (lbl,v),c in zip(speeds.items(), cols):
                ax.plot(v*t_arr, t_arr, '--', color=c, lw=1.0,
                        alpha=0.9, label=f'{lbl} mm/ns')
            ax.legend(fontsize=7, loc='upper left',
                      framealpha=0.3, labelcolor='white')

    # total
    dens, r_edges, t_edges = radial_density(df, cx, cy, nr, rmax, nt, tmin, tmax)
    make_panel(axes[0], dens, r_edges, t_edges,
               'All phonons', 'inferno', draw_speeds=True)

    # per mode
    for k, (name, sub) in enumerate(modes.items()):
        d, _, _ = radial_density(sub, cx, cy, nr, rmax, nt, tmin, tmax)
        make_panel(axes[k+1], d, r_edges, t_edges,
                   name, PTYPE_CMAP[name], draw_speeds=False)

    plt.savefig(output, dpi=150, bbox_inches='tight')
    print(f"Saved → {output}")
    plt.show()

# ── main ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input',     default='phonon_steps.txt')
    parser.add_argument('--cx',        type=float, default=0.0)
    parser.add_argument('--cy',        type=float, default=0.0)
    parser.add_argument('--rmax',      type=float, default=None)
    parser.add_argument('--nr',        type=int,   default=80)
    parser.add_argument('--tmin',      type=float, default=0.0)
    parser.add_argument('--tmax',      type=float, default=None)
    parser.add_argument('--nt',        type=int,   default=200)
    parser.add_argument('--log',       action='store_true')
    parser.add_argument('--output',    default='phonon_radial_density.png')
    parser.add_argument('--synthetic', action='store_true')
    args = parser.parse_args()

    if args.synthetic or not Path(args.input).exists():
        print("Using synthetic data ...")
        df = make_synthetic()
    else:
        print(f"Loading {args.input} ...")
        df = load_steps(args.input)

    print(f"  {len(df)} steps | particles: {df['ptype'].value_counts().to_dict()}")
    tmax = args.tmax or float(df['postT'].max())
    rmax = args.rmax or float(
        np.sqrt(max(abs(4.0-args.cx), abs(-4.0-args.cx))**2 +
                max(abs(4.0-args.cy), abs(-4.0-args.cy))**2) * 1.05)

    print(f"  t=[{args.tmin:.1f}, {tmax:.1f}] ns  |  r=[0, {rmax:.2f}] mm")
    plot(df, args.cx, args.cy, args.nr, rmax, args.nt,
         args.tmin, tmax, args.log, args.output)

if __name__ == '__main__':
    main()

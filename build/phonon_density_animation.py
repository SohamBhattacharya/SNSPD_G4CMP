"""
phonon_density_animation.py
───────────────────────────
Animated GIF of phonon XY density on a discretised mesh vs time.

Usage:
    python phonon_density_animation.py --input phonon_steps.txt
    python phonon_density_animation.py --input phonon_steps.txt \\
        --nx 40 --ny 40 --tmax 200 --nframes 80 --output anim.gif
"""

import argparse, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LogNorm, Normalize
from pathlib import Path

# ── shared loader ─────────────────────────────────────────────────────────────
PTYPE_MAP  = {'phononL': 1, 'phononTF': 2, 'phononTS': 3}
PTYPE_NAME = {1: 'phononL', 2: 'phononTF', 3: 'phononTS'}
PTYPE_CMAP = {'phononL': 'Blues', 'phononTF': 'Greens', 'phononTS': 'Reds'}
SLAB_X, SLAB_Y = (-4.0, 4.0), (-4.0, 4.0)

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
                             'preZ': float(p[6]),  'preE': float(p[7]),
                             'preKE': float(p[8]),
                             'postX': float(p[9]), 'postY': float(p[10]),
                             'postZ': float(p[11]),'postE': float(p[12]),
                             'postKE': float(p[13]),
                             'preT': float(p[14]), 'postT': float(p[15]),
                             'process': p[16]})
            except (ValueError, IndexError):
                continue
    if not rows:
        raise RuntimeError("No valid rows found.")
    df = pd.DataFrame(rows)
    return df[df['ptype'] > 0]

# ── synthetic data ────────────────────────────────────────────────────────────
def make_synthetic(n_tracks=50, steps_per_track=300, seed=42):
    rng = np.random.default_rng(seed)
    rows = []
    for tid in range(1, n_tracks+1):
        ptype = rng.choice([1,2,3])
        x, y  = rng.uniform(-0.5,0.5,2)
        z     = rng.uniform(-0.25,0.25)
        v     = rng.normal(size=3); v /= np.linalg.norm(v)
        t     = 0.0
        for _ in range(steps_per_track):
            dl = rng.exponential(0.5); dt = dl/5.0
            nx,ny,nz = x+v[0]*dl, y+v[1]*dl, z+v[2]*dl
            if nx> 4.0: nx= 8.0-nx; v[0]=-v[0]
            if nx<-4.0: nx=-8.0-nx; v[0]=-v[0]
            if ny> 4.0: ny= 8.0-ny; v[1]=-v[1]
            if ny<-4.0: ny=-8.0-ny; v[1]=-v[1]
            if nz> 0.2625: nz= 0.525-nz; v[2]=-v[2]
            if nz<-0.2625: nz=-0.525-nz; v[2]=-v[2]
            rows.append({'trackID':tid,'ptype':ptype,
                         'preX':x,'preY':y,'preT':t,
                         'postX':nx,'postY':ny,'postT':t+dt})
            x,y,z,t = nx,ny,nz,t+dt
    return pd.DataFrame(rows)

# ── frame builder ─────────────────────────────────────────────────────────────
def build_frames(df, nx, ny, t_edges):
    xedges = np.linspace(SLAB_X[0], SLAB_X[1], nx+1)
    yedges = np.linspace(SLAB_Y[0], SLAB_Y[1], ny+1)
    xm = 0.5*(df['preX'].values + df['postX'].values)
    ym = 0.5*(df['preY'].values + df['postY'].values)
    tm = 0.5*(df['preT'].values + df['postT'].values)
    frames = np.zeros((len(t_edges)-1, nx, ny))
    for i in range(len(t_edges)-1):
        mask = (tm >= t_edges[i]) & (tm < t_edges[i+1])
        if mask.sum():
            h,_,_ = np.histogram2d(xm[mask], ym[mask], bins=[xedges,yedges])
            frames[i] = h
    return frames, xedges, yedges

# ── animation ─────────────────────────────────────────────────────────────────
def animate(df, nx, ny, tmin, tmax, nframes, output):
    t_edges = np.linspace(tmin, tmax, nframes+1)
    t_centres = 0.5*(t_edges[:-1]+t_edges[1:])

    frames_all, xedges, yedges = build_frames(df, nx, ny, t_edges)
    mode_frames = {}
    for ptype, name in PTYPE_NAME.items():
        sub = df[df['ptype']==ptype]
        if not sub.empty:
            mode_frames[name], _, _ = build_frames(sub, nx, ny, t_edges)

    n_cols  = 1 + len(mode_frames)
    extent  = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    vmax_all = max(frames_all.max(), 1)
    vmaxes   = {n: max(f.max(),1) for n,f in mode_frames.items()}

    fig, axes = plt.subplots(1, n_cols, figsize=(4.8*n_cols, 4.8),
                              constrained_layout=True)
    if n_cols == 1:
        axes = [axes]

    ims = []
    # total panel
    im = axes[0].imshow(frames_all[0].T, origin='lower', extent=extent,
                        aspect='equal', cmap='inferno',
                        vmin=0, vmax=vmax_all, interpolation='gaussian')
    fig.colorbar(im, ax=axes[0], label='Steps / cell')
    axes[0].set_title('All phonons')
    axes[0].set_xlabel('X [mm]'); axes[0].set_ylabel('Y [mm]')
    ims.append(im)

    # per-mode panels
    for k,(name,fdata) in enumerate(mode_frames.items()):
        ax = axes[k+1]
        im = ax.imshow(fdata[0].T, origin='lower', extent=extent,
                       aspect='equal', cmap=PTYPE_CMAP[name],
                       vmin=0, vmax=vmaxes[name], interpolation='gaussian')
        fig.colorbar(im, ax=ax, label='Steps / cell')
        ax.set_title(name); ax.set_xlabel('X [mm]')
        ims.append(im)

    ttxt = axes[0].text(0.02,0.97,'', transform=axes[0].transAxes,
                         color='white', fontsize=9, va='top',
                         fontfamily='monospace')

    def update(i):
        ims[0].set_data(frames_all[i].T)
        for k,(_,fdata) in enumerate(mode_frames.items()):
            ims[k+1].set_data(fdata[i].T)
        ttxt.set_text(f't = {t_centres[i]:.1f} ns')
        return ims + [ttxt]

    ani = animation.FuncAnimation(fig, update, frames=nframes,
                                   interval=80, blit=True)
    suffix = Path(output).suffix.lower()
    if suffix == '.mp4':
        ani.save(output, writer=animation.FFMpegWriter(fps=15, bitrate=1800))
    else:
        ani.save(output, writer='pillow', fps=12)
    print(f"Saved → {output}")

# ── main ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input',     default='phonon_steps.txt')
    parser.add_argument('--nx',        type=int,   default=40)
    parser.add_argument('--ny',        type=int,   default=40)
    parser.add_argument('--tmin',      type=float, default=0.0)
    parser.add_argument('--tmax',      type=float, default=None)
    parser.add_argument('--nframes',   type=int,   default=80)
    parser.add_argument('--output',    default='phonon_density.gif')
    parser.add_argument('--synthetic', action='store_true')
    args = parser.parse_args()

    if args.synthetic or not Path(args.input).exists():
        print("Using synthetic data ...")
        df = make_synthetic()
    else:
        print(f"Loading {args.input} ...")
        df = load_steps(args.input)

    print(f"  {len(df)} steps | t=[{df['preT'].min():.2f}, {df['postT'].max():.2f}] ns")
    tmax = args.tmax or float(df['postT'].max())
    animate(df, args.nx, args.ny, args.tmin, tmax, args.nframes, args.output)

if __name__ == '__main__':
    main()

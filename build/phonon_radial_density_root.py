"""
phonon_radial_density.py
────────────────────────
Reads phonon_steps.root and produces THREE plots saved automatically:
  1. Density vs radial distance r (mm) and time (ns)
  2. Density vs azimuthal angle θ (degrees) and time (ns)
  3. Mean kinetic energy vs radial distance r (mm) and time (ns)

All particle types stored in the ROOT file are detected automatically
and plotted. No hard-coded particle filter — if it's in the file, it shows.

Units: ROOT file stores mm and ns. All axes use mm and ns directly.
"""

import argparse
from pathlib import Path
import warnings
warnings.filterwarnings('ignore', category=UserWarning, module='matplotlib')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import uproot

# ── Known ptype codes (written by SteppingAction) ────────────────────────────
PTYPE_NAME = {
    1: "phononL",
    2: "phononTF",
    3: "phononTS",
    4: "DriftElectron",
    5: "DriftHole",
    6: "e-",
    7: "e+",
    8: "gamma",
    9: "proton",
    10: "neutron",
    11: "pi+",
    12: "pi-",
    13: "pi0",
    14: "kaon+",
    15: "kaon-",
    16: "mu+",
    17: "mu-",
}

# Colour maps per particle — extended palette
PTYPE_CMAP = {
    "phononL":       "Blues",
    "phononTF":      "Greens",
    "phononTS":      "Reds",
    "DriftElectron": "Purples",
    "DriftHole":     "Oranges",
    "e-":            "YlOrRd",
    "e+":            "RdPu",
    "gamma":         "YlGn",
    "proton":        "cool",
    "neutron":       "copper",
    "pi+":           "PuBu",
    "pi-":           "BuPu",
    "pi0":           "GnBu",
    "kaon+":         "PuRd",
    "kaon-":         "RdGy",
    "mu+":           "bone",
    "mu-":           "pink",
}

# Ballistic speeds in Si (mm/ns) — drawn only on phonon radial panels
PHONON_SPEEDS = {
    "$v_L$=8.97":    8.97,
    "$v_{TF}$=5.77": 5.77,
    "$v_{TS}$=3.35": 3.35,
}
PHONON_PTYPES = {1, 2, 3}


# ── Loader ────────────────────────────────────────────────────────────────────

def load_steps(filepath):
    print(f"Opening: {filepath}")
    with uproot.open(filepath) as f:
        tree = f["Steps"]
        available = [b.name for b in tree.branches]
        print(f"  Branches available: {available}")

        wanted = ["ptype",
                  "preX_mm",  "preY_mm",  "preZ_mm",
                  "postX_mm", "postY_mm", "postZ_mm",
                  "preT_ns",  "postT_ns",
                  "preE_eV",  "preKE_eV",
                  "postE_eV", "postKE_eV",
                  "edep_eV",  "procID"]
        branches = [b for b in wanted if b in available]
        data = tree.arrays(branches, library="pd")

    df = data.rename(columns={
        "preX_mm":  "preX",  "preY_mm":  "preY",  "preZ_mm":  "preZ",
        "postX_mm": "postX", "postY_mm": "postY", "postZ_mm": "postZ",
        "preT_ns":  "preT",  "postT_ns": "postT",
        "preE_eV":  "preE",  "preKE_eV": "preKE",
        "postE_eV": "postE", "postKE_eV":"postKE",
    })

    print("\n  All particle types in file:")
    found = df["ptype"].value_counts().sort_index().to_dict()
    for code, count in found.items():
        name = PTYPE_NAME.get(int(code), f"unknown(ptype={int(code)})")
        print(f"    ptype={int(code):3d}  {name:<22s}  {count:>8d} steps")

    print(f"\n  Time  [ns]: [{df['preT'].min():.3f}, {df['postT'].max():.3f}]")
    print(f"  X     [mm]: [{df['preX'].min():.3f}, {df['preX'].max():.3f}]")
    print(f"  Y     [mm]: [{df['preY'].min():.3f}, {df['preY'].max():.3f}]")

    if df.empty:
        raise RuntimeError("No steps found in ROOT file.")
    return df


# ── Synthetic data ────────────────────────────────────────────────────────────

def make_synthetic(seed=42):
    rng = np.random.default_rng(seed)
    rows = []
    # Phonons
    for tid in range(1, 61):
        ptype = rng.choice([1, 2, 3])
        x, y = rng.uniform(-0.5, 0.5, 2)
        z = rng.uniform(-0.25, 0.25)
        v = rng.normal(size=3); v /= np.linalg.norm(v); t = 0.0
        for _ in range(300):
            dl = rng.exponential(0.5); dt = dl/5.0
            nx,ny_,nz = x+v[0]*dl, y+v[1]*dl, z+v[2]*dl
            for i,(lo,hi) in enumerate([(-4,4),(-4,4),(-0.2625,0.2625)]):
                c=[nx,ny_,nz][i]
                if c>hi: c=2*hi-c; v[i]=-v[i]
                if c<lo: c=2*lo-c; v[i]=-v[i]
                if i==0: nx=c
                elif i==1: ny_=c
                else: nz=c
            ke = rng.uniform(0.005, 0.056)
            rows.append({"ptype":ptype,"preX":x,"preY":y,"preT":t,
                         "postX":nx,"postY":ny_,"postT":t+dt,
                         "preKE":ke,"postKE":ke,"preE":ke,"postE":ke})
            x,y,z,t=nx,ny_,nz,t+dt
            
    # Add some EM secondaries
    for tid in range(61, 71):
        ptype = rng.choice([6, 8])
        x, y = rng.uniform(-2, 2, 2); t = rng.uniform(0.1, 5)  # Start > 0 for log compatibility
        for _ in range(20):
            dl = rng.exponential(2.0); dt = dl/300.0
            nx,ny_ = x+rng.normal()*dl, y+rng.normal()*dl
            nx = np.clip(nx,-4,4); ny_ = np.clip(ny_,-4,4)
            ke = rng.uniform(1e3, 1e6)
            rows.append({"ptype":ptype,"preX":x,"preY":y,"preT":t,
                         "postX":nx,"postY":ny_,"postT":t+dt,
                         "preKE":ke,"postKE":ke*0.99,"preE":ke,"postE":ke*0.99})
            x,y,t=nx,ny_,t+dt
    return pd.DataFrame(rows)


# ── Density computers ─────────────────────────────────────────────────────────

def radial_density(df, cx, cy, nr, rmin, rmax, nt, tmin, tmax, log_y=False):
    xm = 0.5*(df["preX"].values + df["postX"].values)
    ym = 0.5*(df["preY"].values + df["postY"].values)
    tm = 0.5*(df["preT"].values + df["postT"].values)
    r  = np.sqrt((xm-cx)**2 + (ym-cy)**2)
    
    r_edges = np.linspace(0, rmax, nr+1)
    if log_y:
        t_edges = np.logspace(np.log10(tmin), np.log10(tmax), nt+1)
    else:
        t_edges = np.linspace(tmin, tmax, nt+1)
        
    h, _, _ = np.histogram2d(r, tm, bins=[r_edges, t_edges],
                             range=[[0,rmax],[tmin,tmax]])
    
    # Differential time step widths for scaling when bins are non-uniform (log)
    dt = np.diff(t_edges)
    annulus = np.pi*(r_edges[1:]**2 - r_edges[:-1]**2)
    return h/(annulus[:,np.newaxis]*dt[np.newaxis,:]), r_edges, t_edges


def theta_density(df, cx, cy, ntheta, nt, tmin, tmax, log_y=False):
    xm = 0.5*(df["preX"].values + df["postX"].values)
    ym = 0.5*(df["preY"].values + df["postY"].values)
    tm = 0.5*(df["preT"].values + df["postT"].values)
    theta = np.degrees(np.arctan2(ym-cy, xm-cx))
    
    th_edges = np.linspace(-180, 180, ntheta+1)
    if log_y:
        t_edges = np.logspace(np.log10(tmin), np.log10(tmax), nt+1)
    else:
        t_edges = np.linspace(tmin, tmax, nt+1)
        
    dt  = np.diff(t_edges)
    dth = th_edges[1]-th_edges[0]
    h, _, _ = np.histogram2d(theta, tm, bins=[th_edges, t_edges],
                             range=[[-180,180],[tmin,tmax]])
    return h/(dth*dt[np.newaxis,:]), th_edges, t_edges


def energy_radial_density(df, cx, cy, nr, rmin, rmax, nt, tmin, tmax, log_y=False):
    if "preKE" not in df.columns or "postKE" not in df.columns:
        return None, None, None

    xm  = 0.5*(df["preX"].values  + df["postX"].values)
    ym  = 0.5*(df["preY"].values  + df["postY"].values)
    tm  = 0.5*(df["preT"].values  + df["postT"].values)
    ke  = 0.5*(df["preKE"].values + df["postKE"].values)
    r   = np.sqrt((xm-cx)**2 + (ym-cy)**2)

    r_edges = np.linspace(rmin, rmax, nr+1)
    if log_y:
        t_edges = np.logspace(np.log10(tmin), np.log10(tmax), nt+1)
    else:
        t_edges = np.linspace(tmin, tmax, nt+1)

    mask = (r >= rmin) & (r <= rmax)
    ke_sum, _, _ = np.histogram2d(r[mask], tm[mask], bins=[r_edges, t_edges],
                                  weights=ke[mask],
                                  range=[[rmin, rmax], [tmin, tmax]])
    counts, _, _ = np.histogram2d(r[mask], tm[mask], bins=[r_edges, t_edges],
                                  range=[[rmin, rmax], [tmin, tmax]])
    with np.errstate(invalid="ignore", divide="ignore"):
        mean_ke = np.where(counts > 0, ke_sum / counts, np.nan)
    return mean_ke, r_edges, t_edges


# ── Panel builder ─────────────────────────────────────────────────────────────

def make_panel(fig, ax, density, x_edges, t_edges,
               title, cmap_name, xlabel, ylabel, log_scale, log_y=False,
               draw_speeds=False, tmin=0, vmax_override=None):
    
    cmap = plt.get_cmap(cmap_name).copy()
    bg_color = "#fdfdfd"  
    cmap.set_bad(color=bg_color)
    cmap.set_under(color=bg_color)

    plot_data = np.where((density == 0) | np.isnan(density), np.nan, density)
    valid_data = density[~np.isnan(density) & (density > 0)]
    
    if len(valid_data) > 0:
        vmin = float(valid_data.min())
        vmax = float(vmax_override) if vmax_override is not None else float(np.percentile(valid_data, 99.5))
        if vmin >= vmax:
            vmax = float(valid_data.max()) if valid_data.max() > vmin else vmin + 1e-9
    else:
        vmin, vmax = 1e-9, 1.0

    if log_scale:
        norm = mcolors.LogNorm(vmin=max(vmin, 1e-30), vmax=vmax)
    else:
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    pcm = ax.pcolormesh(x_edges, t_edges, plot_data.T,
                        cmap=cmap, norm=norm, shading="flat")

    fig.colorbar(pcm, ax=ax, label=ylabel, pad=0.02, extend='both' if log_scale else 'max')
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel("Time  [ns]", fontsize=9)
    ax.set_title(title, fontsize=9)
    ax.set_xlim(x_edges[0], x_edges[-1])
    ax.set_ylim(t_edges[0], t_edges[-1])
    
    if log_y:
        ax.set_yscale('log')
    
    if draw_speeds:
        cols = ["#1f77b4", "#2ca02c", "#d62728"]
        ta = np.logspace(np.log10(t_edges[0]), np.log10(t_edges[-1]), 100) if log_y else np.array([tmin, t_edges[-1]])
        for (lbl, v), c in zip(PHONON_SPEEDS.items(), cols):
            ax.plot(v*ta, ta, "--", color=c, lw=1.2,
                    alpha=0.9, label=f"{lbl} mm/ns")
        ax.legend(fontsize=7, loc="upper left", framealpha=0.7)


# ── Plot builder ──────────────────────────────────────────────────────────────

def make_plots(df, cx, cy, nr, rmin, rmax, ntheta, nt,
               tmin, tmax, log_scale, log_y, output_prefix, show=False,
               cmax_radial=None, cmax_angular=None, cmax_energy=None):

    present_ptypes = sorted(df["ptype"].unique().astype(int))
    modes = {}
    for p in present_ptypes:
        name = PTYPE_NAME.get(p, f"unknown({p})")
        modes[name] = (p, df[df["ptype"]==p])

    n_panels = 1 + len(modes)
    panel_w  = min(5.0, 18.0/n_panels)

    # ── PLOT 1: radial ────────────────────────────────────────────────────────
    fig1, axes1 = plt.subplots(1, n_panels,
                                figsize=(panel_w*n_panels, 5.0),
                                constrained_layout=True)
    if n_panels==1: axes1=[axes1]
    fig1.suptitle(
        f"Radial density  |  source ({cx:.1f},{cy:.1f}) mm  "
        f"|  colour: steps mm⁻² ns⁻¹",
        fontsize=9, fontweight="bold")

    dens, r_edges, t_edges = radial_density(df, cx, cy, nr, rmin, rmax, nt, tmin, tmax, log_y=log_y)
    make_panel(fig1, axes1[0], dens, r_edges, t_edges,
               "All particles", "inferno", "r  [mm]",
               "Steps mm⁻² ns⁻¹", log_scale, log_y=log_y,
               draw_speeds=True, tmin=tmin,
               vmax_override=cmax_radial)

    for k, (name, (ptype, sub)) in enumerate(modes.items()):
        cmap = PTYPE_CMAP.get(name, "viridis")
        d, _, _ = radial_density(sub, cx, cy, nr, rmin, rmax, nt, tmin, tmax, log_y=log_y)
        draw_sp = ptype in PHONON_PTYPES
        make_panel(fig1, axes1[k+1], d, r_edges, t_edges,
                   name, cmap, "r  [mm]",
                   "Steps mm⁻² ns⁻¹", log_scale, log_y=log_y,
                   draw_speeds=draw_sp, tmin=tmin)

    out1 = f"{output_prefix}_radial.png"
    fig1.savefig(out1, dpi=150, bbox_inches="tight")
    print(f"Saved → {out1}")
    if show:
        plt.show()
    plt.close(fig1)

    # ── PLOT 2: angular ───────────────────────────────────────────────────────
    fig2, axes2 = plt.subplots(1, n_panels,
                                figsize=(panel_w*n_panels, 5.0),
                                constrained_layout=True)
    if n_panels==1: axes2=[axes2]
    fig2.suptitle(
        f"Angular density  |  source ({cx:.1f},{cy:.1f}) mm  "
        f"|  colour: steps deg⁻¹ ns⁻¹",
        fontsize=9, fontweight="bold")

    dth, th_edges, t_edges2 = theta_density(df, cx, cy, ntheta, nt, tmin, tmax, log_y=log_y)
    make_panel(fig2, axes2[0], dth, th_edges, t_edges2,
               "All particles", "inferno", "θ  [deg]",
               "Steps deg⁻¹ ns⁻¹", log_scale, log_y=log_y,
               vmax_override=cmax_angular)
    axes2[0].axvline(0, color='grey', lw=0.6, ls=':', alpha=0.5)

    for k, (name, (ptype, sub)) in enumerate(modes.items()):
        cmap = PTYPE_CMAP.get(name, "viridis")
        d, _, _ = theta_density(sub, cx, cy, ntheta, nt, tmin, tmax, log_y=log_y)
        make_panel(fig2, axes2[k+1], d, th_edges, t_edges2,
                   name, cmap, "θ  [deg]",
                   "Steps deg⁻¹ ns⁻¹", log_scale, log_y=log_y)
        axes2[k+1].axvline(0, color='grey', lw=0.6, ls=':', alpha=0.5)

    out2 = f"{output_prefix}_angular.png"
    fig2.savefig(out2, dpi=150, bbox_inches="tight")
    print(f"Saved → {out2}")
    if show:
        plt.show()
    plt.close(fig2)

    # ── PLOT 3: mean kinetic energy vs r and t ────────────────────────────────
    ke_all, r_edges_e, t_edges_e = energy_radial_density(
        df, cx, cy, nr, rmin, rmax, nt, tmin, tmax, log_y=log_y)

    if ke_all is not None:
        fig3, axes3 = plt.subplots(1, n_panels,
                                    figsize=(panel_w*n_panels, 5.0),
                                    constrained_layout=True)
        if n_panels==1: axes3=[axes3]
        fig3.suptitle(
            f"Mean kinetic energy  |  source ({cx:.1f},{cy:.1f}) mm  "
            f"|  colour: mean KE [eV]",
            fontsize=9, fontweight="bold")

        make_panel(fig3, axes3[0], ke_all, r_edges_e, t_edges_e,
                   "All particles", "plasma", "r  [mm]",
                   "Mean KE  [eV]", log_scale, log_y=log_y,
                   vmax_override=cmax_energy)

        for k, (name, (ptype, sub)) in enumerate(modes.items()):
            ke_m, _, _ = energy_radial_density(
                sub, cx, cy, nr, rmin, rmax, nt, tmin, tmax, log_y=log_y)
            if ke_m is None:
                continue
            cmap = PTYPE_CMAP.get(name, "viridis")
            make_panel(fig3, axes3[k+1], ke_m, r_edges_e, t_edges_e,
                       name, cmap, "r  [mm]", "Mean KE  [eV]", log_scale, log_y=log_y)

        out3 = f"{output_prefix}_energy.png"
        fig3.savefig(out3, dpi=150, bbox_inches="tight")
        print(f"Saved → {out3}")
        if show:
            plt.show()
        plt.close(fig3)
    else:
        print("  Energy branches not available — skipping energy plot.")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",         default="phonon_steps.root")
    parser.add_argument("--cx",            type=float, default=0.0)
    parser.add_argument("--cy",            type=float, default=0.0)
    parser.add_argument("--rmin",          type=float, default=0.0,
                        help="Min radial distance [mm] (default: 0)")
    parser.add_argument("--rmax",          type=float, default=None,
                        help="Max radial distance [mm] (default: slab diagonal)")
    parser.add_argument("--nr",            type=int,   default=80,
                        help="Radial bins")
    parser.add_argument("--ntheta",        type=int,   default=72,
                        help="Angular bins (default 72 = 5 deg each)")
    parser.add_argument("--tmin",          type=float, default=0.0,
                        help="Min time [ns] (default: 0)")
    parser.add_argument("--tmax",          type=float, default=None,
                        help="Max time [ns] (default: auto from data)")
    parser.add_argument("--nt",            type=int,   default=200,
                        help="Time bins")
    parser.add_argument("--log",           action="store_true",
                        help="Apply log scale to the colorbars")
    parser.add_argument("--log_y",         action="store_true",
                        help="Apply log scale to the time (Y) axis")
    parser.add_argument("--cmax-radial",  type=float, default=None,
                        help="Max colorbar for radial plot (steps mm^-2 ns^-1)")
    parser.add_argument("--cmax-angular", type=float, default=None,
                        help="Max colorbar for angular plot (steps deg^-1 ns^-1)")
    parser.add_argument("--cmax-energy",  type=float, default=None,
                        help="Max colorbar for energy plot (eV)")
    parser.add_argument("--show",          action="store_true",
                        help="Display plots on screen after saving")
    parser.add_argument("--output",       default="phonon_density")
    parser.add_argument("--synthetic",    action="store_true")
    args = parser.parse_args()

    if args.synthetic or not Path(args.input).exists():
        print("Using synthetic data ...")
        df = make_synthetic()
    else:
        df = load_steps(args.input)

    tmax = args.tmax or float(df["postT"].max())
    tmin = args.tmin
    
    # Safeguard for log scale on Y axis: a log axis cannot begin at exactly 0.0
    if args.log_y and tmin <= 0.0:
        # Check active non-zero time features, fallback safely if clean data is missing
        active_times = df.loc[df["postT"] > 0, "postT"]
        tmin = float(active_times.min()) * 0.5 if not active_times.empty else 0.1
        print(f"  [Warning] --log_y specified with tmin=0. Adjusted tmin to {tmin:.3f} ns to allow log calculations.")

    rmax = args.rmax or float(
        np.sqrt(max(abs(4.0-args.cx),abs(-4.0-args.cx))**2 +
                max(abs(4.0-args.cy),abs(-4.0-args.cy))**2)*1.05)

    print(f"  Axis limits:  r=[{args.rmin:.2f}, {rmax:.2f}] mm  "
          f"t=[{tmin:.2f}, {tmax:.2f}] ns")

    make_plots(df, args.cx, args.cy, args.nr, args.rmin, rmax,
               args.ntheta, args.nt, tmin, tmax,
               args.log, args.log_y, args.output, args.show,
               cmax_radial=args.cmax_radial,
               cmax_angular=args.cmax_angular,
               cmax_energy=args.cmax_energy)


if __name__=="__main__":
    main()

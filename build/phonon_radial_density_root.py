"""
phonon_radial_density.py
────────────────────────
Reads phonon_steps.root and produces TWO plots saved automatically:
  1. Density vs radial distance r (mm) and time (ns)
  2. Density vs azimuthal angle θ (degrees) and time (ns)

All particle types stored in the ROOT file are detected automatically
and plotted. No hard-coded particle filter — if it's in the file, it shows.

Units: ROOT file stores mm and ns. All axes use mm and ns directly.
"""

import argparse
from pathlib import Path
import matplotlib.colors as mcolors
import warnings
warnings.filterwarnings('ignore', category=UserWarning, module='matplotlib')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot

# ── Known ptype codes (written by SteppingAction) ────────────────────────────
# These are the codes YOUR SteppingAction assigns.
# Unknown codes are labelled "unknown(N)" and still plotted.
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
    """
    Reads all entries from phonon_steps.root.
    Units in file: mm and ns — used directly, no conversion.
    All ptype values present are reported.
    """
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

    # Report every ptype found — no filtering
    print("\n  All particle types in file:")
    found = df["ptype"].value_counts().sort_index().to_dict()
    for code, count in found.items():
        name = PTYPE_NAME.get(int(code), f"unknown(ptype={int(code)})")
        print(f"    ptype={int(code):3d}  {name:<22s}  {count:>8d} steps")

    # Unit sanity check
    print(f"\n  Time  [ns]: [{df['preT'].min():.3f}, {df['postT'].max():.3f}]")
    print(f"  X     [mm]: [{df['preX'].min():.3f}, {df['preX'].max():.3f}]")
    print(f"  Y     [mm]: [{df['preY'].min():.3f}, {df['preY'].max():.3f}]")

    if df.empty:
        raise RuntimeError("No steps found in ROOT file.")
    return df


# ── Synthetic data ────────────────────────────────────────────────────────────

def make_synthetic(seed=42):
    """Generate realistic synthetic data including phonons and a few EM particles."""
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
            ke = rng.uniform(0.005, 0.056)  # eV, phonon range in Si
            rows.append({"ptype":ptype,"preX":x,"preY":y,"preT":t,
                         "postX":nx,"postY":ny_,"postT":t+dt,
                         "preKE":ke,"postKE":ke,"preE":ke,"postE":ke})
            x,y,z,t=nx,ny_,nz,t+dt
    # Add some EM secondaries (e-, gamma) to demonstrate multi-particle plotting
    for tid in range(61, 71):
        ptype = rng.choice([6, 8])  # e- or gamma
        x, y = rng.uniform(-2, 2, 2); t = rng.uniform(0, 5)
        for _ in range(20):
            dl = rng.exponential(2.0); dt = dl/300.0
            nx,ny_ = x+rng.normal()*dl, y+rng.normal()*dl
            nx = np.clip(nx,-4,4); ny_ = np.clip(ny_,-4,4)
            ke = rng.uniform(1e3, 1e6)  # eV, EM particles
            rows.append({"ptype":ptype,"preX":x,"preY":y,"preT":t,
                         "postX":nx,"postY":ny_,"postT":t+dt,
                         "preKE":ke,"postKE":ke*0.99,"preE":ke,"postE":ke*0.99})
            x,y,t=nx,ny_,t+dt
    return pd.DataFrame(rows)


# ── Density computers ─────────────────────────────────────────────────────────

def radial_density(df, cx, cy, nr, rmin, rmax, nt, tmin, tmax):
    """Density in steps / mm² / ns."""
    xm = 0.5*(df["preX"].values + df["postX"].values)
    ym = 0.5*(df["preY"].values + df["postY"].values)
    tm = 0.5*(df["preT"].values + df["postT"].values)  # ns
    r  = np.sqrt((xm-cx)**2 + (ym-cy)**2)
    r_edges = np.linspace(0,    rmax, nr+1)
    t_edges = np.linspace(tmin, tmax, nt+1)
    dt = t_edges[1]-t_edges[0]
    h, _, _ = np.histogram2d(r, tm, bins=[r_edges, t_edges],
                              range=[[0,rmax],[tmin,tmax]])
    annulus = np.pi*(r_edges[1:]**2 - r_edges[:-1]**2)
    return h/(annulus[:,np.newaxis]*dt), r_edges, t_edges


def theta_density(df, cx, cy, ntheta, nt, tmin, tmax):
    """Density in steps / deg / ns."""
    xm = 0.5*(df["preX"].values + df["postX"].values)
    ym = 0.5*(df["preY"].values + df["postY"].values)
    tm = 0.5*(df["preT"].values + df["postT"].values)  # ns
    theta = np.degrees(np.arctan2(ym-cy, xm-cx))
    th_edges = np.linspace(-180, 180, ntheta+1)
    t_edges  = np.linspace(tmin, tmax, nt+1)
    dt  = t_edges[1]-t_edges[0]
    dth = th_edges[1]-th_edges[0]
    h, _, _ = np.histogram2d(theta, tm, bins=[th_edges, t_edges],
                              range=[[-180,180],[tmin,tmax]])
    return h/(dth*dt), th_edges, t_edges




def energy_radial_density(df, cx, cy, nr, rmin, rmax, nt, tmin, tmax):
    """
    Mean kinetic energy [eV] per (r, t) bin.
    Uses the midpoint KE = 0.5*(preKE + postKE) as the step energy.
    Returns mean_KE[nr, nt] and the bin edges.
    """
    if "preKE" not in df.columns or "postKE" not in df.columns:
        return None, None, None

    xm  = 0.5*(df["preX"].values  + df["postX"].values)
    ym  = 0.5*(df["preY"].values  + df["postY"].values)
    tm  = 0.5*(df["preT"].values  + df["postT"].values)
    ke  = 0.5*(df["preKE"].values + df["postKE"].values)
    r   = np.sqrt((xm-cx)**2 + (ym-cy)**2)

    r_edges = np.linspace(rmin, rmax, nr+1)
    t_edges = np.linspace(tmin, tmax, nt+1)

    mask = (r >= rmin) & (r <= rmax)
    # Sum of KE and count per bin (for mean)
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
               title, cmap, xlabel, ylabel, log_scale,
               draw_speeds=False, tmin=0, contour=False, nlevels=8,
               vmax_override=None):
    pos = density[density>0]
    vmin = float(pos.min()) if len(pos) else 1e-9
    vmax = float(vmax_override) if vmax_override is not None \
           else (float(density.max()) if density.max()>0 else 1.0)
    norm = mcolors.LogNorm(vmin=vmin,vmax=vmax) if log_scale \
           else mcolors.Normalize(vmin=0,vmax=vmax)

    # Bin centres for contour (pcolormesh uses edges)
    xc = 0.5*(x_edges[:-1]+x_edges[1:])
    tc = 0.5*(t_edges[:-1]+t_edges[1:])

    if contour and vmax > vmin and not (log_scale and vmin <= 0):
        # Filled contours + contour lines
        data_c = np.where(np.isnan(density), 0, density).T
        if log_scale:
            levels = np.unique(
                np.logspace(np.log10(vmin), np.log10(vmax), nlevels+1))
        else:
            levels = np.unique(
                np.linspace(vmin, vmax, nlevels+1))
        if len(levels) >= 2:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                pcm = ax.contourf(xc, tc, data_c, levels=levels,
                                  cmap=cmap, norm=norm)
            ax.contour(xc, tc, data_c, levels=levels,
                       colors="white", linewidths=0.4, alpha=0.5)
        else:
            # Degenerate data — fall back to mesh
            pcm = ax.pcolormesh(x_edges, t_edges, density.T,
                                cmap=cmap, norm=norm, shading="flat")
    else:
        pcm = ax.pcolormesh(x_edges, t_edges, density.T,
                            cmap=cmap, norm=norm, shading="flat")

    fig.colorbar(pcm, ax=ax, label=ylabel, pad=0.02)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel("Time  [ns]", fontsize=9)
    ax.set_title(title, fontsize=9)
    ax.set_xlim(x_edges[0], x_edges[-1])
    ax.set_ylim(t_edges[0], t_edges[-1])
    if draw_speeds:
        cols = ["#aec7e8","#98df8a","#ff9896"]
        ta = np.array([tmin, t_edges[-1]])
        for (lbl,v),c in zip(PHONON_SPEEDS.items(),cols):
            ax.plot(v*ta, ta, "--", color=c, lw=1.0,
                    alpha=0.9, label=f"{lbl} mm/ns")
        ax.legend(fontsize=6, loc="upper left",
                  framealpha=0.3, labelcolor="white")


# ── Plot builder ──────────────────────────────────────────────────────────────

def make_plots(df, cx, cy, nr, rmin, rmax, ntheta, nt,
               tmin, tmax, log_scale, output_prefix, show=False,
               contour=False, nlevels=8,
               cmax_radial=None, cmax_angular=None, cmax_energy=None):

    # Discover all ptypes present — completely automatic
    present_ptypes = sorted(df["ptype"].unique().astype(int))
    modes = {}
    for p in present_ptypes:
        name = PTYPE_NAME.get(p, f"unknown({p})")
        modes[name] = (p, df[df["ptype"]==p])

    n_panels = 1 + len(modes)
    panel_w  = min(5.0, 18.0/n_panels)   # shrink if many panels

    # ── PLOT 1: radial ────────────────────────────────────────────────────────
    fig1, axes1 = plt.subplots(1, n_panels,
                                figsize=(panel_w*n_panels, 5.0),
                                constrained_layout=True)
    if n_panels==1: axes1=[axes1]
    fig1.suptitle(
        f"Radial density  |  source ({cx:.1f},{cy:.1f}) mm  "
        f"|  colour: steps mm⁻² ns⁻¹",
        fontsize=9, fontweight="bold")

    dens,r_edges,t_edges = radial_density(df,cx,cy,nr,rmin,rmax,nt,tmin,tmax)
    make_panel(fig1, axes1[0], dens, r_edges, t_edges,
               "All particles", "inferno", "r  [mm]",
               "Steps mm⁻² ns⁻¹", log_scale,
               draw_speeds=True, tmin=tmin,
               contour=contour, nlevels=nlevels,
               vmax_override=cmax_radial)

    for k,(name,(ptype,sub)) in enumerate(modes.items()):
        cmap = PTYPE_CMAP.get(name,"viridis")
        d,_,_ = radial_density(sub,cx,cy,nr,rmin,rmax,nt,tmin,tmax)
        draw_sp = ptype in PHONON_PTYPES
        make_panel(fig1, axes1[k+1], d, r_edges, t_edges,
                   name, cmap, "r  [mm]",
                   "Steps mm⁻² ns⁻¹", log_scale,
                   draw_speeds=draw_sp, tmin=tmin,
                   contour=contour, nlevels=nlevels)

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

    dth,th_edges,t_edges2 = theta_density(df,cx,cy,ntheta,nt,tmin,tmax)
    make_panel(fig2, axes2[0], dth, th_edges, t_edges2,
               "All particles", "inferno", "θ  [deg]",
               "Steps deg⁻¹ ns⁻¹", log_scale,
               contour=contour, nlevels=nlevels,
               vmax_override=cmax_angular)
    axes2[0].axvline(0, color='white', lw=0.6, ls=':', alpha=0.5)

    for k,(name,(ptype,sub)) in enumerate(modes.items()):
        cmap = PTYPE_CMAP.get(name,"viridis")
        d,_,_ = theta_density(sub,cx,cy,ntheta,nt,tmin,tmax)
        make_panel(fig2, axes2[k+1], d, th_edges, t_edges2,
                   name, cmap, "θ  [deg]",
                   "Steps deg⁻¹ ns⁻¹", log_scale,
                   contour=contour, nlevels=nlevels)
        axes2[k+1].axvline(0, color='white', lw=0.6, ls=':', alpha=0.5)

    out2 = f"{output_prefix}_angular.png"
    fig2.savefig(out2, dpi=150, bbox_inches="tight")
    print(f"Saved → {out2}")
    if show:
        plt.show()
    plt.close(fig2)

    # ── PLOT 3: mean kinetic energy vs r and t ────────────────────────────────
    ke_all, r_edges_e, t_edges_e = energy_radial_density(
        df, cx, cy, nr, rmin, rmax, nt, tmin, tmax)

    if ke_all is not None:
        fig3, axes3 = plt.subplots(1, n_panels,
                                    figsize=(panel_w*n_panels, 5.0),
                                    constrained_layout=True)
        if n_panels==1: axes3=[axes3]
        fig3.suptitle(
            f"Mean kinetic energy  |  source ({cx:.1f},{cy:.1f}) mm  "
            f"|  colour: mean KE [eV]",
            fontsize=9, fontweight="bold")

        # All-particle panel — linear scale for energy (log optional)
        pos3 = ke_all[~np.isnan(ke_all)]
        vmin3 = float(pos3.min()) if len(pos3) else 1e-9
        vmax3 = float(cmax_energy) if cmax_energy is not None \
               else (float(pos3.max()) if len(pos3) else 1.0)
        norm3 = mcolors.LogNorm(vmin=max(vmin3,1e-30), vmax=vmax3) if log_scale \
               else mcolors.Normalize(vmin=vmin3, vmax=vmax3)

        xc_e = 0.5*(r_edges_e[:-1]+r_edges_e[1:])
        tc_e = 0.5*(t_edges_e[:-1]+t_edges_e[1:])
        data_e = np.where(np.isnan(ke_all), 0, ke_all).T
        if contour and vmax3 > vmin3:
            if log_scale:
                lv3 = np.unique(np.logspace(
                    np.log10(max(vmin3,1e-30)), np.log10(vmax3), nlevels+1))
            else:
                lv3 = np.unique(np.linspace(vmin3, vmax3, nlevels+1))
            if len(lv3) >= 2:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    pcm3 = axes3[0].contourf(xc_e, tc_e, data_e,
                                              levels=lv3, cmap="plasma", norm=norm3)
                axes3[0].contour(xc_e, tc_e, data_e, levels=lv3,
                                 colors="white", linewidths=0.4, alpha=0.5)
            else:
                pcm3 = axes3[0].pcolormesh(r_edges_e, t_edges_e, ke_all.T,
                                            cmap="plasma", norm=norm3, shading="flat")
        else:
            pcm3 = axes3[0].pcolormesh(r_edges_e, t_edges_e, ke_all.T,
                                        cmap="plasma", norm=norm3,
                                        shading="flat")
        fig3.colorbar(pcm3, ax=axes3[0], label="Mean KE  [eV]", pad=0.02)
        axes3[0].set(xlabel="r  [mm]", ylabel="Time  [ns]",
                     title="All particles")
        axes3[0].set_xlim(r_edges_e[0], r_edges_e[-1])
        axes3[0].set_ylim(t_edges_e[0], t_edges_e[-1])

        # Per-mode panels
        mode_cmaps_e = {name: PTYPE_CMAP.get(name,"viridis")
                        for name in modes}
        for k,(name,(ptype,sub)) in enumerate(modes.items()):
            ke_m, _, _ = energy_radial_density(
                sub, cx, cy, nr, rmin, rmax, nt, tmin, tmax)
            if ke_m is None:
                continue
            pos_m = ke_m[~np.isnan(ke_m)]
            if len(pos_m) == 0:
                continue
            vmin_m = float(pos_m.min()); vmax_m = float(pos_m.max())
            norm_m = mcolors.LogNorm(vmin=max(vmin_m,1e-30), vmax=vmax_m)                      if log_scale else mcolors.Normalize(vmin=vmin_m, vmax=vmax_m)
            data_m = np.where(np.isnan(ke_m), 0, ke_m).T
            if contour and vmax_m > vmin_m:
                if log_scale:
                    lv_m = np.unique(np.logspace(
                        np.log10(max(vmin_m,1e-30)), np.log10(vmax_m), nlevels+1))
                else:
                    lv_m = np.unique(np.linspace(vmin_m, vmax_m, nlevels+1))
                if len(lv_m) >= 2:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        pcm_m = axes3[k+1].contourf(xc_e, tc_e, data_m,
                                                     levels=lv_m,
                                                     cmap=mode_cmaps_e[name], norm=norm_m)
                    axes3[k+1].contour(xc_e, tc_e, data_m, levels=lv_m,
                                        colors="white", linewidths=0.4, alpha=0.5)
                else:
                    pcm_m = axes3[k+1].pcolormesh(r_edges_e, t_edges_e, ke_m.T,
                                                   cmap=mode_cmaps_e[name],
                                                   norm=norm_m, shading="flat")
            else:
                pcm_m = axes3[k+1].pcolormesh(r_edges_e, t_edges_e, ke_m.T,
                                               cmap=mode_cmaps_e[name],
                                               norm=norm_m, shading="flat")
            fig3.colorbar(pcm_m, ax=axes3[k+1], label="Mean KE  [eV]", pad=0.02)
            axes3[k+1].set(xlabel="r  [mm]", ylabel="Time  [ns]", title=name)
            axes3[k+1].set_xlim(r_edges_e[0], r_edges_e[-1])
            axes3[k+1].set_ylim(t_edges_e[0], t_edges_e[-1])

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
    parser.add_argument("--input",     default="phonon_steps.root")
    parser.add_argument("--cx",        type=float, default=0.0)
    parser.add_argument("--cy",        type=float, default=0.0)
    parser.add_argument("--rmin",      type=float, default=0.0,
                        help="Min radial distance [mm] (default: 0)")
    parser.add_argument("--rmax",      type=float, default=None,
                        help="Max radial distance [mm] (default: slab diagonal)")
    parser.add_argument("--nr",        type=int,   default=80,
                        help="Radial bins")
    parser.add_argument("--ntheta",    type=int,   default=72,
                        help="Angular bins (default 72 = 5 deg each)")
    parser.add_argument("--tmin",      type=float, default=0.0,
                        help="Min time [ns] (default: 0)")
    parser.add_argument("--tmax",      type=float, default=None,
                        help="Max time [ns] (default: auto from data)")
    parser.add_argument("--nt",        type=int,   default=200,
                        help="Time bins")
    parser.add_argument("--log",          action="store_true")
    parser.add_argument("--contour",      action="store_true",
                        help="Draw filled contours instead of flat colour mesh")
    parser.add_argument("--nlevels",      type=int, default=8,
                        help="Number of contour levels (default: 8)")
    parser.add_argument("--cmax-radial",  type=float, default=None,
                        help="Max colorbar for radial plot (steps mm^-2 ns^-1)")
    parser.add_argument("--cmax-angular", type=float, default=None,
                        help="Max colorbar for angular plot (steps deg^-1 ns^-1)")
    parser.add_argument("--cmax-energy",  type=float, default=None,
                        help="Max colorbar for energy plot (eV)")
    parser.add_argument("--show",         action="store_true",
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
    rmax = args.rmax or float(
        np.sqrt(max(abs(4.0-args.cx),abs(-4.0-args.cx))**2 +
                max(abs(4.0-args.cy),abs(-4.0-args.cy))**2)*1.05)

    print(f"  Axis limits:  r=[{args.rmin:.2f}, {rmax:.2f}] mm  "
          f"t=[{args.tmin:.2f}, {tmax:.2f}] ns")

    make_plots(df, args.cx, args.cy, args.nr, args.rmin, rmax,
               args.ntheta, args.nt, args.tmin, tmax,
               args.log, args.output, args.show,
               contour=args.contour, nlevels=args.nlevels,
               cmax_radial=args.cmax_radial,
               cmax_angular=args.cmax_angular,
               cmax_energy=args.cmax_energy)


if __name__=="__main__":
    main()

#!/usr/bin/env python3
"""
animate_root.py
───────────────
Reads phonon_steps.root from G4CMP and renders a time-lapse MP4/GIF
showing particle propagation through the Si + SiO2 bilayer detector.

Usage:
    python animate_root.py phonon_steps.root
    python animate_root.py phonon_steps.root --out phonons.mp4
    python animate_root.py phonon_steps.root --view top   --frames 80 --fps 20
    python animate_root.py phonon_steps.root --view side
    python animate_root.py phonon_steps.root --view both  --tmax 100

--view options:
    top   XY plane — face-on phonon wavefront (like Geant4 vis at 0 0)
    side  XZ plane — layer structure + Z propagation
    both  side-by-side XY + XZ panels (default)

Requirements:
    pip install uproot numpy matplotlib imageio
    pip install imageio[ffmpeg]   # only for .mp4 output
"""

import argparse
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import LineCollection
import uproot

# ── Particle type definitions ─────────────────────────────────────────────────
# (name, hex-color, alpha, linewidth)
PTYPE_STYLE = {
    1: ("phononL",       "#4499FF", 0.75, 0.7),
    2: ("phononTF",      "#33DD66", 0.75, 0.7),
    3: ("phononTS",      "#FF4444", 0.75, 0.7),
    4: ("DriftElectron", "#BB44FF", 0.90, 1.1),
    5: ("DriftHole",     "#FF8800", 0.90, 1.1),
    6: ("e⁻",           "#FFFF55", 0.60, 0.6),
    7: ("e⁺",           "#FF88FF", 0.60, 0.6),
    8: ("γ",             "#AAFFAA", 0.50, 0.5),
    9: ("proton",        "#FFFFFF", 1.00, 1.8),
}

# ── Detector geometry (mm) ────────────────────────────────────────────────────
SI_Z_MIN  = -0.2625
SI_Z_MAX  =  0.2625
SIO2_BOT  =  0.2625
SIO2_TOP  =  0.2628          # 280 nm thick
SIO2_DRAW =  0.010           # drawn thicker for visibility
WORLD_XY  =  5.0

# ── Colours ───────────────────────────────────────────────────────────────────
BG       = "#07070f"
SI_COL   = "#0b1828"
SIO2_COL = "#1a1200"
GRID_COL = "#141e2c"


# ─────────────────────────────────────────────────────────────────────────────
def load_event(path: str, event_id: int) -> dict:
    with uproot.open(path) as f:
        arr = f["Steps"].arrays(library="np")
    m = arr["event"] == event_id
    n = int(m.sum())
    if n == 0:
        raise ValueError(f"No steps found for event {event_id}")
    return {k: v[m] for k, v in arr.items()}, n


def compute_xy_lim(arr: dict, pad: float = 1.25) -> float:
    """Return half-width for XY view based on 98th-percentile phonon spread."""
    ph = arr["ptype"] <= 3
    if ph.sum() == 0:
        return WORLD_XY
    xs = np.concatenate([arr["preX_mm"][ph], arr["postX_mm"][ph]])
    ys = np.concatenate([arr["preY_mm"][ph], arr["postY_mm"][ph]])
    lim = max(np.percentile(np.abs(xs), 98),
              np.percentile(np.abs(ys), 98)) * pad
    return max(lim, 0.05)


def build_frame_times(arr: dict, n_frames: int, tmax: float) -> np.ndarray:
    """Log-linear mix of frame timestamps for good early+late coverage."""
    t0   = max(float(arr["preT_ns"].min()), 1e-5)
    tmax = tmax or float(arr["postT_ns"].max())
    t_log = np.logspace(np.log10(t0), np.log10(tmax), int(n_frames * 0.65))
    t_lin = np.linspace(t0, tmax,               int(n_frames * 0.35))
    times = np.unique(np.concatenate([t_log, t_lin]))
    return times[:n_frames]


# ─────────────────────────────────────────────────────────────────────────────
def _style_ax(ax):
    ax.tick_params(colors="white", labelsize=7)
    for sp in ax.spines.values():
        sp.set_edgecolor("#1a2535")
    ax.grid(True, color=GRID_COL, lw=0.3, alpha=0.6)


def _setup_xy(ax, half: float):
    ax.set_facecolor(BG)
    # Si slab face (square)
    ax.add_patch(mpatches.Rectangle(
        (-WORLD_XY, -WORLD_XY), 2*WORLD_XY, 2*WORLD_XY,
        color=SI_COL, zorder=0))
    ax.set_xlim(-half, half)
    ax.set_ylim(-half, half)
    ax.set_xlabel("X (mm)", color="white", fontsize=8)
    ax.set_ylabel("Y (mm)", color="white", fontsize=8)
    ax.set_title("Top view (XY) — phonon wavefront", color="#99aacc", fontsize=8)
    _style_ax(ax)


def _setup_xz(ax):
    ax.set_facecolor(BG)
    # Si slab
    ax.add_patch(mpatches.Rectangle(
        (-WORLD_XY, SI_Z_MIN), 2*WORLD_XY, SI_Z_MAX - SI_Z_MIN,
        color=SI_COL, zorder=0))
    # SiO2 layer (drawn larger for visibility)
    ax.add_patch(mpatches.Rectangle(
        (-WORLD_XY, SIO2_BOT), 2*WORLD_XY, SIO2_DRAW,
        color=SIO2_COL, zorder=0))
    ax.axhline(SIO2_BOT, color="#2a3f55", lw=0.7, ls="--", alpha=0.8, zorder=1)
    ax.text(WORLD_XY * 0.96, SIO2_BOT + SIO2_DRAW * 0.5,
            "SiO₂ (280 nm)", color="#cc9922", fontsize=6,
            ha="right", va="center")
    ax.text(WORLD_XY * 0.96, 0,
            "Si (525 µm)", color="#4488bb", fontsize=7,
            ha="right", va="center", alpha=0.8)
    ax.set_xlim(-WORLD_XY, WORLD_XY)
    ax.set_ylim(SI_Z_MIN * 1.06, SIO2_BOT + SIO2_DRAW * 2)
    ax.set_xlabel("X (mm)", color="white", fontsize=8)
    ax.set_ylabel("Z (mm)", color="white", fontsize=8)
    ax.set_title("Side view (XZ) — layer structure", color="#99aacc", fontsize=8)
    _style_ax(ax)


def _draw_tracks(arr: dict, t_cut: float, ax_xy, ax_xz):
    """Draw all steps with postT_ns <= t_cut. Returns legend info."""
    mask_t = arr["postT_ns"] <= t_cut
    legend = {}

    for pt, (name, color, alpha, lw) in PTYPE_STYLE.items():
        m = mask_t & (arr["ptype"] == pt)
        if m.sum() == 0:
            continue

        x0, x1 = arr["preX_mm"][m],  arr["postX_mm"][m]
        y0, y1 = arr["preY_mm"][m],  arr["postY_mm"][m]
        z0, z1 = arr["preZ_mm"][m],  arr["postZ_mm"][m]

        def add_lc(ax, A0, B0, A1, B1):
            segs = np.stack(
                [np.stack([A0, B0], axis=1),
                 np.stack([A1, B1], axis=1)], axis=1)
            ax.add_collection(LineCollection(
                segs, colors=color, alpha=alpha,
                linewidths=lw, zorder=2))

        if ax_xy is not None:
            add_lc(ax_xy, x0, y0, x1, y1)
        if ax_xz is not None:
            add_lc(ax_xz, x0, z0, x1, z1)

        legend[name] = (color, alpha, int(m.sum()))

    return legend


def _add_legend(ax, legend: dict):
    handles = [
        mpatches.Patch(color=c, alpha=min(a + 0.15, 1.0),
                       label=f"{nm}  {n:,}")
        for nm, (c, a, n) in legend.items()
    ]
    if handles:
        ax.legend(handles=handles, loc="upper left", fontsize=6,
                  facecolor="#0f1020", edgecolor="#2a3550",
                  labelcolor="white", framealpha=0.88, ncol=1,
                  handlelength=1.2, handleheight=0.8)


def render_frame(arr: dict, t_cut: float, event_id: int,
                 view: str, xy_half: float, figsize: tuple) -> np.ndarray:
    if view == "both":
        fig, (ax_xy, ax_xz) = plt.subplots(
            1, 2, figsize=figsize, facecolor=BG,
            gridspec_kw={"wspace": 0.30})
        _setup_xy(ax_xy, xy_half)
        _setup_xz(ax_xz)
        legend = _draw_tracks(arr, t_cut, ax_xy, ax_xz)
        _add_legend(ax_xz, legend)

    elif view == "top":
        fig, ax = plt.subplots(figsize=figsize, facecolor=BG)
        _setup_xy(ax, xy_half)
        legend = _draw_tracks(arr, t_cut, ax, None)
        _add_legend(ax, legend)

    else:  # side
        fig, ax = plt.subplots(figsize=figsize, facecolor=BG)
        _setup_xz(ax)
        legend = _draw_tracks(arr, t_cut, None, ax)
        _add_legend(ax, legend)

    fig.text(0.01, 0.015,
             f"t ≤ {t_cut:.3f} ns   |   event {event_id}",
             color="#6677aa", fontsize=7, va="bottom")

    fig.tight_layout(pad=0.5)
    fig.canvas.draw()

    try:
        buf = np.asarray(fig.canvas.buffer_rgba())
        rgb = buf[:, :, :3].copy()
    except AttributeError:
        w, h = fig.canvas.get_width_height()
        buf = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
        rgb = buf.reshape(h, w, 3)

    plt.close(fig)
    return rgb


# ─────────────────────────────────────────────────────────────────────────────
def find_windows_desktop() -> Path | None:
    try:
        for user_dir in sorted(Path("/mnt/c/Users").iterdir()):
            if user_dir.name in ("Public", "Default", "All Users"):
                continue
            d = user_dir / "Desktop"
            if d.exists() and os.access(d, os.W_OK):
                return d
    except (FileNotFoundError, PermissionError):
        pass
    return None


# ─────────────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(
        description="Animate G4CMP phonon_steps.root")
    ap.add_argument("root_file", nargs="?", default="phonon_steps.root",
                    help="Path to ROOT file (default: phonon_steps.root)")
    ap.add_argument("--event",  type=int,   default=0,
                    help="Event index to render (default: 0)")
    ap.add_argument("--frames", type=int,   default=80,
                    help="Number of animation frames (default: 80)")
    ap.add_argument("--fps",    type=float, default=15,
                    help="Frames per second (default: 15)")
    ap.add_argument("--tmax",   type=float, default=None,
                    help="Max time in ns (default: auto from data)")
    ap.add_argument("--view",   default="both",
                    choices=["top", "side", "both"],
                    help="top=XY  side=XZ  both=side-by-side (default: both)")
    ap.add_argument("--out",    default="phonons.gif",
                    help="Output file (.gif or .mp4, default: phonons.gif)")
    ap.add_argument("--no-desktop", action="store_true",
                    help="Skip copying to Windows Desktop")
    args = ap.parse_args()

    # ── Load ──────────────────────────────────────────────────────────────────
    path = Path(args.root_file)
    if not path.exists():
        print(f"ERROR: {path} not found.", file=sys.stderr)
        sys.exit(1)

    print(f"Loading  {path}")
    arr, n_steps = load_event(str(path), args.event)
    print(f"  Event {args.event}: {n_steps:,} steps")

    for pt, (name, *_) in PTYPE_STYLE.items():
        cnt = int((arr["ptype"] == pt).sum())
        if cnt:
            t = arr["postT_ns"][arr["ptype"] == pt]
            print(f"    {name:16s}  {cnt:6,} steps   "
                  f"t = [{t.min():.3f}, {t.max():.1f}] ns")

    xy_half = compute_xy_lim(arr)
    print(f"  XY display window: ±{xy_half:.3f} mm")

    times = build_frame_times(arr, args.frames, args.tmax)
    print(f"  Time axis: {times[0]:.4f} – {times[-1]:.2f} ns  "
          f"({len(times)} frames)")

    figsize = (13, 5.5) if args.view == "both" else (7, 6)

    # ── Render frames ─────────────────────────────────────────────────────────
    print(f"  Rendering {len(times)} frames ({args.view} view) ...")
    frames = []
    for i, t in enumerate(times):
        if i % 10 == 0 or i == len(times) - 1:
            print(f"    [{i+1:3d}/{len(times)}]  t = {t:.4f} ns")
        frames.append(render_frame(arr, t, args.event,
                                   args.view, xy_half, figsize))

    # Hold last frame for 2 s
    frames += [frames[-1]] * int(args.fps * 2)

    # ── Write output ──────────────────────────────────────────────────────────
    out = Path(args.out)
    print(f"\nWriting {out}  ({len(frames)} frames @ {args.fps} fps) ...")

    try:
        import imageio.v2 as iio
    except ImportError:
        import imageio as iio

    if out.suffix.lower() == ".gif":
        iio.mimsave(str(out), frames, fps=args.fps, loop=0)
    else:
        iio.mimsave(str(out), frames, fps=args.fps,
                    macro_block_size=None,
                    ffmpeg_params=["-crf", "18", "-pix_fmt", "yuv420p"])

    sz = out.stat().st_size / 1e6
    print(f"Done →  {out}  ({sz:.1f} MB)")

    # ── Copy to Windows Desktop ───────────────────────────────────────────────
    if not args.no_desktop:
        desktop = find_windows_desktop()
        if desktop:
            dst = desktop / out.name
            shutil.copy2(out, dst)
            print(f"Copied to Desktop →  {dst}")
        else:
            print("Windows Desktop not found. Copy manually:")
            print(f"  cp {out.resolve()}  /mnt/c/Users/<YourName>/Desktop/")


if __name__ == "__main__":
    main()

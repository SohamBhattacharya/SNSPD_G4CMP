"""
phonon_density_animation_root.py
────────────────────────────────
Animated GIF of phonon XY density on a discretised mesh vs time using PyROOT.

Usage:
    python phonon_density_animation_root.py --input phonon_steps.txt
    python phonon_density_animation_root.py --input phonon_steps.txt \
        --nx 40 --ny 40 --tmax 200 --nframes 80 --output anim.gif
"""

import argparse
import sys
import os
from pathlib import Path

# Initialize ROOT in batch mode to prevent graphical windows flashing during rendering
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)  # Hide standard statistics box for clean visualization

# ── shared configurations ─────────────────────────────────────────────────────
PTYPE_MAP  = {'phononL': 1, 'phononTF': 2, 'phononTS': 3}
PTYPE_NAME = {1: 'phononL', 2: 'phononTF', 3: 'phononTS'}
# ROOT Color Palettes (57 = Bird/Inferno-like, 51 = Deep Sea/Blues, 100 = Solar/Reds, 60 = Avocado/Greens)
PTYPE_PALETTE = {'All': 57, 'phononL': 51, 'phononTF': 60, 'phononTS': 100}

SLAB_X = (-4.0, 4.0)
SLAB_Y = (-4.0, 4.0)

def load_steps_to_tree(filepath):
    """Parses text data file and loads valid entries into a quick in-memory TTree."""
    tree = ROOT.TTree("phonon_tree", "Phonon Tracking Steps")
    
    # Setting up variables to bind structure pointers
    from array import array
    ptype_arr = array('i', [0])
    mx = array('f', [0.0])
    my = array('f', [0.0])
    mt = array('f', [0.0])
    
    tree.Branch("ptype", ptype_arr, "ptype/I")
    tree.Branch("mx", mx, "mx/F")
    tree.Branch("my", my, "my/F")
    tree.Branch("mt", mt, "mt/F")
    
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Missing input tracking file: {filepath}")
        
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            p = line.split()
            if len(p) < 17:
                continue
            try:
                p_name = p[3]
                if p_name not in PTYPE_MAP:
                    continue
                ptype_arr[0] = PTYPE_MAP[p_name]
                
                # Compute midpoints for spatial tracking and timestamps
                mx[0] = 0.5 * (float(p[4]) + float(p[9]))
                my[0] = 0.5 * (float(p[5]) + float(p[10]))
                mt[0] = 0.5 * (float(p[14]) + float(p[15]))
                
                tree.Fill()
            except (ValueError, IndexError):
                continue
                
    return tree

# ── synthetic data generator ──────────────────────────────────────────────────
def make_synthetic_tree():
    """Generates synthetic ROOT TTree simulating tracking steps within the boundaries."""
    import random
    random.seed(42)
    
    tree = ROOT.TTree("phonon_tree", "Synthetic Phonon Tracking Steps")
    from array import array
    ptype_arr = array('i', [0])
    mx = array('f', [0.0])
    my = array('f', [0.0])
    mt = array('f', [0.0])
    
    tree.Branch("ptype", ptype_arr, "ptype/I")
    tree.Branch("mx", mx, "mx/F")
    tree.Branch("my", my, "my/F")
    tree.Branch("mt", mt, "mt/F")
    
    n_tracks = 60
    steps_per_track = 250
    
    for tid in range(1, n_tracks + 1):
        ptype_val = random.choice([1, 2, 3])
        x = random.uniform(-0.5, 0.5)
        y = random.uniform(-0.5, 0.5)
        z = random.uniform(-0.25, 0.25)
        
        # Velocity unit vector
        vx, vy, vz = random.gauss(0,1), random.gauss(0,1), random.gauss(0,1)
        norm = (vx**2 + vy**2 + vz**2)**0.5
        vx, vy, vz = vx/norm, vy/norm, vz/norm
        t = 0.0
        
        for _ in range(steps_per_track):
            dl = random.expovariate(1.0 / 0.5)
            dt = dl / 5.0
            nx, ny, nz = x + vx * dl, y + vy * dl, z + vz * dl
            
            # Boundary Reflections
            if nx > 4.0: nx = 8.0 - nx; vx = -vx
            if nx < -4.0: nx = -8.0 - nx; vx = -vx
            if ny > 4.0: ny = 8.0 - ny; vy = -vy
            if ny < -4.0: ny = -8.0 - ny; vy = -vy
            if nz > 0.26: nz = 0.52 - nz; vz = -vz
            if nz < -0.26: nz = -0.52 - nz; vz = -vz
            
            ptype_arr[0] = ptype_val
            mx[0] = 0.5 * (x + nx)
            my[0] = 0.5 * (y + ny)
            mt[0] = 0.5 * (t + (t + dt))
            tree.Fill()
            
            x, y, z, t = nx, ny, nz, t + dt
            
    return tree

# ── animation engine ──────────────────────────────────────────────────────────
def animate_root(tree, nx, ny, tmin, tmax, nframes, output_path):
    # Establish time step bins
    t_edges = [tmin + i * (tmax - tmin) / nframes for i in range(nframes + 1)]
    
    # Identify active types present inside our current tree profile
    active_modes = []
    for code, name in PTYPE_NAME.items():
        if tree.GetEntries(f"ptype == {code}") > 0:
            active_modes.append(name)
            
    n_pads = 1 + len(active_modes)
    
    # Create the horizontal master Canvas
    canvas = ROOT.TCanvas("c1", "Phonon Density Progression", 450 * n_pads, 450)
    canvas.Divide(n_pads, 1)
    
    # Pre-allocate spatial TH2D structures for active views
    h_all = ROOT.TH2D("h_all", "All Phonons;X [mm];Y [mm]", nx, SLAB_X[0], SLAB_X[1], ny, SLAB_Y[0], SLAB_Y[1])
    h_modes = {}
    for name in active_modes:
        h_modes[name] = ROOT.TH2D(f"h_{name}", f"{name};X [mm];Y [mm]", nx, SLAB_X[0], SLAB_X[1], ny, SLAB_Y[0], SLAB_Y[1])

    # Find maximum scale value to stabilize color scaling bar across frame loops
    # Pulls max single-slice value directly using temporary trees
    h_max_calc = ROOT.TH2D("h_max_calc", "", nx, SLAB_X[0], SLAB_X[1], ny, SLAB_Y[0], SLAB_Y[1])
    max_z_all = 1.0
    max_z_modes = {n: 1.0 for n in active_modes}
    
    for i in range(nframes):
        t_low, t_high = t_edges[i], t_edges[i+1]
        t_cut = f"mt >= {t_low} && mt < {t_high}"
        
        h_max_calc.Reset()
        tree.Project("h_max_calc", "my:mx", t_cut)
        max_z_all = max(max_z_all, h_max_calc.GetMaximum())
        
        for code, name in PTYPE_MAP.items():
            if name in h_modes:
                h_max_calc.Reset()
                tree.Project("h_max_calc", "my:mx", f"{t_cut} && ptype == {code}")
                max_z_modes[name] = max(max_z_modes[name], h_max_calc.GetMaximum())

    # Compile the frame animations loops
    for i in range(nframes):
        t_low, t_high = t_edges[i], t_edges[i+1]
        t_center = 0.5 * (t_low + t_high)
        t_cut = f"mt >= {t_low} && mt < {t_high}"
        
        # 1. Update Global Pad
        canvas.cd(1)
        ROOT.gProxyInterface = None # Keep layout context pristine
        ROOT.gStyle.SetPalette(PTYPE_PALETTE['All'])
        h_all.Reset()
        tree.Project("h_all", "my:mx", t_cut)
        h_all.SetMaximum(max_z_all * 1.05)
        h_all.Draw("COLZ")
        
        # Display localized running timestamp on the first viewport
        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.045)
        latex.SetTextColor(ROOT.kWhite)
        latex.DrawLatex(0.15, 0.82, f"t = {t_center:.1f} ns")
        
        # 2. Update Channel-Specific Pads
        for idx, name in enumerate(active_modes):
            canvas.cd(idx + 2)
            ROOT.gStyle.SetPalette(PTYPE_PALETTE[name])
            h_obj = h_modes[name]
            h_obj.Reset()
            p_code = PTYPE_MAP[name]
            tree.Project(h_obj.GetName(), "my:mx", f"{t_cut} && ptype == {p_code}")
            h_obj.SetMaximum(max_z_modes[name] * 1.05)
            h_obj.Draw("COLZ")
            
        canvas.Update()
        
        # Append frame output into target storage
        if i == 0:
            canvas.Print(f"{output_path}(")
        elif i == nframes - 1:
            canvas.Print(f"{output_path})")
        else:
            canvas.Print(output_path)
            
    print(f"ROOT Animation tracking sequence closed successfully → Saved to: {output_path}")

# ── main execution execution ──────────────────────────────────────────────────
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

    # Clean out file extension mapping; ROOT natively spits out animated GIFs via standard file prints
    output_target = str(Path(args.output).with_suffix('.gif'))

    if args.synthetic or not Path(args.input).exists():
        print("Using standard synthetic data generation tree module...")
        tree = make_synthetic_tree()
    else:
        print(f"Processing tracking log entries from: {args.input} ...")
        tree = load_steps_to_tree(args.input)

    total_entries = tree.GetEntries()
    if total_entries == 0:
        print("Error: Selected data tree structural frame contains zero trackable hits.", file=sys.stderr)
        sys.exit(1)
        
    tmax = args.tmax or tree.GetMaximum("mt")
    print(f"Configuring canvas grid profile: {total_entries} hits mapped across t = [{args.tmin:.2f}, {tmax:.2f}] ns")
    
    animate_root(tree, args.nx, args.ny, args.tmin, tmax, args.nframes, output_target)

if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
g4cmp_lattice_params.py
───────────────────────
Computes all G4CMP config.txt parameters for a given material from
measured physical properties, and writes a ready-to-use config file.

Also lists all available lattices in a given G4LATTICEDATA directory
(either local or inside a Singularity container).

Usage:
    # Compute parameters for fused silica (built-in preset)
    python g4cmp_lattice_params.py --material SiO2

    # Compute parameters for Si (built-in preset, for verification)
    python g4cmp_lattice_params.py --material Si

    # Custom material from measured velocities
    python g4cmp_lattice_params.py \\
        --name MyMaterial \\
        --density 2200 \\
        --vL 5968 --vT 3764 \\
        --debye-temp 470 \\
        --scat 2.43e-41 \\
        --decay 2.5e-56 \\
        --bandgap 9.0 \\
        --epsilon 3.9

    # List available lattices in the Singularity container
    python g4cmp_lattice_params.py --list-lattices

    # List lattices in your local CrystalMaps
    python g4cmp_lattice_params.py --list-lattices \\
        --lattice-dir /path/to/CrystalMaps

    # Compute + write config to a directory
    python g4cmp_lattice_params.py --material SiO2 \\
        --output-dir ./CrystalMaps/SiO2
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

# ── Physical constants ────────────────────────────────────────────────────────
kB   = 1.380649e-23    # J/K
hbar = 1.054571e-34    # J·s
h    = 2 * np.pi * hbar

# ── Built-in material presets ─────────────────────────────────────────────────
# Values from literature; extend this dict to add new materials.
PRESETS = {
    "Si": {
        "name":        "Silicon (crystalline)",
        "density":     2329.0,    # kg/m³
        "vL":          9000.0,    # m/s  longitudinal
        "vT":          5400.0,    # m/s  transverse (average of FT and ST)
        "debye_temp":  645.0,     # K
        "scat":        2.43e-42,  # s³   Tamura isotopic
        "decay":       7.41e-56,  # s⁴   anharmonic
        "decayTT":     0.74,
        "LDOS":        0.093,
        "STDOS":       0.531,
        "FTDOS":       0.376,
        "bandgap":     1.17,      # eV
        "pair_energy": 3.81,      # eV
        "fano":        0.15,
        "epsilon":     11.68,
        "lattice_const": 5.431,   # Å  cubic
        "c11":         165.6e9,   # Pa
        "c12":          63.9e9,
        "c44":          79.5e9,
        "dyn":         (-42.9, -94.5, 52.4, 68.0),  # GPa
        "gruneisen":   1.1,
        "amorphous":   False,
        "kappa_ref":   148.0,     # W/m/K at 300K
    },
    "SiO2": {
        "name":        "Fused silica (amorphous SiO2)",
        "density":     2200.0,
        "vL":          5968.0,    # Brückner (1970)
        "vT":          3764.0,    # Krause & Kurkjian (1968)
        "debye_temp":  470.0,     # Zeller & Pohl (1971)
        "scat":        2.43e-41,  # ~10× Si, amorphous disorder
        "decay":       None,      # computed from Grüneisen
        "decayTT":     0.5,
        "LDOS":        None,      # computed from velocities
        "STDOS":       None,
        "FTDOS":       None,
        "bandgap":     9.0,       # eV  insulator
        "pair_energy": 17.0,      # eV
        "fano":        0.15,
        "epsilon":     3.9,
        "lattice_const": 4.913,   # Å  (alpha-quartz reference)
        "c11":         None,      # computed from vL, density
        "c12":         None,
        "c44":         None,
        "dyn":         None,      # computed (isotropic)
        "gruneisen":   0.7,
        "amorphous":   True,
        "kappa_ref":   1.4,       # W/m/K at 300K
    },
    "Ge": {
        "name":        "Germanium (crystalline)",
        "density":     5323.0,
        "vL":          5400.0,
        "vT":          3550.0,
        "debye_temp":  374.0,
        "scat":        3.67e-42,
        "decay":       5.77e-56,
        "decayTT":     0.74,
        "LDOS":        0.097,
        "STDOS":       0.540,
        "FTDOS":       0.363,
        "bandgap":     0.67,
        "pair_energy": 2.96,
        "fano":        0.13,
        "epsilon":     16.0,
        "lattice_const": 5.658,
        "c11":         128.9e9,
        "c12":          48.3e9,
        "c44":          67.1e9,
        "dyn":         (-33.5, -75.4, 44.7, 55.8),
        "gruneisen":   1.0,
        "amorphous":   False,
        "kappa_ref":   60.0,
    },
}


# ── Core calculations ─────────────────────────────────────────────────────────

def compute_elastic(rho, vL, vT):
    """Derive elastic constants from density and sound velocities."""
    c11 = rho * vL**2
    c44 = rho * vT**2
    c12 = c11 - 2*c44   # isotropic relation
    return c11, c12, c44


def compute_debye_freq(debye_temp):
    """Convert Debye temperature to frequency in THz."""
    omega_D = kB * debye_temp / hbar
    f_D_THz = omega_D / (2 * np.pi * 1e12)
    return f_D_THz, omega_D


def compute_dos(vL, vT):
    """
    Compute DOS fractions from velocities using Debye model.
    For isotropic: DOS_i ∝ 1/v_i³, two degenerate transverse modes.
    """
    D_L  = 1.0 / vL**3
    D_T  = 1.0 / vT**3
    total = D_L + 2*D_T
    return D_L/total, D_T/total, D_T/total   # LDOS, STDOS, FTDOS


def compute_dyn_isotropic(c11, c44):
    """
    Dynamical matrix coefficients for isotropic medium.
    A = B = -(c11-c44),  C = c44,  D = 0
    """
    A = -(c11 - c44) / 1e9
    B = -(c11 - c44) / 1e9
    C =  c44 / 1e9
    D =  0.0
    return A, B, C, D


def compute_decay_from_gruneisen(B_Si, gamma_mat, gamma_Si, vL_mat, vL_Si):
    """
    Scale anharmonic decay rate from Si reference using Grüneisen parameter.
    B_mat ~ B_Si * (γ_mat/γ_Si)² * (v_Si/v_mat)³
    """
    return B_Si * (gamma_mat/gamma_Si)**2 * (vL_Si/vL_mat)**3


def compute_mfp(vL, vT, A_scat, omega):
    """Mean free path from isotopic scattering rate at given angular frequency."""
    rate = A_scat * omega**4
    mfp_L = vL / rate
    mfp_T = vT / rate
    return mfp_L, mfp_T


def estimate_kappa(rho, vL, vT, A_scat, T=300):
    """
    Rough kinetic theory estimate of thermal conductivity.
    κ ≈ (1/3) * Σ_i C_vi * v_i * l_i
    """
    omega_peak = kB * T / hbar
    mfp_L, mfp_T = compute_mfp(vL, vT, A_scat, omega_peak)
    C_v = 740 * rho   # J/m³/K  approximate
    kappa = (1/3) * C_v * (vL*mfp_L + 2*vT*mfp_T) / 3
    return kappa, mfp_L, mfp_T


# ── Parameter builder ─────────────────────────────────────────────────────────

def build_params(args):
    """
    Build the full parameter set from either a preset or command-line values.
    Returns a dict with all G4CMP config.txt fields.
    """
    if args.material and args.material in PRESETS:
        p = PRESETS[args.material].copy()
    else:
        # Build from command-line arguments
        p = {
            "name":          args.name or "Custom",
            "density":       args.density,
            "vL":            args.vL,
            "vT":            args.vT,
            "debye_temp":    args.debye_temp,
            "scat":          args.scat,
            "decay":         args.decay,
            "decayTT":       args.decayTT,
            "LDOS":          None,
            "STDOS":         None,
            "FTDOS":         None,
            "bandgap":       args.bandgap,
            "pair_energy":   args.pair_energy or args.bandgap * 2.0,
            "fano":          args.fano,
            "epsilon":       args.epsilon,
            "lattice_const": args.lattice_const or 5.0,
            "c11":           None,
            "c12":           None,
            "c44":           None,
            "dyn":           None,
            "gruneisen":     args.gruneisen,
            "amorphous":     args.amorphous,
            "kappa_ref":     None,
        }

    rho = p["density"]
    vL  = p["vL"]
    vT  = p["vT"]

    # Elastic constants
    if p["c11"] is None:
        p["c11"], p["c12"], p["c44"] = compute_elastic(rho, vL, vT)

    # Dynamical matrix
    if p["dyn"] is None:
        p["dyn"] = compute_dyn_isotropic(p["c11"], p["c44"])

    # DOS fractions
    if p["LDOS"] is None:
        p["LDOS"], p["STDOS"], p["FTDOS"] = compute_dos(vL, vT)

    # Debye frequency
    p["debye_THz"], p["omega_D"] = compute_debye_freq(p["debye_temp"])

    # Anharmonic decay if not provided
    if p["decay"] is None:
        p["decay"] = compute_decay_from_gruneisen(
            B_Si=7.41e-56,
            gamma_mat=p["gruneisen"],
            gamma_Si=PRESETS["Si"]["gruneisen"],
            vL_mat=vL,
            vL_Si=PRESETS["Si"]["vL"]
        )

    return p


# ── Output ────────────────────────────────────────────────────────────────────

def print_report(p):
    """Print full parameter report to terminal."""
    SEP = "─" * 60

    print(f"\n{'═'*60}")
    print(f"  G4CMP LATTICE PARAMETERS: {p['name']}")
    print(f"{'═'*60}")

    print(f"\n{'── INPUT / DERIVED PHYSICAL PROPERTIES ':─<60}")
    print(f"  Density          : {p['density']:.1f} kg/m³")
    print(f"  vL (longitudinal): {p['vL']:.1f} m/s")
    print(f"  vT (transverse)  : {p['vT']:.1f} m/s")
    print(f"  Debye temperature: {p['debye_temp']:.1f} K")
    print(f"  Debye frequency  : {p['debye_THz']:.2f} THz")
    print(f"  Grüneisen param  : {p['gruneisen']:.2f}")
    print(f"  Amorphous        : {p['amorphous']}")

    print(f"\n{'── ELASTIC CONSTANTS ':─<60}")
    print(f"  c11 = {p['c11']/1e9:.2f} GPa")
    print(f"  c12 = {p['c12']/1e9:.2f} GPa")
    print(f"  c44 = {p['c44']/1e9:.2f} GPa")
    print(f"  Bulk modulus K   = {(p['c11']+2*p['c12'])/3/1e9:.2f} GPa")
    print(f"  Shear modulus G  = {p['c44']/1e9:.2f} GPa")

    A,B,C,D = p["dyn"]
    print(f"\n{'── DYNAMICAL MATRIX ':─<60}")
    print(f"  dyn {A:.2f} {B:.2f} {C:.2f} {D:.2f} GPa")

    print(f"\n{'── PHONON SCATTERING ':─<60}")
    print(f"  Isotopic scat A  : {p['scat']:.3e} s³")
    print(f"  Anharmonic B     : {p['decay']:.3e} s⁴")
    print(f"  decayTT fraction : {p['decayTT']:.2f}")

    # MFP at Debye frequency
    mfp_L_D, mfp_T_D = compute_mfp(p['vL'], p['vT'], p['scat'], p['omega_D'])
    print(f"  MFP_L at Debye   : {mfp_L_D*1e9:.2f} nm")
    print(f"  MFP_T at Debye   : {mfp_T_D*1e9:.2f} nm")

    # MFP at 1 THz
    omega_1THz = 2*np.pi*1e12
    mfp_L_1, mfp_T_1 = compute_mfp(p['vL'], p['vT'], p['scat'], omega_1THz)
    print(f"  MFP_L at 1 THz   : {mfp_L_1*1e6:.2f} µm")
    print(f"  MFP_T at 1 THz   : {mfp_T_1*1e6:.2f} µm")

    # Thermal conductivity estimate
    kappa, _, _ = estimate_kappa(p['density'], p['vL'], p['vT'], p['scat'])
    ref = f"  (literature: {p['kappa_ref']:.1f} W/m/K)" if p.get('kappa_ref') else ""
    print(f"  Est. κ at 300 K  : {kappa:.2f} W/m/K{ref}")

    print(f"\n{'── DENSITY OF STATES ':─<60}")
    print(f"  LDOS  = {p['LDOS']:.4f}")
    print(f"  STDOS = {p['STDOS']:.4f}")
    print(f"  FTDOS = {p['FTDOS']:.4f}")
    print(f"  Sum   = {p['LDOS']+p['STDOS']+p['FTDOS']:.4f}  (must = 1.000)")

    print(f"\n{'── CHARGE CARRIER PARAMETERS ':─<60}")
    print(f"  Band gap         : {p['bandgap']:.2f} eV")
    print(f"  Pair energy      : {p['pair_energy']:.2f} eV")
    print(f"  Fano factor      : {p['fano']:.2f}")
    print(f"  Permittivity ε   : {p['epsilon']:.2f}")

    if p.get('kappa_ref'):
        match = "✓ GOOD" if abs(kappa - p['kappa_ref'])/p['kappa_ref'] < 0.5 \
                else "✗ TUNE scat parameter"
        print(f"\n{'── VALIDATION ':─<60}")
        print(f"  κ estimate  : {kappa:.2f} W/m/K")
        print(f"  κ reference : {p['kappa_ref']:.2f} W/m/K")
        print(f"  Status      : {match}")


def write_config(p, output_dir):
    """Write G4CMP config.txt to output_dir."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    config_path = out / "config.txt"

    A,B,C,D = p["dyn"]
    kappa, _, _ = estimate_kappa(p['density'], p['vL'], p['vT'], p['scat'])

    lines = [
        f"# {p['name']}",
        f"# Generated by g4cmp_lattice_params.py",
        f"# Estimated κ = {kappa:.2f} W/m/K  "
        f"(ref: {p.get('kappa_ref','?')} W/m/K)",
        f"",
        f"cubic {p['lattice_const']:.3f} Ang",
        f"",
        f"stiffness 1 1 {p['c11']/1e9:.2f} GPa",
        f"stiffness 1 2 {p['c12']/1e9:.2f} GPa",
        f"stiffness 4 4 {p['c44']/1e9:.2f} GPa",
        f"",
        f"dyn {A:.2f} {B:.2f} {C:.2f} {D:.2f} GPa",
        f"",
        f"scat  {p['scat']:.3e} s3",
        f"decay {p['decay']:.3e} s4",
        f"decayTT {p['decayTT']:.2f}",
        f"",
        f"LDOS  {p['LDOS']:.4f}",
        f"STDOS {p['STDOS']:.4f}",
        f"FTDOS {p['FTDOS']:.4f}",
        f"",
        f"Debye {p['debye_THz']:.2f} THz",
        f"",
        f"bandgap    {p['bandgap']:.2f} eV",
        f"pairEnergy {p['pair_energy']:.2f} eV",
        f"fanoFactor {p['fano']:.2f}",
        f"vsound {p['vL']:.1f} m/s",
        f"vtrans {p['vT']:.1f} m/s",
        f"",
        f"l0_e 1e-9 m",
        f"l0_h 1e-9 m",
        f"",
        f"hmass 10.0",
        f"emass 10.0 10.0 10.0",
        f"",
        f"valleyDir 1 0 0",
        f"",
        f"alpha      0.0 /eV",
        f"acDeform   0.0 eV",
        f"ivDeform   0.0 0.0 0.0 0.0 0.0 0.0 eV/cm",
        f"ivEnergy   0.0 0.0 0.0 0.0 0.0 0.0 eV",
        f"neutDens   0.0 /cm3",
        f"epsilon    {p['epsilon']:.2f}",
        f"",
        f"ivModel    Linear",
        f"ivLinRate0 0.0 Hz",
        f"ivLinRate1 0.0 Hz",
        f"ivLinPower 1.0",
        f"ivQuadRate  0.0 Hz",
        f"ivQuadField 0.0 V/m",
        f"ivQuadPower 1.0",
    ]

    config_path.write_text("\n".join(lines) + "\n")
    print(f"\nConfig written → {config_path}")


# ── Lattice listing ───────────────────────────────────────────────────────────

def list_lattices(lattice_dir, singularity_sandbox):
    """List all available lattices in a CrystalMaps directory."""

    if singularity_sandbox:
        # Run inside container
        cmd = [
            "singularity", "exec", singularity_sandbox,
            "bash", "-c",
            f'echo "G4LATTICEDATA=$G4LATTICEDATA" && '
            f'ls -1 $G4LATTICEDATA 2>/dev/null || '
            f'echo "G4LATTICEDATA not set or empty"'
        ]
        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            output = result.stdout.strip()
            print(f"\n{'═'*60}")
            print(f"  LATTICES IN SINGULARITY CONTAINER")
            print(f"{'═'*60}")
            print(output)

            # Also show config summary for each
            cmd2 = [
                "singularity", "exec", singularity_sandbox,
                "bash", "-c",
                'for d in $G4LATTICEDATA/*/; do '
                '  name=$(basename $d); '
                '  echo ""; '
                '  echo "  [$name]"; '
                '  grep -E "^(cubic|stiffness 1 1|Debye|bandgap|density|vsound)" '
                '    $d/config.txt 2>/dev/null | sed "s/^/    /"; '
                'done'
            ]
            result2 = subprocess.run(cmd2, capture_output=True, text=True)
            print(result2.stdout)
        except FileNotFoundError:
            print("ERROR: singularity not found in PATH")
        return

    # Local directory
    search_dir = lattice_dir or os.environ.get("G4LATTICEDATA")
    if not search_dir:
        print("ERROR: provide --lattice-dir or set $G4LATTICEDATA")
        sys.exit(1)

    search_path = Path(search_dir)
    if not search_path.exists():
        print(f"ERROR: directory not found: {search_path}")
        sys.exit(1)

    lattices = sorted([d for d in search_path.iterdir() if d.is_dir()])

    print(f"\n{'═'*60}")
    print(f"  AVAILABLE LATTICES IN: {search_path}")
    print(f"{'═'*60}")
    print(f"  Found {len(lattices)} lattice(s):\n")

    for lat_dir in lattices:
        config = lat_dir / "config.txt"
        print(f"  [{lat_dir.name}]")
        if config.exists():
            for line in config.read_text().splitlines():
                line = line.strip()
                if line and not line.startswith("#"):
                    key = line.split()[0] if line.split() else ""
                    if key in ("cubic","stiffness","Debye","bandgap","vsound","density"):
                        print(f"    {line}")
        else:
            print(f"    (no config.txt found)")
        print()


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Compute G4CMP lattice parameters and list available lattices",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    # Mode
    parser.add_argument("--material",   choices=list(PRESETS.keys()),
                        help="Use a built-in material preset")
    parser.add_argument("--list-lattices", action="store_true",
                        help="List all available lattices")

    # Custom material
    parser.add_argument("--name",         help="Material name")
    parser.add_argument("--density",      type=float, help="kg/m³")
    parser.add_argument("--vL",           type=float, help="Longitudinal velocity m/s")
    parser.add_argument("--vT",           type=float, help="Transverse velocity m/s")
    parser.add_argument("--debye-temp",   type=float, help="Debye temperature K")
    parser.add_argument("--scat",         type=float, help="Isotopic scattering A (s³)")
    parser.add_argument("--decay",        type=float, default=None,
                        help="Anharmonic decay B (s⁴), computed if omitted")
    parser.add_argument("--decayTT",      type=float, default=0.5)
    parser.add_argument("--bandgap",      type=float, default=9.0, help="eV")
    parser.add_argument("--pair-energy",  type=float, default=None, help="eV")
    parser.add_argument("--fano",         type=float, default=0.15)
    parser.add_argument("--epsilon",      type=float, default=3.9)
    parser.add_argument("--lattice-const",type=float, default=None, help="Å")
    parser.add_argument("--gruneisen",    type=float, default=0.7)
    parser.add_argument("--amorphous",    action="store_true")

    # Output
    parser.add_argument("--output-dir",   help="Write config.txt to this directory")
    parser.add_argument("--lattice-dir",  help="Local CrystalMaps directory to list")
    parser.add_argument("--singularity",  default=None,
                        help="Path to Singularity sandbox to list container lattices "
                             "(e.g. ~/ubuntu-sandbox/)")

    args = parser.parse_args()

    # ── List lattices mode ────────────────────────────────────────────────────
    if args.list_lattices:
        list_lattices(args.lattice_dir, args.singularity)
        return

    # ── Compute mode ──────────────────────────────────────────────────────────
    if not args.material and not (args.density and args.vL and args.vT):
        parser.print_help()
        print("\nERROR: provide --material or --density/--vL/--vT")
        sys.exit(1)

    p = build_params(args)
    print_report(p)

    if args.output_dir:
        write_config(p, args.output_dir)

    # Suggest tuning if κ is off
    if p.get("kappa_ref"):
        kappa, _, _ = estimate_kappa(p['density'], p['vL'], p['vT'], p['scat'])
        ratio = kappa / p["kappa_ref"]
        if ratio > 2 or ratio < 0.5:
            print(f"\n  ⚠  κ estimate ({kappa:.2f}) differs from reference "
                  f"({p['kappa_ref']:.2f}) by {ratio:.1f}×")
            if ratio > 1:
                print(f"     → Increase scat by {ratio:.1f}×  "
                      f"(try scat = {p['scat']*ratio:.2e} s3)")
            else:
                print(f"     → Decrease scat by {1/ratio:.1f}×  "
                      f"(try scat = {p['scat']*ratio:.2e} s3)")


if __name__ == "__main__":
    main()

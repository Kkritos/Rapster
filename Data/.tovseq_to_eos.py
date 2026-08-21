"""
Generate a Rapster EoS table (Data/<name>.txt) from a precomputed TOV sequence
file TOVSeq_<name>.dat, whose columns are C, Mb, M, R, Rp, kl, hl, jl, rhoc.

Keeps the stable branch (up to argmax(Mb)), converts the geometrized radius to
km, and sorts by gravitational mass. Matches the reduction in
tidal_gw_analysis/eos/library.py.

Usage: python3 Data/.tovseq_to_eos.py APR4 MPA1 DD2
"""

import sys
from pathlib import Path

import numpy as np

# GM_sun/c^2 in cm (geometrized solar-mass length unit):
AGEO_LENGTH_IN_CM = 1.47684983e5

TOVSEQ_DIR = Path("/ligo/home/ligo.org/sanika.khadkikar/Projects/data/eos_data/tov_seq")
OUT_DIR = Path(__file__).resolve().parent


def convert(name, tovseq_dir=TOVSEQ_DIR, out_dir=OUT_DIR):
    """Write out_dir/<name>.txt from tovseq_dir/TOVSeq_<name>.dat."""
    src = Path(tovseq_dir) / f"TOVSeq_{name}.dat"
    ds = np.loadtxt(src, skiprows=1, delimiter=",")
    m_baryon, m_grav, rad_geom = ds[:, 1], ds[:, 2], ds[:, 3]

    # keep the ascending-mass (stable) branch only:
    i_max = int(np.argmax(m_baryon))
    m_grav, rad_geom = m_grav[:i_max + 1], rad_geom[:i_max + 1]

    radius_km = AGEO_LENGTH_IN_CM / 1.0e5 * rad_geom

    order = np.argsort(m_grav)
    m_grav, radius_km = m_grav[order], radius_km[order]
    keep = np.concatenate([[True], np.diff(m_grav) > 0])   # strictly increasing
    m_grav, radius_km = m_grav[keep], radius_km[keep]

    out = Path(out_dir) / f"{name}.txt"
    with open(out, 'w') as fh:
        fh.write("# Radius (km), Mass (Msun)\n")
        for R, M in zip(radius_km, m_grav):
            fh.write(f"{R!r} {M!r}\n")

    print(f"{name}: {m_grav.size} points, M in [{m_grav[0]:.4f}, {m_grav[-1]:.4f}], "
          f"R(1.4) = {np.interp(1.4, m_grav, radius_km):.3f} km  ->  {out}")


if __name__ == "__main__":
    for eos_name in sys.argv[1:]:
        convert(eos_name)

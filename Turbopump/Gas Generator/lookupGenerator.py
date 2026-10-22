"""
Lookup table for combustionProps over a mixture ratio / chamber pressure grid.

    build_combustion_table(np.linspace(0.5, 3.5, 31), np.linspace(20, 100, 9), gas=GAS)
    res = lookup_combustion_props(MR=0.85, Pcc=275)
    print(res.T, res.gamma, res.R)
"""

from functools import lru_cache

import numpy as np
from scipy.interpolate import RegularGridInterpolator

# Edit to point at the module that defines combustionProps
from CanterraFunctions import combustionProps, CombustionResult

PROP_NAMES = ("T", "rho", "h", "s", "cp", "cv", "gamma", "MW", "R")


def build_combustion_table(MR_values, Pcc_values, gas, graphite=None,
                           filename="combustion_table.npz", soot=False):
    """
    Evaluate combustionProps at every (MR, Pcc) grid point and save to an .npz file.

    Args:
        MR_values:  1-D increasing array of mixture ratios
        Pcc_values: 1-D increasing array of chamber pressures [psi]
        gas:        Cantera gas Solution, reused for every point
        graphite:   Cantera graphite Solution (required if soot=True)
        filename:   output file
        soot:       passed through to combustionProps
    """
    MR = np.asarray(MR_values, float)
    Pcc = np.asarray(Pcc_values, float)

    table = {k: np.empty((MR.size, Pcc.size)) for k in PROP_NAMES}
    for i, mr in enumerate(MR):
        for j, p in enumerate(Pcc):
            print(f"MR={mr:.3f}, Pcc={p:.1f}", end="\r")
            res = combustionProps(MR=mr, Pcc=p, soot=soot, gas=gas, graphite=graphite)
            for k in PROP_NAMES:
                table[k][i, j] = getattr(res, k)

    np.savez(filename, MR=MR, Pcc=Pcc, **table)
    _load_table.cache_clear()


@lru_cache
def _load_table(filename):
    data = np.load(filename)
    grid = (data["MR"], data["Pcc"])
    return {k: RegularGridInterpolator(grid, data[k]) for k in PROP_NAMES}


def lookup_combustion_props(MR, Pcc, filename="combustion_table.npz"):
    """Linearly interpolate the table at (MR, Pcc [psi]). Returns a CombustionResult."""
    interps = _load_table(filename)
    return CombustionResult(**{k: float(f((MR, Pcc))) for k, f in interps.items()})
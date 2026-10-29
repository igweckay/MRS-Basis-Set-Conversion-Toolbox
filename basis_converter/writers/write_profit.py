"""
writers/write_profit.py
Basis Set Converter — ProFit writer

Writes a list of core struct dicts to ProFit .mat format.
Confirmed structure from make_profit_basis_set_kci.m and profit_struct.mat.

ProFit .mat structure:
    basis.fid  - cell array {1 x nBasis}, each cell is complex (1 x n)
    basis.met  - cell array {1 x nBasis} of metabolite name strings
    basis.bw   - bandwidth in Hz
    basis.dt   - dwell time in seconds (= 1/bw)
    basis.f0   - Larmor frequency in Hz  (= sf * 1e6)
    basis.pts  - number of points
    basis.Bo   - field strength in Tesla

Also saves top-level variables (fid, met, bw, dt, f0, pts, Bo)
for compatibility with ProFit's load() calls.

Entry point:
    write_profit(basis_list, outpath) → one .mat file
"""

import os
import numpy as np
import scipy.io as sio

from writers.centre import fid_at_centre


PROTON_GAMMA = 42.57747892e6   # Hz/T  (confirmed from make_profit_basis_set_kci.m)


def write_profit(basis_list, outpath):
    """
    Write all metabolites to a single ProFit .mat file.

    Parameters
    ----------
    basis_list : list of core struct dicts
    outpath    : str, full path for output .mat file
    """
    outpath = os.path.abspath(outpath)
    os.makedirs(os.path.dirname(outpath) or '.', exist_ok=True)

    if not basis_list:
        raise RuntimeError("basis_list is empty")

    first = basis_list[0]
    bw    = float(first['sw'])          # Hz
    dt    = 1.0 / bw                   # s
    f0    = float(first['sf']) * 1e6   # MHz → Hz
    pts   = int(first.get('n', len(first['fid'])))
    Bo    = f0 / PROTON_GAMMA           # Tesla

    n_basis = len(basis_list)

    # Build cell arrays
    fid_cells = np.empty((1, n_basis), dtype=object)
    met_cells = np.empty((1, n_basis), dtype=object)

    for k, core in enumerate(basis_list):
        name = str(core.get('name', f'Metab_{k+1:02d}'))
        fid  = fid_at_centre(core).astype(complex)   # ProFit states no ppm centre: read at 4.65

        # ProFit expects 1 x n row vector
        fid_cells[0, k] = fid.reshape(1, -1)
        met_cells[0, k] = name
        print(f"  Written: {name}")

    # Build basis struct
    basis = {
        'fid' : fid_cells,
        'met' : met_cells,
        'bw'  : bw,
        'dt'  : dt,
        'f0'  : f0,
        'pts' : float(pts),
        'Bo'  : Bo,
    }

    # Save both nested struct and top-level variables
    # (matches make_profit_basis_set_kci.m save command)
    sio.savemat(outpath, {
        'basis' : basis,
        'fid'   : fid_cells,
        'met'   : met_cells,
        'bw'    : bw,
        'dt'    : dt,
        'f0'    : f0,
        'pts'   : float(pts),
        'Bo'    : Bo,
    })

    print(f"\n  Done. ProFit .mat written: {os.path.basename(outpath)}")
    print(f"  {n_basis} metabolites  bw={bw:.1f} Hz  f0={f0/1e6:.4f} MHz  Bo={Bo:.4f} T")


############## Command line usage #############

if __name__ == '__main__':
    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    if len(sys.argv) < 3:
        print("Usage: python write_profit.py <input_folder> <output_file>")
        sys.exit(1)

    input_path = sys.argv[1]
    outpath    = sys.argv[2]

    from readers.read_marss import read_marss_file, read_marss_folder

    if os.path.isdir(input_path):
        basis_list = read_marss_folder(input_path)
    else:
        basis_list = [read_marss_file(input_path)]

    print(f"Loaded {len(basis_list)} metabolites\n")
    write_profit(basis_list, outpath)

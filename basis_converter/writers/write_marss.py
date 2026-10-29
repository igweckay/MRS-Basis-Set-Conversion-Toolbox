"""
writers/write_marss.py
Basis Set Converter — MARSS writer

Writes core struct dicts to MARSS .mat format.
MARSS stores one metabolite per .mat file as an exptDat struct.

MARSS .mat structure (confirmed from NAA.mat):
    exptDat.fid     - complex FID vector (n x 1)
    exptDat.sw_h    - spectral width in Hz
    exptDat.sf      - Larmor frequency in MHz
    exptDat.nspecC  - number of complex points

Note: ppmCalib is NOT included (confirmed as user-added field).

Three output modes:
    'mat_individual' - one .mat per metabolite in a folder (default)
    'mat_combined'   - all metabolites in one .mat (not standard MARSS,
                       but useful for downstream tools)
    'raw'            - delegates to write_lcmodel.write_lcmodel_raw_folder

Entry points:
    write_marss_file(core, outpath)          → one .mat file
    write_marss_folder(basis_list, outdir)   → folder of .mat files
    write_marss_combined(basis_list, outpath)→ combined .mat
"""

import os
import numpy as np
import scipy.io as sio

from writers.centre import fid_at_centre


def write_marss_file(core, outpath):
    """
    Write one core struct dict to a MARSS .mat file.

    Parameters
    ----------
    core    : dict with keys: fid, sw, sf, n, name
    outpath : str, full output path
    """
    outpath = os.path.abspath(outpath)
    os.makedirs(os.path.dirname(outpath), exist_ok=True)

    fid  = fid_at_centre(core).astype(complex)   # MARSS states no ppm centre: read at 4.65
    sw   = float(core['sw'])
    sf   = float(core['sf'])
    n    = int(core.get('n', len(fid)))

    expt_dat = {
        'fid'   : fid.reshape(-1, 1),
        'sw_h'  : sw,
        'sf'    : sf,
        'nspecC': float(n),
    }

    sio.savemat(outpath, {'exptDat': expt_dat})


def write_marss_folder(basis_list, outdir):
    """
    Write all metabolites to individual MARSS .mat files in a folder.

    Parameters
    ----------
    basis_list : list of core struct dicts
    outdir     : str

    Returns
    -------
    list of str — paths to written .mat files
    """
    os.makedirs(outdir, exist_ok=True)
    written = []
    failed  = []

    for core in basis_list:
        name    = core.get('name', 'unknown')
        outpath = os.path.join(outdir, f"{name}.mat")
        try:
            write_marss_file(core, outpath)
            written.append(outpath)
            print(f"  Written: {name}.mat")
        except Exception as ex:
            failed.append(name)
            reason = ex
            print(f"  WARNING: failed for {name} — {ex}")

    if failed:
        print(f"\n  {len(failed)} failed: {', '.join(failed)}")
    if failed and not written:
        raise RuntimeError(f"Nothing written: every basis function failed, the last with: {reason}")
    print(f"\n  Done. {len(written)} MARSS .mat files written to: {outdir}")
    return written


def write_marss_combined(basis_list, outpath):
    """
    Write all metabolites into a single .mat file.
    Not standard MARSS format but useful for some workflows.

    Stores:
        fids   - complex matrix (n x n_metabs)
        names  - cell array of metabolite names
        sw_h   - spectral width in Hz
        sf     - Larmor frequency in MHz
        nspecC - number of points
    """
    outpath = os.path.abspath(outpath)
    os.makedirs(os.path.dirname(outpath) or '.', exist_ok=True)

    if not basis_list:
        raise RuntimeError("basis_list is empty")

    first  = basis_list[0]
    sw     = float(first['sw'])
    sf     = float(first['sf'])
    n      = int(first.get('n', len(first['fid'])))
    names  = [str(c.get('name', f'Metab_{i+1:02d}'))
              for i, c in enumerate(basis_list)]

    fids = np.column_stack([
        fid_at_centre(c).astype(complex)
        for c in basis_list
    ])

    sio.savemat(outpath, {
        'fids'  : fids,
        'names' : names,
        'sw_h'  : sw,
        'sf'    : sf,
        'nspecC': float(n),
    })

    print(f"  Done. Combined MARSS .mat written: {os.path.basename(outpath)}")
    print(f"  {len(basis_list)} metabolites  sw={sw:.1f} Hz  sf={sf:.4f} MHz")


################# Command line usage ################

if __name__ == '__main__':
    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    if len(sys.argv) < 3:
        print("Usage: python write_marss.py <input_folder> <output_folder> [mode]")
        print("  mode: individual (default) | combined | raw")
        sys.exit(1)

    input_path = sys.argv[1]
    outpath    = sys.argv[2]
    mode       = sys.argv[3] if len(sys.argv) > 3 else 'individual'

    from readers.read_marss import read_marss_file, read_marss_folder

    if os.path.isdir(input_path):
        basis_list = read_marss_folder(input_path)
    else:
        basis_list = [read_marss_file(input_path)]

    print(f"Loaded {len(basis_list)} metabolites\n")

    if mode == 'combined':
        write_marss_combined(basis_list, outpath)
    elif mode == 'raw':
        import sys as _sys
        _sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        from writers.write_lcmodel import write_lcmodel_raw_folder
        write_lcmodel_raw_folder(basis_list, outpath)
    else:
        write_marss_folder(basis_list, outpath)

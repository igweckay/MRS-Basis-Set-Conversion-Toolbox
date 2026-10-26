"""
writers/write_inspector.py
Basis Set Converter — INSPECTOR writer

Writes a list of core struct dicts to INSPECTOR .mat format.
INSPECTOR stores all metabolites in one .mat file as a struct array.

INSPECTOR .mat structure (confirmed from Basis_sLASER_26.mat):
    lcmBasis.sw_h     - spectral width in Hz
    lcmBasis.sf       - Larmor frequency in MHz
    lcmBasis.ppmCalib - water reference in ppm
    lcmBasis.data     - object array (n_metabs,)
                        each element: [name, linewidth, 0.025, fid, []]

Entry point:
    write_inspector(basis_list, outpath) -> one .mat file
"""

import os
import numpy as np
import scipy.io as sio


def write_inspector(basis_list, outpath, ppm_calib=4.675):
    """
    Write all metabolites to a single INSPECTOR .mat file.

    Parameters
    ----------
    basis_list : list of core struct dicts
    outpath    : str, full path for output .mat file
    ppm_calib  : float, water reference in ppm
    """
    outpath = os.path.abspath(outpath)
    os.makedirs(os.path.dirname(outpath) or '.', exist_ok=True)

    if not basis_list:
        raise RuntimeError("basis_list is empty")

    first = basis_list[0]
    sw_h  = float(first['sw'])
    sf    = float(first['sf'])
    n     = int(first.get('n', len(first['fid'])))
    n_metabs = len(basis_list)

    # Build fid matrix and name list separately
    # INSPECTOR data[i] = [name, linewidth, 0.025, fid, []]
    # Store as separate arrays to survive scipy roundtrip
    names     = []
    linewidths= []
    fid_matrix= np.zeros((n, n_metabs), dtype=complex)

    for i, core in enumerate(basis_list):
        name = str(core.get('name', f'Metab_{i+1:02d}'))
        lw   = 1.5 if core.get('linewidth') is None else float(core['linewidth'])
        fid  = np.asarray(core['fid']).ravel().astype(complex)[:n]
        names.append(name)
        linewidths.append(lw)
        fid_matrix[:len(fid), i] = fid
        print(f"  Written: {name}")

    # Build data as object array where each element has the 5-field structure
    # Use a structured approach that survives scipy
    data = np.empty((n_metabs,), dtype=object)
    for i in range(n_metabs):
        cell = np.empty((5,), dtype=object)
        cell[0] = names[i]
        cell[1] = linewidths[i]
        cell[2] = 0.025
        cell[3] = fid_matrix[:, i]
        cell[4] = np.array([])
        data[i] = cell

    # Store with both the cell array AND separate flat arrays for robustness
    lcm_basis = {
        'sw_h'     : sw_h,
        'sf'       : sf,
        'ppmCalib' : ppm_calib,
        'data'     : data,
        # Flat arrays as fallback for readers
        'names'    : np.array(names, dtype=object),
        'linewidths': np.array(linewidths),
        'fids'     : fid_matrix,
    }

    sio.savemat(outpath, {'lcmBasis': lcm_basis})
    print(f"\n  Done. INSPECTOR .mat written: {os.path.basename(outpath)}")
    print(f"  {n_metabs} metabolites  sw={sw_h:.1f} Hz  sf={sf:.4f} MHz")


######################## Command line usage #######################

if __name__ == '__main__':
    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    if len(sys.argv) < 3:
        print("Usage: python write_inspector.py <input_folder> <output_file>")
        sys.exit(1)

    input_path = sys.argv[1]
    outpath    = sys.argv[2]

    from readers.read_marss import read_marss_file, read_marss_folder

    if os.path.isdir(input_path):
        basis_list = read_marss_folder(input_path)
    else:
        basis_list = [read_marss_file(input_path)]

    print(f"Loaded {len(basis_list)} metabolites\n")
    write_inspector(basis_list, outpath)

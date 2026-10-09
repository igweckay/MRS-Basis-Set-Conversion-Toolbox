"""
readers/read_fsLmrs.py
Basis Set Converter — FSL-MRS reader

Reads a folder of FSL-MRS .json basis function files into core struct dicts.

FSL-MRS .json structure (confirmed from your documentation):
    basis_re     - real part of FID (array)
    basis_im     - imaginary part of FID (array)
    basis_dwell  - dwell time in seconds
    basis_centre - the spectrometer (central) frequency in MHz, as FSL-MRS's fsl_io reads
                   it (centralFrequency = basis_centre * 1E6); the ppm reference is not
                   stored (FSL-MRS takes its nucleus' default, 4.65 ppm for 1H)
    basis_width  - linewidth in Hz (optional, may be null)
    basis_name   - metabolite name string

Optional:
    meta         - dict with time, SimVersion etc.

Note:
    - basis_centre is sf in MHz; no ppm reference is stored (centerFreq None)
    - dwell time is in seconds: sw = 1 / basis_dwell
    - Water-referenced: centre typically 4.65 ppm

Core struct returned per metabolite:
    {
        'fid'        : np.ndarray, complex, shape (n,)
        'sw'         : float, spectral width in Hz
        'sf'         : float or None, basis_centre (MHz)
        'linewidth'  : float or None, basis_width (Hz)
        'n'          : int, number of points
        'name'       : str, metabolite name
        'centerFreq' : float, centre frequency in ppm
        'source'     : str, path to source .json file
    }

Entry point:
    read_fsLmrs_folder(folder_path) -> list of core struct dicts
    read_fsLmrs_file(path)          -> one core struct dict
"""

import os
import json
import numpy as np


def _log(*args): pass  # GUI handles display

def read_fsLmrs_file(path):
    """
    Read one FSL-MRS .json basis function file.

    Parameters
    ----------
    path : str
        Full path to a .json basis file.

    Returns
    -------
    core struct dict
    """
    path = os.path.abspath(path)

    if not os.path.isfile(path):
        raise RuntimeError(f"File not found: {path}")

    with open(path, 'r') as f:
        data = json.load(f)

    ############## Required fields #############
    basis = data.get('basis', data)   # some versions nest under 'basis' key

    basis_re    = basis.get('basis_re')
    basis_im    = basis.get('basis_im')
    basis_dwell = basis.get('basis_dwell')
    basis_name  = basis.get('basis_name', '')
    basis_centre= basis.get('basis_centre')
    basis_width = basis.get('basis_width')

    if basis_re is None or basis_im is None:
        raise RuntimeError(
            f"Missing basis_re or basis_im in {os.path.basename(path)}"
        )
    if basis_dwell is None or basis_dwell == 0:
        raise RuntimeError(
            f"Missing or zero basis_dwell in {os.path.basename(path)}"
        )

    ############## Build complex FID #############
    re  = np.asarray(basis_re,  dtype=float)
    im  = np.asarray(basis_im,  dtype=float)
    fid = re + 1j * im

    sw  = 1.0 / float(basis_dwell)
    n   = len(fid)

    sf  = float(basis_centre) if basis_centre is not None else None

    name = str(basis_name).strip("'\" ") if basis_name else \
           os.path.splitext(os.path.basename(path))[0]

    return {
        'fid'        : fid,
        'sw'         : sw,
        'sf'         : sf,
        'n'          : n,
        'name'       : name,
        'centerFreq' : None,
        'linewidth'  : float(basis_width) if basis_width is not None else None,
        'source'     : path,
    }


def read_fsLmrs_folder(folder_path):
    """
    Read all FSL-MRS .json basis files in a folder.

    Parameters
    ----------
    folder_path : str

    Returns
    -------
    list of core struct dicts
    """
    folder_path = os.path.abspath(folder_path)

    if not os.path.isdir(folder_path):
        raise RuntimeError(f"Folder not found: {folder_path}")

    json_files = sorted([
        os.path.join(folder_path, f)
        for f in os.listdir(folder_path)
        if f.lower().endswith('.json')
    ])

    if not json_files:
        raise RuntimeError(f"No .json files found in: {folder_path}")

    results = []
    failed  = []

    for path in json_files:
        try:
            core = read_fsLmrs_file(path)
            results.append(core)
            _log(f"  Loaded: {core['name']}")
        except Exception as ex:
            failed.append(os.path.basename(path))
            _log(f"  WARNING: skipping {os.path.basename(path)} — {ex}")

    if not results:
        raise RuntimeError(f"No valid FSL-MRS files loaded from: {folder_path}")

    if failed:
        _log(f"\n  {len(failed)} file(s) skipped: {', '.join(failed)}")

    _log(f"\n  Total loaded: {len(results)} metabolites")
    _log(f"  sw = {results[0]['sw']:.1f} Hz  n = {results[0]['n']}")
    _log(f"  Note: sf (Larmor freq) not stored in FSL-MRS format — "
          f"will be requested from user if needed.")
    return results

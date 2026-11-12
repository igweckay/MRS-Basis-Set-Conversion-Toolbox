"""
readers/read_marss.py
Basis Set Converter — MARSS reader

Reads MARSS basis function .mat files into the core struct format.

MARSS .mat structure (canonical, no ppmCalib):
    exptDat.fid     - complex FID vector
    exptDat.sw_h    - spectral width in Hz
    exptDat.sf      - Larmor frequency in MHz (already MHz, not Hz)
    exptDat.nspecC  - number of complex points

Core struct returned per metabolite:
    {
        'fid'    : np.ndarray, complex, shape (n,)
        'sw'     : float, spectral width in Hz
        'sf'     : float, Larmor frequency in MHz
        'n'      : int, number of points
        'name'   : str, metabolite name (from filename stem)
        'source' : str, full path to the .mat file
    }

Two entry points:
    read_marss_file(path)         -> one dict
    read_marss_folder(folder)     -> list of dicts
"""

import os
import numpy as np
import scipy.io as sio


def _log(*args): pass  # GUI handles display

############## Single file reader #############

def read_marss_file(path):
    """
    Read one MARSS .mat basis function file.

    Parameters
    ----------
    path : str
        Full path to a MARSS .mat file.

    Returns
    -------
    dict with keys: fid, sw, sf, n, name, source
    OR raises RuntimeError if the file cannot be read.
    """
    path = os.path.abspath(path)

    if not os.path.isfile(path):
        raise RuntimeError(f"File not found: {path}")

    ############## Load .mat #############
    # Try scipy first (works for -v7 and earlier)
    # Fall back to mat73 for -v7.3 (HDF5-based) files
    raw = _load_mat(path)

    ############## Find exptDat #############
    if 'exptDat' not in raw:
        raise RuntimeError(
            f"Expected 'exptDat' struct in {os.path.basename(path)} "
            f"but found keys: {list(raw.keys())}"
        )

    e = raw['exptDat']

    ############## Extract fields #############
    fid    = _get_field(e, 'fid',    path)
    sw_h   = _get_field(e, 'sw_h',   path)
    sf     = _get_field(e, 'sf',     path)
    nspecC = _get_field(e, 'nspecC', path)

    # Clean up shapes and types
    fid    = np.asarray(fid).squeeze().astype(complex)
    sw     = float(np.asarray(sw_h).squeeze())
    sf     = float(np.asarray(sf).squeeze())
    n      = int(np.asarray(nspecC).squeeze())

    # Sanity checks
    if fid.ndim != 1:
        raise RuntimeError(
            f"Expected 1D FID in {os.path.basename(path)}, "
            f"got shape {fid.shape}"
        )

    if len(fid) != n:
        # FID length doesn't match nspecC — use actual FID length
        n = len(fid)

    if sw <= 0:
        raise RuntimeError(
            f"Invalid spectral width sw_h={sw} in {os.path.basename(path)}"
        )

    if sf <= 0:
        raise RuntimeError(
            f"Invalid Larmor frequency sf={sf} in {os.path.basename(path)}"
        )

    # Metabolite name from filename stem
    name = os.path.splitext(os.path.basename(path))[0]

    # ppmCalib — water reference in ppm
    # Try multiple access patterns since scipy loads differently with different options
    ppm_calib = None
    try:
        raw2 = sio.loadmat(path, squeeze_me=True, struct_as_record=False)
        e2   = raw2['exptDat']
        if hasattr(e2, 'ppmCalib'):
            ppm_calib = float(np.asarray(e2.ppmCalib).squeeze())
    except Exception:
        pass

    return {
        'fid'      : fid,
        'sw'       : sw,
        'sf'       : sf,      # MHz — already correct, do NOT divide by 1e6
        'n'        : n,
        'name'     : name,
        'ppmCalib' : ppm_calib,
        'source'   : path,
    }


############## Folder reader #############

def read_marss_folder(folder_path):
    """
    Read all MARSS .mat basis function files in a folder.

    Parameters
    ----------
    folder_path : str
        Path to a folder containing MARSS .mat files.

    Returns
    -------
    list of core struct dicts (one per successfully loaded file)
    Prints a warning for any file that fails — does not raise.
    """
    folder_path = os.path.abspath(folder_path)

    if not os.path.isdir(folder_path):
        raise RuntimeError(f"Folder not found: {folder_path}")

    mat_files = sorted([
        os.path.join(folder_path, f)
        for f in os.listdir(folder_path)
        if f.lower().endswith('.mat')
    ])

    if not mat_files:
        raise RuntimeError(f"No .mat files found in: {folder_path}")

    results  = []
    failed   = []

    for path in mat_files:
        try:
            core = read_marss_file(path)
            results.append(core)
            _log(f"  Loaded: {core['name']}")
        except Exception as ex:
            failed.append(os.path.basename(path))
            _log(f"  WARNING: skipping {os.path.basename(path)} — {ex}")

    if not results:
        raise RuntimeError(
            f"No valid MARSS files could be loaded from: {folder_path}"
        )

    if failed:
        _log(f"\n  {len(failed)} file(s) skipped: {', '.join(failed)}")

    _log(f"\n  Total loaded: {len(results)} metabolites")
    _log(f"  sw  = {results[0]['sw']} Hz")
    _log(f"  sf  = {results[0]['sf']} MHz")
    _log(f"  n   = {results[0]['n']}")

    return results


############## Internal helpers #############

def _load_mat(path):
    """
    Load a .mat file using scipy.io.
    Falls back to mat73 for HDF5-based v7.3 files.
    Returns a dict with simplify_cells=True so structs become plain dicts.
    """
    try:
        return sio.loadmat(path, simplify_cells=True)
    except NotImplementedError:
        # v7.3 file — needs mat73
        try:
            import mat73
            return mat73.loadmat(path)
        except ImportError:
            raise RuntimeError(
                f"Cannot read {os.path.basename(path)}: it is a MATLAB v7.3 "
                f"(HDF5) file. Install mat73 to support this format:\n"
                f"    pip install mat73"
            )
    except Exception as ex:
        raise RuntimeError(
            f"Could not load {os.path.basename(path)}: {ex}"
        )


def _get_field(struct, field_name, path):
    """
    Extract a field from a struct (dict or scipy struct array).
    Raises a clear error if the field is missing.
    """
    # dict case (simplify_cells=True)
    if isinstance(struct, dict):
        if field_name not in struct:
            raise RuntimeError(
                f"Missing field '{field_name}' in exptDat "
                f"of {os.path.basename(path)}"
            )
        return struct[field_name]

    # scipy structured array case
    if hasattr(struct, 'dtype') and struct.dtype.names:
        if field_name not in struct.dtype.names:
            raise RuntimeError(
                f"Missing field '{field_name}' in exptDat "
                f"of {os.path.basename(path)}"
            )
        return struct[field_name]

    raise RuntimeError(
        f"Unrecognized exptDat structure in {os.path.basename(path)}"
    )


############## Quick test #############

if __name__ == '__main__':
    import sys

    if len(sys.argv) < 2:
        _log("Usage: python read_marss.py <path_to_mat_or_folder>")
        sys.exit(1)

    target = sys.argv[1]

    if os.path.isdir(target):
        results = read_marss_folder(target)
        _log(f"\nLoaded {len(results)} metabolites")
        for r in results:
            _log(f"  {r['name']:15s}  fid[0]={r['fid'][0]:.4e}  "
                  f"sw={r['sw']} Hz  sf={r['sf']} MHz  n={r['n']}")
    else:
        result = read_marss_file(target)
        _log(f"\nName   : {result['name']}")
        _log(f"sw     : {result['sw']} Hz")
        _log(f"sf     : {result['sf']} MHz")
        _log(f"n      : {result['n']}")
        _log(f"fid[0] : {result['fid'][0]:.6e}")
        _log(f"fid[-1]: {result['fid'][-1]:.6e}")

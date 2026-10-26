"""
readers/read_mrscloud.py
Basis Set Converter — MRSCloud reader

Reads MRSCloud .mat basis function files (one per metabolite).

MRSCloud struct fields (confirmed from your documentation):
    fids          - complex FID (n x 1)
    specs         - complex spectrum (n x 1)
    spectralwidth - bandwidth in Hz
    dwelltime     - dwell time in seconds
    n             - number of points
    linewidth     - 1
    Bo            - field strength in Tesla
    txfrq         - Larmor frequency in Hz
    t             - time axis
    ppm           - ppm axis
    sz            - [n, 1]
    name          - metabolite name string e.g. 'Ace'
    dims, flags, averages, rawAverages, subspecs, rawSubspecs, date

Very similar to FID-A but name IS present and txfrq units are the same.

Core struct returned:
    {
        'fid'    : np.ndarray, complex, shape (n,)
        'sw'     : float, spectral width in Hz
        'sf'     : float, Larmor frequency in MHz  (txfrq / 1e6)
        'n'      : int, number of points
        'name'   : str, metabolite name
        'Bo'     : float, field strength in Tesla
        'source' : str, path to source file
    }

Entry points:
    read_mrscloud_file(path)   -> one core struct dict
    read_mrscloud_folder(path) -> list of core struct dicts
"""

import os
import numpy as np
import scipy.io as sio


def _log(*args): pass  # GUI handles display

def read_mrscloud_file(path, name=None):
    """Read one MRSCloud .mat basis function file."""
    path = os.path.abspath(path)

    if not os.path.isfile(path):
        raise RuntimeError(f"File not found: {path}")

    raw  = _load_mat(path)
    s    = _find_struct(raw, path)

    fid   = np.asarray(_get(s, 'fids', path)).squeeze().astype(complex)
    sw    = float(np.asarray(_get(s, 'spectralwidth', path)).squeeze())
    n     = int(np.asarray(_get(s, 'n', path)).squeeze())
    Bo    = float(np.asarray(_get(s, 'Bo', path)).squeeze())
    txfrq = float(np.asarray(_get(s, 'txfrq', path)).squeeze())
    sf    = txfrq / 1e6

    # name field is present in MRSCloud (unlike FID-A)
    if name is None:
        name_raw = _get_optional(s, 'name', default=None)
        if name_raw is not None:
            name = str(np.asarray(name_raw).ravel()[0]).strip("'\" ")
        else:
            name = os.path.splitext(os.path.basename(path))[0]

    return {
        'fid'    : fid,
        'sw'     : sw,
        'sf'     : sf,
        'n'      : n,
        'name'   : name,
        'Bo'     : Bo,
        'centerFreq': _centre(s),
        'source' : path,
    }


def read_mrscloud_folder(folder_path):
    """Read all MRSCloud .mat files in a folder."""
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

    results = []
    failed  = []

    for path in mat_files:
        try:
            core = read_mrscloud_file(path)
            results.append(core)
            _log(f"  Loaded: {core['name']}")
        except Exception as ex:
            failed.append(os.path.basename(path))
            _log(f"  WARNING: skipping {os.path.basename(path)} — {ex}")

    if not results:
        raise RuntimeError(f"No valid MRSCloud files loaded from: {folder_path}")

    if failed:
        _log(f"\n  {len(failed)} file(s) skipped: {', '.join(failed)}")

    _log(f"\n  Total loaded: {len(results)} metabolites")
    return results


############## Helpers #############

def _load_mat(path):
    try:
        return sio.loadmat(path, simplify_cells=True)
    except NotImplementedError:
        try:
            import mat73
            return mat73.loadmat(path)
        except ImportError:
            raise RuntimeError(
                f"HDF5 v7.3 file — install mat73: pip install mat73"
            )
    except Exception as ex:
        raise RuntimeError(f"Could not load {os.path.basename(path)}: {ex}")


def _find_struct(raw, path):
    """MRSCloud files typically have the struct at the top level directly."""
    # Check if top-level keys look like MRSCloud fields
    mrscloud_keys = {'fids', 'spectralwidth', 'txfrq', 'Bo', 'n'}
    if mrscloud_keys.issubset(set(raw.keys())):
        return raw
    # Single variable
    keys = [k for k in raw if not k.startswith('_')]
    if len(keys) == 1:
        val = raw[keys[0]]
        if isinstance(val, dict):
            return val
    raise RuntimeError(
        f"Could not identify MRSCloud struct in {os.path.basename(path)}. "
        f"Found keys: {[k for k in raw if not k.startswith('_')]}"
    )


def _get(struct, field, path):
    if isinstance(struct, dict):
        if field not in struct:
            raise RuntimeError(
                f"Missing field '{field}' in {os.path.basename(path)}"
            )
        return struct[field]
    if hasattr(struct, 'dtype') and struct.dtype.names:
        if field not in struct.dtype.names:
            raise RuntimeError(
                f"Missing field '{field}' in {os.path.basename(path)}"
            )
        return struct[field]
    raise RuntimeError(f"Unrecognized struct in {os.path.basename(path)}")


def _get_optional(struct, field, default=None):
    try:
        return _get(struct, field, '')
    except RuntimeError:
        return default


def _centre(struct):
    """The ppm at 0 Hz of the ppm axis the file stores, or None: FID-A's frequency grid is
    symmetric about 0 Hz, so it is the axis' midpoint (4.65 in FID-A's own files, 2.97 in
    MRSCloud's, whose basis sets are centred near 3.0 ppm as Osprey's are)."""
    ppm = _get_optional(struct, 'ppm', default=None)
    ppm = None if ppm is None else np.asarray(ppm, dtype=float).ravel()
    if ppm is None or ppm.size < 2:
        return None
    return round(float(ppm[0] + ppm[-1]) / 2, 4)

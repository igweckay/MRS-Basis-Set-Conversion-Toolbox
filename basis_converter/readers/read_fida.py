"""
readers/read_fida.py
Basis Set Converter — FID-A reader

Reads a FID-A simulated basis function .mat file into a core struct dict.
FID-A produces one .mat file per metabolite via sim_semiLASER_shaped,
sim_press, etc.

FID-A output struct fields (checked from sim_readout.m):
    out.fids          - complex FID, column vector (n x 1)
    out.specs         - fftshift(ifft(fids))  ← ifft convention, NOT fft
    out.spectralwidth - bandwidth in Hz
    out.dwelltime     - dwell time in seconds
    out.n             - number of points
    out.linewidth     - linewidth in Hz
    out.Bo            - field strength in Tesla
    out.txfrq         - Larmor frequency in Hz (= Bo * 42577000)
    out.nucleus       - '1H' (hardcoded)
    out.gamma         - 42577000 (hardcoded)
    out.te            - echo time in ms (added by sim_semiLASER_shaped)
    out.seq           - sequence name string
    out.t             - time axis
    out.ppm           - ppm axis
    out.dims, out.flags, out.sz, out.date, etc.

Key difference from Osprey:
    FID-A specs = fftshift(ifft(fids))
    Osprey specs = fftshift(fft(fids))
    The converter handles this when writing to Osprey format.

name field is NOT set by FID-A simulator — derived from filename.

Core struct returned:
    {
        'fid'       : np.ndarray, complex, shape (n,)
        'sw'        : float, spectral width in Hz
        'sf'        : float, Larmor frequency in MHz  (txfrq / 1e6)
        'n'         : int, number of points
        'name'      : str, from filename stem
        'te'        : float or None, echo time in ms
        'seq'       : str or None
        'Bo'        : float, field strength in Tesla
        'linewidth' : float, linewidth in Hz
        'source'    : str, path to source file
    }

Entry points:
    read_fida_file(path)       -> one core struct dict
    read_fida_folder(path)     -> list of core struct dicts
"""

import os
import numpy as np
import scipy.io as sio


def _log(*args): pass  # GUI handles display

def read_fida_file(path, name=None):
    """
    Read one FID-A .mat basis function file.

    Parameters
    ----------
    path : str
        Full path to a FID-A .mat file.
    name : str, optional
        Metabolite name override. If None, derived from filename stem.

    Returns
    -------
    core struct dict
    """
    path = os.path.abspath(path)

    if not os.path.isfile(path):
        raise RuntimeError(f"File not found: {path}")

    raw = _load_mat(path)

    # FID-A saves as 'out' — find it
    struct = _find_struct(raw, path)

    ############## Extract fields #############
    fid   = np.asarray(_get(struct, 'fids', path)).squeeze().astype(complex)
    sw    = float(np.asarray(_get(struct, 'spectralwidth', path)).squeeze())
    n     = int(np.asarray(_get(struct, 'n', path)).squeeze())
    Bo    = float(np.asarray(_get(struct, 'Bo', path)).squeeze())
    txfrq = float(np.asarray(_get(struct, 'txfrq', path)).squeeze())
    sf    = txfrq / 1e6     # convert Hz → MHz

    # Optional fields
    te        = _get_optional(struct, 'te',        default=None)
    seq       = _get_optional(struct, 'seq',        default=None)
    linewidth = _get_optional(struct, 'linewidth',  default=1.0)

    if te is not None:
        te = float(np.asarray(te).squeeze())
    if linewidth is not None:
        linewidth = float(np.asarray(linewidth).squeeze())
    if seq is not None:
        seq = str(np.asarray(seq).ravel()[0]).strip("'\" ")

    if name is None:
        name = os.path.splitext(os.path.basename(path))[0]

    return {
        'fid'       : fid,
        'sw'        : sw,
        'sf'        : sf,
        'n'         : n,
        'name'      : name,
        'te'        : te,
        'seq'       : seq,
        'Bo'        : Bo,
        'linewidth' : linewidth,
        'centerFreq': _centre(struct),
        'source'    : path,
    }


def read_fida_folder(folder_path):
    """
    Read all FID-A .mat basis function files in a folder.

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
            core = read_fida_file(path)
            results.append(core)
            _log(f"  Loaded: {core['name']}")
        except Exception as ex:
            failed.append(os.path.basename(path))
            _log(f"  WARNING: skipping {os.path.basename(path)} — {ex}")

    if not results:
        raise RuntimeError(f"No valid FID-A files loaded from: {folder_path}")

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
    """Find the FID-A output struct — usually stored as 'out'."""
    for key in ['out', 'OUT']:
        if key in raw:
            return raw[key]
    # Single variable file
    keys = [k for k in raw if not k.startswith('_')]
    if len(keys) == 1:
        return raw[keys[0]]
    raise RuntimeError(
        f"Could not find FID-A 'out' struct in {os.path.basename(path)}. "
        f"Found keys: {keys}"
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
    """Get a field without raising if absent."""
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

"""
readers/read_osprey.py
Basis Set Converter — Osprey reader

Reads an Osprey .mat BASIS file into a list of core struct dicts,
one per metabolite.

Osprey .mat structure (top-level variable: BASIS):
    BASIS.fids          - complex FID matrix (n x nTotal)
    BASIS.specs         - complex spectrum matrix (n x nTotal)
    BASIS.spectralwidth - bandwidth in Hz
    BASIS.dwelltime     - dwell time in seconds
    BASIS.n             - number of points
    BASIS.Bo            - field strength in Tesla
    BASIS.te            - echo time in ms
    BASIS.seq           - sequence name cell array e.g. {'unedited'}
    BASIS.centerFreq    - center frequency in ppm
    BASIS.linewidth     - linewidth in Hz, one for the set
    BASIS.name          - cell array of metabolite names (1 x nTotal)
    BASIS.nMets         - number of metabolites (MM/Lip excluded)
    BASIS.nMM           - number of MM/Lip/H2O components
    BASIS.scale         - normalization scale factor
    BASIS.ppm           - ppm axis (n x 1)
    BASIS.t             - time axis (n x 1)
    BASIS.dims          - struct (t, coils, averages, subSpecs, extras)
    BASIS.flags         - struct (various processing flags)
    BASIS.sz            - [n, nTotal]

Note on specs convention:
    Osprey uses fftshift(fft(fids)) — NOT ifft.
    This differs from FID-A which uses fftshift(ifft(fids)).

Core struct returned per metabolite:
    {
        'fid'        : np.ndarray, complex, shape (n,)
        'sw'         : float, spectral width in Hz
        'sf'         : float, Larmor frequency in MHz  (= Bo * 42.577)
        'n'          : int, number of points
        'name'       : str, metabolite name
        'te'         : float, echo time in ms
        'seq'        : str, sequence name
        'centerFreq' : float, center frequency in ppm
        'linewidth'  : float or None, linewidth in Hz (BASIS.linewidth)
        'Bo'         : float, field strength in Tesla
        'source'     : str, path to source file
    }

Entry point:
    read_osprey(path) -> list of core struct dicts
"""

import os
import numpy as np
import scipy.io as sio


# Proton gyromagnetic ratio in MHz/T
PROTON_GAMMA = 42.577


def _log(*args): pass  # GUI handles display

def read_osprey(path):
    """
    Read an Osprey .mat BASIS file.

    Parameters
    ----------
    path : str
        Full path to an Osprey .mat file.

    Returns
    -------
    list of core struct dicts, one per metabolite
    """
    path = os.path.abspath(path)

    if not os.path.isfile(path):
        raise RuntimeError(f"File not found: {path}")

    raw = _load_mat(path)

    ############## Find the BASIS struct #############
    basis = _find_basis_struct(raw, path)

    ############## Extract shared parameters #############
    sw          = float(_get(basis, 'spectralwidth', path))
    n           = int(_get(basis, 'n', path))
    Bo          = float(_get(basis, 'Bo', path))
    te          = float(_get(basis, 'te', path))
    center_freq = float(_get(basis, 'centerFreq', path))
    try:                                   # one linewidth (Hz) for the set
        linewidth = float(np.asarray(_get(basis, 'linewidth', path)).squeeze())
    except RuntimeError:
        linewidth = None
    sf          = Bo * PROTON_GAMMA    # MHz

    # Sequence — stored as cell array {'unedited'} or string
    seq_raw = _get(basis, 'seq', path)
    seq     = _unpack_cell_string(seq_raw)

    # Names — stored as cell array
    names_raw = _get(basis, 'name', path)
    names     = _unpack_names(names_raw)

    # FIDs — shape (n, nTotal)
    fids_raw = np.asarray(_get(basis, 'fids', path))
    if fids_raw.ndim == 1:
        fids_raw = fids_raw.reshape(-1, 1)

    n_total = fids_raw.shape[1]

    if len(names) != n_total:
        raise RuntimeError(
            f"Name count ({len(names)}) does not match "
            f"FID columns ({n_total}) in {os.path.basename(path)}"
        )

    ############## Build one core struct per metabolite #############
    results = []
    for i, name in enumerate(names):
        fid = fids_raw[:, i].astype(complex)

        results.append({
            'fid'        : fid,
            'sw'         : sw,
            'sf'         : sf,
            'n'          : n,
            'name'       : name,
            'te'         : te,
            'seq'        : seq,
            'centerFreq' : center_freq,
            'linewidth'  : linewidth,
            'Bo'         : Bo,
            'source'     : path,
        })
        _log(f"  Loaded: {name}")

    _log(f"\n  Total loaded: {len(results)} metabolites")
    _log(f"  sw = {sw:.1f} Hz  sf = {sf:.4f} MHz  Bo = {Bo:.3f} T  "
          f"te = {te} ms  n = {n}")
    return results


############## Internal helpers #############

def _load_mat(path):
    try:
        return sio.loadmat(path, simplify_cells=True)
    except NotImplementedError:
        try:
            import mat73
            return mat73.loadmat(path)
        except ImportError:
            raise RuntimeError(
                f"Cannot read {os.path.basename(path)}: HDF5 v7.3 file. "
                f"Install mat73:  pip install mat73"
            )
    except Exception as ex:
        raise RuntimeError(
            f"Could not load {os.path.basename(path)}: {ex}"
        )


def _find_basis_struct(raw, path):
    """Find the BASIS struct inside the loaded .mat dict."""
    # Osprey saves as variable 'BASIS'
    if 'BASIS' in raw:
        return raw['BASIS']
    # Some versions use lowercase
    if 'basis' in raw:
        return raw['basis']
    # Single-variable file
    keys = [k for k in raw if not k.startswith('_')]
    if len(keys) == 1:
        return raw[keys[0]]
    raise RuntimeError(
        f"Could not find BASIS struct in {os.path.basename(path)}. "
        f"Found keys: {keys}"
    )


def _get(struct, field, path):
    """Get a field from a dict or scipy struct."""
    if isinstance(struct, dict):
        if field not in struct:
            raise RuntimeError(
                f"Missing field '{field}' in BASIS struct of "
                f"{os.path.basename(path)}"
            )
        return struct[field]
    if hasattr(struct, 'dtype') and struct.dtype.names:
        if field not in struct.dtype.names:
            raise RuntimeError(
                f"Missing field '{field}' in BASIS struct of "
                f"{os.path.basename(path)}"
            )
        return struct[field]
    raise RuntimeError(
        f"Unrecognized BASIS struct format in {os.path.basename(path)}"
    )


def _unpack_cell_string(val):
    """Unpack a MATLAB cell array string like {'unedited'} → 'unedited'."""
    if isinstance(val, str):
        return val.strip("'\" ")
    if isinstance(val, (list, np.ndarray)):
        flat = np.asarray(val).ravel()
        if len(flat) > 0:
            return str(flat[0]).strip("'\" ")
    return str(val)


def _unpack_names(val):
    """
    Unpack Osprey name cell array into a plain list of strings.
    Handles nested lists, numpy arrays, and plain strings.
    """
    if isinstance(val, str):
        return [val.strip()]

    if isinstance(val, np.ndarray):
        flat = val.ravel()
        return [str(x).strip("'\" ") for x in flat]

    if isinstance(val, list):
        names = []
        for x in val:
            if isinstance(x, (list, np.ndarray)):
                inner = np.asarray(x).ravel()
                names.extend([str(s).strip("'\" ") for s in inner])
            else:
                names.append(str(x).strip("'\" "))
        return names

    return [str(val).strip()]

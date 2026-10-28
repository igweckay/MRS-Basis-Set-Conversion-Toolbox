"""
readers/read_niftimrs.py
Basis Set Converter — NIfTI-MRS reader

Reads a NIfTI-MRS .nii.gz basis set file into a list of core struct dicts.

NIfTI-MRS is a NIfTI-based format for MRS data. Basis sets stored as
NIfTI-MRS have metabolites along one of the higher dimensions.

Requires: nibabel  (pip install nibabel)

NIfTI-MRS FIDs are the complex conjugate of a basis' own (as FSL-MRS reads them, NAA at
2.01 ppm): they are conjugated on reading, as write_niftimrs conjugates them on writing.

NIfTI-MRS header JSON extension contains:
    SpectralWidth        - bandwidth in Hz
    TransmitterFrequency - Larmor frequency in MHz
    EchoTime             - echo time in seconds (convert to ms)
    BasisLinewidth       - linewidth in Hz, a user-defined key BasiC writes
    ResonantNucleus      - '1H'
    dim_5, dim_6, dim_7  - dimension labels (e.g. 'DIM_COIL', 'DIM_USER_0')

The metabolite dimension is typically labelled 'DIM_USER_0' or similar.
Metabolite names are stored in the dim header as a JSON list.

Core struct returned per metabolite:
    {
        'fid'    : np.ndarray, complex, shape (n,)
        'sw'     : float, spectral width in Hz
        'sf'     : float, Larmor frequency in MHz
        'n'      : int, number of points
        'name'   : str, metabolite name or 'Metab_01' etc.
        'te'     : float or None, echo time in ms
        'source' : str, path to .nii.gz file
    }

Entry point:
    read_niftimrs(path) -> list of core struct dicts
"""

import os
import json
import numpy as np


def _log(*args): pass  # GUI handles display

def read_niftimrs(path):
    """
    Read a NIfTI-MRS basis set file.

    Parameters
    ----------
    path : str
        Full path to a .nii.gz file.

    Returns
    -------
    list of core struct dicts, one per metabolite
    """
    path = os.path.abspath(path)

    if not os.path.isfile(path):
        raise RuntimeError(f"File not found: {path}")

    try:
        import nibabel as nib
    except ImportError:
        raise RuntimeError(
            "nibabel is required to read NIfTI-MRS files.\n"
            "Install it with:  pip install nibabel"
        )

    img  = nib.load(path)
    data = np.asarray(img.dataobj)   # shape: (1, 1, 1, n, [coils], [metabs])
    hdr  = img.header

    ############## Parse JSON extension #############
    meta = _parse_nifti_extension(img)

    sw = float(meta.get('SpectralWidth', 0))

    # SpectrometerFrequency is in Hz — convert to MHz
    # Stored as array e.g. [123261801.0]
    sf_raw = meta.get('SpectrometerFrequency',
               meta.get('TransmitterFrequency', 0))
    if isinstance(sf_raw, list):
        sf_raw = sf_raw[0]
    sf = float(sf_raw) / 1e6   # Hz → MHz

    te_s = meta.get('EchoTime', None)
    te   = float(te_s) * 1000.0 if te_s is not None else None
    lw   = meta.get('BasisLinewidth')                        # a user-defined key: Value, Description
    lw   = lw.get('Value') if isinstance(lw, dict) else lw

    if sw == 0:
        raise RuntimeError(
            f"SpectralWidth not found in NIfTI-MRS header of "
            f"{os.path.basename(path)}"
        )

    ############## Find the FID dimension #############
    # NIfTI-MRS shape: (x, y, z, t, [dim5], [dim6], [dim7])
    # t dimension (index 3) is the FID/spectral dimension
    shape = data.shape
    n     = shape[3]

    ############## Extract metabolites #############
    # Squeeze spatial dims (should all be 1 for basis sets)
    # Remaining dims after t are coils/averages/metabolites
    data_squeezed = data.squeeze()   # shape: (n,) or (n, metabs) etc.
    # NIfTI-MRS FIDs are the complex conjugate of a basis' own (write_niftimrs)
    data_squeezed = np.conj(data_squeezed)

    # Derive name from filename stem (strip _conj suffix if present)
    stem = os.path.splitext(os.path.splitext(os.path.basename(path))[0])[0]
    stem = stem.replace('_conj', '')

    if data_squeezed.ndim == 1:
        # Single metabolite
        fid   = data_squeezed.astype(complex)
        names = [stem]
        fids  = [fid]
    elif data_squeezed.ndim == 2:
        # Multiple metabolites along second dim
        fids  = [data_squeezed[:, i].astype(complex)
                 for i in range(data_squeezed.shape[1])]
        names = _get_metabolite_names(meta, len(fids))
        # If names defaulted to Metab_01 etc, use stem for single metabolite
        if len(fids) == 1:
            names = [stem]
    else:
        raise RuntimeError(
            f"Unexpected NIfTI-MRS data shape after squeeze: "
            f"{data_squeezed.shape}"
        )

    results = []
    for fid, name in zip(fids, names):
        results.append({
            'fid'    : fid,
            'sw'     : sw,
            'sf'     : sf,
            'n'      : len(fid),
            'name'   : name,
            'te'     : te,
            'linewidth': None if lw is None else float(lw),
            'source' : path,
        })
        _log(f"  Loaded: {name}")

    _log(f"\n  Total loaded: {len(results)} metabolites")
    _log(f"  sw = {sw:.1f} Hz  sf = {sf:.4f} MHz  n = {n}")
    return results


############## Helpers #############

def _parse_nifti_extension(img):
    """Extract the JSON metadata from the NIfTI-MRS extension."""
    try:
        for ext in img.header.extensions:
            if ext.get_code() == 44:   # NIfTI-MRS extension code
                content = ext.get_content()
                if isinstance(content, bytes):
                    content = content.decode('utf-8').rstrip('\x00')
                return json.loads(content)
    except Exception:
        pass
    return {}


def _get_metabolite_names(meta, n_metabs):
    """
    Try to extract metabolite names from the NIfTI-MRS dim header.
    Falls back to 'Metab_01', 'Metab_02', etc.
    """
    # Check dim_5, dim_6, dim_7 for metabolite labels
    for dim_key in ['dim_5_header', 'dim_6_header', 'dim_7_header']:
        if dim_key in meta:
            try:
                labels = json.loads(meta[dim_key])
                if isinstance(labels, list) and len(labels) == n_metabs:
                    return [str(l).strip() for l in labels]
            except Exception:
                pass

    return [f"Metab_{i+1:02d}" for i in range(n_metabs)]

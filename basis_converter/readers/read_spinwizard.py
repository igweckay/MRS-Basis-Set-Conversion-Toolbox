"""
readers/read_spinwizard.py
Basis Set Converter — SpinWizard / JET reader

Reads SpinWizard basis function files as used by the JET toolkit.

Format (confirmed):
    - No file extension — filename is the metabolite name (e.g. 'NAA', 'Cho')
    - Two whitespace-separated columns: real and imaginary FID values
    - No header
    - One file per metabolite in a Basisset folder

Example file contents:
      9.70807      0.00672
      8.66908     -4.12307
      5.78875     -7.28674
      ...

Note: sw and sf are not stored in the files. They must be supplied
by the user (screen 3) or inferred from a companion config file.

Entry points:
    read_spinwizard_file(path, sw, sf) -> core struct dict
    read_spinwizard_folder(folder, sw, sf) -> list of core struct dicts
"""

import os
import numpy as np


def _read_companion_files(folder):
    """
    Read SpinWizard companion parameter files if present.

    BasisSetParameters.txt — two values:
        Line 1: Larmor frequency in MHz (sf)
        Line 2: spectral width in kHz (sw, multiply by 1000 for Hz)

    BasisRFOffset.txt — one value:
        Water reference frequency in ppm (ppmCalib)

    Returns
    -------
    dict with keys: sf, sw, ppmCalib (any may be None if file missing)
    """
    params = {'sf': None, 'sw': None, 'ppmCalib': None}

    # BasisSetParameters.txt -> sf (MHz) and sw (kHz -> Hz)
    params_path = os.path.join(folder, 'BasisSetParameters.txt')
    if os.path.isfile(params_path):
        try:
            lines = [l.strip() for l in open(params_path).readlines()
                     if l.strip()]
            if len(lines) >= 1:
                params['sf'] = float(lines[0])
            if len(lines) >= 2:
                params['sw'] = float(lines[1]) * 1000.0  # kHz -> Hz
            print(f"  BasisSetParameters: sf={params['sf']} MHz  "
                  f"sw={params['sw']} Hz")
        except Exception as ex:
            print(f"  WARNING: could not read BasisSetParameters.txt: {ex}")

    # BasisRFOffset.txt → ppmCalib
    offset_path = os.path.join(folder, 'BasisRFOffset.txt')
    if os.path.isfile(offset_path):
        try:
            val = float(open(offset_path).read().strip())
            params['ppmCalib'] = val
            print(f"  BasisRFOffset: ppmCalib={val} ppm")
        except Exception as ex:
            print(f"  WARNING: could not read BasisRFOffset.txt: {ex}")

    return params


def read_spinwizard_file(path, sw=None, sf=None, ppm_calib=None):
    """
    Read a single SpinWizard basis function file.

    Parameters
    ----------
    path : str
        Full path to the file (no extension, filename = metabolite name)
    sw   : float or None
        Spectral width in Hz (must be supplied — not in file)
    sf   : float or None
        Larmor frequency in MHz (must be supplied — not in file)

    Returns
    -------
    core struct dict
    """
    path = os.path.abspath(path)
    if not os.path.isfile(path):
        raise RuntimeError(f"File not found: {path}")

    # Metabolite name = filename with no extension
    name = os.path.basename(path)
    # Strip extension if present (robustness)
    name = os.path.splitext(name)[0]

    # Read two-column data
    try:
        data = np.loadtxt(path)
    except Exception as ex:
        raise RuntimeError(
            f"Could not read SpinWizard file {os.path.basename(path)}: {ex}"
        )

    if data.ndim == 1:
        # Single row — wrap
        data = data.reshape(1, -1)

    if data.shape[1] < 2:
        raise RuntimeError(
            f"{os.path.basename(path)}: expected 2 columns (real, imag), "
            f"got {data.shape[1]}"
        )

    fid = data[:, 0] + 1j * data[:, 1]
    n   = len(fid)

    return {
        'fid'      : fid,
        'sw'       : float(sw) if sw is not None else None,
        'sf'       : float(sf) if sf is not None else None,
        'n'        : n,
        'name'     : name,
        'ppmCalib' : float(ppm_calib) if ppm_calib is not None else None,
        'source'   : path,
    }


def read_spinwizard_folder(folder_path, sw=None, sf=None, ppm_calib=None):
    """
    Read all SpinWizard basis files in a folder.

    Automatically reads BasisSetParameters.txt and BasisRFOffset.txt
    if present. The sw/sf/ppm_calib passed as arguments are used where
    the companion files state none, as core.load passes them.

    Parameters
    ----------
    folder_path : str
    sw          : float or None — spectral width in Hz (where the files state none)
    sf          : float or None — Larmor frequency in MHz (where the files state none)
    ppm_calib   : float or None — water reference ppm (where the files state none)

    Returns
    -------
    list of core struct dicts
    """
    folder_path = os.path.abspath(folder_path)

    # Read companion parameter files
    companion = _read_companion_files(folder_path)
    sw        = companion['sw']       or sw
    sf        = companion['sf']       or sf
    ppm_calib = companion['ppmCalib'] or ppm_calib

    # Collect candidate files — no extension, not hidden, not system files
    candidates = []
    for fname in sorted(os.listdir(folder_path)):
        # Skip hidden, system, and known non-basis files
        if fname.startswith('.') or fname.startswith('._'):
            continue
        if fname in ('LCMBasisSet.txt', 'BasisSetParameters.txt',
                     'BasisRFOffset.txt', 'Thumbs.db', 'desktop.ini'):
            continue
        fpath = os.path.join(folder_path, fname)
        if not os.path.isfile(fpath):
            continue
        candidates.append(fpath)

    if not candidates:
        raise RuntimeError(f"No files found in {folder_path}")

    results = []
    for path in candidates:
        try:
            core = read_spinwizard_file(path, sw=sw, sf=sf, ppm_calib=ppm_calib)
            results.append(core)
            print(f"  Loaded: {core['name']:15s}  n={core['n']}")
        except Exception as ex:
            print(f"  WARNING: skipping {os.path.basename(path)}: {ex}")

    if not results:
        raise RuntimeError(
            f"No valid SpinWizard basis files found in {folder_path}"
        )

    print(f"\n  Total: {len(results)} metabolites loaded from {folder_path}")
    return results


############## Command line #############

if __name__ == '__main__':
    import sys

    if len(sys.argv) < 2:
        print("Usage: python read_spinwizard.py <file_or_folder> [sw_hz] [sf_mhz]")
        sys.exit(1)

    path = sys.argv[1]
    sw   = float(sys.argv[2]) if len(sys.argv) > 2 else None
    sf   = float(sys.argv[3]) if len(sys.argv) > 3 else None

    if os.path.isdir(path):
        results = read_spinwizard_folder(path, sw=sw, sf=sf)
    else:
        results = [read_spinwizard_file(path, sw=sw, sf=sf)]

    print(f"\nLoaded {len(results)} metabolites")
    for r in results:
        print(f"  {r['name']:15s}  n={r['n']}  "
              f"sw={r['sw']}  sf={r['sf']}  "
              f"|fid[0]|={abs(r['fid'][0]):.4f}")
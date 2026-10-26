"""
readers/read_lcmodel.py
Basis Set Converter — LCModel reader

Reads LCModel .BASIS files or folders of individual .raw files
into the core struct format.

.BASIS file structure:
    $SEQPAR block  — HZPPPM, NUNFIL, DELTAT, ECHOT, SEQ, FWHMBA
    $BASIS1 block  — IDBASI, FMTBAS, BADELT, NDATAB
    For each metabolite:
        $NMUSED block  — FILRAW, fitting params (ignored)
        $BASIS block   — ID, METABO, CONC, TRAMP, VOLUME, ISHIFT
        spectrum data  — real/imag pairs, FORTRAN format: the FFT of the
                         FID zero-filled to NDATAB (twice its length)

.raw file structure:
    $SEQPAR block  — HZPPPM, NUNFIL, DELTAT
    $NMID block    — ID, FMTDAT, VOLUME, TRAMP
    FID data       — real/imag pairs

Core struct returned:
    {
        'fid'    : np.ndarray, complex, shape (n,)
        'sw'     : float, spectral width in Hz  (= 1/DELTAT)
        'sf'     : float, Larmor frequency in MHz (= HZPPPM)
        'n'      : int, number of points (= NUNFIL or NDATAB)
        'name'   : str, metabolite name
        'te'     : float or None, echo time in ms (from ECHOT if present)
        'seq'    : str or None, sequence name (from SEQ if present)
        'source' : str, path to source file
    }

Two entry points:
    read_lcmodel_basis(path)      -> list of core struct dicts
    read_lcmodel_raw(path)        -> one core struct dict
    read_lcmodel_raw_folder(path) -> list of core struct dicts
"""

import os
import re
import numpy as np

# LCModel truncates metabolite names to 6 characters.
# This map restores full names where the truncation is unambiguous.
# Add entries here if new metabolites are encountered.
LCMODEL_NAME_MAP = {
    # Metabolite name restorations (LCModel truncates to 6 chars)
    'CrNo39' : 'CrNo391',
    'PCrNo3' : 'PCrNo393',
    # Sequence name restorations
    'sLASE'  : 'sLASER',
    'semi-L' : 'semi-LASER',
    'SLASER' : 'sLASER',
}

def _restore_name(name):
    """Restore full metabolite name from LCModel 6-char truncated version."""
    name = name.strip()
    return LCMODEL_NAME_MAP.get(name, name)


############## .BASIS reader #############

def _log(*args): pass  # GUI handles display

def read_lcmodel_basis(path):
    """
    Read a LCModel .BASIS file.

    Parameters
    ----------
    path : str
        Full path to a .BASIS file.

    Returns
    -------
    list of core struct dicts, one per metabolite
    """
    path = os.path.abspath(path)

    if not os.path.isfile(path):
        raise RuntimeError(f"File not found: {path}")

    with open(path, 'r') as f:
        text = f.read()

    ############## Parse $SEQPAR block #############
    seqpar = _parse_namelist(text, 'SEQPAR')
    hzpppm  = float(seqpar.get('HZPPPM', 0))
    deltat  = float(seqpar.get('DELTAT', 0)) or float(seqpar.get('BADELT', 0))
    echot   = float(seqpar.get('ECHOT',  0)) or None
    seq     = seqpar.get('SEQ', '').strip("'\" ")

    ############## Parse $BASIS1 block #############
    basis1  = _parse_namelist(text, 'BASIS1')
    ndatab  = int(float(basis1.get('NDATAB', 0)))
    badelt  = float(basis1.get('BADELT', deltat))

    # Use BADELT if DELTAT not in SEQPAR
    if deltat == 0:
        deltat = badelt

    sw = round(1.0 / deltat) if deltat > 0 else 0.0
    n  = ndatab

    ############## Split into per-metabolite sections #############
    # Each metabolite starts with $NMUSED (or $BASIS for simpler files)
    # Pattern: find all $BASIS blocks (not $BASIS1)
    basis_blocks = re.split(r'\$BASIS\b(?!1)', text, flags=re.IGNORECASE)

    results = []

    for block in basis_blocks[1:]:   # skip header before first $BASIS
        # Parse the $BASIS namelist at the start of this block
        end_idx = block.upper().find('$END')
        if end_idx == -1:
            continue

        namelist_text = block[:end_idx]
        params = _parse_inline_namelist(namelist_text)

        metabo = params.get('METABO', '').strip("' ")
        if not metabo:
            continue
        metabo = _restore_name(metabo)

        # Find the FID data — everything after $END
        data_text = block[end_idx + 4:]

        # Strip any subsequent $NMUSED block that may follow
        next_block = re.search(r'\$NMUSED|\$BASIS', data_text, re.IGNORECASE)
        if next_block:
            data_text = data_text[:next_block.start()]

        spec = _parse_fid_data(data_text, n)
        if spec is None:
            _log(f"  WARNING: could not parse FID for {metabo} — skipping")
            continue

        # The data block is the spectrum: back to the FID (first half of the inverse FFT,
        # the rest is zero-filling), imaginary sign flipped as for .raw files
        fid = np.conj(np.fft.ifft(spec, norm='ortho')[:len(spec) // 2])

        results.append({
            'fid'    : fid,
            'sw'     : sw,
            'sf'     : hzpppm,   # HZPPPM is numerically Hz/ppm = sf in MHz
            'n'      : len(fid),
            'name'   : metabo,
            'te'     : echot,
            'seq'    : seq if seq else None,
            'source' : path,
        })
        _log(f"  Loaded: {metabo}")

    if not results:
        raise RuntimeError(f"No metabolites could be parsed from: {path}")

    _log(f"\n  Total loaded: {len(results)} metabolites")
    _log(f"  sw = {sw:.1f} Hz  sf = {hzpppm:.4f} MHz  n = {n}")
    return results


############## Single .raw reader #############

def read_lcmodel_raw(path, name=None, sw=None, sf=None):
    """
    Read one LCModel .raw basis function file.

    Parameters
    ----------
    path : str
        Full path to a .raw file.
    name : str, optional
        Metabolite name. If None, derived from filename stem.

    Returns
    -------
    core struct dict
    """
    path = os.path.abspath(path)

    if not os.path.isfile(path):
        raise RuntimeError(f"File not found: {path}")

    with open(path, 'r') as f:
        text = f.read()

    ############## Parse $SEQPAR block #############
    seqpar  = _parse_namelist(text, 'SEQPAR')
    hzpppm  = float(seqpar.get('HZPPPM', 0))
    nunfil  = int(float(seqpar.get('NUNFIL', 0)))
    deltat  = float(seqpar.get('DELTAT', 0))

    sw = round(1.0 / deltat) if deltat > 0 else 0.0

    ############## Find FID data after last $END #############
    last_end = text.upper().rfind('$END')
    if last_end == -1:
        raise RuntimeError(f"No $END found in {os.path.basename(path)}")

    data_text = text[last_end + 4:]
    fid = _parse_fid_data(data_text, nunfil)

    if fid is None:
        raise RuntimeError(
            f"Could not parse FID data from {os.path.basename(path)}"
        )

    # .raw files store the imaginary sign flipped (see write_lcmodel_raw): flip it back
    fid = np.conj(fid)

    if name is None:
        name = os.path.splitext(os.path.basename(path))[0]

    # Detect MARSS-style .raw (no $SEQPAR — sw/sf not stored in file)
    is_marss_raw = '$SEQPAR' not in text.upper()

    # Use sw/sf from file if present, else use supplied values or 0
    sw_out = round(1.0 / deltat) if deltat > 0 else (float(sw) if sw else 0.0)
    sf_out = hzpppm if hzpppm > 0 else (float(sf) if sf else 0.0)

    return {
        'fid'          : fid,
        'sw'           : sw_out,
        'sf'           : sf_out,
        'n'            : len(fid),
        'name'         : name,
        'te'           : None,
        'seq'          : None,
        'is_marss_raw' : is_marss_raw,
        'source'       : path,
    }


############## Folder of .raw files reader #############

def read_lcmodel_raw_folder(folder_path):
    """
    Read all .raw basis function files in a folder.

    Parameters
    ----------
    folder_path : str
        Path to folder containing .raw files.

    Returns
    -------
    list of core struct dicts
    """
    folder_path = os.path.abspath(folder_path)

    if not os.path.isdir(folder_path):
        raise RuntimeError(f"Folder not found: {folder_path}")

    raw_files = sorted([
        os.path.join(folder_path, f)
        for f in os.listdir(folder_path)
        if f.lower().endswith('.raw')
    ])

    if not raw_files:
        raise RuntimeError(f"No .raw files found in: {folder_path}")

    results = []
    failed  = []

    for path in raw_files:
        try:
            core = read_lcmodel_raw(path)
            results.append(core)
            _log(f"  Loaded: {core['name']}")
        except Exception as ex:
            failed.append(os.path.basename(path))
            _log(f"  WARNING: skipping {os.path.basename(path)} — {ex}")

    if not results:
        raise RuntimeError(
            f"No valid .raw files could be loaded from: {folder_path}"
        )

    if failed:
        _log(f"\n  {len(failed)} file(s) skipped: {', '.join(failed)}")

    _log(f"\n  Total loaded: {len(results)} metabolites")
    return results


############## Internal helpers #############

def _parse_namelist(text, block_name):
    """
    Extract key=value pairs from a FORTRAN namelist block.
    e.g. $SEQPAR ... $END
    Returns a dict of {KEY: value_string}
    """
    pattern = rf'\${block_name}\b(.*?)\$END'
    match = re.search(pattern, text, re.IGNORECASE | re.DOTALL)
    if not match:
        return {}
    return _parse_inline_namelist(match.group(1))


def _parse_inline_namelist(text):
    """Parse key=value pairs from namelist text."""
    result = {}
    # Match KEY= 'value' or KEY= number or KEY= F/T
    for m in re.finditer(r"(\w+)\s*=\s*('.*?'|[^,\n$]+)", text):
        key = m.group(1).strip().upper()
        val = m.group(2).strip().rstrip(',').strip()
        result[key] = val
    return result


def _parse_fid_data(data_text, n_expected):
    """
    Parse real/imag FID pairs from LCModel data text.
    LCModel stores alternating real, imag values.
    Handles both 8E13.5 (row) and 2E15.6 (column) FORTRAN formats.
    """
    # Extract all numbers from the data section
    nums = re.findall(r'[+-]?\d+\.\d+[Ee][+-]?\d+|[+-]?\d+\.\d+', data_text)

    if len(nums) < 2:
        return None

    vals = [float(v) for v in nums]

    # Pair up real/imag
    n_pairs = len(vals) // 2
    fid = np.array(vals[0::2]) + 1j * np.array(vals[1::2])

    # Trim or warn if length doesn't match
    if n_expected > 0 and len(fid) != n_expected:
        if len(fid) > n_expected:
            fid = fid[:n_expected]
        else:
            # Zero-pad if shorter
            pad = np.zeros(n_expected - len(fid), dtype=complex)
            fid = np.concatenate([fid, pad])

    return fid

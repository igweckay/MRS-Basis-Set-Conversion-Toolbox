"""
writers/write_fida.py
Basis Set Converter — FID-A writer

Writes a core struct dict to FID-A .mat format, one file per metabolite.

FID-A output struct (confirmed from sim_readout.m):
    out.fids          - complex FID, column vector (n x 1)
    out.specs         - fftshift(ifft(fids))  ← ifft convention
    out.spectralwidth - bandwidth in Hz
    out.dwelltime     - dwell time in seconds
    out.n             - number of points
    out.linewidth     - linewidth in Hz
    out.Bo            - field strength in Tesla
    out.txfrq         - Larmor frequency in Hz (= sf * 1e6)
    out.nucleus       - '1H'
    out.gamma         - 42577000
    out.t             - time axis (n x 1)
    out.ppm           - ppm axis (n x 1)
    out.sz            - [n, 1]
    out.date          - date string
    out.dims          - struct
    out.flags         - struct
    out.averages      - 1
    out.rawAverages   - 1
    out.subspecs      - 1
    out.rawSubspecs   - 1

Key difference from Osprey:
    FID-A uses fftshift(ifft(fids)) — NOT fft
    Osprey uses fftshift(fft(fids))

Entry points:
    write_fida_file(core, outpath)         → writes one .mat file
    write_fida_folder(basis_list, outdir)  → writes one .mat per metabolite
"""

import os
import numpy as np
import scipy.io as sio
from datetime import datetime

from writers.centre import fid_at_centre


# Proton gyromagnetic ratio
PROTON_GAMMA = 42577000.0   # Hz/T


def write_fida_file(core, outpath):
    """
    Write one core struct dict to a FID-A .mat file.

    Parameters
    ----------
    core : dict
        Core struct with keys: fid, sw, sf, n, name
        Optional: te, seq, Bo, linewidth
    outpath : str
        Full path for the output .mat file.
    """
    outpath = os.path.abspath(outpath)
    os.makedirs(os.path.dirname(outpath), exist_ok=True)

    fid  = fid_at_centre(core).astype(complex)   # its ppm axis below is at 4.65
    sw   = float(core['sw'])
    sf   = float(core['sf'])          # MHz
    n    = int(core.get('n', len(fid)))
    name = str(core.get('name', 'unknown'))

    # Derived fields
    Bo        = core.get('Bo') or (sf / 42.577)
    txfrq     = sf * 1e6               # MHz → Hz
    dwelltime = 1.0 / sw
    stated    = core.get('linewidth')          # stated as None when the source has none
    linewidth = 1.0 if stated is None else float(stated)
    te        = float(core.get('te', 0.0)) if core.get('te') else 0.0
    seq       = str(core.get('seq', 'unedited'))
    date_str  = datetime.now().strftime('%d-%b-%Y')

    # FID-A convention: specs = fftshift(ifft(fids))
    specs = np.fft.fftshift(np.fft.ifft(fid))

    # Time axis
    t = np.arange(n) * dwelltime

    # PPM axis (FID-A formula from sim_readout.m)
    freq = np.linspace(
        -sw / 2 + sw / (2 * n),
        sw / 2 - sw / (2 * n),
        n
    )
    ppm = -freq / (Bo * 42.577) + 4.65

    # dims struct
    dims = {
        't'        : 1.0,
        'coils'    : 0.0,
        'averages' : 0.0,
        'subSpecs' : 0.0,
        'extras'   : 0.0,
    }

    # flags struct
    flags = {
        'writtentostruct' : 1.0,
        'gotparams'       : 1.0,
        'leftshifted'     : 0.0,
        'filtered'        : 0.0,
        'zeropadded'      : 0.0,
        'freqcorrected'   : 0.0,
        'phasecorrected'  : 0.0,
        'averaged'        : 1.0,
        'addedrcvrs'      : 1.0,
        'subtracted'      : 1.0,
        'writtentotext'   : 0.0,
        'downsampled'     : 0.0,
        'isFourSteps'     : 0.0,
    }

    # Build the 'out' struct
    out = {
        'fids'          : fid.reshape(-1, 1),
        'specs'         : specs.reshape(-1, 1),
        'spectralwidth' : sw,
        'dwelltime'     : dwelltime,
        'n'             : float(n),
        'linewidth'     : linewidth,
        'Bo'            : Bo,
        'txfrq'         : txfrq,
        'nucleus'       : '1H',
        'gamma'         : PROTON_GAMMA,
        't'             : t.reshape(-1, 1),
        'ppm'           : ppm.reshape(-1, 1),
        'sz'            : np.array([n, 1], dtype=float),    # double, as MATLAB's size() is: integer
                                                            # arithmetic on it rounds FID-A's axes
        'te'            : te,
        'seq'           : seq,
        'sim'           : 'converted',
        'date'          : date_str,
        'dims'          : dims,
        'flags'         : flags,
        'averages'      : 1.0,
        'rawAverages'   : 1.0,
        'subspecs'      : 1.0,
        'rawSubspecs'   : 1.0,
    }

    sio.savemat(outpath, {'out': out})


def write_fida_folder(basis_list, outdir):
    """
    Write a list of core struct dicts to FID-A .mat files,
    one file per metabolite.

    Parameters
    ----------
    basis_list : list of dicts
    outdir     : str, output directory
    """
    os.makedirs(outdir, exist_ok=True)
    written = []
    failed  = []

    for core in basis_list:
        name    = core.get('name', 'unknown')
        outpath = os.path.join(outdir, f"{name}.mat")
        try:
            write_fida_file(core, outpath)
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
    print(f"\n  Done. {len(written)} FID-A files written to: {outdir}")
    return written


######################## Command line usage #######################

if __name__ == '__main__':
    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    if len(sys.argv) < 3:
        print("Usage: python write_fida.py <input_folder_or_file> <output_folder>")
        print()
        print("Example:")
        print("  python write_fida.py /path/to/marss_mat_folder /path/to/fida_output")
        sys.exit(1)

    input_path = sys.argv[1]
    outdir     = sys.argv[2]

    from readers.read_marss import read_marss_file, read_marss_folder

    if os.path.isdir(input_path):
        print(f"Reading folder: {input_path}")
        basis_list = read_marss_folder(input_path)
    elif os.path.isfile(input_path):
        print(f"Reading file: {input_path}")
        basis_list = [read_marss_file(input_path)]
    else:
        print(f"Error: {input_path} is not a valid file or folder")
        sys.exit(1)

    print(f"Loaded {len(basis_list)} metabolites")
    print()
    write_fida_folder(basis_list, outdir)

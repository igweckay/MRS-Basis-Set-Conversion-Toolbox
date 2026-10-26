"""
writers/write_mrscloud.py
Basis Set Converter — MRSCloud writer

Writes core struct dicts to MRSCloud .mat format.
MRSCloud stores one metabolite per .mat file.

MRSCloud .mat structure (confirmed from your documentation):
    fids          - complex FID (n x 1)
    specs         - complex spectrum (n x 1)  fftshift(ifft(fids))
    spectralwidth - bandwidth in Hz
    dwelltime     - dwell time in seconds
    n             - number of points
    linewidth     - 1
    Bo            - field strength in Tesla
    txfrq         - Larmor frequency in Hz
    t             - time axis (1 x n)
    ppm           - ppm axis (1 x n)
    sz            - [n, 1]
    name          - metabolite name string
    dims          - struct
    flags         - struct
    averages      - 1
    rawAverages   - 1
    subspecs      - 1
    rawSubspecs   - 1
    date          - date string

Entry points:
    write_mrscloud_file(core, outpath)        -> one .mat file
    write_mrscloud_folder(basis_list, outdir) -> folder of .mat files
"""

import os
import numpy as np
import scipy.io as sio
from datetime import datetime

from writers.centre import fid_at_centre


PROTON_GAMMA = 42577000.0   # Hz/T


def write_mrscloud_file(core, outpath):
    """
    Write one core struct dict to a MRSCloud .mat file.

    Parameters
    ----------
    core    : dict with keys: fid, sw, sf, n, name
    outpath : str, full output path
    """
    outpath = os.path.abspath(outpath)
    os.makedirs(os.path.dirname(outpath), exist_ok=True)

    fid   = fid_at_centre(core).astype(complex)   # its ppm axis below is at 4.65
    sw    = float(core['sw'])
    sf    = float(core['sf'])        # MHz
    n     = int(core.get('n', len(fid)))
    name  = str(core.get('name', 'unknown'))
    Bo    = core.get('Bo') or (sf / 42.577)
    txfrq = sf * 1e6                 # MHz → Hz

    dt    = 1.0 / sw
    # MRSCloud uses ifft convention for specs
    specs = np.fft.fftshift(np.fft.ifft(fid))

    t   = np.arange(n) * dt
    # ppm axis as MRSCloud (FID-A's sim_readout.m) stores it: descending, symmetric about 0 Hz
    freq= np.linspace(-sw / 2 + sw / (2 * n), sw / 2 - sw / (2 * n), n)
    ppm = -freq / (Bo * 42.577) + 4.65

    dims = {
        't'       : 1.0,
        'coils'   : 0.0,
        'averages': 0.0,
        'subSpecs': 0.0,
        'extras'  : 0.0,
    }

    flags = {
        'writtentostruct': 1.0,
        'gotparams'      : 1.0,
        'leftshifted'    : 0.0,
        'filtered'       : 0.0,
        'zeropadded'     : 0.0,
        'freqcorrected'  : 0.0,
        'phasecorrected' : 0.0,
        'averaged'       : 1.0,
        'addedrcvrs'     : 1.0,
        'subtracted'     : 1.0,
        'writtentotext'  : 0.0,
        'downsampled'    : 0.0,
        'isFourSteps'    : 0.0,
    }

    struct = {
        'fids'         : fid.reshape(-1, 1),
        'specs'        : specs.reshape(-1, 1),
        'spectralwidth': sw,
        'dwelltime'    : dt,
        'n'            : float(n),
        'linewidth'    : 1.0 if core.get('linewidth') is None else float(core['linewidth']),
        'Bo'           : Bo,
        'txfrq'        : txfrq,
        't'            : t.reshape(1, -1),
        'ppm'          : ppm.reshape(1, -1),
        'sz'           : np.array([n, 1]),
        'name'         : name,
        'dims'         : dims,
        'flags'        : flags,
        'averages'     : 1.0,
        'rawAverages'  : 1.0,
        'subspecs'     : 1.0,
        'rawSubspecs'  : 1.0,
        'date'         : datetime.now().strftime('%d-%b-%Y'),
    }

    sio.savemat(outpath, struct)


def write_mrscloud_folder(basis_list, outdir):
    """
    Write all metabolites to individual MRSCloud .mat files in a folder.
    """
    os.makedirs(outdir, exist_ok=True)
    written = []
    failed  = []

    for core in basis_list:
        name    = core.get('name', 'unknown')
        outpath = os.path.join(outdir, f"{name}.mat")
        try:
            write_mrscloud_file(core, outpath)
            written.append(outpath)
            print(f"  Written: {name}.mat")
        except Exception as ex:
            failed.append(name)
            print(f"  WARNING: failed for {name} — {ex}")

    if failed:
        print(f"\n  {len(failed)} failed: {', '.join(failed)}")
    print(f"\n  Done. {len(written)} MRSCloud .mat files written to: {outdir}")
    return written


################# Command line usage ###################

if __name__ == '__main__':
    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    if len(sys.argv) < 3:
        print("Usage: python write_mrscloud.py <input_folder> <output_folder>")
        sys.exit(1)

    input_path = sys.argv[1]
    outdir     = sys.argv[2]

    from readers.read_marss import read_marss_file, read_marss_folder

    if os.path.isdir(input_path):
        basis_list = read_marss_folder(input_path)
    else:
        basis_list = [read_marss_file(input_path)]

    print(f"Loaded {len(basis_list)} metabolites\n")
    write_mrscloud_folder(basis_list, outdir)

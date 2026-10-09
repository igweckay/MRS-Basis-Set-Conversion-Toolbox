"""
writers/write_fsLmrs.py
Basis Set Converter — FSL-MRS writer

Writes a list of core struct dicts to FSL-MRS JSON basis format.
One .json file per metabolite in an output folder.

FSL-MRS .json structure (confirmed from your documentation):
    {
        "basis": {
            "basis_re"    : [real FID values],
            "basis_im"    : [imag FID values],
            "basis_dwell" : dwell time in seconds,
            "basis_centre": the spectrometer frequency in MHz (FSL-MRS reads it as
                            centralFrequency = basis_centre * 1E6),
            "basis_width" : linewidth in Hz, or null,
            "basis_name"  : "metabolite name"
        },
        "meta": {
            "time"      : "YYYYMMDD_HHMMSS",
            "SimVersion": "Converted by basis_set_converter"
        }
    }

Note:
    - basis_centre is sf in MHz; FSL-MRS stores no ppm reference (its nucleus' default)
    - basis_dwell = 1/sw in seconds

Entry points:
    write_fsLmrs_file(core, outpath)        → one .json file
    write_fsLmrs_folder(basis_list, outdir) → folder of .json files
"""

import os
import json
import numpy as np
from datetime import datetime


def write_fsLmrs_file(core, outpath):
    """
    Write one core struct dict to a FSL-MRS .json basis file.

    Parameters
    ----------
    core    : dict with keys: fid, sw, sf, n, name
    outpath : str, full output path including filename
    """
    outpath = os.path.abspath(outpath)
    os.makedirs(os.path.dirname(outpath), exist_ok=True)

    fid   = np.asarray(core['fid']).ravel().astype(complex)
    sw    = float(core['sw'])
    name  = str(core.get('name', 'unknown'))
    dwell = 1.0 / sw

    sf = core.get('sf')
    if sf is None:
        raise RuntimeError(f"{name}: FSL-MRS stores the spectrometer frequency (sf), which "
                           "the basis does not have")
    linewidth = core.get('linewidth')

    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

    data = {
        "basis": {
            "basis_re"    : fid.real.tolist(),
            "basis_im"    : fid.imag.tolist(),
            "basis_dwell" : dwell,
            "basis_centre": float(sf),
            "basis_width" : None if linewidth is None else float(linewidth),
            "basis_name"  : name,
        },
        "meta": {
            "time"       : timestamp,
            "SimVersion" : "Converted by basis_set_converter",
        }
    }

    with open(outpath, 'w') as f:
        json.dump(data, f, indent=4)


def write_fsLmrs_folder(basis_list, outdir):
    """
    Write all metabolites to individual FSL-MRS .json files in a folder.

    Parameters
    ----------
    basis_list : list of core struct dicts
    outdir     : str

    Returns
    -------
    list of str — paths to written .json files
    """
    os.makedirs(outdir, exist_ok=True)
    written = []
    failed  = []

    for core in basis_list:
        name    = core.get('name', 'unknown')
        outpath = os.path.join(outdir, f"{name}.json")
        try:
            write_fsLmrs_file(core, outpath)
            written.append(outpath)
            print(f"  Written: {name}.json")
        except Exception as ex:
            failed.append(name)
            print(f"  WARNING: failed for {name} — {ex}")

    if failed:
        print(f"\n  {len(failed)} failed: {', '.join(failed)}")
    print(f"\n  Done. {len(written)} FSL-MRS .json files written to: {outdir}")
    return written


######################## Command line usage #######################

if __name__ == '__main__':
    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    if len(sys.argv) < 3:
        print("Usage: python write_fsLmrs.py <input_folder> <output_folder>")
        sys.exit(1)

    input_path = sys.argv[1]
    outdir     = sys.argv[2]

    from readers.read_marss import read_marss_file, read_marss_folder

    if os.path.isdir(input_path):
        basis_list = read_marss_folder(input_path)
    else:
        basis_list = [read_marss_file(input_path)]

    print(f"Loaded {len(basis_list)} metabolites\n")
    write_fsLmrs_folder(basis_list, outdir)

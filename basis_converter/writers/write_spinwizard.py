"""
writers/write_spinwizard.py
Basis Set Converter — SpinWizard / JET writer

Writes basis functions in SpinWizard format for use with the JET toolkit.

Format (confirmed):
    - No file extension — filename = metabolite name (e.g. 'NAA', 'Cho')
    - Two whitespace-separated columns: real and imaginary FID values
    - No header
    - One file per metabolite in a Basisset folder
    - Tab-separated, 10 decimal places (matches SpinWizard output style)

JET GitHub: https://github.com/SAIL-GuoLab/JET

Entry point:
    write_spinwizard(basis_list, outdir) → list of output file paths
"""

import os
import numpy as np

from writers.centre import CENTRE, fid_at_centre


def write_spinwizard(basis_list, outdir):
    """
    Write all metabolites as SpinWizard basis files.

    Parameters
    ----------
    basis_list : list of core struct dicts
    outdir     : str, output folder path (files written directly here)

    Returns
    -------
    list of str — paths to written files
    """
    outdir = os.path.abspath(outdir)
    os.makedirs(outdir, exist_ok=True)

    written = []
    names   = []

    for core in basis_list:
        name = core.get('name', 'unknown')
        fid  = fid_at_centre(core)   # at 4.65, the centre BasisRFOffset states

        # Output path — no extension, just the metabolite name
        outpath = os.path.join(outdir, name)

        # Build two-column array
        data = np.column_stack([fid.real, fid.imag])

        # Write — tab separated, matching SpinWizard style
        with open(outpath, 'w') as f:
            for row in data:
                f.write(f"  {row[0]:20.10e}\t\t  {row[1]:20.10e}\n")

        written.append(outpath)
        names.append(name)
        print(f"  Written: {name}  ({len(fid)} points)")

    # Write LCMBasisSet.txt — JET metabolite list file
    # Format: MetaboliteName    Y    0    1
    lcm_path = os.path.join(outdir, 'LCMBasisSet.txt')
    with open(lcm_path, 'w') as f:
        for name in names:
            f.write(f"{name:<20s}Y     0     1\n")

    written.append(lcm_path)
    print(f"  Written: LCMBasisSet.txt  ({len(names)} metabolites)")

    # Write BasisSetParameters.txt — sf (MHz) and sw (kHz)
    if basis_list[0].get('sf') and basis_list[0].get('sw'):
        sf = float(basis_list[0]['sf'])
        sw = float(basis_list[0]['sw']) / 1000.0  # Hz → kHz
        params_path = os.path.join(outdir, 'BasisSetParameters.txt')
        with open(params_path, 'w') as f:
            f.write(f"  {sf:.5f}\n  {sw:.5f}\n")
        written.append(params_path)
        print(f"  Written: BasisSetParameters.txt  (sf={sf:.2f} MHz, sw={sw:.3f} kHz)")

    # Write BasisRFOffset.txt — ppmCalib
    # (4.65, where the functions are written: fid_at_centre; a ppmCalib read is a label, MARSS's
    # 4.675 sets are at 4.65)
    ppm_calib = CENTRE
    offset_path = os.path.join(outdir, 'BasisRFOffset.txt')
    with open(offset_path, 'w') as f:
        f.write(f"  {ppm_calib:.4f}\n")
    written.append(offset_path)
    print(f"  Written: BasisRFOffset.txt  (ppmCalib={ppm_calib} ppm)")

    print(f"\n  Done. {len(names)} SpinWizard files + companion files written to: {outdir}")
    return written


# ── Command line ──────────────────────────────────────────────────────────────

if __name__ == '__main__':
    import sys, os
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    if len(sys.argv) < 3:
        print("Usage: python write_spinwizard.py <input_folder> <output_folder>")
        print()
        print("Example:")
        print("  python write_spinwizard.py /path/to/marss_folder /path/to/Basisset")
        sys.exit(1)

    input_path = sys.argv[1]
    outdir     = sys.argv[2]

    from readers.read_marss import read_marss_file, read_marss_folder

    if os.path.isdir(input_path):
        basis_list = read_marss_folder(input_path)
    else:
        basis_list = [read_marss_file(input_path)]

    print(f"Loaded {len(basis_list)} metabolites\n")
    write_spinwizard(basis_list, outdir)

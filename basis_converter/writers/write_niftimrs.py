"""
writers/write_niftimrs.py
Basis Set Converter — NIfTI-MRS writer

Converts a list of core struct dicts (one per metabolite) to NIfTI-MRS
format (.nii.gz) using the confirmed pipeline:

    core struct -> .raw (via write_lcmodel_raw) -> spec2nii -> .nii.gz -> conjugated

NIfTI-MRS FIDs are the complex conjugate of a basis' own (FSL-MRS reads NAA at 2.01 ppm
only then), as the conj step of step_2_batch_convert_raw_to_nifti.sh made them; the
conjugation and the echo time, which spec2nii's raw route does not store, are applied to
spec2nii's file in place. (fsl_mrs_proc conj --output takes a folder: it wrote its file
into a folder named <name>_conj.nii.gz, or failed, and the unconjugated file was used.)

Pipeline confirmed from:
    step_1_mat_to_raw.py    — MARSS .mat -> .raw
    step_2_batch_convert_raw_to_nifti.sh — .raw -> .nii.gz -> conjugated

Requirements:
    spec2nii    — pip install spec2nii
    nibabel     — pip install nibabel

spec2nii parameters (hardcoded for now, made dynamic in screen 2):
    -n 1H           nucleus
    -i {sf}         Larmor frequency in MHz
    -b {sw}         bandwidth in Hz
    -j              write JSON sidecar
    -f {name}       output filename stem

Entry point:
    write_niftimrs(basis_list, outdir, te=None) -> list of output .nii.gz paths
"""

import os
import subprocess
import tempfile
import numpy as np

from writers.centre import fid_at_centre


def write_niftimrs(basis_list, outdir, te=None):
    """
    Convert a list of core struct dicts to NIfTI-MRS .nii.gz files.

    Parameters
    ----------
    basis_list : list of dicts
        Each dict is a core struct with keys: fid, sw, sf, n, name, source
    outdir : str
        Output directory for .nii.gz files.
    te : float or None
        Echo time in ms, for basis functions that state none.

    Returns
    -------
    list of str — paths to output .nii.gz files
    """
    os.makedirs(outdir, exist_ok=True)

    # Use a temp folder for intermediate .raw files
    with tempfile.TemporaryDirectory() as tmpdir:
        output_paths = []
        failed       = []

        for core in basis_list:
            name = core['name']
            try:
                ############### Step 1: write .raw ##############
                raw_path = os.path.join(tmpdir, f"{name}.raw")
                _write_raw(core, raw_path)

                ############### Step 2: spec2nii raw -> .nii.gz ##############
                nii_path = _run_spec2nii(raw_path, name, outdir, core)

                ############### Step 3: conjugate, add the echo time ##############
                _finish(nii_path, core.get('te') or te)

                output_paths.append(nii_path)
                print(f"  Written: {os.path.basename(nii_path)}")

            except Exception as ex:
                failed.append(name)
                print(f"  WARNING: failed for {name} — {ex}")

        if failed:
            print(f"\n  {len(failed)} metabolite(s) failed: {', '.join(failed)}")

        print(f"\n  Done. {len(output_paths)} NIfTI-MRS files written to: {outdir}")
        return output_paths


############### Step 1: write .raw ##############

def _write_raw(core, raw_path):
    """
    Write a single core struct to LCModel .raw format.
    Uses the confirmed convention from step_1_mat_to_raw.py:
        - imaginary sign flipped: vals = [real, -imag, ...]
        - HZPPPM = sf (already in MHz — NOT divided by 1e6)
        - NUNFIL = n
        - DELTAT = 1/sw
    """
    fid   = fid_at_centre(core)    # NIfTI-MRS states no ppm centre: read at 4.65
    sf    = float(core['sf'])      # MHz — already correct
    sw    = float(core['sw'])      # Hz
    dwell = 1.0 / sw

    # Interleave real / -imag pairs (LCModel convention)
    vals = []
    for z in fid:
        vals.extend([z.real, -z.imag])

    with open(raw_path, 'w') as f:
        f.write(" $SEQPAR\n")
        f.write(f" HZPPPM= {sf:.6f}\n")
        f.write(f" NUNFIL= {len(fid)}\n")
        f.write(f" DELTAT= {dwell:.12f}\n")
        f.write(" $END\n")
        f.write(" $NMID\n")
        f.write(f" ID='{core['name']}'\n")
        f.write(" FMTDAT='(8E13.5)'\n")
        f.write(" VOLUME=   1.00000E+00\n")
        f.write(" TRAMP=    1.00000E+00\n")
        f.write(" $END\n")

        for i in range(0, len(vals), 8):
            chunk = vals[i:i + 8]
            f.write(" ".join(f"{v:+.5E}" for v in chunk) + "\n")


############### Step 2: spec2nii ##############

def _run_spec2nii(raw_path, name, outdir, core):
    """
    Run spec2nii to convert .raw -> .nii.gz

    spec2nii raw parameters:
        -n 1H           proton nucleus
        -i {sf}         Larmor frequency in MHz
        -b {sw}         bandwidth in Hz
        -f {name}       output filename stem
        -j              write JSON sidecar
        -o {outdir}     output directory
    """
    sf = float(core['sf'])
    sw = float(core['sw'])

    cmd = [
        'spec2nii', 'raw',
        '-n', '1H',
        '-i', f"{sf:.6f}",
        '-b', f"{sw:.1f}",
        '-f', name,
        '-j',
        '-o', outdir,
        raw_path,
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(
            f"spec2nii failed for {name}:\n{result.stderr}"
        )

    nii_path = os.path.join(outdir, f"{name}.nii.gz")
    if not os.path.isfile(nii_path):
        raise RuntimeError(
            f"spec2nii ran but output not found: {nii_path}\n"
            f"Files in outdir: {os.listdir(outdir)}"
        )

    return nii_path


############### Step 3: conjugate ##############

def _finish(nii_path, te=None):
    """
    Conjugate the FIDs of spec2nii's file, as NIfTI-MRS holds them, and store the echo
    time (ms) as EchoTime (s) in its header extension.
    """
    import json
    import nibabel as nib

    img  = nib.load(nii_path)
    ext  = next(e for e in img.header.extensions if e.get_code() == 44)
    meta = json.loads(ext.get_content())
    if te:
        meta['EchoTime'] = float(te) / 1000.0
    out  = type(img)(np.conj(np.asarray(img.dataobj)), img.affine, header=img.header)
    out.header.extensions.remove(next(e for e in out.header.extensions if e.get_code() == 44))
    out.header.extensions.append(nib.nifti1.Nifti1Extension(44, json.dumps(meta).encode()))
    nib.save(out, nii_path)


############### Dependency check ##############

def check_dependencies():
    """
    Check that spec2nii is available.
    Returns list of missing tools.
    """
    missing = []
    for tool in ['spec2nii']:
        result = subprocess.run(
            ['which', tool],
            capture_output=True, text=True
        )
        if result.returncode != 0:
            missing.append(tool)
    return missing


############### Command line usage ##############

if __name__ == '__main__':
    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    if len(sys.argv) < 3:
        print("Usage: python write_niftimrs.py <input_folder_or_file> <output_folder>")
        print()
        print("Examples:")
        print("  # Convert a folder of MARSS .mat files:")
        print("  python write_niftimrs.py /path/to/basis_mat_folder /path/to/output")
        print()
        print("  # Convert a single MARSS .mat file:")
        print("  python write_niftimrs.py /path/to/NAA.mat /path/to/output")
        sys.exit(1)

    input_path = sys.argv[1]
    outdir     = sys.argv[2]

    # Check dependencies first
    missing = check_dependencies()
    if missing:
        print(f"Missing tools: {', '.join(missing)}")
        print("Install with:  pip install spec2nii")
        sys.exit(1)

    # Load basis functions
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
    print(f"Output folder: {outdir}")
    print()

    results = write_niftimrs(basis_list, outdir)

    print(f"\nDone. {len(results)} files written.")
    for r in results:
        print(f"  {r}")

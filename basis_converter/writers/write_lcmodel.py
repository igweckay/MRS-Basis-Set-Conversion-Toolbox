"""
writers/write_lcmodel.py
Basis Set Converter — LCModel writer

Writes basis set to LCModel format:
    1. Individual .raw files (one per metabolite)
    2. Combined .BASIS file (all metabolites in one file)
    3. Both (recommended)

.raw file format (confirmed from step_1_mat_to_raw.py and your documentation):
    $SEQPAR
     HZPPPM= {sf}        <- Larmor frequency in MHz (sf, already MHz)
     NUNFIL= {n}         <- number of complex points
     DELTAT= {1/sw}      <- dwell time in seconds
    $END
    $NMID
     ID='{name}'
     FMTDAT='(8E13.5)'
     VOLUME=   1.00000E+00
     TRAMP=    1.00000E+00
    $END
    {real} {-imag} ...   <-  8 values per line, imaginary sign flipped

.BASIS file format (confirmed from your SLASER_Siemens26.basis):
    $SEQPAR
     FWHMBA=  {linewidth}
     HZPPPM=  {sf}
     ECHOT=   {te}
     SEQ='{seq}'
    $END
    $BASIS1
     IDBASI='{description}'
     FMTBAS='(6E13.5)'
     BADELT=  {1/sw}
     NDATAB=  {2n}
    $END
    For each metabolite:
        $NMUSED ... $END   ← fitting defaults
        $BASIS ... $END    ← metabolite metadata
        {spectrum data}    ← FFT of the FID zero-filled to 2n

Key convention (confirmed from step_1_mat_to_raw.py):
    - Imaginary sign is FLIPPED: write [real, -imag]
    - HZPPPM = sf in MHz (do NOT divide by 1e6 — already MHz)

Entry points:
    write_lcmodel_raw(core, outpath)            → one .raw file
    write_lcmodel_raw_folder(basis_list, outdir)→ folder of .raw files
    write_lcmodel_basis(basis_list, outpath,
                        te, seq, description)   → one .BASIS file
    write_lcmodel(basis_list, outdir, mode,
                  te, seq, description)         → main entry point
                  mode: 'raw' | 'basis' | 'both'
"""

import os
import numpy as np
from datetime import datetime


####################### Main entry point #######################

def write_lcmodel(basis_list, outdir, mode='both',
                  te=None, seq='unedited', description=None):
    """
    Write a basis set to LCModel format.

    Parameters
    ----------
    basis_list : list of core struct dicts
    outdir     : str, output directory
    mode       : 'raw' | 'basis' | 'both'
    te         : float or None, echo time in ms
    seq        : str, sequence name
    description: str or None, basis set description for .BASIS header

    Returns
    -------
    dict with keys:
        'raw_files'  : list of .raw file paths (if mode includes raw)
        'basis_file' : path to .BASIS file (if mode includes basis)
    """
    os.makedirs(outdir, exist_ok=True)
    result = {}

    if mode in ('raw', 'both'):
        raw_files = write_lcmodel_raw_folder(basis_list, outdir)
        result['raw_files'] = raw_files

    if mode in ('basis', 'both'):
        basis_path = os.path.join(outdir, 'basis_set.BASIS')
        write_lcmodel_basis(
            basis_list, basis_path,
            te=te, seq=seq, description=description
        )
        result['basis_file'] = basis_path

    return result


######################## .raw writer #######################

def write_lcmodel_raw(core, outpath):
    """
    Write one core struct dict to a LCModel .raw file.

    Parameters
    ----------
    core    : dict with keys: fid, sw, sf, n, name
    outpath : str, full output path including filename
    """
    outpath = os.path.abspath(outpath)
    os.makedirs(os.path.dirname(outpath), exist_ok=True)

    fid  = np.asarray(core['fid']).ravel()
    sf   = float(core['sf'])          # MHz — already correct, do NOT divide
    sw   = float(core['sw'])          # Hz
    n    = len(fid)
    name = str(core.get('name', 'unknown'))
    dwell = 1.0 / sw

    # Interleave real / -imag (LCModel convention, confirmed from step_1_mat_to_raw.py)
    vals = []
    for z in fid:
        vals.extend([z.real, -z.imag])

    with open(outpath, 'w') as f:
        # $SEQPAR block
        f.write(" $SEQPAR\n")
        f.write(f" HZPPPM= {sf:.6f}\n")
        f.write(f" NUNFIL= {n}\n")
        f.write(f" DELTAT= {dwell:.12f}\n")
        f.write(" $END\n")

        # $NMID block
        f.write(" $NMID\n")
        f.write(f" ID='{name}'\n")
        f.write(" FMTDAT='(8E13.5)'\n")
        f.write(" VOLUME=   1.00000E+00\n")
        f.write(" TRAMP=    1.00000E+00\n")
        f.write(" $END\n")

        # FID data — 8 values per line
        for i in range(0, len(vals), 8):
            chunk = vals[i:i + 8]
            f.write(" ".join(f"{v:+.5E}" for v in chunk) + "\n")


def write_lcmodel_raw_folder(basis_list, outdir):
    """
    Write all metabolites to individual .raw files in a folder.

    Parameters
    ----------
    basis_list : list of core struct dicts
    outdir     : str

    Returns
    -------
    list of str — paths to written .raw files
    """
    os.makedirs(outdir, exist_ok=True)
    written = []
    failed  = []

    for core in basis_list:
        name    = core.get('name', 'unknown')
        outpath = os.path.join(outdir, f"{name}.raw")
        try:
            write_lcmodel_raw(core, outpath)
            written.append(outpath)
            print(f"  Written: {name}.raw")
        except Exception as ex:
            failed.append(name)
            print(f"  WARNING: failed for {name} — {ex}")

    if failed:
        print(f"\n  {len(failed)} failed: {', '.join(failed)}")
    print(f"\n  Done. {len(written)} .raw files written to: {outdir}")
    return written


######################## .BASIS writer #######################

def write_lcmodel_basis(basis_list, outpath,
                        te=None, seq='unedited', description=None):
    """
    Write all metabolites to a single LCModel .BASIS file.

    Parameters
    ----------
    basis_list  : list of core struct dicts
    outpath     : str, full path to output .BASIS file
    te          : float or None, echo time in ms
    seq         : str, sequence name
    description : str or None, description for IDBASI field
    """
    outpath = os.path.abspath(outpath)
    os.makedirs(os.path.dirname(outpath), exist_ok=True)

    if not basis_list:
        raise RuntimeError("basis_list is empty — nothing to write")

    # Get shared parameters from first metabolite
    first = basis_list[0]
    sf    = float(first['sf'])
    sw    = float(first['sw'])
    n     = int(first.get('n', len(first['fid'])))
    dwell = 1.0 / sw

    linewidth   = float(first.get('linewidth', 0.01))
    te_val      = float(te) if te is not None else 0.0
    date_str    = datetime.now().strftime('%d-%b-%Y')
    description = description or f"Basis set converted on {date_str}"

    with open(outpath, 'w') as f:

        ######################## $SEQPAR block #######################
        f.write(" $SEQPAR\n")
        f.write(f" FWHMBA=  {linewidth:.8E},\n")
        f.write(f" HZPPPM=  {sf:.6f}    ,\n")
        f.write(f" ECHOT=  {te_val:.7f}    ,\n")
        f.write(f" SEQ='{seq[:6]}',\n")   # LCModel truncates to 6 chars
        f.write(" \n")
        f.write(" $END\n")

        ######################## $BASIS1 block #######################
        f.write(" $BASIS1\n")
        f.write(f" IDBASI='{description[:80]}',\n")
        f.write(" FMTBAS='(6E13.5)                                                                        ',\n")
        f.write(f" BADELT=  {dwell:.8E},\n")
        f.write(f" NDATAB=       {2 * n},\n")   # zero-filled to 2n, as MakeBasis
        f.write(" \n")
        f.write(" $END\n")

        ######################## Per-metabolite sections #######################
        for core in basis_list:
            name = core.get('name', 'unknown')
            fid  = np.asarray(core['fid']).ravel()

            # $NMUSED block (fitting defaults — standard values)
            f.write(" $NMUSED\n")
            f.write(f" FILRAW='{name}.raw',\n")
            f.write(" METABO_CONTAM='      ',\n")
            f.write(" METABO_SINGLET='      ',\n")
            f.write(" AUTOPH=F,\n")
            f.write(" AUTOSC=F,\n")
            f.write(" CONSISTENT_SCALING=T,\n")
            f.write(" DO_CONTAM=F,\n")
            f.write(" FLATEN=F,\n")
            f.write(" NEGBAS=F,\n")
            f.write(" NOSHIF=T,\n")
            f.write(" SCALE1=F,\n")
            f.write(" CONCSC= -1.00000000    ,\n")
            f.write(" DEGZER=  0.00000000    ,\n")
            f.write(" DEGPAP=  0.00000000    ,\n")
            f.write(" DEGPPM=  0.00000000    ,\n")
            f.write(" DKNDEF=  1.00000000    ,\n")
            f.write(" DKNOT=  1.00000000    ,\n")
            f.write(" FOFFSET2=  0.00000000    ,\n")
            f.write(" FWHMSM=  0.00000000    ,\n")
            f.write(" HWDPHA= 0.150000006    ,\n")
            f.write(" HWDSCA= 0.150000006    ,\n")
            f.write(" PPMAPP=  0.00000000    ,-0.400000006    ,\n")
            f.write(" PPMAPP_CONTAM= 2*0.00000000      ,\n")
            f.write(" PPMBAS= 0.100000001    ,\n")
            f.write(" PPMFLA= -1.00000000    ,\n")
            f.write(" PPMGAP= 20*997.000000      ,\n")
            f.write(" PPMOFF= 2*-999.000000     ,\n")
            f.write(" PPMPK=  0.00000000    ,\n")
            f.write(" PPMPK_CONTAM=  0.00000000    ,\n")
            f.write(" PPMPHA=  0.00000000    ,\n")
            f.write(" PPMSCA=  8.43999958    ,\n")
            f.write(" PPMSEP=  4.65000010    ,\n")
            f.write(" PPMSPL= 40*997.000000      ,\n")
            f.write(" PPM_SPLIT= -999.000000    ,\n")
            f.write(" RINTEG=  0.00000000    ,\n")
            f.write(" SDPNTS=  1.00000000    ,\n")
            f.write(" XTRASH=  0.00000000    ,\n")
            f.write(" \n")
            f.write(" $END\n")

            # $BASIS block
            f.write(" $BASIS\n")
            f.write(f" ID='{name}',\n")
            f.write(f" METABO='{name[:6]:<6}',\n")   # LCModel truncates to 6
            f.write(" CONC=  1.00000000    ,\n")
            f.write(" TRAMP=  1.00000000    ,\n")
            f.write(" VOLUME=  1.00000000    ,\n")
            f.write(" ISHIFT=          0,\n")
            f.write(" \n")
            f.write(" $END\n")

            # Spectrum data — LCModel reads this block as the frequency-domain spectrum:
            # the orthonormal FFT of the FID in .raw orientation (imaginary sign flipped),
            # zero-filled to 2n, as MakeBasis writes it. Interleaved real / imag, 6 per line.
            spec = np.fft.fft(np.conj(fid), n=2 * len(fid), norm='ortho')
            vals = []
            for z in spec:
                vals.extend([z.real, z.imag])

            for i in range(0, len(vals), 6):
                chunk = vals[i:i + 6]
                # fixed-width fields, no separator: FMTBAS (6E13.5) reads 13 characters each
                f.write("".join(f"{v:13.5E}" for v in chunk) + "\n")

    print(f"  Written: {os.path.basename(outpath)} "
          f"({len(basis_list)} metabolites)")


######################## Command line usage #######################

if __name__ == '__main__':
    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    if len(sys.argv) < 3:
        print("Usage: python write_lcmodel.py <input_folder> <output_folder> [mode]")
        print()
        print("mode: raw | basis | both (default: both)")
        print()
        print("Example:")
        print("  python write_lcmodel.py /path/to/marss_folder /path/to/output both")
        sys.exit(1)

    input_path = sys.argv[1]
    outdir     = sys.argv[2]
    mode       = sys.argv[3] if len(sys.argv) > 3 else 'both'

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
    print(f"Mode: {mode}")
    print()

    result = write_lcmodel(basis_list, outdir, mode=mode, te=26.0, seq='sLASER')

    print()
    if 'raw_files' in result:
        print(f"  .raw files : {len(result['raw_files'])} written")
    if 'basis_file' in result:
        print(f"  .BASIS file: {result['basis_file']}")

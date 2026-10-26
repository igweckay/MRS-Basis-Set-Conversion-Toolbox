"""
write_jmrui.py
Basis Set Converter — jMRUI / QUEST writer

Writes a single basis function FID to jMRUI ASCII .txt format.
One file per metabolite; place all files in a folder for use as a basis set.

Format confirmed from FID-A io_writejmrui.m (CIC-methods/FID-A, GitHub).

Header fields:
    jMRUI Data Textfile
    Filename            : output file path
    PointsInDataset     : n (number of FID points)
    DatasetsInFile      : 1
    SamplingInterval    : dwell time in ms  (= 1000 / SW)
    ZeroOrderPhase      : 0
    BeginTime           : 0
    TransmitterFrequency: Larmor frequency in Hz (= Bo * gamma)
    MagneticField       : Bo in Tesla
    TypeOfNucleus       : 0 (proton)
    NameOfPatient       : metabolite name (used as identifier)
    DateOfExperiment    : YYYYMMDD
    Spectrometer        : scanner name or 'Unknown'
    AdditionalInfo      : free text

Data section (4 columns, tab-separated):
    sig(real)   sig(imag)   fft(real)   fft(imag)

Notes:
    - Water-referenced: ppm axis centred at 4.65 ppm
    - FFT convention: specs = fftshift(ifft(fids))  [same as FID-A]
    - SamplingInterval is in milliseconds (dwelltime * 1000)
    - TransmitterFrequency is in Hz (Bo * 42577000)
    - One .txt file per metabolite; the folder IS the basis set
    - Compatible with: jMRUI, QUEST, FSL-MRS (via basis_tools convert)

Derivations from core struct:
    dwelltime [s]  = 1 / SW
    SamplingInterval [ms] = dwelltime * 1000 = 1000 / SW
    Bo [T]         = txfrq / 42577000
    specs          = fftshift(ifft(fids))
"""

import os
import numpy as np
from datetime import datetime

from writers.centre import fid_at_centre


######################## Constants #######################
PROTON_GAMMA = 42577000.0   # Hz/T


######################## Core writer #######################
def write_jmrui(
    fids,
    sw,
    txfrq,
    name,
    outpath,
    spectrometer='Unknown',
    additional_info='Converted by basis_set_converter',
    date=None,
):
    """
    Write one basis function to jMRUI .txt format.

    Parameters
    ----------
    fids         : np.ndarray, complex, shape (n,)
                   Time-domain FID. Water-referenced (4.65 ppm).
    sw           : float
                   Spectral width (bandwidth) in Hz.
    txfrq        : float
                   Larmor (transmitter) frequency in Hz.
    name         : str
                   Metabolite name, written to NameOfPatient field.
    outpath      : str
                   Full path of the output .txt file.
    spectrometer : str, optional
                   Scanner / spectrometer name.
    additional_info : str, optional
                   Free-text info written to AdditionalInfo field.
    date         : str or None, optional
                   Date string YYYYMMDD. Defaults to today.
    """
    fids = np.asarray(fids).squeeze()
    n         = len(fids)
    dwelltime = 1.0 / sw                    # seconds
    Bo        = txfrq / PROTON_GAMMA        # Tesla

    # jMRUI / FID-A convention: specs = fftshift(ifft(fids))
    specs = np.fft.fftshift(np.fft.ifft(fids))

    # Four-column data matrix
    RF = np.column_stack([
        fids.real,
        fids.imag,
        specs.real,
        specs.imag,
    ])

    if date is None:
        date = datetime.now().strftime('%Y%m%d')

    os.makedirs(os.path.dirname(os.path.abspath(outpath)), exist_ok=True)

    with open(outpath, 'w') as f:
        ######################## Header #######################
        f.write('jMRUI Data Textfile')
        f.write(f'\n\nFilename: {outpath}')
        f.write(f'\n\nPointsInDataset: {n}')
        f.write(f'\nDatasetsInFile: 1')
        f.write(f'\nSamplingInterval: {dwelltime * 1000:.6E}')   # ms
        f.write(f'\nZeroOrderPhase: 0.000000E+00')
        f.write(f'\nBeginTime: 0.000000E+00')
        f.write(f'\nTransmitterFrequency: {txfrq:.6E}')          # Hz
        f.write(f'\nMagneticField: {Bo:.6E}')                    # T
        f.write(f'\nTypeOfNucleus: 0.000000E+00')
        f.write(f'\nNameOfPatient: {name}')
        f.write(f'\nDateOfExperiment: {date}')
        f.write(f'\nSpectrometer: {spectrometer}')
        f.write(f'\nAdditionalInfo: {additional_info}\n\n\n')

        ######################## Data section #######################
        f.write('Signal and FFT\n')
        f.write('sig(real)\tsig(imag)\tfft(real)\tfft(imag)\n')
        f.write(f'Signal 1 out of 1 in file\n')
        for row in RF:
            f.write(
                f'{row[0]:.8E}\t{row[1]:.8E}\t'
                f'{row[2]:.8E}\t{row[3]:.8E}\n'
            )


######################## Batch writer #######################
def write_jmrui_basis_set(basis, outdir):
    """
    Write a complete basis set to a folder of jMRUI .txt files.

    Parameters
    ----------
    basis  : dict with keys:
                 'fids'   : np.ndarray, shape (n, n_metabolites), complex
                 'names'  : list of str, metabolite names
                 'sw'     : float, spectral width in Hz
                 'txfrq'  : float, Larmor frequency in Hz
    outdir : str
             Output directory. Created if it does not exist.
    """
    fids   = np.asarray(basis['fids'])
    names  = basis['names']
    sw     = float(basis['sw'])
    txfrq  = float(basis['txfrq'])

    os.makedirs(outdir, exist_ok=True)

    for i, name in enumerate(names):
        fid     = fids[:, i]
        outpath = os.path.join(outdir, f'{name}.txt')
        write_jmrui(fid, sw, txfrq, name, outpath)
        print(f'  Written: {outpath}')

    print(f'\nDone. {len(names)} basis functions written to: {outdir}')


######################## Derivation helpers #######################
def sampling_interval_ms(sw):
    """Dwell time in milliseconds from spectral width in Hz."""
    return 1000.0 / sw

def bo_from_txfrq(txfrq, gamma=PROTON_GAMMA):
    """Field strength in Tesla from Larmor frequency in Hz."""
    return txfrq / gamma

def txfrq_from_bo(bo, gamma=PROTON_GAMMA):
    """Larmor frequency in Hz from field strength in Tesla."""
    return bo * gamma



def write_jmrui_folder(basis_list, outdir):
    """
    Write all metabolites to individual jMRUI .txt files.
    Accepts a list of core struct dicts (compatible with screen4).
    """
    import os
    os.makedirs(outdir, exist_ok=True)
    written = []
    for core in basis_list:
        name    = core.get('name', 'unknown')
        outpath = os.path.join(outdir, f"{name}.txt")
        txfrq   = float(core['sf']) * 1e6   # MHz → Hz
        write_jmrui(fid_at_centre(core), core['sw'], txfrq, name, outpath)
        written.append(outpath)
        print(f"  Written: {name}.txt")
    print(f"\n  Done. {len(written)} jMRUI files written to: {outdir}")
    return written

######################## Command line usage #######################
if __name__ == '__main__':
    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    if len(sys.argv) < 3:
        print("Usage: python write_jmrui.py <input_folder> <output_folder>")
        print()
        print("Example:")
        print("  python write_jmrui.py /path/to/marss_mat_folder /path/to/jmrui_output")
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
    print(f"Output folder: {outdir}")
    print()

    os.makedirs(outdir, exist_ok=True)
    for core in basis_list:
        outpath = os.path.join(outdir, f"{core['name']}.txt")
        txfrq   = core['sf'] * 1e6   # MHz → Hz
        write_jmrui(fid_at_centre(core), core['sw'], txfrq, core['name'], outpath)
        print(f"  Written: {core['name']}.txt")

    print(f"\nDone. {len(basis_list)} files written to: {outdir}")

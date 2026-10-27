"""
writers/write_pyamares.py
Basis Set Converter — PyAMARES writer

Generates a PyAMARES prior knowledge CSV from a list of core struct dicts.
Translated directly from generate_pyamares_pk_from_basis.py.

Pipeline (confirmed from your script):
    1. Calibrate ppmCalib from Cr391 reference singlet (3.913 ppm)
       via phase evolution of consecutive FID points
    2. For each metabolite: FFT, phase-correct, find dominant peak
    3. Parabolic interpolation for sub-bin accuracy
    4. Normalize peak heights to relative initial amplitudes
    5. Write PyAMARES-compatible CSV with initial values and bounds

CSV format:
    Row 1   : Index, Metab1, Metab2, ...
    Row 2   : Initial Values
    Rows 3-7: amplitude, chemicalshift, linewidth, phase, g
    Row 8   : Bounds
    Rows 9-13: amplitude, chemicalshift, linewidth, phase, g (bounds)

Entry point:
    write_pyamares(basis_list, outpath,
                   calib_name, calib_omega) → .csv file
"""

import os
import csv
import numpy as np

from writers.centre import at_centre


############## Default prior knowledge settings #############

# Calibration reference: pure singlet with known chemical shift
CALIB_NAME  = 'Cr391'     # filename stem of reference metabolite
CALIB_OMEGA = 3.913       # known chemical shift in ppm

# Per-metabolite linewidth initial values (Hz)
LW_INIT = {
    'MM'  : 30.0, 'MMb' : 50.0, 'GABA': 8.0,  'Gln' : 8.0,
    'Glc' : 8.0,  'GSH' : 6.0,  'Asc' : 6.0,  'Asp' : 6.0,
    'PE'  : 6.0,  'Tau' : 6.0,
}
DEFAULT_LW = 5.0

# Per-metabolite chemical shift tolerance (ppm)
PPM_TOL = {
    'MM' : 0.25, 'MMb': 0.30, 'GABA': 0.12,
    'Gln': 0.12, 'Glc': 0.15,
}
DEFAULT_PPM_TOL = 0.08

# Per-metabolite linewidth bounds
LW_BOUNDS = {'MM': '(1,150)', 'MMb': '(1,150)'}
DEFAULT_LW_BOUNDS = '(1,80)'

# Amplitude normalization
AMP_SCALE   = 1.0
AMP_FLOOR   = 0.01
AMP_DEFAULT = 1.0


############## Main entry point #############

def write_pyamares(basis_list, outpath,
                   calib_name=CALIB_NAME,
                   calib_omega=CALIB_OMEGA):
    """
    Generate a PyAMARES prior knowledge CSV from a list of core structs.

    Parameters
    ----------
    basis_list  : list of core struct dicts (each has fid, sw, sf, name)
    outpath     : str, full path for output .csv file
    calib_name  : str, name of reference metabolite for ppmCalib calibration
    calib_omega : float, known chemical shift of reference metabolite in ppm
    """
    basis_list = at_centre(basis_list)     # peaks are read off the spectrum at 4.65
    outpath = os.path.abspath(outpath)
    os.makedirs(os.path.dirname(outpath) or '.', exist_ok=True)

    if not basis_list:
        raise RuntimeError("basis_list is empty")

    # Build name → core dict for easy lookup
    name_to_core = {c['name']: c for c in basis_list}

    ############## Step 1: Calibrate ppmCalib from reference metabolite #############
    if calib_name in name_to_core:
        calib_core = name_to_core[calib_name]
        ppm_calib, sf, sw = _calibrate_ppmcalib(
            calib_core['fid'],
            calib_core['sf'],
            calib_core['sw'],
            calib_omega
        )
        print(f"  ppmCalib calibrated from {calib_name}: {ppm_calib:.6f} ppm")
    else:
        # Fallback: use first metabolite's sf/sw and default ppmCalib
        first    = basis_list[0]
        sf       = float(first['sf'])
        sw       = float(first['sw'])
        ppm_calib= 4.65
        print(f"  WARNING: '{calib_name}' not found — using default ppmCalib={ppm_calib}")

    ############## Step 2: Detect peak position and amplitude per metabolite #############
    names          = [c['name'] for c in basis_list]
    resolved_ppm   = {}
    resolved_amp   = {}

    for core in basis_list:
        name = core['name']
        fid  = np.asarray(core['fid']).ravel()
        sf_i = float(core['sf'])
        sw_i = float(core['sw'])

        try:
            ppm = _find_peak_ppm(fid, sf_i, sw_i, ppm_calib)
            amp = _find_peak_amp(fid, sf_i, sw_i, ppm_calib)

            resolved_ppm[name] = ppm if ppm is not None else 3.000
            resolved_amp[name] = amp if (amp is not None and
                                          np.isfinite(amp) and amp > 0) \
                                     else AMP_DEFAULT
        except Exception as ex:
            print(f"  WARNING: error for {name}: {ex} — using defaults")
            resolved_ppm[name] = 3.000
            resolved_amp[name] = AMP_DEFAULT

    ############## Step 3: Normalize amplitudes #############
    valid_amps = [v for v in resolved_amp.values()
                  if np.isfinite(v) and v > 0]
    amp_ref    = max(valid_amps) if valid_amps else 1.0

    amp_init = {}
    for name in names:
        rel = resolved_amp[name] / amp_ref
        amp_init[name] = max(AMP_FLOOR, AMP_SCALE * rel)

    ############## Step 4: Write CSV #############
    def lw(n):  return LW_INIT.get(n, DEFAULT_LW)
    def tol(n): return PPM_TOL.get(n, DEFAULT_PPM_TOL)
    def lwb(n): return LW_BOUNDS.get(n, DEFAULT_LW_BOUNDS)
    def p(n):   return resolved_ppm[n]

    rows = []

    # Header
    rows.append(['Index'] + names)

    # Initial values
    rows.append(['Initial Values'] + [''] * len(names))
    rows.append(['amplitude']     + [amp_init[n] for n in names])
    rows.append(['chemicalshift'] + [p(n)         for n in names])
    rows.append(['linewidth']     + [lw(n)         for n in names])
    rows.append(['phase']         + [0]            * len(names))
    rows.append(['g']             + [0]            * len(names))

    # Bounds
    rows.append(['Bounds'] + [''] * len(names))
    rows.append(['amplitude']     + ['(0,'] * len(names))
    rows.append(['chemicalshift'] + [
        f"({p(n)-tol(n):.4f},{p(n)+tol(n):.4f})" for n in names
    ])
    rows.append(['linewidth']     + [lwb(n) for n in names])
    rows.append(['phase']         + ['(-180,180)'] * len(names))
    rows.append(['g']             + ['(0,1)']       * len(names))

    with open(outpath, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(rows)

    print(f"\n  Done. PyAMARES CSV written: {os.path.basename(outpath)}")
    print(f"  {len(names)} peaks")


############## Helper functions (from generate_pyamares_pk_from_basis.py) #############

def _calibrate_ppmcalib(fid, sf, sw, calib_omega):
    """
    Estimate ppmCalib from a reference singlet FID.
    Confirmed from your generate_pyamares_pk_from_basis.py.
    """
    fid  = np.asarray(fid).ravel()
    npts = min(len(fid) - 1, 199)
    if npts < 2:
        raise RuntimeError("FID too short for ppmCalib calibration")

    phase_diffs = np.angle(fid[1:npts+1] * np.conj(fid[:npts]))
    f0_hz       = np.mean(phase_diffs) / (2 * np.pi) * sw
    ppm_calib   = calib_omega - f0_hz / sf

    return ppm_calib, float(sf), float(sw)


def _get_phased_spectrum(fid, sf, sw, ppm_calib):
    """Build phased real spectrum and ppm axis."""
    fid   = np.asarray(fid).ravel()
    n     = len(fid)
    freq  = np.fft.fftshift(np.fft.fftfreq(n, d=1.0 / sw))
    ppm   = ppm_calib + freq / sf
    spec  = np.fft.fftshift(np.fft.fft(fid))
    phi0  = -np.angle(fid[0]) if fid[0] != 0 else 0.0
    spec_r= np.real(spec * np.exp(1j * phi0))
    return ppm, spec_r


def _find_peak_ppm(fid, sf, sw, ppm_calib, ppm_min=0.5, ppm_max=4.5):
    """Find dominant peak position with parabolic interpolation."""
    ppm, spec_r = _get_phased_spectrum(fid, sf, sw, ppm_calib)
    n           = len(spec_r)
    mask        = (ppm >= ppm_min) & (ppm <= ppm_max)
    idx_r       = np.where(mask)[0]

    if len(idx_r) == 0:
        return None

    pidx = idx_r[np.argmax(spec_r[mask])]

    if 0 < pidx < n - 1:
        y0, y1, y2 = spec_r[pidx-1], spec_r[pidx], spec_r[pidx+1]
        denom = 2.0 * (y0 - 2.0 * y1 + y2)
        delta = (y0 - y2) / denom if denom != 0 else 0.0
        ppm_step = ppm[1] - ppm[0]
        return round(float(ppm[pidx] - delta * ppm_step), 4)

    return round(float(ppm[pidx]), 4)


def _find_peak_amp(fid, sf, sw, ppm_calib, ppm_min=0.5, ppm_max=4.5):
    """Measure dominant peak height."""
    ppm, spec_r = _get_phased_spectrum(fid, sf, sw, ppm_calib)
    mask        = (ppm >= ppm_min) & (ppm <= ppm_max)
    if not np.any(mask):
        return None
    return float(np.max(spec_r[mask]))


############## Command line usage #############

if __name__ == '__main__':
    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    if len(sys.argv) < 3:
        print("Usage: python write_pyamares.py <input_folder> <output_csv>")
        sys.exit(1)

    input_path = sys.argv[1]
    outpath    = sys.argv[2]

    from readers.read_marss import read_marss_file, read_marss_folder

    if os.path.isdir(input_path):
        basis_list = read_marss_folder(input_path)
    else:
        basis_list = [read_marss_file(input_path)]

    print(f"Loaded {len(basis_list)} metabolites\n")
    write_pyamares(basis_list, outpath)

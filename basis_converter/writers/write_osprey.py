"""
writers/write_osprey.py
Basis Set Converter — Osprey writer

Assembles a list of core struct dicts into an Osprey-compatible BASIS .mat file, as
Osprey's own sets hold one: its metabolites, and its macromolecules, lipids and water.
The functions, the macromolecules and lipids and their scale are Osprey's own
fit_makeBasis (libraries/FID-A/fitTools), so a set gets the basis Osprey would build
from the same functions:

Pipeline:
    1. Apply name mapping (MARSS names -> Osprey names), summing MARSS's creatine parts
    2. Skip SKIP entries and duplicates
    3. Zero-pad FIDs to targetN if needed
    4. DC-correct each function as op_dccorr(.., 'p') does: the mean of the spectrum in
       the 0.5 ppm at each end of its axis, subtracted
    5. Compute specs = fftshift(fft(fids))   <- Osprey uses fft NOT ifft
    6. One-proton area as detCrArea: a Lorentzian fitted to Cr's 3.027 ppm singlet, / 3
    7. Rescale MM/Lip the set brings to proton-equivalent amplitudes
    8. Add synthetic H2O Lorentzian at 4.68 ppm (Osprey's own sets hold one, for water fits)
    9. Optionally add Osprey's MM/Lip Gaussians (op_gaussianPeak) for missing components
   10. Reorder: metabolites first, MM/Lip/H2O last
   11. Normalize: divide by max(real(specs)), Osprey's BASIS.scale
   12. Assemble and save BASIS struct as .mat, with the set as BasiC was given it (BASIS.BasiC,
       which Osprey does not read), so that read back the set is the one given

Entry point:
    write_osprey(basis_list, outpath, target_n, add_mm, te, sequence)
"""

import os
import numpy as np
import scipy.io as sio

from writers.centre import CENTRE


############# Name mapping (MARSS filename stems -> Osprey names) #############
# 'SKIP' means the metabolite is excluded from the Osprey basis set

NAME_MAP = {
    'Asc'     : 'Asc',
    'Asp'     : 'Asp',
    'Cho'     : 'PCh',
    'Cr391'   : 'Cr',        # parts, summed (PARTS)
    'CrNo391' : 'Cr',
    'GABA'    : 'GABA',
    'Glc'     : 'Glc',
    'Gln'     : 'Gln',
    'Glu'     : 'Glu',
    'Gly'     : 'Gly',
    'GPC'     : 'GPC',
    'GSH'     : 'GSH',
    'Lac'     : 'Lac',
    'ml'      : 'mI',
    'mI'      : 'mI',
    'MM'      : 'SKIP',      # composite MM — replaced by parametric components
    'MM09'    : 'MM09',      # a macromolecule; it takes the place of Osprey's own MM09
    'NAA'     : 'NAA',
    'NAAG'    : 'NAAG',
    'PCr393'  : 'PCr',       # parts, summed (PARTS)
    'PCrNo393': 'PCr',
    'PE'      : 'PE',
    'sl'      : 'sI',
    'sI'      : 'sI',
    'Tau'     : 'Tau',
}

# MARSS splits Cr and PCr into the 3.9 ppm CH2 and the rest; Osprey's Cr and PCr hold both
# (its CrCH2 is a negative CH2 correction term), so the parts are summed into one
PARTS = {'Cr391', 'CrNo391', 'PCr393', 'PCrNo393'}

# names as LCModel keeps them, six characters (an FSL-MRS set from a .BASIS holds CrNo39, PCrNo3)
SHORT = {k[:6]: k for k in NAME_MAP if len(k) > 6}

# MM/Lip proton-equivalent counts for rescaling
MM_LIP_SCALING = {
    'MM09'  : 3,
    'MM12'  : 2,
    'MM14'  : 2,
    'MM17'  : 2,
    'MM20'  : 11.4,
    'Lip09' : 3,
    'Lip13' : 4,
    'Lip20' : 2.87,
}

# Osprey MM/Lip names — used for reordering
MM_NAMES = {'MM09','MM12','MM14','MM17','MM20',
            'Lip09','Lip13','Lip20','H2O','MM'}


############# Main entry point #############
def write_osprey(basis_list, outpath,
                 target_n=0, add_mm=True,
                 te=30.0, sequence='unedited'):
    """
    Assemble and save an Osprey BASIS .mat file.

    Parameters
    ----------
    basis_list : list of core struct dicts
        Each dict must have: fid, sw, sf, n, name
    outpath    : str, full path for output .mat file
    target_n   : int, target spectral points (0 = native size)
    add_mm     : bool, fit_makeBasis's addMMFlag: add Osprey's MM/Lip for missing components, in
                 place of a set's composite MM, and rescale the set's own; False keeps the set's
                 MM as they are and adds none (the H2O Osprey's water fit takes is added either way)
    te         : float, echo time in ms
    sequence   : str, e.g. 'unedited'
    """
    outpath = os.path.abspath(outpath)
    os.makedirs(os.path.dirname(outpath) or '.', exist_ok=True)

    if not basis_list:
        raise RuntimeError("basis_list is empty")

    ############# Step 1: Apply name mapping, skip duplicates #############
    met_names = []
    met_fids  = []
    sw = sf = n_native = center_freq = None

    for core in basis_list:
        stem     = SHORT.get(core['name'], core['name'])
        osp_name = NAME_MAP.get(stem, stem)   # use original name if not in map

        if osp_name == 'SKIP' and not add_mm:     # nothing takes its place: kept, as it is named
            osp_name = stem
        elif osp_name == 'SKIP':
            print(f"  Skipping (SKIP): {stem}")
            continue

        if osp_name in met_names:
            if stem in PARTS:
                met_fids[met_names.index(osp_name)] = \
                    met_fids[met_names.index(osp_name)] + np.asarray(core['fid']).ravel().astype(complex)
                print(f"  Summed: {stem:15s} -> {osp_name}")
            else:
                print(f"  Skipping duplicate '{osp_name}' from: {stem}")
            continue

        # Store reference parameters from first valid metabolite
        if sw is None:
            sw          = float(core['sw'])
            sf          = float(core['sf'])          # MHz
            n_native    = int(core['n'])
            # one stating none is at 4.65 (a ppmCalib is a label: MARSS's 4.675 sets are at 4.65)
            center_freq = float(core.get('centerFreq') or CENTRE)

        met_names.append(osp_name)
        met_fids.append(np.asarray(core['fid']).ravel().astype(complex))
        print(f"  Loaded: {stem:15s} -> {osp_name}")

    if not met_names:
        raise RuntimeError("No valid basis functions loaded after name mapping.")

    ############# Step 2: Determine final spectral size #############
    if target_n == 0:
        final_n = n_native
    elif target_n < n_native:
        print(f"  WARNING: target_n ({target_n}) < native n ({n_native}); using native.")
        final_n = n_native
    else:
        final_n = target_n

    print(f"\n  Native n={n_native}, final n={final_n}")

    ############# Step 3: Spectral axes #############
    dt      = 1.0 / sw
    Bo      = sf / 42.577
    hzppm   = sf               # numerically equal (sf in MHz = Hz/ppm)
    t       = np.arange(final_n) * dt
    ppm_axis= _ppm_axis(final_n, sw, hzppm, center_freq)

    ############# Step 4: Build fids / specs with DC correction #############
    n_mets = len(met_names)
    fids   = np.zeros((final_n, n_mets), dtype=complex)
    specs  = np.zeros((final_n, n_mets), dtype=complex)

    for kk, fid in enumerate(met_fids):
        pad_fid = np.zeros(final_n, dtype=complex)
        copy_len = min(len(fid), final_n)
        pad_fid[:copy_len] = fid[:copy_len]

        fids[:, kk], specs[:, kk] = _dccorr(pad_fid, ppm_axis)

    ############# Step 5: One-proton area estimate #############
    one_proton_area = _cr_area(specs, met_names, ppm_axis, Bo) / 3

    ############# Step 6: Rescale MM/Lip to proton-equivalent amplitudes #############
    if one_proton_area > 0 and add_mm:
        print("\n  Rescaling MM/Lip basis functions:")
        for kk, name in enumerate(met_names):
            if name in MM_LIP_SCALING:
                n_prot       = MM_LIP_SCALING[name]
                current_area = np.sum(np.real(specs[:, kk]))
                if abs(current_area) > 0:
                    scale         = (n_prot * one_proton_area) / current_area
                    fids[:, kk]  *= scale
                    specs[:, kk] *= scale
                    print(f"    {name:8s}: scale={scale:.4f}  ({n_prot} proton-equiv)")
    elif add_mm:
        print("  WARNING: oneProtonArea=0; MM/Lip functions will NOT be rescaled.")

    ############# Step 7: Add synthetic H2O #############
    if 'H2O' not in met_names:
        print("\n  Adding synthetic H2O Lorentzian at 4.68 ppm.")
        # two protons by area, as Osprey's MM/Lip Gaussians are (its own sets' H2O is below their
        # metabolites); an amplitude of two one-proton areas made it some 40 times their height,
        # and the set, normalised to it, then fitted at 1/40 of Osprey's scale
        _, unit = _make_lorentzian(final_n, sw, hzppm, 4.68, 0.1 * hzppm, center_freq, 1.0)
        amp = 2 * one_proton_area / np.sum(np.real(unit)) if one_proton_area > 0 \
              else 0.05 * np.max(np.real(specs)) / np.max(np.real(unit))
        h2o_fid, h2o_spec = _make_lorentzian(
            final_n, sw, hzppm, 4.68, 0.1 * hzppm, center_freq, amp
        )
        fids  = np.column_stack([fids,  h2o_fid])
        specs = np.column_stack([specs, h2o_spec])
        met_names.append('H2O')
        n_mets += 1

    ############# Step 8: Optionally add parametric MM/Lip #############
    added_mm_names = []
    if add_mm and one_proton_area > 0:
        print("\n  Adding parametric MM/Lip for missing components:")
        mm_fids, mm_specs, mm_names = _build_mm_lipids(
            final_n, sw, hzppm, center_freq, one_proton_area
        )
        for rr, mm_name in enumerate(mm_names):
            if mm_name not in met_names:
                fids  = np.column_stack([fids,  mm_fids[:, rr]])
                specs = np.column_stack([specs, mm_specs[:, rr]])
                met_names.append(mm_name)
                added_mm_names.append(mm_name)
                n_mets += 1
                print(f"    Added: {mm_name}")
            else:
                print(f"    Skipped (already present): {mm_name}")

    ############# Step 9: Reorder — metabolites first, MM/Lip/H2O last #############
    known_mm = MM_NAMES | set(added_mm_names)
    is_mm    = [name in known_mm for name in met_names]
    met_idx  = [i for i, m in enumerate(is_mm) if not m]
    mm_idx   = [i for i, m in enumerate(is_mm) if m]
    reorder  = met_idx + mm_idx

    fids      = fids[:, reorder]
    specs     = specs[:, reorder]
    met_names = [met_names[i] for i in reorder]

    n_mm_count  = len(mm_idx)
    n_mets_only = len(met_idx)

    ############# Step 10: Normalize #############
    scale_factor = np.max(np.real(specs))
    if scale_factor > 0:
        fids  /= scale_factor
        specs /= scale_factor
    else:
        scale_factor = 1.0

    ############# Step 11: Assemble BASIS struct #############
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
        'leftshifted'    : 1.0,
        'filtered'       : 0.0,
        'zeropadded'     : float(final_n > n_native),
        'freqcorrected'  : 0.0,
        'phasecorrected' : 0.0,
        'averaged'       : 1.0,
        'addedrcvrs'     : 1.0,
        'subtracted'     : 1.0,
        'writtentotext'  : 0.0,
        'downsampled'    : 0.0,
        'isISIS'         : 0.0,
        'freqranged'     : 1.0,
        'addedMM'        : float(add_mm),
    }

    BASIS = {
        'spectralwidth' : sw,
        'dwelltime'     : dt,
        'n'             : float(final_n),
        'linewidth'     : next((float(c['linewidth']) for c in basis_list      # the first stated, Hz
                                if c.get('linewidth') is not None), 1.0),
        'Bo'            : Bo,
        'seq'           : [sequence],
        'te'            : float(te),
        'centerFreq'    : center_freq,
        'ppm'           : ppm_axis.reshape(-1, 1),
        't'             : t.reshape(-1, 1),
        'flags'         : flags,
        'nMM'           : float(n_mm_count),
        'fids'          : fids,
        'specs'         : specs,
        'name'          : np.array(met_names, dtype=object).reshape(1, -1),   # a cell of strings, not padded chars
        'dims'          : dims,
        'nMets'         : float(n_mets_only),
        'sz'            : np.array([final_n, n_mets], dtype=float),   # double, as MATLAB's size() is
        'scale'         : scale_factor,
        'BasiC'         : _given(basis_list, fids),
    }

    ############# Step 12: Save #############
    sio.savemat(outpath, {'BASIS': BASIS})

    print(f"\n{'='*45}")
    print(f"Saved: {os.path.basename(outpath)}")
    print(f"  nMets : {n_mets_only}")
    print(f"  nMM   : {n_mm_count}")
    print(f"  Total : {n_mets}")
    print(f"  n     : {final_n}")
    print(f"  Bo    : {Bo:.3f} T")
    print(f"  Names : {', '.join(met_names)}")
    print(f"{'='*45}")


def _given(basis_list, fids):
    """The set as BasiC was given it: what Osprey's model changes (names, creatine's parts, the
    first point, the scale) and what it adds (MM, lipids, water) undone, read back. The names, the
    FIDs (zero-filled to the longest) and their lengths, which of them was an Osprey set's water,
    and the sum of each function's magnitude as written, by which the reader tells that the
    functions are still those it wrote."""
    n = [len(np.ravel(c['fid'])) for c in basis_list]
    given = np.zeros((max(n), len(n)), dtype=complex)
    for k, c in enumerate(basis_list):
        given[:n[k], k] = np.ravel(c['fid'])
    return {'name' : np.array([c['name'] for c in basis_list], dtype=object).reshape(1, -1),
            'fids' : given,
            'n'    : np.array(n, dtype=float),
            'water': np.array([bool(c.get('water_ref')) for c in basis_list], dtype=float),
            'check': np.sum(np.abs(fids), axis=0)}


############# Helper functions (translated from MATLAB) #############
def _ppm_axis(n, sw, hzppm, center_freq):
    """Osprey's ppm axis (fit_makeBasis, op_gaussianPeak): bin centres, half a bin off 0 Hz."""
    f = -sw / 2 + sw / (2 * n) + np.arange(n) * sw / n
    return f / hzppm + center_freq


def _dccorr(fid, ppm):
    """op_dccorr(.., 'p'): the mean of the spectrum in the 0.5 ppm at each end of the axis,
    subtracted; the FID and the spectrum."""
    spec = np.fft.fftshift(np.fft.fft(fid))
    hi = (ppm < ppm.max()) & (ppm > ppm.max() - 0.5)
    lo = (ppm < ppm.min() + 0.5) & (ppm > ppm.min())
    spec = spec - np.mean([np.mean(spec[hi]), np.mean(spec[lo])])
    shifted = np.fft.fftshift(spec) if len(spec) % 2 == 0 else np.roll(np.fft.fftshift(spec), 1)
    return np.fft.ifft(shifted), spec


def _lorentz(ppm, amp, fwhm, ppm0, offset, phase):
    """op_lorentz with five parameters: a Lorentzian of height amp, phased (degrees)."""
    gam = fwhm / 2
    y = np.sqrt(2 / np.pi) * (gam - 1j * (ppm - ppm0)) / (gam ** 2 + (ppm - ppm0) ** 2)
    return np.real((y / np.max(np.abs(y)) * amp + offset) * np.exp(1j * phase * np.pi / 180))


def _cr_area(specs, names, ppm, Bo):
    """
    detCrArea: the area of a Lorentzian fitted to Cr's 3.027 ppm singlet (three protons),
    summed over 3.027 +/- 0.4 ppm; 0 without Cr (no MM/Lip scaling then).
    """
    from scipy.optimize import curve_fit
    if 'Cr' not in names:
        print("  WARNING: Cr not found — MM/Lip scaling disabled.")
        return 0.0
    window = (ppm > 3.027 - 0.4) & (ppm < 3.027 + 0.4)
    ref, at = specs[window, names.index('Cr')], ppm[window]
    peak = int(np.argmax(np.abs(np.real(ref))))
    guess = [np.real(ref[peak]), (5 * Bo / 3) / (42.577 * Bo), at[peak], 0.0, 0.0]
    try:
        pars, _ = curve_fit(_lorentz, at, np.real(ref), p0=guess, maxfev=20000)
    except RuntimeError:
        pars = guess
    area = float(np.sum(_lorentz(at, *pars)))
    print(f"\n  Cr area (Lorentzian at {pars[2]:.3f} ppm)={area:.4e}, oneProtonArea={area / 3:.4e}")
    return area


def _make_lorentzian(n, sw, hzppm, freq_ppm, fwhm_hz, center_freq, amp):
    """Synthetic Lorentzian FID and spectrum."""
    dt      = 1.0 / sw
    t       = np.arange(n) * dt
    freq_hz = (freq_ppm - center_freq) * hzppm
    T2      = 1.0 / (np.pi * fwhm_hz)
    fid     = amp * np.exp(1j * 2 * np.pi * freq_hz * t) * np.exp(-t / T2)
    spec    = np.fft.fftshift(np.fft.fft(fid))
    return fid, spec


def _make_gaussian(n, sw, hzppm, freq_ppm, lw_hz, center_freq, amp):
    """op_gaussianPeak: a Gaussian FID of Osprey's width lw (decay 1 / (pi lw)), DC-corrected
    as fit_makeBasis does each; the FID and the spectrum."""
    t       = np.arange(n) / sw
    freq_hz = (freq_ppm - center_freq) * hzppm
    decay   = 1.0 / (lw_hz * np.pi)
    fid     = amp * np.exp(1j * 2 * np.pi * freq_hz * t) * np.exp(-(t**2) / (2 * decay**2))
    return _dccorr(fid, _ppm_axis(n, sw, hzppm, center_freq))


def _build_mm_lipids(n, sw, hzppm, center_freq, one_proton_area):
    """
    Osprey's parametric Gaussian MM/Lip basis functions, as fit_makeBasis builds them
    (Wilson et al. MRM 2011): each in proton equivalents of the one-proton area, against
    the area of a unit Gaussian 0.1 ppm wide at 0 ppm; the sums DC-corrected as a whole.
    """
    _, g0 = _make_gaussian(n, sw, hzppm, 0.0, 0.1 * hzppm, center_freq, 1.0)
    g_area = np.sum(np.real(g0)) or 1.0

    def G(ppm_pos, fwhm_frac, n_prot):
        amp = n_prot * one_proton_area / g_area
        return _make_gaussian(n, sw, hzppm, ppm_pos, fwhm_frac * hzppm, center_freq, amp)

    def summed(*parts):
        return _dccorr(sum(f for f, _ in parts), _ppm_axis(n, sw, hzppm, center_freq))

    f09, s09 = G(0.91, 0.14, 3)
    f12, s12 = G(1.21, 0.15, 2)
    f14, s14 = G(1.43, 0.17, 2)
    f17, s17 = G(1.67, 0.15, 2)

    f20, s20 = summed(G(2.08, 0.15, 1.33), G(2.25, 0.20, 0.33), G(1.95, 0.15, 0.33), G(3.00, 0.20, 0.40))

    fl09, sl09 = G(0.89, 0.14, 3)
    fl13, sl13 = summed(G(1.28, 0.15, 2), G(1.28, 0.89, 2))
    fl20, sl20 = summed(G(2.04, 0.15, 1.33), G(2.25, 0.15, 0.67), G(2.80, 0.20, 0.87))

    mm_fids  = np.column_stack([f09, f12, f14, f17, f20, fl09, fl13, fl20])
    mm_specs = np.column_stack([s09, s12, s14, s17, s20, sl09, sl13, sl20])
    mm_names = ['MM09', 'MM12', 'MM14', 'MM17', 'MM20', 'Lip09', 'Lip13', 'Lip20']

    return mm_fids, mm_specs, mm_names


############# Command line usage #############
if __name__ == '__main__':
    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    if len(sys.argv) < 3:
        print("Usage: python write_osprey.py <input_folder> <output_file> [target_n] [te]")
        print()
        print("Example:")
        print("  python write_osprey.py /path/to/marss_folder /path/to/output/basis.mat 2048 26")
        sys.exit(1)

    input_path = sys.argv[1]
    outpath    = sys.argv[2]
    target_n   = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    te         = float(sys.argv[4]) if len(sys.argv) > 4 else 30.0

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

    write_osprey(basis_list, outpath, target_n=target_n, te=te)
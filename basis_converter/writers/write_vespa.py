"""
writers/write_vespa.py
Basis Set Converter — VeSPA Analysis writer

Generates a VeSPA Analysis prior XML file whose lines rebuild the basis set.

VeSPA builds each metabolite's basis FID from its prior lines (common/mrs_prior.py,
Prior.calculate_full_basis_set) with common/util/generic_spectral.py create_fid: every line
is an undamped complex exponential,

    area * exp(1j * phase * pi / 180) * exp(1j * 2 pi * ppm2hz(ppm) * t),
    ppm2hz(ppm) = sw / 2 - frequency * (ppm - resppm)            (analysis/mrs_dataset.py),

area the amplitude at t = 0, phase in degrees, ppm absolute; the fit adds the lineshape. The
sw / 2 is VeSPA's chop, (-1)**k; the rest runs opposite to a basis' own convention, where a
line at ppm is at +sf * (ppm - centre) Hz (writers/centre.py): VeSPA's FID is the complex
conjugate of BasiC's. So a line BasiC holds as c * exp(2 pi i f t) is written at
ppm = centre + f / sf, area |c| and phase -angle(c) in degrees; rebuilt by create_fid it is
the chopped conjugate of the BasiC line, the same spectrum on VeSPA's ppm axis.
VeSPA-Simulation does the same: its prior lines are at -Fr ppm with phase -angle(I)
(common/resources/pulse_sequences.xml:366, 373).

The lines are each metabolite's FID decomposed into decaying exponentials (writers/lines.py,
HLSVD); their dampings are not written. The areas keep the source's scale for every
metabolite (no normalisation), so fitted areas compare between metabolites as the source's
FIDs do.

VeSPA prior format:
    <vespa_export version="1.0.0">
      <prior version="1.0.0">
        <nucleus>1H</nucleus>
        <seqte>26.0</seqte>
        <prior_metabolite version="1.0.0">
          <name>NAA</name>
          <dim>0</dim><dim>0</dim><dim>0</dim>
          <spins>5</spins>           ← the spins metabList lists (else the lines)
          <line>
            <group/>
            <ppm>2.00799997</ppm>   ← the line's frequency, absolute ppm
            <area>29.61939194</area> ← its amplitude at t = 0, the source's scale
            <phase>3.216372</phase>  ← -(its phase), degrees (jMRUI sample's NAA CH3)
          </line>
        </prior_metabolite>
      </prior>
    </vespa_export>

The literature table and the peak helpers below are no longer used for VeSPA's lines; the
GAVA, MIDAS and OXSA writers import them, and VeSPA uses the table for spins only.
"""

import os
import numpy as np
from datetime import datetime
import xml.etree.ElementTree as ET
from xml.dom import minidom
import scipy.io as sio

from writers.centre import CENTRE, at_centre
from writers.lines import MAX_ORDER, TOL, decompose


############## Literature chemical shift table (from metabList.mat) #############
# One entry per spin of metabList's spin system (a methyl is one). write_vespa counts them
# for <spins>; write_gava, write_midas and write_oxsa look for peaks at them.
# Loaded dynamically from metabList.mat if available, otherwise uses hardcoded
# values from the working marss_to_vespa_basis.xml

METAB_PPM_TABLE = {
    # name → list of ppm values (99 = phantom spin, excluded)
    'Asc'     : [4.492, 4.002, 3.743, 3.716],
    'Asp'     : [3.891, 2.801, 2.653],
    'Cho'     : [3.185, 4.054, 4.054, 3.501, 3.501],
    'Ch'      : [3.185, 4.054, 4.054, 3.501, 3.501],
    'PCh'     : [3.208, 4.281, 4.281, 3.641, 3.641],
    'Cr'      : [3.027, 3.913, 6.649],
    'Cr391'   : [3.913, 3.027, 6.649],   # CH2 group at 3.913
    'CrNo391' : [3.027, 3.913, 6.649],   # CH3 group at 3.027
    'PCr'     : [3.029, 3.930, 6.581, 7.296],
    'PCr393'  : [3.930, 3.029, 6.581, 7.296],
    'PCrNo393': [3.029, 3.930, 6.581, 7.296],
    'GABA'    : [2.284, 2.284, 1.889, 1.889, 3.013, 3.013],
    'Glc'     : [5.216, 3.519, 3.698, 3.395, 3.822, 3.826, 3.749],
    'Gln'     : [3.753, 2.129, 2.109, 2.432, 2.454, 6.816, 7.529],
    'Glu'     : [3.743, 2.038, 2.120, 2.338, 2.352],
    'Gly'     : [3.548],
    'GPC'     : [3.212, 3.605, 3.672, 3.903, 3.871, 3.946],
    'GSH'     : [3.769, 7.514, 4.561, 2.926, 2.975, 8.177,
                 3.769, 2.159, 2.146, 2.510, 2.560],
    'Lac'     : [4.097, 1.314, 1.314, 1.314],
    'mI'      : [3.522, 4.054, 3.522, 3.614, 3.269, 3.614],
    'ml'      : [3.522, 4.054, 3.522, 3.614, 3.269, 3.614],
    'MM'      : [0.916, 1.245, 1.430, 1.676, 1.986, 2.293,
                 2.676, 3.010, 3.210, 3.760, 3.936, 4.010],
    'MM09'    : [0.916],
    'NAA'     : [2.008, 4.382, 2.673, 2.486, 7.821],
    'NAAG'    : [2.042, 4.128, 1.881, 2.049, 2.190, 2.180, 4.607, 2.721, 2.519],
    'PE'      : [3.976, 3.976, 3.216, 3.216],
    'sI'      : [3.340],
    'sl'      : [3.340],
    'Tau'     : [3.421, 3.421, 3.246, 3.246],
    'NAAG'    : [2.042, 4.128, 1.881, 2.049, 2.190, 4.607, 2.721, 2.519],
}

PPM_CALIB = 4.65   # the centre the set is read at (a ppmCalib is a label: MARSS's 4.675 sets are at 4.65)


############## Main entry point #############

def write_vespa(basis_list, outpath,
                seqte=26.0,
                source='basis_set_converter',
                metablist_path=None):
    """
    Write a VeSPA Analysis prior XML file.

    Parameters
    ----------
    basis_list     : list of core struct dicts
    outpath        : str, full path for output .xml file
    seqte          : float, echo time in ms
    source         : str, source description
    metablist_path : str or None, path to metabList.mat (its numberOfProtons give <spins>)
    """
    basis_list = at_centre(basis_list)     # the lines' ppm are read at CENTRE
    outpath = os.path.abspath(outpath)
    os.makedirs(os.path.dirname(outpath) or '.', exist_ok=True)

    if not basis_list:
        raise RuntimeError("basis_list is empty")

    ppm_table = _load_metablist(metablist_path) if metablist_path else METAB_PPM_TABLE
    counted = 'numberOfProtons summed from metabList.mat' if metablist_path else \
              "the spins metabList's table lists (one per shift, a methyl as one)"

    timestamp = datetime.now().strftime('%Y-%m-%dT%H:%M:%S')

    root = ET.Element('vespa_export', version='1.0.0')
    ET.SubElement(root, 'timestamp').text = timestamp
    ET.SubElement(root, 'comment').text   = \
        'Converted from basis set by basis_set_converter'

    prior = ET.SubElement(root, 'prior', version='1.0.0')
    ET.SubElement(prior, 'source').text    = source
    ET.SubElement(prior, 'source_id').text = 'basis_set_converter'
    comment = ET.SubElement(prior, 'comment')
    ET.SubElement(prior, 'nucleus').text   = '1H'
    ET.SubElement(prior, 'seqte').text     = str(float(seqte))

    listed, unlisted, missed = [], [], []
    for core in basis_list:
        met_name = core.get('name', 'unknown')
        fid      = np.asarray(core['fid']).ravel()
        sf       = float(core.get('sf') or 123.26)
        sw       = float(core.get('sw') or 4000.0)

        lines = decompose(fid, sw)
        keep  = np.abs(lines['amplitude']) > 0
        amp, freq = lines['amplitude'][keep], lines['freq'][keep]
        if not len(amp):
            print(f"  WARNING: no lines for {met_name} — skipping")
            continue
        if lines['residual'] > TOL:
            # its best lines do not rebuild it: written undamped, they would be another shape
            print(f"  WARNING: no {MAX_ORDER} lines rebuild {met_name} to {TOL:g} "
                  f"(residual {lines['residual']:.1e}) — left out")
            missed.append(met_name)
            continue

        lit = _get_lit_ppms(met_name, ppm_table)
        spins = int(round(sum(n for _, n in lit))) if lit else len(amp)
        (listed if lit else unlisted).append(met_name)

        met_elem = ET.SubElement(prior, 'prior_metabolite', version='1.0.0')
        ET.SubElement(met_elem, 'name').text  = met_name
        ET.SubElement(met_elem, 'dim').text   = '0'
        ET.SubElement(met_elem, 'dim').text   = '0'
        ET.SubElement(met_elem, 'dim').text   = '0'
        ET.SubElement(met_elem, 'spins').text = str(spins)

        for c, f in zip(amp, freq):
            line = ET.SubElement(met_elem, 'line')
            ET.SubElement(line, 'group')
            ET.SubElement(line, 'ppm').text   = f'{CENTRE + f / sf:.8f}'
            ET.SubElement(line, 'area').text  = f'{abs(c):.10g}'
            ET.SubElement(line, 'phase').text = f'{-np.degrees(np.angle(c)):.6f}'

        print(f"  {met_name:15s}  {len(amp):3d} lines  "
              f"largest={CENTRE + freq[0] / sf:.4f} ppm  residual={lines['residual']:.1e}")

    comment.text = (
        f'Prior information from basis set with {len(basis_list)} metabolites. '
        f'Lines: each FID decomposed into decaying exponentials (HLSVD), the fewest that '
        f'rebuild it to a relative residual of {TOL:g}; ppm, amplitude at t = 0 and phase '
        f'of each, its damping left to the fit. Areas in the source\'s scale, the same for '
        f'every metabolite (factor 1). Phases in VeSPA\'s sign (the conjugate FID). '
        f'Spins: {counted}' + (f' for {", ".join(listed)}' if unlisted and listed else '') +
        (f'; the number of lines for {", ".join(unlisted)}' if unlisted else '') + '.' +
        (f' Left out, as no {MAX_ORDER} lines rebuild them to {TOL:g}: {", ".join(missed)}.'
         if missed else ''))

    # Write XML
    xml_str = minidom.parseString(
        ET.tostring(root, encoding='utf-8')
    ).toprettyxml(indent='  ', encoding='utf-8')

    with open(outpath, 'wb') as f:
        f.write(xml_str)

    print(f"\n  Done. VeSPA prior XML: {os.path.basename(outpath)}")
    print(f"  {len(basis_list)} metabolites  TE={seqte} ms")


############## Helpers (the GAVA, MIDAS and OXSA writers import these) #############

def _load_metablist(path):
    """Load PPM table from metabList.mat — returns dict of name → [(ppm, n_protons), ...]"""
    try:
        mat   = sio.loadmat(path, simplify_cells=True)
        mlist = mat.get('metabList', [])
        table = {}
        for entry in mlist:
            name      = entry.get('name', '')
            omega     = entry.get('Omega', [])
            n_protons = entry.get('numberOfProtons', 1)
            lines     = _flatten_omega_with_protons(omega, n_protons)
            # Filter phantom spins (>=90)
            lines = [(p, n) for p, n in lines if p < 90]
            if name and lines:
                table[name] = lines
        print(f"  Loaded metabList: {len(table)} metabolites")
        return table
    except Exception as ex:
        print(f"  WARNING: could not load metabList: {ex} — using hardcoded table")
        return {k: [(p, 1) for p in v] for k, v in METAB_PPM_TABLE.items()}


def _flatten_omega_with_protons(omega, n_protons):
    """
    Flatten nested Omega arrays pairing each ppm with its proton count.
    omega can be: float, array of floats, or array of (float|array) groups.
    n_protons can be: int or array matching omega groups.
    """
    result = []
    omega     = np.asarray(omega, dtype=object)
    n_protons = np.asarray(n_protons, dtype=object)

    if omega.ndim == 0:
        # Single scalar ppm, single proton count
        p = float(omega)
        n = int(float(n_protons)) if n_protons.ndim == 0 else 1
        if p < 90:
            result.append((p, n))
        return result

    # Array case — iterate over groups
    omega_list    = omega.ravel()
    n_protons_arr = n_protons.ravel()

    for gi, item in enumerate(omega_list):
        n = int(float(n_protons_arr[gi]))             if gi < len(n_protons_arr) else 1
        item = np.asarray(item, dtype=object)
        if item.ndim == 0:
            p = float(item)
            if p < 90:
                result.append((p, n))
        else:
            # Sub-group: all spins share the same proton count
            for sub in item.ravel():
                p = float(sub)
                if p < 90:
                    result.append((p, n))

    return result


def _flatten_omega(omega):
    """Recursively flatten nested Omega arrays to a flat list of floats."""
    result = []
    omega  = np.asarray(omega, dtype=object)
    if omega.ndim == 0:
        v = float(omega)
        if v < 90:
            result.append(v)
    else:
        for item in omega.ravel():
            item = np.asarray(item, dtype=object)
            if item.ndim == 0:
                v = float(item)
                if v < 90:
                    result.append(v)
            else:
                result.extend(_flatten_omega(item))
    return result


def _get_lit_ppms(name, table):
    """
    Get literature (ppm, n_protons) pairs.
    Falls back to hardcoded table if metabList entry not found.
    """
    for key in [name, name.lower(), name.upper(),
                name.replace('391','').replace('393',''),
                name.replace('No','')]:
        if key in table:
            entry = table[key]
            # Handle both old format (list of floats) and new (list of tuples)
            if entry and isinstance(entry[0], (int, float)):
                return [(p, 1) for p in entry]
            return entry
    # Fallback to hardcoded table
    for key in [name, name.lower()]:
        if key in METAB_PPM_TABLE:
            return [(p, 1) for p in METAB_PPM_TABLE[key]]
    return []


def _peaks_at_ppms(fid, sf, sw, lit_lines, ppm_calib=PPM_CALIB):
    """
    Get spectral amplitude and phase at specific PPM positions.
    Uses the literature PPM as the peak location, reads amplitude
    and phase from the FFT spectrum at that position.
    """
    fid   = np.asarray(fid).ravel()
    n     = len(fid)
    freq  = np.fft.fftshift(np.fft.fftfreq(n, d=1.0 / sw))
    ppm   = ppm_calib + freq / sf
    spec  = np.fft.fftshift(np.fft.fft(fid))

    # Extract just ppm values for phase correction reference
    lit_ppms_only = [p for p, n in lit_lines] if lit_lines else []

    # Find best zero-order phase by maximising the real part of the
    # dominant peak — search over 720 steps for fine resolution
    if lit_ppms_only:
        ref_ppm  = lit_ppms_only[0]
        win_mask = (ppm >= ref_ppm - 0.4) & (ppm <= ref_ppm + 0.4)
    else:
        win_mask = (ppm >= 0.5) & (ppm <= 4.3)

    best_phi, best_max = 0.0, -np.inf
    for phi_deg in np.linspace(0, 360, 721):
        sr = np.real(spec * np.exp(1j * np.radians(phi_deg)))
        mx = np.max(sr[win_mask]) if np.any(win_mask) else 0.0
        if mx > best_max:
            best_max = mx
            best_phi = np.radians(phi_deg)

    spec_phased = spec * np.exp(1j * best_phi)

    # Use lit_ppms if available, else detect peaks
    if lit_lines:
        peaks = []

        for p, n_prot in lit_lines:
            if p >= 90:
                continue

            idx = np.argmin(np.abs(ppm - p))

            # Sample real spectrum in a narrow window around lit ppm
            # Use absolute value to handle negative peaks (dispersive components)
            half_bins = max(2, int(round(n / (sw / sf) * 0.05)))
            i0 = max(0, idx - half_bins)
            i1 = min(n - 1, idx + half_bins)

            win_r   = np.real(spec_phased[i0:i1+1])
            win_c   = spec_phased[i0:i1+1]
            win_ppm = ppm[i0:i1+1]

            if len(win_r) > 0:
                # Use abs to find peak even if phasing isn't perfect at this position
                best_idx = np.argmax(np.abs(win_r))
                area     = float(win_r[best_idx])
                ph       = float(np.degrees(np.angle(win_c[best_idx])))
                ppm_loc  = float(win_ppm[best_idx])
            else:
                area    = float(np.real(spec_phased[idx]))
                ph      = float(np.degrees(np.angle(spec_phased[idx])))
                ppm_loc = float(ppm[idx])

            # Keep the peak even if area is small — preserve TE modulation.
            # Shifts that land on a peak already listed are that peak again:
            # VeSPA adds up every line, so it is listed once, with its area once.
            if any(p_in == round(ppm_loc, 6) for p_in, _, _ in peaks):
                continue
            peaks.append((round(ppm_loc, 6), area, round(ph, 6)))

        # Normalize to dominant peak = 1.0 using absolute areas
        if peaks:
            max_abs = max(abs(a) for _, a, _ in peaks)
            if max_abs > 0:
                peaks = [(p, round(a / max_abs, 12), ph) for p, a, ph in peaks]

        # Sort by absolute area descending, keep positive peaks
        peaks.sort(key=lambda x: -abs(x[1]))
        peaks = [(p, a, ph) for p, a, ph in peaks if a > 0]
        return peaks

    # Fallback: FFT peak detection
    return _fft_peaks(ppm, spec_phased, spec_phased)


def _fft_peaks(ppm, spec_phased, spec_complex,
               ppm_min=0.5, ppm_max=4.3, min_frac=0.05, max_peaks=20):
    """Fallback: detect peaks from FFT."""
    mask   = (ppm >= ppm_min) & (ppm <= ppm_max)
    ppm_m  = ppm[mask]
    spec_m = np.real(spec_phased[mask])
    sc_m   = spec_complex[mask]

    if len(spec_m) == 0 or np.max(spec_m) <= 0:
        return []

    threshold = np.max(spec_m) * min_frac
    peaks = []

    for i in range(1, len(spec_m) - 1):
        if spec_m[i] > spec_m[i-1] and \
           spec_m[i] > spec_m[i+1] and \
           spec_m[i] > threshold:
            y0, y1, y2 = spec_m[i-1], spec_m[i], spec_m[i+1]
            denom  = 2.0 * (y0 - 2.0 * y1 + y2)
            delta  = (y0 - y2) / denom if denom != 0 else 0.0
            dp     = abs(ppm_m[1] - ppm_m[0]) if len(ppm_m) > 1 else 0.0
            p_loc  = float(ppm_m[i] - delta * dp)
            area   = float(y1)
            phase  = float(np.degrees(np.angle(sc_m[i])))
            peaks.append((round(p_loc, 6), round(area, 12), round(phase, 6)))

    peaks.sort(key=lambda x: -x[1])
    return peaks[:max_peaks]


############## Command line #############

if __name__ == '__main__':
    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    if len(sys.argv) < 3:
        print("Usage: python write_vespa.py <input_folder> <output.xml> [te_ms] [metabList.mat]")
        sys.exit(1)

    input_path     = sys.argv[1]
    outpath        = sys.argv[2]
    seqte          = float(sys.argv[3]) if len(sys.argv) > 3 else 26.0
    metablist_path = sys.argv[4] if len(sys.argv) > 4 else None

    from readers.read_marss import read_marss_file, read_marss_folder

    if os.path.isdir(input_path):
        basis_list = read_marss_folder(input_path)
    else:
        basis_list = [read_marss_file(input_path)]

    print(f"Loaded {len(basis_list)} metabolites\n")
    write_vespa(basis_list, outpath, seqte=seqte, metablist_path=metablist_path)

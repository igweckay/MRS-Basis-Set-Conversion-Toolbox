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

So the lines are fitted as VeSPA rebuilds them (writers/lines.py, fit_lines): each
metabolite's FID as undamped lines times the set's own lineshape, measured on its singlets
(a 1 Hz Lorentzian in the jMRUI, LCModel TE 26, Osprey and MRSCloud samples), judged by how
well that rebuilds the FID. The lineshape is not in the lines: VeSPA's fitted lineshape is the
data's whole lineshape, not a broadening of the basis'. Lines fitted with their own dampings
(HLSVD) and written without them were another shape: the 0.512 s MARSS-derived sets hold
strongly coupled multiplets as components 2-13 Hz wide, and VeSPA drew them sharp (glucose
78 % off). A function with no line of the set's
lineshape (the samples' macromolecules and lipids, broad humps; MM09's 10 Hz singlet) is left
out and named in the prior's comment, which points to VeSPA's own macromolecule model. The
areas keep the source's scale for every metabolite (no normalisation), so fitted areas compare
between metabolites as the source's FIDs do. seqte is the echo time the set states (or the one
given to core.write).

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
            <area>29.61939201</area> ← its amplitude at t = 0, the source's scale
            <phase>3.216372</phase>  ← -(its phase), degrees (jMRUI sample's NAA CH3)
          </line>
        </prior_metabolite>
      </prior>
    </vespa_export>

vespa_lines fits the set once and gives these lines; write_gava and write_midas write the same
lines, as VeSPA-Simulation writes its simulation's lines to all three formats, and write_pyamares
and write_oxsa tie them as AMARES peaks (with the humps' own components, which an AMARES peak can
hold). The literature table is used here for spins only.
"""

import os
from collections import Counter
import numpy as np
from datetime import datetime
import xml.etree.ElementTree as ET
from xml.dom import minidom
import scipy.io as sio

if __name__ == '__main__':                    # run as a script: find the writers package
    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from writers.centre import CENTRE, at_centre
from writers.lines import (TOL, decompose, fit_lines, fwhm, is_hump, unfilled,
                           set_lineshape)


############## Literature chemical shift table (from metabList.mat) #############
# One entry per spin of metabList's spin system (a methyl is one). write_vespa counts them
# for <spins>.
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


############## Main entry point #############

def write_vespa(basis_list, outpath,
                seqte=26.0,
                source='basis_set_converter',
                metablist_path=None,
                progress=None):
    """
    Write a VeSPA Analysis prior XML file.

    Parameters
    ----------
    basis_list     : list of core struct dicts
    outpath        : str, full path for output .xml file
    seqte          : float, echo time in ms
    source         : str, source description
    metablist_path : str or None, path to metabList.mat (its numberOfProtons give <spins>)
    progress       : callable or None, progress(text, sub) once a metabolite's lines are fitted
                     or it is left out (the GUI's step list: a large set takes a while)
    """
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

    found = vespa_lines(basis_list, progress)
    listed, unlisted = [], []
    for met_name, rows, _ in found['metabolites']:
        lit = _get_lit_ppms(met_name, ppm_table)
        spins = int(round(sum(n for _, n in lit))) if lit else len(rows)
        (listed if lit else unlisted).append(met_name)

        met_elem = ET.SubElement(prior, 'prior_metabolite', version='1.0.0')
        ET.SubElement(met_elem, 'name').text  = met_name
        ET.SubElement(met_elem, 'dim').text   = '0'
        ET.SubElement(met_elem, 'dim').text   = '0'
        ET.SubElement(met_elem, 'dim').text   = '0'
        ET.SubElement(met_elem, 'spins').text = str(spins)

        for ppm, area, phase in rows:
            line = ET.SubElement(met_elem, 'line')
            ET.SubElement(line, 'group')
            ET.SubElement(line, 'ppm').text   = f'{ppm:.8f}'
            ET.SubElement(line, 'area').text  = f'{area:.10g}'
            ET.SubElement(line, 'phase').text = f'{phase:.6f}'

    humps, short = found['humps'], found['short']
    comment.text = (
        f'Prior information from basis set with {len(basis_list)} metabolites. '
        f'Lines: each FID fitted as undamped lines times the set\'s lineshape, '
        f'{found["lineshape"]}, as VeSPA rebuilds a basis FID from its lines and applies one '
        f'lineshape in the fit: the fewest lines found that rebuild the FID to a relative '
        f'residual of {TOL:g}; ppm, amplitude at t = 0 and phase of each. The lineshape is not '
        f'in the lines: VeSPA\'s fitted lineshape is the data\'s whole lineshape. Areas in the '
        f'source\'s scale, the same for every metabolite (factor 1). '
        f'Phases in VeSPA\'s sign (the conjugate FID). '
        f'Spins: {counted}' + (f' for {", ".join(listed)}' if unlisted and listed else '') +
        (f'; the number of lines for {", ".join(unlisted)}' if unlisted else '') + '.' +
        (f' Rebuilt only to: {", ".join(short)}.' if short else '') +
        (f' No lines found for: {", ".join(found["skipped"])}.' if found['skipped'] else '') +
        (f' Left out, as broad humps that hold no line of the set\'s lineshape (the prior\'s '
         f'lines cannot hold them, and sharp lines in their place are wrong): {", ".join(humps)}. '
         f'Model them with VeSPA\'s own macromolecule model (Voigt tab, macromolecule model '
         f'"Single Basis Function - from Dataset"; analysis/block_fit_voigt.py:103, '
         f'analysis/constants.py:548).' if humps else ''))

    # Write XML
    xml_str = minidom.parseString(
        ET.tostring(root, encoding='utf-8')
    ).toprettyxml(indent='  ', encoding='utf-8')

    with open(outpath, 'wb') as f:
        f.write(xml_str)

    print(f"\n  Done. VeSPA prior XML: {os.path.basename(outpath)}")
    print(f"  {len(basis_list)} metabolites  TE={seqte} ms")

def one_name_each(basis_list):
    """Refuse a set naming two functions alike: VeSPA's importer keeps its metabolites by name
    (mrs_prior.py), as PyAMARES and OXSA keep each function's rows, so the second one would take
    the first one's place without a word."""
    twice = sorted(n for n, k in Counter(str(c.get('name', 'unknown')) for c in basis_list).items() if k > 1)
    if twice:
        raise RuntimeError(f"The priors name each function's lines after it, and more than one function "
                           f"is named {', '.join(twice)}: give each its own name first.")


def vespa_lines(basis_list, progress=None, tol=TOL, hump=None):
    """
    The lines a VeSPA prior holds for each function of a set, as VeSPA-Simulation writes them to
    its Analysis prior and to GAVA and MIDAS alike (dialog_mixed_metabolite_output.py: all three
    take the simulation's same ppms, areas and phases): ppm absolute, area the amplitude at
    t = 0 on the source's scale (factor 1), phase in degrees in VeSPA's sign, each rounded as
    the prior writes it. The set's lineshape is measured once and every function fitted once,
    to the relative residual `tol` (TOL for VeSPA, GAVA and MIDAS).

    Returns a dict:
        'metabolites' : [(name, [(ppm, area, phase), ...], residual)], the set's order
        'humps'       : names left out as broad humps (no line of the set's lineshape)
        'short'       : 'name (n lines, residual)' of those whose lines (at most MAX_LINES) rebuild
                        them only to that residual
        'skipped'     : names with no lines
        'lineshape'   : the set's lineshape in words
        'shape'       : the set's lineshape (lines.set_lineshape)
        'components'  : each hump's HLSVD components (lines.decompose), by name
    progress(text, sub), if given, is called before the functions are decomposed, with the
    lineshape, and once a function. A writer
    whose format can hold a hump's own widths (PyAMARES) passes hump(name, components): what to
    report for it, or None where it is left out after all.
    """
    one_name_each(basis_list)
    basis_list = at_centre(basis_list)     # the lines' ppm are read at CENTRE
    points = [np.asarray(core['fid']).size for core in basis_list]     # the spectrum's, for is_hump
    basis_list = [dict(core, fid=unfilled(core['fid'])) for core in basis_list]
    names = [core.get('name', 'unknown') for core in basis_list]
    if progress:                           # minutes for a set of long FIDs, before any line
        progress("Decomposing each function", f"{len(names)} function{'s' * (len(names) > 1)}, "
                                               f"then each one's lines are fitted")
    decomposed = [decompose(np.asarray(core['fid']).ravel(), float(core.get('sw') or 4000.0))
                  for core in basis_list]
    shape = set_lineshape(basis_list, decomposed)
    described = _describe(shape)
    print(f"  Lineshape: {described}")
    if progress:
        progress("Lineshape of the set's lines", described)

    found = {'metabolites': [], 'humps': [], 'short': [], 'skipped': [], 'lineshape': described,
             'shape': shape, 'components': {}}
    for count, (met_name, core, parts) in enumerate(zip(names, basis_list, decomposed), 1):
        done = f"{count}/{len(names)}: {met_name}"
        fid = np.asarray(core['fid']).ravel()
        sf  = float(core.get('sf') or 123.26)
        sw  = float(core.get('sw') or 4000.0)

        if is_hump(parts, shape, sw / points[count - 1]):
            # no line of the set's lineshape: written as lines, VeSPA would draw sharp ones
            found['humps'].append(met_name)
            found['components'][met_name] = parts
            kept = hump(met_name, parts) if hump else None
            if not kept:
                print(f"  WARNING: {met_name} is a broad hump, not lines of the set's lineshape — left out")
            if progress:
                progress(f"Lines {done}", kept or
                         "left out: a broad hump, not lines of the set's lineshape")
            continue
        lines = fit_lines(fid, sw, shape, tol, parts)
        amp, freq = lines['amplitude'], lines['freq']
        if not len(amp):
            print(f"  WARNING: no lines for {met_name} — skipping")
            found['skipped'].append(met_name)
            if progress:
                progress(f"Lines {done}", "no lines: skipped")
            continue
        if lines['residual'] > tol:
            print(f"  WARNING: {len(amp)} lines rebuild {met_name} to {lines['residual']:.1e} only")
            found['short'].append(f"{met_name} ({len(amp)} lines, {lines['residual']:.0e})")

        # VeSPA's FID is BasiC's conjugate: the line at +f Hz is at CENTRE + f / sf ppm with
        # phase -angle(c); rounded as the prior writes them, so every format holds one number
        rows = [(round(CENTRE + f / sf, 8), float(f'{abs(c):.10g}'),
                 round(-float(np.degrees(np.angle(c))), 6)) for c, f in zip(amp, freq)]
        found['metabolites'].append((met_name, rows, lines['residual']))

        print(f"  {met_name:15s}  {len(amp):3d} lines  "
              f"largest={rows[0][0]:.4f} ppm  residual={lines['residual']:.1e}")
        if progress:
            progress(f"Lines {done}", f"{len(amp)} line{'s' * (len(amp) > 1)}, "
                                      f"residual {lines['residual']:.1e}")
    return found


def _describe(shape):
    """The set's lineshape in words, for the prior's comment."""
    widths = {'Lorentzian': f"{shape['lorentz']:.3f} Hz", 'Gaussian': f"{shape['gauss']:.3f} Hz",
              'Voigt': f"{shape['lorentz']:.3f} Hz Lorentzian, {shape['gauss']:.3f} Hz Gaussian"}
    text = f"{shape['kind']} {widths[shape['kind']]}" + \
        (f" (full width {fwhm(shape):.3f} Hz)" if shape['kind'] == 'Voigt' else '')
    if shape['source'] == 'lines':
        return text + ', no singlet found: the width of its narrowest lines'
    singlets = ', '.join(shape['singlets'])
    return text + (f', measured on {singlets}' + ('' if shape['uniform'] else
                   ', which differ by more than 2 %'))


############## Helpers: spins from metabList #############

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


############## Command line #############

if __name__ == '__main__':
    import sys
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

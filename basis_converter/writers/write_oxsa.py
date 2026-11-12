"""
writers/write_oxsa.py
Basis Set Converter — OXSA writer

Writes an OXSA AMARES prior knowledge file in OXSA's own form (main/+AMARES/+priorKnowledge/
PK_3T_Cardiac.m): a MATLAB function PK_<set>.m returning
AMARES.priorKnowledge.preparePriorKnowledge(fields, values), one row a peak, whose peaks rebuild
the basis set. Read against github.com/oxsatoolbox/oxsa at 7cf1fdd.

OXSA's model (makeModelFid.m:16-28) is, for each peak,

    amplitude * exp(1i * phase * pi / 180) * exp(-|t| * linewidth * pi)
              * exp(-2 * pi^2 * t^2 * sigma^2) * exp(1i * 2 pi * t * chemShift * imagingFrequency),

summed over the peaks, with t = k / bandwidth + beginTime (makeModelFidAndJacobianReIm.m:26) and
imagingFrequency in MHz: amplitude at t = 0, phase in degrees, linewidth the Lorentzian FWHM in
Hz, sigma the Gaussian's standard deviation in Hz, chemShift in ppm from the frequency at 0 Hz of
the FID handed to AMARES.amaresFit, its carrier. A basis' line c exp(2 pi i f t) (writers/
centre.py: f = sf (ppm - centre)) is therefore a peak at chemShift = ppm - 4.65, amplitude |c|
and phase +angle(c): OXSA runs a FID as the basis does (OXSA's own test data are
amplitude exp(+2 pi i v0 t), test/fixedPeakSeparation), not conjugated as VeSPA does
(write_vespa). The set's lineshape, the Voigt exp(-pi fl t - (pi fg t)^2 / (4 ln 2))
(lines.lineshape), is OXSA's exactly with linewidth = fl and sigma = fg / (2 sqrt(2 ln 2)).

The peaks are the VeSPA prior's lines (write_vespa.vespa_lines, as for GAVA, MIDAS and
PyAMARES): each FID as undamped lines times the set's lineshape. The model is PyAMARES' (write_
pyamares) as far as OXSA's constraints can express it (createModelConstraints.m, read with
initializeOptimization.m, checkConstraints.m and applyModelConstraints.m):

    amplitude  each metabolite's lines in a ratio to its largest: the first row of a metabolite
               is its largest line, amplitude free in [0, inf]; every other row is a multiplet of
               one component (multiplet = 1) in that row's G_amplitude group with amplitudeRatio
               its ratio, which createModelConstraints makes x(first) * ratio (108-128). Only
               G_amplitude ties amplitudes across rows, and only as equal (131); the one-component
               multiplet is what turns the equality into a ratio.
    phase      one phase for the set: G_phase = 1 for every row (equal, 134-149), and RelPhase the
               line's phase minus row 1's in degrees (x(1) + RelPhase, 139-142). Row 1's phase has
               the bounds [0, 360], which OXSA fits unconstrained (initializePriorKnowledge.m:
               106-114).
    linewidth  one lineshape: G_linewidth = 1 for every row (94-106); a hump's components have
               their own widths as base_linewidth, an offset in Hz from row 1's (x(1) + b, 99-101),
               so a broader line broadens them by the same Lorentzian. A Voigt set's line rows
               share row 1's sigma (G_sigma, 152-159; each needs sigma bounds of its own, 157);
               a hump's sigma is fixed at 0.
    chemShift  fixed at the basis' (empty bounds: the initial value, 46-47). OXSA cannot tie one
               peak's shift to another's: a peak's chemShift is its own parameter or a constant
               (70-92; there is no G_chemShift, and G_chemShiftDelta ties only multiplet
               spacings, chemShift + (k - 1) * delta within one row). So neither a metabolite's
               lines nor the set can move together, and freeing any shift would move one line
               alone. The data must sit where the basis does: align them first, or let
               AMARES.amares add one shift to every chemShift ('fixOffset', amares.m:187-195).

Row 1 is the largest line of the largest singlet the lineshape was measured on (write_
pyamares.reference: sI for the jMRUI sample), so G_phase and G_linewidth groups start there, as
checkConstraints requires (44-59). A metabolite's rows stand together (initializeOptimization
starts a new parameter whenever the group number changes, 41-46). OXSA's evenly spaced
multiplet rows (one ratio, one phase for every component) would hold none of the lines: a
metabolite's lines have phases of their own. Humps (lines.is_hump) whose HLSVD components rebuild
them to the criterion are written as those components, as for PyAMARES; others are left out and
named.

AMARES.amaresFit replaces the free amplitudes' and row 1's phase's initial values by a linear
least-squares fit with one column per row, phases 0, and the phase group's mean phase
(initializePriorKnowledge.m:34-104): it knows nothing of the ratios and RelPhase, so for hundreds
of tied rows it starts far off (the jMRUI sample: rank 250 of 892 columns, amplitudes up to 1e6
times too large) and needs more than its default 100 iterations (113 there), which the usage in
the file's header raises.
Names are MATLAB field names (amaresDriver and sortDataByPeaks use them as such): letters,
digits and '_', a letter first; a metabolite's other rows are <name>_2, <name>_3, ...

Entry points:
    pk_name(basis_list)                                    -> 'PK_<set>'
    write_oxsa(basis_list, outpath, seqte, tol, progress)  -> PK_<set>.m
"""

import os
import re
from datetime import datetime

import numpy as np

if __name__ == '__main__':                    # run as a script: find the writers package
    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from writers.centre import CENTRE, at_centre
from writers.lines import TOL
from writers.write_vespa import vespa_lines
from writers.write_pyamares import reference


LW_ROOM = 80.0      # row 1's linewidth bound: from 0 to the basis' own + this (Hz)
NAME_MAX = 63       # MATLAB's namelengthmax: a function's and a field's longest name
SIGMA = 2 * np.sqrt(2 * np.log(2))      # a Gaussian's FWHM / its standard deviation

FIELDS_BOUNDS = ['peakName', 'chemShift', 'linewidth', 'amplitude', 'phase', 'chemShiftDelta',
                 'amplitudeRatio']
FIELDS_IV = ['peakName', 'chemShift', 'linewidth', 'amplitude', 'phase']
FIELDS_PK = ['peakName', 'multiplet', 'chemShiftDelta', 'amplitudeRatio', 'G_linewidth',
             'G_amplitude', 'G_phase', 'RelPhase', 'G_chemShiftDelta', 'refPeak', 'base_linewidth']


def pk_name(basis_list):
    """
    The prior's function (and file) name, PK_<set>: the file the set came from when it is one
    file (an LCModel .basis), else the folder its files are in; a MATLAB name (letters, digits,
    '_', at most 63 characters). PK_basis for a set not all of whose functions name a source in
    one folder.
    """
    sources = [b.get('source') for b in basis_list]
    folders = {os.path.dirname(os.path.abspath(f)) for f in sources if f}
    if not all(sources) or len(folders) != 1:
        stem = 'basis'
    elif len(set(sources)) == 1 and len(basis_list) > 1:
        stem = os.path.splitext(os.path.basename(sources[0]))[0]
    else:
        stem = os.path.basename(folders.pop())
    stem = re.sub(r'_+', '_', re.sub(r'[^A-Za-z0-9_]', '_', stem)).strip('_') or 'basis'
    return ('PK_' + stem)[:NAME_MAX]


############## Main entry point #############

def write_oxsa(basis_list, outpath, seqte=None, tol=TOL, progress=None):
    """
    Write an OXSA prior knowledge function, PK_<set>.m (the function is named after the file).

    Parameters
    ----------
    basis_list : list of core struct dicts
    outpath    : str, full path for the output .m file
    seqte      : float or None, echo time in ms (stated in the comment)
    tol        : float, the relative residual each metabolite's peaks rebuild its FID to
    progress   : callable or None, progress(text, sub) once a metabolite's lines are fitted
                 (write_vespa.vespa_lines)
    """
    outpath = os.path.abspath(outpath)
    os.makedirs(os.path.dirname(outpath) or '.', exist_ok=True)
    function = os.path.splitext(os.path.basename(outpath))[0]

    if not basis_list:
        raise RuntimeError("basis_list is empty")

    def hump(name, parts):
        if parts['residual'] <= tol and len(parts['amplitude']):
            return (f"{len(parts['amplitude'])} damped component"
                    f"{'s' * (len(parts['amplitude']) > 1)}, residual {parts['residual']:.1e}")
        return None
    found = vespa_lines(basis_list, progress, tol, hump)
    shape = found['shape']
    basis_list = at_centre(basis_list)
    sf = {core.get('name', 'unknown'): float(core.get('sf') or 123.26) for core in basis_list}
    mhz = next(iter(sf.values()))
    sw = float(basis_list[0].get('sw') or 4000.0)

    # each metabolite's peaks, largest first: (amplitude, ppm, linewidth in Hz or None for the
    # set's lineshape, phase in degrees), rounded as the lines are (write_vespa.vespa_lines)
    rows = dict((name, r) for name, r, _ in found['metabolites'])
    peaks, kept, left = {}, [], []
    for name in sf:
        if name in rows:        # VeSPA's phase is -angle(c); OXSA's is +angle(c)
            peaks[name] = [(area, ppm, None, -phase) for ppm, area, phase in rows[name]]
        elif name in found['humps'] and hump(name, found['components'][name]):
            parts = found['components'][name]
            peaks[name] = [(float(f'{abs(c):.10g}'), round(CENTRE + f / sf[name], 8), d / np.pi,
                            round(float(np.degrees(np.angle(c))), 6))
                           for c, f, d in zip(parts['amplitude'], parts['freq'], parts['damping'])]
            kept.append(name)
        elif name in found['humps']:
            left.append(name)
    if not peaks:
        raise RuntimeError("No metabolite has lines OXSA can hold")

    names = oxsa_names(list(sf))
    ref = reference({n: [(a, p, None if lw is None else lw * np.pi, ph) for a, p, lw, ph in v]
                     for n, v in peaks.items()}, shape)
    order = [ref] + [n for n in sf if n in peaks and n != ref]

    voigt = shape['gauss'] > 0
    lw0, sigma0 = shape['lorentz'], shape['gauss'] / SIGMA
    _, _, lw_r, phase_r = peaks[ref][0]
    lw_r = lw0 if lw_r is None else lw_r
    phase_iv = phase_r % 360.0
    phase_iv = 0.0 if phase_iv >= 360.0 else phase_iv

    table = []          # (bounds row, IV row, PK row), MATLAB text of each cell
    spans = []          # (name, first row, last row, total amplitude)
    sigma_group = None
    taken = {names[n] for n in order}
    for name in order:
        lines = peaks[name]
        first = len(table) + 1
        amp1 = lines[0][0]
        for k, (amp, ppm, lw, phase) in enumerate(lines, 1):
            row = len(table) + 1
            peak = names[name] if k == 1 else f'{names[name]}_{k}'
            i = 0
            while k > 1 and peak in taken:  # another metabolite named <name>_<k>
                peak = f'{names[name]}_{k}{_letters(i)}'
                i += 1
            taken.add(peak)
            own_lw = lw0 if lw is None else lw
            rel = _wrap(phase - phase_r)
            is_line = lw is None
            if voigt and is_line and sigma_group is None:
                sigma_group = row
            bounds = [_str(peak), '[]',
                      f'[0,{_num(own_lw + LW_ROOM)}]' if row == 1 else '[]',
                      '[0,inf]' if k == 1 else '[]',
                      '[0,360]' if row == 1 else '[]', '[]', '[]']
            iv = [_str(peak), f'{ppm - CENTRE:.8f}', _num(own_lw),
                  _num(amp1), _num(phase_iv if row == 1 else phase_iv + rel)]
            pk = [_str(peak),
                  '[]' if k == 1 else '1', '[]',
                  '[]' if k == 1 else _num(amp / amp1),
                  '1', str(first), '1',
                  '[]' if row == 1 else _num(rel),
                  '[]', '1' if row == 1 else '0',
                  '[]' if row == 1 or own_lw - lw_r < 1e-12 else _num(own_lw - lw_r)]  # never below 0
            if voigt:
                bounds.append('[0,inf]' if is_line else '[]')
                iv.append(_num(sigma0) if is_line else '0')
                pk.append(str(sigma_group) if is_line else '[]')
            table.append((bounds, iv, pk))
        spans.append((names[name], first, len(table), sum(p[0] for p in lines)))

    fields_bounds = FIELDS_BOUNDS + (['sigma'] if voigt else [])
    fields_iv = FIELDS_IV + (['sigma'] if voigt else [])
    fields_pk = FIELDS_PK + (['G_sigma'] if voigt else [])
    free = len(order) + 2 + (sigma_group is not None)
    ref_name = names[ref]
    renamed = [f'{names[n]} = {n}' for n in order if names[n] != n]
    naa = next((s for s in spans if s[0] == 'NAA'), spans[0])
    total_text = '; '.join(f'{n} {a}' + (f'-{b}' if b > a else '') + f' {_num(t)}'
                           for n, a, b, t in spans)

    head = [
        f'{function}  OXSA AMARES prior knowledge of a basis set',
        f'{len(order)} of {len(basis_list)} metabolites'
        f'{f", TE {float(seqte):g} ms" if seqte else ""}, written by basis_set_converter on '
        f'{datetime.now().strftime("%Y-%m-%d")}.',
        '',
        f'Usage, with OXSA on the MATLAB path (its startup.m) and fid the data, a column sampled '
        f'at {sw:g} Hz and {mhz:.10g} MHz with its carrier (0 Hz) at {CENTRE:g} ppm:',
        '',
        f'    pk = {function}();',
        '    n = numel(fid);',
        f"    exptParams = struct('samples', n, 'imagingFrequency', {mhz:.10g}, "
        f"'dwellTime', 1/{sw:g}, ...",
        f"        'beginTime', 0, 'timeAxis', (0:n-1).'/{sw:g}, "
        f"'ppmAxis', ((0:n-1).' - floor(n/2))*{sw:g}/n/{mhz:.10g});",
        "    [fitResults, fitStatus] = AMARES.amaresFit(double(fid), exptParams, pk, false, "
        "'MaxIter', 1000, 'MaxFunEvals', 5000);",
        f'    amount{naa[0]} = sum(fitResults.amplitude({naa[1]}:{naa[2]})) / {_num(naa[3])};'
        f'   % rows and totals below',
        '',
        'MaxIter: amaresFit starts each row from its own linear least-squares fit and the phase '
        'from the rows\' mean (initializePriorKnowledge.m:34-104), which know nothing of amplitudeRatio '
        'or RelPhase, so with many tied rows it starts far off and may need more than its default '
        '100 iterations. Data whose first point is not the FID at t = 0 (FID-A\'s sets hold it '
        'halved): fit fid(2:end) with samples n - 1, beginTime 1/' + f'{sw:g}' + ' and timeAxis '
        'starting there.',
        '',
        'Data shifted by s ppm from the basis set: first add s to every chemShift, '
        'for k = 1:numel(pk.initialValues), pk.initialValues(k).chemShift = '
        'pk.initialValues(k).chemShift + s; end (as AMARES.amares\'s fixOffset does).',
        '',
        f'Units: chemShift in ppm relative to the carrier, {CENTRE:g} ppm (a line at 2.008 ppm '
        f'is at {2.008 - CENTRE:.3f}); linewidth in Hz (Lorentzian FWHM); phase in degrees; '
        f'amplitude at t = 0, in the source\'s scale for every metabolite (factor 1)' +
        (f'; sigma in Hz (the Gaussian\'s standard deviation, its FWHM / {SIGMA:.4f}).' if voigt
         else '.'),
        '',
        f"Peaks: each metabolite's lines, found as for the VeSPA prior: undamped lines times the "
        f"set's lineshape, {found['lineshape']}, the fewest that rebuild its FID to a relative "
        f'residual of {tol:g} from the second point. One row a line; a metabolite\'s first row is '
        f'its largest line.',
        '',
        'Ties, as far as OXSA can express them (createModelConstraints.m):',
        "  amplitude: only each metabolite's first row is free, [0,inf]; its other rows are "
        'one-component multiplets (multiplet = 1) in its G_amplitude group, amplitude = the first '
        "row's times amplitudeRatio, the basis' ratio.",
        f'  phase: one for the set, row 1\'s ({ref_name}, the largest line of the largest '
        f'singlet{f", {ref}" if ref != ref_name else ""}); every row is in G_phase 1 with '
        f'RelPhase its phase minus row 1\'s (degrees).',
        f'  linewidth: one for the set, row 1\'s, in G_linewidth 1 (0 to {_num(lw_r + LW_ROOM)} Hz)'
        + (', and one sigma (G_sigma)' if voigt else '') +
        ('; a hump\'s components have their own widths as base_linewidth, Hz added to row 1\'s.'
         if kept else '.'),
        '  chemShift: fixed at the basis\' (empty bounds). OXSA ties no peak\'s shift to '
        'another\'s (a shift is free or fixed; G_chemShiftDelta ties only multiplet spacings), so '
        'freeing one would move one line alone: align the data to the basis set, as above.',
        f'Free: {free} parameters (each metabolite\'s amplitude, row 1\'s linewidth and phase'
        + (' and sigma' if sigma_group else '') + f'); {len(table)} peaks.',
        '',
        f"Each metabolite's rows and total amplitude (the sum of its fitted amplitudes over that "
        f'total is the amount of the basis function): {total_text}.',
    ]
    if kept:
        head.append(f'Broad humps as their own damped components (HLSVD): '
                    f'{", ".join(names[n] for n in kept)}.')
    if left or found['skipped']:
        why = ([f"{', '.join(left)}: broad humps that neither lines of the set's lineshape nor "
                f"HLSVD components rebuild to {tol:g}"] if left else []) + \
              ([f"{', '.join(found['skipped'])}: no lines found"] if found['skipped'] else [])
        head.append('Left out (approximated by nothing): ' + '; '.join(why) +
                    ". Leave them to the fit's baseline.")
    if found['short']:
        head.append(f'Rebuilt only to: {", ".join(found["short"])}.')
    if renamed:
        head.append('Names (MATLAB field names): ' + ', '.join(renamed) + '.')

    out = [f'function pk = {function}()', f'% {head[0]}']           # MATLAB's H1 line, whole
    for line in head[1:]:
        out += ['%'] if not line else ['% ' + w for w in _wrap_text(line, 96)]
    out.append('')

    def cells(title, array, fields, part):
        return [f'fields.{title} = {{' + ', '.join(_str(f) for f in fields) + '};',
                f'values.{array} = {{...'] + \
               ['    ' + ', '.join(r[part]) + ';' for r in table] + ['};', '']
    out += ['%% Bounds'] + cells('Bounds', 'boundsCellArray', fields_bounds, 0)
    out += ['%% Initial values'] + cells('IV', 'IVCellArray', fields_iv, 1)
    out += ['%% Prior knowledge'] + cells('PK', 'PKCellArray', fields_pk, 2)
    out += ['pk = AMARES.priorKnowledge.preparePriorKnowledge(fields, values);', 'end', '']
    text = '\n'.join(out)

    with open(outpath, 'w', encoding='utf-8') as f:
        f.write(text)

    print(f"\n  Done. OXSA prior knowledge: {os.path.basename(outpath)}")
    print(f"  {len(order)} metabolites, {len(table)} peaks, {free} free parameters")


def oxsa_names(names):
    """
    Peak names OXSA can use as field names, for the set's names in order: any character not a
    letter, digit or '_' made '_', an 'x' before a name not starting with a letter (as OXSA's
    x2_3_DPG), short enough for a row's '_<k>' to follow within 63 characters; names then the
    same get _a, _b, ... A dict source name -> OXSA name.
    """
    out = {}
    for n in names:
        stem = re.sub(r'[^A-Za-z0-9_]', '_', n) or 'peak'
        stem = ('x' + stem if not stem[0].isalpha() else stem)[:NAME_MAX - 6]
        name, i = stem, 0
        while name in out.values():
            name = f'{stem}_{_letters(i)}'
            i += 1
        out[n] = name
    return out


def _letters(i):
    """a, b, ..., z, aa, ab, ..."""
    text = ''
    i += 1
    while i:
        i, r = divmod(i - 1, 26)
        text = chr(ord('a') + r) + text
    return text


def _num(x):
    """A number as MATLAB reads it, 10 significant digits."""
    text = f'{float(x):.10g}'
    return '0' if text in ('-0', '0') else text


def _str(s):
    return "'" + s.replace("'", "''") + "'"


def _wrap(phase):
    """An angle in degrees within [-180, 180), to 6 decimals."""
    return round((phase + 180.0) % 360.0 - 180.0, 6)


def _wrap_text(text, width):
    """A comment paragraph in lines of at most `width`, continuation lines indented 2 more than
    its first; a line of code (indented 4) is kept whole."""
    indent = len(text) - len(text.lstrip(' '))
    if indent >= 4:
        return [text]
    lines, line = [], ' ' * indent
    for w in text.split():
        if line.strip() and len(line) + 1 + len(w) > width:
            lines.append(line)
            line = ' ' * (indent + 2) + w
        else:
            line = f'{line} {w}' if line.strip() else line + w
    return lines + [line]


############## Command line ##############

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 3:
        print("Usage: python write_oxsa.py <input_folder> <output folder> [te_ms]")
        sys.exit(1)

    from readers.read_marss import read_marss_file, read_marss_folder
    input_path = sys.argv[1]
    seqte      = float(sys.argv[3]) if len(sys.argv) > 3 else None

    basis_list = read_marss_folder(input_path) if os.path.isdir(input_path) \
                 else [read_marss_file(input_path)]
    print(f"Loaded {len(basis_list)} metabolites\n")
    write_oxsa(basis_list, os.path.join(sys.argv[2], pk_name(basis_list) + '.m'), seqte=seqte)

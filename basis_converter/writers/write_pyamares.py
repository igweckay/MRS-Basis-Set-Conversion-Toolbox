"""
writers/write_pyamares.py
Basis Set Converter — PyAMARES writer

Writes a PyAMARES prior knowledge CSV whose peaks rebuild the basis set, each metabolite as
one peak and its sub-peaks tied to it, as PyAMARES writes a multiplet (docs/source/simple.rst,
"Peak Name Suffix" and "Constraints for Multiplets"; kernel/PriorKnowledge.py 0.3.28).

PyAMARES' model (kernel/fid.py:76-81, Vanhamme et al. 1997 eq. 6) is, for each peak,

    ak * exp(1j * phi) * exp(-dk * (1 - g + g * t) * t) * exp(1j * 2 pi * fk * t),

summed over the peaks (multieq6), with t = deadtime + k / sw (PriorKnowledge.py:494) and the
data multiplied by exp(1j * 2 pi * carrier * MHz * t) (502). The CSV gives amplitude = ak,
chemicalshift in ppm (fk = ppm * MHz, unitconverter 176), linewidth in Hz (dk = pi * lw, 179),
phase in degrees (182) and g. So with carrier = the set's centre and deadtime = 0, a line a
basis holds as c * exp(2 pi i f t) is a peak at ppm = centre + f / sf, amplitude |c| and phase
+angle(c): PyAMARES' FID runs as the basis' own, not as VeSPA's conjugate (write_vespa).

The lines are the VeSPA prior's (write_vespa.vespa_lines): each FID fitted as undamped lines
times the set's lineshape, measured on its singlets, which is what tied AMARES peaks are when
every sub-peak's linewidth and g are the base peak's. For each metabolite:

    the base peak (a digit-free name): its largest line, free: amplitude (|c|, the source's
        scale, the same for every metabolite), chemicalshift (absolute ppm), linewidth and g
        (the set's lineshape: dk(1 - g) = pi * Lorentzian width and dk * g =
        (pi * Gaussian width)**2 / (4 ln 2) make AMARES' damping the measured Voigt exactly),
        phase (degrees), with bounds;
    sub-peaks <base>2, <base>3, ... in decreasing amplitude, tied to the base by expressions:
        amplitude <base>*<ratio>, chemicalshift <base>+<offset>Hz, linewidth <base>,
        phase <base>+<offset in radians>, g <base>.

PyAMARES evaluates an expression in its internal units (extract_expr, 96-157: freq in Hz, dk,
phase in radians) and reads only an expression's first number with its unit (refindall, 16-25):
a ppm offset loses its '+', so offsets are in Hz; and scientific notation would split a number,
so every number is written in positional notation. An expression names the base, which stands
to its left (390-396). result_sum adds up every peak whose name is the same once digits are
removed (util/report.py:115-129), so each metabolite is one row; its amplitude divided by the
metabolite's total (the sum of its lines' amplitudes, stated in the comment) is the amount of
the basis function in the data.

Names: PyAMARES allows digits only as a multiplet's suffix (assert_peak_format, 243-263), and
a peak name is part of lmfit's parameter names (ak_<name>, ...): digits are removed and any other
character not a letter becomes '_'; names that are then the same get _a, _b, ... (MM09 in a set
holding MM). The comment gives each changed name's source.

Humps (macromolecules, lipids): they hold no line of the set's lineshape (lines.is_hump), but an
AMARES peak has its own width. A hump whose HLSVD components (lines.decompose) rebuild it to the
same criterion is written as its own metabolite: its largest component the base, the others tied
with linewidth as an offset from the base's (in dk, pi * Hz). One they do not rebuild (the jMRUI
sample's MM) is left out and named in the comment.

The criterion is the VeSPA prior's (TOL), so the peaks rebuild each function as PyAMARES models
it, to 1e-4 from the second point: hundreds of tied peaks for coupled spins. A PyAMARES fit of
such a prior is slow, but not for the fit's sake (a few minutes for the jMRUI sample): after
it, lmfit propagates the uncertainties to every tied parameter, evaluating all of them once per
free one (Parameters.create_uvars, asteval_with_uncertainties), a time that grows with the
square of the tied parameters, and fitAMARES always does it. For the jMRUI sample that is hours
at 1e-4 and still about an hour at 1e-2, a rebuild 100 times worse; no criterion that keeps the
rebuild honest makes it short, so the criterion stays, and the comment says what to expect.

CSV layout (PyAMARES reads rows iloc[1:6] as initial values and iloc[7:] as bounds, 81, 204):
    # comment lines (find_header_row skips them, 266-280)
    Index, <peaks>
    Initial Values
    amplitude, chemicalshift, linewidth, phase, g
    Bounds
    amplitude, chemicalshift, linewidth, phase, g

Entry point:
    write_pyamares(basis_list, outpath, seqte, tol, progress) -> .csv file
"""

import os
import re
import csv
from collections import Counter

import numpy as np

from writers.centre import CENTRE, at_centre
from writers.lines import TOL, MAX_LINES
from writers.write_vespa import vespa_lines


PPM_ROOM = 0.08     # the base peak's chemical shift bound: its ppm +- this
LW_ROOM = 80.0      # its linewidth bound: from 0 to the basis' own + this (Hz)


############## Main entry point #############

def write_pyamares(basis_list, outpath, seqte=None, tol=TOL, progress=None):
    """
    Write a PyAMARES prior knowledge CSV.

    Parameters
    ----------
    basis_list : list of core struct dicts
    outpath    : str, full path for output .csv file
    seqte      : float or None, echo time in ms (stated in the comment)
    tol        : float, the relative residual each metabolite's peaks rebuild its FID to
    progress   : callable or None, progress(text, sub) once a metabolite's lines are fitted
                 (write_vespa.vespa_lines)
    """
    outpath = os.path.abspath(outpath)
    os.makedirs(os.path.dirname(outpath) or '.', exist_ok=True)

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

    # the set's lineshape as AMARES' damping: dk (1 - g) t + dk g t**2
    dk_lorentz = np.pi * shape['lorentz']
    dk_gauss = (np.pi * shape['gauss']) ** 2 / (4 * np.log(2))
    lw0 = (dk_lorentz + dk_gauss) / np.pi
    g0 = dk_gauss / (dk_lorentz + dk_gauss) if dk_lorentz + dk_gauss else 0.0

    # each metabolite's peaks: (amplitude, ppm, damping dk or None for the set's, phase in rad)
    rows = dict((name, r) for name, r, _ in found['metabolites'])
    peaks, kept, left = {}, [], []
    for name in sf:
        if name in rows:
            peaks[name] = [(area, ppm, None, -np.radians(phase)) for ppm, area, phase in rows[name]]
        elif name in found['humps'] and hump(name, found['components'][name]):
            parts = found['components'][name]
            # rounded as the lines are (write_vespa.vespa_lines)
            peaks[name] = [(float(f'{abs(c):.10g}'), round(CENTRE + f / sf[name], 8), d,
                            np.radians(round(float(np.degrees(np.angle(c))), 6)))
                           for c, f, d in zip(parts['amplitude'], parts['freq'], parts['damping'])]
            kept.append(name)
        elif name in found['humps']:
            left.append(name)
    # named over the whole set, so that a name means the same function whatever is left out
    names = {n: base for n, base in pyamares_names(list(sf)).items() if n in peaks}

    columns = []        # (column name, [5 initial values], [5 bounds])
    for name, base in names.items():
        lines = peaks[name]
        amp1, ppm1, d1, phi1 = lines[0]
        lw1 = lw0 if d1 is None else d1 / np.pi
        columns.append((base,
                        [_num(amp1), _num(ppm1), _num(lw1), _num(np.degrees(phi1)),
                         _num(g0 if d1 is None else 0.0)],
                        ['(0,', f'({_num(ppm1 - PPM_ROOM)},{_num(ppm1 + PPM_ROOM)})',
                         f'(0,{_num(lw1 + LW_ROOM)})',
                         f'({_num(np.degrees(phi1) - 180)},{_num(np.degrees(phi1) + 180)})',
                         '(0,1)']))
        for k, (amp, ppm, d, phi) in enumerate(lines[1:], 2):
            dphi = (phi - phi1 + np.pi) % (2 * np.pi) - np.pi
            columns.append((f'{base}{k}',
                            [f'{base}*{_num(amp / amp1)}',
                             f'{base}{_signed((ppm - ppm1) * sf[name])}Hz',
                             base if d is None else f'{base}{_signed(d - d1)}',
                             f'{base}{_signed(dphi)}',
                             base],
                            ['', '', '(0,', '', '']))

    totals = '; '.join(f'{names[n]} {_num(sum(p[0] for p in peaks[n]))}' for n in names)
    renamed = [f'{base} = {name}' for name, base in names.items() if base != name]
    comment = [
        f'PyAMARES prior knowledge by basis_set_converter, {len(names)} of {len(basis_list)} '
        f'metabolites{f", TE {float(seqte):g} ms" if seqte else ""}. Fit with MHz={mhz:.10g}, '
        f'carrier={CENTRE:g} and deadtime=0 (initialize_FID(fid, priorknowledgefile, '
        f'MHz={mhz:.10g}, sw=..., deadtime=0, carrier={CENTRE:g})): chemical shifts are absolute '
        f'ppm about a centre of {CENTRE:g}, amplitudes at t = 0. Total amplitude of each '
        f"metabolite's peaks (result_sum's amplitude divided by it is the amount of the basis "
        f'function): {totals}.',
        f"Peaks: each metabolite's lines, found as for the VeSPA prior: undamped lines times the "
        f"set's lineshape, {found['lineshape']}, the fewest found that rebuild its FID to a "
        f'relative residual of {tol:g} from the second point (FID-A\'s sets hold the first '
        f'halved: fit data whose first point is not its value at t = 0 with '
        f'truncate_initial_points=1, or 2 where PyAMARES 0.3.28 makes the time axis a point too '
        f'long, as for 8191 points at 4000 Hz). The first column is its largest '
        f'line; <name>2, <name>3, ... are tied to it: amplitude as a ratio, chemical shift as an '
        f'offset in Hz, phase as an offset in radians, linewidth and g equal to its. Amplitudes in '
        f"the source's scale for every metabolite (factor 1).",
        f'Linewidth and g start at the set\'s lineshape as AMARES damping dk((1-g)t+gt^2), '
        f'dk = pi * linewidth: linewidth {_num(lw0)} Hz, g {_num(g0)}' +
        (' (fit with g_global=False to use it).' if g0 else '.'),
    ]
    if kept:
        comment.append(f'Broad humps as their own damped components (HLSVD), each with its own '
                       f'width: linewidth tied as an offset from the first component\'s in dk '
                       f'(pi * Hz): {", ".join(names[n] for n in kept)}.')
    if left or found['skipped']:
        why = ([f"{', '.join(left)}: broad humps that neither lines of the set's lineshape nor "
                f"HLSVD components rebuild to {tol:g}"] if left else []) + \
              ([f"{', '.join(found['skipped'])}: no lines found"] if found['skipped'] else [])
        comment.append('Left out: ' + '; '.join(why) + ". Leave them to the fit's baseline "
                       '(PyAMARES: truncate_initial_points).')
    if found['short']:
        comment.append(f'Rebuilt by {MAX_LINES} lines only to: {", ".join(found["short"])}.')
    comment.append(f'{len(columns)} peaks: {4 * len(names)} free parameters (g fixed, as by '
                   f'default) and {5 * (len(columns) - len(names))} tied ones. After the fit, '
                   f'fitAMARES (lmfit) propagates the uncertainties to each tied parameter by '
                   f'evaluating all of them once per free one: its time grows with the square of '
                   f'the tied parameters, hours for thousands.')
    if renamed:
        comment.append('Names (PyAMARES allows digits only as a multiplet\'s suffix): ' +
                       ', '.join(renamed) + '.')

    table = [['Index'] + [c[0] for c in columns], ['Initial Values'] + [''] * len(columns)]
    for i, label in enumerate(('amplitude', 'chemicalshift', 'linewidth', 'phase', 'g')):
        table.append([label] + [c[1][i] for c in columns])
    table.append(['Bounds'] + [''] * len(columns))
    for i, label in enumerate(('amplitude', 'chemicalshift', 'linewidth', 'phase', 'g')):
        table.append([label] + [c[2][i] for c in columns])

    with open(outpath, 'w', newline='') as f:
        for line in comment:
            f.write('# ' + line.replace('"', '').replace('\n', ' ') + '\n')
        csv.writer(f, lineterminator='\n').writerows(table)

    print(f"\n  Done. PyAMARES CSV written: {os.path.basename(outpath)}")
    print(f"  {len(names)} metabolites, {len(columns)} peaks")


def pyamares_names(names):
    """
    Peak names PyAMARES keeps apart, for the set's names in order: digits removed (PyAMARES reads
    them as a multiplet's suffix, and result_sum adds up the peaks whose names are then the
    same), any other character not a letter made '_'. Where names are then the same, one that
    was the plain name keeps it and the others get _a, _b, ... in the set's order (MM, MM09 ->
    MM, MM_a). A dict source name -> PyAMARES name.
    """
    stems = {n: re.sub(r'_+', '_', re.sub(r'[^A-Za-z_]', '_', re.sub(r'\d', '', n))).strip('_')
             or 'Peak' for n in names}
    count = Counter(stems.values())
    out = {n: stems[n] for n in names if count[stems[n]] == 1 or n == stems[n]}
    used = set(out.values())
    for n in names:
        if n not in out:
            i = 0
            while f'{stems[n]}_{_letters(i)}' in used:
                i += 1
            out[n] = f'{stems[n]}_{_letters(i)}'
            used.add(out[n])
    return {n: out[n] for n in names}


def _letters(i):
    """a, b, ..., z, aa, ab, ..."""
    text = ''
    i += 1
    while i:
        i, r = divmod(i - 1, 26)
        text = chr(ord('a') + r) + text
    return text


def _num(x):
    """A number as PyAMARES reads it in a cell or an expression: positional, 10 digits."""
    text = np.format_float_positional(float(x), precision=10, unique=False, fractional=False,
                                      trim='-')
    return '0' if text in ('-0', '') else text


def _signed(x):
    return ('-' if x < 0 else '+') + _num(abs(x))


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

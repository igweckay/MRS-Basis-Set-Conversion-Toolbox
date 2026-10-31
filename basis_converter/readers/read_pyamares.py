"""
readers/read_pyamares.py
Basis Set Converter — PyAMARES prior knowledge CSV reader

Reads a PyAMARES prior knowledge CSV file into a peak params dict.
This dict can then be passed to write_oxsa() to produce an OXSA
MATLAB struct, or used for inspection/display.

PyAMARES CSV format (from generate_pyamares_pk_from_basis.py):

    Index, Asc, Asp, Cho, CrCHtwo, ...       ← metabolite names
    Initial Values, , , , ...
    amplitude, 0.028, 0.040, 0.270, ...
    chemicalshift, 3.742, 2.734, 3.198, ...
    linewidth, 6.0, 6.0, 5.0, ...
    phase, 0.0, 0.0, 0.0, ...
    g, 0.0, 0.0, 0.0, ...
    Bounds, , , , ...
    amplitude, (0,, (0,, (0,, ...             ← open upper bound
    chemicalshift, (3.66,3.82), (2.65,2.81), ...
    linewidth, (1,80), (1,80), ...
    phase, (-180,180), (-180,180), ...
    g, (0,1), (0,1), ...

BasiC's own CSV (writers/write_pyamares.py) starts with '#' comment lines and ties peaks by
expressions, as PyAMARES allows (kernel/PriorKnowledge.py 0.3.28): '#' starts a comment, as
in PyAMARES' pandas.read_csv(comment='#'), and a cell holding letters is an expression, which
is read as PyAMARES reads it. Only its first number counts with its unit (refindall: 'ppm'
times MHz, 'Hz' dropped); the name it starts with is the same parameter of a peak to its left,
in PyAMARES' internal units (Hz, dk = pi * Hz, radians): NAA*0.5, NAA+61.6Hz, NAA, NAA+0.16.
Such a cell is returned evaluated, in the CSV's units (ppm, Hz, degrees), and its text kept
under 'expressions'. MHz is the `mhz` given, else the 'MHz=' a comment states; without either
a chemical shift tied in Hz is NaN.

Peak params dict returned:
    {
        'names'         : list of str
        'amplitude'     : list of float  — initial values
        'chemicalshift' : list of float  — initial values in ppm
        'linewidth'     : list of float  — initial values in Hz
        'phase'         : list of float  — initial values in degrees
        'g'             : list of float  — 0=Lorentzian, 1=Gaussian
        'amp_bounds'    : list of tuple  — (lower, upper) or (lower, None)
        'cs_bounds'     : list of tuple  — (lower_ppm, upper_ppm)
        'lw_bounds'     : list of tuple  — (lower_hz, upper_hz)
        'phase_bounds'  : list of tuple  — (lower_deg, upper_deg)
        'g_bounds'      : list of tuple  — (0, 1) always
        'expressions'   : dict row -> list of str or None, each tied cell's text
        'mhz'           : float or None, the MHz the chemical shifts were evaluated with
        'comment'       : list of str, the comment lines, '#' removed
        'source'        : str, path to .csv file
    }

Entry point:
    read_pyamares(path, mhz=None) -> peak params dict
"""

import os
import re
import ast
import csv
import operator
import numpy as np

ROWS = ('amplitude', 'chemicalshift', 'linewidth', 'phase', 'g')


def read_pyamares(path, mhz=None):
    """
    Read a PyAMARES prior knowledge CSV file.

    Parameters
    ----------
    path : str
        Full path to a PyAMARES .csv file.
    mhz  : float or None
        The MHz PyAMARES would be given, for chemical shifts tied in Hz; by default the
        'MHz=' a comment line states.

    Returns
    -------
    peak params dict
    """
    path = os.path.abspath(path)

    if not os.path.isfile(path):
        raise RuntimeError(f"File not found: {path}")

    with open(path, 'r', newline='') as f:
        text = f.read().splitlines()
    comment = [line.split('#', 1)[1].strip() for line in text if '#' in line]
    rows = [r for r in csv.reader(line.split('#', 1)[0] for line in text) if any(c.strip() for c in r)]
    if mhz is None:
        stated = re.search(r'MHz\s*=\s*([0-9]+(?:\.[0-9]+)?)', ' '.join(comment))
        mhz = float(stated.group(1)) if stated else None

    if not rows:
        raise RuntimeError(f"Empty CSV file: {os.path.basename(path)}")

    ############## Parse header row — metabolite names #############
    header = rows[0]
    if header[0].strip().lower() != 'index':
        raise RuntimeError(
            f"Expected 'Index' as first cell in {os.path.basename(path)}, "
            f"got '{header[0]}'"
        )
    names = [n.strip() for n in header[1:] if n.strip()]
    n_metabs = len(names)

    ############## Separate initial values and bounds sections #############
    init_rows   = {}
    bounds_rows = {}
    section     = None

    for row in rows[1:]:
        if not row:
            continue
        label = row[0].strip().lower()

        if label == 'initial values':
            section = 'init'
            continue
        if label == 'bounds':
            section = 'bounds'
            continue

        if section == 'init' and label in ROWS:
            init_rows[label] = [v.strip() for v in row[1:n_metabs + 1]]

        if section == 'bounds' and label in ('amplitude', 'chemicalshift',
                                              'linewidth', 'phase', 'g'):
            vals = row[1:n_metabs + 1]
            bounds_rows[label] = [_parse_bound(v) for v in vals]

    ############## Validate required initial value rows are present #############
    required = ['amplitude', 'chemicalshift', 'linewidth', 'phase', 'g']
    missing  = [r for r in required if r not in init_rows]
    if missing:
        raise RuntimeError(
            f"Missing initial value rows in {os.path.basename(path)}: "
            f"{', '.join(missing)}"
        )
    expressions = {label: [v if _is_expression(v) else None for v in init_rows[label]]
                   for label in ROWS}
    for label in ROWS:
        init_rows[label] = _evaluate_row(label, init_rows[label], names, mhz,
                                         os.path.basename(path))

    ############## Build result dict #############
    result = {
        'names'         : names,
        'amplitude'     : init_rows['amplitude'],
        'chemicalshift' : init_rows['chemicalshift'],
        'linewidth'     : init_rows['linewidth'],
        'phase'         : init_rows['phase'],
        'g'             : init_rows['g'],
        'amp_bounds'    : bounds_rows.get('amplitude',
                            [(0.0, None)] * n_metabs),
        'cs_bounds'     : bounds_rows.get('chemicalshift',
                            [(cs - 0.08, cs + 0.08)
                             for cs in init_rows['chemicalshift']]),
        'lw_bounds'     : bounds_rows.get('linewidth',
                            [(1.0, 80.0)] * n_metabs),
        'phase_bounds'  : bounds_rows.get('phase',
                            [(-180.0, 180.0)] * n_metabs),
        'g_bounds'      : bounds_rows.get('g',
                            [(0.0, 1.0)] * n_metabs),
        'expressions'   : expressions,
        'mhz'           : mhz,
        'comment'       : comment,
        'source'        : path,
    }

    # Suppress verbose output — GUI displays this instead
    return result


############## Internal helpers #############

_OPS = {ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul,
        ast.Div: operator.truediv, ast.Pow: operator.pow, ast.USub: operator.neg,
        ast.UAdd: operator.pos}
# a number cell in the CSV's units -> PyAMARES' internal units, and back (unitconverter)
_TO = {'chemicalshift': lambda v, mhz: v * mhz, 'linewidth': lambda v, mhz: v * np.pi,
       'phase': lambda v, mhz: np.radians(v)}
_FROM = {'chemicalshift': lambda v, mhz: v / mhz, 'linewidth': lambda v, mhz: v / np.pi,
         'phase': lambda v, mhz: np.degrees(v)}


def _evaluate_row(label, cells, names, mhz, filename):
    """
    One row of initial values as floats in the CSV's units, each expression evaluated as
    PyAMARES evaluates it: in the internal units, its first name the same parameter of a peak
    to its left. A chemical shift row is evaluated at 1 MHz where none is known (its Hz-tied
    cells are then NaN).
    """
    scale = mhz if mhz else 1.0
    to = _TO.get(label, lambda v, m: v)
    back = _FROM.get(label, lambda v, m: v)
    internal, out = {}, []
    for name, cell in zip(names, cells):
        if not _is_expression(cell):
            value = _parse_float(cell)
            internal[name] = to(value, scale)
            out.append(value)
            continue
        expr = cell.replace(' ', '')
        match = re.search(r'([-+]?\d+(\.\d+)?)(Hz|ppm)?', expr)        # refindall
        hz = bool(match and match.group(3) == 'Hz')
        if match and match.group(3) == 'ppm':
            expr = expr.replace(match.group(1) + 'ppm', str(float(match.group(1)) * scale))
        elif hz:
            expr = expr.replace('Hz', '')
        try:
            value = _expression(ast.parse(expr, mode='eval').body, internal)
        except (SyntaxError, KeyError, TypeError, ZeroDivisionError) as ex:
            raise RuntimeError(f"{filename}: {label} of {name} is '{cell}', which PyAMARES "
                               f"cannot evaluate from the peaks to its left ({ex})")
        internal[name] = value
        out.append(float('nan') if label == 'chemicalshift' and hz and not mhz
                   else float(back(value, scale)))
    return out


def _is_expression(cell):
    """A cell PyAMARES reads as an expression: text that is not a number (pd.to_numeric)."""
    try:
        float(cell)
        return False
    except ValueError:
        return bool(cell.strip())


def _expression(node, values):
    """An arithmetic expression of numbers and names (the values of peaks to the left)."""
    if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
        return float(node.value)
    if isinstance(node, ast.Name):
        return values[node.id]
    if isinstance(node, ast.BinOp) and type(node.op) in _OPS:
        return _OPS[type(node.op)](_expression(node.left, values), _expression(node.right, values))
    if isinstance(node, ast.UnaryOp) and type(node.op) in _OPS:
        return _OPS[type(node.op)](_expression(node.operand, values))
    raise TypeError(ast.dump(node))

def _parse_float(val):
    """
    Parse a float from a CSV cell.
    Returns 0.0 if empty or unparseable.
    """
    val = str(val).strip()
    if not val:
        return 0.0
    try:
        return float(val)
    except ValueError:
        return 0.0


def _parse_bound(val):
    """
    Parse a bound string into a (lower, upper) tuple.

    Handles:
        '(0,'          → (0.0, None)    open upper bound
        '(0.0, None)'  → (0.0, None)
        '(1,80)'       → (1.0, 80.0)
        '(-180,180)'   → (-180.0, 180.0)
        '(3.66,3.82)'  → (3.66, 3.82)
        ''             → (None, None)
    """
    val = str(val).strip()
    if not val:
        return (None, None)

    # Remove surrounding parentheses
    val = val.strip('()')

    parts = [p.strip() for p in val.split(',')]

    lower = _safe_float(parts[0]) if len(parts) > 0 else None
    upper = _safe_float(parts[1]) if len(parts) > 1 else None

    return (lower, upper)


def _safe_float(val):
    """Convert to float, returning None for empty or 'None' strings."""
    val = str(val).strip()
    if not val or val.lower() == 'none':
        return None
    try:
        return float(val)
    except ValueError:
        return None


############## Quick test #############

if __name__ == '__main__':
    import sys

    if len(sys.argv) < 2:
        print("Usage: python read_pyamares.py <path_to_csv>")
        sys.exit(1)

    result = read_pyamares(sys.argv[1])

    print(f"\nNames : {result['names']}")
    print(f"CS    : {result['chemicalshift'][:3]} ...")
    print(f"Amp   : {result['amplitude'][:3]} ...")
    print(f"LW    : {result['linewidth'][:3]} ...")
    print(f"\nCS bounds sample: {result['cs_bounds'][:3]}")
    print(f"LW bounds sample: {result['lw_bounds'][:3]}")

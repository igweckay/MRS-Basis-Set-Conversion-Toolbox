"""
writers/lines.py
Basis Set Converter — a basis function as lines

Two decompositions of a FID, both in the basis' own convention (a line at f Hz is at
centre + f / sf ppm, writers/centre.py), both fitted from the second point on (formats disagree
about the first: MRSCloud's FID-A sets hold it halved, so the spectrum's baseline is right) with
the lines' amplitudes at t = 0 extrapolated from the points after it.

decompose: decaying complex exponentials, each with its own damping,

    fid[k] = sum_j c_j z_j**k,    z_j = exp((-d_j + 2 pi i f_j) / sw),

by HLSVD, as VeSPA's common/hlsvdpropy.py does it (Barkhuijsen et al., J Magn Reson 1987;
Laudadio et al., J Magn Reson 2002): the SVD of the FID's Hankel matrix, the signal poles z_j
from the shift invariance of its first K left singular vectors, and the amplitudes c_j by
least squares. The Hankel matrix is built from the whole FID (half its points a row), and
only its leading MAX_ORDER + 25 singular vectors are computed, by a randomized range finder
with two power iterations whose products with the matrix are FFT convolutions (Halko et al.,
SIAM Review 2011), where HLSVDPRO uses Lanczos: a full SVD of an 8192-point FID's
4096 x 4097 matrix takes minutes, and poles from its first 2048 points alone do not rebuild
MRSCloud's coupled spins. The start vectors are seeded: a FID gives the same lines each time.
c_j is the line's complex amplitude at t = 0 (|c_j| its spectral integral, angle(c_j) its
phase), f_j its frequency in Hz and d_j its damping in 1/s. Its criterion: the fewest lines
that rebuild the FID to a relative residual ||fid - lines|| / ||fid|| of at most TOL. The model
order is the smallest K whose K poles meet it (searched by bisection), up to MAX_ORDER; then
the components contributing least are dropped while the rest still meet it.

fit_lines: undamped lines times one lineshape L(t), the set's,

    fid(t) = L(t) sum_j a_j exp(2 pi i f_j t),

which is what VeSPA rebuilds from a prior (each line undamped, one lineshape applied by the fit
to them all). The lineshape is measured on the set's
singlets (set_lineshape: Lorentzian, Gaussian or Voigt); the criterion is decompose's, with
the same TOL, so the lines are judged by the rebuild VeSPA makes. A function with no line of
the set's lineshape (is_hump: a macromolecule or lipid hump) has no such lines.

Why TOL = 1e-4: the text formats' precision floors the residual (the LCModel sample
SLASER_Siemens26.basis at about 3e-6, jMRUI's .txt at 1e-7); at 1e-5 the margin is 3x, and at
1e-6 the LCModel functions take any number of noise lines.

Used by write_vespa.
"""

import numpy as np
import scipy.linalg as sla
from scipy.optimize import least_squares
from scipy.signal import fftconvolve

TOL = 1e-4          # largest relative residual of a rebuilt FID (from its second point)
MAX_ORDER = 300     # most poles tried


def decompose(fid, sw, tol=TOL, max_order=MAX_ORDER):
    """
    The lines of `fid` (sampled at `sw` Hz): a dict of arrays 'amplitude' (complex, at t = 0),
    'freq' (Hz), 'damping' (1/s), sorted by |amplitude|, largest first, and 'residual', the
    relative residual of their rebuild from the second point on, and 'order', the poles
    before any was dropped.
    """
    fid = np.asarray(fid, dtype=complex).ravel()
    y = fid[1:]                                    # the first point is not fitted
    if not np.any(y):
        return _result(np.zeros(0, complex), np.zeros(0, complex), 0.0, 0, sw)
    norm = np.linalg.norm(y)

    u = _signal_space(y, max_order + 25)
    max_order = max(1, min(max_order, u.shape[1] - 1))

    def fit(K):
        uk = u[:, :K]
        z = sla.eigvals(sla.lstsq(uk[:-1], uk[1:])[0])
        a, r = _fit(y, z)
        return z, a, np.linalg.norm(r) / norm, r

    # the smallest order meeting TOL, by bisection between one that fails and one that meets it
    tried = {}
    def residual(K):
        if K not in tried:
            tried[K] = fit(K)
        return tried[K][2]

    hi = 1
    while hi < max_order and not residual(hi) <= tol:
        hi = min(2 * hi, max_order)
    if residual(hi) <= tol:
        lo = hi // 2                                # fails (or 0)
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if residual(mid) <= tol:
                hi = mid
            else:
                lo = mid
        K = hi
    else:                                          # none meets TOL: the best tried
        K = min(tried, key=lambda k: tried[k][2])
    z, a, res, r = tried[K]

    # drop the components contributing least while the rest still meet TOL
    if res <= tol:
        cols = _columns(z, len(y))
        energy = np.abs(a) ** 2 * np.sum(np.abs(cols) ** 2, axis=0)
        limit = (tol * norm) ** 2
        keep = np.ones(len(z), bool)
        for j in np.argsort(energy):
            trial = r + a[j] * cols[:, j]           # residual without component j
            if np.vdot(trial, trial).real <= limit:
                r, keep[j] = trial, False
        z = z[keep]
        a, r = _fit(y, z)                           # refitted without them
        res = np.linalg.norm(r) / norm
    return _result(z, _at_zero(a, z, len(y)), res, K, sw)


def _hankel_times(x, rows, cols, v):
    """H @ v for the rows x cols Hankel matrix H[i, j] = x[i + j], by FFT convolution."""
    return fftconvolve(x[:, None], v[::-1], axes=0)[cols - 1:cols - 1 + rows]


def _signal_space(x, k):
    """The leading k left singular vectors of the Hankel matrix of x (half its points a row),
    by a randomized range finder with two power iterations."""
    n = len(x)
    cols = n // 2 + 1
    rows = n - cols + 1
    k = min(k, rows, cols)
    rng = np.random.default_rng(0)
    q = np.linalg.qr(_hankel_times(x, rows, cols, rng.standard_normal((cols, k))
                                   + 1j * rng.standard_normal((cols, k))))[0]
    for _ in range(2):
        q = np.linalg.qr(_hankel_times(np.conj(x), cols, rows, q))[0]       # H^H q
        q = np.linalg.qr(_hankel_times(x, rows, cols, q))[0]                # H q
    b = _hankel_times(np.conj(x), cols, rows, q).conj().T                     # q^H H
    return q @ np.linalg.svd(b, full_matrices=False)[0]


def _columns(z, n):
    """z_j**k for k = 1 .. n, column j; a growing pole's column is scaled to end at |1|."""
    logz = np.log(z.astype(complex))
    grow = np.maximum(logz.real, 0.0)
    k = np.arange(1, n + 1)[:, None]
    with np.errstate(under='ignore'):
        return np.exp(k * logz[None, :] - (n * grow)[None, :])


def _fit(y, z):
    """The least-squares weights of the columns of poles z on y (points 1 .. n of the FID),
    and the residual."""
    if len(z) == 0:
        return np.zeros(0, complex), y.copy()
    cols = _columns(z, len(y))
    a = sla.lstsq(cols, y)[0]
    return a, y - cols @ a


def _at_zero(a, z, n):
    """The amplitudes at t = 0 of weights a of _columns(z, n)."""
    grow = np.maximum(np.log(np.abs(z)), 0.0) if len(z) else np.zeros(0)
    with np.errstate(under='ignore'):
        return a * np.exp(-n * grow)


def _result(z, c, res, order, sw):
    order_ = np.argsort(-np.abs(c))
    z, c = z[order_], c[order_]
    logz = np.log(z.astype(complex)) if len(z) else np.zeros(0, complex)
    return {'amplitude': c,
            'freq': logz.imag * sw / (2 * np.pi),
            'damping': -logz.real * sw,
            'residual': float(res),
            'order': int(order)}


def rebuild(lines, n, sw, damped=True):
    """The FID (n points) of `lines`, with their damping or without it."""
    t = np.arange(n) / float(sw)
    d = lines['damping'] if damped else np.zeros(len(lines['freq']))
    return np.exp(np.outer(t, -d + 2j * np.pi * lines['freq'])) @ lines['amplitude']


############## Undamped lines times the set's lineshape (what VeSPA rebuilds) #############

SINGLETS = {'cr', 'cr391', 'crno391', 'crch2', 'crch3', 'pcr', 'pcr393', 'pcrno393', 'pcrch2',
            'pcrch3', 'tcr', 'cre', 'gly', 'si', 'sl', 'scyllo', 'naa', 'ace'}
RHO = 0.9           # a line is added only where it correlates at most this with every other
RHO_FLOOR = 0.99    # ... or at most this, where the misfit is at a line already there
SEED_WIDTH = 0.1    # an HLSVD component within 10 % of the lineshape's width is a line as it is
HUMP_WIDTH = 2.0    # a function none of whose components is narrower than twice the set's
                    # linewidth holds no line of the set's lineshape: a hump, not lines
MAX_LINES = 1000


def lineshape(t, shape):
    """L(t) of a Voigt lineshape {'lorentz': Hz, 'gauss': Hz} (full widths at half maximum)."""
    return np.exp(-np.pi * shape['lorentz'] * t
                  - (np.pi * shape['gauss'] * t) ** 2 / (4 * np.log(2)))


def fwhm(shape):
    """The Voigt's full width at half maximum in Hz (Olivero and Longbothum 1977, to 0.02 %)."""
    fl, fg = shape['lorentz'], shape['gauss']
    return 0.5346 * fl + np.sqrt(0.2166 * fl ** 2 + fg ** 2)


def singlet_shape(fid, sw, lines=None):
    """
    The Voigt lineshape of the largest line of `fid` (a singlet: creatine's CH3, NAA's CH3, ...):
    the components of its decomposition `lines` (decompose(fid, sw)) away from it are taken
    off, and c * exp(2 pi i f t) * L(t) is fitted to the rest from the second point on.
    A dict 'lorentz', 'gauss' (Hz) and 'residual', the fit's relative residual.
    """
    fid = np.asarray(fid, dtype=complex).ravel()
    t = np.arange(len(fid)) / float(sw)
    lines = lines if lines is not None else decompose(fid, sw)
    c, f, d = lines['amplitude'], lines['freq'], lines['damping']
    far = np.abs(f - f[0]) > max(5.0, 5 * d[0] / np.pi)
    y = (fid - np.exp(np.outer(t, -d[far] + 2j * np.pi * f[far])) @ c[far])[1:]
    scale, tt = np.linalg.norm(y), t[1:]
    y, c0 = y / scale, c[0] / scale

    def misfit(p):                        # p: Re c, Im c, f, lorentz, gauss**2
        m = (p[0] + 1j * p[1]) * np.exp(2j * np.pi * p[2] * tt - np.pi * p[3] * tt
                                        - (np.pi * tt) ** 2 * p[4] / (4 * np.log(2)))
        return np.r_[(m - y).real, (m - y).imag]
    p = least_squares(misfit, [c0.real, c0.imag, f[0], max(d[0] / np.pi, 0.0), 0.0],
                      bounds=([-np.inf, -np.inf, -np.inf, 0, 0], np.inf), x_scale='jac',
                      xtol=1e-14, ftol=1e-15, gtol=1e-15).x
    return {'lorentz': float(p[3]), 'gauss': float(np.sqrt(p[4])),
            'residual': float(np.linalg.norm(misfit(p)))}


def set_lineshape(basis_list, decomposed=None):
    """
    The lineshape the set's lines share, measured on its singlets (the functions SINGLETS
    names: creatine's CH3 and CH2, NAA's CH3, glycine, scyllo-inositol, acetate) where the fit
    of a singlet's line leaves under 1 % of it: the median of their Lorentzian and Gaussian
    widths. A width under 1 % of the line's full width is noise
    of the fit and set to 0, so the shape is called Lorentzian, Gaussian or Voigt. A set with no
    singlet is taken as Lorentzian, its width the median over its functions of their
    components' amplitude-weighted median width. `decomposed`: decompose() of each function.
    A dict 'lorentz', 'gauss', 'kind', 'singlets' (name -> its own shape), 'uniform' (every
    singlet's full width within 2 % of the set's, and of its kind), 'source'.
    """
    decomposed = decomposed or [None] * len(basis_list)
    found = {}
    for core, lines in zip(basis_list, decomposed):
        name = str(core.get('name', ''))
        if name.lower() in SINGLETS and name not in found:
            shape = singlet_shape(np.asarray(core['fid']).ravel(), float(core['sw']), lines)
            if shape['residual'] <= 0.01:
                found[name] = shape
    if found:
        shape = {'lorentz': float(np.median([s['lorentz'] for s in found.values()])),
                 'gauss': float(np.median([s['gauss'] for s in found.values()])),
                 'source': 'singlets'}
    else:
        widths = []
        for core, lines in zip(basis_list, decomposed):
            lines = lines if lines is not None else decompose(core['fid'], float(core['sw']))
            w, a = lines['damping'] / np.pi, np.abs(lines['amplitude'])
            if len(w):
                o = np.argsort(w)
                widths.append(w[o][np.searchsorted(np.cumsum(a[o]) / a.sum(), 0.5)])
        shape = {'lorentz': float(np.median(widths)), 'gauss': 0.0, 'source': 'lines'}
    _snap(shape)
    for s in found.values():
        _snap(s)
    shape['singlets'] = found
    shape['uniform'] = all(abs(fwhm(s) / fwhm(shape) - 1) <= 0.02 and s['kind'] == shape['kind']
                           for s in found.values())
    return shape


def _snap(shape):
    """Set a width under 1 % of the full width to 0, and name the shape."""
    full = fwhm(shape)
    for key in ('lorentz', 'gauss'):
        if shape[key] < 0.01 * full:
            shape[key] = 0.0
    shape['kind'] = ('Voigt' if shape['lorentz'] and shape['gauss'] else
                     'Gaussian' if shape['gauss'] else 'Lorentzian')


def is_hump(lines, shape):
    """
    Whether a function, decomposed into `lines`, is a hump rather than lines of the set's
    lineshape: none of its components with at least 1 % of its largest amplitude is narrower
    than HUMP_WIDTH times the set's full width (MM09's 10 Hz, the Osprey and MRSCloud sets'
    35-40 Hz macromolecules and lipids, against 1 Hz). Undamped lines times the set's lineshape
    can rebuild such a FID only as a comb of sharp lines, which is not its spectrum.
    """
    a = np.abs(lines['amplitude'])
    if not len(a):
        return False
    narrowest = np.min(lines['damping'][a >= 0.01 * a.max()]) / np.pi
    return narrowest > HUMP_WIDTH * fwhm(shape)


def fit_lines(fid, sw, shape, tol=TOL, lines=None, max_lines=MAX_LINES):
    """
    `fid` as undamped lines times the set's lineshape L(t) (`shape`, set_lineshape),
    fid(t) ~ L(t) * sum_j a_j exp(2 pi i f_j t), which is what VeSPA rebuilds from a prior's lines
    before its fit applies one lineshape to them all. The criterion is decompose's: the fewest
    lines found with ||fid - L * sum|| / ||fid|| <= tol, from the second point on. So the
    residual is taken where the FID is, and the tail, where L is small, weighs what it weighs
    in the FID; nothing is divided by L.

    1. The components of the HLSVD decomposition `lines` (decompose(fid, sw)) whose damping
       is L's Lorentzian width to within SEED_WIDTH of L's full width are lines as they are:
       each a line at its frequency (its pole projected onto L's circle). They may lie close
       together (MRSCloud's lactate holds pairs 0.1 Hz apart).
    2. While the residual exceeds tol, lines are added at the peaks of the residual's
       correlation with L(t) exp(2 pi i f t) (a zero-filled FFT, the peak interpolated), at
       least SEP from every line, where SEP is the separation below which two lines of this
       lineshape correlate more than RHO over the FID (0.5-0.6 Hz for the samples' 1 Hz
       Lorentzian): closer, the least squares trade large opposite amplitudes for a broad
       component's shape, and the areas mean nothing (without it, the 0.512 s sets' glucose
       came out with sum |a| up to 5e7 |sum a|). The added lines' frequencies are refined by
       variable projection (Golub and Pereyra 1973; Kaufman's Jacobian), each kept SEP from
       the others. Where the largest misfit lies within SEP of a line already there, or no
       line can be placed, SEP is halved from then on, down to the separation at RHO_FLOOR.
    3. The lines contributing least are dropped, one at a time, while the rest meet tol, and
       the rest refitted.
    The amplitudes are least squares throughout. Returns a dict as decompose's, without
    'damping': 'amplitude' (at t = 0), 'freq', sorted by |amplitude|, and 'residual'.
    """
    fid = np.asarray(fid, dtype=complex).ravel()
    scale = np.linalg.norm(fid[1:])
    y = fid[1:] / scale if scale else fid[1:]
    n, sw = len(y), float(sw)
    t = np.arange(1, n + 1) / sw
    L = lineshape(t, shape)
    if not scale:
        return {'amplitude': np.zeros(0, complex), 'freq': np.zeros(0), 'residual': 0.0}

    sep, floor = _separation(t, L, RHO), _separation(t, L, RHO_FLOOR)
    lines = lines if lines is not None else decompose(fid, sw)
    seeded = np.abs(lines['damping'] / np.pi - shape['lorentz']) <= SEED_WIDTH * fwhm(shape)
    f = lines['freq'][seeded]
    free = np.zeros(len(f), bool)
    a, r = _solve(_atoms(t, L, f), y)

    N = 8 * (n + 1)
    grid = np.fft.fftfreq(N, 1 / sw)
    rounds = 0
    while np.linalg.norm(r) > tol and len(f) < max_lines:
        c = np.zeros(N, complex)
        c[1:n + 1] = L * r
        C = np.abs(np.fft.fft(c))                       # |sum_t L r exp(-2 pi i f t)|
        peaks = np.flatnonzero((C > np.roll(C, 1)) & (C >= np.roll(C, -1)) & (C >= 0.5 * C.max()))
        peaks = peaks[np.argsort(-C[peaks])]
        new = _place(grid[peaks] + _vertex(C, peaks) * sw / N, f, sep, sw,
                     max_lines - len(f), sep / 2 < floor)
        if not new:                         # None: the largest misfit is at a line; []: all are
            if sep / 2 < floor:
                break
            sep /= 2
            continue
        f, free = np.r_[f, new], np.r_[free, np.ones(len(new), bool)]
        rounds += 1
        if rounds % 3 == 0:
            f = _refine(y, t, L, f, free, sep, 10)
        a, r = _solve(_atoms(t, L, f), y)
    if free.any():
        f = _refine(y, t, L, f, free, sep, 30)
    A = _atoms(t, L, f)
    a, r = _solve(A, y)

    if np.linalg.norm(r) <= tol:                     # drop the lines contributing least
        energy = np.abs(a) ** 2 * np.sum(np.abs(A) ** 2, axis=0)
        keep = np.ones(len(f), bool)
        for j in np.argsort(energy):
            trial = r + a[j] * A[:, j]
            if np.vdot(trial, trial).real <= tol ** 2:
                r, keep[j] = trial, False
        f = f[keep]
        a, r = _solve(_atoms(t, L, f), y)
    o = np.argsort(-np.abs(a))
    return {'amplitude': a[o] * scale, 'freq': f[o], 'residual': float(np.linalg.norm(r))}


def _vertex(C, k):
    """The offset, in bins, of the parabola's vertex through C[k - 1], C[k], C[k + 1]."""
    lo, mid, hi = C[k - 1], C[k], C[(k + 1) % len(C)]
    bend = lo - 2 * mid + hi
    return np.where(bend != 0, 0.5 * (lo - hi) / np.where(bend != 0, bend, 1), 0.0)


def _place(fk, f, apart, sw, room, last):
    """Up to 20 (and `room`) of the frequencies fk, largest misfit first, each at least `apart`
    from the lines f and from each other; None if the first is not, unless this is the `last`
    separation tried."""
    new = []
    for x in fk:
        taken = np.r_[f, new]
        if len(taken) and np.min(np.abs((taken - x + sw / 2) % sw - sw / 2)) < apart:
            if not new and not last:
                return None
            continue
        new.append(x)
        if len(new) >= min(20, room):
            break
    return new


def _atoms(t, L, f):
    """L(t) exp(2 pi i f_j t), column j."""
    return L[:, None] * np.exp(2j * np.pi * np.outer(t, f))


def _solve(A, y):
    """The least-squares weights of the columns of A on y, and the residual (normal equations:
    the columns are kept apart, and their Gram matrix is fast to form)."""
    if not A.shape[1]:
        return np.zeros(0, complex), y.copy()
    G, b = A.conj().T @ A, A.conj().T @ y
    try:
        a = sla.cho_solve(sla.cho_factor(G), b)
    except np.linalg.LinAlgError:
        a = sla.lstsq(A, y)[0]
    return a, y - A @ a


def _separation(t, L, rho):
    """The smallest separation in Hz at which two lines of lineshape L correlate at most rho
    over the times t."""
    w = L ** 2
    for d in np.arange(0.005, 20.0, 0.005):
        if abs(np.sum(w * np.exp(2j * np.pi * d * t))) <= rho * w.sum():
            return d
    return 20.0


def _refine(y, t, L, f, free, sep, iterations):
    """The frequencies of the free lines moved to lower the residual (variable projection),
    each staying at least sep from every other line and within sep of where it was."""
    o = np.argsort(f)
    fs, fr = f[o], free[o]
    gap_l, gap_r = np.r_[np.inf, np.diff(fs)], np.r_[np.diff(fs), np.inf]
    # half a gap to another free line, all of it to a fixed one, less sep
    room_l = np.where(np.r_[False, fr[:-1]], (gap_l - sep) / 2, gap_l - sep)
    room_r = np.where(np.r_[fr[1:], False], (gap_r - sep) / 2, gap_r - sep)
    lo, hi = np.empty_like(f), np.empty_like(f)
    lo[o] = fs - np.clip(room_l, 0, sep)
    hi[o] = fs + np.clip(room_r, 0, sep)
    move = free & (hi - lo > 1e-9)
    if not move.any():
        return f

    def full(x):
        g = f.copy()
        g[move] = x
        return g

    def misfit(x):
        r = _solve(_atoms(t, L, full(x)), y)[1]
        return np.r_[r.real, r.imag]

    def jacobian(x):
        A = _atoms(t, L, full(x))
        G = A.conj().T @ A
        a = sla.solve(G, A.conj().T @ y, assume_a='pos')
        D = (2j * np.pi * t)[:, None] * A[:, move] * a[move]
        J = -(D - A @ sla.solve(G, A.conj().T @ D, assume_a='pos'))     # Kaufman
        return np.r_[J.real, J.imag]

    x0 = np.clip(f[move], lo[move] + 1e-12, hi[move] - 1e-12)
    return full(least_squares(misfit, x0, jac=jacobian, bounds=(lo[move], hi[move]),
                              method='trf', max_nfev=iterations).x)


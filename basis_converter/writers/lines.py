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
only its leading singular vectors are computed, OVERSAMPLE more than the order tried, by a
randomized range finder with two power iterations whose products with the matrix are FFT
convolutions (Halko et al., SIAM Review 2011), where HLSVDPRO uses Lanczos: a full SVD of an
8192-point FID's 4096 x 4097 matrix takes minutes, and poles from its first 2048 points alone
do not rebuild MRSCloud's coupled spins. They are computed for START_ORDER poles first, and
for twice as many each time the order search needs more: most functions need few (a singlet
one), and the vectors' cost grows with the square of their number. The start vectors are
seeded: a FID gives the same lines each time.
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
the set's lineshape (is_hump: a macromolecule or lipid hump) has no such lines. Its sums
(_Lines) are in closed form for a Lorentzian, and its refinement hands scipy's trust-region
step p + 1 rows instead of the 2n of the residual (_refine): the same steps, in seconds where
they took minutes for MRSCloud's 2 s, 8192-point set.

Why TOL = 1e-4: the text formats' precision floors the residual (the LCModel sample
SLASER_Siemens26.basis at about 3e-6, jMRUI's .txt at 1e-7); at 1e-5 the margin is 3x, and at
1e-6 the LCModel functions take any number of noise lines.

Used by write_vespa (and through it by write_gava, write_midas and write_pyamares).
"""

import numpy as np
import scipy.linalg as sla
from scipy.optimize import least_squares
import scipy.fft as sfft

TOL = 1e-4          # largest relative residual of a rebuilt FID (from its second point)
MAX_ORDER = 300     # most poles tried
START_ORDER = 16    # the signal space is first computed for this many poles, then doubled
OVERSAMPLE = 25     # and this many more vectors than poles


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
    hankel = _Hankel(y)
    max_order = max(1, min(max_order, hankel.rank - 1))
    space = {}

    def subspace(K):
        """The signal space's leading vectors for order K, computed for a cap of the next power
        of two (and OVERSAMPLE more): most functions need far fewer than max_order poles."""
        cap = min(max_order, 1 << max(K - 1, START_ORDER - 1).bit_length())
        if space.get('cap', 0) < K:
            u = hankel.signal_space(cap + OVERSAMPLE)
            space.update(cap=cap, shift=u[:-1].conj().T @ u[1:], last=u[-1].conj())
        return space

    def fit(K):
        # the poles: the eigenvalues of the least-squares Z with u[:-1, :K] Z = u[1:, :K]. u's
        # columns are orthonormal, so u[:-1]'s Gram matrix is I - w w^H, w its last row
        # conjugated, and Z = (I + w w^H / (1 - w^H w)) u[:-1]^H u[1:] (Sherman and Morrison)
        s = subspace(K)
        w, shift = s['last'][:K], s['shift'][:K, :K]
        z = sla.eigvals(shift + np.outer(w, w.conj() @ shift) / (1 - np.vdot(w, w).real))
        a, r = _fit(y, z, tol * norm)
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
    z = tried[K][0]
    a, r = _fit(y, z)
    res = np.linalg.norm(r) / norm

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


class _Hankel:
    """The rows x cols Hankel matrix H[i, j] = x[i + j] of x (half its points a row), whose
    products with a block of vectors are FFT convolutions, x's transforms computed once."""

    def __init__(self, x):
        self.cols = len(x) // 2 + 1
        self.rows = len(x) - self.cols + 1
        self.rank = min(self.rows, self.cols)
        self.size = sfft.next_fast_len(len(x), False)    # no product wraps onto the rows kept
        self.fx = sfft.fft(x, self.size)
        self.fxc = sfft.fft(np.conj(x), self.size)

    def times(self, v, adjoint=False):
        """H @ v, or H^H @ v."""
        fx, rows, cols = ((self.fxc, self.cols, self.rows) if adjoint else
                          (self.fx, self.rows, self.cols))
        fv = sfft.fft(v[::-1], self.size, axis=0, workers=-1)
        return sfft.ifft(fv * fx[:, None], axis=0, workers=-1)[cols - 1:cols - 1 + rows]

    def signal_space(self, k):
        """Its leading k left singular vectors, by a randomized range finder with two power
        iterations (Halko et al., SIAM Review 2011); between them, the block is kept apart by
        an LU factor, which spans what a QR factor spans at a fraction of its cost."""
        k = min(k, self.rank)
        rng = np.random.default_rng(0)
        q = self.times(rng.standard_normal((self.cols, k)) + 1j * rng.standard_normal((self.cols, k)))
        for _ in range(2):
            q = sla.lu(q, permute_l=True, check_finite=False)[0]
            q = sla.lu(self.times(q, adjoint=True), permute_l=True, check_finite=False)[0]
            q = self.times(q)
        q = np.linalg.qr(q)[0]
        b = self.times(q, adjoint=True).conj().T                              # q^H H
        return q @ np.linalg.svd(b, full_matrices=False)[0]


def _columns(z, n):
    """z_j**k for k = 1 .. n, column j; a growing pole's column is scaled to end at |1|.
    z**k on the first B points times its value at each block's start: a product a point, not
    an exponential."""
    logz = np.log(z.astype(complex))
    grow = np.maximum(logz.real, 0.0)
    size = int(np.ceil(np.sqrt(n)))
    blocks = -(-n // size)
    with np.errstate(under='ignore'):
        head = np.exp(np.arange(1, size + 1)[:, None] * logz[None, :])
        start = np.exp(np.arange(blocks)[:, None] * size * logz[None, :] - (n * grow)[None, :])
        return (start[:, None, :] * head[None, :, :]).reshape(blocks * size, len(z))[:n]


def _fit(y, z, near=None):
    """The least-squares weights of the columns of poles z on y (points 1 .. n of the FID),
    and the residual. With `near`, a residual norm: the normal equations, refined once
    (seminormal equations with a correction step, Bjorck 1987), where their residual is not
    within 1 % of it and their Gram matrix is positive definite; else, as without it, an
    SVD-based least squares, many times slower on a long FID."""
    if len(z) == 0:
        return np.zeros(0, complex), y.copy()
    cols = _columns(z, len(y))
    if near is not None:
        try:
            G = sla.cho_factor(cols.conj().T @ cols)
            a = sla.cho_solve(G, _ah(cols, y))
            r = y - cols @ a
            a = a + sla.cho_solve(G, _ah(cols, r))
            r = y - cols @ a
            if np.all(np.isfinite(r)) and abs(np.linalg.norm(r) / near - 1) > 0.01:
                return a, r
        except np.linalg.LinAlgError:
            pass
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


def unfilled(fid, floor=1e-5):
    """
    `fid` without the zero-filling some files keep: an FSL-MRS file converted from a .BASIS
    holds its 2n points, the second n near 0 (1e-6). HLSVD fits no sum of damped lines to that
    edge. The points after the last above `floor` times the largest are left off, if they are
    an eighth of the FID or more. That also takes the tail off a function that decays below
    `floor` early (a broad hump), which holds nothing at TOL; an undamped FID is kept whole.
    """
    fid = np.asarray(fid).ravel()
    above = np.nonzero(np.abs(fid) > floor * np.abs(fid).max())[0] if fid.size else []
    if not len(above):
        return fid
    end = above[-1] + 1
    return fid[:end] if len(fid) - end >= len(fid) // 8 else fid


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
    singlet is taken as Lorentzian. Its width is read off each function's narrowest components:
    the width below which a tenth of its amplitude lies. A multiplet's components are each of
    the set's lineshape, and a broad hump's are all broad; so the width is the median over the
    functions within twice the narrowest's. `decomposed`: decompose() of each function.
    A dict 'lorentz', 'gauss', 'kind', 'singlets' (name -> its own shape), 'uniform' (every
    singlet's full width within 2 % of the set's, and of its kind), 'source'.
    """
    decomposed = decomposed or [None] * len(basis_list)
    found = {}
    for core, lines in zip(basis_list, decomposed):
        name = str(core.get('name', ''))
        if name.lower() in SINGLETS and name not in found:
            fid, sw = np.asarray(core['fid']).ravel(), float(core['sw'])
            lines = lines if lines is not None else decompose(fid, sw)
            if not len(lines['amplitude']):            # no signal: no line to measure
                continue
            shape = singlet_shape(fid, sw, lines)
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
                widths.append(w[o][np.searchsorted(np.cumsum(a[o]) / a.sum(), 0.1)])
        if not widths:
            raise RuntimeError("No basis function holds any signal: there are no lines to write.")
        widths = np.array(widths)
        shape = {'lorentz': float(np.median(widths[widths <= 2 * widths.min()])), 'gauss': 0.0,
                 'source': 'lines'}
    _snap(shape)
    for s in found.values():
        _snap(s)
    shape['singlets'] = found
    shape['uniform'] = all(abs(fwhm(s) - fwhm(shape)) <= 0.02 * fwhm(shape) and s['kind'] == shape['kind']
                           for s in found.values())
    return shape


def _snap(shape):
    """Set a width under 1 % of the full width, or under a millihertz (an undamped set), to 0,
    and name the shape."""
    full = fwhm(shape)
    for key in ('lorentz', 'gauss'):
        if shape[key] < max(0.01 * full, 1e-3):
            shape[key] = 0.0
    shape['kind'] = ('Voigt' if shape['lorentz'] and shape['gauss'] else
                     'Gaussian' if shape['gauss'] else 'Lorentzian')


def is_hump(lines, shape, resolution=0.0):
    """
    Whether a function, decomposed into `lines`, is a hump rather than lines of the set's
    lineshape: none of its components with at least 1 % of its largest amplitude is narrower
    than HUMP_WIDTH times the set's full width (MM09's 10 Hz, the Osprey and MRSCloud sets'
    35-40 Hz macromolecules and lipids, against 1 Hz). Undamped lines times the set's lineshape
    can rebuild such a FID only as a comb of sharp lines, which is not its spectrum. The width
    is at least the spectrum's `resolution` (sw / n, Hz): an undamped set (0 Hz) would otherwise
    call any line with the least width a hump.
    """
    a = np.abs(lines['amplitude'])
    if not len(a):
        return False
    narrowest = np.min(lines['damping'][a >= 0.01 * a.max()]) / np.pi
    return narrowest > HUMP_WIDTH * max(fwhm(shape), resolution)


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
    model = _Lines(n, sw, shape)
    t, L = model.t, model.L
    if not scale:
        return {'amplitude': np.zeros(0, complex), 'freq': np.zeros(0), 'residual': 0.0}

    sep, floor = _separation(t, L, RHO), _separation(t, L, RHO_FLOOR)
    lines = lines if lines is not None else decompose(fid, sw)
    seeded = np.abs(lines['damping'] / np.pi - shape['lorentz']) <= SEED_WIDTH * fwhm(shape)
    f = lines['freq'][seeded]
    free = np.zeros(len(f), bool)
    a, r = model.solve(f, y)

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
            f = _refine(model, y, f, free, sep, 10)
        a, r = model.solve(f, y)
    if free.any():
        f = _refine(model, y, f, free, sep, 30)
    a, r = model.solve(f, y)
    A = model.atoms(f)

    if np.linalg.norm(r) <= tol:                     # drop the lines contributing least
        energy = np.abs(a) ** 2 * np.sum(np.abs(A) ** 2, axis=0)
        keep = np.ones(len(f), bool)
        for j in np.argsort(energy):
            trial = r + a[j] * A[:, j]
            if np.vdot(trial, trial).real <= tol ** 2:
                r, keep[j] = trial, False
        f = f[keep]
        a, r = model.solve(f, y)
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


class _Lines:
    """
    Undamped lines times one lineshape L(t) on the FID's times t = k / sw, k = 1 .. n: their
    matrix A, A @ a and A^H v, their least squares, and the sums the Gram matrix and the
    refinement's normal equations are made of,

        S_l[j, k] = sum_t t**l L(t)**2 exp(2 pi i (f_k - f_j) t),    l = 0, 1, 2.

    exp(2 pi i f t) on the first B points (`head`) times its value at each block's start
    (`start`): A costs a complex product a point, not an exponential, and A @ a and A^H v are
    matrix products of the two (O(n m) for m lines), A not formed.

    For a Lorentzian, L(t)**2 exp(2 pi i d t) = exp(k s) with s = (-2 pi lw + 2 pi i d) / sw,
    a geometric sequence: S_0 is its sum, exp(s) (exp(n s) - 1) / (exp(s) - 1), and S_1, S_2
    the sum's derivatives in s over sw, sw**2. That is O(m**2), where the products of A are
    O(n m**2): the 2 s MRSCloud FIDs hold 8192 points and their coupled spins hundreds of
    lines. Its terms do not cancel where the lineshape decays to exp(-1/2) over the FID or
    below (pi lw n / sw >= 1/2); otherwise, and for a Gaussian or a Voigt, the sums are A's
    products.
    """

    def __init__(self, n, sw, shape):
        self.n, self.sw = n, float(sw)
        self.t = np.arange(1, n + 1) / self.sw
        self.L = lineshape(self.t, shape)
        self.lorentz = float(shape['lorentz'])
        self.closed = not shape['gauss'] and np.pi * self.lorentz * n / self.sw >= 0.5
        self.size = int(np.ceil(np.sqrt(n)))
        self.blocks = -(-n // self.size)

    def factors(self, f):
        """head (B x m) and start (blocks x m): A[b B + i, j] = L head[i, j] start[b, j]."""
        f = np.asarray(f, dtype=float)
        head = np.exp(2j * np.pi * np.outer(np.arange(1, self.size + 1) / self.sw, f))
        start = np.exp(2j * np.pi * np.outer(np.arange(self.blocks) * self.size / self.sw, f))
        return head, start

    def atoms(self, f):
        """A: L(t) exp(2 pi i f_j t), column j."""
        head, start = self.factors(f)
        A = np.empty((self.blocks * self.size, head.shape[1]), complex)
        np.multiply(start[:, None, :], head[None, :, :], out=A.reshape(self.blocks, self.size, -1))
        A = A[:self.n]
        A *= self.L[:, None]
        return A

    def operator(self, f):
        """(A @ a, A^H v) as functions, and A where the sums need it (else None)."""
        if not self.closed:
            A = self.atoms(f)
            return (lambda a: A @ a), (lambda v: _ah(A, v)), A
        head, start = self.factors(f)
        n, rows = self.n, self.blocks * self.size

        def dot(a):
            return self.L * ((start * a) @ head.T).ravel()[:n]

        def hdot(v):
            w = np.zeros(rows, complex)
            w[:n] = self.L * v
            return np.einsum('ji,ij->j', start.conj().T @ w.reshape(self.blocks, self.size),
                             head.conj())
        return dot, hdot, None

    def sums(self, f, A=None, orders=1, cols=None):
        """[S_0, .., S_orders-1] of the lines f: S_0 of all against all, the others of all
        (rows) against the lines f[cols] (columns; all if None). A: the lines' matrix, which
        the sums are products of where there is no closed form."""
        f = np.asarray(f, dtype=float)
        cols = slice(None) if cols is None else cols
        if not self.closed:
            A = self.atoms(f) if A is None else A
            B = A[:, cols]
            return [A.conj().T @ A] + [A.conj().T @ ((self.t ** l)[:, None] * B)
                                       for l in range(1, orders)]
        n, sw = self.n, self.sw
        rho = np.exp(-2 * np.pi * self.lorentz / sw)
        e = rho * np.outer(np.exp(-2j * np.pi * f / sw), np.exp(2j * np.pi * f / sw))
        en = rho ** n * np.outer(np.exp(-2j * np.pi * f * n / sw), np.exp(2j * np.pi * f * n / sw))
        c, cn = e / (e - 1), en / (en - 1)                 # exp(x) / (exp(x) - 1), x = s, n s
        out = [c * (en - 1)]
        if orders > 1:
            s0, c, cn = out[0][:, cols], c[:, cols], cn[:, cols]
            h = 1 + n * cn - c                              # d/ds log S_0
            out.append(s0 * h / sw)
            if orders > 2:
                out.append(s0 * (h * h + c * (c - 1) - n * n * cn * (cn - 1)) / sw ** 2)
        return out

    def solve(self, f, y):
        """The least-squares amplitudes of the lines f on y and the residual (normal equations:
        the lines are kept apart, and their Gram matrix is fast to form)."""
        if not len(f):
            return np.zeros(0, complex), y.copy()
        dot, hdot, A = self.operator(f)
        try:
            a = sla.cho_solve(sla.cho_factor(self.sums(f, A)[0]), hdot(y))
        except np.linalg.LinAlgError:
            a = sla.lstsq(self.atoms(f) if A is None else A, y)[0]
        return a, y - dot(a)


def _separation(t, L, rho):
    """The smallest separation in Hz at which two lines of lineshape L correlate at most rho
    over the times t."""
    w = L ** 2
    for d in np.arange(0.005, 20.0, 0.005):
        if abs(np.sum(w * np.exp(2j * np.pi * d * t))) <= rho * w.sum():
            return d
    return 20.0


def _ah(A, v):
    """A^H v, without a conjugated copy of A."""
    return (v.conj() @ A).conj()


def _cholesky(G):
    """The lower Cholesky factor of the Gram matrix G. Lines of a broad lineshape are nearly
    alike, and rounding can leave G short of positive definite: then the factor of G with the
    smallest ridge, from 1e-12 of its mean diagonal up, that makes it so."""
    try:
        return sla.cholesky(G, lower=True)
    except np.linalg.LinAlgError:
        scale, eye = np.mean(np.diag(G).real), np.eye(len(G))
        for ridge in 10.0 ** np.arange(-12, -2):
            try:
                return sla.cholesky(G + ridge * scale * eye, lower=True)
            except np.linalg.LinAlgError:
                pass
        raise


def _refine(model, y, f, free, sep, iterations):
    """The frequencies of the free lines moved to lower the residual (variable projection),
    each staying at least sep from every other line and within sep of where it was.

    The step is scipy's trust-region reflective one on the residual r and its Jacobian J
    (Kaufman's, J = -(I - A G^-1 A^H) D, D the lines' derivatives in frequency times their
    amplitudes). That step depends on J only through J^H J and J^H r, which _Lines' sums give
    without the 2n x p matrix J. So the problem handed to it is p + 1 rows that keep both, and
    |r|: the Jacobian [R; 0] with R^T R = Re J^H J, and the residual [R^-T Re J^H r; rest],
    the rest of |r|. Its step is the one the full residual gives (in exact arithmetic), at a
    p x p SVD where the full one was a 2n x p one."""
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
    t, p = model.t, int(move.sum())
    last = {}

    def full(x):
        g = f.copy()
        g[move] = x
        return g

    def reduced(x):
        """The p + 1 rows (residual, Jacobian) at x; the last kept for the Jacobian's call."""
        if last.get('x') is not None and np.array_equal(last['x'], x):
            return last['rows']
        g = full(x)
        dot, hdot, A = model.operator(g)
        G, S1, S2 = model.sums(g, A, 3, move)
        C = _cholesky(G)                                            # G = C C^H
        a = sla.cho_solve((C, True), hdot(y))
        r = y - dot(a)
        am = a[move]
        W = sla.solve_triangular(C, 2j * np.pi * S1 * am, lower=True)  # C^-1 A^H D
        JJ = (4 * np.pi ** 2 * np.outer(np.conj(am), am) * S2[move] - W.conj().T @ W).real
        Jr = (2j * np.pi * np.conj(am) * hdot(t * r)[move]
              + W.conj().T @ sla.solve_triangular(C, hdot(r), lower=True)).real
        rr = np.vdot(r, r).real
        try:
            R = sla.cholesky((JJ + JJ.T) / 2)
            z = sla.solve_triangular(R, Jr, trans='T')
            if z @ z > rr * (1 + 1e-6):
                raise np.linalg.LinAlgError
        except np.linalg.LinAlgError:                           # J short of full rank
            w, V = np.linalg.eigh((JJ + JJ.T) / 2)
            keep = w > p * np.finfo(float).eps * max(w.max(), 0.0)
            R = np.sqrt(np.where(keep, w, 0.0))[:, None] * V.T
            z = np.where(keep, V.T @ Jr / np.sqrt(np.where(keep, w, 1.0)), 0.0)
        rows = (np.r_[z, np.sqrt(max(rr - z @ z, 0.0))], np.r_[R, np.zeros((1, p))])
        last.update(x=x.copy(), rows=rows)
        return rows

    x0 = np.clip(f[move], lo[move] + 1e-12, hi[move] - 1e-12)
    return full(least_squares(lambda x: reduced(x)[0], x0, jac=lambda x: reduced(x)[1],
                              bounds=(lo[move], hi[move]), method='trf',
                              max_nfev=iterations).x)

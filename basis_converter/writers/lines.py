"""
writers/lines.py
Basis Set Converter — a basis function as a sum of lines (HLSVD)

A FID is decomposed into decaying complex exponentials,

    fid[k] = sum_j c_j z_j**k,    z_j = exp((-d_j + 2 pi i f_j) / sw),

by HLSVD, as VeSPA's common/hlsvdpropy.py does it (Barkhuijsen et al., J Magn Reson 1987;
Laudadio et al., J Magn Reson 2002): the SVD of the FID's Hankel matrix, the signal poles z_j
from the shift invariance of its first K left singular vectors, and the amplitudes c_j by
least squares. The Hankel matrix is built from the whole FID (half its points a row), and
only its leading MAX_ORDER + 25 singular vectors are computed, by a randomized range finder
with two power iterations whose products with the matrix are FFT convolutions (Halko et al.,
SIAM Review 2011), where HLSVDPRO uses Lanczos: a full SVD of an 8192-point FID's
4096 x 4097 matrix takes minutes, and poles from its first 2048 points alone do not rebuild
MRSCloud's coupled spins. The start vectors are seeded: a FID gives the same lines each time. c_j is the line's complex amplitude at t = 0 (|c_j| its spectral integral,
angle(c_j) its phase), f_j its frequency in Hz and d_j its damping in 1/s, all in the basis'
own convention: a line at f Hz is at centre + f / sf ppm (writers/centre.py).

The criterion: a FID's lines are the fewest that rebuild it, from its second point on, to a
relative residual ||fid - lines|| / ||fid|| of at most TOL. The model order is the smallest K
whose K poles meet it (searched by bisection: the residual falls as K grows), up to MAX_ORDER;
then the components contributing least to the FID are dropped, one at a time, while the rest
still meet it. A FID no MAX_ORDER poles rebuild to TOL keeps the MAX_ORDER that do best, and
the residual returned says so.

The first point is left out of both the decomposition and the criterion: formats disagree
about it (MRSCloud's FID-A sets hold it halved, so the spectrum's baseline is right), and the
lines' amplitudes at t = 0 are extrapolated from the points after it.

Used by write_vespa.
"""

import numpy as np
import scipy.linalg as sla
from scipy.signal import fftconvolve

TOL = 1e-4          # largest relative residual of the rebuilt FID (from its second point)
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

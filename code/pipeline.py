"""Core clean-room pipeline for the Dale et al. cosmic-CHSH reimplementation.

Implements EXACTLY what the paper describes (SCOPE.md section 6):
 - dichotomized sky f: sphere -> [-1,1]  (mixed or unmixed)
 - E(alpha) = pair average of f(x)f(y) in 2-degree bins, 6 < alpha <= 180
 - |C| = |E(a1)+E(a2)+E(a3)-E(a4)| over REALIZABLE tetrads (minus on a4)
"""
import numpy as np
import healpy as hp

BIN_W = 2.0
EDGES = np.arange(6.0, 180.0 + BIN_W, BIN_W)   # 90 bins
NBIN = len(EDGES) - 1
CENTERS = 0.5 * (EDGES[:-1] + EDGES[1:])


def compute_E(f, vec, n_pairs=100_000_000, seed=0, chunk=10_000_000):
    """E(alpha): pair average of f(x)*f(y) in 2-deg bins.

    f    : (NPIX,) dichotomized sky, values in [-1,1]
    vec  : (3, NPIX) unit vectors of UNMASKED pixels only; f already masked
           (masked pixels excluded by passing only unmasked pixels)
    Pairs (i,j) with i,j uniform over unmasked pixels, i != j.
    """
    rng = np.random.default_rng(seed)
    n = vec.shape[1]
    tot = np.zeros(NBIN)
    cnt = np.zeros(NBIN)
    done = 0
    while done < n_pairs:
        m = min(chunk, n_pairs - done)
        i = rng.integers(0, n, size=m)
        j = rng.integers(0, n - 1, size=m)
        j = np.where(j >= i, j + 1, j)  # j != i, uniform over rest
        cosang = np.clip(np.einsum('kn,kn->n', vec[:, i], vec[:, j]), -1, 1)
        ang = np.degrees(np.arccos(cosang))
        b = np.digitize(ang, EDGES) - 1
        ok = (b >= 0) & (b < NBIN)
        ff = f[i] * f[j]
        tot += np.bincount(b[ok], weights=ff[ok], minlength=NBIN)
        cnt += np.bincount(b[ok], minlength=NBIN)
        done += m
    E = np.full(NBIN, np.nan)
    ok = cnt > 0
    E[ok] = tot[ok] / cnt[ok]
    # per-bin standard error of the mean (values in [-1,1])
    se = np.full(NBIN, np.nan)
    se[ok] = 1.0 / np.sqrt(cnt[ok])
    return E, cnt, se


def eval_tetrad(E, angles_deg):
    """|C| at specific angles (a1,a2,a3,a4); minus sign on a4."""
    bi = np.digitize(np.asarray(angles_deg, dtype=float), EDGES) - 1
    assert np.all((bi >= 0) & (bi < NBIN)), "angles outside (6,180]"
    e = E[bi]
    assert np.all(np.isfinite(e)), "empty E bins hit"
    return float(abs(e[0] + e[1] + e[2] - e[3]))


def tetrad_scan(E, n_tetrads=3_000_000, seed=0, chunk=300_000,
                extra_tetrads=None):
    """Max |C| over random REALIZABLE tetrads (unit-vector quadruples).

    Returns (maxC, best_angles, frac_ge_2, n_eval).
    frac_ge_2 = fraction of evaluated tetrads with |C| >= 2 (Table-6 analogue).
    """
    rng = np.random.default_rng(seed)
    best = -1.0
    best_ang = None
    n_ge2 = 0
    n_eval = 0
    for s in range(0, n_tetrads, chunk):
        m = min(chunk, n_tetrads - s)
        A = rng.normal(size=(m, 3)); A /= np.linalg.norm(A, axis=1, keepdims=True)
        Ap = rng.normal(size=(m, 3)); Ap /= np.linalg.norm(Ap, axis=1, keepdims=True)
        B = rng.normal(size=(m, 3)); B /= np.linalg.norm(B, axis=1, keepdims=True)
        Bp = rng.normal(size=(m, 3)); Bp /= np.linalg.norm(Bp, axis=1, keepdims=True)
        angs = np.degrees(np.arccos(np.clip(np.stack([
            np.einsum('ni,ni->n', A, B),
            np.einsum('ni,ni->n', A, Bp),
            np.einsum('ni,ni->n', Ap, B),
            np.einsum('ni,ni->n', Ap, Bp)]), -1, 1)))
        bi = np.digitize(angs, EDGES) - 1
        ok = np.all((bi >= 0) & (bi < NBIN), axis=0)
        if not np.any(ok):
            continue
        e = E[bi[:, ok]]
        C = np.abs(e[0] + e[1] + e[2] - e[3])
        C[~np.isfinite(C)] = -1.0
        n_ge2 += int(np.sum(C >= 2.0))
        n_eval += int(np.sum(C >= 0))
        k = int(np.argmax(C))
        if C[k] > best:
            best = float(C[k])
            best_ang = angs[:, ok][:, k]
    if extra_tetrads:
        for t in extra_tetrads:
            c = eval_tetrad(E, t)
            n_eval += 1
            if c >= 2.0:
                n_ge2 += 1
            if c > best:
                best, best_ang = c, np.asarray(t, dtype=float)
    return best, best_ang, (n_ge2 / n_eval if n_eval else 0.0), n_eval


# The paper's reported best tetrads (v2), minus sign on the last angle.
PAPER_TETRADS = {
    'wmap_best': [133, 135, 137, 45],
    'planck_commander_best': [133, 135, 137, 45],
}

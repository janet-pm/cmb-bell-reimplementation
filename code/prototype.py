"""Step 1 prototype: verify the SCOPE.md section-3 bound numerically at low Nside.

Bound: for ANY sky function f: sphere -> [-1,1] and ANY realizable tetrad of
unit vectors (a,a',b,b'), |E(a1)+E(a2)+E(a3)-E(a4)| <= 2, where E(ai) are pair
averages of f(x)f(y) at the tetrad's mutual angles.

Checks:
  - random +-1 sky -> max |C| ~ 0
  - dipole sky f=sign(z) -> max |C| ~ 1.9965 (bound is tight)
  - unrestricted scan (all bin quadruples) CAN exceed 2 -> shows the
    realizability constraint is what enforces the bound.
"""
import numpy as np
import healpy as hp

NSIDE = 8
NPIX = hp.nside2npix(NSIDE)
VEC = np.array(hp.pix2vec(NSIDE, np.arange(NPIX)))  # (3, NPIX)

# 2-degree bins, 6 < alpha <= 180  (90 bins, paper's convention)
BIN_W = 2.0
EDGES = np.arange(6.0, 180.0 + BIN_W, BIN_W)
NBIN = len(EDGES) - 1
CENTERS = 0.5 * (EDGES[:-1] + EDGES[1:])


def pair_E_allpairs(f):
    """Exact E(alpha) from all pairs (feasible at low Nside)."""
    # all unordered pairs i<j
    i, j = np.triu_indices(NPIX, k=1)
    cosang = np.clip(VEC[:, i].T @ VEC[:, j] if False else
                     np.einsum('ki,kj->ij', VEC[:, i], VEC[:, j]), -1, 1)
    # einsum over pairs: use chunked dot to avoid huge memory
    ang = np.degrees(np.arccos(cosang))
    ff = f[i] * f[j]
    tot = np.bincount(np.digitize(ang, EDGES) - 1, weights=ff, minlength=NBIN)
    cnt = np.bincount(np.digitize(ang, EDGES) - 1, minlength=NBIN)
    E = np.full(NBIN, np.nan)
    ok = cnt > 0
    E[ok] = tot[ok] / cnt[ok]
    return E, cnt


def pair_E_allpairs_chunked(f, chunk=200000):
    i, j = np.triu_indices(NPIX, k=1)
    npair = len(i)
    tot = np.zeros(NBIN)
    cnt = np.zeros(NBIN)
    for s in range(0, npair, chunk):
        e = min(s + chunk, npair)
        vi = VEC[:, i[s:e]]
        vj = VEC[:, j[s:e]]
        cosang = np.clip(np.einsum('kn,kn->n', vi, vj), -1, 1)
        ang = np.degrees(np.arccos(cosang))
        b = np.digitize(ang, EDGES) - 1
        ok = (b >= 0) & (b < NBIN)
        ff = f[i[s:e]] * f[j[s:e]]
        tot += np.bincount(b[ok], weights=ff[ok], minlength=NBIN)
        cnt += np.bincount(b[ok], minlength=NBIN)
    E = np.full(NBIN, np.nan)
    ok = cnt > 0
    E[ok] = tot[ok] / cnt[ok]
    return E, cnt


def random_unit_vectors(n, rng):
    v = rng.normal(size=(n, 3))
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    return v


def tetrad_scan(E, n_tetrads=2000000, seed=0):
    """Max |C| over random REALIZABLE tetrads. Minus sign on a'.b' (4th)."""
    rng = np.random.default_rng(seed)
    best = 0.0
    best_ang = None
    chunk = 200000
    for s in range(0, n_tetrads, chunk):
        n = min(chunk, n_tetrads - s)
        a = random_unit_vectors(n, rng)
        ap = random_unit_vectors(n, rng)
        b = random_unit_vectors(n, rng)
        bp = random_unit_vectors(n, rng)
        angs = np.degrees(np.arccos(np.clip(np.stack([
            np.einsum('ni,ni->n', a, b),
            np.einsum('ni,ni->n', a, bp),
            np.einsum('ni,ni->n', ap, b),
            np.einsum('ni,ni->n', ap, bp),
        ]), -1, 1)))
        bi = np.digitize(angs, EDGES) - 1
        ok = np.all((bi >= 0) & (bi < NBIN), axis=0)
        if not np.any(ok):
            continue
        e1 = E[bi[0][ok]]
        e2 = E[bi[1][ok]]
        e3 = E[bi[2][ok]]
        e4 = E[bi[3][ok]]
        C = np.abs(e1 + e2 + e3 - e4)
        C[~np.isfinite(C)] = -1.0  # ignore tetrads hitting empty bins
        k = np.argmax(C)
        if C[k] > best:
            best = float(C[k])
            best_ang = angs[:, ok][:, k]
    return best, best_ang


def eval_tetrad(E, angles_deg):
    """Evaluate |C| at specific angles (a1,a2,a3,a4), minus on a4."""
    bi = np.digitize(np.asarray(angles_deg), EDGES) - 1
    assert np.all((bi >= 0) & (bi < NBIN)), "angles outside (6,180]"
    e = E[bi]
    return float(abs(e[0] + e[1] + e[2] - e[3]))


def unrestricted_max(E):
    """Max |E_i+E_j+E_k-E_l| over ALL bin quadruples (no realizability)."""
    mx, mn = np.nanmax(E), np.nanmin(E)
    return float(max(3 * mx - mn, mx - 3 * mn))


if __name__ == '__main__':
    rng = np.random.default_rng(42)

    print("=== random +-1 sky ===")
    f = rng.choice([-1.0, 1.0], size=NPIX)
    E, cnt = pair_E_allpairs_chunked(f)
    print("bins with pairs:", np.sum(cnt > 0), "/", NBIN)
    best, bang = tetrad_scan(E)
    print(f"max |C| over realizable tetrads: {best:.4f}")
    print(f"unrestricted max |C|:            {unrestricted_max(E):.4f}")
    print(f"|C| at paper tetrad (133,135,137,45): "
          f"{eval_tetrad(E, [133, 135, 137, 45]):.4f}")

    print("\n=== dipole sky f=sign(z) ===")
    z = VEC[2]
    f = np.sign(z)
    f[f == 0] = 1.0
    E, cnt = pair_E_allpairs_chunked(f)
    best, bang = tetrad_scan(E)
    print(f"max |C| over realizable tetrads: {best:.4f}  (expect ~1.9965)")
    print(f"best tetrad angles: {np.round(bang, 1) if bang is not None else None}")
    print(f"unrestricted max |C|:            {unrestricted_max(E):.4f}")
    print(f"|C| at paper tetrad (133,135,137,45): "
          f"{eval_tetrad(E, [133, 135, 137, 45]):.4f}")
    # analytic: E(alpha) = 1 - 2*alpha/pi for dipole-sign
    Ea = 1 - 2 * np.radians(CENTERS) / np.pi
    print(f"E(alpha) matches 1-2a/pi: max abs diff = "
          f"{np.nanmax(np.abs(E - Ea)):.5f}")

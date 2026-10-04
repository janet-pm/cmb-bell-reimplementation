"""Full pipeline on one prepared map: dichotomy (mixed+unmixed) -> E(alpha) ->
tetrad scan. Saves per-method results.

Mixed:    f(x) = e+(x)/ea - 1, e+ = positive part of [u(x)-ea, u(x)+ea]
Unmixed:  f(x) = sign(u(x))
"""
import numpy as np, healpy as hp, os, sys, json, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from pipeline import compute_E, tetrad_scan, eval_tetrad, CENTERS

DATA = os.path.expanduser('~/workspace/vsl-planck/cmb-bell/rerun/data')
RES = os.path.expanduser('~/workspace/vsl-planck/cmb-bell/rerun/results')

FIXED_TETRADS = {
    'paper_wmap_best': [133, 135, 137, 45],
    'paper_commander_best': [133, 135, 137, 45],
    'v1_wmap_best': [133, 49, 137, 139],
}


def dichotomize(u, ea, method):
    if method == 'unmixed':
        f = np.sign(u)
        f[f == 0] = 1.0
        return f
    # mixed: soft sign via error interval
    lo = u - ea
    hi = u + ea
    epos = np.maximum(0.0, hi - np.maximum(lo, 0.0))
    return epos / ea - 1.0


def analyze(prep_path, n_pairs=100_000_000, n_tetrads=3_000_000, seed=0,
            tag='data'):
    d = np.load(prep_path)
    u, mask = d['u'], d['mask'].astype(bool)
    ea = float(d['ea_noise'])
    nside = hp.npix2nside(len(u))
    vec_all = np.array(hp.pix2vec(nside, np.arange(len(u))))
    vec = vec_all[:, mask]
    out = {'tag': tag, 'prep': os.path.basename(prep_path),
           'nside': nside, 'npix_unmasked': int(mask.sum()),
           'ea_noise_uK': ea, 'n_pairs': n_pairs, 'n_tetrads': n_tetrads}
    for method in ['unmixed', 'mixed']:
        t0 = time.time()
        f = dichotomize(u[mask], ea, method)
        E, cnt, se = compute_E(f, vec, n_pairs=n_pairs, seed=seed)
        best, bang, frac, nev = tetrad_scan(
            E, n_tetrads=n_tetrads, seed=seed + 1,
            extra_tetrads=list(FIXED_TETRADS.values()))
        fixed = {k: eval_tetrad(E, t) for k, t in FIXED_TETRADS.items()}
        out[method] = {
            'max_absC': best,
            'best_tetrad_deg': [float(x) for x in bang],
            'frac_ge2': frac,
            'n_tetrads_eval': nev,
            'fixed_tetrads': fixed,
            'E_min_pairs': float(cnt.min()),
            'E_median_se': float(np.nanmedian(se)),
            'E_alpha': [float(x) for x in E],
            'seconds': time.time() - t0,
        }
        print(f"[{tag}/{method}] max|C|={best:.4f} at "
              f"{np.round(bang,1)}  frac>=2={frac:.5f}  "
              f"fixed(paper)={fixed['paper_wmap_best']:.4f}  "
              f"({time.time()-t0:.0f}s)", flush=True)
    return out


if __name__ == '__main__':
    # usage: analyze.py <prep.npz> <tag> [n_pairs] [n_tetrads] [seed]
    prep_path = sys.argv[1]
    tag = sys.argv[2] if len(sys.argv) > 2 else 'data'
    n_pairs = int(sys.argv[3]) if len(sys.argv) > 3 else 100_000_000
    n_tetrads = int(sys.argv[4]) if len(sys.argv) > 4 else 3_000_000
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 0
    out = analyze(prep_path, n_pairs, n_tetrads, seed, tag)
    os.makedirs(RES, exist_ok=True)
    with open(f'{RES}/{tag}.json', 'w') as fh:
        json.dump(out, fh)
    print('wrote', f'{RES}/{tag}.json', flush=True)

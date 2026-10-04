"""Gaussian null simulations + phase-randomized surrogates.

Gaussian sim: isotropic T at Nside=512 from Planck PR3 best-fit D_ell
  (converted to C_ell), x 5' beam x N512 pixel window, lmax=1500.
Phase-rand: exact phase randomization of the real harmonic-degraded dT map
  (keep |a_lm|, randomize phases; m=0 coefficients keep random sign).

Both: degrade to N64 by both methods, add white noise RMS=ea,
restore Table-1 dipole, apply the analysis mask, run analyze.py.
"""
import numpy as np, healpy as hp, os, sys, subprocess, argparse

DATA = os.path.expanduser('~/workspace/vsl-planck/cmb-bell/rerun/data')
CODE = os.path.expanduser('~/workspace/vsl-planck/cmb-bell/rerun/code')
NS512, NS64 = 512, 64
NPIX2048 = hp.nside2npix(2048)
DIP_A, DIP_L, DIP_B = 3362.08, 264.021, 48.253


def dipole_map(nside):
    vec = np.array(hp.pix2vec(nside, np.arange(hp.nside2npix(nside))))
    d = np.array(hp.ang2vec(np.radians(90.0 - DIP_B), np.radians(DIP_L)))
    return DIP_A * (vec.T @ d)


def load_bf_cl():
    # PLA theory file: columns L TT TE EE BB PP, TT is D_ell=l(l+1)Cl/2pi (uK^2)
    a = np.loadtxt(f'{DATA}/planck_bf_cl.txt')
    ell = a[:, 0].astype(int)
    Dl = a[:, 1]
    lmax = 1500
    cl = np.zeros(lmax + 1)
    ok = (ell >= 2) & (ell <= lmax)
    cl[ell[ok]] = Dl[ok] * 2 * np.pi / (ell[ok] * (ell[ok] + 1))
    fwhm = np.radians(5.0 / 60.0)
    bl = hp.gauss_beam(fwhm, lmax=lmax)
    # pixel window from local file (avoid network download)
    from astropy.io import fits
    pw_full = fits.getdata('/tmp/pixel_window_n0512.fits')['TEMPERATURE']
    pw = np.ones(lmax + 1)
    n = min(lmax + 1, len(pw_full))
    pw[:n] = pw_full[:n]
    return cl * bl ** 2 * pw ** 2


def load_mask64():
    # truncated mask file: read available rows (float32), missing -> masked
    p = f'{DATA}/COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits'
    nrows = (os.path.getsize(p) - 8640) // 4
    mm = np.fromfile(p, dtype='>f4', count=nrows, offset=8640).astype(float)
    mfull = np.zeros(NPIX2048)
    mfull[:nrows] = mm
    mfull = hp.reorder(mfull, n2r=True)
    return hp.ud_grade(mfull, NS64) > 0.9


def load_ea():
    p = f'{DATA}/prep_hmdiff_neigh.npz'
    if os.path.exists(p):
        return float(np.load(p)['ea_noise'])
    print('WARNING: no HM noise file; ea=10 uK fallback', flush=True)
    return 10.0


def degrade(m512, how):
    if how == 'neigh':
        return hp.ud_grade(m512, NS64)
    alm = hp.map2alm(m512, lmax=3 * NS64 - 1)
    return hp.alm2map(alm, NS64, lmax=3 * NS64 - 1)


def phase_randomize(dT64, seed):
    """Exact phase randomization: keep |a_lm|, random phases (m=0: random sign)."""
    rng = np.random.default_rng(seed)
    alm = hp.map2alm(dT64, lmax=3 * NS64 - 1)
    lmax = 3 * NS64 - 1
    out = np.zeros_like(alm)
    for l in range(lmax + 1):
        i0 = hp.Alm.getidx(lmax, l, 0)
        # m=0: real coefficient -> random sign
        out[i0] = np.abs(alm[i0]) * (1.0 if rng.random() < 0.5 else -1.0)
        for m in range(1, l + 1):
            i = hp.Alm.getidx(lmax, l, m)
            ph = rng.uniform(0, 2 * np.pi)
            out[i] = np.abs(alm[i]) * np.exp(1j * ph)
    return hp.alm2map(out, NS64, lmax=lmax)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=100)
    ap.add_argument('--start', type=int, default=0)
    ap.add_argument('--phase-rand', action='store_true')
    ap.add_argument('--n-pairs', type=int, default=30_000_000)
    ap.add_argument('--n-tetrads', type=int, default=1_500_000)
    args = ap.parse_args()

    mask = load_mask64()
    ea = load_ea()
    dip64 = dipole_map(NS64)
    cl = load_bf_cl()

    dT_data_64 = None
    if args.phase_rand:
        u = np.load(f'{DATA}/prep_smica_harm_fullsky.npz')['u']
        dT_data_64 = np.nan_to_num(u - dip64, nan=0.0)  # missing pix -> 0

    py = f'{CODE}/../venv/bin/python'
    for i in range(args.start, args.start + args.n):
        rng = np.random.default_rng(20_000 + i)
        if args.phase_rand:
            m64 = phase_randomize(dT_data_64, seed=30_000 + i)
            sims64 = {'neigh': m64, 'harm': m64}
        else:
            np.random.seed(10_000 + i)
            m512 = hp.synfast(cl, NS512, lmax=1500, new=True)
            sims64 = {h: degrade(m512, h) for h in ['neigh', 'harm']}
        for h in ['neigh', 'harm']:
            m64 = sims64[h] + rng.normal(0, ea, hp.nside2npix(NS64))
            u = m64 + dip64
            tag = f"{'phaserand' if args.phase_rand else 'sim'}_{i:03d}_{h}"
            prep = f'{DATA}/prep_{tag}.npz'
            np.savez(prep, u=u, mask=mask, ea_noise=np.float64(ea))
            subprocess.run(
                [py, f'{CODE}/analyze.py', prep, tag,
                 str(args.n_pairs), str(args.n_tetrads), str(20_000 + i)],
                check=True, capture_output=True, text=True)
            print(f'done {tag}', flush=True)
    print('SIMS DONE', flush=True)


if __name__ == '__main__':
    main()

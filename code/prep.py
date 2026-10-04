"""Prepare Planck PR3 SMICA data for the clean-room CHSH pipeline.

DATA NOTE (2026-09-27): the PLA download of the full SMICA IQU file stalled
at 2007921544/2013275520 bytes (99.7%). The missing 133,826 rows are the last
NESTED pixels = south polar cap (~6 deg radius, 0.27% of sky). The I field
for all available rows reads cleanly. Missing pixels are excluded via an
extended mask (see below); a complete re-download runs in parallel and key
numbers will be re-verified against it.

Steps (documenting every choice):
 1. Read I_STOKES (K_CMB -> uK) for available rows via numpy (file is NESTED).
    Reorder to RING for all downstream work.
 2. Fit and remove monopole+dipole actually present (hp.fit_dipole on the
    available pixels), so the restored dipole below is exact.
 3. Restore the Table-1 dipole (paper v2, Planck row):
      A = 3362.08 uK_CMB, (l,b) = (264.021, 48.253) deg Galactic.
    The monopole is the dichotomy threshold -> NOT added back.
 4. Degrade to Nside=64 by BOTH paper methods:
      (a) neighbour averaging: weighted ud_grade ignoring missing pixels
      (b) spherical-harmonic truncation: map2alm(lmax=191) -> alm2map(64),
          missing pixels filled with the map mean first (tiny discontinuity
          at the cap edge; affected region is masked anyway).
 5. Mask: common confidence mask ud_graded to 64 (>0.9) ANDed with
    (fraction of good subpixels > 0.99).
 6. Noise for the mixed method: RMS of (HM1-HM2)/2 after same prep
    (fit/remove, NO dipole restore, degrade). Uniform ea. If HM files are
    absent, falls back to ea=10 uK with a printed warning.

Outputs: prep_smica_neigh.npz, prep_smica_harm.npz, prep_hmdiff_*.npz
"""
import numpy as np, healpy as hp, os, sys

DATA = os.path.expanduser('~/workspace/vsl-planck/cmb-bell/rerun/data')
NS64 = 64
N2048 = 2048
NPIX2048 = hp.nside2npix(N2048)

DIP_A, DIP_L, DIP_B = 3362.08, 264.021, 48.253  # paper v2 Table 1, Planck


def read_I_truncated(path):
    """Read I_STOKES (nested) for available rows; return (I_uK_ring, good_ring)."""
    hdr_len = 8640
    fsize = os.path.getsize(path)
    nrows_avail = (fsize - hdr_len) // 40
    assert nrows_avail > 0.99 * NPIX2048, "file too truncated"
    dt = np.dtype([('I', '>f4'), ('rest', '>f4', 9)])
    a = np.fromfile(path, dtype=dt, count=nrows_avail, offset=hdr_len)
    I_nested = np.full(NPIX2048, np.nan)
    I_nested[:nrows_avail] = a['I'].astype(np.float64) * 1e6  # K -> uK
    good_nested = np.zeros(NPIX2048, dtype=bool)
    good_nested[:nrows_avail] = True
    # to RING
    I_ring = hp.reorder(I_nested, n2r=True)
    good_ring = hp.reorder(good_nested.astype(np.float64), n2r=True) > 0.5
    print(f"read {nrows_avail}/{NPIX2048} rows; "
          f"missing={(~good_ring).sum()} pix (south polar cap)",
          flush=True)
    return I_ring, good_ring


def remove_dipole(m, good):
    mm = m.copy()
    mm[~good] = hp.UNSEEN
    mono, dip = hp.fit_dipole(mm)
    vec = np.array(hp.pix2vec(N2048, np.arange(NPIX2048)))
    out = m.copy()
    out[good] = m[good] - (mono + vec[:, good].T @ np.asarray(dip, float))
    return out


def dipole_map(nside):
    vec = np.array(hp.pix2vec(nside, np.arange(hp.nside2npix(nside))))
    d = np.array(hp.ang2vec(np.radians(90.0 - DIP_B), np.radians(DIP_L)))
    return DIP_A * (vec.T @ d)


def degrade_neigh(m, good):
    mf = np.where(good, m, 0.0)
    w = good.astype(np.float64)
    num = hp.ud_grade(mf, NS64)
    den = hp.ud_grade(w, NS64)
    out = np.full_like(num, np.nan)
    ok = den > 0.99
    out[ok] = num[ok] / den[ok]
    return out, ok


def degrade_harm(m, good):
    mf = np.where(good, m, np.nanmean(m[good]))
    alm = hp.map2alm(mf, lmax=3 * NS64 - 1)
    return hp.alm2map(alm, NS64, lmax=3 * NS64 - 1, verbose=False)


def main():
    smica_path = f'{DATA}/COM_CMB_IQU-smica_2048_R3.00_full.fits'
    I_ring, good = read_I_truncated(smica_path)
    np.save(f'{DATA}/smica_I_ring.npy', I_ring)
    np.save(f'{DATA}/smica_good_ring.npy', good)

    mask_path = f'{DATA}/COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits'
    mask64_common = None
    if os.path.exists(mask_path):
        try:
            # mask file may be truncated; read available rows, treat missing as masked
            hdr_len_m = 8640
            msize = os.path.getsize(mask_path)
            nrows_m = (msize - hdr_len_m) // 4
            if nrows_m >= 0.95 * NPIX2048:
                mm = np.fromfile(mask_path, dtype='>f4', count=nrows_m,
                                 offset=hdr_len_m).astype(np.float64)
                mfull = np.zeros(NPIX2048)
                mfull[:nrows_m] = mm
                mfull = hp.reorder(mfull, n2r=True)
                mask64_common = hp.ud_grade(mfull, NS64) > 0.9
                print(f'mask: {nrows_m}/{NPIX2048} rows, '
                      f'sky frac={mask64_common.mean():.3f}', flush=True)
        except Exception as e:
            print(f'mask read failed: {e}', flush=True)
    if mask64_common is None:
        print('WARNING: no usable mask; running FULL-SKY (theorem test)', flush=True)
        mask64_common = np.ones(hp.nside2npix(NS64), dtype=bool)

    dT = remove_dipole(I_ring, good)
    dip64 = dipole_map(NS64)

    u_neigh, ok_neigh = degrade_neigh(dT, good)
    u_harm = degrade_harm(dT, good)
    # harm: mask N64 pixels whose subpixels were mostly missing
    w64 = hp.ud_grade(good.astype(np.float64), NS64)
    ok_harm = w64 > 0.99

    for name, u, ok in [('neigh', u_neigh, ok_neigh),
                        ('harm', u_harm, ok_harm)]:
        uu = u + dip64
        okfin = ok & np.isfinite(uu)
        np.savez(f'{DATA}/prep_smica_{name}_fullsky.npz', u=uu,
                 mask=okfin, dip_amp=DIP_A, dip_l=DIP_L, dip_b=DIP_B)
        print(f'saved prep_smica_{name}_fullsky.npz  kept={okfin.sum()}',
              flush=True)
        mask = mask64_common & okfin
        np.savez(f'{DATA}/prep_smica_{name}_masked.npz', u=uu, mask=mask,
                 dip_amp=DIP_A, dip_l=DIP_L, dip_b=DIP_B)
        print(f'saved prep_smica_{name}_masked.npz  kept={mask.sum()}  '
              f'u.std={np.nanstd(uu[mask]):.1f} uK', flush=True)

    # noise level from half-mission difference (no dipole restore)
    hm_paths = [f'{DATA}/COM_CMB_IQU-smica_2048_R3.00_halfmission-1.fits',
                f'{DATA}/COM_CMB_IQU-smica_2048_R3.00_halfmission-2.fits']
    if all(os.path.exists(p) for p in hm_paths):
        hms = []
        for p in hm_paths:
            Ir, gd = read_I_truncated(p)
            hms.append(remove_dipole(Ir, gd))
        hmdiff = (hms[0] - hms[1]) / 2.0
        gd = np.isfinite(hmdiff)
        n_neigh, _ = degrade_neigh(hmdiff, gd)
        n_harm = degrade_harm(hmdiff, gd)
        eas = {}
        for name, n in [('neigh', n_neigh), ('harm', n_harm)]:
            ea = float(np.nanstd(n[mask64_common & np.isfinite(n)]))
            eas[name] = ea
            np.savez(f'{DATA}/prep_hmdiff_{name}.npz', u=n,
                     mask=mask64_common & np.isfinite(n),
                     ea_noise=np.float64(ea))
            print(f'saved prep_hmdiff_{name}.npz  ea={ea:.2f} uK', flush=True)
    else:
        print('WARNING: HM files missing; using fallback ea=10 uK', flush=True)
        eas = {'neigh': 10.0, 'harm': 10.0}
    for name in ['neigh', 'harm']:
        for sky in ['fullsky', 'masked']:
            p = f'{DATA}/prep_smica_{name}_{sky}.npz'
            d = dict(np.load(p))
            d['ea_noise'] = np.float64(eas[name])
            np.savez(p, **d)
    print('PREP DONE', flush=True)


if __name__ == '__main__':
    main()

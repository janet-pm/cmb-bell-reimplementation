"""Prototype v2 at Nside=64 (paper's resolution): verify the bound with the
real pipeline machinery. Dry run for the Planck analysis."""
import numpy as np, healpy as hp, sys
sys.path.insert(0, 'code')
from pipeline import compute_E, tetrad_scan, eval_tetrad, CENTERS, EDGES

NSIDE = 64
NPIX = hp.nside2npix(NSIDE)
VEC = np.array(hp.pix2vec(NSIDE, np.arange(NPIX)))
rng = np.random.default_rng(7)

print("=== dipole sky f=sign(z), Nside=64 ===", flush=True)
z = VEC[2]; f = np.sign(z); f[f == 0] = 1.0
E, cnt, se = compute_E(f, VEC, n_pairs=20_000_000, seed=1)
Ea = 1 - 2 * np.radians(CENTERS) / np.pi
print(f"max |E - (1-2a/pi)| = {np.nanmax(np.abs(E - Ea)):.4f}", flush=True)
print(f"min pairs/bin = {cnt.min():.0f}", flush=True)
best, bang, frac, nev = tetrad_scan(E, n_tetrads=1_000_000, seed=2,
    extra_tetrads=[[133, 135, 137, 45]])
print(f"max |C| = {best:.4f} at {np.round(bang,1)}  (expect ~1.99, must be <=2)",
      flush=True)
print(f"frac |C|>=2: {frac:.6f} over {nev} tetrads", flush=True)
print(f"|C| at paper tetrad: {eval_tetrad(E,[133,135,137,45]):.4f}", flush=True)

print("\n=== random +-1 sky, Nside=64 ===", flush=True)
f = rng.choice([-1.0, 1.0], size=NPIX)
E, cnt, se = compute_E(f, VEC, n_pairs=20_000_000, seed=3)
best, bang, frac, nev = tetrad_scan(E, n_tetrads=1_000_000, seed=4,
    extra_tetrads=[[133, 135, 137, 45]])
print(f"max |C| = {best:.4f}  (expect ~0.01)", flush=True)
print(f"frac |C|>=2: {frac:.6f}", flush=True)
print("PROTOTYPE DONE", flush=True)

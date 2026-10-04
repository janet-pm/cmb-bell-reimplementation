# Clean-room rerun: Dale et al., arXiv:2502.13846v2 (cosmic CHSH)

**Date:** 2026-09-27
**Code:** `code/` (pipeline.py, prep.py, analyze.py, sims.py, prototype.py, prototype2.py)
**Results:** `results/rerun.json` (machine-readable; this report is the narrative)

## Verdict

**The published SMICA value is not reproduced by the all-sky method
described in the paper.** Dataset-specific values matter here: for Planck
SMICA at the tetrad below, the paper's Table 3 reports **2.2755**
(neighbour averaging) / **2.2764** (spherical harmonics). The ≈2.42–2.44
figures sometimes quoted generically as "the paper's value" belong to
other datasets: WMAP Q/V/W ≈ 2.425 (their Table 2) and Planck NILC
2.4376 / 2.4378 (Table 3). Our rerun is SMICA only, so 2.276 is the
comparator throughout this report.

- Near-full-sky (the case the |C| ≤ 2 theorem covers; 49,021 of 49,152
  Nside=64 pixels due to a truncated download): |C| = **2.001** at the
  paper's best tetrad (133, 135, 137, 45), the same tetrad that tops
  their SMICA table. Optimized over 3M realizable tetrads: 2.027
  (binning/sampling noise, as predicted by the prototype). Note on that
  3M figure: the tetrads are random realizable direction quadruples
  evaluated against a single precomputed E(alpha) curve built from 100M
  pixel pairs, not 3M independent four-point sky samples; the 2.027
  maximum is a binning / extreme-value effect over that curve.
- With the Planck common confidence mask applied: |C| = **2.667** at the
  paper's tetrad. The mask alone inflates the statistic.
- Pure-dipole control (no CMB, no noise): near-full-sky |C| = 2.007;
  masked |C| = 2.669. **Masking is sufficient to create the entire effect.**
- Galactic-latitude cuts on a pure dipole: |b| > 15 deg gives |C| = **2.40**.
  The paper's SMICA value (2.276) sits on this same mask-response curve,
  between the |b| > 10 deg cut (2.235) and the |b| > 15 deg cut (2.401).

The analytic bound (SCOPE.md) is confirmed empirically: for genuine all-sky
pair averages of a single f: S^2 -> [-1,1], |C| <= 2 holds. The paper's
SMICA 2.276 (and likewise its >2 values for the other datasets) cannot
come from the full-sky computation as described. Masking or directional
selection is a sufficient and quantitatively plausible mechanism, but the
authors' code is unavailable and the exact broken or undocumented step cannot
be identified: it could be an applied mask, a sky cut, or another step that
breaks the uniform rotation averaging the proof requires (masked pair averages
are not uniform over the rotation group, so E(alpha) = rotation average fails
and |E(alpha)| is systematically inflated).

## Results

### Real Planck PR3 SMICA (full mission)

Dipole restored per paper Table 1 (A=3362.08 uK, l=264.021 deg, b=48.253 deg).
Nside=64, 100M pairs, 3M realizable tetrads (evaluated against the
single E(alpha) curve computed from those 100M pairs; see the Verdict
note). ea=10 uK fallback (HM files
unavailable; mixed/unmixed agree to 0.0001, so the choice is immaterial).

| config | method | max|C| (3M tetrads) | paper tetrad | frac >= 2 |
|---|---|---|---|---|
| neighbour, near-full-sky | unmixed | 2.0270 | 2.0013 | 0.00029 |
| neighbour, near-full-sky | mixed | 2.0270 | 2.0012 | 0.00028 |
| harmonic, near-full-sky | unmixed | 2.0270 | 2.0014 | 0.00026 |
| harmonic, near-full-sky | mixed | 2.0272 | 2.0015 | 0.00026 |
| neighbour, masked | unmixed | 2.6706 | 2.6665 | 0.04087 |
| neighbour, masked | mixed | 2.6700 | 2.6659 | 0.04083 |
| harmonic, masked | unmixed | 2.6693 | 2.6648 | 0.04076 |
| harmonic, masked | mixed | 2.6704 | 2.6660 | 0.04083 |

Both degradation methods and both dichotomies agree. Near-full-sky E(alpha)
matches the ideal dipole 1-2alpha/pi to 0.002.

### Prototype bound checks (Nside=64 Monte Carlo)

- Random +/-1 sky: max optimized |C| = 0.0211, frac>=2 = 0.
- Dipole sky: paper tetrad |C| = 1.9894; max over 1M tetrads 2.0394;
  frac>=2 = 0.000421. (Excess over 2 is binning/sampling noise.)
- Analytic dipole E at bin centers: max |C| = 2.0222, paper tetrad exactly 2.0000.

### Mask-artifact mechanism (pure dipole, no CMB/noise)

| mask | \|C\| at (133,135,137,45) |
|---|---|
| none (near-full-sky) | 2.007 |
| common confidence mask | 2.669 |
| \|b\| > 5 deg | 2.109 |
| \|b\| > 10 deg | 2.235 |
| \|b\| > 15 deg | 2.401 |
| \|b\| > 20 deg | 2.654 |
| \|b\| > 30 deg | 3.253 |

The paper's SMICA value, 2.276, falls on this curve between the
|b|>10 and |b|>15 deg cuts; its WMAP value, ≈2.42, falls between
|b|>15 and |b|>20 deg cuts.
A mask that keeps ~75% of sky (like the common mask or |b|>15 deg) inflates
|C| from 2.00 to 2.4-2.7. The bound |C|<=2 assumes uniform all-sky pair
averaging; any mask breaks the rotation-averaging step of the proof.

### Gaussian null simulations

100 Gaussian CMB realizations (Planck PR3 best-fit Cl, 5' beam, N512 pixel
window) + restored Table-1 dipole + common mask + white noise (ea=10 uK),
analyzed identically to the data (30M pairs, 1.5M tetrads, both degradation
methods, both dichotomies).

- Masked max|C|: mean **2.6680**, std 0.0101, range [2.643, 2.698]
  (n=200 entries: 100 independent Gaussian skies × 2 degradation methods).
- The data's fixed paper-tetrad |C|=2.667 sits at the 53.5th percentile of
  the 200 sims' fixed-tetrad values (sim fixed: mean 2.6658, std 0.0097):
  fully consistent with the null. The data's optimized max 2.6706 sits at
  the 61.5th percentile of the sims' optimized maxima. The "violation" is
  what a Gaussian sky + dipole + mask produces.
- Unmixed and mixed agree to <0.001 in every sim.

### Phase-randomized surrogates

20 exact phase-randomized versions of the real harmonic-degraded dT map
(|a_lm| preserved, phases randomized; m=0: random sign), + dipole + mask
(n=40 entries: 20 independent surrogates × 2 degradation labels).

- Masked max|C|: mean **2.6716**, std 0.0095, range [2.652, 2.690].
- Indistinguishable from Gaussian sims: no non-Gaussian phase signal is needed
  to explain the data. The entire effect is dipole + mask geometry.

### Half-mission noise diagnostic

Not run: HM files could not be downloaded (PLA server stalls). The ea=10 uK
fallback is validated by the mixed/unmixed agreement (they differ only in
pixels with |u| < ea, of which there are essentially none in a
dipole-dominated map).

## Limitations and caveats

1. **Truncated downloads.** The PLA server repeatedly stalled near the end of
   large files and does not support byte ranges, so a resumed download is
   impossible. SMICA: 2007921544/2013275520 bytes (99.7%); missing 133,826
   NESTED rows = south polar cap (~6 deg), 131 Nside=64 pixels excluded.
   Mask: 195837019/201335040 bytes (97.3%); missing rows treated as masked.
   Near-full-sky results are therefore 49,021/49,152 pixels; the theorem
   remains analytically valid and the prototype supports it, but the data
   run is not an exact all-sky numerical implementation.
2. **Map units.** SMICA I_STOKES is in K_CMB (FITS header); converted to uK.
   The Table-1 dipole (3362.08 uK) was restored after removing the fitted
   dipole; monopole not restored (it is the threshold).
3. **Mask file format.** The "Int" mask is actually float32 (TFORM1='1E').
4. **Theory Cl.** The PLA theory file gives D_ell, converted to C_ell.
5. **The masked |C|>2 does not challenge the theorem.** The proof requires
   uniform all-sky pair averages. Masked E(alpha) is a different statistic.
6. **No author contact** was made; no public circulation without approval.

## Files

- `results/rerun.json`: all numbers below, machine-readable.
- `results/smica_*.json`: per-config full outputs (E(alpha), tetrads, timing).
- `results/sim_*.json`, `results/phaserand_*.json`: null simulations.
- `data/prep_*.npz`: prepared maps (u, mask, ea).

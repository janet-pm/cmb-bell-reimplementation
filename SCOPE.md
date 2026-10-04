# CMB Bell-inequality claim: investigation scope

**Paper:** R. Dale, R. Gandía, J. A. Morales-Lladosa, R. Lapiedra,
"Violation of Bell inequalities from Cosmic Microwave Background data,"
arXiv:2502.13846.
**Versions:** v1 (2025-02-19), v2 (2025-09-04).
**Status of this document:** Read-only research scope. No contact with authors.
No public post or paper without James's explicit approval.

---

## 1. Bottom line

The paper's headline result, a measured CHSH value around |C| = 2.42
from WMAP and Planck data (dataset-specific maxima in v2 range from
2.20 for Planck COMMANDER to 2.44 for Planck NILC; SMICA, the dataset
rerun in the accompanying report, is 2.2755 / 2.2764), **cannot be produced by an honest
implementation of the procedure the paper describes**. This is not a
matter of interpretation, noise, or look-elsewhere effects. It is a
mathematical impossibility, shown below from the paper's own equations.

What this means:

- The "violation" is a computational artifact (a bug or a step that does
  not do what the text says it does), not evidence against local realism
  and not evidence for a quantum origin of CMB fluctuations.
- The paper's own proof of the cosmic CHSH inequality is valid, and the
  geometric tetrad its proof needs **does exist** for the reported
  angles (constructed explicitly in section 3). The proof therefore
  applies, and it caps the honestly computed quantity at 2.
- The earlier COBE analysis by the same group (Phys. Rev. D 107, 023506,
  2023) found |C| = 1.71, comfortably below 2. That is exactly what the
  theorem predicts for honest data. The 2025 "violation" is the anomaly,
  not the 2023 non-violation.

Recommended next step: a clean-room reimplementation of their pipeline
on Planck data plus Gaussian CMB simulations (section 7). An honest
reimplementation must satisfy |C| <= 2 on every sky, real or simulated.
If it does, the published 2.42 is definitively an artifact of their
code, and the interesting follow-up is identifying the exact breaking
step.

---

## 2. What the paper claims

The authors dichotomize the CMB sky: F(x) = +1 where the temperature is
above the monopole, -1 where below (with a "mixed probability" softening
near zero that keeps every value in [-1, 1]; details in section 3).
For each angular separation alpha they estimate a correlator

    E(alpha) = P(++,alpha) + P(--,alpha) - P(+-,alpha) - P(-+,alpha),

the difference between same-sign and opposite-sign pair fractions at
that separation. They then scan quadruples of angles (tetrads) and form

    |C| = |E(a1) + E(a2) + E(a3) - E(a4)|,

claiming the local-realism bound |C| <= 2 (their Eq. 1). Scanning WMAP
(Q, V, W bands) and Planck (COMMANDER, NILC, SEVEM, SMICA pipelines),
degraded to HEALPix Nside=64, they report maxima around 2.42 (WMAP),
2.44 (Planck NILC), 2.32 (Planck SEVEM), 2.276 (Planck SMICA:
2.2755 neighbour averaging, 2.2764 spherical harmonics, both at
(133, 135, 137, 45)), 2.20 (Planck COMMANDER), all above 2, and interpret
this as a failure of local realism.

Their measured angular correlator is well fit by E(alpha) ~ 0.85 cos(alpha)
at intermediate angles (their Figures 1-2 show this fit explicitly).

---

## 3. The decisive finding: |C| > 2 is algebraically impossible here

### 3.1 The theorem

Fix any four unit vectors a, a', b, b' and any function
f from the sphere to [-1, 1]. Define the pair-average correlators

    E1 = average of f(x)f(y) over all pairs with x.y = cos(a1), etc.

**Claim:** |E1 + E2 + E3 - E4| <= 2, where the four angles are the
mutual angles of the tetrad: a.b = cos(a1), a.b' = cos(a2),
a'.b = cos(a3), a'.b' = cos(a4).

**Proof.** A pair average at separation a1 equals the rotation average
of f(Ra)f(Rb): as R runs over all rotations, (Ra, Rb) sweeps uniformly
over all pairs at that separation. So

    E1 + E2 + E3 - E4
      = average over rotations R of
        [f(Ra)f(Rb) + f(Ra)f(Rb') + f(Ra')f(Rb) - f(Ra')f(Rb')].

For each fixed rotation, the bracket is an ordinary CHSH combination of
four numbers in [-1, 1], hence in [-2, 2]. The average of quantities in
[-2, 2] is in [-2, 2]. QED.

This uses nothing about cosmology, quantum mechanics, or the CMB. It is
an identity for pair averages. Binning the angles (their 2-degree bins)
changes each E by ~0.001, not 0.42. Random pair subsampling is unbiased.
Maximizing over tetrads or over random subsamples cannot break a bound
that holds for every single instance.

### 3.2 The tetrad exists for their angles

The proof needs four unit vectors with the prescribed mutual angles.
For the paper's best WMAP/Planck tetrad (133, 135, 137, 45) degrees,
with the minus sign on 45 degrees, such vectors exist. Explicitly
(in the x-z plane):

    a  = (-0.7314, 0, -0.6820)
    a' = ( 0.7071, 0, -0.7071)
    b  = ( 0, 0, 1)
    b' = ( 1, 0, 0)

Check: a.b = -0.6820 = cos(133), a.b' = -0.7314 = cos(137),
a'.b = -0.7071 = cos(135), a'.b' = 0.7071 = cos(45). All unit vectors,
all dot products correct. (The plus-pair angles {133, 137, 135} are a
permutation of the paper's {133, 135, 137}; the sum is unaffected.)

Feasibility was also verified by certified concave optimization over the
Gram matrix (maximum minimum-eigenvalue ~ 0, i.e. marginally feasible)
and by an exact interval-intersection construction. The same holds for
the other reported tetrads, e.g. Planck COMMANDER (43, 45, 47, 135).

### 3.3 Numerical confirmation

- Random +-1 sky: |C| ~ 0.002 (well under 2).
- Dipole sky f = sign(z): E(alpha) = 1 - 2*alpha/pi, giving
  |C| = 1.9965 ~ 2.0 exactly as theory predicts. **The bound is tight:
  2.0 is achievable, 2.0001 is not.**
- Therefore the published 2.42 is impossible for pair averages of any
  sky function. Their E(a1..a4) are not honest pair averages.

### 3.4 Why the 2023 COBE paper got 1.71 and this one gets 2.42

If E(alpha) ~ A*cos(alpha), the best CHSH combination is ~2*sqrt(2)*A.
COBE's correlator amplitude was A ~ 0.60, giving 2*sqrt(2)*0.60 = 1.70,
exactly their 2023 maximum (1.71 +- 0.03, no violation). WMAP/Planck give
A ~ 0.85, giving 2*sqrt(2)*0.85 = 2.40, exactly the 2025 "violation."
The "violation" turns on precisely when the fitted amplitude crosses
1/sqrt(2) ~ 0.707. But A*cos(alpha) with A > 0.707 **cannot be the pair
correlation of any +-1 (or [-1,1]) function on the sphere**: it violates
the tetrad constraint proved above. So the E(alpha) curve they plot
cannot be the pair average they claim to compute. Something in the chain
from pixels to E(alpha) to |C| does not do what the text describes.

### 3.5 What this does and does not show

- It shows the 2.42 cannot come from the described computation. The
  "sound violation" is an artifact.
- It does not identify the exact bug (their code is not public). Likely
  candidates: an error in the pair-selection or angle-binning code, a
  normalization error in P(X,Y,alpha), or the tetrad scan combining
  E values computed under inconsistent conventions.
- It does not rescue or refute the broader idea of cosmic Bell tests;
  it only shows this measurement cannot be what it claims to be.
- The authors' inequality proof itself (rotation-averaged CHSH) is
  correct. What fails is the measurement, not the theorem.

---

## 4. v1 -> v2 changes (2025-02-19 to 2025-09-04)

### 4.1 Numerical changes in the headline tables

**WMAP best tetrads (Table 2).** The reported maxima moved and the
angle sets changed completely:

- v1 best (neighbour averaging): (133, 49, 137, 139),
  Q/V/W = 2.4208 / 2.4222 / 2.4224, error 0.0009.
- v2 best (neighbour averaging): (133, 135, 137, 45),
  Q/V/W = 2.4248 / 2.4262 / 2.4263, error 0.0009.

Note: with E(alpha) ~ 0.85*cos(alpha), the v1-style tetrad
(133, 49, 137, 139) gives |C| ~ 0.85*(cos133 + cos49 + cos137 - cos139)
~ 0.85*(-0.002) ~ 0, not 2.42. The v2-style tetrad (133, 135, 137, 45)
gives ~2.40, matching the reported value. So v1's table either listed
the angles in an order inconsistent with Eq. (1) or the scan itself
changed; v2's tetrads are the ones that maximize |C| for a
cosine-shaped correlator. Either way, both versions report |C| > 2 and
both are subject to the impossibility in section 3.

**Planck (Table 3).**

- v1 COMMANDER top: ~2.1949 at (43, 135, 45, 47).
- v2 COMMANDER top: 2.1960 (neighbour averaging) / 2.1966 (spherical
  harmonics) at (133, 135, 137, 45).
- v1 NILC top: ~2.4353. v2 NILC top: 2.4376 / 2.4378.
- SEVEM and SMICA shift comparably; the full side-by-side table is in
  the working notes (/tmp/cmbbell/, from the extracted v1.txt/v2.txt).

**Claimed violation-angle ranges** (angles appearing in violating
tetrads) narrowed slightly:

- WMAP: v1 [45,51] U [129,141] -> v2 [45,49] U [131,141].
- Planck: v1 [39,49] U [131,141] -> v2 [41,49] U [133,141].

### 4.2 New in v2

- **Appendix B, "Unmixed Probabilities."** The v1 method ("mixed
  probabilities") assigns per-pixel probabilities from the temperature
  error interval; the unmixed method uses the raw sign of dT. v2 reports
  unmixed maxima ~1.5% larger (WMAP Q 2.4637, Planck COMMANDER 2.2321,
  at the same tetrad shape) and calls the mixed method "conservative."
  Both methods exceed 2, so both are subject to section 3.
- **Table 6: violation fractions.** Fraction of all scanned tetrads
  with |C| >= 2: WMAP ~3.3-3.4%, Planck NILC ~3.4%, SEVEM ~2.3%,
  SMICA ~1.8%, COMMANDER ~1.1%, for both degradation methods. These
  fractions are useful for the null test in section 7: an honest
  pipeline must give 0% on every sky.
- The random-pair sampling is now described as 10^7 pairs per bin,
  repeated 100 times, reporting the "most unfavorable" (lowest-maximum)
  repetition. Error bars on maxima remain ~0.0009-0.0010; these are
  propagated sampling errors, not global significances, and they do not
  account for maximizing over the tetrad scan.

### 4.3 What did not change

Section structure, the Eq. (1) inequality, the Nside=64 degradation, the
two degradation methods, the 2-degree bins, and the interpretation
(violation = failure of local realism) are essentially unchanged.

---

## 5. Prior literature

### 5.1 Dale, Lapiedra, Morales-Lladosa, Phys. Rev. D 107, 023506 (2023)
(arXiv:2302.05125)

- Proved the cosmic CHSH inequality used by the 2025 paper and applied
  it to COBE data with the hard dichotomy f = sign(dT).
- Scanned ~95 billion angle quadruples (871 angle values); maximum
  |C| = 1.71 +- 0.03 at (125, 39.2, 135, 139) degrees. No violation.
- This is exactly what section 3 predicts for honest data: the maximum
  over any scan is still <= 2. The 2023 paper is consistent with the
  theorem; the 2025 paper is not.
- The 2023 paper also derives new inequalities from temperature-
  polarization cross-correlations for WMAP (not pursued in 2025).

### 5.2 Martin and Vennin, Phys. Rev. D 96, 063501 (2017)
(arXiv:1706.05001), "Obstructions to Bell CMB experiments"

- Different approach: pseudo-spin operators on the inflationary
  two-mode squeezed state. Shows formal violations exist but practical
  obstructions (non-commuting observables that cannot actually be
  measured on the CMB, decoherence, cosmic variance, the impossibility
  of spacelike-separated measurement choices) are probably
  insurmountable.
- Relevant as context: the field's default stance is deep skepticism of
  cosmic Bell tests. But Martin-Vennin does not address Dale et al.'s
  specific dichotomized-temperature method; the refutation in section 3
  is independent and stronger for this paper.

### 5.3 What the v2 cites and does not answer

v2 cites the Bell-test literature but contains no Gaussian-simulation
null test: at no point do the authors run their full tetrad scan on
simulated Lambda-CDM skies to show the maximum stays below 2. That is
the missing control, and section 7 scopes it.

---

## 6. Method as described (for the reimplementation)

- **Data:** WMAP Q/V/W and Planck COMMANDER, NILC, SEVEM, SMICA
  temperature maps; monopole/dipole restored per Table 1, then maps
  degraded to Nside=64 (49,152 pixels) by neighbour averaging and by
  spherical-harmonic truncation.
- **Dichotomy (mixed):** per-pixel p+(x) = e+/(2*ea), p-(x) = 1 - p+(x),
  where [dT - ea, dT + ea] is the temperature error interval and e+ its
  positive part. f(x) = p+(x) - p-(x) in [-1, 1].
- **Correlator:** P(X,Y,alpha) = (1/N_alpha) sum over pairs of
  pX(x)*pY(y); E(alpha) = sum_{X,Y} X*Y*P(X,Y,alpha).
  Note: this equals the pair average of f(x)*f(y), so section 3 applies
  exactly as written.
- **Scan:** 90 bins of width 2 deg (6 deg < alpha <= 180 deg); 10^7
  random pairs per bin; 100 repetitions, least-favorable reported;
  tetrad scan maximizing |C|; errors by standard propagation.

---

## 7. Reproduction and null-test plan

### 7.1 Goal

Determine whether an honest implementation of the described pipeline
reproduces |C| > 2 on Planck data, and characterize the maximum-|C|
distribution under the Gaussian null.

### 7.2 Predicted outcome (from section 3)

An honest implementation **must** give max |C| <= 2 (up to ~0.001
binning noise) on real Planck maps, on Gaussian simulations, and on any
other sky. If it does, the published 2.42 is confirmed as an artifact of
the authors' code, and the residual task is identifying the breaking
step by differential testing (see 7.5).

### 7.3 Steps

1. **Data:** public Planck PR3 foreground-cleaned maps (start with
   SMICA), common mask, Nside=64 degradation by both methods
   (healpy: ud_grade for neighbour averaging; harmonic truncation via
   map2alm/alm2map for the spectral method).
2. **Dichotomy:** implement mixed probabilities exactly per section 3.5
   (needs per-pixel temperature errors; approximate first with uniform
   or hit-count-based errors, then refine) and unmixed (sign of dT).
3. **Correlator:** E(alpha) in 2-deg bins from pair averages. Use all
   pairs where feasible (49k pixels -> ~1.2B pairs; use random 10^7/bin
   as they do, plus one all-pairs run at lower Nside for validation).
4. **Tetrad scan:** replicate their scan over angle quadruples;
   record max |C|, the maximizing tetrad, and the Table-6 violation
   fraction (must be 0% if honest).
5. **Gaussian null:** 100+ Gaussian isotropic simulations from the
   Planck best-fit Cl (beam + pixel window + noise + mask +
   degradation), full pipeline on each; histogram of max |C| and of
   violation fractions.
6. **Diagnostics:** phase-randomized skies (same Cl, destroyed phases);
   sign-shuffled maps; half-mission difference maps (pure noise, should
   give |C| ~ small); pre-registered fixed tetrad evaluated without
   optimization.

### 7.4 Tractability

High. All data and tools are public (Planck Legacy Archive, healpy,
numpy). The main work is matching their preprocessing (monopole/dipole
handling, mask, degradation) and pair-sampling details. A low-Nside
prototype can validate the section-3 bound numerically before the full
run.

### 7.5 If the honest pipeline stays <= 2 (expected): finding their bug

Differential tests to isolate the breaking step: (a) compare our E(alpha)
curve to their Fig. 1 (if ours matches but our |C| <= 2, the bug is in
their tetrad scan/combination, not in E); (b) test whether their E(alpha)
satisfies the tetrad constraint at their own angles using digitized
values from their figures; (c) check sensitivity to pair-sampling
implementation (with/without replacement, per-bin vs global sampling).

---

## 8. Open questions and caveats

- The exact computational step that produces 2.42 is unidentified (code
  not public). Section 3 proves such a step exists; it does not name it.
- The analysis above takes the paper's equations at face value. If the
  code computes something materially different from Eqs. (3)-(6), the
  "bug" may be a description-code mismatch rather than a numerical slip.
- No claim is made here about cosmic Bell tests in general, only about
  this measurement.
- All arXiv PDFs and extracted text are in /tmp/cmbbell/ (ephemeral);
  the key extracted numbers are summarized above.

---

## 9. Files

- This scope: ~/workspace/vsl-planck/cmb-bell/SCOPE.md
- Ephemeral working material (PDFs, extracted text, check scripts):
  /tmp/cmbbell/ (v1.pdf, v2.pdf, prd2023.pdf, v1.txt, v2.txt,
  prd2023.txt). Re-extract from arXiv if needed; do not treat /tmp as
  durable.

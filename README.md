# Independent reimplementation of the cosmic CHSH analysis in arXiv:2502.13846

A clean-room reimplementation of the analysis in Dale et al., "Violation of
Bell inequalities from Cosmic Microwave Background data" (arXiv:2502.13846),
built to independently check the published CHSH values on Planck PR3 SMICA
data. Shared publicly at the authors' request so the two implementations
can be compared directly.

## Headline results (Planck PR3 SMICA, Nside=64)

- Near-full-sky |C| = 2.001 at the paper's best tetrad (133, 135, 137, 45),
  where the paper's Table 3 reports 2.2755 (neighbour averaging) /
  2.2764 (spherical harmonics).
- Optimized over ~3M realizable tetrads evaluated against a single
  precomputed correlation curve: maximum 2.027.
- With the Planck common confidence mask applied: |C| = 2.667 at the same
  tetrad. A pure-dipole control (no CMB, no noise) gives 2.669 masked and
  2.007 unmasked: masking alone is sufficient to produce the effect.

See `REPORT.md` for the full narrative and `SCOPE.md` for the analytic
argument (the |C| <= 2 bound for all-sky pair averages, and why masked pair
averages fall outside it). Machine-readable headline numbers are in
`results/rerun.json`.

## Contents

- `code/` — prep.py, pipeline.py, analyze.py, sims.py, consolidate.py,
  run_all.py, prototype.py, prototype2.py
- `REPORT.md` — full report of the rerun
- `SCOPE.md` — scope and analytic bound
- `results/rerun.json` — headline results

## Caveats (also documented in REPORT.md)

1. The Planck downloads used were slightly truncated: the SMICA file at
   99.7% of pixels (south polar cap excluded and treated as masked), the
   mask file at 97.3% (missing rows treated as masked).
2. The half-mission noise diagnostic was not run; the noise level uses a
   fallback (ea = 10 uK) validated by mixed/unmixed dichotomy agreement.
3. `prep.py` and `run_all.py` contain local path constants that must be
   adjusted before running. The Planck data files themselves are public
   Planck Legacy Archive downloads and are not included here.

Contact: Janet, AI research assistant (janetmariepm@gmail.com).

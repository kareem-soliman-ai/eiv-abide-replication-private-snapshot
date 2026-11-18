---
status: current
last_verified: 2026-11-12
verdict_evidence: "Plan written and hashed on 2026-11-12 before any ABIDE-II outcome (FD or IQ-FD relation) was read. Only the bucket layout and the header plus two rows of one participants.tsv had been inspected. Not externally registered."
---

# EXP-06 analysis plan: ABIDE-II replication (frozen before data access)

## Plain-English summary

ABIDE-I gave us three results: head-motion scores are reliable within a scan, the corrected method agrees with ordinary regression once noise is measured, and higher IQ goes with less motion after controls, though that last result was borderline. All three were found by exploring one dataset. This plan fixes, in advance, how the same questions will be tested on ABIDE-II, a separate set of sites and people, and what counts as a pass. It was written before looking at any ABIDE-II motion data. It is not registered with an outside registry, so it is weaker than a formal preregistration, but the file hash and the run receipt record the order of events.

## Sample

- ABIDE-II cross-sectional collections on the public bucket (`data/Projects/ABIDE2/RawData/ABIDEII-*`), excluding the two longitudinal collections (`UCLA_Long`, `UPSM_Long`) from the primary analysis because they include ABIDE-I participants.
- Include a participant when full-scale IQ is present and above 40, age, sex and diagnosis are present, and the fMRIPrep output has at least one resting-state run with 40 or more frames of framewise displacement (FD).
- ABIDE-II participant tables carry no rater quality flag, so none is applied. ABIDE-I checks showed the flag changes the slope by under 7%.
- Sites with fewer than 5 included participants are dropped.
- Outcome: mean FD of the first run of the first session. Primary scale: log FD. Raw FD is reported alongside.
- IQ error: the same age bands as ABIDE-I (4.0, 3.4, 3.0 points).

## Confirmatory tests and pass criteria

| ID | Claim from ABIDE-I | Test in ABIDE-II | Pass if |
|---|---|---|---|
| H1 | Run-mean FD is reliable within a scan | G coefficient of log FD, two contiguous halves, within site | G ≥ 0.80 |
| H2 | Higher IQ, less motion, after controls | OLS slope of log FD on IQ with site, age, age², sex and diagnosis partialled out; wild cluster bootstrap by site (Rademacher, null imposed, 4,999 draws, two-sided) | slope negative and p < 0.05 |
| H2b | Size is similar | Same slope | estimate lies inside the ABIDE-I 95% site-bootstrap interval [−0.0091, −0.0010] |
| H3 | With measured noise the corrected method agrees with OLS | Structural PCR-2 with per-subject noise (8-block batch means) against OLS, log scale, with controls | ratio of slopes between 0.90 and 1.10 |
| H4 | The v1 result is a weighting artefact that reproduces | Structural PCR-2 with site-SD noise against OLS, raw scale, no controls | OLS slope ÷ PCR-2 slope > 1.5 |
| H5 | A fit without a scatter term fails with realistic noise | Corrected York with measured noise, log scale, with controls | reduced χ² > 5 |

H2 is the primary test. If the OLS slope in H3 or H4 is within its own bootstrap noise of zero, the ratio is reported as not interpretable instead of pass or fail.

## Secondary, reported without pass criteria

- Diagnosis main effect on log FD with controls.
- IQ × diagnosis interaction on log FD (ABIDE-I: steeper in the autism group, p = 0.035, uncorrected).
- Reliability across runs and across sessions wherever a participant has more than one, including the longitudinal collections.
- Variance components for K = 2, 4, 8 and per-site reliability.
- Pooled ABIDE-I and ABIDE-II slope with controls.

## Fixed analysis details

Code: the functions in `experiments/EXP-04_FrameLevelGStudy_2026-11-12/scripts/run_gstudy.py` and `experiments/EXP-05_RobustnessChecks_2026-11-12/scripts/run_robustness.py`, unchanged. Seed 20261112. Site-cluster bootstrap with 1,000 replicates for intervals. No outlier removal. No transformation other than the natural log. Any departure from this plan is listed in the report under "Deviations".

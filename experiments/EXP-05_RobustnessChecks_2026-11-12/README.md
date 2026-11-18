---
status: current
last_verified: 2026-11-12
verdict_evidence: "Run 2026-11-12 (27.7 s, exit 0): eight robustness checks on ABIDE-I; 4,999 wild cluster draws; 300 simulations. results/summary.json holds hashes and seed."
---

# EXP-05: robustness checks before write-up

## Plain-English summary

Before writing up, we ran the checks a sceptical reviewer would ask for. Most results held. One did not hold as firmly as reported earlier: with a significance test suited to having only 19 sites, the IQ-motion link after controls is borderline (p about 0.04 to 0.05 on the log scale, 0.07 on the raw scale). The noise estimate does not depend on how the scan is cut up, excluding or including failed-quality scans changes little, and no single site drives the result. A simulation shows the corrected method is unbiased on the log scale and slightly biased on the raw scale when noise grows with motion.

## Objective

Close open questions from EXP-01 to EXP-04: sensitivity of the noise estimate, selection and truncation, linearity, inference with few clusters, single-site influence, raw-scale behaviour of PCR-2, and added controls.

## Methodology

ABIDE-I, n=935 (phenotypic FD) or 853 (frame-level FD). Controls throughout: site, age, age², sex, diagnosis. Wild cluster bootstrap by site with Rademacher weights and the null imposed, 4,999 draws. Simulation: 300 datasets of n=853, true log slope −0.006, person SD of log FD 0.60, noise proportional to level (coefficient of variation 0.12). Run: `python experiments/EXP-05_RobustnessChecks_2026-11-12/scripts/run_robustness.py`.

## Results

| Check | Result |
|---|---|
| R1 Noise from 4, 8, 16 blocks | PCR-2 log slope with controls −0.00357, −0.00363, −0.00371; raw −0.00073, −0.00071, −0.00072 |
| R2 Include quality-fail scans (n=1,029) | Log slope with controls −0.00425 against −0.00419; raw −0.00075 against −0.00070 |
| R3 Exclude mean FD above 0.5 mm (n=907) | Log −0.00374; raw −0.00056 |
| R3 Exclude mean FD above 0.2 mm (n=794) | Log −0.00185; raw −0.00015 |
| R4 Quadratic IQ term | p = 0.79 |
| R4 Slope by IQ third (log, controls) | low −0.0055, middle +0.0061, high −0.0015 |
| R5 Log slope with controls | −0.00419; cluster-robust p = 0.041; wild cluster p = 0.052; random-slope mixed model −0.00489, p = 0.041 |
| R5 Raw slope with controls | −0.00070; cluster-robust p = 0.070; wild cluster p = 0.070; mixed model −0.00091, p = 0.063 |
| R6 Leave one site out | Log −0.0053 to −0.0033; raw −0.00086 to −0.00047; negative in all 19 |
| R7 Simulation, log scale | PCR-2 recovers 99.9% of the true slope; OLS 94.2% |
| R7 Simulation, raw scale | PCR-2 recovers 88.1% of the target; OLS 95.5% |
| R8 Add IQ test type as a control | Log −0.00392; raw −0.00065 |
| R8 Site scan length against site mean FD | Spearman 0.16 (median frames per site 81 to 299) |

## Discussion

1. **The controlled association is borderline, not established.** Site-bootstrap percentile intervals in EXP-03 excluded zero, but tests built for few clusters put p between 0.04 and 0.07. EXP-06 tests it on new data. [EVIDENCE: R5]
2. **Excluding high-motion participants shrinks the slope sharply.** Cutting at 0.2 mm halves it or worse. Exclusion on the outcome is expected to do this, so studies that scrub high movers will see a weaker association. [EVIDENCE: R3]
3. **The slope is not smooth across the IQ range.** The quadratic term is null, but the thirds are erratic, which fits a weak effect. [EVIDENCE: R4]
4. **PCR-2 belongs on the log scale here.** Proportional noise breaks the constant-scatter assumption on the raw scale and costs about 12%; that is smaller than the raw-scale gap seen in EXP-04, so slope differences between high and low movers contribute too. [EVIDENCE: R7]

## Limits

Simulation parameters were chosen to resemble EXP-04, not fitted. No correction for the number of checks.

## Files

`results/summary.json`, `r2_r3_selection_truncation.csv`, `r6_leave_one_site_out.csv`, `r7_simulation.csv`. No subject identifiers; the script asserts this.

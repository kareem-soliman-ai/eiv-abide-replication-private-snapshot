---
status: current
last_verified: 2026-11-12
verdict_evidence: "Run 2026-11-12T04:05Z (61.7 s, exit 0) against analysis_plan.md, SHA-256 e27ce85f…a4cd, hashed 2026-11-12T03:58Z before ABIDE-II outcomes were read. n=460, 11 sites; 1,000 site-cluster bootstraps; 4,999 wild cluster draws. summary.json holds manifest and script hashes."
---

# EXP-06: ABIDE-II replication against a frozen plan

## Plain-English summary

We tested the ABIDE-I results on a separate dataset, ABIDE-II, using a plan and pass marks written down beforehand. The measurement results held: motion scores are again reliable within a scan (0.88), and a fit without a scatter term again fails badly. The headline association did not hold: after controlling for site, age, sex and diagnosis, the link between IQ and motion in ABIDE-II is small and statistically indistinguishable from zero (p = 0.68). Autistic participants move more than controls in both datasets, strongly. The apparent ABIDE-I finding that the IQ slope is steeper in the autism group reversed direction here, so it should be dropped. Only about half of the eligible ABIDE-II participants have the frame-level files publicly released, so this test has limited power.

## Objective

Confirmatory test of five claims derived from ABIDE-I (EXP-01 to EXP-05). Plan: [analysis_plan.md](analysis_plan.md).

## Methodology

As planned. Data: 16 ABIDE-II participant tables and 699 fMRIPrep confounds tables (65 MB) from the public bucket, stored on `D:` (see `data_pointer.txt`). Sample flow: 1,007 listed, 909 with IQ above 40, age, sex and diagnosis, 461 with frame data, 460 after dropping one site with fewer than 5. Eleven sites; KKI contributes 192 (42%). Mean age 12.8 years, mean IQ 111, 215 ASD and 245 TD. Code: the EXP-04 and EXP-05 functions, imported unchanged (hashes in `results/summary.json`).

## Results

### Confirmatory

| ID | Test | Result | Criterion | Verdict |
|---|---|---|---|---|
| H1 | G coefficient, log FD, two halves | 0.881 | ≥ 0.80 | **Pass** |
| H2 | Log slope with controls, wild cluster bootstrap | −0.00176, p = 0.68 | negative and p < 0.05 | **Fail** |
| H2b | Same slope inside the ABIDE-I interval | −0.00176 | within [−0.0091, −0.0010] | Pass (uninformative given H2) |
| H3 | PCR-2 with measured noise ÷ OLS | 1.013 | 0.90–1.10 | Not interpretable by the plan's rule (OLS interval includes zero); the ratio itself meets the criterion |
| H4 | OLS ÷ PCR-2 with site-SD noise, raw, no controls | 7.02 | > 1.5 | Not interpretable by the plan's rule |
| H5 | York without scatter, reduced χ² | 13.4 | > 5 | **Pass** |

### Variance components

G coefficient: raw 0.79, 0.85, 0.89 and log 0.88, 0.91, 0.93 for K = 2, 4, 8. Noise share of within-site variance: 11–21% raw, 7–12% log. Site share of total variance: 8% raw, 10% log.

### Slopes (95% site-bootstrap intervals)

| Estimator | Log, no controls | Log, controls |
|---|---|---|
| OLS | −0.0020 [−0.0082, +0.0037] | −0.0018 [−0.0105, +0.0033] |
| PCR-2, measured noise | −0.0022 [−0.0084, +0.0036] | −0.0018 [−0.0110, +0.0038] |
| PCR-2, site-SD noise | −0.0012 [−0.0064, +0.0038] | −0.0007 [−0.0071, +0.0032] |

Raw-scale slopes are in `results/slopes_with_intervals.csv`; none excludes zero.

### Secondary

- Diagnosis: TD participants have 0.37 lower log FD than ASD participants with controls (p < 10⁻⁸). ABIDE-I gave 0.23 (p < 0.001).
- IQ × diagnosis: coefficient −0.0056 (p = 0.045), opposite in sign to ABIDE-I (+0.0055, p = 0.035).
- Between-run reliability, 57 participants with 3 runs each at two sites (IP_1, OHSU_1), same session: ICC 0.42 raw, 0.60 log.
- Pooled ABIDE-I and ABIDE-II, 1,313 participants, 30 sites, with controls: log slope −0.0031, wild cluster p = 0.072; raw −0.00081, p = 0.20.

## Deviations from the plan

1. Participant tables were read as Latin-1 because one file is not valid UTF-8. No values are affected.
2. The longitudinal collections and ETHZ have no participant table on the bucket, so they were never listed; between-session reliability could not be estimated. Between-run reliability within a session is reported.

## Discussion

1. The measurement findings replicate: within-scan reliability and the failure of scatter-free fits. [EVIDENCE: H1, H5]
2. The IQ-motion association with controls does not replicate. The ABIDE-II interval is wide and contains the ABIDE-I estimate, so this is absence of confirmation, not evidence of no effect. Pooled over both datasets it is not significant (p = 0.072). [EVIDENCE: H2, pooled]
3. Agreement between PCR-2 with measured noise and OLS replicates numerically (ratio 1.013), though the plan's rule classes it as not interpretable because the slope itself is not distinguishable from zero. [EVIDENCE: H3]
4. The diagnosis effect on motion is the robust person-level result. [EVIDENCE: secondary]
5. The IQ × diagnosis interaction is unreliable and should not be reported as a finding. [EVIDENCE: sign reversal]

## Limits

Only 460 of 909 eligible participants have public fMRIPrep outputs, and which ones were processed is not random by site; KKI dominates; the sample is younger than ABIDE-I; 11 clusters give low power. No quality-rating filter exists for ABIDE-II. The plan was not externally registered.

## Files

`analysis_plan.md`, `analysis_plan.sha256`, `results/summary.json`, `g_variance_components.csv`, `g_per_site.csv`, `slopes_point.csv`, `slopes_with_intervals.csv`, `data_pointer.txt`. No subject identifiers; the script asserts this.

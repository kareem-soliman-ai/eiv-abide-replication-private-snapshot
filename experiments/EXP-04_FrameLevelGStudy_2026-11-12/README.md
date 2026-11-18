---
status: current
last_verified: 2026-11-12
verdict_evidence: "Full run 2026-11-12 (97.7 s, exit 0) on 853 of 935 subjects with frame-level FD from public fMRIPrep outputs; variance components for K=2,4,8 segments, per-site components, between-run ICC (23 subjects), 1,000 site-cluster bootstraps (1,000/1,000 successful). summary.json holds input, manifest and source hashes."
---

# EXP-04: G study of head motion from frame-level data, and reanalysis with measured noise

## Plain-English summary

Using frame-by-frame motion records, we measured how much of the spread in head motion is noise. Within a scan, the answer is 5–16%: the scan-average motion score is 84–95% reliable. The paper treated 100% of each site's spread as noise, and the measured noise is about a tenth of the value it used. Sites with a wide spread of motion have a wide spread of real differences between people, not noisier measurement. When the measured noise is fed into the updated PCR-2, the IQ-motion slope on the log scale is the same as ordinary least squares. The paper's claim that measurement uncertainty inflates the slope is therefore contradicted by direct measurement. One caveat: stability across separate runs is much lower (about 0.4 in the one site that has repeat runs), so a single scan is a weak measure of a person's general tendency to move.

## Objective

Run the G study (generalisability study: a split of score variance into person, occasion and error) that EXP-03 could not, and rerun the IQ-motion slope with measured, per-subject motion noise in place of the paper's assumed site spread.

## Methodology

- **Frame data:** fMRIPrep confounds tables from the public ABIDE bucket, one per run, holding framewise displacement (FD) at every frame. 894 files, 127 MB, stored on `D:` (see `data_pointer.txt`). `scripts/download_confounds.py` fetches them.
- **Coverage:** 853 of the 935 analysis subjects. USM loses 70 of 71 subjects (the release lacks them); 12 others are missing across six sites. Runs have 81–319 frames (median 179).
- **Outcome:** the fMRIPrep run-mean FD, so that outcome and noise come from the same frames. It correlates 0.956 with the phenotypic FD used in the paper. The two use different FD formulas, so levels differ (means 0.204 vs 0.127 mm).
- **G study:** design (persons within sites) × segments. Each run is cut into K contiguous segments, K = 2, 4, 8. Variance components come from the ANOVA mean squares with site × segment means removed; the site component from a REML model on person means. The G coefficient for the run mean is var(person) / (var(person) + var(person × segment) / K).
- **Between-run reliability:** ICC(1), the intraclass correlation for a single run, in the 23 OHSU subjects with 2–3 runs (64 runs).
- **Per-subject noise:** SD of 8 block means ÷ √8 (batch means, which allows for correlation between neighbouring frames). On the log scale, noise SD ÷ mean.
- **Reanalysis:** OLS; regression calibration (RC, corrects IQ error only); corrected York (functional, no scatter term); hazy York; structural PCR-2. Each is run with measured noise and with the paper's site-SD noise, on the same 853 subjects, with no controls and with site, age, sex and diagnosis partialled out, on raw and log FD. Site-cluster bootstrap, 1,000 replicates. USM's one remaining subject is given the median site SD for the site-SD variants.
- **Run:** `python experiments/EXP-04_FrameLevelGStudy_2026-11-12/scripts/run_gstudy.py --bootstrap 1000 --n-jobs 8`

## Results

### Variance components (pooled within site)

| Scale | Segments K | Person variance | Noise variance of run mean | G coefficient | Noise share of within-site variance |
|---|---:|---:|---:|---:|---:|
| Raw | 2 | 0.0433 | 0.0083 | 0.839 | 16.1% |
| Raw | 4 | 0.0454 | 0.0063 | 0.879 | 12.1% |
| Raw | 8 | 0.0470 | 0.0049 | 0.906 | 9.4% |
| Log | 2 | 0.351 | 0.0377 | 0.903 | 9.7% |
| Log | 4 | 0.342 | 0.0271 | 0.927 | 7.3% |
| Log | 8 | 0.331 | 0.0188 | 0.946 | 5.4% |

Split-half reliability with the Spearman-Brown correction: 0.87 raw, 0.91 log. Site accounts for 6% of total variance on the raw scale and 11–13% on the log scale. The G coefficient rises with K because motion drifts slowly within a run; K = 2 is the conservative figure.

### Per site (K = 4)

G coefficients run from 0.65 (OHSU) to 0.95 (SBL) on the raw scale and from 0.80 to 0.96 on the log scale. In every site the true person SD is close to the observed site SD. For example, UM_1 has site SD 0.431 with true person SD 0.409 and noise SD 0.135; NYU has 0.069, 0.064 and 0.027. Full table: `results/g2_per_site.csv`.

### Between runs

Single-run ICC in OHSU: 0.39 raw, 0.45 log (23 subjects, one site).

### Measured noise against the paper's assumption

The median measured noise SD is 0.013 mm, and the mean is 0.034 mm. The paper-style site SD averages 0.188 mm on this FD measure. The median ratio is 0.098. Measured noise rises with the subject's motion level (Spearman 0.84).

### IQ slope, 853 subjects

| Estimator | Raw, no controls | Raw, controls | Log, no controls | Log, controls |
|---|---:|---:|---:|---:|
| OLS | −0.00168 | −0.00083 | −0.0061 | −0.0037 |
| RC | −0.00179 | −0.00089 | −0.0065 | −0.0039 |
| PCR-2, measured noise | −0.00073 | −0.00071 | −0.0060 | −0.0036 |
| Hazy York, measured noise | −0.00092 | −0.00073 | −0.0060 | −0.0036 |
| York (no scatter), measured noise | −0.0118 | −0.0235 | −0.138 | −0.129 |
| PCR-2, paper-style site SD | −0.00032 | −0.00015 | −0.0039 | −0.0019 |
| OLS on the phenotypic FD, same subjects | −0.00119 | −0.00070 | −0.0068 | −0.0043 |

95% site-bootstrap intervals for PCR-2 with measured noise: raw [−0.00171, −0.00016] and [−0.00157, −0.00013] with controls; log [−0.0118, −0.0017] and [−0.0085, −0.0003] with controls. All exclude zero. The OLS interval with controls on the raw scale is [−0.00218, +0.00006]. Paper-style intervals cross zero in three of four cells. York without a scatter term has reduced χ² of 15–17, so it does not fit.

## Discussion

1. **The paper's noise assumption is wrong by about a factor of ten.** Measured within-scan noise is 5–16% of within-site variance. Site differences in spread are differences in real between-person variation. [GROUNDED]
2. **With measured noise, PCR-2 agrees with OLS on the log scale**, with and without controls. No inflation from measurement uncertainty exists there. [GROUNDED]
3. **On the raw scale PCR-2 still gives a flatter slope without controls** (−0.00073 against −0.00168). Noise grows with motion level, so high movers get less weight. This is a scale and weighting effect, not attenuation; with controls the gap nearly closes (−0.00071 against −0.00083). The log scale, where noise is roughly proportional, is the appropriate one. [GROUNDED for the numbers; PLAUSIBLE for the interpretation]
4. **A fit with no scatter term fails outright with realistic noise.** York's slope is 10–30 times too steep with reduced χ² near 16. This is the flaw the PCR update addressed, shown on real data: the structural PCR-2 handles it and the functional form does not. [GROUNDED]
5. **Two different reliabilities.** Within a scan, motion is measured well (0.84–0.95). As a measure of a person's general tendency to move, one scan is weak (between-run ICC about 0.4, one site, 23 subjects). That affects studies that treat motion as a trait. It does not bias a slope in which motion is the outcome. [GROUNDED for the estimates; the between-run figure is imprecise]

## Limits

- 82 subjects lack frame data, almost all from USM, so results describe 18 effective sites.
- fMRIPrep FD is not the same formula as the phenotypic FD (correlation 0.956).
- Segment-based noise captures variation over minutes within one session. It does not capture session-to-session variation or error in the realignment estimates themselves.
- The batch-means noise uses 8 blocks; block count affects it, as the K table shows.
- Between-run reliability rests on one site.
- Rows are treated as independent within site; no multiplicity correction.

## Next steps

A between-session reliability estimate needs a dataset with repeat sessions (ABIDE-II has longitudinal subsets). Not started.

## Files

`results/g2_variance_components.csv`, `g2_per_site.csv`, `g5_slopes.csv`, `g5_point_details.csv`, `g5_bootstrap_replicates.csv`, `summary.json`, and `data_pointer.txt`. None contains subject identifiers; the script asserts this.

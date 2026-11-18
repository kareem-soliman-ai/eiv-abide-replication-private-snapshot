---
status: current
last_verified: 2026-11-12
verdict_evidence: "Full run 2026-11-12 (117.4 s, exit 0): variance partition, noise bounds, IQ part-score reliability, five-rung control ladder x estimators on raw and log FD, 1,000 site-cluster bootstraps (1,000/1,000 successful where bootstrapped). summary.json holds input and source hashes."
---

# EXP-03: variance partition of head motion, noise bounds and a controlled reanalysis

## Plain-English summary

We asked how much of the spread in head motion is measurement noise and how much is real differences between people, then reran the IQ-motion analysis with proper controls. A full G study (generalisability study, which splits variance into person, occasion and error) is not possible here: each person has one scan, so person and error cannot be told apart directly. What the data do give is a ceiling on noise. Image-quality metrics computed independently of the motion estimate predict about two-thirds of the within-site motion variance at sites the model has not seen, so at most about a third can be estimation noise. The paper assumed all of it was. With site, age, sex and diagnosis controlled, the IQ-motion slope is still negative, but about 40% smaller than the pooled figure; diagnosis accounts for most of the drop after site. The measurement error that matters most turns out to be in IQ, not motion: IQ may be much less reliable than the paper assumed, which would make the true slope steeper, not flatter.

## Objective

1. Separate measurement noise from genuine between-person variance in mean framewise displacement (FD), as far as one scan per person allows.
2. Partition FD variance among site, person traits and residual.
3. Rerun the IQ-motion slope under controls, with noise levels set by the bounds instead of by site spread.

## Methodology

- **Data:** as EXP-01 and EXP-02: n=935, 19 sites, input SHA-256 `009f01a8…0865`. Raw FD and log FD.
- **P1, variance partition:** R² for site and for person traits (age, age², sex, diagnosis) in both entry orders; the IQ increment after all of them; a REML site variance component. This is the one-facet nested design "persons within sites", in which person and error share one term.
- **P2, noise bounds within site:** variance in FD that can be predicted cannot be FD estimation noise. Predictors: (a) person traits and IQ; (b) seven image-based quality metrics that are not computed from the realignment parameters (DVARS, outlier fraction, quality, ghost-to-signal, entropy focus, foreground-background energy, smoothness); (c) both. R² is cross-validated with whole sites held out (5 folds).
- **P3, IQ reliability:** verbal and performance IQ treated as two part-scores, combined with the Spearman-Brown formula (the reliability of a composite of two parts). Because the parts measure partly different abilities, this is a lower bound.
- **P4, controlled reanalysis:** five rungs of controls partialled out of both IQ and FD (the Frisch-Waugh method). Estimators: OLS; regression calibration (RC, which corrects for IQ error only); structural PCR-2 with σy² = (1 − reliability) × residual variance, for assumed FD reliabilities 0.5, 0.7 and 0.9, once with the same noise for everyone and once site-specific. Site-cluster bootstrap, 1,000 replicates; PCR-2 bootstrapped on rungs M0 and M3 only.
- **P5, diagnosis:** slopes within the autism (ASD) and typically developing (TD) groups, and the IQ × diagnosis interaction with site-clustered errors.
- **Not identifiable:** session-to-session variation in motion. No literature value was verified, so FD reliability is a swept assumption.
- **Run:** `python experiments/EXP-03_VariancePartitionControls_2026-11-12/scripts/run_partition.py --bootstrap 1000 --n-jobs 8`

## Results

### P1: where the motion variance sits

| Share of FD variance | Raw FD | Log FD |
|---|---:|---:|
| Site alone | 11.0% | 17.8% |
| Person traits alone (age, sex, diagnosis) | 4.9% | 6.4% |
| Person traits, added after site | 6.1% | 8.7% |
| IQ, added after site and person traits | 0.41% | 0.56% |
| Residual (person and error, confounded) | 82.5% | 73.0% |

Site spreads are unequal on both scales (Levene p < 10⁻⁷). The largest site SD is 10.3 times the smallest on the raw scale and 5.1 times on the log scale.

### P2: ceiling on FD measurement noise, within site

| Predictors | Held-out-site R², raw | Noise ceiling, raw | Held-out-site R², log | Noise ceiling, log |
|---|---:|---:|---:|---:|
| Person traits and IQ | 0.024 | 97.6% | 0.069 | 93.1% |
| Image quality metrics | 0.663 | 33.7% | 0.540 | 46.0% |
| Both | 0.659 | 34.1% | 0.547 | 45.3% |

At most 34% (raw) or 46% (log) of within-site FD variance can be estimation noise. In SD terms the raw ceiling is 0.086 mm, against the paper's mean σy of 0.129 mm.

### P3: IQ reliability

Verbal and performance IQ correlate 0.47 within site (n=778). Spearman-Brown gives 0.64 as a lower bound for full-scale IQ reliability, against the 0.945 implied by the paper's σx. The lower bound corresponds to σx ≈ 9.0 IQ points instead of 3–4.

### P4: IQ slope under controls

| Controls | Raw OLS | Raw RC | Raw 95% interval (OLS) | Log OLS | Log RC | Log 95% interval (OLS) |
|---|---:|---:|---|---:|---:|---|
| M0 none | −0.00125 | −0.00132 | [−0.00248, −0.00038] | −0.00695 | −0.00736 | [−0.0124, −0.0025] |
| M1 site | −0.00104 | −0.00111 | [−0.00222, −0.00035] | −0.00612 | −0.00650 | [−0.0123, −0.0023] |
| M2 + age, sex | −0.00096 | −0.00102 | [−0.00196, −0.00029] | −0.00556 | −0.00590 | [−0.0105, −0.0021] |
| M3 + diagnosis | −0.00070 | −0.00075 | [−0.00163, −0.00008] | −0.00419 | −0.00445 | [−0.0091, −0.0010] |
| M4 + eyes open/closed | −0.00070 | −0.00074 | [−0.00163, −0.00007] | −0.00418 | −0.00445 | [−0.0091, −0.0010] |

Structural PCR-2 at M3:

| FD noise model | Raw slope | Raw 95% interval | Log slope | Log 95% interval |
|---|---:|---|---:|---|
| Same noise for everyone, reliability 0.5, 0.7 or 0.9 | −0.00075 | [−0.00178, −0.00009] | −0.00449 | [−0.0099, −0.0010] |
| Site-specific, reliability 0.9 | −0.00064 | [−0.00168, −0.00006] | −0.00423 | [−0.0096, −0.0009] |
| Site-specific, reliability 0.7 | −0.00043 | [−0.00143, +0.00002] | −0.00372 | [−0.0090, −0.0007] |
| Site-specific, reliability 0.5 | −0.00031 | [−0.00119, +0.00008] | −0.00327 | [−0.0085, −0.0005] |

IQ-reliability sweep at M3 (RC): with σx = 3.5, 6 and 9.0 the log slope is −0.0045, −0.0051 and −0.0069; the raw slope is −0.00075, −0.00085 and −0.00116.

### P5: diagnosis

With site, age, sex and IQ held fixed, TD participants have 0.043 mm lower FD than ASD participants (p = 0.005), or 0.23 lower log FD (p < 0.001). The within-group slopes are ASD −0.00064 and TD −0.00053 on the raw scale (interaction p = 0.21), and ASD −0.0048 and TD −0.0022 on the log scale (interaction p = 0.035).

## Discussion

1. **A G study is not identifiable, but the paper's noise assumption is bounded out.** The paper set noise equal to the whole within-site spread. At least 54–66% of that spread is real motion in that scan, because independent image metrics predict it at unseen sites. [GROUNDED]
2. **The level of FD noise does not move the slope.** With the same noise for everyone, PCR-2 gives an identical slope at reliability 0.5, 0.7 and 0.9, and that slope equals RC. Only noise that differs by site moves it, and it does so through weighting. This repeats EXP-01 and EXP-02 under controls. [GROUNDED]
3. **Controls shrink the slope by 40–44%, and it remains negative.** Site removes 12–16% and diagnosis a further 20–27 percentage points. IQ accounts for about half a percent of FD variance once the rest is controlled. [GROUNDED]
4. **The consequential measurement error is in IQ.** If IQ reliability is closer to the part-score lower bound, the controlled slope is up to 1.5 times steeper than OLS. Ordinary errors-in-variables logic says the paper's direction of correction was backwards. [GROUNDED for the bound; PLAUSIBLE for where true reliability sits, since 0.64 is a lower bound and published full-scale reliabilities are higher]
5. **The association may be stronger in the autism group** on the log scale (p = 0.035, one test, not corrected for multiplicity). [PLAUSIBLE]

## Limits

- The image-metric bound concerns estimation noise for this scan. It says nothing about how stable a person's motion is across sessions.
- Image quality metrics are themselves degraded by motion; that is why they predict it, and it is also why the bound applies to estimation noise only.
- Rows are treated as independent within site; diagnosis groups are not matched; no multiplicity correction is applied across P5.
- Structural PCR-2 was bootstrapped on M0 and M3 only.

## Next steps

A real G study needs replicates. The ABIDE-I raw scans are public: recomputing per-frame FD and splitting each run into halves or odd/even segments would give a persons × segments design within sites, and an actual noise estimate for FD. That is a download of tens of gigabytes and a motion-correction run, and has not been started.

## Files

`results/p1_variance_partition.csv`, `p2_noise_bounds.csv`, `p4_controlled_slopes.csv`, `p4_bootstrap_replicates.csv`, `p5_diagnosis.csv`, `summary.json` (IQ reliability, sweep, hashes). None contains subject identifiers; the script asserts this.

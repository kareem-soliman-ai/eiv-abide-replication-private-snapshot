---
status: current
last_verified: 2026-11-12
verdict_evidence: "Full run 2026-11-12 (169.8 s, exit 0): 7 methods on n=935/19 sites, 8x8 grid, sigma_y-structure sensitivity, 1,000 site-cluster bootstraps; receipt.json holds input/source hashes. Legacy pipeline re-run first and its verification block passed."
---

# EXP-01: ABIDE rerun with the corrected EIV method

## Plain-English summary

The paper claims that the usual IQ-motion slope in ABIDE is about 4.7 times too steep once measurement uncertainty is accounted for. We reran it with the corrected method from Rasad. The correction to the PCR formula itself changes nothing here: the corrected fit gives the same slope as the original. But the rerun shows where the 4.7× comes from. It is produced entirely by giving high-motion sites less weight, because they have a larger "σy", the within-site spread of motion. Plain weighted least squares, with no measurement-error correction at all, gives the same answer. When the predictor's measurement error is corrected on its own, the slope gets slightly steeper than OLS, not 4.7× flatter. The headline claim therefore depends on treating real differences between people's motion as instrument noise. The corrected method warns against exactly that assumption.

## Objective

Decide whether the paper's headline (OLS −0.00125 vs PCR −0.00027, "bias factor" 4.67×) survives the corrected errors-in-variables methodology. Errors-in-variables (EIV) regression is regression that allows for measurement error in the predictor. The flaw was found in the 2026-09-22 numerical review in PCR-methodology and is fixed in Rasad (`docs/references/system/pcr-method-integration.md`). The PCR update step used combined-variance weights that are not the correct conditional step, and the model had no term for intrinsic scatter: genuine variation around the line that is not measurement error.

## Methodology

- **Data:** ABIDE-I `Phenotypic_V1_0b_preprocessed1.csv`, downloaded 2026-11-12 from the PCP S3 bucket (SHA-256 `009f01a8…0865`); the paper's filters give n=935 across 19 sites. σx (IQ uncertainty) is set by age band (4.0/3.4/3.0). σy is the within-site SD of mean FD (framewise displacement, the motion measure), exactly as in the paper.
- **Method core:** Rasad `src/eiv_core.py` at commit 302b394, `rasad-eiv/2`, vendored byte-identical as `src/eiv_core.py`. It is preferred to the PCR-methodology `structural_eiv.py` because it is the integrated, test-covered version (1,173 tests) and shares one implementation for the functional and structural fits.
- **Methods, all on identical inputs:** OLS; legacy PCR (the paper's rule); **corrected functional** EIV (fixed latent IQ values, no intrinsic scatter); **structural PCR-2**, which estimates a latent IQ distribution and the intrinsic scatter τ². Three diagnostics isolate the mechanism: WLS with weights 1/σy² (weighting with no EIV correction), regression calibration (the EIV correction without weighting) and site fixed effects (the within-site slope).
- **Checks:** leave-site-out (LOSO) R²; the paper's 8×8 σ-multiplier grid; a new **σy-structure** sensitivity (site-specific vs constant σy); and a **site-cluster bootstrap** (1,000 replicates, whole sites resampled, seed 20261112), which the paper did not report.
- **Run:** `python experiments/EXP-01_CorrectedEIVRerun_2026-11-12/scripts/rerun_corrected.py` from the repo root. Python 3.12.10.

## Results

### Headline, n=935

| Method | Slope (mm FD per IQ point) | OLS ÷ method | LOSO R² |
|---|---:|---:|---:|
| OLS | −0.001248 | 1.00 | −0.012 |
| Legacy PCR (paper) | −0.000265 | 4.70 | −0.074 |
| Corrected functional | −0.000265 | 4.71 | −0.074 |
| Structural PCR-2 (τ² = 0.00063, interior) | −0.000319 | 3.91 | −0.069 |
| WLS 1/σy², **no EIV correction** | −0.000245 | 5.08 | −0.075 |
| Regression calibration, **EIV correction only** (IQ reliability 0.945) | −0.001320 | 0.94 | −0.013 |
| Site fixed effects | −0.001042 | 1.20 | −0.012 |

All fits converged. The functional fit's reduced χ² is 1.21: by construction, the site-SD σy absorbs almost all of the residual spread.

### σy structure decides the answer

| σy assumption | Legacy PCR | Corrected functional | Structural PCR-2 |
|---|---:|---:|---:|
| Site-specific SD (paper) | −0.000265 | −0.000265 | −0.000319 |
| Constant, pooled within-site SD (0.149) | −0.001335 | −0.001336 | −0.001319 |
| Constant, 0.25 × pooled | −0.007359 | −0.008399 | −0.001319 |
| Constant, 0.05 × pooled | −0.051526 (not converged) | −0.087951 | −0.001319 |

When σy is the same at every site, all EIV methods agree with regression calibration: about −0.00132, slightly *steeper* than OLS. The two functional fits blow up as σy shrinks, because they have to explain residual scatter as measurement error. This is the omitted-intrinsic-scatter failure. Structural PCR-2 is stable because it moves that scatter into τ².

### Site-cluster bootstrap, 1,000/1,000 successful for every method

| Method | 95% interval for slope | Share of replicates with negative slope |
|---|---|---:|
| OLS | [−0.00234, −0.00028] | 0.997 |
| Legacy PCR | [−0.00097, +0.00004] | 0.949 |
| Structural PCR-2 | [−0.00114, +0.00003] | 0.961 |
| Site fixed effects | [−0.00210, −0.00027] | 1.000 |

The OLS-minus-PCR difference interval is [−0.00178, −0.00012]. The difference between the two slopes is therefore real under site resampling, but what produces it is the weighting.

### The paper's 8×8 grid

Constant multipliers on σy leave the relative site weights unchanged, so the grid cannot test the assumption that matters. Even so, legacy PCR fails to converge in 4 of the 64 cells, and in 8 cells (σy ≤ 0.75× with σx ≥ 1×) the conclusion reverses: PCR comes out 8–55× *steeper* than OLS. Structural PCR-2 converges in all 64 cells, with OLS ÷ method between 0.58 and 5.08.

### Within tiers

Across all four motion tiers the correction changes little: legacy and functional agree to three decimal places, and structural is within 3%. The paper's tier pattern of steeper slopes in higher-motion tiers is an OLS result and is unaffected.

## Discussion

1. **The formula flaw does not affect ABIDE numerically.** With b ≈ −0.001 and σx ≈ 3, the term b²σx² (≈10⁻⁵) is negligible next to σy² (≈0.017), so the old and corrected update steps coincide.
2. **The headline "measurement-uncertainty" claim is not supported.** The 4.7× ratio is reproduced exactly by WLS with weights 1/σy², which has no measurement-error model. Correcting IQ measurement error alone makes the slope slightly steeper (reliability 0.945, so about 6%). The paper's mechanism is down-weighting of high-motion sites. That holds only if each site's FD spread is instrument noise, rather than the genuine between-subject variation in motion it mostly is. Rasad's integration note states this boundary explicitly: unmodelled variation cannot be reassigned to instrument noise.
3. **What survives:** (a) the pooled OLS slope is fragile under site resampling, though its interval still excludes zero; (b) the site-level descriptive pattern (Figure 1, tier slopes); (c) site fixed effects shrink the slope by about 17%. Point (c) is modest site confounding, not 4.7×.
4. **Limits:** no direct estimate of FD measurement error exists in these data, so the "true" σy is unknown. Every σy here is an assumption, including the paper's. The bootstrap treats sites as exchangeable, ignores within-site dependence beyond clustering, and does not use survey weights.

## Stress tests (added 2026-11-12, `scripts/stress_tests.py`, `results/stress_tests.json`)

Rerun first: the legacy pipeline and EXP-01 reproduced identically, with only the runtime line differing.

| Test | Result | Reading |
|---|---|---|
| S1: simulate from the paper's own model (common slope; σy is pure outcome noise) | Mean OLS ÷ PCR = 0.96 when the true slope is the OLS value, 1.00 when it is the PCR value (500 simulations each) | Under the paper's assumptions OLS cannot come out about 4.7× steeper. Noise in the outcome does not bias OLS. The observed gap needs a different mechanism. |
| S2: slope heterogeneity | Site slopes steepen with site mean FD: Spearman ρ = −0.52 (p = 0.023); FIQ × site-mean-FD interaction p = 0.006 with site-clustered errors | That mechanism is genuine heterogeneity, which the weighting then treats as noise |
| S3: log(FD) | OLS ÷ PCR falls from 4.7 to 1.76 | Much of the heterogeneity is scale: FD skewness is 3.87 (0.72 after logging), and site SD tracks site mean (ρ = 0.81, p < 0.001; S6) |
| S4: robust and multilevel fits | Huber −0.00054; median −0.00044; random-slope mixed model −0.00123 (p = 0.005); log-scale mixed −0.0072 (p = 0.002) | The association survives multilevel modelling. Robust fits flatten it about 2.5×, which is a tail effect, not measurement error |
| S5: site jackknife | OLS ÷ PCR between 2.56 (NYU dropped) and 6.13 (YALE dropped) | Unstable, and dependent on single sites |

Verdict from the stress tests: the title's mechanism ("measurement uncertainty") is falsified by S1 together with the regression-calibration result. What survives is a modest, defensible claim. The IQ-motion association in ABIDE-I survives measurement-error correction and multilevel modelling (about 0.7% lower FD per IQ point on the log scale). On the raw scale, its size is dominated by right-skewed FD and high-motion sites.

## Recommendation

Withdraw or substantially correct the headline. The corrected method does not reproduce the claim that measurement uncertainty inflates the confound. It reproduces a heteroscedastic weighting effect under an untested σy assumption. A defensible revision would describe the result as a site-weighting sensitivity rather than an EIV correction, report the site-bootstrap intervals and the σy-structure table, and drop the "bias factor" framing. This is a public-claim decision for Kareem; no README, paper text or citation metadata was changed by this experiment.

## Files

`results/headline_by_method.csv`, `tier_slopes_by_method.csv`, `sensitivity_grid_by_method.csv`, `sigma_y_structure_sensitivity.csv`, `site_bootstrap_intervals.csv`, `site_bootstrap_replicates.csv`, `receipt.json` (hashes, versions, runtime). None contains subject identifiers; the script asserts this.

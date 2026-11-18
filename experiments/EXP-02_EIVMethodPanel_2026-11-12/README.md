---
status: current
last_verified: 2026-11-12
verdict_evidence: "Full run 2026-11-12 (141.8 s, exit 0): 13 estimators x 3 sigma specifications on n=935/19 sites; LOSO R^2; 1,000 site-cluster bootstraps per spec, 1,000/1,000 successful for every bootstrapped method; receipt.json hashes inputs and vendored sources."
---

# EXP-02: full EIV method panel on ABIDE-I

## Plain-English summary

We fitted the updated PCR and eleven other estimators to the same ABIDE data, including the standard errors-in-variables (EIV) methods. EIV methods are regression methods that allow for measurement error in the predictor. The methods do not really disagree: they split into groups according to a single input assumption. That assumption is whether each site's spread in head motion counts as measurement noise. If it does, every method that uses it gives the paper's flat slope. If it doesn't, all twelve agree within 8% on a slope slightly steeper than ordinary least squares. The updated PCR gives the same numbers as York regression, and PCR-2 gives the same numbers as an independent structural likelihood fit. On this dataset, PCR adds no information beyond the established methods. The paper's 4.7× is therefore a consequence of one assumption, not something any method found.

## Objective

Compare the updated PCR (Rasad `eiv_core`, functional and structural PCR-2) against the established EIV estimators on identical inputs. Find out whether the paper's conclusion depends on the method or on the inputs.

## Methodology

- **Data:** as in EXP-01: ABIDE-I phenotypic, n=935, 19 sites, input SHA-256 `009f01a8…0865`.
- **Estimators (13):**
  - naive: OLS;
  - correcting IQ measurement error only: regression calibration and SIMEX;
  - Deming, with constant error ratio δ = mean σy² / mean σx²;
  - WLS with weights 1/σy² (weighting only, not EIV);
  - functional fits without intrinsic scatter: corrected York, the objective-coherent coordinate descent, legacy PCR and updated PCR functional;
  - fits with intrinsic scatter: hazy York (York with τ² set so that reduced χ² = 1), the structural ML reference, and updated PCR-2 structural;
  - a random-slope mixed model as a non-EIV reference.
- **Sources:** `eiv_methods.py` and `pcr2_reference.py` are vendored byte-identical from the PCR-methodology working tree. That repo's suite of 13 tests passes, and these files are uncommitted there at e858adc. `src/eiv_core.py` comes from Rasad 302b394.
- **Specifications:**
  - **A**, the paper's: σx by age band, σy = within-site SD of mean FD;
  - **B**: constant σy = pooled within-site SD;
  - **C**: log(FD) with σy = within-site SD of log FD.
- **Checks:** leave-site-out R², with σ recomputed on the training sites; site-cluster bootstrap with 1,000 replicates per spec and SIMEX at 50 simulations per bootstrap replicate. The mixed model is not bootstrapped.
- **Run:** `python experiments/EXP-02_EIVMethodPanel_2026-11-12/scripts/run_panel.py --bootstrap 1000 --n-jobs 8`

## Results

### Spec A (paper σ): three clusters

| Cluster | Methods | Slope | OLS ÷ method | LOSO R² | 95% site-bootstrap interval |
|---|---|---:|---:|---:|---|
| σy not used as a per-point weight | OLS, RC, SIMEX, Deming, mixed model | −0.00123 to −0.00133 | 0.94–1.02 | −0.004 to −0.013 | about [−0.0026, −0.0004]; 99.6% negative |
| Site σy as the whole residual variance | WLS, York, coordinate descent, legacy PCR, **updated PCR** | −0.000245 to −0.000265 | 4.7–5.1 | −0.074 | about [−0.0010, +0.00002]; 96% negative |
| Site σy plus a common τ² | hazy York, structural ML, **PCR-2** | −0.000319 to −0.000367 | 3.4–3.9 | −0.061 to −0.072 | about [−0.0012, ±0.00001]; 97–98% negative |

### Spec B (constant σy): agreement

All twelve EIV and naive estimators give between −0.001248 and −0.001336, a spread of 7%, with intervals of about [−0.0027, −0.0004]. Every method that corrects for IQ measurement error lands 5–7% steeper than OLS. That is the size of correction implied by the IQ reliability of 0.945.

### Spec C (log FD): smaller gap, same structure

The cluster that ignores σy gives −0.0070 to −0.0074 per IQ point (about 0.7% lower FD per point). The weighted cluster gives −0.0036 to −0.0039, and the scatter-aware cluster −0.0047. The gap shrinks from 4.7× to at most 1.9× and comes from the same assumption.

### Method agreement

- Updated PCR functional, corrected York and coordinate descent agree to six significant figures in every spec. Legacy PCR differs in the fourth significant figure.
- PCR-2 structural and the independent structural ML reference agree to three or four significant figures, including τ².
- No fit failed, and every bootstrapped method succeeded in 1,000/1,000 replicates.

## Discussion

1. **The method makes no difference; the σy assumption decides everything.** Within any one specification the estimators agree. Between specifications A and B, only the methods that use site-specific σy as per-point noise change. [GROUNDED]
2. **The genuine measurement-error correction is small and goes the other way.** Correcting IQ error steepens the slope by about 6%; nothing in the panel flattens it except down-weighting high-motion sites. [GROUNDED]
3. **Updated PCR is York, and PCR-2 is standard structural ML, on these data.** That checks the Rasad implementation independently. It also means ABIDE provides no evidence for a PCR-specific advantage. [GROUNDED]
4. **Prediction prefers the methods that don't down-weight.** Out of site, R² is −0.013 for OLS against −0.074 for PCR and −0.004 for the mixed model. All values are negative because between-site differences in average motion dominate, so IQ adds little prediction at a new site. [GROUNDED]
5. **Scatter-aware fits sit in between, not with OLS.** τ² is modelled as constant, so the site-specific σy still drives the weights. A defensible PCR-2 analysis needs either an external estimate of FD measurement error or a site-varying scatter model. [PLAUSIBLE]

## Limits

The true measurement error in FD is unknown; every σy here, including the paper's, is an assumption. The bootstrap resamples sites only. No survey weights are used, and diagnosis groups are pooled, as in the paper.

## Files

`results/panel_point_estimates.csv`, `panel_bootstrap_intervals.csv`, `panel_bootstrap_replicates.csv`, `receipt.json`. None contains subject identifiers; the script asserts this.

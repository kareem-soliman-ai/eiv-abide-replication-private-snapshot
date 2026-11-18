"""EXP-03: variance partition of head motion, reliability bounds, controlled reanalysis.

P1 Variance partition of FD (raw, log): site / person covariates / residual,
   in both entry orders, plus REML site variance component (one-facet nested
   G-study design p:s; person and error are confounded in the residual).
P2 Bounds on FD measurement noise within sites:
   a) share predicted by stable person traits (age, sex, DX, FIQ)
   b) share predicted (5-fold CV) by image-based QC metrics that are not
      computed from the realignment parameters (DVARS, outlier fraction,
      quality, GSR, EFC, FBER, FWHM)
   c) both together. Predictable variance cannot be FD estimation noise.
   Occasion (session-to-session) variance is NOT identifiable: one scan each.
P3 IQ reliability check from part-scores (VIQ, PIQ): Spearman-Brown.
P4 Controlled reanalysis: IQ slope under a ladder of controls x estimator
   (OLS, regression calibration, structural PCR-2 with sigma_y from an assumed
   reliability, homoscedastic and site-specific), raw and log FD, with
   site-cluster bootstrap. Controls are partialled out (Frisch-Waugh).
P5 Diagnosis: slope within ASD and TD, and IQ x DX interaction.

Run from repo root:
  python experiments/EXP-03_VariancePartitionControls_2026-11-12/scripts/run_partition.py --bootstrap 1000 --n-jobs 8
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from joblib import Parallel, delayed
from scipy.stats import levene

EXP = Path(__file__).resolve().parents[1]
ROOT = EXP.parents[1]
sys.path.insert(0, str(ROOT / "src"))
warnings.filterwarnings("ignore")
from eiv_core import fit_eiv, measurement_covariances  # noqa: E402
from utils import filter_abide_sample, load_abide_csv, resolve_input_path, sigma_x_from_age  # noqa: E402

QC = ["func_dvars", "func_outlier", "func_quality", "func_gsr", "func_efc", "func_fber", "func_fwhm"]
LADDER = {
    "M0 none": [],
    "M1 site": ["site"],
    "M2 site+age+sex": ["site", "age", "sex"],
    "M3 site+age+sex+DX": ["site", "age", "sex", "dx"],
    "M4 M3+eyes": ["site", "age", "sex", "dx", "eyes"],
}
RELIABILITIES = [0.5, 0.7, 0.9]


def design(df: pd.DataFrame, controls: list[str]) -> np.ndarray:
    cols = [np.ones(len(df))]
    if "site" in controls:
        cols.append(pd.get_dummies(df["SITE_ID"], drop_first=True).to_numpy(float))
    if "age" in controls:
        a = df["AGE_AT_SCAN"].to_numpy(float)
        cols += [a, a**2]
    if "sex" in controls:
        cols.append((df["SEX"] == 2).to_numpy(float))
    if "dx" in controls:
        cols.append((df["DX_GROUP"] == 1).to_numpy(float))
    if "eyes" in controls:
        cols.append((df["EYE_STATUS_AT_SCAN"] == 2).to_numpy(float))
    return np.column_stack(cols)


def resid(v: np.ndarray, Z: np.ndarray) -> np.ndarray:
    return v - Z @ np.linalg.lstsq(Z, v, rcond=None)[0]


def r2(y: np.ndarray, Z: np.ndarray) -> float:
    return 1 - resid(y, Z).var() / y.var()


def slopes(df: pd.DataFrame, ycol: str, controls: list[str], pcr2: bool = True) -> dict[str, float]:
    """IQ slope after partialling controls from both IQ and FD (Frisch-Waugh)."""
    Z = design(df, controls)
    n, k = Z.shape
    x, y = resid(df["FIQ"].to_numpy(float), Z), resid(df[ycol].to_numpy(float), Z)
    sx = df["sigma_x"].to_numpy(float)
    sxx, sxy = x @ x / (n - k), x @ y / (n - k)
    out = {"ols": sxy / sxx, "rc": sxy / (sxx - np.mean(sx**2)),
           "iq_reliability_after_controls": 1 - np.mean(sx**2) / sxx}
    if pcr2:
        res_var = float(np.var(y - out["ols"] * x))
        site_var = pd.Series(y).groupby(df["SITE_ID"].to_numpy()).transform("var").to_numpy()
        for rel in RELIABILITIES:
            for label, sy2 in (("homo", np.full(n, (1 - rel) * res_var)), ("site", (1 - rel) * site_var)):
                try:
                    X, yy, C = measurement_covariances(x, y, sigma_X=sx, sigma_y=np.sqrt(sy2))
                    f = fit_eiv(X, yy, C, method="structural", max_iter=2000)
                    out[f"pcr2_{label}_rel{rel}"] = float(f.coef[0]) if f.converged else np.nan
                except Exception:
                    out[f"pcr2_{label}_rel{rel}"] = np.nan
    return out


def boot_rep(df: pd.DataFrame, sites: np.ndarray, seed: int) -> list[dict]:
    rng = np.random.default_rng(seed)
    parts = []
    for j, s in enumerate(rng.choice(sites, len(sites), replace=True)):
        g = df[df.SITE_ID == s].copy()
        g["SITE_ID"] = f"{s}#{j}"
        parts.append(g)
    b = pd.concat(parts, ignore_index=True)
    rows = []
    for ycol in ("fd", "log_fd"):
        for name, controls in LADDER.items():
            # PCR-2 only on the two rungs that matter, to keep the bootstrap tractable
            s = slopes(b, ycol, controls, pcr2=name in ("M0 none", "M3 site+age+sex+DX"))
            rows.append({"seed": seed, "scale": ycol, "model": name, **s})
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input")
    ap.add_argument("--bootstrap", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=20261112)
    ap.add_argument("--n-jobs", type=int, default=1, help="parallel workers (8 recommended here)")
    args = ap.parse_args()
    t0 = time.time()
    out = EXP / "results"
    path = resolve_input_path(args.input)
    df = filter_abide_sample(load_abide_csv(path)).reset_index(drop=True)
    assert len(df) == 935 and df.SITE_ID.nunique() == 19
    df["fd"] = df["func_mean_fd"]
    df["log_fd"] = np.log(df["fd"])
    df["sigma_x"] = sigma_x_from_age(df["AGE_AT_SCAN"])
    for c in ("VIQ", "PIQ"):
        v = pd.to_numeric(df[c], errors="coerce")
        df[c] = v.where(v > 40)
    res: dict = {"n": len(df), "sites": int(df.SITE_ID.nunique())}

    # ---- P1 variance partition
    person = ["age", "sex", "dx"]
    p1 = []
    for ycol in ("fd", "log_fd"):
        y = df[ycol].to_numpy(float)
        r_site, r_person = r2(y, design(df, ["site"])), r2(y, design(df, person))
        r_both = r2(y, design(df, ["site"] + person))
        Zx = np.column_stack([design(df, ["site"] + person), df["FIQ"].to_numpy(float)])
        r_all = r2(y, Zx)
        mm = smf.mixedlm(f"{ycol} ~ 1", df, groups=df["SITE_ID"]).fit(reml=True)
        vs, ve = float(mm.cov_re.iloc[0, 0]), float(mm.scale)
        groups = [g[ycol].to_numpy() for _, g in df.groupby("SITE_ID")]
        p1.append({"scale": ycol, "site_alone": r_site, "person_traits_alone": r_person,
                   "site_then_person_increment": r_both - r_site,
                   "person_then_site_increment": r_both - r_person,
                   "iq_increment_after_all": r_all - r_both, "residual": 1 - r_all,
                   "reml_site_variance_share": vs / (vs + ve),
                   "reml_var_site": vs, "reml_var_person_within_site_plus_error": ve,
                   "levene_p_equal_site_variances": float(levene(*groups, center="median").pvalue),
                   "site_sd_ratio_max_over_min": float(max(g.std() for g in groups) / min(g.std() for g in groups))})
    p1 = pd.DataFrame(p1)
    p1.to_csv(out / "p1_variance_partition.csv", index=False)

    # ---- P2 noise bounds, within site
    from sklearn.linear_model import LinearRegression
    from sklearn.model_selection import GroupKFold, cross_val_predict
    p2 = []
    Zs = design(df, ["site"])
    feats = {
        "a_person_traits": np.column_stack([design(df, person)[:, 1:], df["FIQ"]]),
        "b_image_qc": df[QC].to_numpy(float),
    }
    feats["c_both"] = np.column_stack([feats["a_person_traits"], feats["b_image_qc"]])
    for ycol in ("fd", "log_fd"):
        yw = resid(df[ycol].to_numpy(float), Zs)
        for name, F in feats.items():
            Fw = np.column_stack([resid(F[:, j].astype(float), Zs) for j in range(F.shape[1])])
            insample = 1 - resid(yw, np.column_stack([np.ones(len(yw)), Fw])).var() / yw.var()
            # leave-sites-out CV so the bound is not optimistic
            pred = cross_val_predict(LinearRegression(), Fw, yw, groups=df["SITE_ID"], cv=GroupKFold(5))
            cv = 1 - np.mean((yw - pred) ** 2) / yw.var()
            p2.append({"scale": ycol, "predictors": name, "within_site_r2_insample": insample,
                       "within_site_r2_cv_grouped": cv,
                       "noise_share_upper_bound": 1 - cv,
                       "noise_sd_upper_bound": float(np.sqrt((1 - cv) * yw.var()))})
    p2 = pd.DataFrame(p2)
    p2.to_csv(out / "p2_noise_bounds.csv", index=False)

    # ---- P3 IQ reliability from part scores
    iq = df.dropna(subset=["VIQ", "PIQ"])
    Zi = design(iq, ["site"])
    r_within = float(np.corrcoef(resid(iq.VIQ.to_numpy(float), Zi), resid(iq.PIQ.to_numpy(float), Zi))[0, 1])
    r_pooled = float(iq[["VIQ", "PIQ"]].corr().iloc[0, 1])
    res["P3_iq"] = {"n_with_viq_piq": len(iq), "viq_piq_r_pooled": r_pooled, "viq_piq_r_within_site": r_within,
                    "spearman_brown_pooled": 2 * r_pooled / (1 + r_pooled),
                    "spearman_brown_within_site": 2 * r_within / (1 + r_within),
                    "assumed_reliability_from_sigma_x": float(1 - np.mean(df.sigma_x**2) / df.FIQ.var()),
                    "fiq_sd": float(df.FIQ.std()),
                    "implied_sigma_x_if_sb_within": float(df.FIQ.std() * np.sqrt(1 - 2 * r_within / (1 + r_within)))}

    # ---- P4 controlled reanalysis: point estimates
    p4 = []
    for ycol in ("fd", "log_fd"):
        for name, controls in LADDER.items():
            p4.append({"scale": ycol, "model": name, **slopes(df, ycol, controls)})
    p4 = pd.DataFrame(p4)

    sites = np.array(sorted(df.SITE_ID.unique()))
    seeds = np.random.SeedSequence(args.seed).generate_state(args.bootstrap)
    keep = ["SITE_ID", "FIQ", "fd", "log_fd", "sigma_x", "AGE_AT_SCAN", "SEX", "DX_GROUP", "EYE_STATUS_AT_SCAN"]
    reps = Parallel(n_jobs=args.n_jobs, backend="loky")(delayed(boot_rep)(df[keep], sites, int(s)) for s in seeds)
    boot = pd.DataFrame([r for rep in reps for r in rep])
    boot.to_csv(out / "p4_bootstrap_replicates.csv", index=False)
    est_cols = [c for c in p4.columns if c not in ("scale", "model", "iq_reliability_after_controls")]
    long = []
    for _, row in p4.iterrows():
        b = boot[(boot.scale == row.scale) & (boot.model == row.model)]
        for e in est_cols:
            s = b[e].dropna() if e in b else pd.Series(dtype=float)
            long.append({"scale": row.scale, "model": row.model, "estimator": e, "slope": row[e],
                         "iq_reliability_after_controls": row.iq_reliability_after_controls,
                         "ci_low": s.quantile(.025) if len(s) else np.nan,
                         "ci_high": s.quantile(.975) if len(s) else np.nan,
                         "share_negative": float((s < 0).mean()) if len(s) else np.nan,
                         "n_boot_ok": int(len(s)), "n_boot": int(len(b))})
    p4long = pd.DataFrame(long)
    p4long.to_csv(out / "p4_controlled_slopes.csv", index=False)

    # ---- P5 diagnosis
    p5 = []
    for ycol in ("fd", "log_fd"):
        for dx, lab in ((1, "ASD"), (2, "TD")):
            g = df[df.DX_GROUP == dx]
            s = slopes(g, ycol, ["site", "age", "sex"], pcr2=False)
            p5.append({"scale": ycol, "group": lab, "n": len(g), "ols": s["ols"], "rc": s["rc"]})
        m = smf.ols(f"{ycol} ~ FIQ * C(DX_GROUP) + AGE_AT_SCAN + I(AGE_AT_SCAN**2) + C(SEX) + C(SITE_ID)", df).fit(
            cov_type="cluster", cov_kwds={"groups": df.SITE_ID})
        k = [c for c in m.params.index if c.startswith("FIQ:")][0]
        dxk = [c for c in m.params.index if c.startswith("C(DX_GROUP)") and "FIQ" not in c][0]
        p5.append({"scale": ycol, "group": "interaction FIQ x DX (TD vs ASD)", "n": len(df),
                   "ols": float(m.params[k]), "p_cluster": float(m.pvalues[k])})
        m2 = smf.ols(f"{ycol} ~ FIQ + C(DX_GROUP) + AGE_AT_SCAN + I(AGE_AT_SCAN**2) + C(SEX) + C(SITE_ID)", df).fit(
            cov_type="cluster", cov_kwds={"groups": df.SITE_ID})
        dxk = [c for c in m2.params.index if c.startswith("C(DX_GROUP)")][0]
        p5.append({"scale": ycol, "group": "DX main effect (TD minus ASD), M3", "n": len(df),
                   "ols": float(m2.params[dxk]), "p_cluster": float(m2.pvalues[dxk])})
    p5 = pd.DataFrame(p5)
    p5.to_csv(out / "p5_diagnosis.csv", index=False)

    for p in out.glob("*.csv"):
        assert not set(pd.read_csv(p, nrows=0).columns) & {"SUB_ID", "FILE_ID", "subject"}
    h = lambda q: hashlib.sha256(Path(q).read_bytes()).hexdigest()  # noqa: E731
    res["receipt"] = {"input_csv": h(path), "eiv_core.py": h(ROOT / "src/eiv_core.py"), "script": h(__file__),
                      "bootstrap": {"reps": args.bootstrap, "seed": args.seed, "unit": "site"},
                      "runtime_s": round(time.time() - t0, 1)}
    (out / "summary.json").write_text(json.dumps(res, indent=2))
    with pd.option_context("display.width", 250, "display.max_rows", 200, "display.max_columns", 30,
                           "display.float_format", "{:.6g}".format):
        print(p1.T.to_string()); print(p2.to_string(index=False)); print(json.dumps(res["P3_iq"], indent=1))
        print(p4long.drop(columns=["n_boot"]).to_string(index=False)); print(p5.to_string(index=False))
    print(f"runtime {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()

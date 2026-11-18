r"""EXP-05: robustness checks on the ABIDE-I programme before write-up.

R1 Noise-estimate sensitivity: per-subject noise from 4, 8, 16 blocks -> PCR-2 slope.
R2 Selection: include QC-fail subjects (phenotypic FD, n=1,029).
R3 Outcome truncation: exclude mean FD > 0.5 and > 0.2 mm.
R4 Linearity: quadratic IQ term; slope by IQ tercile.
R5 Inference with 19 clusters: mixed model with controls, cluster-robust SE,
   wild cluster bootstrap (Rademacher, null imposed).
R6 Leave-one-site-out range of the controlled slope.
R7 Simulation: is raw-scale PCR-2 biased when noise is proportional to level?
R8 Scan length by site, and IQ test type as an added control.

Run from repo root:
  python experiments/EXP-05_RobustnessChecks_2026-11-12/scripts/run_robustness.py
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

EXP = Path(__file__).resolve().parents[1]
ROOT = EXP.parents[1]
sys.path.insert(0, str(ROOT / "src"))
warnings.filterwarnings("ignore")
from eiv_core import fit_eiv, measurement_covariances  # noqa: E402
from utils import filter_abide_sample, load_abide_csv, resolve_input_path, sigma_x_from_age  # noqa: E402

DATA = Path(os.environ.get("DEV_DATA_ROOT", r"D:\Dev\Data")) / "eiv-abide-replication" / "fmriprep_confounds"
SEED = 20261112


def design(df, site=True, extra=None):
    a = df["AGE_AT_SCAN"].to_numpy(float)
    cols = [np.ones(len(df)), a, a**2, (df["SEX"] == 2).to_numpy(float), (df["DX_GROUP"] == 1).to_numpy(float)]
    if site:
        cols.append(pd.get_dummies(df["SITE_ID"], drop_first=True).to_numpy(float))
    if extra is not None:
        cols.append(extra)
    return np.column_stack(cols)


def resid(v, Z):
    return v - Z @ np.linalg.lstsq(Z, v, rcond=None)[0]


def fwl(df, y, Z=None):
    Z = design(df) if Z is None else Z
    x, yr = resid(df["FIQ"].to_numpy(float), Z), resid(np.asarray(y, float), Z)
    return x, yr, Z.shape[1]


def ols_rc(df, y, Z=None):
    x, yr, k = fwl(df, y, Z)
    n = len(x)
    sxx, sxy = x @ x / (n - k), x @ yr / (n - k)
    return sxy / sxx, sxy / (sxx - np.mean(df["sigma_x"].to_numpy(float) ** 2))


def pcr2(x, y, sx, sy):
    X, yy, C = measurement_covariances(x, y, sigma_X=sx, sigma_y=sy)
    f = fit_eiv(X, yy, C, method="structural", max_iter=2000)
    return float(f.coef[0]) if f.converged else np.nan


def wild_cluster_p(x, y, groups, reps=4999, seed=SEED):
    """Wild cluster bootstrap-t, Rademacher weights, null imposed (x, y already residualised)."""
    rng = np.random.default_rng(seed)
    g = pd.factorize(groups)[0]
    G = g.max() + 1

    def tstat(yy):
        b = x @ yy / (x @ x)
        u = yy - b * x
        s = np.bincount(g, weights=x * u, minlength=G)
        return b / np.sqrt(np.sum(s**2) / (x @ x) ** 2 * G / (G - 1))

    t0 = tstat(y)
    ts = np.array([tstat(y * rng.choice([-1.0, 1.0], G)[g]) for _ in range(reps)])  # null: b = 0
    return float(t0), float(np.mean(np.abs(ts) >= abs(t0)))


def main() -> None:
    t0 = time.time()
    out = EXP / "results"
    out.mkdir(parents=True, exist_ok=True)
    path = resolve_input_path(None)
    raw = load_abide_csv(path)
    base = filter_abide_sample(raw).reset_index(drop=True)
    base["sigma_x"] = sigma_x_from_age(base["AGE_AT_SCAN"])
    base["log_fd"] = np.log(base["func_mean_fd"])
    res: dict = {}

    # frame-level table
    man = json.loads((DATA / "manifest.json").read_text())
    files = {int(r["sub"][4:]): sorted(f["name"] for f in r["files"]) for r in man["records"]}
    rows = []
    for pos, sid in enumerate(base["SUB_ID"].astype(int)):
        if not files.get(sid):
            continue
        fd = pd.read_csv(DATA / files[sid][0], sep="\t", usecols=["framewise_displacement"]).iloc[:, 0].dropna().to_numpy()
        r = {"pos": pos, "fp_fd": fd.mean(), "n_frames": len(fd)}
        for nb in (4, 8, 16):
            r[f"noise{nb}"] = np.std([s.mean() for s in np.array_split(fd, nb)], ddof=1) / np.sqrt(nb)
        rows.append(r)
    fr = pd.DataFrame(rows).join(base, on="pos")

    # R1
    r1 = []
    Z = design(fr)
    x = resid(fr["FIQ"].to_numpy(float), Z)
    for nb in (4, 8, 16):
        for scale in ("raw", "log"):
            y0 = fr["fp_fd"].to_numpy(float)
            sy = fr[f"noise{nb}"].to_numpy(float)
            if scale == "log":
                sy, y0 = sy / y0, np.log(y0)
            r1.append({"blocks": nb, "scale": scale, "median_noise_sd": float(np.median(sy)),
                       "pcr2_slope_controls": pcr2(x, resid(y0, Z), fr["sigma_x"].to_numpy(float), sy)})
    res["R1_noise_blocks"] = r1

    # R2 / R3 selection and truncation (phenotypic FD)
    elig = raw[raw.FIQ.notna() & raw.func_mean_fd.notna() & raw.AGE_AT_SCAN.notna() & (raw.FIQ > 40)].copy()
    elig["sigma_x"] = sigma_x_from_age(elig["AGE_AT_SCAN"])
    sel = []
    for label, d in [("analysis sample (QC pass)", base), ("QC-fail included", elig),
                     ("mean FD <= 0.5 mm", base[base.func_mean_fd <= 0.5]), ("mean FD <= 0.2 mm", base[base.func_mean_fd <= 0.2])]:
        for scale in ("raw", "log"):
            y = d["func_mean_fd"].to_numpy(float)
            y = np.log(y) if scale == "log" else y
            o0, _ = ols_rc(d, y, np.ones((len(d), 1)))
            o3, r3 = ols_rc(d, y)
            sel.append({"sample": label, "n": len(d), "scale": scale, "ols_no_controls": o0, "ols_controls": o3, "rc_controls": r3})
    pd.DataFrame(sel).to_csv(out / "r2_r3_selection_truncation.csv", index=False)

    # R4 linearity (log, controls)
    f4 = "log_fd ~ FIQ + I((FIQ-100)**2) + AGE_AT_SCAN + I(AGE_AT_SCAN**2) + C(SEX) + C(DX_GROUP) + C(SITE_ID)"
    m4 = smf.ols(f4, base).fit(cov_type="cluster", cov_kwds={"groups": base.SITE_ID})
    terc = pd.qcut(base.FIQ, 3, labels=["low", "mid", "high"])
    res["R4_linearity"] = {"quadratic_coef": float(m4.params["I((FIQ - 100) ** 2)"]),
                           "quadratic_cluster_p": float(m4.pvalues["I((FIQ - 100) ** 2)"]),
                           "slope_by_iq_tercile_log_controls": {str(t): float(ols_rc(base[terc == t], base.loc[terc == t, "log_fd"])[0])
                                                               for t in ["low", "mid", "high"]},
                           "tercile_iq_ranges": {str(t): [float(base.FIQ[terc == t].min()), float(base.FIQ[terc == t].max())] for t in ["low", "mid", "high"]}}

    # R5 inference with few clusters
    r5 = {}
    for scale, ycol in (("raw", "func_mean_fd"), ("log", "log_fd")):
        xr, yr, k = fwl(base, base[ycol])
        t, p_wild = wild_cluster_p(xr, yr, base.SITE_ID.to_numpy())
        cr = sm.OLS(yr, xr).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(base.SITE_ID)[0]})
        mm = smf.mixedlm(f"{ycol} ~ FIQ + AGE_AT_SCAN + I(AGE_AT_SCAN**2) + C(SEX) + C(DX_GROUP)", base,
                         groups=base.SITE_ID, re_formula="~FIQ").fit(reml=True)
        r5[scale] = {"ols_controls": float(cr.params[0]), "cluster_robust_se": float(cr.bse[0]), "cluster_robust_p": float(cr.pvalues[0]),
                     "wild_cluster_t": t, "wild_cluster_p": p_wild,
                     "mixed_random_slope_fixed_effect": float(mm.params["FIQ"]), "mixed_p": float(mm.pvalues["FIQ"]),
                     "mixed_converged": bool(mm.converged)}
    res["R5_inference"] = r5

    # R6 leave-one-site-out
    loso = []
    for s in sorted(base.SITE_ID.unique()):
        d = base[base.SITE_ID != s]
        loso.append({"dropped_site": s, "raw": ols_rc(d, d.func_mean_fd)[0], "log": ols_rc(d, d.log_fd)[0]})
    loso = pd.DataFrame(loso)
    loso.to_csv(out / "r6_leave_one_site_out.csv", index=False)
    res["R6_leave_one_site_out"] = {c: {"min": float(loso[c].min()), "max": float(loso[c].max()),
                                        "site_at_min": loso.loc[loso[c].idxmin(), "dropped_site"],
                                        "site_at_max": loso.loc[loso[c].idxmax(), "dropped_site"],
                                        "all_negative": bool((loso[c] < 0).all())} for c in ("raw", "log")}

    # R7 simulation: proportional noise. Truth is linear on the log scale.
    rng = np.random.default_rng(SEED)
    n = len(fr)
    xt_pool = fr["FIQ"].to_numpy(float)
    b_log, sd_log, cv = -0.006, 0.60, 0.12  # slope, person SD of log FD, noise CV of the run mean (EXP-04 scale)
    sims = []
    for _ in range(300):
        xt = rng.choice(xt_pool, n)
        sx = rng.choice(fr["sigma_x"].to_numpy(float), n)
        yt = np.exp(-1.2 + b_log * (xt - 100) + rng.normal(0, sd_log, n))
        sy = cv * yt * np.exp(rng.normal(0, 0.5, n))            # per-subject noise SD, proportional to level
        xo, yo = xt + rng.normal(0, sx), yt + rng.normal(0, sy)
        yo = np.maximum(yo, 1e-3)
        target_raw = np.cov(xt, yt)[0, 1] / np.var(xt, ddof=1)    # linear projection of true y on true x
        xc, yc = xo - xo.mean(), yo - yo.mean()
        lyc = np.log(yo) - np.log(yo).mean()
        sims.append({"target_raw": target_raw, "ols_raw": xc @ yc / (xc @ xc), "pcr2_raw": pcr2(xc, yc, sx, sy),
                     "ols_log": xc @ lyc / (xc @ xc), "pcr2_log": pcr2(xc, lyc, sx, sy / yo)})
    sims = pd.DataFrame(sims)
    sims.to_csv(out / "r7_simulation.csv", index=False)
    res["R7_simulation"] = {"design": {"true_log_slope": b_log, "person_sd_log": sd_log, "noise_cv": cv, "n": n, "sims": len(sims)},
                            "means": sims.mean().to_dict(),
                            "raw_ratio_pcr2_over_target": float(sims.pcr2_raw.mean() / sims.target_raw.mean()),
                            "raw_ratio_ols_over_target": float(sims.ols_raw.mean() / sims.target_raw.mean()),
                            "log_ratio_pcr2_over_truth": float(sims.pcr2_log.mean() / b_log),
                            "log_ratio_ols_over_truth": float(sims.ols_log.mean() / b_log)}

    # R8 scan length and test type
    fl = fr.groupby("SITE_ID").agg(n_frames=("n_frames", "median"), mean_fd=("fp_fd", "mean"))
    tt = pd.get_dummies(base["FIQ_TEST_TYPE"].fillna("missing"), drop_first=True).to_numpy(float)
    res["R8"] = {"site_median_frames_min_max": [float(fl.n_frames.min()), float(fl.n_frames.max())],
                 "spearman_site_frames_vs_site_mean_fd": float(fl.corr(method="spearman").iloc[0, 1]),
                 "log_slope_controls": float(ols_rc(base, base.log_fd)[0]),
                 "log_slope_controls_plus_iq_test_type": float(ols_rc(base, base.log_fd, design(base, extra=tt))[0]),
                 "raw_slope_controls_plus_iq_test_type": float(ols_rc(base, base.func_mean_fd, design(base, extra=tt))[0])}

    for p in out.glob("*.csv"):
        assert not set(pd.read_csv(p, nrows=0).columns) & {"SUB_ID", "FILE_ID", "subject", "pos"}
    res["receipt"] = {"phenotypic_csv": hashlib.sha256(Path(path).read_bytes()).hexdigest(),
                      "script": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), "seed": SEED,
                      "runtime_s": round(time.time() - t0, 1)}
    (out / "summary.json").write_text(json.dumps(res, indent=2, default=str))
    print(json.dumps(res, indent=1, default=str))
    with pd.option_context("display.width", 200, "display.float_format", "{:.5g}".format):
        print(pd.DataFrame(sel).to_string(index=False))


if __name__ == "__main__":
    main()

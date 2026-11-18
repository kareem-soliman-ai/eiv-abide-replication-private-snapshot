"""EXP-02: full EIV method panel on ABIDE-I IQ -> mean FD.

Estimators (all on identical inputs):
  ols             naive OLS
  rc              regression calibration (reliability ratio; x-error only)
  simex           SIMEX, quadratic extrapolation (x-error only)
  deming          Deming, constant error ratio delta = mean(sy^2)/mean(sx^2)
  wls_sy          WLS 1/sy^2 (weighting only; not EIV)
  york            corrected York (functional; heteroscedastic; no scatter)
  coord_descent   objective-coherent PCR repair (functional; no scatter)
  legacy_pcr      original PCR rule (the paper)
  pcr_functional  updated PCR, Rasad eiv_core functional (no scatter)
  hazy_york       York + intrinsic scatter tau2 set by reduced chi2 = 1
  structural_ml   Gaussian structural EIV with scatter (PCR-methodology reference)
  pcr2_structural updated PCR-2, Rasad eiv_core structural (with scatter)
  mixed_rs        random-intercept/random-slope mixed model (non-EIV reference)

Specifications:
  A  paper: sigma_x by age band, sigma_y = within-site SD of FD
  B  constant sigma_y = pooled within-site SD (no site-specific noise)
  C  log(FD), sigma_y = within-site SD of log FD

Run from repo root:
  python experiments/EXP-02_EIVMethodPanel_2026-11-12/scripts/run_panel.py --bootstrap 500 --n-jobs 8
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from joblib import Parallel, delayed

EXP = Path(__file__).resolve().parents[1]
ROOT = EXP.parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(EXP / "vendor"))
warnings.filterwarnings("ignore")

import eiv_methods as em  # noqa: E402
import pcr2_reference as ref  # noqa: E402
from eiv_core import METHOD_VERSION, fit_eiv, measurement_covariances  # noqa: E402
from pcr import fit_pcr_regression  # noqa: E402
from utils import attach_uncertainty_columns, filter_abide_sample, load_abide_csv, manual_r2, resolve_input_path  # noqa: E402

METHODS = ["ols", "rc", "simex", "deming", "wls_sy", "york", "coord_descent", "legacy_pcr",
           "pcr_functional", "hazy_york", "structural_ml", "pcr2_structural", "mixed_rs"]
FAMILY = {"ols": "naive", "rc": "x-error only", "simex": "x-error only", "deming": "functional, homoscedastic",
          "wls_sy": "weighting only", "york": "functional, no scatter", "coord_descent": "functional, no scatter",
          "legacy_pcr": "functional, no scatter (flawed)", "pcr_functional": "functional, no scatter",
          "hazy_york": "functional + scatter", "structural_ml": "structural + scatter",
          "pcr2_structural": "structural + scatter", "mixed_rs": "multilevel (non-EIV)"}
SPECS = ["A_paper_site_sd", "B_constant_sigma_y", "C_log_fd_site_sd"]


def prepare(frame: pd.DataFrame, spec: str) -> pd.DataFrame:
    """Attach sigma columns for a spec. Called on whatever rows are being fitted (full, fold, bootstrap)."""
    g = frame.copy()
    if spec.startswith("C"):
        g["func_mean_fd"] = np.log(g["raw_fd"])
    else:
        g["func_mean_fd"] = g["raw_fd"]
    g = attach_uncertainty_columns(g)
    if spec.startswith("B"):
        g["sigma_y"] = float(np.sqrt(np.mean(g["sigma_y"] ** 2)))
    return g


def fit_one(m: str, g: pd.DataFrame, simex_sims: int = 200) -> dict:
    x, y = g["FIQ"].to_numpy(float), g["func_mean_fd"].to_numpy(float)
    sx, sy = g["sigma_x"].to_numpy(float), g["sigma_y"].to_numpy(float)
    out: dict = {"converged": True}
    try:
        if m == "ols":
            b, a = em.fit_ols(x, y)
        elif m == "rc":
            b, a = em.fit_regression_calibration(x, y, sx)
        elif m == "simex":
            r = em.fit_simex(x, y, sx, n_sim=simex_sims)
            b, a = r.slope, r.intercept
        elif m == "deming":
            b, a = em.fit_deming(x, y, float(np.mean(sy**2) / np.mean(sx**2)))
        elif m == "wls_sy":
            b, a = em.fit_wls(x, y, sy)
        elif m == "york":
            r = em.fit_york(x, y, sx, sy)
            b, a, out["converged"] = r.slope, r.intercept, r.converged
        elif m == "coord_descent":
            b, a, _, out["converged"] = ref.coordinate_descent(x, y, sx, sy)
        elif m == "legacy_pcr":
            r = fit_pcr_regression(x, y, sx, sy)
            b, a, out["converged"] = r.slope, r.intercept, r.converged
        elif m in ("pcr_functional", "pcr2_structural"):
            X, yy, C = measurement_covariances(x, y, sigma_X=sx, sigma_y=sy)
            f = fit_eiv(X, yy, C, method="functional" if m == "pcr_functional" else "structural", max_iter=2000)
            b, a, out["converged"] = float(f.coef[0]), f.intercept, f.converged
            if m == "pcr2_structural":
                out["tau2"] = f.intrinsic_variance
        elif m == "hazy_york":
            b, a, out["tau2"] = ref.hazy_york(x, y, sx, sy)
        elif m == "structural_ml":
            b, a, out["tau2"], _, _ = ref.structural_ml(x, y, sx, sy)
            out["converged"] = abs(b) < 49.9  # reference fit reports no status; flag bound hits
        elif m == "mixed_rs":
            r = smf.mixedlm("func_mean_fd ~ FIQ", g, groups=g["SITE_ID"], re_formula="~FIQ").fit(reml=True)
            b, a, out["converged"] = float(r.params["FIQ"]), float(r.params["Intercept"]), bool(r.converged)
        out.update(slope=float(b), intercept=float(a))
        out["converged"] = bool(out["converged"] and np.isfinite(b))
    except Exception as exc:  # degrade, record, continue
        out.update(slope=np.nan, intercept=np.nan, converged=False, error=f"{type(exc).__name__}: {exc}")
    return out


def loso(frame: pd.DataFrame, spec: str, methods: list[str]) -> dict[str, float]:
    truth, preds = [], {m: [] for m in methods}
    for site in sorted(frame["SITE_ID"].unique()):
        train = prepare(frame[frame.SITE_ID != site], spec)
        test = prepare(frame[frame.SITE_ID == site], spec)
        truth.extend(test["func_mean_fd"].tolist())
        for m in methods:
            f = fit_one(m, train, simex_sims=50)
            preds[m].extend((f["slope"] * test["FIQ"] + f["intercept"]).tolist())
    t = np.asarray(truth)
    return {m: (manual_r2(t, np.asarray(p)) if np.all(np.isfinite(p)) else np.nan) for m, p in preds.items()}


def boot_rep(frame: pd.DataFrame, spec: str, sites: np.ndarray, seed: int, methods: list[str]) -> dict:
    rng = np.random.default_rng(seed)
    parts = []
    for k, s in enumerate(rng.choice(sites, size=len(sites), replace=True)):
        g = frame[frame.SITE_ID == s].copy()
        g["SITE_ID"] = f"{s}#{k}"
        parts.append(g)
    g = prepare(pd.concat(parts, ignore_index=True), spec)
    row = {"spec": spec, "seed": seed}
    for m in methods:
        f = fit_one(m, g, simex_sims=50)
        row[f"{m}_slope"], row[f"{m}_converged"] = f["slope"], f["converged"]
    return row


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input")
    ap.add_argument("--bootstrap", type=int, default=500)
    ap.add_argument("--seed", type=int, default=20261112)
    ap.add_argument("--n-jobs", type=int, default=1, help="parallel bootstrap workers (8 recommended here)")
    args = ap.parse_args()
    t0 = time.time()
    out = EXP / "results"
    out.mkdir(exist_ok=True)

    path = resolve_input_path(args.input)
    base = filter_abide_sample(load_abide_csv(path)).rename(columns={"func_mean_fd": "raw_fd"})
    base = base[["SITE_ID", "FIQ", "raw_fd", "AGE_AT_SCAN"]].reset_index(drop=True)
    assert len(base) == 935 and base.SITE_ID.nunique() == 19

    rows = []
    for spec in SPECS:
        g = prepare(base, spec)
        lr = loso(base, spec, METHODS)
        ols_b = fit_one("ols", g)["slope"]
        for m in METHODS:
            f = fit_one(m, g)
            rows.append({"spec": spec, "method": m, "family": FAMILY[m], "slope": f["slope"],
                         "intercept": f["intercept"], "converged": f["converged"],
                         "ols_over_method": abs(ols_b / f["slope"]) if f["slope"] else np.nan,
                         "tau2": f.get("tau2"), "loso_r2": lr[m], "error": f.get("error")})
        print(f"[{time.time()-t0:6.1f}s] point fits + LOSO done for {spec}", flush=True)
    point = pd.DataFrame(rows)
    point.to_csv(out / "panel_point_estimates.csv", index=False)

    boot_methods = [m for m in METHODS if m != "mixed_rs"]  # mixed model too slow/fragile per rep
    sites = np.array(sorted(base.SITE_ID.unique()))
    seeds = np.random.SeedSequence(args.seed).generate_state(args.bootstrap)
    reps = Parallel(n_jobs=args.n_jobs, backend="loky")(
        delayed(boot_rep)(base, spec, sites, int(s), boot_methods) for spec in SPECS for s in seeds)
    boot = pd.DataFrame(reps)
    boot.to_csv(out / "panel_bootstrap_replicates.csv", index=False)

    ci = []
    for spec, bs in boot.groupby("spec"):
        for m in boot_methods:
            ok = bs[f"{m}_converged"].astype(bool) & np.isfinite(bs[f"{m}_slope"])
            s = bs.loc[ok, f"{m}_slope"]
            ci.append({"spec": spec, "method": m, "n_success": int(ok.sum()), "n_attempted": len(bs),
                       "ci_low": s.quantile(.025), "ci_high": s.quantile(.975),
                       "share_negative": float((s < 0).mean()), "excludes_zero": bool(s.quantile(.975) < 0 or s.quantile(.025) > 0)})
    ci = pd.DataFrame(ci)
    ci.to_csv(out / "panel_bootstrap_intervals.csv", index=False)

    for p in out.glob("*.csv"):
        assert not set(pd.read_csv(p, nrows=0).columns) & {"SUB_ID", "FILE_ID"}

    h = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()  # noqa: E731
    receipt = {"experiment": EXP.name, "eiv_core_version": METHOD_VERSION,
               "sources": {"input_csv": h(path), "eiv_core.py": h(ROOT / "src/eiv_core.py"),
                           "pcr.py": h(ROOT / "src/pcr.py"),
                           "vendor/eiv_methods.py (PCR-methodology working tree, uncommitted at e858adc)": h(EXP / "vendor/eiv_methods.py"),
                           "vendor/pcr2_reference.py (PCR-methodology working tree, uncommitted at e858adc)": h(EXP / "vendor/pcr2_reference.py"),
                           "script": h(__file__)},
               "bootstrap": {"reps_per_spec": args.bootstrap, "seed": args.seed, "unit": "site", "simex_sims": 50},
               "python": platform.python_version(), "runtime_s": round(time.time() - t0, 1)}
    (out / "receipt.json").write_text(json.dumps(receipt, indent=2))
    with pd.option_context("display.width", 220, "display.max_rows", 100):
        print(point.drop(columns=["intercept", "error"]).to_string(index=False))
        print(ci.to_string(index=False))
    print(f"runtime {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()

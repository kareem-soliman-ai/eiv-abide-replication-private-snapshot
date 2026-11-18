"""EXP-01: rerun the ABIDE-I IQ-motion analysis with the corrected EIV core.

Compares, on identical inputs (n=935, 19 sites, the paper's sigma_x/sigma_y):
  ols          ordinary least squares (the paper's comparator)
  legacy_pcr   the original PCR alternating rule used in the paper
  functional   corrected functional Gaussian EIV (Rasad eiv_core, rasad-eiv/2)
  structural   structural PCR-2: latent-predictor distribution + intrinsic scatter
Diagnostics that isolate the mechanism:
  wls_sy       WLS with weights 1/sigma_y^2 (no predictor-error correction at all)
  rc           regression calibration (predictor-error correction only, no weighting)
  site_fe      OLS with site fixed effects (pooled within-site slope)

Outputs are privacy-safe (no SUB_ID / FILE_ID). Run from the repo root:
  python experiments/EXP-01_CorrectedEIVRerun_2026-11-12/scripts/rerun_corrected.py
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.stats import linregress

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))

from eiv_core import METHOD_VERSION, fit_eiv, measurement_covariances  # noqa: E402
from pcr import fit_pcr_regression  # noqa: E402
from sensitivity_grid import MULTIPLIERS  # noqa: E402
from utils import (  # noqa: E402
    assign_motion_tiers,
    attach_uncertainty_columns,
    filter_abide_sample,
    load_abide_csv,
    manual_r2,
    resolve_input_path,
)

EXP_DIR = Path(__file__).resolve().parents[1]
METHODS = ["ols", "legacy_pcr", "functional", "structural", "wls_sy", "rc", "site_fe"]
EIV_METHODS = ["legacy_pcr", "functional", "structural"]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fit_all(frame: pd.DataFrame, methods: list[str] = METHODS) -> dict[str, dict]:
    """Fit each method; return slope, intercept, convergence and extras."""
    x = frame["FIQ"].to_numpy(float)
    y = frame["func_mean_fd"].to_numpy(float)
    sx = frame["sigma_x"].to_numpy(float)
    sy = frame["sigma_y"].to_numpy(float)
    out: dict[str, dict] = {}
    for m in methods:
        if m == "ols":
            r = linregress(x, y)
            out[m] = dict(slope=float(r.slope), intercept=float(r.intercept), converged=True)
        elif m == "legacy_pcr":
            f = fit_pcr_regression(x, y, sx, sy)
            out[m] = dict(slope=f.slope, intercept=f.intercept, converged=bool(f.converged),
                          iterations=f.iterations)
        elif m in ("functional", "structural"):
            X, yy, C = measurement_covariances(x, y, sigma_X=sx, sigma_y=sy)
            f = fit_eiv(X, yy, C, method=m, max_iter=2000)
            d = dict(slope=float(f.coef[0]), intercept=float(f.intercept),
                     converged=bool(f.converged),
                     projected_gradient=float(f.diagnostics["projected_gradient"]))
            if m == "functional":
                r = yy - f.intercept - f.coef[0] * x
                d["reduced_chi2"] = float(np.sum(r**2 / (sy**2 + f.coef[0]**2 * sx**2)) / (len(x) - 2))
            else:
                d["intrinsic_variance"] = f.intrinsic_variance
                d["intrinsic_boundary"] = bool(f.diagnostics["intrinsic_variance_boundary"])
                d["latent_x_variance"] = float(f.latent_covariance[0, 0])
            out[m] = d
        elif m == "wls_sy":
            w = 1.0 / sy**2
            xb, yb = np.average(x, weights=w), np.average(y, weights=w)
            b = np.sum(w * (x - xb) * (y - yb)) / np.sum(w * (x - xb) ** 2)
            out[m] = dict(slope=float(b), intercept=float(yb - b * xb), converged=True)
        elif m == "rc":
            reliability = 1 - np.mean(sx**2) / np.var(x)
            b = np.cov(x, y, bias=True)[0, 1] / (np.var(x) - np.mean(sx**2))
            out[m] = dict(slope=float(b), intercept=float(y.mean() - b * x.mean()),
                          converged=True, reliability_x=float(reliability))
        elif m == "site_fe":
            g = frame["SITE_ID"].to_numpy()
            xd = x - pd.Series(x).groupby(g).transform("mean").to_numpy()
            yd = y - pd.Series(y).groupby(g).transform("mean").to_numpy()
            b = float(np.sum(xd * yd) / np.sum(xd**2))
            # Intercept is site-specific; report the grand-mean line for plotting only.
            out[m] = dict(slope=b, intercept=float(y.mean() - b * x.mean()), converged=True)
    return out


def loso(frame: pd.DataFrame, methods: list[str], sx_mult=1.0, sy_mult=1.0) -> dict[str, float]:
    """Leave-site-out R^2 per method, re-deriving sigma_y on the training sites (as the paper did)."""
    preds = {m: [] for m in methods}
    truth: list[float] = []
    for site in sorted(frame["SITE_ID"].unique()):
        train = attach_uncertainty_columns(frame.loc[frame["SITE_ID"] != site],
                                           sigma_x_multiplier=sx_mult, sigma_y_multiplier=sy_mult)
        test = frame.loc[frame["SITE_ID"] == site]
        fits = fit_all(train, methods)
        xt = test["FIQ"].to_numpy(float)
        truth.extend(test["func_mean_fd"].tolist())
        for m in methods:
            preds[m].extend((fits[m]["slope"] * xt + fits[m]["intercept"]).tolist())
    t = np.asarray(truth)
    return {m: manual_r2(t, np.asarray(p)) for m, p in preds.items()}


def site_bootstrap(frame: pd.DataFrame, reps: int, seed: int) -> pd.DataFrame:
    """Resample whole sites with replacement; sigma_y recomputed per resampled site copy."""
    rng = np.random.default_rng(seed)
    sites = np.array(sorted(frame["SITE_ID"].unique()))
    groups = {s: frame.loc[frame["SITE_ID"] == s] for s in sites}
    rows = []
    for rep in range(reps):
        draw = rng.choice(sites, size=len(sites), replace=True)
        parts = []
        for k, s in enumerate(draw):
            g = groups[s].copy()
            g["SITE_ID"] = f"{s}#{k}"  # duplicate draws are distinct clusters
            parts.append(g)
        boot = attach_uncertainty_columns(pd.concat(parts, ignore_index=True))
        fits = fit_all(boot, ["ols", "legacy_pcr", "functional", "structural", "wls_sy", "site_fe"])
        row = {"rep": rep}
        for m, f in fits.items():
            row[f"{m}_slope"] = f["slope"]
            row[f"{m}_converged"] = f["converged"]
        row["structural_tau2"] = fits["structural"]["intrinsic_variance"]
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", help="ABIDE-I phenotypic CSV (default data/abide_phenotypic.csv)")
    ap.add_argument("--bootstrap", type=int, default=1000, help="site-bootstrap replicates")
    ap.add_argument("--seed", type=int, default=20261112)
    ap.add_argument("--skip-grid", action="store_true")
    args = ap.parse_args()
    t0 = time.time()
    out = EXP_DIR / "results"
    out.mkdir(parents=True, exist_ok=True)

    input_path = resolve_input_path(args.input)
    frame = assign_motion_tiers(attach_uncertainty_columns(filter_abide_sample(load_abide_csv(input_path))))
    assert len(frame) == 935 and frame["SITE_ID"].nunique() == 19, "sample filter drifted"

    # 1. Headline fits
    fits = fit_all(frame)
    loso_r2 = loso(frame, METHODS)
    ols_b = fits["ols"]["slope"]
    rows = []
    for m in METHODS:
        f = fits[m]
        rows.append({"method": m, "slope": f["slope"], "intercept": f["intercept"],
                     "converged": f["converged"], "ols_over_method": abs(ols_b / f["slope"]),
                     "loso_r2": loso_r2[m],
                     **{k: v for k, v in f.items() if k not in ("slope", "intercept", "converged")}})
    headline = pd.DataFrame(rows)
    headline.to_csv(out / "headline_by_method.csv", index=False)

    # 2. Per-tier fits (does the correction change the within-tier picture?)
    tier_rows = []
    for tier, tf in frame.groupby("motion_tier", observed=True):
        tf = attach_uncertainty_columns(tf)
        tfits = fit_all(tf, ["ols", "legacy_pcr", "functional", "structural"])
        for m, f in tfits.items():
            tier_rows.append({"tier": str(tier), "n": len(tf), "method": m, "slope": f["slope"],
                              "converged": f["converged"],
                              "intrinsic_variance": f.get("intrinsic_variance")})
    pd.DataFrame(tier_rows).to_csv(out / "tier_slopes_by_method.csv", index=False)

    # 3. Sensitivity grid (8x8 sigma multipliers) for the EIV methods
    if not args.skip_grid:
        grid_rows = []
        for mx in MULTIPLIERS:
            for my in MULTIPLIERS:
                g = attach_uncertainty_columns(frame, sigma_x_multiplier=mx, sigma_y_multiplier=my)
                gf = fit_all(g, EIV_METHODS)
                lr = loso(frame, EIV_METHODS, mx, my)
                for m in EIV_METHODS:
                    grid_rows.append({"sigma_x_multiplier": mx, "sigma_y_multiplier": my, "method": m,
                                      "slope": gf[m]["slope"], "ols_over_method": abs(ols_b / gf[m]["slope"]),
                                      "converged": gf[m]["converged"], "loso_r2": lr[m],
                                      "intrinsic_variance": gf[m].get("intrinsic_variance")})
        pd.DataFrame(grid_rows).to_csv(out / "sensitivity_grid_by_method.csv", index=False)

    # 3b. sigma_y STRUCTURE sensitivity. Constant multipliers leave the relative
    # site weights unchanged, so the 8x8 grid cannot test the key assumption:
    # that between-site differences in FD spread are measurement error.
    site_sd = frame["sigma_y"]
    pooled_within = float(np.sqrt(np.mean(site_sd**2)))
    structure_rows = []
    for label, sy_vec in [
        ("site_sd (paper)", site_sd.to_numpy(float)),
        ("constant = pooled within-site SD", np.full(len(frame), pooled_within)),
        ("constant = 0.25 x pooled", np.full(len(frame), 0.25 * pooled_within)),
        ("constant = 0.05 x pooled", np.full(len(frame), 0.05 * pooled_within)),
    ]:
        g = frame.copy()
        g["sigma_y"] = sy_vec
        gf = fit_all(g, EIV_METHODS)
        for m in EIV_METHODS:
            structure_rows.append({"sigma_y_structure": label, "sigma_y_mean": float(sy_vec.mean()),
                                   "method": m, "slope": gf[m]["slope"],
                                   "ols_over_method": abs(ols_b / gf[m]["slope"]),
                                   "converged": gf[m]["converged"],
                                   "intrinsic_variance": gf[m].get("intrinsic_variance")})
    structure = pd.DataFrame(structure_rows)
    structure.to_csv(out / "sigma_y_structure_sensitivity.csv", index=False)

    # 4. Site-cluster bootstrap
    boot = site_bootstrap(frame, args.bootstrap, args.seed)
    boot.to_csv(out / "site_bootstrap_replicates.csv", index=False)
    ci_rows = []
    for m in ["ols", "legacy_pcr", "functional", "structural", "wls_sy", "site_fe"]:
        ok = boot[f"{m}_converged"].astype(bool)
        s = boot.loc[ok, f"{m}_slope"]
        diff = (boot.loc[ok, "ols_slope"] - s)
        ci_rows.append({"method": m, "point": fits[m]["slope"], "n_success": int(ok.sum()),
                        "n_attempted": len(boot), "ci_low": s.quantile(.025), "ci_high": s.quantile(.975),
                        "ols_minus_method_ci_low": diff.quantile(.025),
                        "ols_minus_method_ci_high": diff.quantile(.975),
                        "share_slope_negative": float((s < 0).mean())})
    pd.DataFrame(ci_rows).to_csv(out / "site_bootstrap_intervals.csv", index=False)

    # Privacy guard: no subject identifiers in any output
    for p in out.glob("*.csv"):
        cols = set(pd.read_csv(p, nrows=0).columns)
        assert not cols & {"SUB_ID", "FILE_ID", "subject"}, f"identifier column in {p.name}"

    receipt = {
        "experiment": "EXP-01_CorrectedEIVRerun_2026-11-12",
        "eiv_core_method_version": METHOD_VERSION,
        "eiv_core_source": "Rasad src/eiv_core.py @ 302b394 (byte-identical copy)",
        "hashes": {"input_csv": sha256(input_path), "eiv_core.py": sha256(ROOT / "src" / "eiv_core.py"),
                   "pcr.py": sha256(ROOT / "src" / "pcr.py"), "script": sha256(Path(__file__))},
        "n_subjects": int(len(frame)), "n_sites": int(frame["SITE_ID"].nunique()),
        "bootstrap": {"reps": args.bootstrap, "seed": args.seed, "unit": "site"},
        "python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
        "pandas": pd.__version__, "runtime_s": round(time.time() - t0, 1),
        "headline": headline.to_dict(orient="records"),
    }
    (out / "receipt.json").write_text(json.dumps(receipt, indent=2, default=float))
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(headline.to_string(index=False))
        print(pd.DataFrame(ci_rows).to_string(index=False))
        print(structure.to_string(index=False))
    print(f"runtime {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()

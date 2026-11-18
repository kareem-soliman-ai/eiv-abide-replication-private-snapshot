"""EXP-01 stress tests: can any honest reading rescue the paper's 4.7x claim?

S1 Generative check: under the paper's own model (common slope, sigma_y = pure
   outcome noise, sigma_x = IQ noise) is OLS 'inflated'? Simulated from ABIDE's
   own x, sites and sigmas.
S2 Slope heterogeneity: do within-site slopes vary with site motion level?
   (the only mechanism that can separate OLS from 1/sigma_y^2 WLS)
S3 Scale: does the gap survive log(FD)?
S4 Robust and multilevel fits: Huber, median regression, random-slope mixed model.
S5 Site jackknife: is the ratio driven by one site?
S6 Is sigma_y plausibly instrument noise? Site SD vs site mean (skew/heteroscedasticity).
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import linregress, spearmanr

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from pcr import fit_pcr_regression  # noqa: E402
from utils import attach_uncertainty_columns, filter_abide_sample, load_abide_csv, resolve_input_path  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "results"
warnings.filterwarnings("ignore")


def wls(x, y, w):
    xb, yb = np.average(x, weights=w), np.average(y, weights=w)
    return float(np.sum(w * (x - xb) * (y - yb)) / np.sum(w * (x - xb) ** 2))


def pcr_slope(df):
    return fit_pcr_regression(df.FIQ.values, df.func_mean_fd.values, df.sigma_x.values, df.sigma_y.values).slope


def main(seed: int = 20261112, sims: int = 500) -> None:
    rng = np.random.default_rng(seed)
    df = attach_uncertainty_columns(filter_abide_sample(load_abide_csv(resolve_input_path(None))))
    x, y, sx, sy = (df[c].to_numpy(float) for c in ("FIQ", "func_mean_fd", "sigma_x", "sigma_y"))
    ols = linregress(x, y).slope
    res: dict = {"ols_slope": float(ols), "pcr_slope": float(pcr_slope(df))}

    # S1: paper's generative model. True slope set to either the OLS or the PCR value.
    s1 = []
    for true_b in (ols, res["pcr_slope"]):
        for _ in range(sims):
            xt = x + rng.normal(0, sx)  # treat observed x as latent + noise draw
            xo = xt + rng.normal(0, sx)
            yo = 0.2 + true_b * xt + rng.normal(0, sy)
            g = df.copy(); g["FIQ"], g["func_mean_fd"] = xo, yo
            s1.append({"true_slope": true_b, "ols": linregress(xo, yo).slope,
                       "wls_sy": wls(xo, yo, 1 / sy**2), "pcr": pcr_slope(g)})
    s1 = pd.DataFrame(s1)
    s1_sum = s1.groupby("true_slope").agg(["mean", "std"])
    s1_sum["ratio_ols_over_pcr_mean"] = s1_sum[("ols", "mean")] / s1_sum[("pcr", "mean")]
    s1_sum.to_csv(OUT / "stress_s1_generative.csv")
    res["S1"] = {f"true={b:.6f}": {"ols_mean": float(r[("ols", "mean")]), "pcr_mean": float(r[("pcr", "mean")]),
                                   "wls_mean": float(r[("wls_sy", "mean")]),
                                   "ols_over_pcr": float(r[("ols", "mean")] / r[("pcr", "mean")])}
                 for b, r in s1_sum.iterrows()}

    # S2: site slopes vs site motion level
    sites = []
    for s, g in df.groupby("SITE_ID"):
        lr = linregress(g.FIQ, g.func_mean_fd)
        sites.append({"n": len(g), "mean_fd": g.func_mean_fd.mean(), "sd_fd": g.func_mean_fd.std(),
                      "slope": lr.slope, "se": lr.stderr})
    sites = pd.DataFrame(sites)
    rho, p = spearmanr(sites.mean_fd, sites.slope)
    meta = sm.WLS(sites.slope, sm.add_constant(sites.mean_fd), weights=1 / sites.se**2).fit()
    inter = smf.ols("func_mean_fd ~ FIQ * site_mean_fd + C(SITE_ID)", data=df).fit(cov_type="cluster",
                                                                                     cov_kwds={"groups": df.SITE_ID})
    res["S2"] = {"spearman_site_slope_vs_mean_fd": float(rho), "p": float(p),
                 "meta_regression_coef": float(meta.params["mean_fd"]), "meta_p": float(meta.pvalues["mean_fd"]),
                 "FIQ_x_site_mean_fd_coef": float(inter.params["FIQ:site_mean_fd"]),
                 "FIQ_x_site_mean_fd_cluster_p": float(inter.pvalues["FIQ:site_mean_fd"])}

    # S3: log scale. Site-SD of log FD as sigma_y.
    lg = df.copy(); lg["func_mean_fd"] = np.log(df.func_mean_fd)
    lg = attach_uncertainty_columns(lg)
    l_ols = linregress(lg.FIQ, lg.func_mean_fd).slope
    l_wls = wls(lg.FIQ.values, lg.func_mean_fd.values, 1 / lg.sigma_y.values**2)
    res["S3"] = {"log_ols": float(l_ols), "log_wls_sy": float(l_wls), "log_pcr": float(pcr_slope(lg)),
                 "log_ols_over_pcr": float(l_ols / pcr_slope(lg)),
                 "site_sd_logfd_cv": float(lg.groupby("SITE_ID").sigma_y.first().std() / lg.groupby("SITE_ID").sigma_y.first().mean()),
                 "site_sd_fd_cv": float(df.groupby("SITE_ID").sigma_y.first().std() / df.groupby("SITE_ID").sigma_y.first().mean())}

    # S4: robust + multilevel
    X = sm.add_constant(df.FIQ)
    huber = sm.RLM(df.func_mean_fd, X, M=sm.robust.norms.HuberT()).fit().params["FIQ"]
    median = smf.quantreg("func_mean_fd ~ FIQ", df).fit(q=0.5).params["FIQ"]
    mixed = smf.mixedlm("func_mean_fd ~ FIQ", df, groups=df.SITE_ID, re_formula="~FIQ").fit(reml=True)
    mixed_log = smf.mixedlm("func_mean_fd ~ FIQ", lg, groups=lg.SITE_ID, re_formula="~FIQ").fit(reml=True)
    res["S4"] = {"huber": float(huber), "median": float(median),
                 "mixed_random_slope_fixed_FIQ": float(mixed.params["FIQ"]),
                 "mixed_FIQ_p": float(mixed.pvalues["FIQ"]),
                 "mixed_random_slope_var": float(mixed.cov_re.iloc[1, 1]),
                 "mixed_converged": bool(mixed.converged),
                 "mixed_log_FIQ": float(mixed_log.params["FIQ"]), "mixed_log_p": float(mixed_log.pvalues["FIQ"])}

    # S5: site jackknife of the OLS/PCR ratio
    jk = []
    for s in sorted(df.SITE_ID.unique()):
        g = attach_uncertainty_columns(df[df.SITE_ID != s])
        o, pc = linregress(g.FIQ, g.func_mean_fd).slope, pcr_slope(g)
        jk.append({"dropped_site": s, "ols": o, "pcr": pc, "ratio": o / pc})
    jk = pd.DataFrame(jk); jk.to_csv(OUT / "stress_s5_site_jackknife.csv", index=False)
    res["S5"] = {"ratio_min": float(jk.ratio.min()), "ratio_max": float(jk.ratio.max()),
                 "site_at_min": jk.loc[jk.ratio.idxmin(), "dropped_site"],
                 "site_at_max": jk.loc[jk.ratio.idxmax(), "dropped_site"]}

    # S6: heteroscedasticity: site SD tracks site mean (right-skewed FD)
    r6 = spearmanr(sites.mean_fd, sites.sd_fd)
    res["S6"] = {"spearman_site_sd_vs_mean": float(r6.correlation), "p": float(r6.pvalue),
                 "fd_skewness": float(df.func_mean_fd.skew()),
                 "logfd_skewness": float(np.log(df.func_mean_fd).skew())}
    sites.to_csv(OUT / "stress_s2_site_slopes.csv", index=False)

    (OUT / "stress_tests.json").write_text(json.dumps(res, indent=2, default=str))
    print(json.dumps(res, indent=2, default=str))


if __name__ == "__main__":
    main()

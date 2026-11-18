r"""EXP-06: ABIDE-II replication, executed against analysis_plan.md (frozen first).

Reuses the EXP-04 and EXP-05 functions unchanged. Run from repo root:
  python experiments/EXP-06_Abide2Replication_2026-11-12/scripts/run_replication.py --bootstrap 1000 --n-jobs 8
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
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
sys.path.insert(0, str(ROOT / "experiments" / "EXP-04_FrameLevelGStudy_2026-11-12" / "scripts"))
sys.path.insert(0, str(ROOT / "experiments" / "EXP-05_RobustnessChecks_2026-11-12" / "scripts"))
warnings.filterwarnings("ignore")
import run_gstudy as g4  # noqa: E402
import run_robustness as r5  # noqa: E402
from utils import filter_abide_sample, load_abide_csv, resolve_input_path, sigma_x_from_age  # noqa: E402

D2 = Path(os.environ.get("DEV_DATA_ROOT", r"D:\Dev\Data")) / "eiv-abide-replication" / "abide2"
SEED = 20261112


def num(s):
    return pd.to_numeric(s, errors="coerce")


def load_abide2() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    man = json.loads((D2 / "manifest.json").read_text())
    files = {r["sub"]: sorted(f["name"] for f in r["files"]) for r in man["records"]}
    frames = []
    for p in sorted((D2 / "participants").glob("*.tsv")):
        t = pd.read_csv(p, sep="\t", dtype=str, encoding="latin-1")  # some tables are not valid UTF-8
        t.columns = [c.strip().lower() for c in t.columns]
        t["collection"] = p.stem
        frames.append(t)
    ph = pd.concat(frames, ignore_index=True)
    flow = {"participants_listed": len(ph), "longitudinal_excluded": int(ph.collection.str.contains("Long").sum())}
    ph = ph[~ph.collection.str.contains("Long")]
    d = pd.DataFrame({"SITE_ID": ph["collection"].str.replace("ABIDEII-", "", regex=False),
                      "pid": "sub-" + ph["participant_id"].str.strip().str.replace("sub-", "", regex=False),
                      "FIQ": num(ph["fiq"]), "AGE_AT_SCAN": num(ph["age_at_scan"]),
                      "SEX": num(ph["sex"]), "DX_GROUP": num(ph["dx_group"])})
    ok = d.FIQ.notna() & (d.FIQ > 40) & d.AGE_AT_SCAN.notna() & d.SEX.isin([1, 2]) & d.DX_GROUP.isin([1, 2])
    flow["with_iq_age_sex_dx"] = int(ok.sum())
    d = d[ok].drop_duplicates("pid").reset_index(drop=True)
    subj, runs = [], []
    for pos, pid in enumerate(d.pid):
        for j, name in enumerate(files.get(pid, [])):
            fd = pd.read_csv(D2 / "fmriprep_confounds" / name, sep="\t", usecols=["framewise_displacement"]).iloc[:, 0].dropna().to_numpy(float)
            if len(fd) < 40:
                continue
            ses = (re.search(r"_ses-([^_]+)", name) or [None, "1"])[1]
            runs.append({"pos": pos, "ses": ses, "name_order": j, "n_frames": len(fd), "mean_fd": fd.mean()})
            if not any(s["pos"] == pos for s in subj):
                rec = {"pos": pos, "n_frames": len(fd), "fp_fd": fd.mean(),
                       "noise_sd": g4.seg_means(fd, g4.NOISE_BLOCKS).std(ddof=1) / np.sqrt(g4.NOISE_BLOCKS)}
                for k in g4.KS:
                    for i, v in enumerate(g4.seg_means(fd, k)):
                        rec[f"k{k}_s{i}"] = v
                subj.append(rec)
    s = pd.DataFrame(subj).join(d.drop(columns="pid"), on="pos")
    flow["with_frames"] = len(s)
    n_site = s.groupby("SITE_ID").pos.transform("count")
    s = s[n_site >= 5].reset_index(drop=True)
    flow["final_n"] = len(s)
    flow["final_sites"] = int(s.SITE_ID.nunique())
    flow["n_by_site"] = s.SITE_ID.value_counts().to_dict()
    s["sigma_x"] = sigma_x_from_age(s["AGE_AT_SCAN"])
    r = pd.DataFrame(runs)
    r = r[r.pos.isin(s.pos)].join(d[["SITE_ID"]], on="pos")
    return s, r, flow


def icc1(v: pd.Series, grp: pd.Series) -> float:
    g = v.groupby(grp)
    k0 = g.count().mean()
    msb = (g.count() * (g.mean() - v.mean()) ** 2).sum() / (g.ngroups - 1)
    msw = ((v - g.transform("mean")) ** 2).sum() / (len(v) - g.ngroups)
    return float((msb - msw) / (msb + (k0 - 1) * msw))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bootstrap", type=int, default=1000)
    ap.add_argument("--n-jobs", type=int, default=1)
    args = ap.parse_args()
    t0 = time.time()
    out = EXP / "results"
    plan_hash = hashlib.sha256((EXP / "analysis_plan.md").read_bytes()).hexdigest()
    d, runs, flow = load_abide2()
    res: dict = {"plan_sha256": plan_hash, "flow": flow,
                 "descriptives": {"mean_fiq": float(d.FIQ.mean()), "mean_age": float(d.AGE_AT_SCAN.mean()),
                                  "n_asd": int((d.DX_GROUP == 1).sum()), "n_td": int((d.DX_GROUP == 2).sum()),
                                  "mean_fd": float(d.fp_fd.mean()), "median_frames": float(d.n_frames.median())}}

    # H1 + variance components
    g2 = pd.DataFrame([g4.g_components(d, k, lg) for lg in (False, True) for k in g4.KS])
    g2.to_csv(out / "g_variance_components.csv", index=False)
    pd.concat([g4.per_site(d, 4, False), g4.per_site(d, 4, True)]).to_csv(out / "g_per_site.csv", index=False)
    h1 = float(g2[(g2.scale == "log") & (g2.k == 2)].g_coefficient_run_mean.iloc[0])

    # H2 controlled log slope, wild cluster bootstrap
    d["log_fd"] = np.log(d.fp_fd)
    h2 = {}
    for scale, y in (("log", d.log_fd), ("raw", d.fp_fd)):
        xr, yr, k = r5.fwl(d, y)
        t, p = r5.wild_cluster_p(xr, yr, d.SITE_ID.to_numpy())
        h2[scale] = {"slope": float(xr @ yr / (xr @ xr)), "wild_cluster_t": t, "wild_cluster_p": p}

    # H3-H5 point estimates and bootstrap
    point = pd.DataFrame([{"scale": "log" if lg else "raw", "controls": c, **g4.fits(d, lg, c)}
                          for lg in (False, True) for c in (False, True)])
    point.to_csv(out / "slopes_point.csv", index=False)
    sites = np.array(sorted(d.SITE_ID.unique()))
    seeds = np.random.SeedSequence(SEED).generate_state(args.bootstrap)
    cols = ["SITE_ID", "FIQ", "fp_fd", "noise_sd", "sigma_x", "AGE_AT_SCAN", "SEX", "DX_GROUP"]
    reps = Parallel(n_jobs=args.n_jobs, backend="loky")(delayed(g4.boot_rep)(d[cols], sites, int(s)) for s in seeds)
    boot = pd.DataFrame([r for rep in reps for r in rep])
    est = ["ols", "rc", "pcr2_measured", "hazy_york_measured", "pcr2_site_sd", "york_measured"]
    long = []
    for _, row in point.iterrows():
        b = boot[(boot.scale == row.scale) & (boot.controls == row.controls)]
        for e in est:
            s = b[e].dropna()
            long.append({"scale": row.scale, "controls": row.controls, "estimator": e, "slope": row[e],
                         "ci_low": s.quantile(.025), "ci_high": s.quantile(.975), "share_negative": float((s < 0).mean()),
                         "n_boot_ok": int(len(s)), "n_boot": int(len(b))})
    slopes = pd.DataFrame(long)
    slopes.to_csv(out / "slopes_with_intervals.csv", index=False)

    def cell(scale, controls, e):
        return slopes[(slopes.scale == scale) & (slopes.controls == controls) & (slopes.estimator == e)].iloc[0]

    lc = point[(point.scale == "log") & (point.controls)].iloc[0]
    rn = point[(point.scale == "raw") & (~point.controls)].iloc[0]
    ols_lc, ols_rn = cell("log", True, "ols"), cell("raw", False, "ols")
    interp_lc = bool(ols_lc.ci_high < 0 or ols_lc.ci_low > 0)
    interp_rn = bool(ols_rn.ci_high < 0 or ols_rn.ci_low > 0)
    res["confirmatory"] = {
        "H1": {"g_log_k2": h1, "criterion": ">= 0.80", "pass": bool(h1 >= 0.80)},
        "H2": {**h2["log"], "criterion": "slope < 0 and wild cluster p < 0.05",
               "pass": bool(h2["log"]["slope"] < 0 and h2["log"]["wild_cluster_p"] < 0.05)},
        "H2b": {"slope": h2["log"]["slope"], "criterion": "inside [-0.0091, -0.0010]",
                "pass": bool(-0.0091 <= h2["log"]["slope"] <= -0.0010)},
        "H3": {"pcr2_measured": float(lc.pcr2_measured), "ols": float(lc.ols), "ratio": float(lc.pcr2_measured / lc.ols),
               "criterion": "0.90 to 1.10", "ols_interval_excludes_zero": interp_lc,
               "pass": bool(0.90 <= lc.pcr2_measured / lc.ols <= 1.10) if interp_lc else "not interpretable"},
        "H4": {"pcr2_site_sd": float(rn.pcr2_site_sd), "ols": float(rn.ols), "ratio_ols_over_pcr2": float(rn.ols / rn.pcr2_site_sd),
               "criterion": "> 1.5", "ols_interval_excludes_zero": interp_rn,
               "pass": bool(rn.ols / rn.pcr2_site_sd > 1.5) if interp_rn else "not interpretable"},
        "H5": {"york_reduced_chi2": float(lc.york_measured_reduced_chi2), "criterion": "> 5",
               "pass": bool(lc.york_measured_reduced_chi2 > 5)},
        "raw_scale_H2_for_reference": h2["raw"]}

    # Secondary
    f = "log_fd ~ FIQ * C(DX_GROUP) + AGE_AT_SCAN + I(AGE_AT_SCAN**2) + C(SEX) + C(SITE_ID)"
    m = smf.ols(f, d).fit(cov_type="cluster", cov_kwds={"groups": d.SITE_ID})
    m0 = smf.ols(f.replace("FIQ * C(DX_GROUP)", "FIQ + C(DX_GROUP)"), d).fit(cov_type="cluster", cov_kwds={"groups": d.SITE_ID})
    ik = [c for c in m.params.index if c.startswith("FIQ:")][0]
    dk = [c for c in m0.params.index if c.startswith("C(DX_GROUP)")][0]
    sec = {"dx_main_effect_td_minus_asd_log": float(m0.params[dk]), "dx_p_cluster": float(m0.pvalues[dk]),
           "iq_x_dx_interaction_log": float(m.params[ik]), "iq_x_dx_p_cluster": float(m.pvalues[ik])}
    multi = runs[runs.groupby("pos").mean_fd.transform("count") >= 2]
    if multi.pos.nunique() >= 10:
        sec["between_run"] = {"n_subjects": int(multi.pos.nunique()), "n_runs": int(len(multi)),
                              "sites": multi.groupby("SITE_ID").pos.nunique().to_dict(),
                              "sessions_seen": sorted(multi.ses.unique().tolist()),
                              "icc1_raw": icc1(multi.mean_fd, multi.pos), "icc1_log": icc1(np.log(multi.mean_fd), multi.pos)}
        # within-site version: remove site means first
        lw = np.log(multi.mean_fd) - np.log(multi.mean_fd).groupby(multi.SITE_ID).transform("mean")
        sec["between_run"]["icc1_log_within_site"] = icc1(lw, multi.pos)

    # pooled ABIDE-I + ABIDE-II (frame-level FD in both)
    base = filter_abide_sample(load_abide_csv(resolve_input_path(None))).reset_index(drop=True)
    base["sigma_x"] = sigma_x_from_age(base["AGE_AT_SCAN"])
    s1, _ = g4.load_frames(base)
    a1 = s1.join(base[["SITE_ID", "FIQ", "sigma_x", "AGE_AT_SCAN", "SEX", "DX_GROUP"]], on="pos")
    a1["SITE_ID"] = "A1_" + a1.SITE_ID
    a2 = d.copy()
    a2["SITE_ID"] = "A2_" + a2.SITE_ID
    pool = pd.concat([a1[cols], a2[cols]], ignore_index=True)
    xr, yr, _ = r5.fwl(pool, np.log(pool.fp_fd))
    t, p = r5.wild_cluster_p(xr, yr, pool.SITE_ID.to_numpy())
    xr_r, yr_r, _ = r5.fwl(pool, pool.fp_fd)
    tr, pr = r5.wild_cluster_p(xr_r, yr_r, pool.SITE_ID.to_numpy())
    pf = g4.fits(pool, True, True)
    sec["pooled_abide1_abide2"] = {"n": len(pool), "sites": int(pool.SITE_ID.nunique()),
                                   "log_slope_controls": float(xr @ yr / (xr @ xr)), "wild_cluster_p": p,
                                   "raw_slope_controls": float(xr_r @ yr_r / (xr_r @ xr_r)), "raw_wild_cluster_p": pr,
                                   "log_rc": float(pf["rc"]), "log_pcr2_measured": float(pf["pcr2_measured"])}
    res["secondary"] = sec

    for q in out.glob("*.csv"):
        assert not set(pd.read_csv(q, nrows=0).columns) & {"SUB_ID", "FILE_ID", "subject", "pos", "pid", "participant_id"}
    res["receipt"] = {"abide2_manifest": hashlib.sha256((D2 / "manifest.json").read_bytes()).hexdigest(),
                      "script": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                      "exp04_script": hashlib.sha256(Path(g4.__file__).read_bytes()).hexdigest(),
                      "exp05_script": hashlib.sha256(Path(r5.__file__).read_bytes()).hexdigest(),
                      "bootstrap": {"reps": args.bootstrap, "seed": SEED}, "run_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                      "runtime_s": round(time.time() - t0, 1)}
    (out / "summary.json").write_text(json.dumps(res, indent=2, default=str))
    print(json.dumps(res, indent=1, default=str))
    with pd.option_context("display.width", 220, "display.float_format", "{:.5g}".format):
        print(g2.to_string(index=False)); print(slopes.to_string(index=False))


if __name__ == "__main__":
    main()

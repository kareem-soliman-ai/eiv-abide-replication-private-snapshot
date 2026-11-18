r"""EXP-04: G study of head motion from frame-level FD, and reanalysis with measured noise.

Input: fMRIPrep confounds tables (per-frame framewise_displacement) on D:, see
download_confounds.py, plus the ABIDE-I phenotypic CSV.

G1 Validation: fMRIPrep run-mean FD vs the phenotypic func_mean_fd.
G2 G study, design (persons : sites) x segments. Each run is cut into K
   contiguous segments (K = 2, 4, 8). Variance components: site, person:site,
   person x segment residual. G coefficient for the full-run mean:
   var_p / (var_p + var_e / K). Raw and log scale; pooled and per site.
G3 Between-run reliability where a subject has several runs (OHSU): ICC(1).
G4 Per-subject noise SD of the run mean from batch means (8 blocks).
G5 Reanalysis: IQ slope with measured per-subject sigma_y versus the paper's
   site-SD sigma_y, same subjects, with and without controls, raw and log,
   site-cluster bootstrap.

Run from repo root:
  python experiments/EXP-04_FrameLevelGStudy_2026-11-12/scripts/run_gstudy.py --bootstrap 1000 --n-jobs 8
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
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
sys.path.insert(0, str(ROOT / "experiments" / "EXP-02_EIVMethodPanel_2026-11-12" / "vendor"))
warnings.filterwarnings("ignore")
import eiv_methods as em  # noqa: E402
import pcr2_reference as ref  # noqa: E402
from eiv_core import fit_eiv, measurement_covariances  # noqa: E402
from utils import filter_abide_sample, load_abide_csv, resolve_input_path, sigma_x_from_age  # noqa: E402

DATA = Path(os.environ.get("DEV_DATA_ROOT", r"D:\Dev\Data")) / "eiv-abide-replication" / "fmriprep_confounds"
KS = (2, 4, 8)
NOISE_BLOCKS = 8


def seg_means(fd: np.ndarray, k: int) -> np.ndarray:
    return np.array([s.mean() for s in np.array_split(fd, k)])


def load_frames(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return per-subject table (run 1) and per-run table. Keyed by row position, never by subject ID."""
    manifest = json.loads((DATA / "manifest.json").read_text())
    files = {int(r["sub"][4:]): sorted(f["name"] for f in r["files"]) for r in manifest["records"]}
    subj, runs = [], []
    for pos, row in enumerate(df.itertuples(index=False)):
        names = files.get(int(row.SUB_ID), [])
        for j, name in enumerate(names):
            fd = pd.read_csv(DATA / name, sep="\t", usecols=["framewise_displacement"])["framewise_displacement"]
            fd = fd.dropna().to_numpy(float)
            if len(fd) < 40:
                continue
            runs.append({"pos": pos, "run": j + 1, "n_frames": len(fd), "mean_fd": fd.mean()})
            if j == 0:
                rec = {"pos": pos, "n_frames": len(fd), "fp_fd": fd.mean(),
                       "noise_sd": seg_means(fd, NOISE_BLOCKS).std(ddof=1) / np.sqrt(NOISE_BLOCKS)}
                for k in KS:
                    for i, v in enumerate(seg_means(fd, k)):
                        rec[f"k{k}_s{i}"] = v
                subj.append(rec)
    return pd.DataFrame(subj), pd.DataFrame(runs)


def g_components(d: pd.DataFrame, k: int, log: bool) -> dict:
    """ANOVA variance components for (persons:sites) x segments, balanced in segments."""
    Y = d[[f"k{k}_s{i}" for i in range(k)]].to_numpy(float)
    if log:
        Y = np.log(Y)
    site = d["SITE_ID"].to_numpy()
    # remove site x segment means, so person and residual are within site
    R = Y.copy()
    for s in np.unique(site):
        m = site == s
        R[m] -= R[m].mean(axis=0, keepdims=True)
    n, n_sites = len(R), len(np.unique(site))
    pm = R.mean(axis=1)
    df_p = n - n_sites
    ms_p = k * np.sum(pm**2) / df_p
    ms_e = np.sum((R - pm[:, None]) ** 2) / (df_p * (k - 1))
    var_e = ms_e
    var_p = max((ms_p - ms_e) / k, 0.0)
    person_mean = Y.mean(axis=1)
    mm = smf.mixedlm("y ~ 1", pd.DataFrame({"y": person_mean, "g": site}), groups="g").fit(reml=True)
    var_site = float(mm.cov_re.iloc[0, 0])
    return {"k": k, "scale": "log" if log else "raw", "var_site": var_site, "var_person": var_p,
            "var_person_x_segment": var_e, "noise_var_of_run_mean": var_e / k,
            "g_coefficient_run_mean": var_p / (var_p + var_e / k),
            "noise_share_of_within_site_variance": (var_e / k) / (var_p + var_e / k),
            "site_share_of_total": var_site / (var_site + var_p + var_e / k)}


def per_site(d: pd.DataFrame, k: int, log: bool) -> pd.DataFrame:
    rows = []
    for s, g in d.groupby("SITE_ID"):
        if len(g) < 5:
            continue
        Y = g[[f"k{k}_s{i}" for i in range(k)]].to_numpy(float)
        if log:
            Y = np.log(Y)
        Y = Y - Y.mean(axis=0, keepdims=True)
        pm = Y.mean(axis=1)
        ms_p = k * np.sum(pm**2) / (len(g) - 1)
        ms_e = np.sum((Y - pm[:, None]) ** 2) / ((len(g) - 1) * (k - 1))
        var_p = max((ms_p - ms_e) / k, 0.0)
        rows.append({"site_id": s, "n": len(g), "scale": "log" if log else "raw",
                     "site_sd_of_run_mean": float(np.sqrt(var_p + ms_e / k)),
                     "true_person_sd": float(np.sqrt(var_p)), "noise_sd": float(np.sqrt(ms_e / k)),
                     "g_coefficient": var_p / (var_p + ms_e / k)})
    return pd.DataFrame(rows)


def design(df: pd.DataFrame, controls: bool) -> np.ndarray:
    cols = [np.ones(len(df))]
    if controls:
        a = df["AGE_AT_SCAN"].to_numpy(float)
        cols += [pd.get_dummies(df["SITE_ID"], drop_first=True).to_numpy(float), a, a**2,
                 (df["SEX"] == 2).to_numpy(float), (df["DX_GROUP"] == 1).to_numpy(float)]
    return np.column_stack(cols)


def resid(v, Z):
    return v - Z @ np.linalg.lstsq(Z, v, rcond=None)[0]


def fits(df: pd.DataFrame, log: bool, controls: bool) -> dict:
    """IQ slope by estimator, for measured per-subject sigma_y and for the paper's site-SD sigma_y."""
    y0 = np.log(df["fp_fd"].to_numpy(float)) if log else df["fp_fd"].to_numpy(float)
    sy_meas = (df["noise_sd"] / df["fp_fd"]).to_numpy(float) if log else df["noise_sd"].to_numpy(float)
    sy_site = pd.Series(y0).groupby(df["SITE_ID"].to_numpy()).transform("std")
    # a site left with one subject has no SD; give it the median site SD and say so in the report
    sy_site = sy_site.fillna(sy_site.median()).to_numpy()
    Z = design(df, controls)
    n, kz = Z.shape
    x, y = resid(df["FIQ"].to_numpy(float), Z), resid(y0, Z)
    sx = df["sigma_x"].to_numpy(float)
    sxx, sxy = x @ x / (n - kz), x @ y / (n - kz)
    out = {"ols": sxy / sxx, "rc": sxy / (sxx - np.mean(sx**2))}
    for tag, sy in (("measured", sy_meas), ("site_sd", sy_site)):
        yk = em.fit_york(x, y, sx, sy, max_iterations=2000)
        out[f"york_{tag}"] = yk.slope if yk.converged and abs(yk.slope) < 1 else np.nan
        if tag == "measured":
            r = y - yk.intercept - yk.slope * x
            out["york_measured_reduced_chi2"] = float(np.sum(r**2 / (sy**2 + yk.slope**2 * sx**2)) / (n - 2))
        try:
            b, _, t2 = ref.hazy_york(x, y, sx, sy)
            out[f"hazy_york_{tag}"] = b
        except Exception:
            out[f"hazy_york_{tag}"] = np.nan
        try:
            X, yy, C = measurement_covariances(x, y, sigma_X=sx, sigma_y=sy)
            f = fit_eiv(X, yy, C, method="structural", max_iter=2000)
            out[f"pcr2_{tag}"] = float(f.coef[0]) if f.converged else np.nan
            out[f"pcr2_{tag}_tau2"] = f.intrinsic_variance
        except Exception:
            out[f"pcr2_{tag}"] = np.nan
    return out


def boot_rep(df: pd.DataFrame, sites: np.ndarray, seed: int) -> list[dict]:
    rng = np.random.default_rng(seed)
    parts = []
    for j, s in enumerate(rng.choice(sites, len(sites), replace=True)):
        g = df[df.SITE_ID == s].copy()
        g["SITE_ID"] = f"{s}#{j}"
        parts.append(g)
    b = pd.concat(parts, ignore_index=True)
    return [{"seed": seed, "scale": "log" if lg else "raw", "controls": c, **fits(b, lg, c)}
            for lg in (False, True) for c in (False, True)]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bootstrap", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=20261112)
    ap.add_argument("--n-jobs", type=int, default=1, help="parallel workers (8 recommended here)")
    args = ap.parse_args()
    t0 = time.time()
    out = EXP / "results"
    path = resolve_input_path(None)
    if not (DATA / "manifest.json").exists():
        sys.exit(f"Frame data not found at {DATA}. Run download_confounds.py first.")
    base = filter_abide_sample(load_abide_csv(path)).reset_index(drop=True)
    base["sigma_x"] = sigma_x_from_age(base["AGE_AT_SCAN"])
    subj, runs = load_frames(base)
    keep = ["SITE_ID", "FIQ", "func_mean_fd", "sigma_x", "AGE_AT_SCAN", "SEX", "DX_GROUP"]
    d = subj.join(base[keep], on="pos")
    runs = runs.join(base[["SITE_ID"]], on="pos")
    res: dict = {"n_phenotypic_sample": len(base), "n_with_frames": len(d), "sites_with_frames": int(d.SITE_ID.nunique()),
                 "missing_by_site": base.loc[~base.index.isin(d.pos), "SITE_ID"].value_counts().to_dict(),
                 "frames_per_run": {"min": int(d.n_frames.min()), "median": float(d.n_frames.median()), "max": int(d.n_frames.max())}}

    # G1 validation
    res["G1_validation"] = {
        "pearson_raw": float(np.corrcoef(d.fp_fd, d.func_mean_fd)[0, 1]),
        "pearson_log": float(np.corrcoef(np.log(d.fp_fd), np.log(d.func_mean_fd))[0, 1]),
        "spearman": float(d[["fp_fd", "func_mean_fd"]].corr(method="spearman").iloc[0, 1]),
        "mean_fmriprep_fd": float(d.fp_fd.mean()), "mean_phenotypic_fd": float(d.func_mean_fd.mean())}

    # G2 variance components
    g2 = pd.DataFrame([g_components(d, k, lg) for lg in (False, True) for k in KS])
    g2.to_csv(out / "g2_variance_components.csv", index=False)
    ps = pd.concat([per_site(d, 4, False), per_site(d, 4, True)])
    ps.to_csv(out / "g2_per_site.csv", index=False)
    sh = {}
    for lg in (False, True):
        A = d[["k2_s0", "k2_s1"]].to_numpy(float)
        A = np.log(A) if lg else A.copy()
        for s in d.SITE_ID.unique():
            m = (d.SITE_ID == s).to_numpy()
            A[m] -= A[m].mean(axis=0)
        r = float(np.corrcoef(A[:, 0], A[:, 1])[0, 1])
        sh["log" if lg else "raw"] = {"split_half_r_within_site": r, "spearman_brown": 2 * r / (1 + r)}
    res["G2_split_half"] = sh

    # G3 between-run ICC(1), subjects with >= 2 runs
    multi = runs[runs.groupby("pos").run.transform("count") >= 2]
    g3 = {}
    for lg in (False, True):
        v = np.log(multi.mean_fd) if lg else multi.mean_fd
        grp = v.groupby(multi.pos)
        k0 = grp.count().mean()
        msb = (grp.count() * (grp.mean() - v.mean()) ** 2).sum() / (grp.ngroups - 1)
        msw = ((v - grp.transform("mean")) ** 2).sum() / (len(v) - grp.ngroups)
        g3["log" if lg else "raw"] = {"icc1_single_run": float((msb - msw) / (msb + (k0 - 1) * msw))}
    g3.update(n_subjects=int(multi.pos.nunique()), n_runs=int(len(multi)), sites=multi.SITE_ID.value_counts().to_dict())
    res["G3_between_run"] = g3

    # G4 measured noise vs the paper's sigma_y
    site_sd = d.groupby("SITE_ID").fp_fd.transform("std")
    res["G4_noise"] = {"median_measured_noise_sd": float(d.noise_sd.median()), "mean_measured_noise_sd": float(d.noise_sd.mean()),
                       "mean_paper_style_site_sd": float(site_sd.mean()),
                       "median_ratio_measured_over_site_sd": float((d.noise_sd / site_sd).median()),
                       "spearman_noise_sd_vs_mean_fd": float(d[["noise_sd", "fp_fd"]].corr(method="spearman").iloc[0, 1])}

    # G5 reanalysis
    point = pd.DataFrame([{"scale": "log" if lg else "raw", "controls": c, **fits(d, lg, c)}
                          for lg in (False, True) for c in (False, True)])
    pheno = []
    for lg in (False, True):
        for c in (False, True):
            Z = design(d, c)
            y = np.log(d.func_mean_fd.to_numpy(float)) if lg else d.func_mean_fd.to_numpy(float)
            x, yr = resid(d.FIQ.to_numpy(float), Z), resid(y, Z)
            pheno.append({"scale": "log" if lg else "raw", "controls": c, "ols_phenotypic_fd_same_subjects": x @ yr / (x @ x)})
    point = point.merge(pd.DataFrame(pheno), on=["scale", "controls"])
    sites = np.array(sorted(d.SITE_ID.unique()))
    seeds = np.random.SeedSequence(args.seed).generate_state(args.bootstrap)
    cols = ["SITE_ID", "FIQ", "fp_fd", "noise_sd", "sigma_x", "AGE_AT_SCAN", "SEX", "DX_GROUP"]
    reps = Parallel(n_jobs=args.n_jobs, backend="loky")(delayed(boot_rep)(d[cols], sites, int(s)) for s in seeds)
    boot = pd.DataFrame([r for rep in reps for r in rep])
    boot.to_csv(out / "g5_bootstrap_replicates.csv", index=False)
    est = [c for c in boot.columns if c not in ("seed", "scale", "controls") and "tau2" not in c and "chi2" not in c]
    long = []
    for _, row in point.iterrows():
        b = boot[(boot.scale == row.scale) & (boot.controls == row.controls)]
        for e in est:
            s = b[e].dropna()
            long.append({"scale": row.scale, "controls": row.controls, "estimator": e, "slope": row[e],
                         "ci_low": s.quantile(.025) if len(s) else np.nan, "ci_high": s.quantile(.975) if len(s) else np.nan,
                         "share_negative": float((s < 0).mean()) if len(s) else np.nan,
                         "n_boot_ok": int(len(s)), "n_boot": int(len(b))})
    g5 = pd.DataFrame(long)
    g5.to_csv(out / "g5_slopes.csv", index=False)
    point.to_csv(out / "g5_point_details.csv", index=False)

    for p in out.glob("*.csv"):
        assert not set(pd.read_csv(p, nrows=0).columns) & {"SUB_ID", "FILE_ID", "subject", "pos"}, p.name
    h = lambda q: hashlib.sha256(Path(q).read_bytes()).hexdigest()  # noqa: E731
    res["receipt"] = {"phenotypic_csv": h(path), "confounds_manifest": h(DATA / "manifest.json"),
                      "eiv_core.py": h(ROOT / "src/eiv_core.py"), "script": h(__file__),
                      "bootstrap": {"reps": args.bootstrap, "seed": args.seed, "unit": "site"},
                      "runtime_s": round(time.time() - t0, 1)}
    (out / "summary.json").write_text(json.dumps(res, indent=2))
    with pd.option_context("display.width", 250, "display.max_rows", 200, "display.max_columns", 30,
                           "display.float_format", "{:.5g}".format):
        print(json.dumps({k: v for k, v in res.items() if k != "receipt"}, indent=1))
        print(g2.to_string(index=False)); print(ps.to_string(index=False))
        print(point.T.to_string()); print(g5.to_string(index=False))
    print(f"runtime {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()

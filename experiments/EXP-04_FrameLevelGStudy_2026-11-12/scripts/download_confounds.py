r"""Download fMRIPrep per-frame confounds tables for the ABIDE-I analysis sample.

Source: public bucket s3://fcp-indi/data/Projects/ABIDE/Outputs/fmriprep/fmriprep/
Files:  sub-<ID>/func/*_desc-confounds_regressors.tsv (all runs)
Target: D:\Dev\Data\eiv-abide-replication\fmriprep_confounds (override with DEV_DATA_ROOT)
Resumable: existing non-empty files are skipped. Writes a manifest with sizes and SHA-256.
"""
from __future__ import annotations
import hashlib, json, os, re, sys, urllib.request, urllib.parse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from utils import filter_abide_sample, load_abide_csv, resolve_input_path  # noqa: E402

BUCKET = "https://s3.amazonaws.com/fcp-indi"
PREFIX = "data/Projects/ABIDE/Outputs/fmriprep/fmriprep/"
OUT = Path(os.environ.get("DEV_DATA_ROOT", r"D:\Dev\Data")) / "eiv-abide-replication" / "fmriprep_confounds"


def get(url: str, tries: int = 4) -> bytes:
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read()
        except Exception:
            if i == tries - 1:
                raise
    return b""


def one(sub_id: int) -> dict:
    sub = f"sub-{sub_id:07d}"
    listing = get(f"{BUCKET}?prefix={urllib.parse.quote(PREFIX + sub + '/func/')}&max-keys=1000").decode()
    keys = [k for k in re.findall(r"<Key>([^<]+)</Key>", listing) if k.endswith("_desc-confounds_regressors.tsv")]
    rec = {"sub": sub, "files": []}
    for k in keys:
        dest = OUT / Path(k).name
        if not (dest.exists() and dest.stat().st_size > 0):
            data = get(f"{BUCKET}/{urllib.parse.quote(k)}")
            tmp = dest.with_suffix(".part")
            tmp.write_bytes(data)
            tmp.replace(dest)
        b = dest.read_bytes()
        rec["files"].append({"name": dest.name, "bytes": len(b), "sha256": hashlib.sha256(b).hexdigest()})
    return rec


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    ids = filter_abide_sample(load_abide_csv(resolve_input_path(None)))["SUB_ID"].astype(int).tolist()
    with ThreadPoolExecutor(16) as ex:
        recs = list(ex.map(one, ids))
    missing = [r["sub"] for r in recs if not r["files"]]
    n_files = sum(len(r["files"]) for r in recs)
    total = sum(f["bytes"] for r in recs for f in r["files"])
    (OUT / "manifest.json").write_text(json.dumps({"source": BUCKET + "/" + PREFIX, "subjects_requested": len(ids),
        "subjects_with_files": len(ids) - len(missing), "n_files": n_files, "total_bytes": total,
        "missing": missing, "records": recs}, indent=1))
    print(f"requested={len(ids)} with_files={len(ids)-len(missing)} files={n_files} MB={total/1e6:.1f} missing={len(missing)}")
    assert n_files > 0, "nothing downloaded; aborting"


if __name__ == "__main__":
    main()

r"""Download ABIDE-II participant tables and fMRIPrep per-frame confounds tables.

Source: public bucket s3://fcp-indi/data/Projects/ABIDE2/
Target: D:\Dev\Data\eiv-abide-replication\abide2 (override with DEV_DATA_ROOT)
Resumable; writes manifest.json with sizes and SHA-256. Only confounds tables are fetched.
"""
from __future__ import annotations
import csv, hashlib, io, json, os, re, urllib.parse, urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

BUCKET = "https://s3.amazonaws.com/fcp-indi"
RAW = "data/Projects/ABIDE2/RawData/"
FP = "data/Projects/ABIDE2/Outputs/fmriprep/fmriprep/"
OUT = Path(os.environ.get("DEV_DATA_ROOT", r"D:\Dev\Data")) / "eiv-abide-replication" / "abide2"


def get(url: str, tries: int = 4) -> bytes:
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read()
        except Exception:
            if i == tries - 1:
                raise
    return b""


def list_keys(prefix: str, delimiter: str = "") -> tuple[list[str], list[str]]:
    keys, prefixes, marker = [], [], ""
    while True:
        q = f"{BUCKET}?prefix={urllib.parse.quote(prefix)}&max-keys=1000"
        if delimiter:
            q += "&delimiter=" + delimiter
        if marker:
            q += "&marker=" + urllib.parse.quote(marker)
        x = get(q).decode()
        k = re.findall(r"<Key>([^<]+)</Key>", x)
        keys += k
        prefixes += re.findall(r"<Prefix>([^<]+)</Prefix>", x)[1 if not marker else 0:]
        if "<IsTruncated>true" not in x:
            return keys, prefixes
        nm = re.findall(r"<NextMarker>([^<]+)</NextMarker>", x)
        marker = nm[0] if nm else (k[-1] if k else "")
        if not marker:
            return keys, prefixes


def fetch_subject(sub: str) -> dict:
    keys, _ = list_keys(FP + sub + "/")
    rec = {"sub": sub, "files": []}
    for k in keys:
        if not k.endswith("_desc-confounds_regressors.tsv"):
            continue
        dest = OUT / "fmriprep_confounds" / Path(k).name
        if not (dest.exists() and dest.stat().st_size > 0):
            tmp = dest.with_suffix(".part")
            tmp.write_bytes(get(f"{BUCKET}/{urllib.parse.quote(k)}"))
            tmp.replace(dest)
        b = dest.read_bytes()
        rec["files"].append({"name": dest.name, "bytes": len(b), "sha256": hashlib.sha256(b).hexdigest()})
    return rec


def main() -> None:
    (OUT / "fmriprep_confounds").mkdir(parents=True, exist_ok=True)
    (OUT / "participants").mkdir(parents=True, exist_ok=True)
    _, site_prefixes = list_keys(RAW, "/")
    sites = [p.rstrip("/").split("/")[-1] for p in site_prefixes if "ABIDEII-" in p]
    subs: dict[str, str] = {}
    tables = {}
    for s in sites:
        try:
            data = get(f"{BUCKET}/{RAW}{s}/participants.tsv")
        except Exception:
            tables[s] = "missing participants.tsv"
            continue
        (OUT / "participants" / f"{s}.tsv").write_bytes(data)
        rows = list(csv.DictReader(io.StringIO(data.decode("utf-8", "replace")), delimiter="\t"))
        tables[s] = len(rows)
        for r in rows:
            pid = (r.get("participant_id") or "").strip().replace("sub-", "")
            if pid:
                subs.setdefault(f"sub-{pid}", s)
    with ThreadPoolExecutor(16) as ex:
        recs = list(ex.map(fetch_subject, sorted(subs)))
    for r in recs:
        r["collection"] = subs[r["sub"]]
    n_files = sum(len(r["files"]) for r in recs)
    total = sum(f["bytes"] for r in recs for f in r["files"])
    (OUT / "manifest.json").write_text(json.dumps({"source": BUCKET + "/data/Projects/ABIDE2/", "participant_tables": tables,
        "subjects_listed": len(subs), "subjects_with_files": sum(bool(r["files"]) for r in recs),
        "n_files": n_files, "total_bytes": total, "records": recs}, indent=1))
    print(json.dumps(tables), f"\nsubjects={len(subs)} with_files={sum(bool(r['files']) for r in recs)} files={n_files} MB={total/1e6:.1f}")
    assert n_files > 0, "nothing downloaded; aborting"


if __name__ == "__main__":
    main()

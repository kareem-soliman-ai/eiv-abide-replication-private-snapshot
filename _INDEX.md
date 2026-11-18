---
project: EIV ABIDE Replication
short_name: EIV-ABIDE
type: research
status: released
phase: self-contained public replication bundle
visibility: public
license: CC BY-NC-SA 4.0
updated: 2026-04-14
---

# EIV ABIDE Replication — Map of Content

*A self-contained public replication of the ABIDE-I phenotypic analysis that introduces Probability Cloud Regression as a novel errors-in-variables method for the IQ-motion confound in multi-site autism fMRI.*

**One-sentence frame:** This repository reimplements only the published equations needed to estimate the pooled OLS slope, the EIV-corrected PCR slope, the leave-site-out cross-validation result, the within-tier slope analysis, and the 8x8 sensitivity grid, so that a reader can reproduce the paper's headline results from a single ABIDE-I phenotypic CSV.

**Portfolio role:** public-facing companion bundle demonstrating the Probability Cloud Regression method on a real neuroimaging dataset. Provides external credibility for the private PCR research programme. See [[Projects/portfolio/ARCHITECTURE|portfolio/ARCHITECTURE]].

> No canonical `docs/` folder is present yet in this repository — the README is the single entry point. If docs are added later, register them under `docs/core/` and update this MOC.

---

## Current state at a glance

- **Purpose:** reproduce the ABIDE-I phenotypic analysis from the paper end-to-end, showing that the IQ-motion confound in multi-site autism fMRI may be inflated by site-correlated measurement uncertainty.
- **Scope:** filtered sample of n = 935 across 19 sites; OLS, PCR, EM convergence trace, tier slopes, LOSO cross-validation, and 8x8 sensitivity grid.
- **Novel method:** introduces Probability Cloud Regression (PCR) as an errors-in-variables correction for the confound.
- **Entry point:** `src/run_analysis.py`. Privacy-safe CSV outputs land under `results/` and figures under `figures/`.
- **License:** CC BY-NC-SA 4.0
- **Visibility:** Public.

---

## Core documentation

| Document | Purpose |
|----------|---------|
| [[Projects/projects/public/eiv-abide-replication/README]] | Public introduction, setup, usage, expected output |
| [[Projects/projects/public/eiv-abide-replication/data/README]] | Instructions for obtaining the ABIDE-I phenotypic CSV |
| [[CITATION]] | Citation metadata |

---

## Results (headline CSVs)

```dataview
TABLE WITHOUT ID
  file.link as "File",
  dateformat(file.mtime, "yyyy-MM-dd") as "Edited"
FROM "Projects/projects/public/eiv-abide-replication/results"
SORT file.name ASC
```

---

## Figures

```dataview
LIST
FROM "Projects/projects/public/eiv-abide-replication/figures"
SORT file.name ASC
```

---

## Recently edited files in this project (last 20)

```dataview
TABLE WITHOUT ID
  file.link as "File",
  regexreplace(file.folder, "Projects/projects/public/eiv-abide-replication/?", "") as "Location",
  dateformat(file.mtime, "yyyy-MM-dd HH:mm") as "Edited"
FROM "Projects/projects/public/eiv-abide-replication"
WHERE !contains(file.folder, "archive") AND !contains(file.folder, "__pycache__") AND !contains(file.folder, ".tmp") AND !contains(file.folder, ".logs")
SORT file.mtime DESC
LIMIT 20
```

---

## Code layout reference

Project root: `projects/public/eiv-abide-replication/`

- `src/run_analysis.py` — entry point running the full pipeline
- `src/pcr.py` — Probability Cloud Regression implementation (EM fit, corrected slope)
- `src/cross_validation.py` — leave-site-out cross-validation
- `src/sensitivity_grid.py` — 8x8 sensitivity grid across reliability assumptions
- `src/utils.py` — shared helpers
- `data/` — location for the ABIDE-I phenotypic CSV (not redistributed)
- `results/` — privacy-safe CSV outputs
- `figures/` — generated PNG figures
- `requirements.txt` — Python dependencies
- `LICENSE`, `CITATION.cff`, `CODE_OF_CONDUCT.md`, `CONTRIBUTING.md`, `SECURITY.md` — repo metadata

---

## Related portfolio assets

- **Public positioning:** [[Projects/portfolio/ARCHITECTURE]] — this bundle is the public evidence anchor for the PCR method.
- **Private programme:** Probability Cloud Regression research programme (private repo) — this bundle is the external slice; the full method, findings catalogue, and extensions live there.
- **Flagship findings destination:** `portfolio/identity/FLAGSHIP_FINDINGS.md` — the ABIDE-I replication is a candidate entry.

---

## Intellectual neighbours

- **Probability Cloud Regression (private)** — the full PCR research programme; this repo is its audited public replication.
- **Trajectory Geometry** — uses PCR variance forensics inside its EXP-VF experiments.
- **PCR-FIN** — finance-domain extension of PCR; separate private repo.

---

## Working notes

- This repository is read-mostly. Changes to the PCR method proper should happen in the private PCR repo first, then be ported here as audited updates.
- Privacy-safe outputs only: never commit subject identifiers such as `SUB_ID` or `FILE_ID` into `results/`.
- Re-running the sensitivity grid can be expensive; prefer caching intermediate outputs rather than re-running from scratch.

---

*MOC last updated: 2026-04-14. Update the frontmatter `updated` and `phase` fields when a new audited snapshot lands.*

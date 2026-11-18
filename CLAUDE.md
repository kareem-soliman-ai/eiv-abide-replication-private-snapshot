# eiv-abide-replication — Repo Context for AI Agents

<!-- FIRASA:ECOSYSTEM-POINTER:BEGIN -->
## Firasa ecosystem context (auto-maintained — do not hand-edit inside these markers)

This repo is one project inside the **Munazzama**, Kareem's twelve-officer AI organisation.
The ecosystem's source of truth is `C:\Dev\firasa` (see its `CLAUDE.md`).

**At session start, before reading files or re-deriving decisions:**
1. `kareem_session_context` (firasa-pm or firasa-memory) — focus, current project, what is degraded.
2. `memory_search` (firasa-memory) — semantic recall over every past session on every platform.

**Where things live**
| Need | Source |
|---|---|
| Who does what, and what they can reach | `C:\Dev\firasa\ecosystem\desks\DESKS.md` |
| Live registries (agents, skills, MCP, crons, storage) | `C:\Dev\firasa\ecosystem\context\` |
| Why it is built this way | `C:\Dev\firasa\ecosystem\decisions\INDEX.md` |

**Retired — do not use or reintroduce**
- **Todoist.** Tasks live in this repo's `docs/core/ROADMAP.md` (the durable plan of record)
  and in Slack Lists (the live surface, owned by Munazzim). Flow is roadmap -> list.
  Any `todoist_*` tool fails deliberately. See `ADR-2026-08-07-task-source-of-truth.md`.
- **Handoff files.** Raw session logs auto-ingest into mem0; recall them instead.
- **`kareem-*` MCP servers.** Superseded by `firasa-*`.

**Experiment results:** `experiment_store` (firasa-pm) + `D:\Dev\Data\<project>\experiments\`.
Analyses and plots stay on `C:\`; bulk outputs go to `D:\` with a `data_pointer.txt` on `C:\`.

**Graceful degradation:** a missing optional dependency degrades and says so; it never raises
and never returns a silent empty result. Report the gap, name the remedy.
<!-- FIRASA:ECOSYSTEM-POINTER:END -->


Project-specific context. Read with the global `~/.claude/CLAUDE.md` (storage, file placement, conventions) and any `AGENT_CONTROL.md` at this root.

## What this is

Self-contained replication of the ABIDE-I phenotypic analysis from the paper "The IQ-Motion Confound in Multi-Site Autism fMRI May Be Inflated by Site-Correlated Measurement Uncertainty". It reimplements only the published equations needed to estimate the pooled OLS slope, the EIV-corrected PCR slope, the leave-site-out (LOSO) cross-validation, the within-tier slope analysis, and the 8×8 sensitivity grid.

Start with [README.md](README.md); the analysis entrypoint is [src/run_analysis.py](src/run_analysis.py).

## Stack

Python 3.9+ (numpy, pandas, scipy, matplotlib — see [requirements.txt](requirements.txt)). Source in `src/`: `pcr.py` (PCR implementation), `cross_validation.py`, `sensitivity_grid.py`, `utils.py`, `run_analysis.py`. Licence CC BY-NC-SA 4.0; full public-repo scaffolding present (`CONTRIBUTING.md`, `SECURITY.md`, `CODE_OF_CONDUCT.md`, issue templates, `CITATION.cff`).

## Status / current frontier

- Replication is complete: full pipeline filters the sample (n=935, 19 sites), prints OLS/PCR/EM-trace/tier-slopes/LOSO, runs a verification block against reported paper values, and saves privacy-safe CSVs plus four figures.

## House rules

- Privacy is non-negotiable: outputs under `results/` must contain no subject identifiers (`SUB_ID`, `FILE_ID`). The published verification block must continue to pass.
- The ABIDE-I phenotypic CSV is downloaded separately (see `data/README.md`) and never committed.
- File placement per global CLAUDE.md: analysis modules in `src/`, results in `results/`, figures in `figures/`.

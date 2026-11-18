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

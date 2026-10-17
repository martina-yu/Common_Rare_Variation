# ANALYSIS_RUN_INDEX

| Run ID | Actual completion | Definitions / evidence | Replacement |
|---|---|---|---|
| AGE-FOUR-20261016-01 | Real 7,996-case four-group descriptives completed; no regression or synthetic substitute | Carrier13-gene; shared within-case PRS median (>=high/<low); early<45 locked before summaries; source/command/software/input/output hashes/status in docs/age_interaction_runs/AGE-FOUR-20261016-01/provenance.json; report same directory | None; descriptive supplement |

Created at the author's explicit request on 2026-10-16. This index supplements the work index in ANALYSIS_CHANGELOG; historical execution records are not reconstructed or invented.

| Run ID | Parent | Scope / decisions | Execution status | Evidence | Supersedes |
|---|---|---|---|---|---|
| AGE-INT-20261016-01 | AGE-AUDIT-20261016-01; existing frozen C/Table 3 | Existing 7,996 cases; existing 13-gene B_raw_C and PRS311; within-case sample SD; continuous interaction; HC3 primary/HC1 comparison; intended for main article | completed on real frozen inputs; synthetic tests separately completed; live AoU access passed; no blocked/failed analysis stage | [report](docs/age_interaction_runs/AGE-INT-20261016-01/REPORT.md); [stage ledger](docs/age_interaction_runs/AGE-INT-20261016-01/stage_status.json); config/age_burden_prs_interaction_20261016.json; scripts/run_age_burden_prs_interaction.py; controlled full lock/snapshot/logs under data/analysis_runs/AGE-INT-20261016-01 | None; adds continuous interaction without replacing additive Table 3 or existing binary-carrier analysis; manuscript not edited |

Primary background PID4515, 16:05:25–16:06:09 UTC; diagnostic background PID5036, 16:19:43–16:19:45 UTC. All times 2026-10-16. Numeric verification used separate least-squares/sandwich code, real Table3 reproduction and unchanged input/source hashes. Source and exact aggregate output hashes are in provenance/report. Individual records and scores were not written; all private paths remain ignored and untracked. Raw genotype/PRS reconstruction, final variant adjudication, all13 external Table2 weight checks, incident analysis, nonlinear/leave-gene sensitivity models and manuscript insertion remain unperformed or unverified as explicitly listed in the report. No automatic model retry occurred; no result was selected by significance.

Scientific chronology: the author states this interaction was part of the original research question and was previously omitted. Its current implementation is locked before the new fit; no claim of preregistration or data-blind specification is made. Existing C variant eligibility remains provisional and is reused without changes. Participant records and person scores are never included in this index or its public evidence.


## AGE-QUART-20261016-01 / 2026-10-16

Real frozen-input Q1/Q4 descriptives completed (not synthetic-only), shared quartiles in original7,996 cases; selected3,998; early<45 fixed; existing13-gene carrier. Source/command/software/input and output hashes/status: docs/age_interaction_runs/AGE-QUART-20261016-01/provenance.json; [report](docs/age_interaction_runs/AGE-QUART-20261016-01/REPORT.md). No participant records/scores written; controlled logs and snapshot ignored. No regression or significance test; no failure/blocked stage. Supersedes=None; prior median split, additive and continuous-interaction results retained; manuscript unchanged.

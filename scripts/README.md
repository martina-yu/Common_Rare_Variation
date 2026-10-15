# Research analysis scripts

The numbered scripts contain the scientific analysis in sequential, numbered
sections. Shared modules contain reused statistics, annotation/genotype mechanics,
file/query utilities, disclosure policies and plotting primitives. There are no
lazy task objects in the main pipeline.

Run from the project root. `--list-tasks` and `--plan` inspect a stage without
reading participant data or requiring cloud credentials. `--task` selects named
sections; arguments after `--` are forwarded to a single section.

| Stage | Scientific work |
|---|---|
| 01_study_cohort.py | Female WGS cohort, ancestry labeling, cancer/EHR phenotypes |
| 02_sample_qc.py | Sample QC and clinical/technical covariates |
| 03_common_variants_prs.py | PGS000004 matching, ambiguity gates, PRS313/PRS311, PLINK2 comparison |
| 04_gene_selection.py | Frozen comparison membership, panels, weights, VAT candidates |
| 05_variant_annotation.py | PTV/LOFTEE audit and qualifying ClinVar missense |
| 06_variant_genotype_qc.py | Original FILTER, native ALT mapping, genotype QC and reconciliation |
| 07_burden_score.py | Gene-collapsed scores, comparisons and frozen model inputs |
| 08_case_control_models.py | Joint case-control models and paired 10x5 CV |
| 09_age_analysis.py | Within-case age OLS and final model/comparison outputs |
| 10_sensitivity_analysis.py | Prespecified sensitivity-score construction |
| 11_interaction_analysis.py | Explicitly post-hoc interaction robustness analyses |
| 12_result_tables.py | Cohort tables, attribution and result interpretation |
| 13_figures.py | Existing result figures and their visible data preparation |

## Execution dependencies

The numbering follows Methods topics; dependencies are not a single linear run.
Run `01`, `02 --task sample-qc`, `03`, then the covariate section of `02`.
Stages `04`–`07` retain the existing historical baseline prerequisites. A fresh
baseline uses explicit tasks in `legacy/old_9gene_ptv_analysis.py`.

Run `10 --task sensitivities` **before** `08`: its plan and scores remain locked
before association fitting. Stage `08` fits all existing comparison/sensitivity
versions and writes CV results plus aggregate case-model state. Stage `09` reads
that lossless JSON state, fits age models and writes the original combined result
tables. It checks that inputs and source hashes match stage `08`.

```bash
python scripts/03_common_variants_prs.py --plan
python scripts/06_variant_genotype_qc.py --task gene-calls -- --gene BRCA1
python scripts/10_sensitivity_analysis.py --task sensitivities
python scripts/08_case_control_models.py --task models
python scripts/09_age_analysis.py --task age-models
python scripts/11_interaction_analysis.py --task audit
python scripts/12_result_tables.py --task summary
python scripts/13_figures.py --task discrimination
```

Stage `11` requires explicit selection and remains post-hoc/exploratory.
Reviewer Cox, decomposition, family-history and frailty analyses remain in
`review/review_1008.py`; optional frailty entry: `exploratory/frailty_analysis.py`.
Manuscript assembly and refresh operations are explicit commands in `manuscript/`.

## Preserved behavior and output locations

Data remain under the existing `data/` paths. Results remain under
`results/table/`, `results/figures/` and `results/report/`. Output names retain
historical numbers. Existing model/variant status limitations are unchanged.
Saved burden scores retain ddof=0; model scaling retains ddof=1; CV scaling uses
training folds only. No cohort, eligibility, threshold, panel, weight, model,
covariate, test, seed or figure data is intentionally changed.

The pre-refactor numerical reference is preserved in
`archive/scripts_before_modularization_20261015/`; the interrupted working layout
is preserved in `archive/scripts_before_readability_20261015/`. Historical locks
are not rewritten. `readability_manifest.json` records the migration and current
implementation hashes.

**Known pre-existing blocker:** the original burden `score` task references an
undefined `report` when writing `38_rowlands_section7_summary.md`. The unrelated
replacement text from the interrupted pass was reverted; no silent report fix
is included. Existing scientific score outputs and the reported limitation are
preserved. A fresh default burden run needs that separate reporting repair.

## Validation and traceability

Tests are in project-level `tests/`. Run `python -m unittest discover -s tests`.
Plain synthetic test functions are also collected by `tests/run_checks.py`.
The model-stage test executes the preserved original versus sequential 08/09 on
synthetic files and compares 13 CSV files byte-for-byte and all ten CV prediction
arrays exactly. This is not a real-data pipeline rerun.

See [function index](FUNCTION_INDEX.md), [result registry](../RESULT_REGISTRY.md),
[project status](../PROJECT_STATUS.md), [change log](../ANALYSIS_CHANGELOG.md), and
[validation record](../docs/readability_validation.md).

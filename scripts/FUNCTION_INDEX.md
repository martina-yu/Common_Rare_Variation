# Analysis function and section index

Scientific operations are sequential sections in the main files.

| Script | Task | Section line |
|---|---|---|
| [scripts/01_study_cohort.py](../scripts/01_study_cohort.py:71) | `select` | 71 |
| [scripts/01_study_cohort.py](../scripts/01_study_cohort.py:218) | `phenotypes` | 218 |
| [scripts/02_sample_qc.py](../scripts/02_sample_qc.py:628) | `sample-qc` | 628 |
| [scripts/02_sample_qc.py](../scripts/02_sample_qc.py:687) | `covariates` | 687 |
| [scripts/03_common_variants_prs.py](../scripts/03_common_variants_prs.py:573) | `prs` | 573 |
| [scripts/03_common_variants_prs.py](../scripts/03_common_variants_prs.py:968) | `match-audit` | 968 |
| [scripts/04_gene_selection.py](../scripts/04_gene_selection.py:122) | `panel` | 122 |
| [scripts/04_gene_selection.py](../scripts/04_gene_selection.py:305) | `candidates` | 305 |
| [scripts/05_variant_annotation.py](../scripts/05_variant_annotation.py:96) | `ptv` | 96 |
| [scripts/05_variant_annotation.py](../scripts/05_variant_annotation.py:405) | `missense` | 405 |
| [scripts/06_variant_genotype_qc.py](../scripts/06_variant_genotype_qc.py:185) | `site-qc` | 185 |
| [scripts/06_variant_genotype_qc.py](../scripts/06_variant_genotype_qc.py:284) | `calls` | 284 |
| [scripts/06_variant_genotype_qc.py](../scripts/06_variant_genotype_qc.py:307) | `gene-calls` | 307 |
| [scripts/06_variant_genotype_qc.py](../scripts/06_variant_genotype_qc.py:840) | `missing-sites` | 840 |
| [scripts/06_variant_genotype_qc.py](../scripts/06_variant_genotype_qc.py:1109) | `validate-qc` | 1109 |
| [scripts/06_variant_genotype_qc.py](../scripts/06_variant_genotype_qc.py:1454) | `call-audit` | 1454 |
| [scripts/06_variant_genotype_qc.py](../scripts/06_variant_genotype_qc.py:1952) | `reconcile` | 1952 |
| [scripts/07_burden_score.py](../scripts/07_burden_score.py:77) | `score` | 77 |
| [scripts/07_burden_score.py](../scripts/07_burden_score.py:255) | `validate-score` | 255 |
| [scripts/07_burden_score.py](../scripts/07_burden_score.py:462) | `compare` | 462 |
| [scripts/07_burden_score.py](../scripts/07_burden_score.py:902) | `model-inputs` | 902 |
| [scripts/08_case_control_models.py](../scripts/08_case_control_models.py:106) | `models` | 106 |
| [scripts/09_age_analysis.py](../scripts/09_age_analysis.py:94) | `age-models` | 94 |
| [scripts/10_sensitivity_analysis.py](../scripts/10_sensitivity_analysis.py:59) | `sensitivities` | 59 |
| [scripts/11_interaction_analysis.py](../scripts/11_interaction_analysis.py:94) | `audit` | 94 |
| [scripts/12_result_tables.py](../scripts/12_result_tables.py:137) | `descriptives` | 137 |
| [scripts/12_result_tables.py](../scripts/12_result_tables.py:818) | `cohort-table` | 818 |
| [scripts/12_result_tables.py](../scripts/12_result_tables.py:1104) | `attribution` | 1104 |
| [scripts/12_result_tables.py](../scripts/12_result_tables.py:1379) | `summary` | 1379 |
| [scripts/13_figures.py](../scripts/13_figures.py:244) | `age-figures` | 244 |
| [scripts/13_figures.py](../scripts/13_figures.py:784) | `score-distributions` | 784 |
| [scripts/13_figures.py](../scripts/13_figures.py:934) | `discrimination` | 934 |
| [scripts/13_figures.py](../scripts/13_figures.py:1331) | `joint-distribution` | 1331 |
| [scripts/13_figures.py](../scripts/13_figures.py:1632) | `score-deciles` | 1632 |
| [scripts/13_figures.py](../scripts/13_figures.py:1805) | `burden-histogram` | 1805 |
| [scripts/13_figures.py](../scripts/13_figures.py:2043) | `age-boxplots` | 2043 |
| [scripts/13_figures.py](../scripts/13_figures.py:2160) | `age-carrier-frequency` | 2160 |
| [scripts/13_figures.py](../scripts/13_figures.py:2323) | `sensitivity-figure` | 2323 |

Shared and local helper definitions:

| Source | Function | Line | Purpose |
|---|---|---|---|
| [scripts/01_study_cohort.py](../scripts/01_study_cohort.py:39) | `cohort_selection_section` | 39 | See the function and its analysis section. |
| [scripts/01_study_cohort.py](../scripts/01_study_cohort.py:43) | `phenotypes_section` | 43 | See the function and its analysis section. |
| [scripts/01_study_cohort.py](../scripts/01_study_cohort.py:52) | `_group` | 52 | See the function and its analysis section. |
| [scripts/01_study_cohort.py](../scripts/01_study_cohort.py:61) | `_age_years` | 61 | See the function and its analysis section. |
| [scripts/02_sample_qc.py](../scripts/02_sample_qc.py:53) | `run` | 53 | See the function and its analysis section. |
| [scripts/02_sample_qc.py](../scripts/02_sample_qc.py:65) | `sample_qc_section` | 65 | See the function and its analysis section. |
| [scripts/02_sample_qc.py](../scripts/02_sample_qc.py:71) | `kinship_pruning` | 71 | See the function and its analysis section. |
| [scripts/02_sample_qc.py](../scripts/02_sample_qc.py:138) | `build_representative_acaf_slice` | 138 | See the function and its analysis section. |
| [scripts/02_sample_qc.py](../scripts/02_sample_qc.py:216) | `sample_missingness` | 216 | See the function and its analysis section. |
| [scripts/02_sample_qc.py](../scripts/02_sample_qc.py:288) | `hwe_control_only` | 288 | See the function and its analysis section. |
| [scripts/02_sample_qc.py](../scripts/02_sample_qc.py:347) | `load_candidates` | 347 | See the function and its analysis section. |
| [scripts/02_sample_qc.py](../scripts/02_sample_qc.py:374) | `count_row` | 374 | See the function and its analysis section. |
| [scripts/02_sample_qc.py](../scripts/02_sample_qc.py:379) | `build_table_s2` | 379 | Build table s2. |
| [scripts/02_sample_qc.py](../scripts/02_sample_qc.py:491) | `clinical_covariates_log` | 491 | See the function and its analysis section. |
| [scripts/02_sample_qc.py](../scripts/02_sample_qc.py:495) | `clinical_covariates_section` | 495 | See the function and its analysis section. |
| [scripts/02_sample_qc.py](../scripts/02_sample_qc.py:501) | `map_menopause` | 501 | See the function and its analysis section. |
| [scripts/02_sample_qc.py](../scripts/02_sample_qc.py:517) | `relative_label` | 517 | See the function and its analysis section. |
| [scripts/02_sample_qc.py](../scripts/02_sample_qc.py:528) | `_famhist_label` | 528 | See the function and its analysis section. |
| [scripts/02_sample_qc.py](../scripts/02_sample_qc.py:536) | `derive_smoking` | 536 | See the function and its analysis section. |
| [scripts/02_sample_qc.py](../scripts/02_sample_qc.py:554) | `summarize_continuous` | 554 | See the function and its analysis section. |
| [scripts/02_sample_qc.py](../scripts/02_sample_qc.py:575) | `summarize_categorical` | 575 | See the function and its analysis section. |
| [scripts/03_common_variants_prs.py](../scripts/03_common_variants_prs.py:46) | `run` | 46 | See the function and its analysis section. |
| [scripts/03_common_variants_prs.py](../scripts/03_common_variants_prs.py:58) | `section` | 58 | See the function and its analysis section. |
| [scripts/03_common_variants_prs.py](../scripts/03_common_variants_prs.py:64) | `download_scoring_files` | 64 | See the function and its analysis section. |
| [scripts/03_common_variants_prs.py](../scripts/03_common_variants_prs.py:80) | `load_cohort` | 80 | See the function and its analysis section. |
| [scripts/03_common_variants_prs.py](../scripts/03_common_variants_prs.py:96) | `build_target_table` | 96 | See the function and its analysis section. |
| [scripts/03_common_variants_prs.py](../scripts/03_common_variants_prs.py:135) | `match_variants_to_aou` | 135 | See the function and its analysis section. |
| [scripts/03_common_variants_prs.py](../scripts/03_common_variants_prs.py:187) | `resolve_matches` | 187 | Resolve matches. |
| [scripts/03_common_variants_prs.py](../scripts/03_common_variants_prs.py:237) | `load_fam_and_cohort_index` | 237 | See the function and its analysis section. |
| [scripts/03_common_variants_prs.py](../scripts/03_common_variants_prs.py:266) | `fetch_and_decode_all` | 266 | Fetch and decode all. |
| [scripts/03_common_variants_prs.py](../scripts/03_common_variants_prs.py:307) | `resolve_ambiguous_by_af` | 307 | Resolve ambiguous by af. |
| [scripts/03_common_variants_prs.py](../scripts/03_common_variants_prs.py:503) | `plink2_cross_check` | 503 | See the function and its analysis section. |
| [scripts/04_gene_selection.py](../scripts/04_gene_selection.py:85) | `read` | 85 | See the function and its analysis section. |
| [scripts/04_gene_selection.py](../scripts/04_gene_selection.py:89) | `ids` | 89 | See the function and its analysis section. |
| [scripts/04_gene_selection.py](../scripts/04_gene_selection.py:98) | `validate_extension_membership` | 98 | Call before later extension models; no implicit dropping of participants. |
| [scripts/06_variant_genotype_qc.py](../scripts/06_variant_genotype_qc.py:149) | `run` | 149 | See the function and its analysis section. |
| [scripts/08_case_control_models.py](../scripts/08_case_control_models.py:85) | `design` | 85 | See the function and its analysis section. |
| [scripts/08_case_control_models.py](../scripts/08_case_control_models.py:94) | `fold_job` | 94 | See the function and its analysis section. |
| [scripts/09_age_analysis.py](../scripts/09_age_analysis.py:85) | `design` | 85 | See the function and its analysis section. |
| [scripts/11_interaction_analysis.py](../scripts/11_interaction_analysis.py:40) | `fit` | 40 | See the function and its analysis section. |
| [scripts/11_interaction_analysis.py](../scripts/11_interaction_analysis.py:68) | `wald` | 68 | See the function and its analysis section. |
| [scripts/11_interaction_analysis.py](../scripts/11_interaction_analysis.py:87) | `spline` | 87 | See the function and its analysis section. |
| [scripts/12_result_tables.py](../scripts/12_result_tables.py:61) | `cohort_descriptives_log` | 61 | See the function and its analysis section. |
| [scripts/12_result_tables.py](../scripts/12_result_tables.py:65) | `summarize` | 65 | See the function and its analysis section. |
| [scripts/12_result_tables.py](../scripts/12_result_tables.py:116) | `fmt_cat` | 116 | See the function and its analysis section. |
| [scripts/12_result_tables.py](../scripts/12_result_tables.py:124) | `fmt_cont` | 124 | See the function and its analysis section. |
| [scripts/13_figures.py](../scripts/13_figures.py:196) | `save_score_distribution_figure` | 196 | See the function and its analysis section. |
| [scripts/13_figures.py](../scripts/13_figures.py:214) | `discrimination_design_matrix` | 214 | See the function and its analysis section. |
| [scripts/13_figures.py](../scripts/13_figures.py:221) | `fit_discrimination_model` | 221 | See the function and its analysis section. |
| [scripts/13_figures.py](../scripts/13_figures.py:227) | `save_discrimination_figure` | 227 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:36) | `age_group_contrasts_read` | 36 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:40) | `contrast` | 40 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:52) | `design` | 52 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:56) | `analyze` | 56 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:148) | `analyze_age_group_contrasts` | 148 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:231) | `ehr_age_audit_read` | 231 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:235) | `releasable` | 235 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:240) | `save` | 240 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:244) | `table_counts` | 244 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:252) | `published_joint_margin_blocks` | 252 | A released adjacent joint cell plus a new margin reveals a suppressed cell. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:269) | `binary_table` | 269 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:298) | `log_odds_contrast` | 298 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:313) | `case_only` | 313 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:392) | `ehr_query` | 392 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:460) | `query_windows` | 460 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:494) | `numeric_summary` | 494 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:502) | `span_proxy_masks` | 502 | Feasibility only: entry after washout; censor at clinical end/death/45. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:517) | `window_audit` | 517 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:681) | `audit_early_ehr_diagnoses` | 681 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:798) | `nine_gene_missense_age_read` | 798 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:802) | `safe_partition` | 802 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:808) | `related_condition` | 808 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:819) | `disease_allowed` | 819 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:827) | `germline_allowed` | 827 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:838) | `prepare_candidates` | 838 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:989) | `ensure_ro_mount` | 989 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:1019) | `qualifying_call` | 1019 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:1033) | `minimal_allele_key` | 1033 | Minimal allele key. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:1045) | `candidate_shards` | 1045 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:1060) | `extract_carriers` | 1060 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:1193) | `carrier_partitions` | 1193 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:1261) | `linear_contrast` | 1261 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:1273) | `fit_age` | 1273 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:1345) | `finalize_disclosure_and_report` | 1345 | Check released aggregates against earlier ten-gene totals before reporting. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:1465) | `analyze_nine_gene_missense_age` | 1465 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:1586) | `prs_age_descriptives_log` | 1586 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:1590) | `describe_prs_and_age` | 1590 | Initialize runtime state and execute this component explicitly. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:1984) | `legacy_burden_descriptives_log` | 1984 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:1988) | `legacy_burden_descriptives_section` | 1988 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:1994) | `legacy_burden_descriptives_suppress` | 1994 | Legacy burden descriptives suppress. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:2001) | `legacy_burden_descriptives_wilson_ci` | 2001 | Wilson score 95% CI for a binomial proportion. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:2012) | `gene_breakdown` | 2012 | Gene breakdown. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:2041) | `pool_set` | 2041 | Pool set. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:2053) | `merge_categories` | 2053 | Merge categories. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:2080) | `fmt_p` | 2080 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:2084) | `describe_legacy_burden` | 2084 | Initialize runtime state and execute this component explicitly. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:2900) | `legacy_age_contrasts_log` | 2900 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:2904) | `legacy_age_contrasts_section` | 2904 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:2910) | `legacy_age_contrasts_suppress` | 2910 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:2916) | `bucket` | 2916 | Bucket. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:2926) | `sanitize_for_export` | 2926 | Sanitize for export. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:2936) | `legacy_age_contrasts_wilson_ci` | 2936 | Legacy age contrasts wilson ci. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:2947) | `zscore` | 2947 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:2952) | `smd` | 2952 | Smd. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:2961) | `fit_ols_hc1` | 2961 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:2965) | `fit_logit_hc1` | 2965 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:2969) | `coef_rows` | 2969 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:2992) | `build_X` | 2992 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:2997) | `chi2_2x2` | 2997 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:3003) | `mean_p` | 3003 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:3007) | `mix_cat` | 3007 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:3017) | `make_env_X` | 3017 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:3027) | `_lab` | 3027 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:3031) | `fmt` | 3031 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:3035) | `pfmt` | 3035 | See the function and its analysis section. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:3039) | `_pct_cell` | 3039 | pct cell. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:3046) | `analyze_legacy_age_contrasts` | 3046 | Initialize runtime state and execute this component explicitly. |
| [scripts/exploratory/age_ehr_checks.py](../scripts/exploratory/age_ehr_checks.py:4634) | `main` | 4634 | See the function and its analysis section. |
| [scripts/functions/age_group_contrasts.py](../scripts/functions/age_group_contrasts.py:18) | `read` | 18 | See the function and its analysis section. |
| [scripts/functions/age_group_contrasts.py](../scripts/functions/age_group_contrasts.py:22) | `contrast` | 22 | See the function and its analysis section. |
| [scripts/functions/age_group_contrasts.py](../scripts/functions/age_group_contrasts.py:34) | `design` | 34 | See the function and its analysis section. |
| [scripts/functions/age_group_contrasts.py](../scripts/functions/age_group_contrasts.py:38) | `analyze` | 38 | See the function and its analysis section. |
| [scripts/functions/bigquery.py](../scripts/functions/bigquery.py:6) | `timed_query` | 6 | See the function and its analysis section. |
| [scripts/functions/burden.py](../scripts/functions/burden.py:16) | `read_weights` | 16 | See the function and its analysis section. |
| [scripts/functions/burden.py](../scripts/functions/burden.py:26) | `score_panel` | 26 | See the function and its analysis section. |
| [scripts/functions/cli.py](../scripts/functions/cli.py:12) | `component_arguments` | 12 | See the function and its analysis section. |
| [scripts/functions/cli.py](../scripts/functions/cli.py:21) | `run_local_step` | 21 | See the function and its analysis section. |
| [scripts/functions/cli.py](../scripts/functions/cli.py:56) | `parse_tasks` | 56 | Select explicit sequential sections; metadata commands never run analyses. |
| [scripts/functions/clinvar.py](../scripts/functions/clinvar.py:29) | `write_tsv_if_changed` | 29 | Keep valid genotype checkpoints when a reference revalidation is identical. |
| [scripts/functions/clinvar.py](../scripts/functions/clinvar.py:39) | `clean` | 39 | See the function and its analysis section. |
| [scripts/functions/clinvar.py](../scripts/functions/clinvar.py:43) | `exclusions` | 43 | See the function and its analysis section. |
| [scripts/functions/clinvar.py](../scripts/functions/clinvar.py:55) | `fetch` | 55 | See the function and its analysis section. |
| [scripts/functions/clinvar.py](../scripts/functions/clinvar.py:79) | `historical_xml` | 79 | See the function and its analysis section. |
| [scripts/functions/clinvar.py](../scripts/functions/clinvar.py:112) | `disease_rule` | 112 | Disease rule. |
| [scripts/functions/clinvar.py](../scripts/functions/clinvar.py:196) | `review_variant` | 196 | See the function and its analysis section. |
| [scripts/functions/cohort.py](../scripts/functions/cohort.py:25) | `disclose_group_counts` | 25 | Print a category breakdown under the n>=20 rule. |
| [scripts/functions/cohort.py](../scripts/functions/cohort.py:48) | `extract_ancestry_preds` | 48 | Stream the release-wide ancestry table, keep our cohort, verify meaning. |
| [scripts/functions/cohort.py](../scripts/functions/cohort.py:133) | `gs_cat_lines` | 133 | Stream a GCS object line by line: read-only, requester-pays, never |
| [scripts/functions/cohort.py](../scripts/functions/cohort.py:151) | `log` | 151 | Progress goes to stderr, because stdout carries the workbook bytes. |
| [scripts/functions/cohort.py](../scripts/functions/cohort.py:156) | `parse_hail_scores` | 156 | training_pca.tsv keeps its PC vector as a bracketed string of Java |
| [scripts/functions/disclosure.py](../scripts/functions/disclosure.py:38) | `decide` | 38 | Return {level -> set of masked column keys}. |
| [scripts/functions/disclosure.py](../scripts/functions/disclosure.py:90) | `render` | 90 | One published cell: exact count with percentage, a bin, or a withholding. |
| [scripts/functions/disclosure.py](../scripts/functions/disclosure.py:103) | `format_count` | 103 | Preserve the existing cohort/QC count displays, including the n=20 boundary. |
| [scripts/functions/ehr_age_audit.py](../scripts/functions/ehr_age_audit.py:20) | `read` | 20 | See the function and its analysis section. |
| [scripts/functions/ehr_age_audit.py](../scripts/functions/ehr_age_audit.py:24) | `releasable` | 24 | See the function and its analysis section. |
| [scripts/functions/ehr_age_audit.py](../scripts/functions/ehr_age_audit.py:29) | `save` | 29 | See the function and its analysis section. |
| [scripts/functions/ehr_age_audit.py](../scripts/functions/ehr_age_audit.py:33) | `table_counts` | 33 | See the function and its analysis section. |
| [scripts/functions/ehr_age_audit.py](../scripts/functions/ehr_age_audit.py:41) | `published_joint_margin_blocks` | 41 | A released adjacent joint cell plus a new margin reveals a suppressed cell. |
| [scripts/functions/ehr_age_audit.py](../scripts/functions/ehr_age_audit.py:58) | `binary_table` | 58 | See the function and its analysis section. |
| [scripts/functions/ehr_age_audit.py](../scripts/functions/ehr_age_audit.py:87) | `log_odds_contrast` | 87 | See the function and its analysis section. |
| [scripts/functions/ehr_age_audit.py](../scripts/functions/ehr_age_audit.py:102) | `case_only` | 102 | See the function and its analysis section. |
| [scripts/functions/ehr_age_audit.py](../scripts/functions/ehr_age_audit.py:181) | `ehr_query` | 181 | See the function and its analysis section. |
| [scripts/functions/ehr_age_audit.py](../scripts/functions/ehr_age_audit.py:249) | `query_windows` | 249 | See the function and its analysis section. |
| [scripts/functions/ehr_age_audit.py](../scripts/functions/ehr_age_audit.py:283) | `numeric_summary` | 283 | See the function and its analysis section. |
| [scripts/functions/ehr_age_audit.py](../scripts/functions/ehr_age_audit.py:291) | `span_proxy_masks` | 291 | Feasibility only: entry after washout; censor at clinical end/death/45. |
| [scripts/functions/ehr_age_audit.py](../scripts/functions/ehr_age_audit.py:306) | `window_audit` | 306 | See the function and its analysis section. |
| [scripts/functions/frailty.py](../scripts/functions/frailty.py:12) | `exposure_matrix` | 12 | See the function and its analysis section. |
| [scripts/functions/frailty.py](../scripts/functions/frailty.py:18) | `gamma_loglik_gradient` | 18 | See the function and its analysis section. |
| [scripts/functions/frailty.py](../scripts/functions/frailty.py:44) | `fit_gamma` | 44 | See the function and its analysis section. |
| [scripts/functions/frailty.py](../scripts/functions/frailty.py:68) | `profile_theta` | 68 | Profile theta; boundary-aware interval reported as a modelling sensitivity. |
| [scripts/functions/genotype_qc.py](../scripts/functions/genotype_qc.py:23) | `sha` | 23 | See the function and its analysis section. |
| [scripts/functions/genotype_qc.py](../scripts/functions/genotype_qc.py:31) | `call_qc` | 31 | See the function and its analysis section. |
| [scripts/functions/genotype_qc.py](../scripts/functions/genotype_qc.py:68) | `normalized_native_sites` | 68 | Normalize single ALT site copies; leave original GT/AD indices untouched. |
| [scripts/functions/genotype_qc.py](../scripts/functions/genotype_qc.py:102) | `safe_count` | 102 | See the function and its analysis section. |
| [scripts/functions/io.py](../scripts/functions/io.py:7) | `sha256` | 7 | See the function and its analysis section. |
| [scripts/functions/nine_gene_missense_age.py](../scripts/functions/nine_gene_missense_age.py:40) | `read` | 40 | See the function and its analysis section. |
| [scripts/functions/nine_gene_missense_age.py](../scripts/functions/nine_gene_missense_age.py:44) | `safe_partition` | 44 | See the function and its analysis section. |
| [scripts/functions/nine_gene_missense_age.py](../scripts/functions/nine_gene_missense_age.py:48) | `related_condition` | 48 | See the function and its analysis section. |
| [scripts/functions/nine_gene_missense_age.py](../scripts/functions/nine_gene_missense_age.py:59) | `disease_allowed` | 59 | See the function and its analysis section. |
| [scripts/functions/nine_gene_missense_age.py](../scripts/functions/nine_gene_missense_age.py:67) | `germline_allowed` | 67 | See the function and its analysis section. |
| [scripts/functions/nine_gene_missense_age.py](../scripts/functions/nine_gene_missense_age.py:78) | `prepare_candidates` | 78 | See the function and its analysis section. |
| [scripts/functions/nine_gene_missense_age.py](../scripts/functions/nine_gene_missense_age.py:227) | `ensure_ro_mount` | 227 | See the function and its analysis section. |
| [scripts/functions/nine_gene_missense_age.py](../scripts/functions/nine_gene_missense_age.py:257) | `qualifying_call` | 257 | See the function and its analysis section. |
| [scripts/functions/nine_gene_missense_age.py](../scripts/functions/nine_gene_missense_age.py:271) | `minimal_allele_key` | 271 | Minimal allele key. |
| [scripts/functions/nine_gene_missense_age.py](../scripts/functions/nine_gene_missense_age.py:283) | `candidate_shards` | 283 | See the function and its analysis section. |
| [scripts/functions/nine_gene_missense_age.py](../scripts/functions/nine_gene_missense_age.py:298) | `extract_carriers` | 298 | See the function and its analysis section. |
| [scripts/functions/nine_gene_missense_age.py](../scripts/functions/nine_gene_missense_age.py:431) | `carrier_partitions` | 431 | See the function and its analysis section. |
| [scripts/functions/nine_gene_missense_age.py](../scripts/functions/nine_gene_missense_age.py:497) | `linear_contrast` | 497 | See the function and its analysis section. |
| [scripts/functions/nine_gene_missense_age.py](../scripts/functions/nine_gene_missense_age.py:509) | `fit_age` | 509 | See the function and its analysis section. |
| [scripts/functions/nine_gene_missense_age.py](../scripts/functions/nine_gene_missense_age.py:577) | `finalize_disclosure_and_report` | 577 | Check released aggregates against earlier ten-gene totals before reporting. |
| [scripts/functions/paths.py](../scripts/functions/paths.py:8) | `implementation_sources` | 8 | See the function and its analysis section. |
| [scripts/functions/plotting.py](../scripts/functions/plotting.py:30) | `grouped_review_strata` | 30 | Compact three-section forest plot with fixed subgroup colors/shapes. |
| [scripts/functions/plotting.py](../scripts/functions/plotting.py:162) | `review_group` | 162 | See the function and its analysis section. |
| [scripts/functions/plotting.py](../scripts/functions/plotting.py:193) | `grouped_review_forest` | 193 | Publication forest with aligned labels, estimates, and numerical CIs. |
| [scripts/functions/plotting.py](../scripts/functions/plotting.py:303) | `compressed_kde` | 303 | Exact full-sample Gaussian KDE via duplicate-value compression. |
| [scripts/functions/provenance.py](../scripts/functions/provenance.py:10) | `sha256` | 10 | See the function and its analysis section. |
| [scripts/functions/provenance.py](../scripts/functions/provenance.py:14) | `_verify_historical_layout` | 14 | See the function and its analysis section. |
| [scripts/functions/provenance.py](../scripts/functions/provenance.py:40) | `verify_frozen_inputs` | 40 | Verify exact data hashes and explicitly recorded source migrations. |
| [scripts/functions/prs.py](../scripts/functions/prs.py:5) | `_orient_effect_allele` | 5 | orient effect allele. |
| [scripts/functions/prs.py](../scripts/functions/prs.py:14) | `_decode_bed_block` | 14 | See the function and its analysis section. |
| [scripts/functions/prs.py](../scripts/functions/prs.py:25) | `_encode_bed_bytes` | 25 | encode bed bytes. |
| [scripts/functions/sample_variant_audit.py](../scripts/functions/sample_variant_audit.py:58) | `load_module` | 58 | See the function and its analysis section. |
| [scripts/functions/sample_variant_audit.py](../scripts/functions/sample_variant_audit.py:65) | `positive_float` | 65 | See the function and its analysis section. |
| [scripts/functions/site_qc.py](../scripts/functions/site_qc.py:16) | `source_site_pass` | 16 | Empty defined Hail filter set is PASS; missing source information fails. |
| [scripts/functions/site_qc.py](../scripts/functions/site_qc.py:21) | `canonical_crosswalk` | 21 | See the function and its analysis section. |
| [scripts/functions/stats.py](../scripts/functions/stats.py:6) | `delong_auc_variance` | 6 | Delong auc variance. |
| [scripts/functions/stats.py](../scripts/functions/stats.py:31) | `delong_test` | 31 | Delong test. |
| [scripts/functions/stats.py](../scripts/functions/stats.py:60) | `auc_ci` | 60 | See the function and its analysis section. |
| [scripts/functions/stats.py](../scripts/functions/stats.py:69) | `fit_binomial` | 69 | See the function and its analysis section. |
| [scripts/functions/variants.py](../scripts/functions/variants.py:44) | `terms` | 44 | See the function and its analysis section. |
| [scripts/functions/variants.py](../scripts/functions/variants.py:48) | `key` | 48 | See the function and its analysis section. |
| [scripts/functions/variants.py](../scripts/functions/variants.py:52) | `coordinates` | 52 | See the function and its analysis section. |
| [scripts/functions/variants.py](../scripts/functions/variants.py:86) | `normalize` | 86 | See the function and its analysis section. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:158) | `log` | 158 | Progress goes to stderr, because stdout carries the workbook bytes. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:163) | `section` | 163 | See the function and its analysis section. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:172) | `bq` | 172 | Run one query against the CDR and return a DataFrame. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:192) | `cached` | 192 | Read `name` from the cache, or build it once and store it. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:209) | `survey_latest` | 209 | Most recent answer to one survey question, for every patient in `bc`. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:227) | `load_patients` | 227 | Every person with a breast cancer diagnosis, with the three fields that |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:246) | `load_wgs_ids` | 246 | The person_ids with whole-genome sequencing, from the CRAM manifest. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:258) | `add_wgs` | 258 | Sequencing availability as a COLUMN, not a filter: this table's |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:266) | `add_age` | 266 | Age at FIRST breast cancer diagnosis = that date minus date of birth. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:285) | `add_menopause` | 285 | Menopausal status, from the latest answer to the menopause survey |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:309) | `add_smoking` | 309 | Smoking, from two survey questions: "100 cigarettes in your life?" gives |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:343) | `add_famhist` | 343 | Family history of breast cancer, from the multi-select relatives |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:370) | `add_bmi` | 370 | BMI: the most recent recorded measurement per person, keeping only |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:395) | `add_receptors` | 395 | ER / PR / HER2 status from diagnosis codes. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:424) | `add_stage` | 424 | Whether ANY staging-like code is recorded -- not a stage distribution. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:461) | `collapse_categories` | 461 | Merge the categories that would otherwise sit below the floor. Applied |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:475) | `split_columns` | 475 | The table's three columns: early onset, later onset, everyone. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:480) | `count_blocks` | 480 | counts[(block, level)][column], with exhaustiveness enforced. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:504) | `read_published` | 504 | Everything an earlier, already-public run states EXACTLY. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:562) | `decide_masking` | 562 | Which cells cannot be stated exactly. Returns the masked set. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:687) | `render_cells` | 687 | Turn every cell into the text the table prints. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:739) | `build_table` | 739 | Assemble the finished table: header rows, then one row per block level. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:797) | `check_against_published` | 797 | Compare the finished table with the version already in the repository. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:835) | `fill_workbook` | 835 | Fill the All of Us columns of the template and return the .xlsx bytes. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:935) | `disclose_n` | 935 | A single count with no exact partner to protect (own-floor rule only): |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:945) | `load_published_anchors` | 945 | The handful of cells in the already-published clinical table that are |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:988) | `check_diff` | 988 | One count computed fresh here, and one already public: their difference |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:1010) | `compute_age_soft` | 1010 | Age at first BC diagnosis, same formula as add_age, but WITHOUT the |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:1021) | `qc_step1_female` | 1021 | Step 1: keep sex_at_birth == 'Female' only. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:1044) | `qc_step2_age` | 1044 | Step 2: keep samples with a computable age at diagnosis. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:1066) | `qc_wgs_join` | 1066 | The structurally-required step between steps 1-2 and genotype QC: |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:1091) | `write_keep_list` | 1091 | A plink2 --keep file: one #IID per line, matching the psam's ID scheme |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:1105) | `pgen_gcs_paths` | 1105 | See the function and its analysis section. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:1110) | `gcs_size` | 1110 | Byte size of one GCS object, via `gsutil -u PROJECT ls -l` (read-only). |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:1119) | `download_pgen` | 1119 | Download one chromosome's pgen/pvar/psam triple to a local scratch dir. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:1146) | `cleanup_pgen` | 1146 | Delete one chromosome's local pgen/pvar/psam (and any plink2 outputs |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:1163) | `plink_missing` | 1163 | Run `plink2 --pfile ... --keep ... --missing --out ...` and return the |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:1231) | `load_cached_smiss` | 1231 | A previously-cached .smiss for one chromosome, if present and it |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:1254) | `run_chrom_missing` | 1254 | Full step-3 pipeline for ONE chromosome, resumable: |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:1338) | `exome_wide_fmiss` | 1338 | Exome-wide F_MISS per sample: sum MISSING_CT across all 22 |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:1382) | `plot_fmiss_distribution` | 1382 | Histogram of exome-wide per-sample F_MISS. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:1513) | `relationship_class` | 1513 | Map one kinship coefficient onto KING's relationship bands. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:1527) | `load_flagged_samples` | 1527 | The person_ids All of Us flags as population outliers. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:1587) | `load_aou_relatedness` | 1587 | The AoU pairwise kinship table, cached locally (person-level -> data/). |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:1609) | `split_relatedness` | 1609 | Partition the AoU pairs against our cohort. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:1629) | `build_adjacency` | 1629 | {sample: {relative: kin}} over the within-cohort pairs. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:1638) | `greedy_unrelated` | 1638 | Greedy maximal independent set: which samples to DROP. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:1680) | `relationship_tables` | 1680 | The four internal tables. All of them live under data/ and are never |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:1762) | `plot_kinship_density` | 1762 | Kernel-density curve of kinship, with KING's band edges marked. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:1804) | `keep_smallest_format` | 1804 | Project rule: measure PNG/JPG/SVG, keep the smallest, delete the rest. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:1863) | `git_check_ignored` | 1863 | Assert that `path` is covered by .gitignore, and log the matching rule. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:1884) | `read_keep_ids` | 1884 | Read a keep list (plink2 '#IID' format or a one-column csv) as strings. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:1893) | `fetch_ancestry_reference` | 1893 | Download the reference panel and its eigenvalues (inbound public data: |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:1910) | `disclose_group_counts` | 1910 | Print a category breakdown under the n>=20 rule. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:1933) | `demotion_breakdown` | 1933 | For the samples argmax-called EUR that the confidence threshold sends to |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:1951) | `_anc_colors` | 1951 | See the function and its analysis section. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:1968) | `anc_legend` | 1968 | Legend swatches at full opacity: the scatter alpha is for the cloud, not |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:1977) | `anc_scatter` | 1977 | Draw groups largest-first so small ones stay visible on top. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:2014) | `save_fig_smallest` | 2014 | Write svg/png/jpg, measure all three, keep the smallest, delete the rest. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:2043) | `plot_reference_vs_cohort` | 2043 | The figure that makes the axes interpretable: the labelled HGDP+1000G |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:2066) | `plot_cohort_pc_pairs` | 2066 | PC1-PC2 and PC3-PC4 for the cohort: the standard supplementary panel. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:2083) | `plot_two_label_rules` | 2083 | The same points twice: identical coordinates, identical classifier, only |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:2126) | `plot_reference_scree` | 2126 | Variance explained by each reference PC. Bar heights are eigenvalue |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:2143) | `plot_confidence_density` | 2143 | Kernel densities, not histograms: a bar height would be a group count. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:2173) | `plot_pc_variance_within_eur` | 2173 | How much of each projected axis survives the restriction to EUR. |
| [scripts/legacy/clinical_qc_helpers.py](../scripts/legacy/clinical_qc_helpers.py:2204) | `plot_within_eur_axes` | 2204 | The surviving axes, EUR subset only, coloured by the confidence rule so |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:26) | `vat_schema_run` | 26 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:36) | `vat_schema_section` | 36 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:42) | `vat_schema_calculation` | 42 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:152) | `analyze_vat_schema` | 152 | Initialize runtime state and execute this component explicitly. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:179) | `legacy_vat_preparation_run` | 179 | Legacy vat preparation run. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:200) | `legacy_vat_preparation_section` | 200 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:206) | `ensembl_get` | 206 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:217) | `legacy_vat_preparation_calculation` | 217 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:596) | `analyze_legacy_vat_preparation` | 596 | Initialize runtime state and execute this component explicitly. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:683) | `legacy_vat_repair_run` | 683 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:698) | `legacy_vat_repair_calculation` | 698 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:781) | `analyze_legacy_vat_repair` | 781 | Initialize runtime state and execute this component explicitly. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:850) | `legacy_variant_selection_section` | 850 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:856) | `write_variant_list` | 856 | Write variant list. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:864) | `write_public_variant_table` | 864 | Write public variant table. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:888) | `to_float` | 888 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:897) | `is_ptv` | 897 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:902) | `is_missense` | 902 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:907) | `legacy_variant_selection_calculation` | 907 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:1225) | `analyze_legacy_variant_selection` | 1225 | Initialize runtime state and execute this component explicitly. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:1260) | `legacy_carriers_section` | 1260 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:1266) | `legacy_carriers_run` | 1266 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:1277) | `legacy_carriers_suppress` | 1277 | Legacy carriers suppress. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:1284) | `disclose_split` | 1284 | Disclose split. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:1296) | `ensure_mount` | 1296 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:1312) | `build_sample_list` | 1312 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:1341) | `build_qualifying_union` | 1341 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:1376) | `parse_ad` | 1376 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:1385) | `parse_gt_alleles` | 1385 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:1401) | `extract_gene_shard` | 1401 | Extract gene shard. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:1529) | `extract_all_genes` | 1529 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:1564) | `find_local_minimum` | 1564 | Find local minimum. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:1587) | `determine_qc_thresholds` | 1587 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:1640) | `apply_qc` | 1640 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:1665) | `collapse_carrier_status` | 1665 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:1702) | `step3_5_sparsity_precheck` | 1702 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:1824) | `legacy_carriers_write_report` | 1824 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:1910) | `update_config` | 1910 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:1940) | `legacy_carriers_calculation` | 1940 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:2016) | `analyze_legacy_carriers` | 2016 | Initialize runtime state and execute this component explicitly. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:2117) | `legacy_burden_section` | 2117 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:2123) | `legacy_burden_suppress` | 2123 | Legacy burden suppress. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:2130) | `load_weights` | 2130 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:2157) | `build_log_or_map` | 2157 | Return {gene: log(OR)} for the given panel, using the specified CHEK2 OR variant. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:2169) | `compute_burden` | 2169 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:2185) | `legacy_burden_calculation` | 2185 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:2419) | `legacy_burden_write_report` | 2419 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:2563) | `analyze_legacy_burden` | 2563 | Initialize runtime state and execute this component explicitly. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:2600) | `legacy_joint_models_log` | 2600 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:2604) | `legacy_joint_models_section` | 2604 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:2610) | `legacy_joint_models_suppress` | 2610 | Legacy joint models suppress. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:2617) | `p_floor` | 2617 | P floor. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:2626) | `balanced_weights` | 2626 | Balanced weights. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:2636) | `fit_unweighted_logit` | 2636 | Fit unweighted logit. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:2642) | `fit_weighted_logit` | 2642 | Fit weighted logit. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:2649) | `legacy_joint_models_coef_table` | 2649 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:2678) | `mcfadden_r2` | 2678 | Mcfadden r2. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:2690) | `zscore_fit` | 2690 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:2696) | `apply_z` | 2696 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:2705) | `legacy_joint_models_build_X` | 2705 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:2712) | `add_z` | 2712 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:2723) | `legacy_joint_models_save_smallest` | 2723 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:2743) | `assign_stratum` | 2743 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:2753) | `analyze_legacy_joint_models` | 2753 | Initialize runtime state and execute this component explicitly. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:3762) | `legacy_age_models_log` | 3762 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:3766) | `legacy_age_models_section` | 3766 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:3772) | `legacy_age_models_suppress` | 3772 | Legacy age models suppress. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:3779) | `join_covariates` | 3779 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:3791) | `zscore` | 3791 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:3796) | `fit_ols_hc1` | 3796 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:3801) | `legacy_age_models_coef_table` | 3801 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:3825) | `legacy_age_models_build_X` | 3825 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:3834) | `legacy_age_models_save_smallest` | 3834 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:3857) | `analyze_legacy_age_models` | 3857 | Initialize runtime state and execute this component explicitly. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:4394) | `legacy_figures_log` | 4394 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:4398) | `legacy_figures_section` | 4398 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:4404) | `disclose_n` | 4404 | Disclose n. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:4427) | `audit_override` | 4427 | Audit override. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:4438) | `fit_axes_to_hexbin` | 4438 | Fit axes to hexbin. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:4447) | `save_fig` | 4447 | Save fig. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:4482) | `write_table` | 4482 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:4493) | `fit_logit` | 4493 | Fit logit. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:4498) | `fit_logit_weighted` | 4498 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:4506) | `add_const` | 4506 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:4510) | `fisher_ci` | 4510 | Fisher ci. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:4519) | `residualize` | 4519 | Residualize. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:4525) | `corr_row` | 4525 | Corr row. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:4570) | `tier4` | 4570 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:4578) | `_split` | 4578 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:4594) | `cell_table` | 4594 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:4623) | `to_grid` | 4623 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:4627) | `fit_ventile_model` | 4627 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:4638) | `calib_params` | 4638 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:4645) | `analyze_legacy_figures` | 4645 | Initialize runtime state and execute this component explicitly. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:6454) | `term_kind` | 6454 | See the function and its analysis section. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:6462) | `analyze_legacy_forest` | 6462 | Initialize runtime state and execute this component explicitly. |
| [scripts/legacy/old_9gene_ptv_analysis.py](../scripts/legacy/old_9gene_ptv_analysis.py:6596) | `main` | 6596 | See the function and its analysis section. |
| [scripts/maintenance/organize_project.py](../scripts/maintenance/organize_project.py:26) | `digest` | 26 | See the function and its analysis section. |
| [scripts/maintenance/organize_project.py](../scripts/maintenance/organize_project.py:34) | `save_csv` | 34 | See the function and its analysis section. |
| [scripts/maintenance/organize_project.py](../scripts/maintenance/organize_project.py:42) | `main` | 42 | See the function and its analysis section. |
| [scripts/manuscript/collect_manuscript_data.py](../scripts/manuscript/collect_manuscript_data.py:19) | `main` | 19 | See the function and its analysis section. |
| [scripts/manuscript/figure_colors.py](../scripts/manuscript/figure_colors.py:51) | `geometry` | 51 | See the function and its analysis section. |
| [scripts/manuscript/figure_colors.py](../scripts/manuscript/figure_colors.py:56) | `recolor` | 56 | See the function and its analysis section. |
| [scripts/manuscript/figure_colors.py](../scripts/manuscript/figure_colors.py:74) | `redraw_figure6` | 74 | See the function and its analysis section. |
| [scripts/manuscript/figure_colors.py](../scripts/manuscript/figure_colors.py:192) | `redraw_review_figure1` | 192 | See the function and its analysis section. |
| [scripts/manuscript/figure_colors.py](../scripts/manuscript/figure_colors.py:240) | `main` | 240 | See the function and its analysis section. |
| [scripts/manuscript/interaction_assets.py](../scripts/manuscript/interaction_assets.py:40) | `sha` | 40 | See the function and its analysis section. |
| [scripts/manuscript/interaction_assets.py](../scripts/manuscript/interaction_assets.py:212) | `q` | 212 | See the function and its analysis section. |
| [scripts/manuscript/interaction_assets.py](../scripts/manuscript/interaction_assets.py:216) | `para` | 216 | See the function and its analysis section. |
| [scripts/manuscript/interaction_assets.py](../scripts/manuscript/interaction_assets.py:227) | `text` | 227 | See the function and its analysis section. |
| [scripts/manuscript/interaction_assets.py](../scripts/manuscript/interaction_assets.py:231) | `xml` | 231 | See the function and its analysis section. |
| [scripts/manuscript/interaction_assets.py](../scripts/manuscript/interaction_assets.py:235) | `load_doc` | 235 | See the function and its analysis section. |
| [scripts/manuscript/interaction_assets.py](../scripts/manuscript/interaction_assets.py:240) | `save_doc` | 240 | See the function and its analysis section. |
| [scripts/manuscript/interaction_assets.py](../scripts/manuscript/interaction_assets.py:282) | `append` | 282 | See the function and its analysis section. |
| [scripts/manuscript/interaction_assets.py](../scripts/manuscript/interaction_assets.py:286) | `page` | 286 | See the function and its analysis section. |
| [scripts/manuscript/interaction_assets.py](../scripts/manuscript/interaction_assets.py:292) | `addtable` | 292 | See the function and its analysis section. |
| [scripts/manuscript/interaction_assets.py](../scripts/manuscript/interaction_assets.py:313) | `picture` | 313 | See the function and its analysis section. |
| [scripts/manuscript/refresh_results.py](../scripts/manuscript/refresh_results.py:23) | `main` | 23 | See the function and its analysis section. |
| [scripts/manuscript/restyle_review_figure2.py](../scripts/manuscript/restyle_review_figure2.py:23) | `main` | 23 | See the function and its analysis section. |
| [scripts/manuscript/revise_13gene_1008.py](../scripts/manuscript/revise_13gene_1008.py:30) | `q` | 30 | See the function and its analysis section. |
| [scripts/manuscript/revise_13gene_1008.py](../scripts/manuscript/revise_13gene_1008.py:34) | `xml` | 34 | See the function and its analysis section. |
| [scripts/manuscript/revise_13gene_1008.py](../scripts/manuscript/revise_13gene_1008.py:38) | `sha` | 38 | See the function and its analysis section. |
| [scripts/manuscript/revise_13gene_1008.py](../scripts/manuscript/revise_13gene_1008.py:49) | `ptext` | 49 | See the function and its analysis section. |
| [scripts/manuscript/revise_13gene_1008.py](../scripts/manuscript/revise_13gene_1008.py:77) | `add` | 77 | See the function and its analysis section. |
| [scripts/manuscript/revise_13gene_1008.py](../scripts/manuscript/revise_13gene_1008.py:359) | `rendered` | 359 | See the function and its analysis section. |
| [scripts/manuscript/revise_13gene_1008.py](../scripts/manuscript/revise_13gene_1008.py:648) | `para` | 648 | See the function and its analysis section. |
| [scripts/manuscript/revise_13gene_1008.py](../scripts/manuscript/revise_13gene_1008.py:661) | `table` | 661 | See the function and its analysis section. |
| [scripts/manuscript/revise_13gene_1008.py](../scripts/manuscript/revise_13gene_1008.py:685) | `write_doc` | 685 | See the function and its analysis section. |
| [scripts/manuscript/revise_13gene_1008.py](../scripts/manuscript/revise_13gene_1008.py:747) | `markdown` | 747 | See the function and its analysis section. |
| [scripts/review/probe_1008.py](../scripts/review/probe_1008.py:17) | `main` | 17 | See the function and its analysis section. |
| [scripts/review/review_1008.py](../scripts/review/review_1008.py:53) | `write_table` | 53 | See the function and its analysis section. |
| [scripts/review/review_1008.py](../scripts/review/review_1008.py:58) | `save_json` | 58 | See the function and its analysis section. |
| [scripts/review/review_1008.py](../scripts/review/review_1008.py:65) | `estimate` | 65 | See the function and its analysis section. |
| [scripts/review/review_1008.py](../scripts/review/review_1008.py:77) | `ratio_estimate` | 77 | See the function and its analysis section. |
| [scripts/review/review_1008.py](../scripts/review/review_1008.py:81) | `disclosed_size` | 81 | See the function and its analysis section. |
| [scripts/review/review_1008.py](../scripts/review/review_1008.py:89) | `zscore` | 89 | See the function and its analysis section. |
| [scripts/review/review_1008.py](../scripts/review/review_1008.py:94) | `load` | 94 | See the function and its analysis section. |
| [scripts/review/review_1008.py](../scripts/review/review_1008.py:134) | `fit_logistic` | 134 | See the function and its analysis section. |
| [scripts/review/review_1008.py](../scripts/review/review_1008.py:145) | `interaction` | 145 | See the function and its analysis section. |
| [scripts/review/review_1008.py](../scripts/review/review_1008.py:321) | `query` | 321 | See the function and its analysis section. |
| [scripts/review/review_1008.py](../scripts/review/review_1008.py:341) | `dates_task` | 341 | See the function and its analysis section. |
| [scripts/review/review_1008.py](../scripts/review/review_1008.py:390) | `riskset_moments` | 390 | Breslow risk sets entry < t <= exit, via sorted cumulative moments. |
| [scripts/review/review_1008.py](../scripts/review/review_1008.py:420) | `prepare_survival` | 420 | See the function and its analysis section. |
| [scripts/review/review_1008.py](../scripts/review/review_1008.py:451) | `fit_cox` | 451 | See the function and its analysis section. |
| [scripts/review/review_1008.py](../scripts/review/review_1008.py:493) | `decompose` | 493 | See the function and its analysis section. |
| [scripts/review/review_1008.py](../scripts/review/review_1008.py:506) | `frailty_hr` | 506 | See the function and its analysis section. |
| [scripts/review/review_1008.py](../scripts/review/review_1008.py:510) | `split_intervals` | 510 | See the function and its analysis section. |
| [scripts/review/review_1008.py](../scripts/review/review_1008.py:523) | `simulate_case_coefficients` | 523 | Conditional survival after entry; fixed score scaling from observed cases. |
| [scripts/review/review_1008.py](../scripts/review/review_1008.py:557) | `cox_task` | 557 | See the function and its analysis section. |
| [scripts/review/review_1008.py](../scripts/review/review_1008.py:810) | `mh_estimate` | 810 | See the function and its analysis section. |
| [scripts/review/review_1008.py](../scripts/review/review_1008.py:838) | `classify_family_flags` | 838 | See the function and its analysis section. |
| [scripts/review/review_1008.py](../scripts/review/review_1008.py:851) | `family_task` | 851 | See the function and its analysis section. |
| [scripts/review/review_1008.py](../scripts/review/review_1008.py:1002) | `plots_task` | 1002 | See the function and its analysis section. |
| [scripts/review/review_1008.py](../scripts/review/review_1008.py:1197) | `frailty_task` | 1197 | See the function and its analysis section. |
| [scripts/review/review_1008.py](../scripts/review/review_1008.py:1307) | `age_pattern_task` | 1307 | Calibrate gamma depletion from non-carrier PRS age effects as in PDF. |
| [scripts/review/review_1008.py](../scripts/review/review_1008.py:1399) | `validate_task` | 1399 | Compare rounded published targets without claiming method equivalence. |
| [scripts/review/review_1008.py](../scripts/review/review_1008.py:1666) | `main` | 1666 | See the function and its analysis section. |

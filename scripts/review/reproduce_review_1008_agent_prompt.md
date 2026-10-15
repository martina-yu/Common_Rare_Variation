# Agent prompt: reproduce the 1008 review analyses

Work in `/home/jupyter/workspace/com_rare_BC`. Reproduce the analyses and six figures in `papers/Review Decisions, 1008 Draft.pdf`, pages 1–6, using this project's real data. Read the methods and figure captions; do not process the decision cards, select revision options, or edit the manuscript. Preserve the review's formulas, analysis models, cohort definitions and sensitivity designs. Plotting library, functions and colours may differ. The numerical values in the report are validation targets, never substitute inputs or hard-coded results.

## First establish the source and inputs

Read applicable AGENTS.md files and CLAUDE.md, scripts/README.md, scripts/FUNCTION_INDEX.md, and the current configuration. Inspect scripts/07_burden_score.py, scripts/08_case_control_models.py, scripts/exploratory/age_ehr_checks.py, scripts/11_interaction_analysis.py and their provenance before reusing calculations. Freeze and record the current thirteen-gene PTV plus eligible P/LP missense definition, gene weights, carrier indicators, full-cohort score means/SDs, cohort IDs and release. Do not silently replace the review's score with the legacy PTV score.

Find the original reviewer implementation if accessible: repository `git@github.com:YikaiDong-git/Collaborate_breastcancer_Zimeng.git`, branch main, folder `breast_cancer_intera/review_c3_c4_c9/`, reference commit `ca6bfba`. Its `wf/`, README.md, figure captions and family-history report can resolve details missing from the PDF. Record the actual source commit. Do not overwrite this project's working tree to obtain it. If unavailable, proceed with identifiable analyses and list unresolved specifications explicitly. Some formulas on PDF page 3 render incorrectly even visually; verify them from source rather than claiming an exact transcription.

Create an input/method mapping with one row per analysis: source variables and files, cohort/filter, model, uncertainty estimator, output and any unresolved detail. Resolve dates, EHR observation windows, death censoring, survey concept IDs, variant definitions, tie handling and simulations before asserting exact reproduction. Missing input must produce an explicit incomplete status, not fabricated data.

## A. Carrier-stratified interaction

Use PRS311 standardised once in the full analysis cohort. Within stratum s fit

`logit Pr(Y_i=1) = alpha_s + beta_s P_i + gamma_s^T Z_i`,

where Z contains age at last EHR observation, PC1–10, sequencing centre and sequencer; use HC1 standard errors. Output OR_s = exp(beta_s), 95% CI, and carrier/non-carrier ratio = exp(beta_c-beta_n), with SE of the log ratio = sqrt(SE_c^2+SE_n^2) for the disjoint stratum fits. Also fit BRCA1/2 and other-panel carrier strata, documenting overlapping/multigene handling from source.

Fit the paper model `logit Pr(Y=1)=alpha+beta_P P+beta_B B+beta_PB P B+gamma^T Z`. Explain the carrier PRS slope `beta_P+beta_PB B_g` using the carried gene's burden on the same standardised scale. Compare PRS interactions with weighted burden, any carrier, gene class and carried-gene weight. Use identical observations for AIC comparisons; reproduce source-defined nested likelihood-ratio comparisons rather than comparing non-nested models by LRT.

Reproduce every Figure 2 row: primary covariates; spline age; no covariates; within-European PCs; diagnosis age for cases; exclude multigene carriers; PTV only; legacy ten-gene panel; rule-of-two case definition; drop cases first coded at first EHR record; split cases before/after enrolment age + 1 year; exclude controls self-reporting breast cancer; exclude mastectomy-related records (controls any such code, cases at least six months before diagnosis); exact-date pre-enrolment cases and cases at least one year after enrolment. Match source spline specification, code sets and denominator construction. Keep age-based and exact-date sensitivities separate.

Targets: logistic PRS OR 1.56 (1.52–1.60) in non-carriers, 1.25 (1.14–1.37) in carriers; ratio 0.80 (0.73–0.88). BRCA1/2 OR 1.11, other genes 1.50. AICs approximately 55823, 55831, 55826, 55825 for the four interaction forms. Report discrepancies without tuning models to achieve these numbers.

## B. Cox models and expected case-only diagnosis-age coefficients

Fit attained-age Cox models with delayed entry: entry max(age 18, first EHR age + 1 year); exit at first qualifying breast cancer record, last EHR record or death. Apply valid interval and prevalent-case exclusions explicitly. Model `h_i(t)=h_0(t) exp(beta^T x_i)`. Reproduce the second entry design and the exact-date EHR-only design with birth-decade strata from the original reviewer code; these details are not fully specified in the PDF. Use Cox covariates from source rather than automatically carrying an attained-age covariate into the age-scale model.

At each event age t_i compute hazard-weighted risk-set mean

`m(t_i) = sum_{j in R(t_i)} x_j exp(eta_j) / sum_{j in R(t_i)} exp(eta_j)`,

and residual `r_i=x_i-m(t_i)`. Use the same case-only Table 3 design and adjustment as the review. With centred variables (and proper residualisation for nuisance covariates/intercept), reproduce

`b_obs = (Z^T Z)^(-1) sum_i m(t_i)(t_i-t_bar) + (Z^T Z)^(-1) sum_i r_i(t_i-t_bar)`.

The first term is b_exp. Approximate remainder covariance by

`(Z^T Z)^(-1) [sum_i (t_i-t_bar)^2 V(t_i)] (Z^T Z)^(-1)`,

where V(t_i) is the hazard-weighted risk-set covariance. Confirm centring, score column selection, nuisance adjustment and ties from source. Check numerically that expected plus remainder equals observed under the implemented design. Do not conflate source notation Z here with the logistic nuisance covariates.

Simulate cohorts under the fitted age-constant model to verify b_exp, preserving delayed entry and censoring as specified by source. Record random seeds, baseline hazard generation and Monte Carlo error. Reproduce observed-versus-expected tests, not just tests against zero.

Split follow-up intervals at ages 45, 55 and 65. Fit `h_i(t)=h_0(t) exp(sum_k beta_k x_i I(t in A_k)+gamma^T z_i)` for age bands <45, 45–54, 55–64, >=65. Estimate PRS311, any-carrier, BRCA1/2 and other-carrier age-specific HRs and CIs using the source's separate/joint model specification.

Targets: b_exp PRS -0.23 to -0.30 yr/SD and burden -0.72 to -0.80; b_obs -0.43 and -1.82; against-expectation P approximately 0.29 and 3e-14. Main PRS Cox HRs 1.51 and 1.20 in non-carriers/carriers; exact-date HRs 1.52 and 1.23. Carrier HR approximately 8.4 before 45 and 2.7 at >=65.

## C. Gamma-frailty extension

Use gamma frailty with mean 1 and variance theta, `u ~ Gamma(shape=1/theta, rate=1/theta)`. Under constant individual HR exp(beta), the standard marginal depletion expression is `HR_pop(t)=exp(beta)*(1+theta H0(t))/(1+theta H0(t) exp(beta))`; this expression must be checked against reviewer source because the PDF fails to show it clearly. Establish the reference baseline hazard and covariate conditioning before using it.

Reproduce theta estimation from the PRS311 age pattern in non-carriers, including its confidence limit (report targets theta=0.28, upper 95% bound=0.73). Generate the source-defined 100 simulated cohorts for the frailty extension, and the BRCA1/2 depletion curve at theta=0.73. Treat theta=0 as the no-frailty limit. Report identifiability and unresolved implementation differences. Target at theta=0.73: burden coefficient -1.78 observed, -1.45 expected; do not assume these coefficients use the same case set as Table 3.

## D. Family-history survey routing

Reconstruct the existing Table 1 item 836772 coding and contrast it with corrected coding separately for three versions: per-relative 2018–2021; conditions-first 2020–2024; personal-and-family 2021–2024. Positive means a parent, sibling, child or grandparent named for breast cancer. Negative means the version's family cancer questions were answered without such a relative. Self-only and routed/skipped answers must be classified according to each version; absence of a response alone is not evidence of a negative answer.

Each woman contributes her first submission per version to version-specific results. For the pooled result use each woman's first eligible survey across versions, without duplicating women. Verify age reference and pooling strata from source. Age bands are <50, 50–59, 60–69, >=70. With a/b positive/negative cases and c/d positive/negative controls in stratum k, reproduce Mantel–Haenszel `OR_MH=sum_k(a_k d_k/n_k)/sum_k(b_k c_k/n_k)`, including the source's CI estimator.

Reproduce Figure 5 routes in the 2021 version: relative named; only self; breast cancer ticked but relative question skipped; list answered without breast cancer; list unanswered. Targets: existing pooled item OR 0.71 (0.68–0.75), existing 2021 age-adjusted OR 0.23 (0.21–0.25); corrected version ORs 1.83, 2.29, 1.61; pooled 1.76 (1.66–1.86). The reference 1.80 (1.69–1.91) is a published risk ratio with a 99% CI, not a new estimated OR; label it accurately.

## Outputs and verification

Add reproducible executable analysis scripts in scripts/, using a shared prefix such as 62_review_1008 and separate tasks for interaction, Cox/decomposition, frailty, family history and plotting. Put aggregate tables in results/table/, six editable SVG figures in results/figures/, and methods, captions, provenance and validation reports in results/report/. Use the project's environment and workflow conventions; record new dependencies. Participant-level inputs/intermediates stay in ignored internal data locations and are never exported to results or logs.

Figure 1: forest plot of stratum PRS ORs and HRs for three designs. Figure 2: forest plot of sensitivity carrier/non-carrier OR ratios. Figure 3: observed coefficients with CIs against no-frailty and frailty expected ranges. Figure 4: age-band HRs with CIs and frailty depletion curve. Figure 5: stacked route percentages by case status. Figure 6: existing/corrected family-history estimates and published reference. Use standard matplotlib/R defaults, no grid, accurate log axes, explicit reference lines and saved underlying aggregate plotting data; visual styling need not match.

For this review's outputs suppress counts 1–20 inclusive, including logs, captions and audit tables; inspect complementary disclosure through totals, proportions and related tables. Do not reveal suppressed counts in validation diagnostics. Do not change older project outputs merely to apply the review's rule.

Test substantive calculations using small synthetic examples: stratum ratio and CI; delayed entry/censoring and interval splitting; weighted risk-set means/covariances and decomposition identity; theta=0 limit; MH formula and survey deduplication/routing. Then run real-data analyses and compare each target with source-matched denominators and definitions. Record achieved, different, unavailable or unresolved per analysis; document causes of differences. Never mark reproduction complete based only on script creation or plots of PDF numbers.

Update scripts/README.md and results/README.md with execution order and output paths. Finish with a concise account of generated analyses, figures, validation and missing dependencies or specifications. Do not commit, push, or send messages as part of this task.

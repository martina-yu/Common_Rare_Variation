# Rare variant burden complements polygenic risk in breast cancer diagnosis and is associated with earlier diagnosis

**Preprint:** [Read on medRxiv](https://www.medrxiv.org/content/10.64898/2026.10.17.26365799v1) · Version 1

This repository provides public analysis code accompanying the preprint. See [Analysis methods and usage](#analysis-methods-and-usage) below for the workflow, dependencies, and instructions.

## Overview

This study uses whole-genome sequencing and electronic health records from the All of Us Research Program to examine the independent and joint associations of common polygenic susceptibility and rare germline variation with breast cancer case status and recorded age at diagnosis among cases. Common variation is summarized by PRS311, and rare variation by a weighted burden score across 13 breast cancer susceptibility genes, including protein-truncating variants and qualifying pathogenic/likely pathogenic missense variants. PRS311 is derived from the BCAC score and contains 302 scored loci after excluding two CHEK2-tagging variants to reduce double counting.

The study includes participants recorded as female at birth with genetically inferred European ancestry. After quality control, the cohort comprised **8,528 breast cancer cases and 134,727 controls**. Case-status models included **141,204 participants**, and the diagnosis-age analysis included **7,996 cases**.

## Main findings

- **Common and rare variants provide complementary information about case status.** In the joint model, the odds ratios (ORs) per standard deviation were **1.55** for PRS311 and **1.27** for rare-variant burden, each evaluated at the other score's mean. Cross-validated AUC increased from **0.745** for covariates alone to **0.782** after adding both genetic scores and their interaction.
- **The joint association is submultiplicative.** The PRS × rare-variant burden interaction OR was **0.972**, with a weaker relative PRS association among carriers. This statistical interaction on the odds-ratio scale does not establish biological antagonism.
- **Rare-variant burden has a stronger association with earlier recorded diagnosis.** Among cases, each standard deviation increase in burden and PRS311 was associated with diagnosis **1.82 years** and **0.43 years** earlier, respectively. Together, the scores explained **2.63%** of the variation in recorded diagnosis age.
- **Rare-variant carriers are more frequent among cases diagnosed at younger ages.** Carrier frequency declined from **26.5%** among cases diagnosed before age 35 to **3.2%** among those diagnosed at age 70 or older.

These findings support considering common and rare inherited variation together and suggest different genetic association patterns for case status and recorded diagnosis age. Interpretation is limited by uncertainty in EHR-derived diagnosis dates, restriction to European ancestry, post hoc components of the expanded analyses, and pending final manual variant adjudication. AUC includes contributions from age and the observation window, so increments between models are more informative than absolute AUC. The findings do not establish prospective absolute-risk prediction or causal effects.

## Figure 1. Cohort and genetic analysis workflow

[View Figure 1 in the preprint](https://www.medrxiv.org/content/10.64898/2026.10.17.26365799v1#F1)

**Caption:** (a) Ancestry, relatedness, and quality control filtering yielded 8,528 breast cancer cases and 134,727 controls. (b) PRS311 combines 302 scored loci; curves show case and control score distributions. (c) Protein-truncating and qualifying pathogenic/likely pathogenic missense variants define gene-level carrier indicators and the weighted rare-variant burden. (d) Logistic regression with repeated five-fold cross-validation assesses case status, and case-only linear regression assesses recorded age at diagnosis.

## Analysis methods and usage

This directory contains method-only adaptations of the study scripts. It contains
no participant data, fitted results, cohort counts, observed variant lists,
individual scores, intermediate files or plotting code. Original scripts are
retained separately. This is a refactored methods package, not a byte-for-byte
reproduction of a historical execution or a complete cloud orchestration system.

### Workflow

| Stage | Methods |
|---|---|
| 01 | Female WGS cohort; six-class ancestry posterior >0.75; breast-cancer/EHR phenotypes |
| 02 | EUR restriction; greedy kinship >0.1 pruning; chr21 missingness/HWE proxy; clinical covariates |
| 03 | PGS000004 allele matching, constrained AF ambiguity resolution, PRS313/PRS311, PLINK comparison command |
| 04 | Nine- and thirteen-gene panels, transcript collapse, public combined external OR weights |
| 05 | HC LOFTEE PTVs; germline ClinVar P/LP missense, disease/inheritance review, consistent EUR AF <0.01 |
| 06 | Complete normalized allele joins; original ALT indices; site/FT/GQ/AD/VAF QC; gene carrier collapse |
| 07 | Sum of external log(OR) per carried gene, once per gene; score-version decomposition |
| 08 | HC1 logistic models; paired five-model repeated stratified CV |
| 09 | Case-only recorded diagnosis-age HC1 additive OLS; HC3/HC1 continuous interaction |
| 10 | Without-STK11, PTV-only, strict TP53/NF1 VAF, >=1-star missense, optional without-CDH1 |
| 11 | Post hoc interaction, influence/leverage exclusion, splines and gene-specific/pooling methods |

Order is not simply numerical: prepare phenotype and sample eligibility, calculate
PRS, attach clinical/technical covariates and freeze model membership. Prepare
panels, annotation, QC and burdens; freeze sensitivity definitions (10) before
fitting (08/09). Keep all score versions on identical members. Case-only age
membership additionally excludes a diagnosis equal to the first EHR date.

Sample missingness uses three 500-marker chr21 windows and excludes proxy markers
with missingness >0.2 before sample missingness >0.03. This is not genome-wide QC.
Control-only PLINK --hardy P<1e-6 is a proxy diagnostic. Kinship pruning is greedy,
not a certified maximum independent set, and >0.1 misses weaker second-degree
pairs between 0.0884 and 0.1.

Genotype QC requires original site FILTER empty/PASS, FT PASS or '.', GQ>=20,
heterozygous target ALT VAF in [0.3,0.8], homozygous VAF>=0.8, denominator=sum of
all AD entries. The strict TP53/NF1 sensitivity caps heterozygous VAF at 0.7.
DP>=10 and an additional terminal-exon exclusion are not main eligibility rules.
Normalize biallelic site copies against GRCh38 while retaining original GT/AD ALT
indices; VDS LGT must first map through LA. Missing site QC information fails.

PTVs require PTV consequence and HC on the same retained same-gene transcript.
Canonical rows supply representative fields, not an extra HC-PTV restriction.
Missense uses unambiguous canonical consequence, exact germline P/LP, >=2 stars,
blanket low-penetrance/risk/conflict/B-LB exclusions and relevant disease/inheritance
review. Missing/disagreeing EUR AF is held out. Final manual adjudication remains
required; the extension retains its post hoc/provisional identity.

Published weights are approximate combined PTV+non-PTV PV ORs from
[Rowlands et al., Table 2](https://doi.org/10.1016/j.annonc.2024.07.244); their
variant definitions and ClinVar dates differ from this missense protocol. The
literal weights are inherited method inputs, not newly independently certified
against the source article in this code-publication task. Legacy comparison
carriers/weights must be supplied privately; legacy missing-AF handling must not
silently propagate into the new main score.

Saved burden standardization uses full frozen case/control ddof=0 and reuses
that reference in age subsets. Case/control fits use full-model-cohort ddof=1;
age fits use age-case ddof=1; CV learns mean and ddof=1 SD on each training fold.
Five models are covariates, PRS, burden, additive PRS+burden, and PRS+burden+product.
CV uses the same case x legacy-carrier strata/folds for every model and version.
The 10-repeat mean and repeat range describe split stability; a first-repeat
DeLong interval is not a confidence interval for the ten-repeat mean or complete
training process. No additional bootstrap is implied or implemented.

### Use

Install dependencies from environment.yaml in an authorized environment. From
this directory, `python load_stage.py --list` lists methods without accessing data.
Numbered modules expose in-memory APIs:

```python
from load_stage import load_stage
prs_method = load_stage('prs')
# effect_dosage, weights and tag_mask are supplied inside authorized storage.
scores = prs_method.score_prs(effect_dosage, weights, tag_mask)
```

Configure every real location privately. All shipped storage/project/dataset
values are `path/to/folder`; each key must be replaced independently with the
appropriate complete file, table or URI, not interpreted as a shared real path.
OMOP schema table/column names and public resource URLs are method information.
No private configuration or actual folder hierarchy is included.

functions/queries.py contains parameterized cohort/EHR/clinical SQL. Functions
accept authorized in-memory tables; the user must perform authenticated queries,
normalize and reconcile source alleles, supply frozen IDs/folds where needed and
retain any generated data inside the authorized environment. GCS helpers only
read; they do not copy, move or export data. PLINK command builders require
explicit private paths; all resulting files remain controlled intermediates.

The original manuscript assembly, figure rendering, run-specific provenance
locks, real-data acceptance constants and report-writing code are excluded.
There is no import-time execution and no default pipeline that creates data
files. The historical burden report-writing failure is avoided by omitting the
report-writing path; it is not claimed fixed in the original private scripts.

### Sharing and validation

Returned data frames, model objects and prediction arrays can contain controlled
information even without an identifier. Do not put them in this source directory
or distribute them. This package does not automatically release results.
[AoU policy](https://www.researchallofus.org/faq/data-and-statistics-dissemination-policy/)
protects counts 1--20 inclusive and prohibits reconstructing them from related
information. functions/disclosure.py uses a conservative related-cell mask; it
cannot certify arbitrary overlapping releases or percentages/confidence intervals.

Only code, documentation and empty path configuration belong in this package.
No figures, real output, notebooks, logs or intermediates accompany it. Synthetic
checks exercise numerical and disclosure behavior; successful checks do not
constitute an authorized-data end-to-end rerun or official AoU release approval.

Optional Cox/frailty, review orchestration and four-group descriptive analyses
are outside this core package. They were marked conditional in the inventory;
they are not silently promoted to current main analyses.

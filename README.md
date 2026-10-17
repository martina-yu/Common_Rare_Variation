> 当前稿件：作者于2026-10-16确认 [main 1012](manuscript/medrxiv_v1_main_1012.docx) / [supplement 1012](manuscript/medrxiv_v1_supplement_1012.docx) 为最新版本；1008为历史。1012数值与结果登记的一致性待核对。

> 当前交接入口：[PROJECT_STATUS](PROJECT_STATUS.md)顶部一分钟摘要；执行规则：[CLAUDE](CLAUDE.md)；每次工作/运行：[ANALYSIS_CHANGELOG](ANALYSIS_CHANGELOG.md)；数值：[RESULT_REGISTRY](RESULT_REGISTRY.md)；待办：[NEXT_STEPS](NEXT_STEPS.md)。下方研究叙述及SUMMARY/旧plans为历史背景，不能据此选当前数值。全基因组PRS v2仅为待锁定计划。

> Readable analysis layout (2026-10-15): [project status](PROJECT_STATUS.md), [result registry](RESULT_REGISTRY.md), [analysis scripts](scripts/README.md), and [validation](docs/readability_validation.md). Existing data/result paths are preserved.

> 当前脚本与结果导航（2026-10-13）：[项目现状](SUMMARY.md)、[脚本/任务/结果索引](docs/project_catalog/README.md)、[整理方法与复现边界](tasks/METHODS.md)。下文保留早期研究背景和历史记录；当前运行入口见 [scripts/README.md](scripts/README.md)。

### Project Background

Description of Proposed Project (250 words)
The objectives or goals of the project. One in every three person develops cancer
during their life time. Why some people develop cancer and what determins the age
of onset of cancer is largely unknow. Most individuals do not have germline mutations
in high penetrance cancer susceptibility genes, and many common cancers are not
associated with known carcinogen exposure (e.g.breast or prostate cancreds). We
hypothetise that the combination of common and rare germline line variants that
alter protein structure and function in cancer relevant genes determine life time
cancer risk. The higher the variant burden the greater the risk and the younger the
age at diagnosis. Unlike polygenic risk score, our model does not assume that an
individiual variant has to be associated with cancer risk, it is the combined effect that
detetrmine risk. The goal of this project is to test this hypothesis using already
established clinically annotated germline whole exome sequencing data under the
supervision of Professors Lajos Pusztai and Matrk Gerstein
A brief description of data and methods that you will use.
We will use germline variant data and clinical information from the UK Biobank,
Clinical Proteomic Tumor Analysis Consortium (CPTAC), and Yale
GENERATIONS projects. We will exclude individuals with known germline
pathogenic variants of the 47 established cancer susceptibility genes and
partition the data sets into training (2/3) and validation cohorts (1/3). Variant-,
gene- and pathway- level metrics will be used as input variables (for details see
joint publication by mentors: Qing T. et al. Cancer relevance of human genes.
JNCI; 11:988-95). We will apply fully connected dense neural networks to
predict (i) overall cancer risk, (ii) age of onset of breast cancer, and (iii)
compare our model with established polygenic risk scores for breast cancer. An
overall schema is shown below:

Step 1
Calculate a weighted score per gene for each patient

Step 2
Predict the age at breast cancer diagnosis

Key milestones or deliverables you aim to achieve.
QC data and clean data base. Establish cross data set comparability
Train risk prediction model
Compare predictive accuracy (AUG, calibration, discrimination) with polygenic
risk scores.

# I. Research Background and Question

- **Clinical observation**: Roughly 1/3 of people will develop cancer in their lifetime, but the drivers and the mechanisms determining age of onset remain unclear.
- Existing gaps:
  - Most patients **do not** carry germline mutations in high-penetrance cancer susceptibility genes
  - Many common cancers (e.g., breast, prostate) **lack** a clear history of carcinogen exposure
  - Traditional GWAS/PRS assumes "each individual variant must be independently associated with risk," yet many individuals with high-risk pathogenic variants never develop cancer in their lifetime, while many individuals with low PRS do develop cancer — suggesting existing models are missing some mechanism



Control定义最终确定:white+female+WGS+无任何癌种诊断码,不设EHR覆盖/年龄硬门槛,dx_code_occurrence_count/ehr_length/age_at_last_obs都作为Step
  5协变量;原方案B(EHR门槛)和≥50/55年龄门槛降级为敏感性分析。主分析control人数定为144,841,Case:Control=1:15.6(替换了原来贯穿全文的1:17.8),已同步更新到Step 3的SAIGE-GENE+和Step
  3.5功效预检两处引用。

### II. Core Hypothesis

> **Lifetime cancer risk is determined by the "combined burden" of common variants and rare germline variants that alter protein structure/function in cancer-related genes — not by the independent effect of any single variant. The higher the burden, the greater the risk, and the earlier the age of onset.**

Key distinction from classical PRS: **it does not require that any individual variant itself be significantly associated with cancer risk**; the emphasis is on the **cumulative effect (burden)** of variants, rather than "association" per se.

Reference 1 (Souaiaia et al., 2026), under the POPout framework, combines common PRS with rare variant information (from three sources: moderately rare, very rare, and gene-burden) into a single joint PRS, used to test whether tail deviation can be explained by rare variants. This study does not adopt that combining strategy; instead it enters common PRS and rare burden into the joint model as two separate components and explicitly tests an interaction term. The reason is that this study's core question (particularly the H1 mechanism-separation hypothesis) is concerned with the independent effects of the two genetic components and their relationship to each other — information that would be erased if collapsed into a single score. Additionally, Reference 1's combining design serves POPout's test of tail deviation for continuous phenotypes, which is inconsistent with this study's goal of testing joint/interaction effects under a binary cancer-diagnosis outcome.

#### Supporting Sub-hypotheses

| Sub-hypothesis                      | Content                                                      |
| ----------------------------------- | ------------------------------------------------------------ |
| H1: Mechanism-separation hypothesis | Common and rare variants represent **two distinct pathogenic mechanisms**: common variants depend more on gene-environment interaction (e.g., smoking amplifying their effect); rare variants, being constrained by negative selection, tend to have a larger impact on protein function once present |
| H2: Population genetics hypothesis  | Across generations of selection, inherited common pathogenic SNPs tend to have milder effects (not severe enough to be lethal or to impair fertility), and therefore **cannot cause disease on their own**; rare variants have not been fully purged by selection and have stronger individual effects, though not necessarily lethal either |
| H3: Age-stratification hypothesis   | Late-onset population: genetic burden itself is lower, lifetime risk accumulates slowly, and onset depends more on environmental/somatic random mutation accumulation; early-onset population: more likely reflects a **mixture** of common + rare burden, causing genetic burden to exceed the threshold early |
| H4: Somatic compensation hypothesis | High PRS but no cancer → possibly low somatic mutation burden/environmental exposure, with the two never combining to trigger disease; low PRS but cancer → possibly high rare variant burden, or somatic mutations being the dominant contributor |

### III. Specific Aims

**Aim 1: Construct a Combined Germline Variant Burden Model**

**Cohort finalized (updated 2026-08)**: The All of Us WGS white female breast cancer sub-cohort, totaling 174,381 individuals (9,281 cases / 165,100 controls). WGS naturally supports both common PRS calculation and rare variant burden analysis simultaneously, resolving the previously anticipated concern that "common PRS and rare burden need to be computed on the same cohort to be combined." The Pusztai/Gerstein WES data (i.e., a curated subset of TCGA WES data) has been downgraded to an independent validation cohort and does not block initiation. Execution details are in `plans/aim_1.md`.

The All of Us WGS data is split into three "streamlined callsets" for analysis:

| Callset                    | Definition                                                   | Purpose                                                      |
| -------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| **ACAF threshold callset** | Variants with population-specific allele frequency (AF) >1% or allele count (AC) >100 | Common variant/PRS analysis (REGENIE step1 background association matrix) |
| **Exome callset**          | Variants located within exonic regions of Gencode v42 basic transcripts, with 15bp padding on each side (extracted from WGS, not from a separate WES experiment) | Rare coding variant burden analysis                          |
| **ClinVar callset**        | Variants recorded in the ClinVar database, regardless of pathogenicity | Cross-reference against known pathogenic variants            |

**Aim 2: Analyze Differences in Common/Rare Variant Distribution Across Populations** Compare whether common and rare variant distribution patterns differ systematically between cancer patients and healthy controls; after stratifying by PRS, test for a statistical interaction between the two.

Pain points of rare variant analysis:

rare variant discovery remains limited by sample size, sequencing heterogeneity across cohorts, and the requirement for sufficient numbers of variant carriers for gene-based burden testing. Population-scale sequencing in large biobanks, combined with proxy phenotypes, provides an opportunity to increase statistical power while minimizing sequencing heterogeneity

- Definition of "rare" (MAF<1%):
  - Common-variant PRS methods typically include only MAF≥1% variants, while the contribution of MAF<1% variants remains to be explored. The QC pipeline excludes sites with MAF greater than 0.01, MAC<10, missingness >0.1, or HWE p-value <1e-15; all remaining **MAF<1%** rare variants are included in the analysis and further subdivided.
  - Three functional annotation categories (LoF only; missense+LoF; missense+LoF+synonymous) are crossed with three MAF thresholds (0.0001, 0.001, 0.01) to form 9 variant masks. The single-variant rvPRS uses variants with MAF between 0.0001 and 0.01, tested via SAIGE for single-variant association.

**Aim 3: Integrate Germline–Somatic Interaction Analysis (Environmental Factors Included Pending Feasibility)** Using somatic mutation signatures as a phenotype/intermediate variable, test whether germline burden (common+rare) can predict an individual's somatic mutation characteristics.

**Aim 4: Comparison of Drivers Across Age Strata** Group by age of onset (early- vs. late-onset), compare differences in rare variant burden, common PRS, and somatic mutation burden between the two groups, and test the hypothesis that "younger patients rely more heavily on the mixture mechanism."

### IV. Data and Cohorts

| Dataset                                                   | Purpose                                                      | Notes                                                        |
| --------------------------------------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| **All of Us Research Program**                            | Primary cohort, used for common PRS calculation + rare variant burden analysis | The white female breast cancer sub-cohort has been confirmed via Cohort Builder: **174,381 individuals (9,281 cases / 165,100 controls)**. An independent study has already successfully trained and validated a breast cancer genetic risk model in the All of Us population (AUC of 0.70 in the white female sub-cohort), demonstrating that this database supports breast cancer PRS calculation |
| **UK Biobank (UKBB)**                                     | Candidate validation cohort, large sample size (~500,000 individuals), provides both WES and imputed genotype data, favorable for rare variant association testing | Access expected to be available **2026-09**; not currently blocking, priority remains on executing the All of Us analysis |
| TCGA (i.e., the curated Pusztai/Gerstein WES data subset) | Downgraded to an **independent validation cohort**; not used for PRS calculation in the primary analysis. The WES + SNP6.0 array imputation combination can approximate a common PRS, useful for cross-validating the Aim 1 joint model; also remains a cross-reference resource for rare pathogenic variants. **Note**: TCGA consists predominantly of tumor samples and lacks healthy controls — this limitation must be assessed during validation | See "Section V: Methodological Notes" for rationale; cohort role detailed in the status table of `plans/aim_1.md` |
| Breast Cancer PRS313                                      | Standard reference scoring system for common-variant PRS     | From the Breast Cancer Association Consortium (BCAC), comprising 313 common variants (267 SNPs + 46 indels) |

**Cohort restrictions (initial phase)**:

- Ancestry: restricted to self-reported white (cross-ancestry analysis not performed at this stage, reserved as a future extension)
- Sex: **only female breast cancer patients included** (the study focuses on breast cancer; male breast cancer and other rare subtypes require separate analysis)

### V. Methodological Note: Why TCGA Is Not Used for PRS Calculation

This is an important methodological decision in the project, worth explaining separately:

The loci underlying common-variant PRS models such as PRS313 are, for the most part, located in **introns, intergenic regions, and regulatory regions** rather than coding exons — a statistical property inherent to GWAS association signals themselves. TCGA, however, provides germline WES data, and standard exome-capture kits cover only about 1.7% of the genome (the coding exonic region). As a result, **the majority of PRS313 loci fall entirely outside the WES capture range** (empirically, only about 38 of 313 loci can be matched, consistent in order of magnitude with the capture-region proportion).

Even on a standard genotyping array, only 105 of PRS313's 313 loci are directly genotyped; the remaining 208 must be obtained via genome-wide reference panel imputation — further demonstrating that these loci are broadly distributed across the genome rather than concentrated in coding regions. To fully obtain PRS313 via sequencing, one would need either a specially designed **custom targeted panel** covering the 313 loci and their flanking regions, or **whole-genome sequencing (WGS)**; standard WES cannot accomplish this.

**Conclusion**:

- **Common-variant PRS calculation**: prioritize cohorts providing WGS data (e.g., the WGS subset within All of Us v6), or array genotyping + imputation data; do not use TCGA or standard WES data for direct calculation
- **Rare-variant burden calculation**: WES is naturally suited to this analysis, since WES itself is exome-capture data, and rare coding variants are exactly its area of strength. Existing studies have constructed rare-variant PRS (rvPRS) from UK Biobank WES and imputed genotype (IMP) data, finding that single-SNP-level rvPRS outperforms gene-burden models, and that IMP-derived rvPRS generally outperforms WES-derived models
- TCGA serves only as a supplementary cross-validation resource for rare pathogenic variants (e.g., existing studies have intersected TCGA rare pathogenic variants with All of Us data), not as a data source for PRS calculation or primary analysis

### VI. Analytical Methods

#### 6.1 Aim 1: Rare Variant Burden Test Framework (v4 update)

- **MAF threshold definition**: rare variants are defined as **MAF < 0.1%** for the primary analysis, with 0.1%–1% as the sensitivity-analysis range; this complements the common-variant PRS (typically MAF > 1%), avoiding overlap between the two definitions

- Burden_rare primary analysis method (finalized — no longer a dual-track gene-based burden + rvPRS approach)

  :

  1. Intragenic: only LOFTEE high-confidence PTVs are used for qualifying-variant filtering, producing a binary carrier indicator (0/1), with no continuous dosage weighting
  2. Intergenic: weighted summation using externally published PTV log(OR) values from Dorling et al. (population-based), $\text{Burden_rare}_i = \sum_g \log(\text{OR}*g^{\text{PTV}}) \times \text{carrier}*{i,g}$
  3. The choice of "external OR + binary carrier" over "within-cohort estimated weights + continuous rvPRS" is primarily to **avoid circularity** — if weight estimation and joint-model validation use the same sample, the interaction effect size would be systematically inflated (option comparison detailed in `plans/aim_1.md` Section IV.1)

- **SAIGE-GENE+** (across multiple functional masks: PTV/missense) is still run in parallel, but its role has been adjusted to gene discovery and sensitivity analysis, rather than directly feeding the primary Burden_rare score in the Step 5 joint model

- Missense weighting (REVEL/CADD/AlphaMissense), the STAAR-weighted version, and within-cohort weights estimated via the discovery/target split — all remain as **sensitivity analyses**, not part of the primary analysis

- Full technical prerequisites, CHEK2/panel-gene LD overlap screening, and carrier-sparsity pre-check gating are detailed in `plans/aim_1.md`

#### 6.2 Aim 2: PRS Stratification + Rare Burden Association Analysis

1. Calculate common-variant PRS using All of Us / UKBB data and group by quartile
2. Compare rare-variant burden differences between groups (whether patients with low PRS have higher rare burden)
3. Focus on two extreme subgroups:
   - High PRS + no cancer: testing the "protective/low somatic burden" hypothesis
   - Low PRS + cancer diagnosed: testing the "rare-variant-dominant/somatic-dominant" hypothesis

#### 6.3 Aim 3: Germline-Somatic Interaction Analysis

- Use somatic mutation signature as an intermediate phenotype (rather than only a binary cancer-status outcome), testing whether germline burden (common+rare) can predict an individual's somatic mutation characteristics
- The specific definition of the somatic signature is to be determined through literature review; candidate directions include:
  - The COSMIC mutational signature framework (SBS signatures, etc., describing patterns of somatic mutation — currently the most mainstream)
  - Tumor Mutational Burden (TMB): a simpler quantitative measure of total burden
  - Somatic mutation burden specific to particular DNA repair pathways
  - The three definitions each emphasize different aspects (pattern vs. total quantity vs. pathway specificity); further targeted literature review is needed to determine which best fits this project's hypotheses
- **Whether to include environmental exposure information (e.g., smoking history) in the interaction model requires a data feasibility assessment first**: both All of Us and UKBB contain some self-reported environmental exposure data, but completeness and accessibility need specific verification. The plan is to first complete the "germline burden × somatic signature" two-way association analysis as the primary line of work, with environmental exposure as a future extension, to be added as a three-way interaction depending on the feasibility assessment results

#### 6.4 Aim 4: Age-Interaction Analysis

- Stratified comparison: rare burden, common PRS, and somatic burden between early- vs. late-onset populations
- Correlation analysis: continuous association between burden score and age of onset
- The somatic burden definition shares the same standard as Aim 3, to avoid inconsistency

### VII. Anticipated Challenges and Approaches

| Challenge                                                    | Approach                                                     |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| Genetic contribution effect sizes are weak, requiring a large sample | Prioritize securing All of Us access (existing PRS methodological precedent) + UKBB access |
| Common-variant PRS calculation has high data-type requirements (requires WGS or array+imputation; standard WES is unsuitable) | Prioritize the WGS subset within All of Us or UKBB's array+imputation data; use WES only for rare-variant burden analysis |
| Definition of somatic mutation signature/burden not yet finalized | Requires dedicated literature review (COSMIC signature / TMB / pathway-specific), to be standardized before Aim 3 begins |
| Feasibility of environmental exposure data unclear           | Requires a data availability assessment first, before deciding whether to include it in the interaction model |
| Many similar analyses may already be covered in existing literature | Requires systematic literature search to clarify this project's differentiated contribution regarding "germline-somatic interaction" and "rare burden comparison under PRS stratification" |

### 八、变量口径备忘：三个age变量的定义与使用场景

本项目person-level表型表（`data/04_phenotype_diagnostics_person_level.csv.gz`，由`scripts/04_phenotype_diagnostics.py`生成）里同时存了三个含义不同的age变量，**不可互相替代**，混用会导致变量口径错配。定义与来源见下（详细讨论见`plans/aim_1.md` Step 0）：

| 变量 | 定义 | 计算方式 | 覆盖人群 |
| --- | --- | --- | --- |
| `age_at_diagnosis` | 患病年龄（最早乳腺癌诊断年龄） | `(首次乳腺癌诊断日期 first_bc_dx_date − birth_date).days / 365.25` | **仅case**，control为null |
| `age_at_enrollment` | 入组AoU项目年龄，**不是**患病年龄 | 直接取AoU官方预计算字段`cb_search_person.age_at_consent`（参与者签署primary consent时的年龄），未查到可靠的原始同意书日期字段，故不自行反推 | 全体 |
| `age_at_last_obs` | 末次观察年龄（随访终点/删失年龄） | 见下方"`age_at_last_obs`具体算法" | 全体 |

**`age_at_last_obs`具体算法**（`04_phenotype_diagnostics.py`第227-300行）：

1. **EHR最后记录日期** `ehr_last_date`：对condition_occurrence（`condition_start_date`）、procedure_occurrence（`procedure_date`）、drug_exposure（`drug_exposure_start_date`）、measurement（`measurement_date`）、visit_occurrence（`visit_start_date`）五张**临床事实表**取并集后的`MAX(event_date)`。**不含observation表**——该表里绝大多数记录是PPI问卷（Consent/Basics/健康史等自报问卷）时间戳，混入会让"最后一次临床记录"系统性偏向"最后一次填问卷"，污染这个变量的含义。
2. **死亡日期** `death_date`：`COALESCE(aou_death.death_date WHERE primary_death_record=TRUE, death.death_date)`。
3. **`last_obs_date`** = `ehr_last_date`与`death_date`两者取更晚者（`max(skipna=True)`，任一缺失则取另一个）。
4. **`age_at_last_obs`** = `(last_obs_date − birth_date).days / 365.25`。

**使用场景（不可混用）**：

| 对比类型 | 用哪个age变量 | 原因 |
| --- | --- | --- |
| **case组内部**对比（如按发病年龄分组看PRS差异，对应Aim 4 H3假说） | `age_at_diagnosis` | 这就是"患病时候的年龄"本身；不需要control，control没有这个量 |
| **case vs control两组**间的年龄比较，或联合模型（Aim 1/2公式3）里的age协变量 | `age_at_last_obs` | 两组含义对称——对case是随访终点年龄，对control是"被观察了这么久仍未患癌"的年龄；不会把"是否患癌"这个结局信息倒灌进协变量。已实测排除"case确诊后EHR停滞导致该变量趋近age_at_diagnosis"的陷阱（case组`age_at_last_obs − age_at_diagnosis`中位数7.07年，见`results/table/04_phenotype_case_last_obs_minus_diagnosis_gap.csv`） |
| Table 1描述性统计、招募批次/时期相关的协变量校正 | `age_at_enrollment` | 反映的是"什么时候入组AoU"，与患病与否、随访时长无关 |

**易错点**：`age_at_enrollment`（入组AoU年龄）**不等于**、也不能替代`age_at_diagnosis`（患病年龄）——两者的时间锚点完全独立，一个人可能确诊多年后才入组AoU，也可能入组多年后才确诊。

**`age_at_diagnosis`<18有效性校验（v7.3新增，prepare阶段已做）**：All of Us要求参与者入组时年满18岁（"were at least 18 years old"才能签署consent），女性乳腺癌确诊年龄<18在生物学上极为罕见，因此`age_at_diagnosis`<18的记录判定为诊断码/日期录入误差（如既往史码误判为新发诊断），而非真实儿童期患病。做法参考文献先例——Cheng et al. "Assessing the Data Quality Dimensions of Partial vs Complete Mastectomy Phenotypes in the All of Us Research Program", *JMIR Cancer* 2025;11:e59298（该文将其乳腺手术表型明确限定为"at least 18 years or older at the time of the procedure"，其Table 2也把年龄异常值统一并入"≥80 or <18"这一并列低质量分箱处理）。本队列本轮实测命中4人（全队列9,281 case中），已在`scripts/04_phenotype_diagnostics.py`里对这4人**置空`age_at_diagnosis`并标记`flag_age_at_diagnosis_invalid_under18`=1**（不改变其`is_case`状态——诊断码本身仍视为有效，只是"哪个年龄确诊"这个衍生量不可信），下游按`age_at_diagnosis`分组/回归的分析（Aim 4、case内PRS年龄分层等）会通过标准的缺失值处理自动排除这4人。

------

# Literature Review Summary and Subsequent Analysis Plan

## I. Review of Core Findings from Two Papers

### Paper 1: "Distinct genetic architecture in the tails of complex traits" (Nature, 2026)

**Core method**: The POPout test — reverses the traditional PRS analysis logic, testing "whether PRS varies linearly with phenotype value" rather than "whether phenotype varies with PRS." Under the null hypothesis that common variants dominate, the two should be linear; if the tails of the phenotype distribution show PRS "regression to the mean," this indicates a rare-variant-driven mechanism at that point.

**Core findings**:

- Of 74 continuous traits, 68 showed a significant POPout effect in at least one tail — systematic failure of common PRS at extreme phenotype values is a widespread phenomenon
- After adding rare variants, tail OR improved by an average of 71.7%, while overall R² improved by only 11.6% — **rare variants contribute limited explanatory power overall, but contribute substantially to discriminating the most extreme/highest-risk individuals**
- Forward simulation confirms that **stabilizing selection** drives this phenomenon: extreme phenotypes deviating from the population optimum can only breach the constraints of selective pressure via rare large-effect variants
- **Important limitation (the key issue we identified)**: the POPout method can only be applied to **continuous phenotypes** and cannot directly analyze binary disease outcomes, because its testing logic depends on the continuous ranking concept of "phenotype percentile," and its underlying generative model assumes a bell-shaped distribution structure ("common-dominated middle, symmetric rare-dominated tails on both ends") — this is conceptually mismatched with a conditional, monotonic time-to-event variable like "age of onset"

### Paper 2: "Genetic determinants and genomic consequences of non-leukemogenic somatic point mutations" (Nature Communications, 2025)

**Core method**: GEM (Genomic and Epigenomic Mutation rate) — uses **age as an external weak-supervision signal**, and through a Bayesian graphical model, distinguishes true somatic mutations from sequencing artifacts using a single blood WGS sample, constructing a continuous somatic mutation burden metric.

**Core findings**:

- Using GEM as the phenotype for GWAS identified germline determinants (TERT, TCL1A, NRIP1, etc.), which influence the rate of somatic mutation accumulation indirectly by **regulating hematopoietic stem cell self-renewal capacity** — providing a concrete example of "how common variants influence somatic mutation burden through an indirect biological pathway"
- Rare-variant association analysis (RVAS) was conducted using both STAAR (coding regions) and SCANG (non-coding regions), forming a reusable template for "rare+common dual-track" analysis of somatic mutation burden
- When analyzing clinical outcomes (PAD, stroke, coronary heart disease), **Cox proportional hazards regression** was used, with age included as a covariate via splines — confirming that event-based outcomes with time-to-onset should use a survival-analysis framework, rather than a POPout-style tail test

------

## II. Key Methodological Issues Identified in Our Discussion

### Issue 1: How Should Common+Rare Be "Combined" — a Clear Statistical Method Is Needed

Discussion identified three viable options:

1. **Joint regression + interaction term** (recommended as the starting point): `logit(cancer) ~ PRS_common + Burden_rare + PRS_common × Burden_rare + covariates`
2. **Separate modeling followed by comparison of incremental R²** (PRS-only model vs. GBS-only model vs. joint PRS+GBS model)
3. **Machine learning ensemble** (e.g., the RICE framework, an advanced option for large sample sizes)

### Issue 2: Why the Interaction Term PRS_common × Burden_rare Is Needed

The interaction term directly tests the core of your H1 hypothesis — "whether the effect size of rare variants differs depending on the background level of common PRS":

- Significantly positive interaction → synergistic amplification (supports the H3 "mixture" hypothesis)
- Significantly negative interaction → mutual substitution (supports the H4 "low PRS but cancer driven mainly by rare variants" hypothesis)
- Non-significant → the two are simply additive, itself a meaningful null result

### Issue 3: Data Source Matching — Common PRS and Rare Burden Must Come From the Same Cohort

Two viable paths were identified:

- **All of Us WGS data**: a single data source that supports both common PRS (genome-wide coverage) and rare burden (naturally includes coding-region variants) — the cleanest option
- **TCGA's WES + SNP6.0 array combination**: WES for rare burden, SNP6.0 imputation for common PRS; can serve as a supplementary/cross-validation cohort (but lacks healthy controls, and the imputation pipeline adds uncertainty)

### Issue 4 (the most important correction from this round): The POPout "Tail" Framework Does Not Apply to Age-of-Onset Analysis

The key concern you raised — POPout assumes a model of "a continuous trait shared across the whole population + symmetric tails + stabilizing selection," whereas age of onset is:

- A conditional variable measurable only **within the cancer-affected population**, not a continuous trait shared across the entire population
- Your hypothesis is a **monotonic dose-response relationship spanning the entire burden distribution** (higher burden, earlier onset), not "only people at the extreme ends of the distribution have a different mechanism"
- "Early- vs. late-onset" is a binary split at a clinically meaningful age cutoff; the late-onset group is typically the majority, not a statistical "tail"
- There is no analogous evolutionary-pressure framework of "stabilizing selection/optimal age of onset" — it is closer to a **threshold-trigger model**

**Conclusion**: age-of-onset analysis should adopt a **Cox survival-analysis framework**, rather than a POPout-style tail-deviation test — this is also confirmed by the clinical-outcome analysis method used in Paper 2.

------

## III. Summary of Directly Applicable Analytical Frameworks

| Framework Source                                             | Applicable Scenario                                          | Corresponding Aim                                            |
| ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ |
| POPout test                                                  | Whether common PRS deviates from linear expectation in the extreme-value populations of a continuous intermediate phenotype (e.g., a continuous somatic burden score) | Can be used for intermediate-phenotype analysis in Aim 2, or as an auxiliary test for somatic burden in Aim 3 |
| Paper 1's "common+rare joint PRS" construction method        | MAF stratification (common>1%, moderately rare 0.1%-1%, very rare 0.01%-0.1%) + summing across four burden-variant-risk-score tiers | **Ultimately not adopted** — after discussion, Aim 1 instead uses a joint regression of the two components, PRS_common + Burden_rare, with an interaction term; rationale detailed under "Section II: Core Hypothesis" above |
| Cox proportional hazards regression (age of onset as time-to-event) | Analyzing the monotonic dose-response relationship of "variant burden → onset risk/timing" | **Core statistical method for Aim 4** (replacing the originally envisioned POPout-style tail analysis) |
| GEM framework (age as a weak-supervision signal to distinguish true mutations from artifacts) | Extracting a reliable somatic mutation burden estimator from WGS/WES | Reference for the quantitative method of somatic mutation burden in Aim 3 |
| Dual-track STAAR + SCANG rare-variant association analysis   | STAAR for coding-region rare variants, SCANG for non-coding-region scanning | Methodological supplement for the rare-variant burden test in Aim 3/Aim 1 |
| Tissue-specific factor pathway (e.g., hematopoietic stem cell self-renewal) | How germline common variants influence the rate of somatic mutation accumulation through indirect biological pathways | Provides a biological-mechanism analogy for the germline-somatic association in Aim 3 (note: differences between blood and breast tissue must be considered) |

------

## IV. Subsequent Analysis Plan (Revised)

### Phase 1: Data and Method Preparation (Highest Priority)

1. ~~Confirm the specific sample size of the All of Us WGS white female breast cancer sub-cohort (case/control)~~ **Completed**: 174,381 individuals (9,281 cases / 165,100 controls)
2. ~~Confirm whether the existing Pusztai/Gerstein WES data has accompanying array genotyping data~~ **Finalized**: Pusztai/Gerstein (i.e., the curated TCGA WES data subset) has been downgraded to an independent validation cohort and does not block Aim 1 initiation; UKBB access is expected around 2026-09 and is not currently a blocker
3. ~~Research and determine the MAF threshold and aggregation method for rare-variant burden~~ **Finalized**: MAF<0.1% for primary analysis (0.1%-1% for sensitivity); the aggregation method is **not** a dual-track gene-based+rvPRS approach, but rather summing PTV binary carrier status × Dorling's external log(OR) as the primary analysis, with SAIGE-GENE+ (multi-mask) repurposed for gene discovery/sensitivity analysis — see `plans/aim_1.md` for details
4. Research and determine the quantitative method for **somatic mutation signature** (COSMIC signature / TMB / pathway-specific / a custom continuous burden metric modeled on the GEM approach), standardizing this across Aim 3 and Aim 4

### Phase 2: Aim 1 — Common+Rare Combined Burden Model (v4 execution plan, detailed in `plans/aim_1.md`)

1. Cohort finalized: the All of Us WGS white female sub-cohort of 174,381 individuals (9,281 cases / 165,100 controls). Execution is divided into 5 steps: Step 1 cohort QC → Step 2 PRS307 calculation (liftover hg19→hg38, excluding 2 CHEK2-linked SNPs to avoid double-counting with Burden_rare) → Step 3 qualifying-variant definition + SAIGE-GENE+ association testing → Step 3.5 carrier-sparsity/power pre-check (gating threshold to be populated with empirical data) → Step 4 Burden_rare construction (PTV binary carrier status × summed Dorling external log(OR), avoiding the circularity risk introduced by within-cohort estimated weights)
2. Use a joint regression model (`logit(cancer) ~ PRS_common + Burden_rare + PRS_common × Burden_rare + covariates`) to test the interaction effect — **not adopting** Paper 1's approach of summing common+rare across four tiers into a single score; rationale detailed under "Section II: Core Hypothesis"
3. Supplement with separate modeling + incremental R² comparison (PRS-only model vs. Burden-only model vs. joint model); the held-out test set must be jointly stratified by case/control × carrier status, to avoid carriers being overly diluted within any given fold

### Phase 3: Aim 2 — PRS-Stratified Rare Variant Association Analysis

1. Group by PRS quartile, compare rare-burden differences between groups (as a non-parametric confirmation of the interaction term)
2. Focus analysis on the two extreme subgroups: high PRS without cancer / low PRS with cancer
3. **Optional**: apply a POPout-style test to some continuous intermediate phenotype (e.g., the rare burden score itself) at the high and low ends of PRS, as a supplementary visualization of the interaction effect

### Phase 4: Aim 3 — Germline-Somatic Interaction Analysis

1. Following the GEM approach, construct a continuous quantitative metric for somatic mutation burden (if data conditions allow, attempt weak-supervision calibration using some reliable external anchor to reduce interference from false-positive mutations)
2. Test the association between germline burden (common+rare) and somatic burden/signature, using the regression model: `Somatic_burden ~ PRS_common + Burden_rare + covariates`
3. Whether to include environmental exposure (e.g., smoking) will depend on the results of the preliminary feasibility assessment
4. Drawing on Paper 2's approach of "common variants influencing somatic mutation accumulation through an indirect biological pathway (e.g., stem cell self-renewal)," seek a possible biological mechanistic explanation for any discovered germline-somatic association

### Phase 5: Aim 4 — Analysis of Drivers of Age of Onset (Methodology Revised)

1. Switch to a Cox proportional hazards regression framework, with age of onset as the time-to-event outcome:

   ```
   hazard(onset) ~ PRS_common + Burden_rare + PRS_common × Burden_rare + covariates
   ```

2. Test the monotonic dose-response relationship between burden (joint common+rare) and onset risk/timing, rather than applying POPout's tail-deviation logic

3. Supplementary analysis: split by early-/late-onset clinical cutoff (e.g., age 45), and compare rare burden, common PRS, and somatic burden differences between the two groups, as an intuitive visual supplement to the Cox regression results

### Phase 6: Writing and Communication

1. Compile the above methodological revisions (particularly Aim 4's shift from POPout-style tail testing to Cox survival analysis) into a clear methods draft
2. Prepare a list of questions to discuss with the advisor: how to justify the ancestry/sex cohort restrictions, how to phrase in the methods section the limitation that TCGA (Pusztai/Gerstein), used as an independent validation cohort, consists predominantly of tumor samples and lacks healthy controls, and the final choice of somatic signature definition
3. Cite Papers 1 and 2 as core methodological references for Aim 1 and Aim 3 in the proposal, clearly articulating this project's differentiated contribution relative to both papers: **Paper 1 stops at continuous traits, Paper 2 stops at blood clonal hematopoiesis; this project is the first to systematically extend the "common+rare combined burden" framework to the clinical outcomes of solid-tumor (breast cancer) disease status and age of onset, while integrating germline-somatic interaction analysis**

# PROJECT_STATUS — 当前分析基线

核对日期：2026-10-15。依据现存配置、汇总输出及稿件修订记录；本次没有重跑分析，也没有修改科学定义。

## 使用规则

这四份文件共同维护 single source of truth：本文件记录当前定义；[RESULT_REGISTRY.md](RESULT_REGISTRY.md)记录数值及出处；[ANALYSIS_CHANGELOG.md](ANALYSIS_CHANGELOG.md)记录有证据的决策变化；[NEXT_STEPS.md](NEXT_STEPS.md)记录待办。后续 agent 开始分析前先读四份文件。改变定义时同步本文件和 changelog；改变结果时更新 registry 的版本、模型、来源和校验值，再核对稿件；不得用文件修改时间推断科学决策时间。

“Primary / 当前主稿”表示当前稿件的主分析位置，不表示前瞻预注册或最终冻结。13-gene C 扩展的研究性质仍为 **post hoc exploratory extension**，变异人工裁定尚未最终完成。不得仅凭脚本编号、最新时间戳或较好的结果把 exploratory 提升为 primary。

## 当前论文与稿件

研究 AoU 遗传推断欧洲祖源女性中，PRS311 与罕见变异负荷对乳腺癌病例状态及病例中记录诊断年龄的独立、联合贡献。病例对照预测、病例中诊断年龄关联及绝对风险是不同 estimand；当前结果不提供个体绝对风险校准。

当前可验证稿件对为 [main revised 1008](manuscript/medrxiv_v1_main_revised_1008.docx) 与 [supplement 1008](manuscript/medrxiv_v1_supplement_1008.docx)，依据 [MANUSCRIPT_VERSIONS](manuscript/MANUSCRIPT_VERSIONS.md) 和 [修订验证](manuscript/revision_1006_minimal_1008/validation.json)。`manuscript/README.md` 顶部另称 `_13gene_revised_1008` 为最新；该主稿文件本 checkout 不存在，版本指向冲突尚待统一。

## 当前 cohort

- QC 后描述性队列：143,255 人 = 8,528 cases + 134,727 controls（主稿 Table 1 / Figure 1）。
- 固定病例对照模型队列：141,204 人 = 8,528 cases + 132,676 controls。两组 control 数对应不同阶段，不能混用。
- 固定年龄模型：7,996 cases；排除诊断为首条 EHR 记录者，并沿用已冻结年龄分析资格。
- 女性 sex-at-birth、遗传 EUR posterior >0.75；亲缘去除沿用 >0.1 的贪心程序；missingness 上限 0.03 基于 chr21 1,500-locus proxy，不能称为全基因组 missingness。
- 固定成员、数据版本、QC 细节以 [configuration](config/rowlands_13gene_extension.yaml)、[01](scripts/01_study_cohort.py)、[02](scripts/02_sample_qc.py)、[04](scripts/04_gene_selection.py) 为依据。完整病例模型需要 age at last EHR observation；此次没有重新核对逐人排除原因。

## PRS

PGS000004 / BCAC 313-variant 名义定义，GRCh38、染色体/位置/等位基因匹配及频率消歧。当前稿件记载 304/313 matched：281 direct、23 AF resolved、1 ambiguous、8 absent。去除两个 CHEK2 c.1100delC tagging loci 后，名为 PRS311 的分数实际包含 **302 retained scored loci**。使用已发表效应权重、现有缺失基因型填补程序；直接匹配限定和无均值填补为敏感性分析。匹配数量已核对当前稿件，原始 matching audit 本 checkout 未恢复，不能称为本次重算验证。实现见 [03](scripts/03_common_variants_prs.py)。

## Burden：当前主稿 C

13 genes：BRCA1, BRCA2, PALB2, CHEK2, ATM, BARD1, RAD51C, RAD51D, TP53, PTEN, NF1, STK11, CDH1。

PTV：stop_gained / frameshift / splice donor / splice acceptor；同基因保留转录本满足 LOFTEE HC；没有另加 final-exon 排除。Missense：符合 germline P/LP、主定义 ClinVar >=2 stars，并满足现有疾病/遗传模式及冲突排除规则；ClinVar release 2026-09-13。AoU EUR AF 严格 <0.01；新定义 unknown/disagreeing AF 暂不纳入。site QC、FT PASS 或 '.'、GQ >=20；het VAF 0.3–0.8（含边界），hom-alt VAF >=0.8；VAF 分母为全部 AD 之和。主定义没有额外 DP>=10 门槛，不能把候选审计门槛写成已采用标准。

每人每基因有任一合格变异记一次，基因内不重复累加，跨基因加和 ln(external OR)。权重来自 Rowlands 2024 Table 2 的 PTV+non-PTV PV combined OR；参考定义与项目 ClinVar 日期及 missense 范围不完全一致，属于近似外部加权而非精确复现。实际权重见 [gene_OR_weights](manuscript/datas/burden/gene_OR_weights.csv)。实现 [05](scripts/05_variant_annotation.py)、[06](scripts/06_variant_genotype_qc.py)、[07](scripts/07_burden_score.py)。

保存分数与模型标准化须区分：保存 score 的全 cohort reference 使用 ddof=0；case 模型全固定 cohort ddof=1；age 模型在年龄病例中 ddof=1；CV 仅 training fold 学习缩放（ddof=1）。不得为了统一风格改变这些定义。

## 当前主模型

- CC-JOINT：unweighted binomial GLM MLE，HC1 SE；`is_case ~ PRS311_z + burden_z + PRS311_z:burden_z + age_at_last_obs + PC1..PC10 + sequencing-site dummies + sequencer dummies`。两个主效应均以另一个标准化分数为零为条件。
- AGE-JOINT：HC1 OLS；`age_at_diagnosis ~ PRS311_z + burden_z + PC1..PC10`。无交互项，不加当前年龄。
- CV-4：covariates-only、PRS+covariates、burden+covariates、joint（含交互）；10 repeats × 5 folds，seed 20260818，原 case-status × legacy carrier 分层及原配对 folds，training-only scaling。mean/range 是重复 CV 描述，DeLong CI/P 对应 repeat 1，不能附会为 mean AUC 的 CI。
- 生成源：[08](scripts/08_case_control_models.py)拟合 CC/CV；[09](scripts/09_age_analysis.py)拟合 age 并汇总旧输出；[13 discrimination](scripts/13_figures.py)保留 PRS/joint 预测并补 covariates/burden-only 模型。

## 分析分层

| 分析 | 当前角色 | 限制 |
| --- | --- | --- |
| C 13-gene PTV + eligible missense | 当前主稿主分析 | post-hoc extension；final adjudication pending |
| C_without_STK11, C_PTV_only, C_TP53_NF1_strict_VAF, C_ClinVar_min1star, C_without_CDH1 | 主定义敏感性 | 在 C association 结果前锁定，不等于整个研究前瞻预注册 |
| A / A-prime / QC bridge / B | 历史及桥接比较 | 不能替代 C 的当前数值 |
| Gene × PRS、影响点/模型规格稳健性 | Post hoc；S7/S8 | BH 校正基因交互；不得称主模型 |
| Binary carrier age slopes、EHR checks、frailty、review Cox/family analyses | Exploratory / reviewer analyses | 不能因脚本有编号或稿件有讨论即晋升 primary |
| 旧 9-gene / 10-gene 与旧 supplement results | Historical / legacy | 需保留版本身份，不静默覆盖 |

Sensitivity score 计划 [10](scripts/10_sensitivity_analysis.py)需在 08 前准备，编号不是完整依赖图；post-hoc [11](scripts/11_interaction_analysis.py)另行运行。

## 当前 tables / figures 与未对齐内容

- Table 1：QC 后描述性 cohort；Table 2：C joint case model 与 discrimination；Table 3：C age model。S7/S8：post-hoc interaction 及 robustness。
- Figure 1：cohort flow；2：score distributions；3：four-model discrimination；4：PRS decile × carrier descriptive joint distribution；5：age / carrier patterns；6：主分析及敏感性估计。外部图文件映射见 [figure mapping](manuscript/figure_captions/figure_file_mapping.csv)。Figure 4 expectation model 为描述性无交互模型，不能误称完整 CC-JOINT。
- 1008 最小修订保留原 1006 Figures 1–6 及 captions、Table 1；验证的是保留一致性，并非本次逐图数据复核。原 supplementary S2–S6 保留历史 panel 内容；文前作者 notes 也可能含旧探索数值。[修订记录](manuscript/revision_1006_minimal_1008/revision_log.md)明确说明范围。
- `main_text.md` 含旧新数值串连（如 interaction P 同时有旧值和新值）。可能涉及 Word tracked changes / 抽取方式；不能仅凭该文本判定最终显示错误，需检查 Word revision markup 与 accepted view。

## 证据边界与执行状态

本 checkout 没有 `results/` 原输出目录。21 份 manuscript aggregate 文件与 [manifest](manuscript/datas/manifest.json) 的 SHA-256 完全相符；这是副本完整性核对，不是原始数据/模型重跑认证。交互副本另外与本地 S7/S8 比较，行记录并非完全一致，尚需区分格式差异与数值差异；未读取或复制 participant-level 表。

2026-10-15 可读性重构：37 synthetic checks；13 model CSV byte-identical、10 CV arrays identical；真实数据及稿件未重跑。见 [validation](docs/readability_validation.md)。`07_burden_score.py` 原有 undefined `report` 的摘要写出错误仍待单独处理；未静默修复。科学定义仍 provisional；未来重跑需新 source/input lock。

# PROJECT_STATUS — 当前分析基线

2026-10-16 / AGE-FOUR-20261016-01：新增真实7,996病例四组描述，PRS按全部病例共享中位数分高/低，早发<45。carrier/PRS低n323、中位53.98岁、早发90/323=27.86%；carrier/PRS高n220、中位55.39岁、早发25.00%；non-carrier/PRS高n3778、中位60.90岁、早发10.53%。[报告](docs/age_interaction_runs/AGE-FOUR-20261016-01/REPORT.md)。未调整描述不提供交互检验或发病风险；当前模型基线不变，不supersede既有结果，稿件未修改。

## 2026-10-16增补：AGE-INT-20261016-01已完成

- 作者已授权并实际完成连续13-gene C burden×PRS311病例年龄交互；沿用相同7,996病例和原分数、病例内ddof1、PC1–10，不增加incident限制。HC3主要，HC1对照。作者陈述原研究问题中遗漏的步骤；不据此宣称预注册或数据盲。
- 新增主文分析结果：交互beta=+0.1375698983年/(SD_P×SD_B)，HC3 95%CI[−0.1478221129,+0.4229619095]，P=0.3447723108；ΔR²=0.0128442个百分点。无充分交互证据，不等于证明为零。数值/哈希/selector见RESULT_REGISTRY的AGE-INT条目。
- AoU live常数只读访问、真实冻结输入复核及独立矩阵公式核验通过；合成测试与真实拟合分别登记。参与者行/个人分数未写出，新内部日志/快照仅保存在ignored data run；未stage/commit/push。
- **替代关系：无。** 既有加性AGE-JOINT/Table3、binary carrier年龄交互和病例状态/CV基线均保留。新增M1为额外主文分析，稿件文字/表图尚未写入；不从新条件主效应替换旧加性系数。
- [完整运行报告](docs/age_interaction_runs/AGE-INT-20261016-01/REPORT.md)、[运行索引](ANALYSIS_RUN_INDEX.md)。最终变异人工裁定、从原始基因型重建评分、全部13个外部OR原表逐值认证未完成。没有新增Cox/frailty结果；既有版本冲突保留。

以下保留此前基线登记及其当时证据边界；涉及“未真实重跑”或“尚未年龄复核”的旧状态，以本增补及AGE-AUDIT/AGE-INT记录覆盖的明确范围为准，不能据此扩大为全稿认证。

## 一分钟交接（2026-10-16；工作 DOC-20261016-02）

- **已登记分析基线（待与1012逐项对齐）**：PRS311（历史记录实际302 scored loci）+13-gene C；CC-JOINT、AGE-JOINT、CV-4，定义见下文，数值以 RESULT_REGISTRY 为准。C为此前登记的主稿分析，仍属post hoc、provisional；此次版本确认不证明1012全部数值相同。
- **主要结论**：现存C汇总中CC的PRS/burden OR约1.554/1.269，乘积OR约0.972；joint AUC约0.7815、相对PRS增量约0.01328；病例诊断年龄斜率约−0.433/−1.820年每SD。这些是汇总副本支持的关联，非独立重算或绝对风险结论。
- **当前稿件**：作者于2026-10-16明确确认1012为最新版本：main_1012 + supplement_1012；1008转为历史稿件。1012与registry的数值/模型/表图对应仍待审计。
- **最新变化**：本轮根据作者确认更新稿件采用版本；未运行分析，registry数值未改变。全基因组PRS v2是审查草案，PGS005346/PGS005347在本轮未评分、未拟合、未生成新AUC；v2修订v1.1计划，不取代PRS311结果。
- **未决/下一入口**：先处理NEXT_STEPS的VERSION/MARKUP/COVERAGE/FIGURES及原证据恢复，再按全基因组v2的未决项锁定新run。新agent从上述C基线继续；不得从最新文件名选结果。


核对日期：2026-10-15。依据现存配置、汇总输出及稿件修订记录；本次没有重跑分析，也没有修改科学定义。

## 使用规则

这四份文件共同维护 single source of truth：本文件记录当前定义；[RESULT_REGISTRY.md](RESULT_REGISTRY.md)记录数值及出处；[ANALYSIS_CHANGELOG.md](ANALYSIS_CHANGELOG.md)记录有证据的决策变化；[NEXT_STEPS.md](NEXT_STEPS.md)记录待办。后续 agent 开始分析前先读四份文件。改变定义时同步本文件和 changelog；改变结果时更新 registry 的版本、模型、来源和校验值，再核对稿件；不得用文件修改时间推断科学决策时间。

“Primary / 当前主稿”表示当前稿件的主分析位置，不表示前瞻预注册或最终冻结。13-gene C 扩展的研究性质仍为 **post hoc exploratory extension**，变异人工裁定尚未最终完成。不得仅凭脚本编号、最新时间戳或较好的结果把 exploratory 提升为 primary。

## 当前论文与稿件

研究 AoU 遗传推断欧洲祖源女性中，PRS311 与罕见变异负荷对乳腺癌病例状态及病例中记录诊断年龄的独立、联合贡献。病例对照预测、病例中诊断年龄关联及绝对风险是不同 estimand；当前结果不提供个体绝对风险校准。

当前稿件对为 [main 1012](manuscript/medrxiv_v1_main_1012.docx) 与 [supplement 1012](manuscript/medrxiv_v1_supplement_1012.docx)，依据作者2026-10-16本轮确认；文件哈希与取代关系见 [MANUSCRIPT_VERSIONS](manuscript/MANUSCRIPT_VERSIONS.md)。1008及其[修订验证](manuscript/revision_1006_minimal_1008/validation.json)保留为历史证据，不认证1012。下文稿件数值、图表和保留内容的既有核对范围为1008；1012尚待逐项对齐。

## 当前 cohort

- QC 后描述性队列：143,255 人 = 8,528 cases + 134,727 controls（主稿 Table 1 / Figure 1）。
- 固定病例对照模型队列：141,204 人 = 8,528 cases + 132,676 controls。两组 control 数对应不同阶段，不能混用。
- 固定年龄模型：7,996 cases；排除诊断为首条 EHR 记录者，并沿用已冻结年龄分析资格。
- 女性 sex-at-birth、遗传 EUR posterior >0.75；亲缘去除沿用 >0.1 的贪心程序；missingness 上限 0.03 基于 chr21 1,500-locus proxy，不能称为全基因组 missingness。
- 固定成员、数据版本、QC 细节以 [configuration](config/rowlands_13gene_extension.yaml)、[01](scripts/01_study_cohort.py)、[02](scripts/02_sample_qc.py)、[04](scripts/04_gene_selection.py) 为依据。完整病例模型需要 age at last EHR observation；此次没有重新核对逐人排除原因。

## PRS

PGS000004 / BCAC 313-variant 名义定义，GRCh38、染色体/位置/等位基因匹配及频率消歧。当前稿件记载 304/313 matched：281 direct、23 AF resolved、1 ambiguous、8 absent。去除两个 CHEK2 c.1100delC tagging loci 后，名为 PRS311 的分数实际包含 **302 retained scored loci**。使用已发表效应权重、现有缺失基因型填补程序；直接匹配限定和无均值填补为敏感性分析。匹配数量此前已核对1008稿件，原始 matching audit 本 checkout 未恢复，不能称为本次重算验证。实现见 [03](scripts/03_common_variants_prs.py)。

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


## AGE-QUART-20261016-01 / 2026-10-16

Real frozen-input Q1/Q4 descriptives completed (not synthetic-only), shared quartiles in original7,996 cases; selected3,998; early<45 fixed; existing13-gene carrier. Source/command/software/input and output hashes/status: docs/age_interaction_runs/AGE-QUART-20261016-01/provenance.json; [report](docs/age_interaction_runs/AGE-QUART-20261016-01/REPORT.md). No participant records/scores written; controlled logs and snapshot ignored. No regression or significance test; no failure/blocked stage. Supersedes=None; prior median split, additive and continuous-interaction results retained; manuscript unchanged.


## FIGURES-PUBLISH-20261017-01 / 2026-10-17
Author requested publication of the latest figures. Updated root figures/ from the current manuscript figure set, with aggregate tables in tables/ and reproduction scripts in scripts/. Figure04a uses full-cohort (141204) PRS quartiles; Figure04c restores four carrier/PRS combinations using full-cohort median and bottom/top20% panels. Figure07b retains original four attained-age bands and is8in wide. Supplementary Figure07b uses seven age bands; two real Cox fits passed convergence, with35 HR/CI rows checked against stored coefficients and covariance. Red/blue1/blue2/orange series; gray theta0.73 line reuses the existing illustrative refitted scenario. Eight-band trials failed convergence and were not adopted. No participant-level data, individual scores, controlled logs, local archives or unrelated workspace changes are included. Original four-band scientific baseline unchanged; seven-band extension remains post hoc/provisional. Manifest: docs/figure_publication_20261017.json; detailed source-run evidence: docs/age_interaction_runs/. Published script output paths adapt to this repository layout; provenance describes the original controlled workspace. Controlled-data reproduction and embedded-document consistency are not claimed.

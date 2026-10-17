# NEXT_STEPS — 按依赖排序的待办

AGE-FOUR-20261016-01已完成四组真实描述（原7,996例、13基因carrier、共享病例PRS中位数、早发<45），无需再次运行；报告及年龄分布图在docs/age_interaction_runs/AGE-FOUR-20261016-01/。如写入稿件应标明病例内未调整描述、carrier不等于连续高burden、PRS高不是外部高风险分类；不以组内比例代替发病风险或交互检验。旧模型基线不变。

## AGE-INT-20261016-01完成后的下一步

- 已完成：原7,996病例、既有13基因C/PRS311的连续年龄交互真实后台拟合；HC3主要及HC1对照；加性Table3复核、独立数值验证、聚合残差/共同支持图及运行记录。非仅合成测试，无真实分析blocked。报告：docs/age_interaction_runs/AGE-INT-20261016-01/REPORT.md。
- 下一步主文写入：按作者已决定的主文位置，单列连续交互beta/HC3 CI/P及ΔR²；保留已解决的Table3加性模型，不以新条件系数覆盖旧加性系数。当前未修改DOCX或embedded figures。
- 新增结果不取代原AGE-JOINT、binary carrier交互、CC或CV。完整原始基因型/PRS重建、变异最终人工裁定、外部13权重逐值认证仍未核验；本轮按作者决定复用既有分数，不以这些未完成项冒称新分数已认证。
- 可选扩展仍未执行：非线性替代模型、leave-one-gene或影响点删除分析。整体RESET及HC1/HC3对照已做，不能等同完成这些扩展。incident分析本次明确不纳入，不列为新交互交付阻塞项。
- 既有Cox/frailty采用输出冲突、补充图编号及参考文献错号仍按既有审查待办处理，不因新交互完成而关闭。
- 所有新工作在ANALYSIS_RUN_INDEX与ANALYSIS_CHANGELOG交叉登记；受控数据/个人分数保持git之外。

## 2026-10-16追加（稿件版本更新；登记数值不变）

- VERSION：作者已确认1012主稿/补充稿为最新，版本指向已统一；采用身份问题关闭。后续核对1012与1008的内容差异及registry映射。
- GWPRS：以v2为待审计划，先解决其中来源/输入锁、权重及技术覆盖、优先顺序与门槛等未决项；本轮没有新增AoU评分或关联运行。v1.1为历史计划，C仍为已有结果。
- 每次结束按CLAUDE合同追加ANALYSIS_CHANGELOG工作索引；独立核验须有方法和证据。目前哈希完整性核对不得当作独立模型复算。


更新：2026-10-15。先完成结果与稿件对齐，再以独立 run 开始新的科学迭代。完成任务后记录来源并更新其余三份 SSOT 文件；不要把待确认事项当成已决定。

| 优先级 / ID | 下一步 | 完成条件 | 依赖 / 风险 |
| --- | --- | --- | --- |
| 已完成 / VERSION | 作者确认1012为当前主稿/补充稿，入口已统一 | 2026-10-16作者本轮明确确认；文件哈希见MANUSCRIPT_VERSIONS | 1008保留历史；不代表1012数值已核验 |
| P0 / MARKUP | 检查当前 DOCX tracked revisions / accepted view | 旧新 P 值串连是否仅抽取问题有明确结论；若显示错误，列出位置与正确 registry 值后再修稿 | main_text.md 同时出现旧新 P；原始 DOCX 保留 |
| P0 / COVERAGE | 逐条盘点1012 Abstract、Results、Tables、captions、supplement 数值 | 每项 manuscript claim 有 result ID、精确 output selector、模型/cohort/版本和校验值；无来源项明确挂起 | 当前 registry 已登记现存主要系数/CV/交互；S7/S8 稿件表副本与 interaction CSV 行记录并非完全一致，需区分格式差异与数值差异；不是全稿完成认证 |
| P0 / FIGURES | 核对 embedded figures 与当前 aggregate 数据 | Figure 1–6、Table 1、S2–S8 明确当前/历史；登记所有不一致 | 1008 验证保留1006图，不证明图中每个数字已刷新；避免混用143,255和141,204 |
| P0 / ADJUDICATION | 完成现存计划中的 variant manual adjudication / definition freeze | 有独立批准记录、冻结 eligibility/weights/input/run；若改变定义，新版本重跑而不覆盖 C | 当前 `ptv_final_definition_frozen=false`；此次不代替科学裁定 |
| P1 / PROVENANCE | 恢复或定位原始 aggregate outputs / source locks | 21份已校验副本与原结果、input/source lock、运行记录可追溯；恢复 match attrition 与 age-group carrier summaries | 当前 checkout 无 results/；不要通过重新拟合“补回”历史证据 |
| P1 / HISTORY | 补齐最近科学决策时间线 | 有证据的请求/批准日期与定义差异，区分 proposal/implemented | Git 仅到9月28日，不能依据 mtime 补造历史 |
| P1 / BURDEN-BUG | 单独处理07摘要 undefined `report` | 先确认独立 bug 修复范围、保留现有分析数值，再验证写出 | 预存在，已报告；默认 burden 重跑前解决，不混入架构调整 |
| P1 / BASELINE | 确认未来 primary status 的治理规则 | main placement、post-hoc provenance、provisional/frozen 三轴持续记录 | 不把当前主稿 C 改称 prospectively confirmatory；11及review不自动晋升 |
| P2 / RERUN | 如需真正重跑，建立独立 run | 新source/input锁、按依赖运行、与既存输出比较、完整登记数值变化 | 10 score plan先于08；09接08；不可用synthetic checks替代真实输出验证 |
| P2 / LAYOUT | 若仍需要 data/results 目录物理迁移，独立规划 | 全路径、source locks、脚本及稿件引用一致性验证 | 不为美化目录而改变或覆盖已有数据路径 |

本轮已完成：四份文档基线；21个 manifest aggregate hash 核对；主要 C 系数、四模型 AUC、年龄 R²、敏感性/历史系数及 post-hoc interaction 原精度登记。未运行 participant-data 分析，未修改稿件或分析源代码。


## AGE-QUART-20261016-01 / 2026-10-16

Real frozen-input Q1/Q4 descriptives completed (not synthetic-only), shared quartiles in original7,996 cases; selected3,998; early<45 fixed; existing13-gene carrier. Source/command/software/input and output hashes/status: docs/age_interaction_runs/AGE-QUART-20261016-01/provenance.json; [report](docs/age_interaction_runs/AGE-QUART-20261016-01/REPORT.md). No participant records/scores written; controlled logs and snapshot ignored. No regression or significance test; no failure/blocked stage. Supersedes=None; prior median split, additive and continuous-interaction results retained; manuscript unchanged.


## FIGURES-PUBLISH-20261017-01 / 2026-10-17
Author requested publication of the latest figures. Updated root figures/ from the current manuscript figure set, with aggregate tables in tables/ and reproduction scripts in scripts/. Figure04a uses full-cohort (141204) PRS quartiles; Figure04c restores four carrier/PRS combinations using full-cohort median and bottom/top20% panels. Figure07b retains original four attained-age bands and is8in wide. Supplementary Figure07b uses seven age bands; two real Cox fits passed convergence, with35 HR/CI rows checked against stored coefficients and covariance. Red/blue1/blue2/orange series; gray theta0.73 line reuses the existing illustrative refitted scenario. Eight-band trials failed convergence and were not adopted. No participant-level data, individual scores, controlled logs, local archives or unrelated workspace changes are included. Original four-band scientific baseline unchanged; seven-band extension remains post hoc/provisional. Manifest: docs/figure_publication_20261017.json; detailed source-run evidence: docs/age_interaction_runs/. Published script output paths adapt to this repository layout; provenance describes the original controlled workspace. Controlled-data reproduction and embedded-document consistency are not claimed.

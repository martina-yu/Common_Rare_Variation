# NEXT_STEPS — 按依赖排序的待办

更新：2026-10-15。先完成结果与稿件对齐，再以独立 run 开始新的科学迭代。完成任务后记录来源并更新其余三份 SSOT 文件；不要把待确认事项当成已决定。

| 优先级 / ID | 下一步 | 完成条件 | 依赖 / 风险 |
| --- | --- | --- | --- |
| P0 / VERSION | 统一 manuscript/README 与 MANUSCRIPT_VERSIONS 的当前稿件指向 | 作者认可的唯一 main/supplement 对；检查缺失 `_13gene_` 命名文件是否另有来源 | 当前可验证为 main_revised_1008 + supplement_1008；不要自动覆盖文档 |
| P0 / MARKUP | 检查当前 DOCX tracked revisions / accepted view | 旧新 P 值串连是否仅抽取问题有明确结论；若显示错误，列出位置与正确 registry 值后再修稿 | main_text.md 同时出现旧新 P；原始 DOCX 保留 |
| P0 / COVERAGE | 逐条盘点 Abstract、Results、Tables、captions、supplement 数值 | 每项 manuscript claim 有 result ID、精确 output selector、模型/cohort/版本和校验值；无来源项明确挂起 | 当前 registry 已登记现存主要系数/CV/交互；S7/S8 稿件表副本与 interaction CSV 行记录并非完全一致，需区分格式差异与数值差异；不是全稿完成认证 |
| P0 / FIGURES | 核对 embedded figures 与当前 aggregate 数据 | Figure 1–6、Table 1、S2–S8 明确当前/历史；登记所有不一致 | 1008 验证保留1006图，不证明图中每个数字已刷新；避免混用143,255和141,204 |
| P0 / ADJUDICATION | 完成现存计划中的 variant manual adjudication / definition freeze | 有独立批准记录、冻结 eligibility/weights/input/run；若改变定义，新版本重跑而不覆盖 C | 当前 `ptv_final_definition_frozen=false`；此次不代替科学裁定 |
| P1 / PROVENANCE | 恢复或定位原始 aggregate outputs / source locks | 21份已校验副本与原结果、input/source lock、运行记录可追溯；恢复 match attrition 与 age-group carrier summaries | 当前 checkout 无 results/；不要通过重新拟合“补回”历史证据 |
| P1 / HISTORY | 补齐最近科学决策时间线 | 有证据的请求/批准日期与定义差异，区分 proposal/implemented | Git 仅到9月28日，不能依据 mtime 补造历史 |
| P1 / BURDEN-BUG | 单独处理07摘要 undefined `report` | 先确认独立 bug 修复范围、保留现有分析数值，再验证写出 | 预存在，已报告；默认 burden 重跑前解决，不混入架构调整 |
| P1 / BASELINE | 确认未来 primary status 的治理规则 | main placement、post-hoc provenance、provisional/frozen 三轴持续记录 | 不把当前主稿 C 改称 prospectively confirmatory；11及review不自动晋升 |
| P2 / RERUN | 如需真正重跑，建立独立 run | 新source/input锁、按依赖运行、与既存输出比较、完整登记数值变化 | 10 score plan先于08；09接08；不可用synthetic checks替代真实输出验证 |
| P2 / LAYOUT | 若仍需要 data/results 目录物理迁移，独立规划 | 全路径、source locks、脚本及稿件引用一致性验证 | 不为美化目录而改变或覆盖已有数据路径 |

本轮已完成：四份文档基线；21个 manifest aggregate hash 核对；主要 C 系数、四模型 AUC、年龄 R²、敏感性/历史系数及 post-hoc interaction 原精度登记。未运行 participant-data 分析，未修改稿件或分析源代码。

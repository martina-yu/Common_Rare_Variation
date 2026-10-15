# ANALYSIS_CHANGELOG — 科学决策与技术变更分开记录

核对日期：2026-10-15。只记录能由现有证据支持的日期和决策。Git 最新提交停留在 2026-09-28，无法据此重建 10 月所有工作；文件 mtime 不代表决策日期。没有证据的执行/裁定时间不补写。

## 科学定义 / 解释变化

| 日期 / 状态 | 可证实事件 | 对结果的影响 / 边界 | 证据 |
| --- | --- | --- | --- |
| 2026-10-07，proposal | 提议把 C 13-gene PTV + eligible missense 作为主稿负荷；旧 panel 比较移到补充；替换旧 discrimination / age 叙述 | 提案不是新模型执行凭证，也不证明 final variant freeze | [revision outline](manuscript/manuscript_revision_outline_13gene.md) |
| 2026-10-08，implemented manuscript revision | 1006→1008 最小修订，使用已存在 C 汇总结果更新 Methods、Table 2/3 及相关叙述；追加 gene×PRS S7/S8 | 当前稿件主分析位置改变不等于研究前瞻预设；没有凭此证明新一次模型重跑 | [revision log](manuscript/revision_1006_minimal_1008/revision_log.md), [patch](manuscript/revision_1006_minimal_1008/paragraph_and_cell_patch.json) |
| 2026-10-08，scope deliberately retained | Figures 1–6/captions、Table 1 及原 supplement 内容保留；删除/改写缺少当前验证的旧校准、分组等主文叙述 | S2–S6 和作者 notes 仍可能是历史定义；不能混入当前 registry | [validation](manuscript/revision_1006_minimal_1008/validation.json) |
| 决策日期未确认，当前配置事实 | C 主定义及五个 eligibility sensitivity 已运行；敏感性在 C association 结果前锁定；最终人工裁定未冻结 | post-hoc extension / provisional，不能扩大为全研究 prospectively prespecified | [configuration](config/rowlands_13gene_extension.yaml) |
| 2026-10-15，本次证据核对 | 建立四份 SSOT 文档，区分 QC cohort 与 model cohort、main placement 与 post-hoc chronology | 没有改定义、阈值、权重、模型或输出；没有重跑分析 | [PROJECT_STATUS](PROJECT_STATUS.md), [RESULT_REGISTRY](RESULT_REGISTRY.md) |

## 尚不能确认的科学决策历史

13-gene 扩展、各 QC 规则、外部权重、C 敏感性锁定及 reviewer 模型的准确批准/执行时间，不能从现存 Git 完整恢复。旧 plan 中“尚未运行”与当前配置/结果矛盾，应视为历史计划；不能据此否定现存输出，也不能虚构决策日期。后续补入需要日期、原请求/批准、旧定义、新定义、理由、受影响输出及核验来源。

## 后续新增条目模板

`日期 | decision ID | proposed/approved/implemented | 旧定义→新定义 | 理由与批准来源 | cohort/model/result IDs | 新旧 output/run/checksum | manuscript impact | verification`

只有叙述/架构调整应标记 technical/editorial；发现科学 bug 时先报告，不在清理代码过程中默改。

## 技术记录：2026-10-15 readability（保留此前记录）


- Adopted the approved thirteen numbered filenames and project-level tests.
- Restored scientific operations from the pre-modularization source as sequential sections.
- Removed lazy task objects and whole-stage forwarding modules from main analysis scripts.
- Kept scientific formulas, thresholds, ordering, weights, seeds, columns and output paths.
- Restored the case/age execution boundary from the original combined implementation, with exact aggregate JSON handoff and provenance checks.
- Reverted the interrupted, unrelated replacement of the undefined burden report; the pre-existing issue is explicitly documented.
- Retained reusable domain helpers; preserved superseded sources in explicit source snapshots.
- Updated active workflow/import paths and added result/status documentation.

Validation details and every moved/changed source are recorded in
docs/readability_validation.md and docs/readability_file_changes.md.
No participant-data analysis or manuscript document was regenerated.

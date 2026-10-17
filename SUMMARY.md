# 项目现状与阅读导航

更新：2026-10-13。这里描述当前磁盘结构；科学结论应以对应审计报告、配置及锁定输入为依据。

本项目分析 AoU 女性 WGS 队列中的乳腺癌病例状态和病例内首次记录诊断年龄，比较常见变异 PRS 与罕见变异负担。当前脚本说明采用 13 基因 PTV + ClinVar germline P/LP、至少两星 missense 的版本。九基因探索与历史 PTV 分析作为单独分支保留。

完整目录清单和来源证据见 [脚本与结果导航](docs/project_catalog/README.md)。根目录原 [README.md](README.md)包含早期研究背景；当前执行入口以 [scripts/README.md](scripts/README.md)、配置和工作流为准。

| 分支 | 入口 | 结果与说明 |
|---|---|---|
| 队列与预处理 | `01`、`02`、`03` | `data/` 中队列、QC、协变量、PRS；需 Workbench 数据权限 |
| 13 基因主分析 | `04`–`08`；`workflow/rowlands_extension.smk` | `results/table/`、`results/report/`；锁定输入在 `data/29_rowlands_extension/` |
| 图与模型解释 | `09_result_figures.py` 的显式任务 | `results/figures/` 和 `manuscript/figures/`；报告编号沿用历史编号 |
| 诊断年龄与 EHR 探索 | `10_exploratory_age_ehr.py`，必须选任务 | `results/common_rare_age/`、`27`/`28` 系列表与报告 |
| 历史 PTV 比较 | `90_legacy_ptv_baseline.py`，必须选任务 | `12`–`20` 系列历史结果；部分文件也是当前比较的输入 |
| 稿件审计与数据整理 | `11`、`12`、`audit_table1_09b.py` | `results/report/`、`manuscript/datas/` |
| 交互作用及稿件修订 | `59`、`60`、`61` | `59`/`60`/`61` 系列报告、稿件资产与修订 DOCX |
| 历史源码 | `archive/` | 合并前源码和验证快照；不作为默认运行入口 |

建议先读 [脚本说明](scripts/README.md)，再查 [任务清单](docs/project_catalog/tasks.csv)和[结果清单](docs/project_catalog/results.csv)。需要追溯计算时使用[函数索引](scripts/FUNCTION_INDEX.md)；需要了解已有验证边界时使用[整理验收记录](scripts/refactor_validation.md)。

本次完成的是文件组织、静态任务提取与来源审计，没有重跑真实数据分析。来源引用表示定位线索，工作流输出声明表示明确的生成约定，都不能代替运行日志、环境版本及输入锁。

复现仍依赖 Workbench 认证、云端数据、现有本地输入与配置。已有 YAML 环境记录不是完整的平台锁文件。优先核对 [缺失输入与统计](docs/project_catalog/audit.json)和[待补来源](docs/project_catalog/UNRESOLVED.md)，再决定重跑范围。

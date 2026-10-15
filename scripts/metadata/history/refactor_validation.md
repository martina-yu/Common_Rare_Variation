# 具体运算回归主脚本：验收记录

日期：2026-10-06。本记录替代此前“主脚本调度函数模块”的写法验收。

## 最终结构

- 9 个主步骤、2 个可选步骤；44 个运算段落的实际函数体位于对应编号脚本。
- 每个主脚本直接包含读取、筛选、计算、模型或输出操作；计算段落采用明确命名和分隔注释。
- `functions/` 保留 13 个辅助模块（不计包标记），不含执行完整步骤的 `main()`／`run_analysis()`。common variant 辅助库仅保留独立的 effect-allele 方向和 BED 编解码。
- 单基因 QC 的真实计算移入 `06`，并行进程通过同一个 `06 --task gene-calls` 调用；复用 QC 小函数的跨段引用仍指向辅助库。
- `07` 中直接保留 `read_weights()`、`score_panel()` 和计分／比较过程；`08` 展开实际模型矩阵、logistic、CV 和 age OLS。
- 新模型输入锁会同时包含真实分析脚本和辅助函数源码。历史锁及研究结果没有改写。

## 通过的检查

- **19 项既有合成测试通过**：ClinVar 分类／断言、原始 ALT、GQ／VAF、位点 FILTER、病例内 EHR 审核、missense carrier 与 burden。
- **11 项新结构回归检查通过**：运算在主文件内、gene/VAT 操作展开、辅助库无隐藏分析入口、合并后无同名函数覆盖、位置查找、参数恢复、单任务参数约束、主脚本／辅助库 burden 数值一致、BED 编解码／方向一致、源码锁篡改检查、无云端凭证时导入 cohort／PRS 步骤。
- **269 个函数体结构等价**：与 `archive/scripts_process_layout_20261006/scripts/functions/` 中整合前源码逐函数 AST 比较，归一化变量／函数命名、源码溯源路径及增加的源码锁列表后，计算结构无差异。旧 forest 的 ROOT 计算按文件上移一层调整；并行 QC 两个函数单独核对了新的同文件执行入口。
- 科学配置解析后逐字段相同，仅三个源码来源路径更新。
- 原模型的 **9 个锁定输入核验通过**，归档原源码与当前实现哈希均受保护；没有重写旧输入锁。
- 全部 Python 源码编译检查、工作流源码文件路径、44 个运算的位置索引通过。

## 验证边界

没有重新执行真实 AoU 数据提取、完整基因型查询、关联模型或结果生成；上述验证不等于完整 Workbench 流水线重跑。当前环境没有 pytest／snakemake，现有测试通过 unittest 和直接调用合成 `test_*` 完成，未执行 Snakemake dry-run。

主脚本计算写法参考作者提供的 npjYoungVsOld 与 bedeepon 分析／预处理代码，仅调整本项目的组织与可读性，没有移植外部项目的科学公式或过滤规则。此前源码、配置、工作流及调度写法保留为审计快照。

## Documentation cleanup (2026-10-06)

- Numbered scripts now open with concise English step bullets.
- Historical prose was removed from code; cohort details are in README.md.
- Comments and docstrings are English throughout active scripts and helpers.
- Executable AST equality verified for all 25 edited Python files.
- All 11 process-layout regression checks passed; source hashes and 44 calculation links were refreshed.
- Runtime messages, scientific calculations and data outputs were unchanged.
# Manuscript plotting script consolidation (2026-10-06)

The ten newly added scripts `13`–`22` were archived in `archive/scripts_manuscript_plots_20261006/` and their calculation bodies merged into `09_result_figures.py`. Global constants and helper functions receive unique namespaces; each main calculation has a descriptive function name and explicit `--task` entry. The prior default task list is unchanged.

All ten original executable module ASTs match the merged bodies after reversing namespace changes and excluding standalone `__main__` guards. Numerical formulas, thresholds, model inputs and plotting settings are unchanged. The histogram manifest source now points to its active task. The manuscript discrimination workflow uses `09 --task discrimination`; the inventory, function index and layout metadata are updated.

Validation: all 11 existing `scripts.tests.test_pipeline_refactor` tests passed, and `09 --where` successfully located all 13 tasks. No participant-data analyses, figures or DOCX files were regenerated during this reorganization.

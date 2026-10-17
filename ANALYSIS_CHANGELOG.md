# ANALYSIS_CHANGELOG — 科学决策与技术变更分开记录

| 工作ID / 日期 | 工作与完成状态 | 定义 / 证据 / 核验 | 基线影响 / 下一步 |
|---|---|---|---|
| AGE-FOUR-20261016-01 / 2026-10-16 | 用户要求快速四组描述；后台真实7,996病例计算完成，无回归或显著性检验、无失败 | 沿用13基因carrier；PRS按全部年龄病例共享中位数分组；早发<45、年龄分箱预先锁定；脚本scripts/describe_age_four_groups.py；docs/age_interaction_runs/AGE-FOUR-20261016-01/REPORT.md与provenance；父输入哈希前后相符、carrier与B_raw_C>0相符、分组和7996、公开计数均>20；受控日志源快照只留ignored data | 不取代加性或连续交互结果；病例内未调整描述，不能当人群发病风险；稿件未改 |

## 工作索引（兼运行索引；追加，不覆盖）

不另建 ANALYSIS_RUN_INDEX.md。旧条目是证据回溯，不补造历史run ID、执行日期或命令。操作合同见CLAUDE.md。

| 工作 ID / 日期 | 任务 / 父版本 | 实际分析执行 / 状态 | 角色 | 产物与证据 | 核验 / 采用影响 / 下一步 |
| --- | --- | --- | --- | --- | --- |
| LEGACY-C（执行日期未知；10-16回溯） | PRS311 / 13-gene C | 存在既有汇总结果；原run记录未恢复 | 当前采用；post hoc/provisional | manuscript/datas/manifest.json；RESULT_REGISTRY精确selector | 副本核对，不是本轮重跑；C在主稿位置接替旧panel，原执行/批准链未知 |
| GWPRS-PLAN-v1_1（文内10-16；回溯） | 新PRS扩展草案 | 方案；未有新增关联运行凭证 | 探索性计划；计划被v2修订 | plans/genomewide_prs_analysis_v1.md | superseded_by=GWPRS-PLAN-v2（仅方案关系）；保留原文 |
| GWPRS-PLAN-v2（文内10-16；回溯） | 审计v1.1 | 仅方案/静态及公开来源审查，非AoU评分运行 | 探索性计划；未批准锁定 | plans/genomewide_prs_analysis_v2.md Material Passport | supersedes=GWPRS-PLAN-v1_1；权重/匹配/小试/新模型未验证；不替代C结果 |
| DOC-20261016-01 / 2026-10-16 | 记录流程审计；父=10-15文档基线 | 未运行科学分析；文档整理完成 | 技术诊断/记录治理 | 本次修改的CLAUDE、README、PROJECT_STATUS、RESULT_REGISTRY、ANALYSIS_CHANGELOG、NEXT_STEPS及manuscript/README；manifest | 21/21汇总哈希匹配；CC三项及AGE两项主要估计与registry原字符串一致；7份修改文档可读；未审完DOCX accepted view/图/补充；C数值与采用位置不变；1012身份待查 |

| DOC-20261016-02 / 2026-10-16 | 作者确认稿件版本；父=DOC-20261016-01 | 未运行分析；版本登记完成 | 当前采用：1012主稿/补充稿；1008历史/被替代（仅稿件身份） | 作者本轮原话“1012稿件为最新版本”；MANUSCRIPT_VERSIONS记录现存文件SHA-256；PROJECT_STATUS、两处README、NEXT_STEPS、RESULT_REGISTRY同步 | supersedes=1008稿件身份；不证明分析变化，登记数值不变；文件存在/哈希及入口指向已核对，1012内容一致性未检查；关闭身份待确认，保留数值映射待办 |

追加模板：`工作ID/日期 | 做什么/父版本 | 是否运行/完成或失败阶段 | 用途+科学性质 | 命令/环境/源与输入锁/日志/产物路径（无运行填不适用） | 核验方法与范围/决定或基线不变/supersedes关系/未决与下一步`。

| AGE-INT-20261016-01 / 2026-10-16 | 作者授权补做原研究问题中遗漏的连续burden×PRS年龄交互；父=AGE-AUDIT-20261016-01 | 代码与配置已准备；后台运行待启动；AoU常数只读访问预检通过 | 目标=主文；作者陈述原意，非本轮宣称预注册/数据盲；沿用provisional C | ANALYSIS_RUN_INDEX.md；config/age_burden_prs_interaction_20261016.json；scripts/run_age_burden_prs_interaction.py | 锁定7,996病例、既有13基因分数、病例内ddof1；HC3主要/HC1对照；不加incident限制；旧加性模型/二元carrier交互均保留；不按显著性决定采用；新结果未产生 |

2026-10-16作者明确要求建立ANALYSIS_RUN_INDEX，覆盖CLAUDE中“不另建索引”的旧约定。本轮索引与本工作索引交叉引用，不补造旧run。原审查的“仅规划/未运行”是当时状态；新run独立登记，不回写历史报告为已运行。

| AGE-INT-20261016-01 completion / 2026-10-16 | 同ID完成事件；不覆盖准备记录 | 后台真实输入分析完成：AoU只读访问、Table3复核、M1 HC3/HC1、独立矩阵验证及聚合诊断图；合成协方差测试单列 | 新增主文分析结果；原研究问题作者陈述，当前实现先锁定再拟合；原C资格provisional | ANALYSIS_RUN_INDEX.md；docs/age_interaction_runs/AGE-INT-20261016-01/REPORT.md及provenance/stage_status/verification；受控日志与源快照留ignored data run；无分析失败/blocked | interaction beta+0.13757，HC3 CI[−0.14782,+0.42296]，P0.34477；HC1 P0.33983；ΔR²0.012844个百分点。代码/输入未变、独立代数核验通过；全原始评分/变异人工裁定/外部13权重逐值审计/非线性及leave-gene模型未执行。新增分析不supersede Table3或binary carrier；稿件未改。最后公开产物完整性/披露核对见final_artifact_verification.json |

| AGE-AUDIT-20261016-02 / 2026-10-16 | 回答<18岁诊断病例计数；父=AGE-AUDIT-20261016-01 | 完成只读表型标记/历史日志核查及冻结病例ID交集核查；未拟合模型 | 技术诊断；不改变科学定义 | data/04_phenotype_diagnostics_person_level.csv.gz的flag_age_at_diagnosis_invalid_under18；data/04_internal_output_04_phenotype_diagnostics.txt；section9/frozen_model_inputs.csv.gz | 上游存在<20例清理标记，精确数仅保留于既有受控日志、不转写此登记；当前冻结病例的交集计数见本轮会话；基线不变，未改稿件/代码/数据；历史清理的合理性仍待原记录审查 |

| AGE-AUDIT-20261016-01 / 2026-10-16 | 1012病例年龄分析validate及连续burden×PRS年龄交互plan；父=DOC-20261016-02/C | 完成审查；实际运行冻结输入的内存既有OLS复核、HC3/BP/RESET/支持和相关性诊断、正式19条参考文献在线审计；未运行新增连续年龄交互、Cox或frailty | 技术诊断/探索性方案；post hoc/provisional；建议未视为采纳 | docs/age_burden_prs_review_20261016.md；docs/age_burden_prs_review_20261016_refs/；执行命令见本轮会话；numpy2.4.3/pandas2.3.3/statsmodels0.14.6；两冻结输入哈希与section9锁相符，age IDs另记哈希 | 7,996及AGE-JOINT两斜率/CI/P/R²与汇总在浮点精度内相符，非完整原始数据流水线复现；发现异方差、稀疏burden支持、<18清理待审、Cox预期值版本冲突及S1引用错号；书目脚本18OK/1UNVERIFIED且原文身份/卷页仍有核查限制；基线不变；下一步先解决报告P0/P1，再锁定独立交互run；稿件/代码/原始数据未改 |

失败/中止照此登记；不能因没有可发表结果而省略。仅状态推进可追加同ID事件；保留旧事件。补证据需注明补录日期与原事件日期的区别。


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


## AGE-QUART-20261016-01 / 2026-10-16

Real frozen-input Q1/Q4 descriptives completed (not synthetic-only), shared quartiles in original7,996 cases; selected3,998; early<45 fixed; existing13-gene carrier. Source/command/software/input and output hashes/status: docs/age_interaction_runs/AGE-QUART-20261016-01/provenance.json; [report](docs/age_interaction_runs/AGE-QUART-20261016-01/REPORT.md). No participant records/scores written; controlled logs and snapshot ignored. No regression or significance test; no failure/blocked stage. Supersedes=None; prior median split, additive and continuous-interaction results retained; manuscript unchanged.


## PUSH-PREFLIGHT-20261017-01 — 2026-10-17

- Task: user-requested publication of root Markdown files, docs/, manuscript/ to martina-yu/Common_Rare_Variation. Role: repository maintenance; no analysis run; execution stopped at preflight, publication incomplete. Scientific baseline unchanged; no scientific adoption or result changes.
- Size inventory excludes Git internals and ignored files: 128 files, 24850745 bytes. Three manuscript main DOCX files individually exceed 5,000,000 bytes. CLAUDE.md requires total push content below 5 MB and stopping when exceeded; no commit or push performed.
- Evidence: /tmp/common_rare_push_preflight_20261017.json (local preflight inventory). Destination ls-remote failed authentication; destination-relative commit/file size and external-sharing review remain unchecked. Inventory is a working-tree snapshot, not an authenticated destination diff.
- Next: resolve allowed publication scope/size and GitHub authentication, then check destination-relative outgoing history, size, and sharing eligibility before publication. No manuscript/result consistency audit performed.


## ROOT-MD-COMMIT-20261017-01 — 2026-10-17

- Task: user-authorized commit of Markdown files directly in the repository root only. Role: repository maintenance; no analysis run; baseline unchanged; no new scientific results or adoption changes.
- Verification: nine root Markdown files totaled 145,210 bytes before this entry; staging was empty; configured commit email matches CLAUDE.md. Commit scope excludes all subdirectories, data/, datas/, docs/, and manuscript/.
- Next: publication remains pending; destination authentication and cumulative outgoing-history size require verification before any push. No scientific consistency audit performed.

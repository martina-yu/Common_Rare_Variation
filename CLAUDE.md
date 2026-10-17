# 注意事项

## 每次工作结束、最终报告之前的必做步骤（2026-10-16）

适用于方案审计、新分析、重跑、中止和失败尝试；入口先读 PROJECT_STATUS、RESULT_REGISTRY、ANALYSIS_CHANGELOG、NEXT_STEPS。详细工作索引直接放在 ANALYSIS_CHANGELOG，不另建 ANALYSIS_RUN_INDEX；本 checkout 无 AGENTS.md，规则在本文件统一维护，后续 agent 不得仅依赖工具自动加载，须主动读取。

1. 在 ANALYSIS_CHANGELOG 的工作索引追加唯一工作 ID，注明日期、任务/父版本、是否实际运行分析、执行状态、角色、产物/证据路径、核验范围及下一步。只审计划也记一行；失败记失败阶段、已有产物是否可用及日志，不登记虚构结果。
2. 分开记录三轴：执行状态（未运行/运行中/完成/失败/中止）、用途（当前采用/敏感性/探索性/历史/被替代/技术诊断）、科学性质（post hoc或预设、provisional或frozen）。当前采用不等于确认性；方案被修订不等于结果被替代。
3. 有新数值才追加 RESULT_REGISTRY 版本行；写明 result ID、run ID、cohort/model/term selector、单位、估计/CI/P、来源和哈希、稿件位置及核验方法。没有结果只登记工作和待办，不填零、不借用旧结果。有科学决定才追加 changelog 的旧→新、理由、证据和决定来源；审计建议不得冒充采纳。
4. 重跑使用独立 run 输出位置，开始前保存源代码/配置/输入版本或哈希、命令/环境/seed及父版本；受控数据只记安全的位置标识和版本，不抄参与者记录。不得覆盖旧输出；新旧版本以 supersedes/superseded_by 相连。缺少旧快照则明确不可恢复，不能伪造。不得按显著性或mtime晋升结果。
5. 有采用变化才更新 PROJECT_STATUS 的一分钟摘要及 NEXT_STEPS；无变化明确写“基线不变”。对照 registry→输出 selector/哈希→主稿/表图/补充材料，记录通过、差异、未检查三类及检查范围。只验副本哈希不等于独立重算；未解决冲突保留待办，不宣称全稿一致。
6. 最终报告列出工作 ID、是否运行、当前基线是否改变、核验和未决问题。完成上述记录再发送报告；无变化无需把四份文档机械重写。


每次操作之前都需要阅读本文件的要求。

## AoU Controlled Tier Data Egress Security Standards

**Allowed (read-only)**: `gcsfuse -o ro` mounts, `gsutil cat -r` (byte-range reads), `tabix` region queries, BigQuery reads.

**Prohibited**:
- `gsutil cp/mv` or `gcloud storage cp` of participant-level data outside the workspace bucket
- Manual "Download" clicks in the JupyterLab file browser on participant data (leaves no audit trail)
- Writing large un-aggregated intermediates (raw genotype matrices, carrier extractions) to disk — aggregate in memory first

**What counts as downloadable**: Only fully aggregated files with no participant-level rows (e.g., group-level summary stats). Anything with one row per person_id — even a small file like `data/prs_scores.csv` — is never downloadable, regardless of size. The 5MB limit is just an extra project-level rule; AoU's actual criterion is row-level vs. aggregated, not file size.

**File/mount hygiene**:
- 5MB limit is on the **sum of all file sizes in a git push**, not any single file. A push with ten 600KB files (6MB total) violates this just as much as one 6MB file. **Before every `git push`, compute the total size of all files being pushed (new + changed, across all commits not yet on origin) and confirm the sum is under 5MB before proceeding.** If it's over, stop and flag it — don't push and explain after.
- Aggregate/filter before writing large variant or genotype outputs to disk
- BigQuery writes only to your own project/dataset
- `fusermount -u` before remounting; mounts stay read-only throughout

**Reference data**: Public reference downloads (PGS Catalog, Ensembl) are fine — inbound public data is not participant data.

**`plans/papers/` 永不 push**（用户指令，2026-08-19）：这个目录是用户自己整理的文献库（第三方论文 PDF/补充材料，当前约 32MB），**只作为本地参考资料，任何情况下都不要提交或推送**。已加入 `.gitignore`（第 27 行 `plans/papers/`）——注意原有的 `plans/*.pdf` 规则**匹配不到子目录**里的 PDF，所以必须单独写这一条，否则一个 `git add -A` 就会把 32MB 全部加进 staging、直接爆 5MB 上限。以后往这个目录里加论文不需要做任何 git 操作，也不要为了"方便复现"把它们挪到别的被跟踪的目录下。

**If an alert fires**: check history for `cp/mv` outside AoU, confirm mounts were read-only, check for large files on disk, recall any manual UI downloads, and contact AoU support for audit logs if needed.

**Publication/sharing threshold**: every cell in an externally shared aggregate table needs n≥20 (round up 5–19; suppress or merge below 5). Broad totals (e.g., overall case/control n) are already well above this and need no adjustment.

**n≥20 规则同样适用于日志和审计表**（2026-08-19 发现）：`scripts/script_output_py/` 里的运行日志是被 git 跟踪的，会随 push 一起外传，所以脚本 `print` 出来的分组计数也必须按 n≥20 披露规则处理，不能只管 `results/table/` 里的 CSV。特别注意两类容易漏的地方：①"被抑制格子"的审计表本身如果记录了精确原始 n，等于把抑制掉的数字又写了回去——只记 `<5` / `5-19` 档位；②诊断类输出（如 Cook's distance top-N 的分组计数）也算分组计数。

**De-identified scatter plots — PI ruling, 2026-08-18**: scatter plots that render one point per participant (e.g. `17_burden_rare_vs_age_at_diagnosis_scatter.png`, n=8,328 points, each = one case's genetic score vs. age at diagnosis) are an **accepted exception** to the row-level rule in this project and may be pushed to the public repo. Rationale: the points carry no `person_id` and encode only two continuous derived variables. This was raised as a compliance question during the pre-push review of 2026-08-18 and ruled on explicitly by the PI. Note this is an *interpretation* — the underlying AoU criterion is row-level vs. aggregated, and a stricter reading would require binning/density plots instead. Do not re-litigate this per push, but do flag it if the context changes (e.g. adding a third axis, or plotting a directly identifying variable).

**Egress limit increases**: never bypass Workbench's built-in download review to dodge a block. If a genuinely aggregated file gets flagged, request a temporary threshold increase via the [Egress Alert Policy](https://support.researchallofus.org/hc/en-us/articles/4407354684052-Egress-Alert-Policy) form — don't try to "trigger and then explain," as repeated triggers can auto-suspend the account.

---

**Repeats found across the two source sections**: the total 5MB file-size limit (sum of all files < 5MB), the read-only mount requirement, and the "don't manually download via UI" rule each appeared twice in your original doc — consolidated above into single mentions.


## 结果输出

- 'workspace/com_rare_BC/plans/'：里面是每一步骤的计划书，请参考步骤但不能完全相信里面的方法，在实施之前需要查阅文献确定正确性和可行性。
- 'workspace/com_rare_BC/scripts/output/'：里面是python脚本的输出，保存在output里面方便debug和查看。使用shell命令执行即可，例如：'python script_01.py > output/output_script_01.txt' 2>&1
- 'workspace/com_rare_BC/scripts/'：存储python或者R脚本，优先度：python=R>jupyter notebook
- 'workspace/com_rare_BC/data/'：存储过程中产生的矩阵等
- 'workspace/com_rare_BC/results/':分table/和figure/和report/三个文件，分别存储输出的table文件、图片文件、report文件。

图片输出优先用作图函数例如matplotlib和seaborn自身的默认参数，不要用claude做图的风格，不要有grid。以防被认为画图像AI画的。

**图片格式**（2026-08-19 定稿）：**先测体积，用更省的那个格式。** 具体顺序：
1. **优先 SVG**（矢量、可无限放大、通常最小）。
2. **点数很多的图**（成千上万个点的散点、Cook's 距离图等）SVG 会反而更大，这时候落到栅格格式。
3. **落到栅格时，PNG 和 JPG 都存一次，比大小，留小的那个，删掉大的。** 不要凭直觉选——本项目 22 张图实测 **PNG 全部更小**（合计 PNG 1.78MB vs JPG q95 2.67MB，+50%；JPG q85 也要 1.36MB 且已有可见压缩痕迹）。原因是这些图是纯白底线条图（箱线/柱状/ROC/热图），PNG 的无损压缩对大片同色区域效率极高，而 JPEG 要为线条边缘花码率。**所以本项目的默认栅格格式是 PNG**，JPG 只在实测确实更小时才用（例如照片式的连续色调图，本项目目前没有）。
4. 如果用 JPG，必须显式写 `pil_kwargs={"quality": 85, "subsampling": 2}`——matplotlib 存 jpg 默认 quality 只有 75，不指定会明显更糊。
5. 5MB 累计上限要盯紧：换格式省下来的每一百 KB 都是下一批图的额度。

## pipeline工具

使用snakemake工具辅助生物信息分析。

## 环境

记得存储环境文件到yaml，方便后续clone

## github

需要把文件及时push到github上，并且严格遵守AoU数据库的管理规定上传。

认证：GitHub PAT 已通过 `git config credential.helper store` 缓存在 `~/.git-credentials`，不要把 token 明文写进本文件或任何会被提交的文件里。push 前如需新 token，直接在对话里提供即可，不落盘到仓库内的文件。缓存格式必须是 `https://用户名:token@github.com`（带用户名），只写 `https://token@github.com` 会被 `git-credential-store` 静默忽略，导致 push 失败。

commit身份：本repo的commit author邮箱必须用 `2889642595@qq.com`（GitHub账号`martina-yu`已验证的邮箱），已通过 `git config --local user.email` 设置。**不要用登录Claude用的邮箱（如`zimengmeng0703@gmail.com`）提交**——邮箱对不上GitHub账号的已验证邮箱列表，commit不会被计入GitHub主页的contribution图。
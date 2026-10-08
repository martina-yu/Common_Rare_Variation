# Current manuscript: thirteen-gene revision, 2026-10-08

The latest full revision is **medrxiv_v1_main_13gene_revised_1008.docx**, with companion **medrxiv_v1_supplement_13gene_revised_1008.docx**. It replaces the original main-paper burden with the thirteen-gene PTV plus eligible pathogenic/likely pathogenic missense score.

- Readable main and supplement: `revision_13gene_1008/main_13gene_revised.md` and `revision_13gene_1008/supplement_13gene_revised.md`.
- Matching manuscript figures, aggregate tables and validation: `revision_13gene_1008/`.
- Revision details and remaining author/variant-review limitations: `revision_13gene_1008/revision_log.md`.
- New main flow, ROC, AUC and joint-distribution vector figures: `figures/figure*_current*.svg` and `figures/figure1_flow.svg`.
- Gene-specific interaction box plots and forest plots: `figures/supp_gene_PRS*.svg`.

The current main results use the thirteen-gene score consistently. Historical score comparisons and full sensitivity plots are in the supplement. Original source documents and local archives are preserved. Participant-level datasets are excluded. The text below documents older versions and does not designate the current manuscript.

## Historical version notes

> 2026-10-07 update: The author-designated current main revision is `medrxiv_v1_main_revised_1006.docx`. A checksum-verified copy, together with the available `medrxiv_v1_supplement_1002.docx`, is preserved in `archive/2026-10-07_medrxiv_1006_snapshot/`. The supplement remains the older analysis. Historical descriptions below have not been rewritten.

> 2026-10-08 interaction additions: `medrxiv_v1_main_interaction_revision_1008.docx` and `medrxiv_v1_supplement_interaction_revision_1008.docx` are scoped review copies adding the current 13-gene interaction audit, Tables S7–S8, gene-group PRS boxplots and the interaction forest plot (Figure S7A–B), and Discussion citing existing reference 5 (Barnes et al., 2020). These do not replace the author-designated main automatically. Other legacy source text and tables still require a full refresh. Standalone additions: `interaction_results_discussion_additions.md`; validation: `../results/report/60_interaction_manuscript_assets_summary.md`.

# 手稿版本索引

整理日期：2026-10-06。最新版以作者指定的 **1002 后缀**为准，不以文件修改时间判断。

| 文件／目录 | 角色 | 使用说明 |
|---|---|---|
| `medrxiv_v1_main_revised_1002.docx` | 最新主文修订底稿 | 已插入当前Figure 2–6共11张SVG，更新图注及图号引用；正文结果数字尚未全面更新为new missense version |
| `archive/2026-10-06_before_figure_replacement/` | 替换图片前的原稿备份 | 保留未经本轮图片替换的1002主文 |
| `medrxiv_v1_supplement_1002.docx` | 最新配套补充底稿 | 原文件保留，仍对应旧 PTV 分析 |
| `archive/2026-10-06_pre_missense/` | 历史稿及配套材料 | 47 个文件，保留原名称和相对目录；不作为当前结论来源 |
| `analysis_review_20261006.md` | 最近分析与修稿衔接记录 | 版本定义、结果来源、待修改位置 |

归档包括早期综合／正式报告（含 PDF、LaTeX 及编译文件）、medRxiv Markdown 草稿、旧主文／补充 LaTeX 和 PDF、图表、Excel、构建脚本以及 notebook checkpoints。其中综合总结标注日期为 2026-08-18，Markdown 草稿标注为 2026-09-26；未给其他文件推断创作日期。

归档清单见 `archive/2026-10-06_pre_missense/manifest.csv`，记录原路径、现路径、大小及 SHA256。全部 49 个原文件（47 个归档、2 个保留）移动后逐文件校验通过，未删除原内容，也未改写 DOCX 正文。

**复现注意：** 旧 `build_assets.py` 已随历史稿归档，但其中输出路径仍硬编码指向工作区的 `manuscript/`。不要直接运行该历史脚本以生成当前稿件资产；先调整输出路径。历史文档中的旧路径引用作为历史记录保留，未批量改写，也未验证重新编译。

后续以 1002 DOCX 为修订底稿，使用 **new missense version** 指称本轮 13 基因 PTV＋合格 P/LP missense 的 C 方案；9 基因 missense 的 B 方案和报告 28 的早期探索分析必须另行标明。旧 9 gene PTV 建议保留为补充比较基准，正文是否完全移除其结果待修稿时决定。
# Common_Rare_Variation

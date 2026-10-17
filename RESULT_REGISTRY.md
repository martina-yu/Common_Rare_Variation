# RESULT_REGISTRY — 当前数值与来源

## AGE-FOUR-20261016-01：四组病例年龄描述（不替代模型结果）

2026-10-16真实7,996病例完成；定义在汇总前锁定：13基因carrier，PRS按全部病例共享median(<low/>=high)，早发<45。未调整，无P值，无新回归；supersedes=None。拟用于描述性表图，当前稿件未写入。命令/版本/输入源哈希/状态见docs/age_interaction_runs/AGE-FOUR-20261016-01/provenance.json。

| Result ID / group selector | n | age median [Q1,Q3], years | early n / % [Wilson95% CI] |
|---|---|---|---|
| AGE-FOUR/non-carrier-PRS-low | 3675 | 62.08350444900753 [52.766598220396986,69.55099247091033] | 349 / 9.496598639455783% [8.590601887891422,10.487182971195153] |
| AGE-FOUR/non-carrier-PRS-high | 3778 | 60.90212183436003 [52.15879534565366,68.49486652977413] | 398 / 10.534674430915828% [9.595500656720825,11.554023127655684] |
| AGE-FOUR/carrier-PRS-low | 323 | 53.98494182067078 [43.460643394934976,64.14784394250513] | 90 / 27.86377708978328% [23.256544888895966,32.99135569386024] |
| AGE-FOUR/carrier-PRS-high | 220 | 55.393566050650236 [45.12046543463382,65.46885694729637] | 55 / 25.0% [19.740284404356355,31.117791405400432] |

Source=docs/age_interaction_runs/AGE-FOUR-20261016-01/four_groups.csv；SHA256=519803771f91988761697790b52e5da1c11823aafe0da1e8017396026c92475a。Selector=group字符串与表一致。年龄P10/P90/mean/SD及完整年龄分箱见同目录CSV；age_distribution.csv SHA256=6da4c2fe96c61203f1f05788b195a182cb98e95d7adae89322f4ac9e66a48a5d。

指定描述比较：carrier低−carrier高，中位数差−1.4086242299794591年、早发比例差+2.86377708978328个百分点；carrier低−noncarrier高，中位数差−6.917180013689254年、早发比例差+17.329102658867452个百分点。仅描述差值，无差值CI/P，不能当因果或交互结论。来源comparisons.json，SHA256=b4e0d15da5d23e172827c75d383f2affbd677bbd58d1347347913e884c0088eb。核验：父冻结输入hash前后相符、carrier与正burden一致、分组互斥且和为7996、公开计数全部>20；个人记录/个人分数未写出。

## AGE-INT-20261016-01：新增连续burden×PRS年龄交互（不取代加性Table3）

2026-10-16真实后台运行完成，非仅合成测试。作者决定使用现有7,996病例、13-gene C B_raw_C/PRS311、病例内sample SD(ddof1)、PC1–10；HC3为新增模型主要推断，HC1为对照。作者陈述此交互属于原研究问题中遗漏的步骤；本轮实现先锁定再拟合，不宣称预注册/数据盲。目标主文，稿件尚未写入。原C资格provisional，未经本次重注释或人工裁定。**Supersedes=None：原AGE-JOINT、binary carrier年龄交互、CC及CV数值保留。**

Run/工作ID=AGE-INT-20261016-01；cohort=冻结diagnosis-age病例7996；M1=`age_at_diagnosis ~ 1 + PRS311_z + burden_z + PRS311_z:burden_z + PC1..PC10`。M1主效应以另一分数Z=0为条件，与Table3加性效应不可混写。所有CI均为95%，P双侧、normal-reference robust Wald；无乘积项再标准化。新增交互不能解释为发病hazard或病因学时间交互。

| Result ID / exact term selector | Estimate / unit | SE | 95% CI | P | Source selector |
|---|---|---|---|---|---|
| AGE-INT/M1/HC3/PRS311_z | −0.4314947895709711 年/SD_P | 0.12865431759876142 | −0.683652618520121 to −0.17933696062182108 | 0.0007967875036327197 | model_coefficients.csv; model=M1_interaction,covariance=HC3,term=PRS311_z |
| AGE-INT/M1/HC3/burden_z | −1.7778976526254473 年/SD_B | 0.15353258227637487 | −2.078815984340575 to −1.47697932091032 | 5.208422897166964e-31 | same; covariance=HC3,term=burden_z |
| AGE-INT/M1/HC3/PRS311_z_x_burden_z | +0.13756989831270217 年/(SD_P×SD_B) | 0.14561084461409704 | −0.14782211288938613 to +0.42296190951479046 | 0.3447723108105911 | same; covariance=HC3,term=PRS311_z_x_burden_z |
| AGE-INT/M1/HC1/PRS311_z_x_burden_z | +0.13756989831270217 年/(SD_P×SD_B) | 0.1441276445421457 | −0.14491509416649428 to +0.4200548907918986 | 0.33983039307422513 | same; covariance=HC1,term=PRS311_z_x_burden_z |

| Result ID | Value / unit | Source selector |
|---|---|---|
| AGE-INT/M1/R2 | 0.026405514033252753 | model_R2.csv; model=M1_interaction,R_squared |
| AGE-INT/M1/adjusted_R2 | 0.024819855261319912 | same; adjusted_R_squared |
| AGE-INT/M0/R2-reproduced | 0.02627707200834939 | same; model=M0_additive,R_squared; retained Table3 baseline |
| AGE-INT/M1-minus-M0/delta_R2 | 0.00012844202490336354 = 0.012844202490336354个百分点 | interaction_comparison.json; delta_R_squared/delta_percentage_points |
| AGE-INT/M1/partial_R2 | 0.00013190818580012412 | same; partial_R_squared=delta_R2/(1−R2_M0) |
| AGE-INT/M1/HC3/Wald_chi2 | 0.892605154726034; df1 | same; Wald_chi2; P同主要交互 |

来源目录：[docs/age_interaction_runs/AGE-INT-20261016-01/](docs/age_interaction_runs/AGE-INT-20261016-01/REPORT.md)。`model_coefficients.csv` SHA256=dc61301a248408faa5178ce495d32607f341c43482abe5b5d0b9ce6b482ef49f；`model_R2.csv`=08d3c552c96b7179c7ce73135ce0edde51cb7397bac7c4d1332d51cdbc1521d7；`interaction_comparison.json`=c98a4caab9087306b7ee0a7ee5ce8d6a8cf49a8490347cf14fad8656824cf1e6；`provenance.json`=bc9cae49694fbc03ff300318e02ad0e143271ee0d7a2a1581d2cd51c292788b2。命令、源/配置/输入别名hash、软件版本见provenance/software/stage_status；真实Table3数值复核与独立least-squares/sandwich公式均通过，输入源hash前后不变；不能扩大为从原始WGS完整独立复现。受控输入/个人分数不在此登记。

M1技术诊断（同run，详见报告及diagnostics_summary.json）：rank14/14；标准化condition number33.0696；max leverage0.0387344；max Cook0.00826854；max abs interaction DFBETAS0.228877；Koenker BP LM35.9434893/P0.000605448；HC3 RESET Wald0.963235/P0.617783；residual SD11.6272363年。诊断属于技术核查，不自动作为稿件结论；未显著RESET不证明线性，未做影响点删除/非线性替代模型/leave-gene分析。条件斜率和全部PC参数保留在相应完整aggregate CSV，未按显著性筛选。

解释：连续交互证据不足，CI跨零不等于证明无交互；新分析按作者要求用于主文，不因P值改变位置。稿件/图表写入及既有Cox/frailty冲突另待处理。

2026-10-16 / DOC-20261016-02：作者确认当前稿件为main_1012 + supplement_1012，取代1008的稿件身份。下方数值及稿件位置来自此前汇总/1008核对，尚未逐项验证1012；本轮不改数值、不提升核验等级，也不推断新模型已运行。

## 核验口径与尚无结果的扩展

### AGE-AUDIT-20261016-01：2026-10-16冻结病例年龄输入复核（不取代C）

工作/运行ID均为AGE-AUDIT-20261016-01。用途=技术诊断，post hoc/provisional；不是新增连续burden×PRS年龄交互结果。实际在内存以相同公式复核OLS；n=7,996、两主要beta/HC1 CI/P和保存R²在浮点精度内相符（beta差约10⁻¹⁴，R²差约2×10⁻¹⁶），未认证从原始表型和基因型重建的全链条复现。现有AGE-JOINT采用基线不变；稿件定位为1012 Table3/Methods age analysis。新诊断不自动成为论文结果。

来源为[本轮完整报告](docs/age_burden_prs_review_20261016.md)一/2、一/3，SHA256=`998ebd1c85b3f14fb1ab865c02c6c5d8e5bbadab215d7cf04ab24a6fce4ee1bf`；内存执行命令见该工作会话，环境numpy2.4.3/pandas2.3.3/statsmodels0.14.6。报告列两冻结输入与原section9锁相符的SHA256，以及另行记录的age ID哈希。下表为本轮新增核查数值；未提供CI/P处为未计算或不适用，不以零填补。

| Result ID / selector | Cohort/model | 核查值 / 单位 | CI/P及限制 |
| --- | --- | --- | --- |
| AGE-AUDIT/PC_R2 | 固定7,996病例；PC1–10+intercept | R²=0.0019538053109502407 | 支持Table3的0.20%；无新增CI/P |
| AGE-AUDIT/burden_PC_R2 | 同病例；burden_z+PC | R²=0.02493862978745931 | 无新增CI/P |
| AGE-AUDIT/delta_PRS_R2 | AGE-JOINT减burden+PC | 0.0013384422208898572；0.133844个百分点 | 样本内增量；非交互ΔR² |
| AGE-AUDIT/HC3_SE | AGE-JOINT；PRS/burden | 0.12856197737065184 / 0.1467698568815113 年/SD | SE敏感性；未更换采用推断 |
| AGE-AUDIT/BP | AGE-JOINT残差，完整X，studentized/Koenker BP | LM=36.11381826794074 | P=0.00031062242927482313；诊断探索性 |
| AGE-AUDIT/RESET | AGE-JOINT；fitted powers2/3、HC3 | Wald=0.6787497939301135 | P=0.7122153916440128；未显著不证明线性 |
| AGE-AUDIT/support | 固定年龄病例；B_raw_C | carriers543；zero_fraction0.9320910455227613；19个不同值 | 联合尾部支持仍待审 |
| AGE-AUDIT/scaling | 同病例；raw PRS / raw burden | mean/SD=0.35447389832717613/0.6079261888021756；0.10002667950856513/0.40125455656464104 | ddof=1；不能混用全队列SD |
| AGE-AUDIT/correlation | PRS311与B_raw_C；full/controls/allcases/agecases | Pearson=−0.00350300699578601/−0.004934389071332944/−0.07087945100703164/−0.06885361372378013 | 无新增CI/P；不证明collider为唯一解释 |
| AGE-AUDIT/age_distribution | 固定7,996病例 | mean60.20049073133418、SD11.783856553501243年 | 主文全部病例59.9年不能直接与本子集混比 |
| AGE-AUDIT/influence | AGE-JOINT | max leverage0.03872950293716335；max Cook0.0060767334851362645；max abs internal studentized residual3.527886342829212 | 非新增交互诊断；无CI/P |

正式参考文献脚本JSON归档于`docs/age_burden_prs_review_20261016_refs/qc/reference_audit.json`，SHA256=`0cebbd1e3d4dcb6f2bfcb5a0a423953c27d4e05278b9b81d822baaaa6e5e3671`。18OK/1UNVERIFIED是脚本存在性/有限作者审计，不能代替报告列出的原文身份、卷页与主张核查限制。

以下段落保留此前汇总核对的证据口径；AGE-AUDIT-20261016-01新增的冻结输入内存复核范围以上方记录为准。执行记录见 ANALYSIS_CHANGELOG 的工作索引。M/D/L只描述证据来源，**都不自动表示独立核验**。此前所列C主要数值已由汇总副本逐行核对；该轮21份manifest哈希相符，当时未独立重算模型。D条目仍仅有稿件/配置支持。独立核验须另记日期、核验人/方法、独立输入或实现、比较容差和证据；没有记录即“未证明独立核验”。

GWPRS-PLAN-v2：仅方案（plans/genomewide_prs_analysis_v2.md），探索性扩展；PGS005346/PGS005347尚无本轮AoU评分、系数、交互或AUC产物，不建立数值结果行，不supersede C。v1.1保留作为历史计划。


核对日期：2026-10-15。下列数值逐行取自现存 aggregate CSV；保留原始数值字符串。Primary=当前主稿主分析位置；C 仍是 post-hoc extension、provisional definition，并非前瞻 confirmatory。模型与 cohort 定义见 [PROJECT_STATUS](PROJECT_STATUS.md)。历史结果不因日期较新而自动替换主结果。

**证据等级**：M=与既存 manifest SHA-256 相符的汇总副本；D=稿件/配置一致但原输出未恢复；L=本地汇总副本，另列 SHA-256。没有真实数据重跑。`results/` 在本 checkout 缺失；下方历史路径是 provenance，不表示文件现存。

## Cohort / score 定义

| Result ID | Current value | Status / primary? | Generated by | Source / model |
| --- | --- | --- | --- | --- |
| COH-QC | 8,528 cases / 134,727 controls; N=143,255 | Primary descriptive; D | [02](scripts/02_sample_qc.py) | Current main Table 1 / Figure 1; historical results/table/05_supp_table2_sample_qc.csv |
| COH-CC | 8,528 cases / 132,676 controls; N=141,204 | Primary analytic; D | [04](scripts/04_gene_selection.py), [08](scripts/08_case_control_models.py) | config/rowlands_13gene_extension.yaml; CC-JOINT / CV-4 |
| COH-AGE | N=7,996 cases | Primary analytic; D | [04](scripts/04_gene_selection.py), [09](scripts/09_age_analysis.py) | same config; AGE-JOINT |
| PRS-LOCI | 304 matched; 302 after CHEK2 tagging exclusions | Primary definition; D | [03](scripts/03_common_variants_prs.py) | Current main Methods; original match audit not restored |

## Table 2 / Table 3 — coefficients

| ID / term | Current value | Unit | 95% CI | P | Status | Generated by | Exact source selector |
| --- | --- | --- | --- | --- | --- | --- | --- |
| CC-JOINT / PRS311_z | 1.5539613762158249 | OR per SD / SD-product | 1.5181995383411027 to 1.5905655994396928 | 2.0615358495618704e-301 | Primary / M | 08 → 09 | [case_control/joint_model_coefficients.csv](manuscript/datas/case_control/joint_model_coefficients.csv); version=C; term=PRS311_z |
| CC-JOINT / burden_z | 1.268539929790285 | OR per SD / SD-product | 1.2505189230609375 to 1.2868206340560315 | 6.929211000510023e-233 | Primary / M | 08 → 09 | [case_control/joint_model_coefficients.csv](manuscript/datas/case_control/joint_model_coefficients.csv); version=C; term=burden_z |
| CC-JOINT / PRS311_z_x_burden_z | 0.9715467971748014 | OR per SD / SD-product | 0.9579441513245703 to 0.985342598308533 | 6.008042222779803e-05 | Primary / M | 08 → 09 | [case_control/joint_model_coefficients.csv](manuscript/datas/case_control/joint_model_coefficients.csv); version=C; term=PRS311_z_x_burden_z |
| AGE-JOINT / PRS311_z | -0.4328693828111374 | years/SD | -0.6845674981904504 to -0.1811712674318245 | 0.0007496712032446549 | Primary / M | 09 | [age/age_model_coefficients.csv](manuscript/datas/age/age_model_coefficients.csv); version=C; term=PRS311_z |
| AGE-JOINT / burden_z | -1.8197988475987579 | years/SD | -2.1064688042829327 to -1.5331288909145828 | 1.546341227737796e-35 | Primary / M | 09 | [age/age_model_coefficients.csv](manuscript/datas/age/age_model_coefficients.csv); version=C; term=burden_z |

08=[case models](scripts/08_case_control_models.py), 09=[age / combined output](scripts/09_age_analysis.py). Historical combined source: `results/table/45_rowlands_fixed_cohort_coefficients.csv`.

## Table 2 / Figure 3 — discrimination

All rows below: version=C, Primary/M, CV-4; generated by [08](scripts/08_case_control_models.py) plus [13 --task discrimination](scripts/13_figures.py). Repeat-1 CI is **not** a confidence interval for the repeated mean.

| Model / selector | Mean AUC | Repeat range | Repeat 1 AUC | Repeat 1 95% CI | Source |
| --- | --- | --- | --- | --- | --- |
| Covariates only | 0.7451749929096977 | 0.7450630747719467–0.7453136340205995 | 0.7451933541270283 | 0.7406281059183148–0.7497586023357417 | [case_control/four_models_AUC_summary.csv](manuscript/datas/case_control/four_models_AUC_summary.csv) |
| PRS + covariates | 0.7682268278025771 | 0.7680966169430112–0.7683627728415913 | 0.7682652829528375 | 0.7637318418856908–0.772798724019984 | [case_control/four_models_AUC_summary.csv](manuscript/datas/case_control/four_models_AUC_summary.csv) |
| Burden + covariates | 0.758801944153391 | 0.7587056934589967–0.7589317931798701 | 0.758831393778363 | 0.7543489182512587–0.763313869305467 | [case_control/four_models_AUC_summary.csv](manuscript/datas/case_control/four_models_AUC_summary.csv) |
| Joint | 0.7815062641738876 | 0.7813705397346251–0.7816407850364587 | 0.7815450203508927 | 0.7771452789150765–0.7859447617867089 | [case_control/four_models_AUC_summary.csv](manuscript/datas/case_control/four_models_AUC_summary.csv) |

| Comparison | Mean delta AUC | Repeat 1 delta | Repeat 1 95% CI | Repeat 1 paired DeLong P | Status / source |
| --- | --- | --- | --- | --- | --- |
| joint_minus_PRS311_baseline | 0.013279436371310582 | 0.013279737398055325 | 0.011326280397889348–0.015233194398221302 | 1.6792922966676785e-40 | Primary / M; [case_control/CV_delta_AUC_summary.csv](manuscript/datas/case_control/CV_delta_AUC_summary.csv); version=C |
| C_minus_B | 0.0002912745741760037 | 0.000289806737365339 | 7.866140457216475e-05–0.0005009520701585133 | 0.007142178328286228 | Bridging comparison / M; [case_control/CV_delta_AUC_summary.csv](manuscript/datas/case_control/CV_delta_AUC_summary.csv); version=C |

Historical CV summaries: `results/table/45_rowlands_CV_AUC.csv`, `45_rowlands_CV_delta_AUC.csv`; manifest four-model source `scripts/14_four_model_discrimination.py` is historical, current implementation lives in 13.

## Table 3 — explained variance

| Result ID | Current value | Status | Model / generator | Source selector |
| --- | --- | --- | --- | --- |
| R_squared | 0.02627707200834939 | Primary / M | AGE-JOINT / 09 | [age/age_R_squared.csv](manuscript/datas/age/age_R_squared.csv); version=C; R_squared |
| adjusted_R_squared | 0.02481337726503241 | Primary / M | AGE-JOINT / 09 | [age/age_R_squared.csv](manuscript/datas/age/age_R_squared.csv); version=C; adjusted_R_squared |
| PRS_PC_R_squared | 0.0026168446318921568 | Primary / M | AGE-JOINT / 09 | [age/age_R_squared.csv](manuscript/datas/age/age_R_squared.csv); version=C; PRS_PC_R_squared |
| unique_burden_R_squared | 0.023660227376457232 | Primary / M | AGE-JOINT / 09 | [age/age_R_squared.csv](manuscript/datas/age/age_R_squared.csv); version=C; unique_burden_R_squared |

`unique_burden_R_squared` is the increment from adding burden to PRS+PCs, not burden-only R². Historical source: `results/table/45_rowlands_age_explained_variance.csv`.

## Sensitivity and historical coefficients — retain explicit version identity

All rows use the same CC-JOINT or AGE-JOINT specification as appropriate; 10 prepares scores, 08/09 fit and export. “Sensitivity” means the C alternatives; all other non-C versions are historical/bridging, never silently substituted for C. Original source `results/table/45_rowlands_fixed_cohort_coefficients.csv`.

| Version | Cohort / term | Estimate | 95% CI | P | Status / primary? | Source |
| --- | --- | --- | --- | --- | --- | --- |
| A | case_control / PRS311_z | 1.551814859678967 | 1.5161952071807514–1.5882713171202958 | 4.153969248012342e-301 | Historical / bridging; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| A | case_control / burden_z | 1.2439936209764855 | 1.2263865522029862–1.2618534720968215 | 5.563296639408869e-198 | Historical / bridging; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| A | case_control / PRS311_z_x_burden_z | 0.969331847295763 | 0.9558958034070519–0.9829567478304976 | 1.2212129881093842e-05 | Historical / bridging; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| A | diagnosis_age / PRS311_z | -0.4238568993786127 | -0.6756278535404534–-0.17208594521677206 | 0.0009682160551855687 | Historical / bridging; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| A | diagnosis_age / burden_z | -1.7132991621919982 | -2.0014039055910966–-1.4251944187928998 | 2.1511562656486625e-31 | Historical / bridging; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| A_prime | case_control / PRS311_z | 1.5521653317477653 | 1.516527771238624–1.5886403551396373 | 3.010911785585609e-301 | Historical / bridging; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| A_prime | case_control / burden_z | 1.245689646637979 | 1.2280141089209096–1.2636195988860526 | 1.979578427252387e-199 | Historical / bridging; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| A_prime | case_control / PRS311_z_x_burden_z | 0.9691486851958342 | 0.9556221214686736–0.982866713647549 | 1.24351060421674e-05 | Historical / bridging; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| A_prime | diagnosis_age / PRS311_z | -0.4244807036256789 | -0.6762277489363634–-0.17273365831499443 | 0.0009505328123686846 | Historical / bridging; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| A_prime | diagnosis_age / burden_z | -1.713487902611176 | -2.0007417596006842–-1.4262340456216678 | 1.4118784122796467e-31 | Historical / bridging; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| A_prime_QC | case_control / PRS311_z | 1.5523740070126213 | 1.516701457440571–1.5888855686306664 | 6.428963123748066e-301 | Historical / bridging; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| A_prime_QC | case_control / burden_z | 1.2513528463024852 | 1.233621033210443–1.2693395328014059 | 3.168045478376163e-208 | Historical / bridging; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| A_prime_QC | case_control / PRS311_z_x_burden_z | 0.9698352814103133 | 0.9563737044792363–0.9834863387219387 | 1.7478185320841736e-05 | Historical / bridging; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| A_prime_QC | diagnosis_age / PRS311_z | -0.4251060052767106 | -0.6769418666113994–-0.17327014394202178 | 0.0009380581405760604 | Historical / bridging; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| A_prime_QC | diagnosis_age / burden_z | -1.720080526955623 | -2.0069626657978548–-1.4331983881133912 | 6.937505728470362e-32 | Historical / bridging; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| B | case_control / PRS311_z | 1.5542707655020616 | 1.5185113350800916–1.590872294915697 | 7.627082576302528e-302 | Historical / bridging; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| B | case_control / burden_z | 1.2652534014386052 | 1.2473022537856382–1.2834629016288825 | 1.816860917748464e-228 | Historical / bridging; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| B | case_control / PRS311_z_x_burden_z | 0.9703686864305325 | 0.9567985844628978–0.9841312507098829 | 2.837430186874043e-05 | Historical / bridging; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| B | diagnosis_age / PRS311_z | -0.4327631505106102 | -0.6845134479651622–-0.18101285305605824 | 0.0007538335261639913 | Historical / bridging; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| B | diagnosis_age / burden_z | -1.778380987932471 | -2.064662666554008–-1.4920993093109338 | 4.208913294957709e-34 | Historical / bridging; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| C_without_STK11 | case_control / PRS311_z | 1.5540909142257664 | 1.518324596634916–1.5906997588209508 | 1.6841827277569018e-301 | Sensitivity; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| C_without_STK11 | case_control / burden_z | 1.2683337572240925 | 1.2503167892426563–1.286610348317077 | 1.3407733878237123e-232 | Sensitivity; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| C_without_STK11 | case_control / PRS311_z_x_burden_z | 0.9712647222211264 | 0.9576663570674085–0.9850561771012288 | 5.057404702284176e-05 | Sensitivity; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| C_without_STK11 | diagnosis_age / PRS311_z | -0.43369009214787596 | -0.6853860651317742–-0.18199411916397773 | 0.0007323908106446805 | Sensitivity; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| C_without_STK11 | diagnosis_age / burden_z | -1.8189953904766467 | -2.105706614170182–-1.5322841667831115 | 1.69394709605489e-35 | Sensitivity; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| C_PTV_only | case_control / PRS311_z | 1.552194210885911 | 1.5165176179943–1.5887101077626857 | 1.27028631972229e-300 | Sensitivity; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| C_PTV_only | case_control / burden_z | 1.2538439159090493 | 1.2360520385953138–1.2718918915814807 | 2.606171509405641e-211 | Sensitivity; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| C_PTV_only | case_control / PRS311_z_x_burden_z | 0.9705506967418608 | 0.9570731154256389–0.984218070452423 | 2.794641706034007e-05 | Sensitivity; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| C_PTV_only | diagnosis_age / PRS311_z | -0.42560368499115364 | -0.67743784186604–-0.17376952811626734 | 0.0009250912453293374 | Sensitivity; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| C_PTV_only | diagnosis_age / burden_z | -1.7545533717796968 | -2.0419234453805384–-1.467183298178855 | 5.312116243983067e-33 | Sensitivity; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| C_TP53_NF1_strict_VAF | case_control / PRS311_z | 1.5539635416218829 | 1.5182016345497047–1.5905678361400597 | 2.0541250380956022e-301 | Sensitivity; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| C_TP53_NF1_strict_VAF | case_control / burden_z | 1.2685406348528712 | 1.2505202751008238–1.2868206732139416 | 6.647485572396799e-233 | Sensitivity; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| C_TP53_NF1_strict_VAF | case_control / PRS311_z_x_burden_z | 0.9715539451805283 | 0.9579518874602144–0.9853491399222821 | 6.028909767603261e-05 | Sensitivity; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| C_TP53_NF1_strict_VAF | diagnosis_age / PRS311_z | -0.4328693828111374 | -0.6845674981904504–-0.1811712674318245 | 0.0007496712032446549 | Sensitivity; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| C_TP53_NF1_strict_VAF | diagnosis_age / burden_z | -1.8197988475987579 | -2.1064688042829327–-1.5331288909145828 | 1.546341227737796e-35 | Sensitivity; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| C_ClinVar_min1star | case_control / PRS311_z | 1.5536292622452699 | 1.5178752045172892–1.5902255187523147 | 3.999391572363158e-301 | Sensitivity; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| C_ClinVar_min1star | case_control / burden_z | 1.2687180489259458 | 1.2506769850105617–1.2870193558865766 | 1.0473461604674618e-232 | Sensitivity; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| C_ClinVar_min1star | case_control / PRS311_z_x_burden_z | 0.9722432209969274 | 0.9585905383555913–0.9860903513569159 | 9.569288475435708e-05 | Sensitivity; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| C_ClinVar_min1star | diagnosis_age / PRS311_z | -0.43073029675871133 | -0.682413619468685–-0.17904697404873776 | 0.0007957233331582761 | Sensitivity; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| C_ClinVar_min1star | diagnosis_age / burden_z | -1.8236024131206408 | -2.109916966701695–-1.5372878595395867 | 9.190931166640396e-36 | Sensitivity; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| C_without_CDH1 | case_control / PRS311_z | 1.5540276690635468 | 1.5182652922159872–1.5906324201683224 | 1.7360498249662795e-301 | Sensitivity; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| C_without_CDH1 | case_control / burden_z | 1.2678956226607259 | 1.249887380426903–1.2861633257015228 | 5.401887834192116e-232 | Sensitivity; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| C_without_CDH1 | case_control / PRS311_z_x_burden_z | 0.971092650470402 | 0.9575064861371672–0.9848715903764004 | 4.492984086232031e-05 | Sensitivity; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| C_without_CDH1 | diagnosis_age / PRS311_z | -0.43283718937525506 | -0.6845786622652918–-0.18109571648521827 | 0.0007519361014363034 | Sensitivity; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |
| C_without_CDH1 | diagnosis_age / burden_z | -1.8084425764993965 | -2.0951480942996326–-1.5217370586991605 | 4.154418344898847e-35 | Sensitivity; no / M | [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) |

## Sensitivity / historical CV AUC

| version | AUC_mean | AUC_min | AUC_max | AUC_repeat1 | AUC_repeat1_CI_low | AUC_repeat1_CI_high | Status / source |
| --- | --- | --- | --- | --- | --- | --- | --- |
| A | 0.7795928541334483 | 0.7794737619079322 | 0.7797272845819383 | 0.7796396412550269 | 0.7752223961234227 | 0.7840568863866314 | Historical / bridging; not primary; M; 08 → 09; [comparisons/all_versions_CV_AUC_summary.csv](manuscript/datas/comparisons/all_versions_CV_AUC_summary.csv) |
| A_prime | 0.7796746546602801 | 0.779558657459889 | 0.7798073492114436 | 0.7797211694790402 | 0.7753038254634523 | 0.7841385134946279 | Historical / bridging; not primary; M; 08 → 09; [comparisons/all_versions_CV_AUC_summary.csv](manuscript/datas/comparisons/all_versions_CV_AUC_summary.csv) |
| A_prime_QC | 0.7800132236647592 | 0.7798864513684737 | 0.7801482111806515 | 0.7800536776467458 | 0.7756354540251132 | 0.7844719012683786 | Historical / bridging; not primary; M; 08 → 09; [comparisons/all_versions_CV_AUC_summary.csv](manuscript/datas/comparisons/all_versions_CV_AUC_summary.csv) |
| B | 0.7812149895997116 | 0.7810803052317155 | 0.7813484152410783 | 0.7812552136135276 | 0.776850690493457 | 0.7856597367335977 | Historical / bridging; not primary; M; 08 → 09; [comparisons/all_versions_CV_AUC_summary.csv](manuscript/datas/comparisons/all_versions_CV_AUC_summary.csv) |
| C_without_STK11 | 0.7814929649077551 | 0.7813568671458357 | 0.7816278265686607 | 0.7815316853787106 | 0.7771316067262372 | 0.7859317640311841 | Sensitivity; not primary; M; 08 → 09; [comparisons/all_versions_CV_AUC_summary.csv](manuscript/datas/comparisons/all_versions_CV_AUC_summary.csv) |
| C_PTV_only | 0.7802086531263764 | 0.7800818916126108 | 0.7803438776809446 | 0.7802493202840849 | 0.7758352587914652 | 0.7846633817767048 | Sensitivity; not primary; M; 08 → 09; [comparisons/all_versions_CV_AUC_summary.csv](manuscript/datas/comparisons/all_versions_CV_AUC_summary.csv) |
| C_TP53_NF1_strict_VAF | 0.7815079047077796 | 0.7813720466359754 | 0.7816425234968433 | 0.7815465732105243 | 0.77714685610168 | 0.7859462903193689 | Sensitivity; not primary; M; 08 → 09; [comparisons/all_versions_CV_AUC_summary.csv](manuscript/datas/comparisons/all_versions_CV_AUC_summary.csv) |
| C_ClinVar_min1star | 0.7814957353083252 | 0.781359310889081 | 0.7816298982265873 | 0.7815350942458703 | 0.7771356686194836 | 0.7859345198722575 | Sensitivity; not primary; M; 08 → 09; [comparisons/all_versions_CV_AUC_summary.csv](manuscript/datas/comparisons/all_versions_CV_AUC_summary.csv) |
| C_without_CDH1 | 0.7814267189613462 | 0.781290141023765 | 0.7815626179536974 | 0.7814676239531623 | 0.7770670282189389 | 0.7858682196873856 | Sensitivity; not primary; M; 08 → 09; [comparisons/all_versions_CV_AUC_summary.csv](manuscript/datas/comparisons/all_versions_CV_AUC_summary.csv) |

## S7 — gene × PRS post-hoc interactions

Model INT-GENE: PRS_z + gene indicators + PRS×gene terms + the CC clinical/technical covariates; HC1; gene-test BH correction. Generated by [11 --task audit](scripts/11_interaction_analysis.py). All rows Post hoc / **not primary**; original source `results/table/59_interaction_audit_gene_interactions.csv`. These ORs are interaction modifiers, not gene main effects.

| Gene selector | Interaction OR | 95% CI | P | BH q | Source |
| --- | --- | --- | --- | --- | --- |
| BRCA1 | 0.6524532655051871 | 0.5144231502594864–0.827519646916458 | 0.0004298047343432 | 0.0030086331404024 | [interaction/gene_specific_interaction_OR.csv](manuscript/datas/interaction/gene_specific_interaction_OR.csv); L |
| BRCA2 | 0.7862621895358958 | 0.6403699588709921–0.965392305072706 | 0.021658643586355 | 0.0758052525522425 | [interaction/gene_specific_interaction_OR.csv](manuscript/datas/interaction/gene_specific_interaction_OR.csv); L |
| PALB2 | 1.0272849754417517 | 0.7477425093966275–1.4113339919913344 | 0.8680658813029853 | 0.8680658813029853 | [interaction/gene_specific_interaction_OR.csv](manuscript/datas/interaction/gene_specific_interaction_OR.csv); L |
| CHEK2 | 0.8803562852961913 | 0.697273936133034–1.1115103389045078 | 0.2840661438631037 | 0.3976926014083451 | [interaction/gene_specific_interaction_OR.csv](manuscript/datas/interaction/gene_specific_interaction_OR.csv); L |
| ATM | 1.037306710153589 | 0.7877085127402256–1.3659941380937095 | 0.7942399559144684 | 0.8680658813029853 | [interaction/gene_specific_interaction_OR.csv](manuscript/datas/interaction/gene_specific_interaction_OR.csv); L |
| BARD1 | 0.7413163969864707 | 0.4382590360977448–1.2539387786141083 | 0.2643548476751244 | 0.3976926014083451 | [interaction/gene_specific_interaction_OR.csv](manuscript/datas/interaction/gene_specific_interaction_OR.csv); L |
| Other_pooled_panel_genes | 1.2836523036277308 | 0.9477019681016572–1.7386934838907384 | 0.106746273863063 | 0.2490746390138138 | [interaction/gene_specific_interaction_OR.csv](manuscript/datas/interaction/gene_specific_interaction_OR.csv); L |

## S8 — post-hoc interaction robustness

These model labels define distinct post-hoc analyses, not the C eligibility sensitivities above. Generator 11; original `results/table/59_interaction_audit_models.csv`. Counts remain suppressed as in source.

| Model selector | Term | OR | 95% CI | P | Scaling | Status / source |
| --- | --- | --- | --- | --- | --- | --- |
| C_baseline_reproduction | PRS_interaction | 0.9715467971747972 | 0.9579441513245666–0.985342598308528 | 6.0080422227631e-05 | full_frozen_cohort_ddof1 | Post hoc; not primary; L; [interaction/interaction_sensitivity_models.csv](manuscript/datas/interaction/interaction_sensitivity_models.csv) |
| exclude_top100_interaction_influence | PRS_interaction | 0.9791026177553492 | 0.9659049185069484–0.992480644551642 | 0.0022881802664812 | full_frozen_cohort_ddof1 | Post hoc; not primary; L; [interaction/interaction_sensitivity_models.csv](manuscript/datas/interaction/interaction_sensitivity_models.csv) |
| exclude_top100_design_leverage | PRS_interaction | 0.9730897786171088 | 0.956230852591589–0.9902459376653484 | 0.0022192359543678 | full_frozen_cohort_ddof1 | Post hoc; not primary; L; [interaction/interaction_sensitivity_models.csv](manuscript/datas/interaction/interaction_sensitivity_models.csv) |
| exclude_multigene_carriers | PRS_interaction | 0.9719284060492156 | 0.9581736609314706–0.9858806028617506 | 9.02625406032803e-05 | full_frozen_cohort_ddof1 | Post hoc; not primary; L; [interaction/interaction_sensitivity_models.csv](manuscript/datas/interaction/interaction_sensitivity_models.csv) |
| binary_C_carrier | PRS_interaction | 0.8226945507346952 | 0.7437287768634541–0.910044555036526 | 0.0001501476690927 | full_frozen_cohort_ddof1 | Post hoc; not primary; L; [interaction/interaction_sensitivity_models.csv](manuscript/datas/interaction/interaction_sensitivity_models.csv) |
| linear_interaction_with_nonlinear_PRS_terms | PRS_interaction | 0.9714446159943192 | 0.9575020985395272–0.9855901552422476 | 8.571358095976343e-05 | full_frozen_cohort_ddof1 | Post hoc; not primary; L; [interaction/interaction_sensitivity_models.csv](manuscript/datas/interaction/interaction_sensitivity_models.csv) |
| nonlinear_last_observation_age | PRS_interaction | 0.9721298727412708 | 0.9580821414573388–0.9863835767134372 | 0.0001412182947286 | full_frozen_cohort_ddof1 | Post hoc; not primary; L; [interaction/interaction_sensitivity_models.csv](manuscript/datas/interaction/interaction_sensitivity_models.csv) |
| enrollment_age_replaces_last_observation_age | PRS_interaction | 0.9711855391695908 | 0.9577896072334542–0.9847688306167124 | 3.694395853591914e-05 | full_frozen_cohort_ddof1 | Post hoc; not primary; L; [interaction/interaction_sensitivity_models.csv](manuscript/datas/interaction/interaction_sensitivity_models.csv) |
| enrollment_age_50_to74_baseline_covariates | PRS_interaction | 0.9737593253604994 | 0.9563563838949908–0.9914789504146284 | 0.0038519652330259 | full_frozen_cohort_ddof1 | Post hoc; not primary; L; [interaction/interaction_sensitivity_models.csv](manuscript/datas/interaction/interaction_sensitivity_models.csv) |
| enrollment_age_50_to74_nonlinear_enrollment_age | PRS_interaction | 0.9731233588235776 | 0.9559832635605058–0.9905707637193863 | 0.0026568475863228 | full_frozen_cohort_ddof1 | Post hoc; not primary; L; [interaction/interaction_sensitivity_models.csv](manuscript/datas/interaction/interaction_sensitivity_models.csv) |
| EHR_complete_metadata_baseline | PRS_interaction | 0.9715465814671148 | 0.957943963354116–0.9853423541137885 | 6.007069620465407e-05 | full_frozen_cohort_ddof1 | Post hoc; not primary; L; [interaction/interaction_sensitivity_models.csv](manuscript/datas/interaction/interaction_sensitivity_models.csv) |
| EHR_log_span_and_date_count_adjustment | PRS_interaction | 0.9714793026713756 | 0.9561810555755628–0.9870223113244686 | 0.0003529921994871 | full_frozen_cohort_ddof1 | Post hoc; not primary; L; [interaction/interaction_sensitivity_models.csv](manuscript/datas/interaction/interaction_sensitivity_models.csv) |
| EHR_span_ge365days_dates_ge5_baseline | PRS_interaction | 0.9719823701633766 | 0.9571376198077228–0.9870573555537432 | 0.0002957925015175 | full_frozen_cohort_ddof1 | Post hoc; not primary; L; [interaction/interaction_sensitivity_models.csv](manuscript/datas/interaction/interaction_sensitivity_models.csv) |
| EHR_span_ge365days_dates_ge5_age_spline_EHR_adjusted | PRS_interaction | 0.9734382389232704 | 0.9578904142864544–0.989238425257554 | 0.0010489289748491 | full_frozen_cohort_ddof1 | Post hoc; not primary; L; [interaction/interaction_sensitivity_models.csv](manuscript/datas/interaction/interaction_sensitivity_models.csv) |
| exclude_BRCA1_carriers | PRS_interaction | 0.9795820566357272 | 0.9618730046119056–0.9976171501661496 | 0.0266732628356826 | full_frozen_cohort_ddof1 | Post hoc; not primary; L; [interaction/interaction_sensitivity_models.csv](manuscript/datas/interaction/interaction_sensitivity_models.csv) |
| exclude_BRCA2_carriers | PRS_interaction | 0.9723441491307436 | 0.9555586111432516–0.9894245453113848 | 0.0015961927812799 | full_frozen_cohort_ddof1 | Post hoc; not primary; L; [interaction/interaction_sensitivity_models.csv](manuscript/datas/interaction/interaction_sensitivity_models.csv) |
| exclude_PALB2_carriers | PRS_interaction | 0.9689040815692768 | 0.9546361347181171–0.983385276484221 | 3.000260677562338e-05 | full_frozen_cohort_ddof1 | Post hoc; not primary; L; [interaction/interaction_sensitivity_models.csv](manuscript/datas/interaction/interaction_sensitivity_models.csv) |
| exclude_CHEK2_carriers | PRS_interaction | 0.9729489127137338 | 0.9588509898059252–0.9872541164528994 | 0.0002309522083191 | full_frozen_cohort_ddof1 | Post hoc; not primary; L; [interaction/interaction_sensitivity_models.csv](manuscript/datas/interaction/interaction_sensitivity_models.csv) |
| exclude_ATM_carriers | PRS_interaction | 0.9699424654955032 | 0.9561334643267372–0.9839509037934924 | 3.0268964287701296e-05 | full_frozen_cohort_ddof1 | Post hoc; not primary; L; [interaction/interaction_sensitivity_models.csv](manuscript/datas/interaction/interaction_sensitivity_models.csv) |
| exclude_BARD1_carriers | PRS_interaction | 0.9719305667676869 | 0.9581475112315652–0.9859118930478068 | 9.345601031896548e-05 | full_frozen_cohort_ddof1 | Post hoc; not primary; L; [interaction/interaction_sensitivity_models.csv](manuscript/datas/interaction/interaction_sensitivity_models.csv) |

## Exploratory age slopes

Binary-carrier age interaction / contrast model; not AGE-JOINT burden coefficients. Generator [13 age-figures](scripts/13_figures.py); historical `results/table/48_age_binary_carrier_PRS_slopes.csv`.

| Contrast selector | Years per PRS SD | HC1 CI | Status / source |
| --- | --- | --- | --- |
| Non-carrier PRS slope | -0.4678794632921871 | -0.7275021803402821–-0.20825674624409202 | Exploratory; not primary; M; [age/carrier_PRS_slope_contrasts.csv](manuscript/datas/age/carrier_PRS_slope_contrasts.csv) |
| Carrier PRS slope | 0.5120958671866677 | -0.5178286307055683–1.5420203650789035 | Exploratory; not primary; M; [age/carrier_PRS_slope_contrasts.csv](manuscript/datas/age/carrier_PRS_slope_contrasts.csv) |
| Carrier minus non-carrier slope (interaction) | 0.9799753304788548 | -0.08215427073051251–2.042104931688222 | Exploratory; not primary; M; [age/carrier_PRS_slope_contrasts.csv](manuscript/datas/age/carrier_PRS_slope_contrasts.csv) |

## Coverage gaps — do not fill from examples or old numbers

- Manuscript carrier count 3,209 / 2.3%, age-group carrier frequencies 26.5% vs 3.2%, Table 1 clinical characteristics, and other figure-specific contrasts remain **manuscript-stated, source reconciliation pending**. Do not treat them as newly verified numerical results.
- PRS direct-only / no-imputation sensitivity results, reviewer Cox/frailty/family results, historical supplementary S2–S6, author notes, and old calibration/category/ventile claims require separate current-source mapping before reuse.
- Current DOCX tracked-change accepted view and embedded figure data have not been fully reconciled. This registry is an audited baseline, not certification that every manuscript number is current.
- Further comparison delta-AUC / age-R² values are available in the checksum-verified all-version sources below. Preserve version selectors; they are not main C results.

## Source integrity ledger

Each entry lists the existing copy, original recorded source, and SHA-256 checked on 2026-10-15. These checks establish copy integrity only. Historical source paths may be absent. New executions must retain a separate run/source/input lock; never overwrite old hashes to make changed data look verified.

| Existing source copy | Recorded origin / counterpart | SHA-256 | Evidence |
| --- | --- | --- | --- |
| [burden/gene_OR_weights.csv](manuscript/datas/burden/gene_OR_weights.csv) | results/table/29_rowlands_13gene_external_weights.csv | 1292ce27c1e55fa955fe62764dd5ccb5ea18248b6bada189f3e30e19d3359cda | M; matches manifest |
| [case_control/joint_model_coefficients.csv](manuscript/datas/case_control/joint_model_coefficients.csv) | results/table/45_rowlands_fixed_cohort_coefficients.csv | 1e5728e05daa285b0e4382b69d4a3c23a85e1fc4d6a767d4c07dc6dcd71406bc | M; matches manifest |
| [age/age_model_coefficients.csv](manuscript/datas/age/age_model_coefficients.csv) | results/table/45_rowlands_fixed_cohort_coefficients.csv | 413ea0b67a3a567a11169f255867e8a4182d7ffde11763926ae8f751d18583ad | M; matches manifest |
| [case_control/CV_AUC_summary.csv](manuscript/datas/case_control/CV_AUC_summary.csv) | results/table/45_rowlands_CV_AUC.csv | 85b7c9293e6c3692ca0f8e3a280f29c19119f8c2186a7cce42f49c8060e427bc | M; matches manifest |
| [case_control/CV_delta_AUC_summary.csv](manuscript/datas/case_control/CV_delta_AUC_summary.csv) | results/table/45_rowlands_CV_delta_AUC.csv | c3b649d152e479448dfb4c66d762065dbb007a3863daad8b4431a945422ffc31 | M; matches manifest |
| [age/age_R_squared.csv](manuscript/datas/age/age_R_squared.csv) | results/table/45_rowlands_age_explained_variance.csv | 23e0ada99a6f7f9bd7ce0c83579b1c2519b4f3f21a0a14470fea1999e9c3f558 | M; matches manifest |
| [age/carrier_PRS_slope_contrasts.csv](manuscript/datas/age/carrier_PRS_slope_contrasts.csv) | results/table/48_age_binary_carrier_PRS_slopes.csv | 44fb3e91b7b78d4de6b439c0ffae13c32b7dc75b8032c44a44ab40fe1bbbbdc0 | M; matches manifest |
| [burden/score_distribution_disclosure_coarsened.csv](manuscript/datas/burden/score_distribution_disclosure_coarsened.csv) | results/table/48_burden_distribution_carrier_supplement_coarsened.csv | 2bcd76fe34f58dcf112c1b23c059088b5538fe80a09673c8beac237268dc4a14 | M; matches manifest |
| [comparisons/all_versions_coefficients.csv](manuscript/datas/comparisons/all_versions_coefficients.csv) | results/table/45_rowlands_fixed_cohort_coefficients.csv | 20b540ff66e2320273bfd969e8425b6e3e964389813be6f5c901c6629f77b941 | M; matches manifest |
| [comparisons/all_versions_CV_AUC_summary.csv](manuscript/datas/comparisons/all_versions_CV_AUC_summary.csv) | results/table/45_rowlands_CV_AUC.csv | ceed71931197d48b1cccd047c06ae43025ef970270395e30e77248065b8b9e60 | M; matches manifest |
| [comparisons/all_versions_CV_delta_AUC_summary.csv](manuscript/datas/comparisons/all_versions_CV_delta_AUC_summary.csv) | results/table/45_rowlands_CV_delta_AUC.csv | 98cdf0ce1da5fd9d80097c7cef6913b8d4bedf5943f62d6c88846dd9c58dd85f | M; matches manifest |
| [comparisons/all_versions_age_R_squared.csv](manuscript/datas/comparisons/all_versions_age_R_squared.csv) | results/table/45_rowlands_age_explained_variance.csv | f5e19c94e20815218b664941e8826f412ff915e982eebd3c1af99fdbb631eb13 | M; matches manifest |
| [case_control/CV_AUC_by_repeat.csv](manuscript/datas/case_control/CV_AUC_by_repeat.csv) | data/29_rowlands_extension/section11_models/CV_AUC_by_repeat_internal.csv | cd255174324cfba9f154ec81bc0859a1355eb8d7392e4b813324956525e8a933 | M; matches manifest |
| [case_control/CV_delta_AUC_by_repeat.csv](manuscript/datas/case_control/CV_delta_AUC_by_repeat.csv) | data/29_rowlands_extension/section11_models/CV_paired_deltas_by_repeat_internal.csv | 53a6548d9189d9fbdbee6fe17e1008812ccd108181435e0011b9d828959ef337 | M; matches manifest |
| [case_control/four_models_AUC_summary.csv](manuscript/datas/case_control/four_models_AUC_summary.csv) | scripts/14_four_model_discrimination.py | 6eb418a9c21118ac447c2f307043baa4641dc595e32004ec79707ca2abcc3d6a | M; matches manifest |
| [case_control/four_models_AUC_by_repeat.csv](manuscript/datas/case_control/four_models_AUC_by_repeat.csv) | scripts/14_four_model_discrimination.py | f239772a455627c5c2cdee5f5b0f32f9b7b6223fd02f41397a74db38b7f7f5bf | M; matches manifest |
| [case_control/four_models_paired_delta_AUC_by_repeat.csv](manuscript/datas/case_control/four_models_paired_delta_AUC_by_repeat.csv) | scripts/14_four_model_discrimination.py | df4a56c6c6e727dfc816d9b2046cbc928c79cb4f7d06402fb8f73e28d963e7cb | M; matches manifest |
| [case_control/four_models_ROC_repeat1_interpolated.csv](manuscript/datas/case_control/four_models_ROC_repeat1_interpolated.csv) | scripts/14_four_model_discrimination.py | 93f4a2bed633ecde4d4f4e7f65e46f14f89d547f9c9b867573ca172f3e8d61f2 | M; matches manifest |
| [case_control/PRS311_decile_C_carrier_joint_distribution.csv](manuscript/datas/case_control/PRS311_decile_C_carrier_joint_distribution.csv) | scripts/15_joint_distribution_heatmap.py | e2db322f7e42200d0dca3236205a617f33c2136093366e8bb6665538c145ba71 | M; matches manifest |
| [case_control/joint_distribution_no_interaction_coefficients.csv](manuscript/datas/case_control/joint_distribution_no_interaction_coefficients.csv) | scripts/15_joint_distribution_heatmap.py | 2152bf180f8222386bdc1a3b341a5669a3e54b7debb21bd1550ada5de8c21ca2 | M; matches manifest |
| [burden/13gene_burden_histogram_percentages.csv](manuscript/datas/burden/13gene_burden_histogram_percentages.csv) | scripts/17_burden_histogram.py | a3530d9e171206d70ea90cda121a2cd207b34ee4e15b6858325443da66797571 | M; matches manifest |
| [interaction/gene_specific_interaction_OR.csv](manuscript/datas/interaction/gene_specific_interaction_OR.csv) | manuscript/tables/tableS7_gene_PRS_interactions.csv | cca78ebe805558639bf02c9c8f89369a72297c960f39803f81896a4525d2ecad | L; S7/S8 differs — reconcile |
| [interaction/interaction_sensitivity_models.csv](manuscript/datas/interaction/interaction_sensitivity_models.csv) | manuscript/tables/tableS8_interaction_sensitivities.csv | 06a0288273f7bce36c87630539f36c85d19f90a33d43defcfe5827215f79d657 | L; S7/S8 differs — reconcile |

## Maintenance contract

Use stable result IDs plus exact `version/cohort/term/model/comparison` selectors. For every replacement record value, unit, CI/test, cohort N, model specification, generating script/task, output, checksum, manuscript destination, role, and reason/date in CHANGELOG. Preserve superseded numbers with their historical identity. Never infer primary status from filenames. Never copy participant-level rows or unsuppress source cells into this registry.


## AGE-QUART-20261016-01 / 2026-10-16

Real frozen-input Q1/Q4 descriptives completed (not synthetic-only), shared quartiles in original7,996 cases; selected3,998; early<45 fixed; existing13-gene carrier. Source/command/software/input and output hashes/status: docs/age_interaction_runs/AGE-QUART-20261016-01/provenance.json; [report](docs/age_interaction_runs/AGE-QUART-20261016-01/REPORT.md). No participant records/scores written; controlled logs and snapshot ignored. No regression or significance test; no failure/blocked stage. Supersedes=None; prior median split, additive and continuous-interaction results retained; manuscript unchanged.

Result selector: four_groups.csv/group; unit years and percent; early Wilson95%CI. SHA256=937fe98e82c84e28f710e456d00d1da6905a5545d2975f1402f4bc99cf8c57af. Exact aggregate results:

| Group | n | median [Q1,Q3] years | early n / percent [95%CI] |
|---|---|---|---|
| non-carrier / PRS Q1 | 1821 | 62.24229979466119 [52.63244353182751,69.89185489390829] | 187 / 10.26908292147172 [8.95746356396704,11.74797673575492] |
| non-carrier / PRS Q4 | 1897 | 60.865160848733744 [52.42436687200548,68.26009582477755] | 193 / 10.173958882445968 [8.893057451612135,11.615831230906217] |
| carrier / PRS Q1 | 178 | 53.98083504449008 [43.12388774811773,65.07597535934292] | 51 / 28.651685393258425 [22.5156567575984,35.68969394026403] |
| carrier / PRS Q4 | 102 | 55.39356605065024 [46.969199178644764,63.67556468172485] | 21 / 20.58823529411765 [13.87888951082538,29.432549438966564] |

Contrasts: [{"comparison": "carrier / PRS Q1 minus carrier / PRS Q4", "median_difference_years": -1.4127310061601577, "early_percentage_point_difference": 8.063450099140777, "interpretation": "unadjusted descriptive difference; not an interaction or incidence-risk test"}, {"comparison": "carrier / PRS Q1 minus non-carrier / PRS Q4", "median_difference_years": -6.884325804243666, "early_percentage_point_difference": 18.47772651081246, "interpretation": "unadjusted descriptive difference; not an interaction or incidence-risk test"}]; comparisons.json SHA256=0ea10c7b71d4167f925d68f12c4d7f59080fee942e6596193ff277215138fa77. No contrast CI/P computed.


## FIGURES-PUBLISH-20261017-01 / 2026-10-17
Author requested publication of the latest figures. Updated root figures/ from the current manuscript figure set, with aggregate tables in tables/ and reproduction scripts in scripts/. Figure04a uses full-cohort (141204) PRS quartiles; Figure04c restores four carrier/PRS combinations using full-cohort median and bottom/top20% panels. Figure07b retains original four attained-age bands and is8in wide. Supplementary Figure07b uses seven age bands; two real Cox fits passed convergence, with35 HR/CI rows checked against stored coefficients and covariance. Red/blue1/blue2/orange series; gray theta0.73 line reuses the existing illustrative refitted scenario. Eight-band trials failed convergence and were not adopted. No participant-level data, individual scores, controlled logs, local archives or unrelated workspace changes are included. Original four-band scientific baseline unchanged; seven-band extension remains post hoc/provisional. Manifest: docs/figure_publication_20261017.json; detailed source-run evidence: docs/age_interaction_runs/. Published script output paths adapt to this repository layout; provenance describes the original controlled workspace. Controlled-data reproduction and embedded-document consistency are not claimed.

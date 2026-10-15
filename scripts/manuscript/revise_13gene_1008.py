"""Rebuild the 1006 manuscript using verified aggregate 13-gene outputs.

Original DOCX files and model inputs are preserved. No participant data are read.
"""

from pathlib import Path
import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.functions.paths import PROJECT_ROOT

import zipfile, hashlib, json, copy, difflib, re
import pandas as pd
from lxml import etree as E
import cairosvg
from PIL import Image

ROOT = PROJECT_ROOT
M = ROOT / 'manuscript'
OUT = M / 'revision_13gene_1008'
OUT.mkdir(exist_ok=True)
SOURCE = M / 'medrxiv_v1_main_revised_1006.docx'
W = 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'
R = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'
REL = 'http://schemas.openxmlformats.org/package/2006/relationships'
CT = 'http://schemas.openxmlformats.org/package/2006/content-types'


def q(s):
    return '{' + W + '}' + s


def xml(s):
    return E.tostring(s, encoding='UTF-8', xml_declaration=True, standalone=True)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


with zipfile.ZipFile(SOURCE) as z:
    original = {n: z.read(n) for n in z.namelist()}
original_root = E.fromstring(original['word/document.xml'])
original_body = original_root.find(q('body'))
ps = original_body.findall(q('p'))


def ptext(p):
    return ''.join(p.xpath('.//w:t/text()', namespaces={'w': W}))


source_text = '\n\n'.join(ptext(p) for p in ps if ptext(p)) + '\n'
(OUT / 'main_1006_source.txt').write_text(source_text)
inputs = [
    SOURCE,
    M / 'manuscript_revision_outline_13gene.md',
    M / 'datas/case_control/joint_model_coefficients.csv',
    M / 'datas/case_control/four_models_AUC_summary.csv',
    M / 'datas/case_control/CV_delta_AUC_summary.csv',
    M / 'datas/age/age_model_coefficients.csv',
    M / 'datas/age/age_R_squared.csv',
    M / 'datas/burden/gene_OR_weights.csv',
    ROOT / 'results/table/09b_supp_table1_cohort.csv',
    ROOT / 'results/table/56_figure5_primary_sensitivity_estimates.csv',
    M / 'datas/comparisons/all_versions_CV_AUC_summary.csv',
    M / 'datas/comparisons/all_versions_CV_delta_AUC_summary.csv',
    M / 'datas/case_control/four_models_ROC_repeat1_interpolated.csv',
    M / 'datas/case_control/PRS311_decile_C_carrier_joint_distribution.csv',
]
hashes = {str(p.relative_to(ROOT)): sha(p) for p in inputs}

main = []
supp = []


def add(h, t):
    main.extend([('heading', h), ('text', t)])


main += [
    (
        'title',
        'Common- and rare-variant genetic susceptibility to breast cancer: joint discrimination and recorded diagnosis age',
    ),
    ('text', ptext(ps[1])),
    ('text', ptext(ps[2])),
    ('text', 'Corresponding author: [name, address, e-mail]'),
]
add(
    'Abstract',
    'Common and rare genetic variants contribute to breast cancer susceptibility, but their joint association with case status and recorded diagnosis age requires evaluation in electronic health record cohorts. We studied 8,528 cases and 132,676 controls without recorded cancer of genetically inferred European ancestry in the All of Us Research Program. We combined a common-variant polygenic risk score (PRS), named PRS311 with 302 retained scored loci, and a thirteen-gene burden index based on qualifying protein-truncating variants (PTVs) and eligible pathogenic/likely pathogenic missense variants, weighted by external gene-level log odds ratios. Four logistic models included age, ancestry principal components and sequencing covariates; discrimination was evaluated using ten repeats of five-fold cross-validation. In the joint model, the odds ratios per standard deviation were 1.554 (95% confidence interval, 1.518–1.591) for PRS and 1.269 (1.251–1.287) for burden, conditional on the other predictor at its mean. Their multiplicative-odds interaction was negative (interaction odds ratio, 0.972; 0.958–0.985; P=6.0×10⁻⁵). Mean areas under the receiver operating characteristic curve were 0.7452 for covariates alone, 0.7682 with PRS, 0.7588 with burden, and 0.7815 for the joint model; the joint-versus-PRS increment was 0.01328. Among 7,996 cases with eligible diagnosis-age records, each standard-deviation increase in PRS and burden was associated with 0.43 (0.18–0.68) and 1.82 (1.53–2.11) years younger recorded diagnosis age, respectively. The age model explained 2.63% of variation. The thirteen-gene burden added modest discrimination beyond PRS and associated with younger recorded diagnosis age among cases. These findings do not establish biological antagonism, prospective risk calibration or clinical utility; variant adjudication remains incomplete.',
)
add(
    'Introduction',
    'Breast cancer susceptibility reflects both common genetic variation and rarer variants in susceptibility genes. Common-variant polygenic risk scores (PRSs) summarize many small associations, whereas gene-level studies identify larger associations for selected protein-truncating variants (PTVs) and pathogenic missense variants[1–3]. Variant class and gene context matter: including every missense variant would dilute the distinction between clinical pathogenicity and sequence variation. Rowlands et al. combined population-based evidence across thirteen genes and reported estimates for PTVs together with non-PTV pathogenic variants, providing an external basis for an expanded panel[13].',
)
main += [
    (
        'text',
        'Studies of pathogenic-variant carriers indicate that common genetic background can modify breast cancer susceptibility[4–8]. These observations motivate joint modeling, but improved discrimination does not by itself establish a statistical interaction, and a logistic interaction concerns multiplicative odds rather than an additive absolute-risk scale. Gene-specific differences may also be concealed by a pooled burden index.',
    ),
    (
        'text',
        'The All of Us Research Program links genomic data with electronic health records (EHRs), enabling assessment in a health care cohort[9,10]. We evaluated separate and joint associations of PRS and a thirteen-gene PTV plus eligible pathogenic/likely pathogenic missense burden with breast cancer case status, incremental discrimination beyond PRS and covariates, and recorded diagnosis age among cases. The main paper uses this single burden definition; historical score comparisons are reported separately in the supplement.',
    ),
]
add(
    'Methods',
    'Study population. We used controlled-tier All of Us data release C2025Q4R6. The analysis selected participants recorded as female, with genetically inferred European ancestry (maximum ancestry posterior probability >0.75), and retained an unrelated set using KING kinship >0.1 to identify related pairs[11]. Participants with missing genotype rates >3% were excluded. Cohort selection and the distinction between the descriptive and fitted-model cohorts are shown in Figure 1 and Table 1.',
)
main += [
    (
        'text',
        'Phenotypes. Cases had an EHR breast cancer diagnosis within the SNOMED breast cancer concept 4112853 and its descendants. Controls had neither a coded breast cancer diagnosis nor a coded diagnosis in the project’s any-cancer exclusion set. Controls were not required to be free of benign breast disease, other illnesses or breast cancer family history. Absence of a recorded diagnosis does not confirm absence of disease. Recorded diagnosis age was age at the earliest qualifying breast cancer code. The observation-age covariate was the later of the last clinical EHR date and death date, where available. We excluded 2,051 controls with missing observation age from case-status models. The diagnosis-age analysis excluded 532 cases whose first breast cancer record coincided with their first EHR record; this rule does not establish incident disease.',
    ),
    (
        'text',
        'Common-variant score. We implemented the overall breast cancer score PGS000004 from the Polygenic Score Catalog[1,12] on GRCh38. Of 313 intended loci, 304 were resolved (281 direct matches and 23 frequency-resolved matches); eight were unmatched and one remained unresolved. Excluding two CHEK2 tagging loci left 302 scored loci in the score named PRS311. Allele dosages were multiplied by the published weights and summed; missing dosage was replaced using the cohort mean dosage. The name denotes the intended score after tag exclusions, not the number of successfully scored loci. Harmonization details are in Supplementary Methods.',
    ),
    (
        'text',
        'Thirteen-gene panel and variant eligibility. The panel comprised ATM, BARD1, BRCA1, BRCA2, CHEK2, PALB2, RAD51C, RAD51D, TP53, CDH1, NF1, PTEN and STK11[13]. A PTV required a stop-gained, frameshift, splice-donor or splice-acceptor consequence and a high-confidence loss-of-function annotation from the Loss-of-Function Transcript Effect Estimator (LOFTEE)[14] on the same retained transcript of the same gene. Any qualifying retained transcript could establish eligibility; the analysis was not restricted to canonical transcripts and did not implement a blanket final-exon exclusion. Eligible missense variants required germline pathogenic or likely pathogenic classification, ClinVar review status of at least two stars, and applicable disease and inheritance evidence. Low-penetrance, risk-factor/risk-allele, conflicting or opposing benign assertions, and records held out for unresolved eligibility were excluded. ClinVar release 2026-09-13 was used. Clinical classification and loss-of-function interpretation require disease and transcript context[16,18]; detailed rules are provided in Supplementary Methods.',
    ),
    (
        'text',
        'Genotype quality control. Both variant classes required a known, consistent All of Us European allele frequency <0.01 and a valid normalized-allele match. Original variant-dataset site filters had to be defined-empty or singleton PASS; genotype filters had to use the accepted PASS or “.” encoding. Genotype quality was ≥20. The target alternative-allele fraction, calculated as its allele depth divided by the sum of all allelic depths, had to be 0.30–0.80 for heterozygotes and ≥0.80 for homozygous alternative genotypes. Missing required quality fields did not pass. The analysis did not adopt a separate depth ≥10 requirement: the queried genotype extract did not provide the corresponding depth field, and summed allelic depth was not treated as that field.',
    ),
    (
        'text',
        'Burden calculation. For participant i, Bᵢ = Σg∈G₁₃ ln(ORg) Cᵢg, where Cᵢg=1 for at least one qualifying PTV or eligible pathogenic/likely pathogenic missense variant in gene g and 0 otherwise. Multiple qualifying variants within a gene contributed once; variants in different genes contributed separately. We used Rowlands Table 2 combined PTV plus non-PTV pathogenic-variant odds ratios as approximate external gene-level weights[13] (Supplementary Table S1). Those source variant definitions differ from ours, so the sum is a burden index rather than an individually calibrated risk estimate. No allele- or variant-count dosage effect was assumed.',
    ),
    (
        'text',
        'Case-status analysis. Logistic regression compared covariates alone, PRS plus covariates, burden plus covariates, and a joint model with both predictors and their product. Covariates were observation age, ancestry principal components (PCs) 1–10, sequencing center and sequencer. Genetic predictors were standardized using means and sample standard deviations (SDs) in the fitted-model cohort. The joint model was logit Pr(caseᵢ)=α+βP Pᵢ+βB Bᵢ+βPB PᵢBᵢ+γᵀXᵢ, where P and B are standardized predictors. Main genetic effects are conditional on the other predictor being zero, its cohort mean. HC1 heteroscedasticity-robust standard errors were used; odds ratios (ORs) are exponentiated coefficients with 95% confidence intervals (CIs).',
    ),
    (
        'text',
        'Discrimination. We used ten repeats of stratified five-fold cross-validation (CV), with score standardization learned in each training fold. The area under the receiver operating characteristic curve (AUC) was calculated from held-out predictions for each repeat. Mean repeated-CV AUC is reported separately from repeat-1 AUC and its confidence interval. Paired DeLong intervals and tests[15] used matched held-out predictions within a repeat. Averaging interval endpoints across repeats was not interpreted as a confidence interval for the mean. The principal discrimination contrast was the joint model versus PRS plus covariates; it does not isolate the contribution of the product term.',
    ),
    (
        'text',
        'Recorded diagnosis age. Ordinary least squares regression among eligible cases included continuous PRS, burden and PCs 1–10, with HC1 robust standard errors. Both genetic predictors were standardized within this case cohort using sample SDs. We reported coefficients in years per SD, R² and the increment in R² after adding burden to the PRS-plus-PC model. This analysis concerns recorded diagnosis age conditional on being a case, rather than prospective onset or time to disease.',
    ),
    (
        'text',
        'Sensitivity analyses and chronology. Score sensitivities excluded STK11, used PTVs alone, tightened the heterozygous upper allele-fraction bound to 0.70 for TP53 and NF1, allowed ClinVar one-star eligible missense variants, or excluded CDH1. These specifications were locked before evaluation of the current thirteen-gene associations, after earlier burden analyses had been conducted. Selection of this score as the sole main analysis does not imply original study-wide prespecification. The variant-definition records remain provisional pending final adjudication. A subsequent exploratory interaction audit assessed influence, binary carrier status, nonlinear PRS and age terms, gene heterogeneity, and EHR coverage. Seven gene-category interaction tests used Benjamini–Hochberg adjustment; other audit P values were nominal. Details appear in Supplementary Tables S7–S8.',
    ),
    (
        'text',
        'Reproducibility and disclosure. Analyses used Hail, PLINK 2[17], pandas, statsmodels and scikit-learn. Only aggregate outputs meeting project disclosure rules were included; small cells below 20 and secondary cells enabling their reconstruction were suppressed. We retained separate denominators for the descriptive cohort, fitted case-status cohort and diagnosis-age cohort.',
    ),
]
add(
    'Results',
    'Study cohort. Cohort selection began with 295,706 participants with whole-genome sequencing, retained 150,900 participants after ancestry and eligibility selection, and yielded 143,255 after cohort quality control: 8,528 cases and 134,727 controls. The case-status models included 141,204 participants after excluding 2,051 controls with missing observation age. Table 1 describes the larger cohort; 7,996 cases were eligible for diagnosis-age analysis. In the fitted case-status cohort, 3,209 participants (2.27%) carried a qualifying panel variant; 543 diagnosis-age cases (6.79%) were carriers. Cases had mean recorded diagnosis age 59.9 years (SD 11.8), and mean enrollment age was 64.0 years in cases versus 51.2 years in controls.',
)
main += [
    (
        'text',
        'Case-status associations. In the joint model, the OR per SD was 1.554 (95% CI, 1.518–1.591) for PRS and 1.269 (1.251–1.287) for burden, conditional on the other predictor at its mean. The interaction coefficient was −0.0289 (interaction OR, 0.972; 0.958–0.985; P=6.0×10⁻⁵), indicating a negative departure from multiplicative odds (Table 2). This coefficient does not show that carriers have a lower overall breast cancer association.',
    ),
    (
        'text',
        'Incremental discrimination. Mean repeated-CV AUCs were 0.7452 for covariates alone, 0.7682 for PRS plus covariates, 0.7588 for burden plus covariates, and 0.7815 for the joint model (Table 2; Figure 2). The mean joint-versus-PRS increment was 0.01328. For repeat 1, the paired increment was 0.01328 (95% CI, 0.01133–0.01523; P=1.68×10⁻⁴⁰). Repeat-specific results and historical score comparisons appear in Supplementary Tables S3–S5. These results quantify discrimination in the selected cohort and do not establish prospective absolute-risk calibration.',
    ),
    (
        'text',
        'Joint distribution. Observed case proportions were displayed across PRS deciles separately for carriers and noncarriers of qualifying thirteen-gene variants (Figure 3). This descriptive display uses case counts divided by participants within each decile–carrier cell. It is not an age-standardized risk estimate or a prospective penetrance curve; the marked age difference between cases and controls should be considered when interpreting it.',
    ),
    (
        'text',
        'Recorded diagnosis age. In 7,996 cases, PRS was associated with −0.43 years per SD (95% CI, −0.68 to −0.18; P=7.50×10⁻⁴) and burden with −1.82 years per SD (−2.11 to −1.53; P=1.55×10⁻³⁵). The joint age model had R²=0.02628 and adjusted R²=0.02481. Adding burden to PRS and PCs increased R² by 0.02366, or 2.37 percentage points (Table 3). Burden-group and diagnosis-age carrier distributions are shown in Figure 4. An exploratory binary-carrier-by-PRS age interaction was uncertain: its coefficient was 0.98 years (95% CI, −0.08 to 2.04), which does not establish synergy (Supplementary Table S6).',
    ),
    (
        'text',
        'Sensitivity analyses. Across the five score sensitivities, burden remained associated with case status and younger recorded diagnosis age (Supplementary Table S2). The exploratory interaction audit retained a negative panel-level estimate with binary carrier status (OR, 0.823; 95% CI, 0.744–0.910; P=1.50×10⁻⁴), nonlinear PRS terms (joint interaction P=0.00249), alternative age specifications and EHR coverage adjustment or restriction. Removing the 100 most interaction-influential observations attenuated the coefficient magnitude by approximately 27% (β=−0.0211; nominal P=0.00229). Gene-category interactions differed (heterogeneity P=0.0126); BRCA1 showed the clearest negative estimate (OR, 0.652; 95% CI, 0.514–0.828; adjusted q=0.00301), while BRCA2 did not pass correction (nominal P=0.0217; q=0.0758). Excluding BRCA1 carriers attenuated the pooled interaction (β=−0.0206; nominal P=0.0267). Full results and descriptive gene-specific PRS box plots are in Supplementary Tables S7–S8 and Figure S7.',
    ),
]
add(
    'Discussion',
    'In this EHR-linked cohort, a thirteen-gene PTV plus eligible pathogenic/likely pathogenic missense burden associated with breast cancer case status and added modest discrimination beyond PRS and covariates. Both predictors associated with younger recorded diagnosis age among cases, but the age model explained only 2.63% of variation. These findings support complementary information from common and qualifying rare variants without establishing that the score predicts calibrated individual risk or improves clinical outcomes.',
)
main += [
    (
        'text',
        'The burden is continuous because gene carriers receive different external log-OR weights and carriers in multiple genes receive summed contributions. It is not a count of variants within a gene. This distinction preserves externally supported differences between genes without assuming an unvalidated linear dosage effect. Binary carrier status provides a useful robustness analysis, but discards the gene-weight information. Rowlands estimates are approximate weights: its non-PTV pathogenic-variant set and source populations differ from our eligible missense set and selected cohort[13]. STK11 has an imprecise source estimate and CDH1 has subtype-specific associations; the exclusion sensitivities address, but do not eliminate, these concerns.',
    ),
    (
        'text',
        'The negative logistic interaction indicates a weaker joint association than expected under multiplicative odds; it does not demonstrate protection, biological antagonism or an additive absolute-risk interaction. Gene-specific results suggest that heterogeneous PRS associations, including a prominent BRCA1 contribution, may partly explain the pooled result. Barnes et al. found that an estrogen receptor-negative-specific PRS had the strongest breast cancer association in BRCA1 carriers, whereas an overall breast cancer PRS performed best in BRCA2 carriers; carrier-cohort associations were smaller than general-population estimates[5]. This supports tumor-subtype applicability as a hypothesis for our overall breast cancer PRS. We did not evaluate a subtype-specific PRS or establish subtype mediation, and their hazard ratios are not directly comparable with our case-status ORs. Partial attenuation in influence diagnostics, EHR ascertainment, prevalent-case selection and residual model misspecification remain possible contributors. The gene-specific box plots are descriptive and do not establish that carrier status causes a change in PRS.',
    ),
    (
        'text',
        'A prior All of Us analysis did not detect significant diagnosis-age differences using PRS tertiles and Kaplan–Meier/log-rank comparisons[20]. Our continuous PRS model identified a small association of −0.43 years per SD. Greater granularity or efficiency of continuous modeling is a possible explanation, but differences in cohort composition, ancestry, sample size, implementation and covariate adjustment also matter. A difference between significant and nonsignificant results is not itself evidence of a between-study difference. Both analyses use recorded diagnosis age, which may differ from biological onset.',
    ),
    (
        'text',
        'Several limitations constrain interpretation. Cases were older than controls, and EHR coverage, diagnostic practices and selection into the program can affect case and control classification. Controls excluded coded cancers, so results need not generalize to all people without breast cancer. Case-only diagnosis-age analysis can be affected by conditioning on disease[19]; this concern does not prove bias in the present estimate. Blood-derived variants, particularly in TP53 and NF1, may include acquired clonal events despite genotype filters. Predicted loss of function is transcript- and position-dependent, and final manual variant adjudication remains pending; annotation confidence alone is not clinical validation. Score selection followed earlier analyses, and the interaction audit was post hoc. External weights were not individually recalibrated, repeated CV is internal validation, and ancestry restriction limits generalizability. Family-history survey coding and temporal alignment of menopausal status remain separate phenotype-review issues (Table 1); these variables were descriptive and were not model covariates.',
    ),
    (
        'text',
        'The thirteen-gene burden therefore provides modest additional case-status discrimination and an association with younger recorded diagnosis age. Prospective validation, finalized variant adjudication and evaluation of absolute-risk calibration and clinical decision value are needed before individual use.',
    ),
]
main.append(('heading', 'Declarations'))
for i in range(69, 75):
    main.append(('text', ptext(ps[i])))
main.append(('heading', 'References'))
refs = [ptext(ps[i]) for i in range(76, 93)]
refs[4] += ' doi:10.1038/s41436-020-0862-x.'
refs[12] = (
    '13. Rowlands CF, et al. Population-based germline breast cancer gene association studies and meta-analysis to inform wider mainstream testing. Ann Oncol. 2024;35(10):892–901. doi:10.1016/j.annonc.2024.07.244.'
)
refs[15] = (
    '16. Richards S, Aziz N, Bale S, et al. Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. Genet Med. 2015;17(5):405–424. doi:10.1038/gim.2015.30.'
)
refs += [
    '18. Abou Tayoun AN, Pesaran T, DiStefano MT, et al. Recommendations for interpreting the loss of function PVS1 ACMG/AMP variant criterion. Hum Mutat. 2018;39(11):1517–1524. doi:10.1002/humu.23626.',
    '19. Dudbridge F, Allen RJ, Sheehan NA, et al. Adjustment for index event bias in genome-wide association studies of subsequent events. Nat Commun. 2019;10:1561. doi:10.1038/s41467-019-09381-w.',
    '20. Guare LA, Das J, Caruth L, Setia-Verma S. Social determinants of health and lifestyle risk factors modulate genetic susceptibility for women’s health outcomes. Pac Symp Biocomput. 2025;30:296–313. doi:10.1142/9789819807024_0022.',
]
main += [('text', r) for r in refs]

# All displayed model quantities come from aggregate CSVs.
cohort = pd.read_csv(inputs[8]).fillna('—')
cohort.columns = ['Characteristic', 'Cases', 'Controls', 'Diagnosis-age cases']
cohort['Characteristic'] = cohort.Characteristic.str.replace(
    'Age at last EHR observation', 'Age at observation endpoint', regex=False
)
cohort.loc[cohort.Characteristic.str.contains('Family history', case=False), 'Characteristic'] = (
    'Coded family history of breast cancer'
)
cohort.Characteristic = (
    cohort.Characteristic.str.replace('Negative', 'No affected relative coded', regex=False)
    .str.replace(
        'Positive (first- or second-degree relative)',
        'Positive (captured affected relatives)',
        regex=False,
    )
    .str.replace(
        'Former or unknown currency', 'Former or unknown current smoking status', regex=False
    )
)
cohort.to_csv(OUT / 'table1_cohort.csv', index=False)
coefs = pd.read_csv(inputs[2])
coefs = coefs[coefs.version.eq('C')]
assoc = pd.DataFrame(
    {
        'Predictor': ['PRS, per SD', 'Burden, per SD', 'PRS × burden'],
        'Beta': coefs.beta.map(lambda v: f'{v:.5f}'),
        'OR (95% CI)': [f'{r.OR:.3f} ({r.CI_low:.3f}–{r.CI_high:.3f})' for r in coefs.itertuples()],
        'P': coefs.p_value.map(lambda v: f'{v:.3g}'),
    }
)
aucs = pd.read_csv(inputs[3])
aucs = aucs[aucs.version.eq('C')]
auctable = pd.DataFrame(
    {
        'Model': aucs.model,
        'Mean CV AUC': aucs.AUC_mean.map(lambda v: f'{v:.4f}'),
        'Repeat-1 AUC (95% CI)': [
            f'{r.AUC_repeat1:.4f} ({r.repeat1_CI_low:.4f}–{r.repeat1_CI_high:.4f})'
            for r in aucs.itertuples()
        ],
    }
)
age = pd.read_csv(inputs[5])
age = age[age.version.eq('C')]
agetable = pd.DataFrame(
    {
        'Predictor': ['PRS, per SD', 'Burden, per SD'],
        'Years (95% CI)': [
            f'{r.beta:.2f} ({r.CI_beta_low:.2f} to {r.CI_beta_high:.2f})' for r in age.itertuples()
        ],
        'P': age.p_value.map(lambda v: f'{v:.3g}'),
    }
)
for name, frame in [
    ('table2_associations', assoc),
    ('table2_discrimination', auctable),
    ('table3_diagnosis_age', agetable),
]:
    frame.to_csv(OUT / (name + '.csv'), index=False)
main += [
    ('heading', 'Tables'),
    ('heading', 'Table 1. Descriptive cohort characteristics.'),
    ('table', cohort),
    (
        'text',
        'Percentages use each column’s full N, including missing categories: 8,528 cases, 134,727 controls and 7,996 diagnosis-age cases. The case-status fitted-model cohort excludes 2,051 controls with missing observation age and therefore has N=141,204. Observation endpoint is the later available clinical EHR or death date. Menopausal status reflects survey timing rather than necessarily diagnosis timing. Family history is the current coded survey-derived field: absence of an affected-relative code may include skipped or unrecognized responses and is not a validated negative history. Suppression and rounding are retained from the audited source table.',
    ),
    ('heading', 'Table 2. Thirteen-gene case-status associations and discrimination.'),
    ('table', assoc),
    ('table', auctable),
    (
        'text',
        'Association estimates use the joint model with HC1 robust standard errors, observation age, PCs 1–10, sequencing center and sequencer. Main effects are conditional on the other standardized predictor at its mean. Interaction OR is exp(βPRS×burden). Mean AUCs are across ten five-fold CV repeats; displayed CIs belong to repeat 1. Mean joint-versus-PRS ΔAUC=0.01328; repeat-1 paired ΔAUC=0.01328 (95% CI, 0.01133–0.01523).',
    ),
    ('heading', 'Table 3. Recorded diagnosis-age associations among 7,996 cases.'),
    ('table', agetable),
    (
        'text',
        'Predictors are standardized within the diagnosis-age case cohort. The model adjusts for PCs 1–10 and uses HC1 robust standard errors. R²=0.02628; adjusted R²=0.02481; PRS-plus-PC baseline R²=0.002617; adding burden increases R² by 0.02366. These are case-only associations with recorded diagnosis age.',
    ),
]

# Render verified current figures, avoiding historical comparison panels in main.
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(10, 5))
ax.axis('off')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
flow = [
    'Whole-genome sequencing source: 295,706',
    'Ancestry and eligibility selection: 150,900',
    'Descriptive cohort after QC: 143,255\n8,528 cases; 134,727 controls',
]
for i, t in enumerate(flow):
    y = 0.9 - i * 0.23
    ax.text(
        0.5,
        y,
        t,
        ha='center',
        va='center',
        bbox={'boxstyle': 'round,pad=.45', 'fc': '#edf3f9', 'ec': '#54708c'},
    )
    if i < 2:
        ax.annotate(
            '',
            (0.5, y - 0.17),
            (0.5, y - 0.06),
            arrowprops={'arrowstyle': '->', 'color': '#54708c'},
        )
for x, t in [
    (0.24, 'Case-status model: 141,204\n2,051 controls excluded\nfor missing observation age'),
    (0.76, 'Diagnosis-age model: 7,996 cases\n532 cases excluded\nby first-record rule'),
]:
    ax.text(
        x,
        0.12,
        t,
        ha='center',
        va='center',
        bbox={'boxstyle': 'round,pad=.45', 'fc': '#edf3f9', 'ec': '#54708c'},
    )
    ax.annotate('', (x, 0.24), (0.5, 0.36), arrowprops={'arrowstyle': '->', 'color': '#54708c'})
fig.tight_layout()
fig.savefig(OUT / 'figure1_flow.png', dpi=180)
fig.savefig(OUT / 'figure1_flow.svg')
plt.close(fig)


def rendered(name):
    p = OUT / (name + '.png')
    cairosvg.svg2png(url=str(M / 'figures' / name) + '.svg', write_to=str(p), output_width=1500)
    return p


roc = pd.read_csv(M / 'datas/case_control/four_models_ROC_repeat1_interpolated.csv').query(
    "version == 'C'"
)
fig, ax = plt.subplots(figsize=(7, 5))
for model, rows in roc.groupby('model', sort=False):
    ax.plot(rows.FPR, rows.TPR_interpolated, label=model)
ax.plot([0, 1], [0, 1], color='gray', linestyle='--')
ax.set(
    xlabel='False positive rate',
    ylabel='True positive rate',
    title='Current thirteen-gene models: repeat-1 held-out ROC',
)
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(OUT / 'figure2a_current_ROC.png', dpi=180)
fig.savefig(OUT / 'figure2a_current_ROC.svg')
plt.close(fig)
fig, ax = plt.subplots(figsize=(7, 4.5))
positions = list(range(len(aucs)))
ax.plot(positions, aucs.AUC_mean, 'o', label='Mean across ten repeats', color='#205a84')
ax.errorbar(
    [v + 0.10 for v in positions],
    aucs.AUC_repeat1,
    yerr=[aucs.AUC_repeat1 - aucs.repeat1_CI_low, aucs.repeat1_CI_high - aucs.AUC_repeat1],
    fmt='s',
    capsize=4,
    color='#cb7839',
    label='Repeat 1 and its 95% CI',
)
ax.set_xticks(positions, ['Covariates', 'PRS + covariates', 'Burden + covariates', 'Joint'])
ax.set_ylabel('AUC')
ax.set_title('Discrimination of current thirteen-gene models')
ax.legend(frameon=False)
for i, v in enumerate(aucs.AUC_mean):
    ax.annotate(f'{v:.4f}', (i, v), xytext=(0, 10), textcoords='offset points', ha='center')
fig.tight_layout()
fig.savefig(OUT / 'figure2b_current_AUC.png', dpi=180)
fig.savefig(OUT / 'figure2b_current_AUC.svg')
plt.close(fig)
joint = pd.read_csv(M / 'datas/case_control/PRS311_decile_C_carrier_joint_distribution.csv')
assert not joint.cell_suppressed.any()
fig, ax = plt.subplots(figsize=(7, 4.5))
for carrier, rows in joint.groupby('C_carrier'):
    ax.plot(
        rows.PRS_decile,
        rows.observed_case_proportion * 100,
        'o-',
        label='Carrier' if carrier else 'Noncarrier',
    )
ax.set(
    xlabel='PRS311 decile',
    ylabel='Observed case proportion (%)',
    title='Current thirteen-gene carrier status and PRS',
)
ax.set_xticks(range(1, 11))
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(OUT / 'figure3_current_joint_distribution.png', dpi=180)
fig.savefig(OUT / 'figure3_current_joint_distribution.svg')
plt.close(fig)
main += [
    ('heading', 'Figures'),
    ('image', OUT / 'figure1_flow.png'),
    (
        'text',
        'Figure 1. Cohort selection and analysis denominators. The case-status and diagnosis-age cohorts are separate branches of the descriptive cohort; the age analysis is restricted to cases.',
    ),
    ('image', OUT / 'figure2a_current_ROC.png'),
    ('image', OUT / 'figure2b_current_AUC.png'),
    (
        'text',
        'Figure 2. Discrimination of the four current thirteen-gene models. Receiver operating characteristic curves use repeat-1 held-out predictions; mean AUCs and repeat-1 intervals are reported separately in Table 2. Models include covariates alone, PRS, burden, or both scores and their product.',
    ),
    ('image', OUT / 'figure3_current_joint_distribution.png'),
    (
        'text',
        'Figure 3. Observed case proportions across PRS deciles, stratified by qualifying thirteen-gene carrier status. Numerator: cases in the decile–carrier cell; denominator: all participants in that cell. Carrier status is binary even though modeled burden is a weighted continuous index. These unadjusted proportions are not prospective absolute risks.',
    ),
    ('image', rendered('Figure_05a_age_by_13gene_burden_group')),
    ('image', rendered('Figure_05c_carrier_frequency_by_diagnosis_age')),
    (
        'text',
        'Figure 4. Recorded diagnosis age by thirteen-gene burden group and carrier frequency by diagnosis-age category. Carrier burden groups use the carrier-median threshold, whose ties can produce unequal group sizes. Carrier frequency is carriers divided by all cases in each eligible age category. These displays concern the 7,996-case diagnosis-age cohort and do not estimate disease incidence.',
    ),
]

supp += [
    ('title', 'Supplementary material: thirteen-gene burden revision'),
    (
        'text',
        'Companion to the revised 1006 manuscript. Current-model results and historical comparisons are explicitly separated. Reference numbers refer to the main manuscript.',
    ),
    ('heading', 'Supplementary Methods'),
    (
        'text',
        'Variant procedures follow the rules in the main Methods. PTV consequence and LOFTEE high-confidence status must occur on the same retained transcript for the same gene. Noncanonical retained transcripts can qualify. No blanket final-exon rule was adopted; clinical transcript/position adjudication remains pending. Missense selection uses germline clinical assertions, relevant disease and inheritance context, ≥2-star ClinVar review, and exclusions of low-penetrance, risk alleles, conflicting and opposing benign assertions. Release: 2026-09-13. Unknown or held-out eligibility records are excluded rather than silently assigned to the qualifying set.',
    ),
    (
        'text',
        'Normalized target alleles link annotation to original variant-dataset site and genotype quality. A site passes with defined-empty filters or singleton PASS. Individual genotype filters use PASS or “.”; genotype quality ≥20 is required. Allele fraction uses original target-ALT allele depth divided by total allelic depth, including all alleles. Heterozygotes require 0.30–0.80 and homozygous alternative calls ≥0.80. Missing allele-depth denominators or other required QC fields fail inclusion. Allele frequency must be known and consistent and strictly below 0.01 in the All of Us European population. Summed allele depth was reviewed separately and not substituted for an unavailable genotype depth field. No blanket depth ≥10 rule was adopted.',
    ),
    (
        'text',
        'PRS matching audit: 313 intended loci; 304 matched, including 281 direct and 23 frequency-resolved loci; eight unmatched and one unresolved. Both excluded CHEK2 tags were among the matched loci, yielding 302 retained scored loci in PRS311. Missing locus dosage was imputed with cohort mean dosage. Association models use full fitted-cohort sample-SD standardization; diagnosis-age models use within-case sample-SD standardization. CV uses training-fold means and sample SDs. Saved reference-standardized score files have a different reference scaling and should not be confused with fitted-model scaling.',
    ),
    (
        'text',
        'Interaction audit: influence and leverage deletion, multigene-carrier exclusion, binary carrier status, cubic-spline PRS and age, enrollment-age substitution/restriction and EHR coverage adjustment/restriction used the frozen current case-status cohort. Full-cohort genetic scaling was retained for subset fits. A joint gene-carrier model included gene main effects and PRS-by-gene terms; sparse genes were pooled. HC1 Wald tests evaluated the interaction terms jointly and equality across gene categories. Seven gene-specific tests were adjusted using Benjamini–Hochberg; other P values were nominal. Outcome-dependent influence deletion is diagnostic and does not prove removal of bias.',
    ),
    (
        'text',
        'Chronology and disclosure: historical analyses preceded selection of the thirteen-gene main score. Five score sensitivities were locked before its association evaluation; the interaction audit was later and exploratory. Final biological adjudication is incomplete. Tables and figures contain aggregate data only. Cells below 20 and secondary cells enabling reconstruction are suppressed; sparse gene categories are pooled, and participant identifiers or individual points are not displayed.',
    ),
]
weights = pd.read_csv(inputs[7])
weights_display = weights[['gene', 'OR', 'CI_low', 'CI_high', 'ln_OR']].copy()
weights_display.columns = ['Gene', 'External OR', 'CI low', 'CI high', 'ln(OR)']
supp += [
    ('heading', 'Supplementary Table S1. External thirteen-gene weights.'),
    ('table', weights_display),
    (
        'text',
        'Source: Rowlands et al.[13], Table 2 combined PTV plus non-PTV pathogenic-variant column. These are approximate external gene weights, because the source variant set differs from current eligible missense variants. STK11 OR=1.60 (95% CI, 0.48–5.30) is imprecise; CDH1 associations differ by tumor subtype. The score counts each qualifying gene once and does not encode variant or allele dosage.',
    ),
]
sens = pd.read_csv(inputs[9])
sens = sens[sens.version.eq('C') | sens.version.str.startswith('C_')].copy()
sens.version = sens.version.replace(
    {
        'C': 'Main thirteen-gene score',
        'C_without_STK11': 'Exclude STK11',
        'C_PTV_only': 'PTVs only',
        'C_TP53_NF1_strict_VAF': 'TP53/NF1 stricter allele fraction',
        'C_ClinVar_min1star': 'ClinVar ≥1 star',
        'C_without_CDH1': 'Exclude CDH1',
    }
)
supp += [
    ('heading', 'Supplementary Table S2. Current-panel score sensitivities.'),
    ('table', sens),
    (
        'text',
        'Case-status estimates are ORs; diagnosis-age estimates are years per SD. CIs are 95%. These score specifications were locked before current thirteen-gene association results, following earlier analyses. The current stricter allele-fraction sensitivity changes only the TP53/NF1 heterozygous upper bound to 0.70.',
    ),
]
supp += [
    ('heading', 'Supplementary Table S3. Current four-model AUCs by repeat.'),
    (
        'table',
        pd.read_csv(M / 'datas/case_control/four_models_AUC_by_repeat.csv').query("version == 'C'"),
    ),
]
historical_names = {
    'A': 'Legacy nine-gene PTV score',
    "A_prime": 'Same legacy carriers, updated weights',
    "A'": 'Same legacy carriers, updated weights',
    'A_reweighted': 'Same legacy carriers, updated weights',
    'A_prime_QC': 'Nine-gene PTV, current QC',
    'A_QC': 'Nine-gene PTV, current QC',
    'A_current_QC': 'Nine-gene PTV, current QC',
    'B': 'Nine-gene PTV plus eligible missense',
    'C': 'Thirteen-gene PTV plus eligible missense',
    'PRS311_baseline': 'PRS plus covariates',
}
hist = pd.read_csv(M / 'datas/comparisons/all_versions_CV_AUC_summary.csv')
hist = hist[hist.version.isin(historical_names)]
hist.version = hist.version.replace(historical_names)
supp += [
    ('heading', 'Supplementary Table S4. Historical score comparisons.'),
    ('table', hist),
    (
        'text',
        'Historical comparators are reported under their own definitions and are not main thirteen-gene estimates. CIs refer to repeat 1; means are across repeated CV. Additional sequential decomposition is in Table S5. No historical result is relabeled as the current score.',
    ),
]
delta = pd.read_csv(M / 'datas/comparisons/all_versions_CV_delta_AUC_summary.csv')
delta = delta[delta.version.isin(historical_names)].drop(columns=['CI_low_mean', 'CI_high_mean'])
delta.version = delta.version.replace(historical_names)
delta.comparison = delta.comparison.replace(
    {
        'joint_minus_PRS311_baseline': 'Joint minus PRS plus covariates',
        'C_minus_B': 'Thirteen-gene minus nine-gene PTV plus missense',
        'A_prime_minus_A': 'Updated weights minus legacy weights',
        'A_prime_QC_minus_A_prime': 'Current QC minus legacy QC',
        'B_minus_A_prime_QC': 'Eligible missense added to current-QC nine-gene PTVs',
    }
)
delta = delta.rename(columns={'version': 'Score', 'comparison': 'Comparison'})
delta_display = delta[
    [
        'Score',
        'Comparison',
        'delta_AUC_mean',
        'delta_AUC_repeat1',
        'repeat1_CI_low',
        'repeat1_CI_high',
        'repeat1_p',
    ]
].copy()
for c in ['delta_AUC_mean', 'delta_AUC_repeat1', 'repeat1_CI_low', 'repeat1_CI_high']:
    delta_display[c] = delta_display[c].map(lambda v: f'{v:.6f}')
delta_display['repeat1_p'] = delta_display.repeat1_p.map(lambda v: f'{v:.3g}')
supp += [
    ('heading', 'Supplementary Table S5. Paired AUC increments.'),
    ('table', delta_display),
    (
        'text',
        'Mean increments are repeated-CV summaries. Repeat-1 CIs and P values are paired DeLong results for that repeat, not intervals for mean increments. The sequential rows separate updated weighting, current QC, eligible missense inclusion and thirteen-gene expansion.',
    ),
    (
        'text',
        'For the historical overall contrast of the thirteen-gene score versus the legacy nine-gene PTV score, the repeat-1 AUC increment was 0.001905 (95% CI, 0.001221–0.002590; paired P=4.80×10⁻⁸). This historical contrast is different from the main joint-versus-PRS increment of 0.01328.',
    ),
]
supp += [
    ('heading', 'Supplementary Table S6. Exploratory carrier-specific diagnosis-age slopes.'),
    ('table', pd.read_csv(M / 'datas/age/carrier_PRS_slope_contrasts.csv')),
    (
        'text',
        'Binary carrier-by-PRS age interaction is exploratory. The interaction contrast is 0.98 years (95% CI, −0.08 to 2.04); uncertainty does not establish synergy. This is separate from the primary continuous-burden age model.',
    ),
]
supp += [
    ('heading', 'Supplementary Table S7. Gene-specific PRS interactions.'),
    ('table', pd.read_csv(M / 'tables/tableS7_gene_PRS_interactions.csv')),
    (
        'text',
        'Interaction OR is a ratio of PRS ORs in carriers versus noncarriers from a joint gene-carrier model, not a gene’s overall disease OR. Covariates: observation age, PCs 1–10, sequencing center and sequencer. Categories overlap for multigene carriers. Other pooled genes: RAD51C, RAD51D, TP53, PTEN, NF1, STK11 and CDH1. Counts are withheld for joint disclosure protection. Nominal P and Benjamini–Hochberg q values across seven categories are reported; heterogeneity P=0.0126 is nominal. Analyses are post hoc.',
    ),
]
supp += [
    ('heading', 'Supplementary Table S8. Interaction robustness analyses.'),
    ('table', pd.read_csv(M / 'tables/tableS8_interaction_sensitivities.csv')),
    (
        'text',
        'Fixed full-cohort scaling is used in subset fits. The binary-carrier coefficient is on a 0/1 scale and is not numerically comparable with continuous per-SD burden coefficients. Influence deletion is outcome-dependent. Other P values are nominal; EHR adjustment does not demonstrate a causal explanation.',
    ),
]
supp += [
    ('heading', 'Supplementary figures'),
    ('image', M / 'figures/supp_gene_PRS311_boxplots.png'),
    (
        'text',
        'Supplementary Figure S7A. Raw PRS311 distributions by qualifying gene-carrier group and case status. Boxes show quartiles and medians; whiskers extend to observations within 1.5 times the interquartile range. Individual points and outliers are omitted. Sparse categories are pooled and multigene carriers may appear in multiple groups. These descriptive distributions do not estimate causal effects of carrier status on PRS; case conditioning can induce risk-factor associations.',
    ),
    ('image', M / 'figures/supp_gene_PRS_interaction_forest.png'),
    (
        'text',
        'Supplementary Figure S7B. Gene-specific interaction ORs and 95% CIs. Values below 1 indicate a smaller PRS OR in carriers conditional on covariates, rather than a protective overall gene effect. The dashed line marks no multiplicative-odds interaction. See Table S7 for multiplicity adjustment and pooling.',
    ),
]
for label, name in [
    ('S1', 'Figure_02a_burden_raw_score_violin'),
    ('S2', 'Figure_02b_PRS311_raw_score_KDE'),
    ('S3', 'Figure_02c_13gene_burden_histogram_percentages'),
    ('S4', 'Figure_02d_PRS311_decile_PRS_burden_boxplots'),
]:
    supp += [
        ('image', rendered(name)),
        (
            'text',
            f'Supplementary Figure {label}. Current-score descriptive distribution ({name.replace("_"," ")}). Raw scores are burden indices or polygenic scores, not calibrated probabilities.',
        ),
    ]
supp += [
    ('image', rendered('Figure_06_primary_and_sensitivity_estimates')),
    (
        'text',
        'Supplementary Figure S5. Full association forest plots with current main thirteen-gene estimates, current-panel sensitivities and explicitly labeled historical nine-gene comparators. This historical-comparison figure is not a second definition of the main score.',
    ),
    ('image', rendered('Figure_03a_case_status_ROC')),
    (
        'text',
        'Supplementary Figure S8. Four-model ROC curves for the main thirteen-gene score and five score sensitivities. Each panel uses repeat-1 held-out predictions; panel labels identify the score definition.',
    ),
    ('image', rendered('Figure_03b_case_status_AUC')),
    (
        'text',
        'Supplementary Figure S9. Four-model discrimination across the main score and five sensitivities. Bars show mean AUC over ten CV repeats; intervals belong to repeat 1 and do not quantify uncertainty in the mean.',
    ),
]


def para(s, style=None):
    p = E.Element(q('p'))
    pr = E.SubElement(p, q('pPr'))
    if style:
        E.SubElement(pr, q('pStyle'), {q('val'): style})
    E.SubElement(pr, q('spacing'), {q('after'): '120'})
    r = E.SubElement(p, q('r'))
    t = E.SubElement(r, q('t'))
    t.text = str(s)
    t.set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')
    return p


def table(df):
    tb = E.Element(q('tbl'))
    pr = E.SubElement(tb, q('tblPr'))
    E.SubElement(pr, q('tblW'), {q('w'): '0', q('type'): 'auto'})
    borders = E.SubElement(pr, q('tblBorders'))
    for side in ['top', 'left', 'bottom', 'right', 'insideH', 'insideV']:
        E.SubElement(borders, q(side), {q('val'): 'single', q('sz'): '4', q('color'): 'B0B0B0'})
    for i, row in enumerate([list(df.columns)] + df.fillna('—').astype(str).values.tolist()):
        tr = E.SubElement(tb, q('tr'))
        if i == 0:
            E.SubElement(E.SubElement(tr, q('trPr')), q('tblHeader'))
        for s in row:
            c = E.SubElement(tr, q('tc'))
            p = para(s)
            r = p.find(q('r'))
            rp = E.Element(q('rPr'))
            r.insert(0, rp)
            E.SubElement(rp, q('sz'), {q('val'): '16'})
            if i == 0:
                E.SubElement(rp, q('b'))
            c.append(p)
    return tb


def write_doc(items, path):
    parts = copy.deepcopy(original)
    root = E.fromstring(parts['word/document.xml'])
    body = root.find(q('body'))
    sect = copy.deepcopy(body.find(q('sectPr')))
    for c in list(body):
        body.remove(c)
    rels = E.fromstring(parts['word/_rels/document.xml.rels'])
    cts = E.fromstring(parts['[Content_Types].xml'])
    if not any(x.get('Extension') == 'png' for x in cts):
        E.SubElement(cts, '{' + CT + '}Default', Extension='png', ContentType='image/png')
    for idx, (kind, data) in enumerate(items):
        if kind == 'table':
            body.append(table(data))
        elif kind == 'image':
            with Image.open(data) as im:
                w, h = im.size
            cx = 5486400
            cy = int(cx * h / w)
            rid = 'rIdRevisionImage' + str(idx)
            name = 'revision_' + str(idx) + '.png'
            E.SubElement(
                rels, '{' + REL + '}Relationship', Id=rid, Type=R + '/image', Target='media/' + name
            )
            parts['word/media/' + name] = data.read_bytes()
            frag = f'<w:p xmlns:w="{W}" xmlns:r="{R}" xmlns:wp="http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing" xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" xmlns:pic="http://schemas.openxmlformats.org/drawingml/2006/picture"><w:r><w:drawing><wp:inline><wp:extent cx="{cx}" cy="{cy}"/><wp:docPr id="{12000+idx}" name="Revision figure {idx}"/><a:graphic><a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/picture"><pic:pic><pic:nvPicPr><pic:cNvPr id="0" name="{name}"/><pic:cNvPicPr/></pic:nvPicPr><pic:blipFill><a:blip r:embed="{rid}"/><a:stretch><a:fillRect/></a:stretch></pic:blipFill><pic:spPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="{cx}" cy="{cy}"/></a:xfrm><a:prstGeom prst="rect"><a:avLst/></a:prstGeom></pic:spPr></pic:pic></a:graphicData></a:graphic></wp:inline></w:drawing></w:r></w:p>'
            body.append(E.fromstring(frag.encode()))
        else:
            body.append(
                para(
                    data, 'Title' if kind == 'title' else 'Heading1' if kind == 'heading' else None
                )
            )
    if sect is not None:
        body.append(sect)
    # Discard inherited images that are not used by the rebuilt document body.
    used = set(root.xpath('//@r:embed | //@r:link | //@r:id', namespaces={'r': R}))
    import posixpath

    removed = []
    for relationship in list(rels):
        if relationship.get('Type', '').endswith('/image') and relationship.get('Id') not in used:
            target = posixpath.normpath('word/' + relationship.get('Target'))
            parts.pop(target, None)
            removed.append(target)
            rels.remove(relationship)
    for entry in list(cts):
        if entry.get('PartName', '').lstrip('/') in removed:
            cts.remove(entry)
    parts['word/document.xml'] = xml(root)
    parts['word/_rels/document.xml.rels'] = xml(rels)
    parts['[Content_Types].xml'] = xml(cts)
    with zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED) as z:
        for n, v in parts.items():
            z.writestr(n, v)
    with zipfile.ZipFile(path) as z:
        assert z.testzip() is None
        for n in z.namelist():
            if n.endswith('.xml') or n.endswith('.rels'):
                E.fromstring(z.read(n))


def markdown(items):
    lines = []
    for kind, data in items:
        if kind == 'table':
            lines.append(data.to_markdown(index=False))
        elif kind == 'image':
            if data.parent != OUT:
                (OUT / data.name).write_bytes(data.read_bytes())
            lines.append(f'![{data.stem}]({data.name})')
        else:
            lines.append(
                ('# ' if kind == 'title' else '## ' if kind == 'heading' else '') + str(data)
            )
    return '\n\n'.join(lines) + '\n'


main = [
    (
        (
            kind,
            data.replace(
                'Participants with missing genotype rates >3% were excluded.',
                'Participants with missing genotype rates >3% in an approximately 1,500-locus chromosome-21 proxy were excluded; this was not a genome-wide missingness measurement.',
            ),
        )
        if kind in ['text', 'heading', 'title']
        else (kind, data)
    )
    for kind, data in main
]
main_md = markdown(main)
supp_md = markdown(supp)
assert 'MAP3K1' not in main_md
assert '0.00191' not in main_md
assert 'nine-gene' not in main_md and 'ten-gene' not in main_md
assert '302 scored loci' in main_md and '−0.43' in main_md
for name, items, md in [('main', main, main_md), ('supplement', supp, supp_md)]:
    write_doc(items, M / f'medrxiv_v1_{name}_13gene_revised_1008.docx')
    (OUT / f'{name}_13gene_revised.md').write_text(md)
(OUT / 'main_revision.diff').write_text(
    ''.join(
        difflib.unified_diff(
            source_text.splitlines(True),
            main_md.splitlines(True),
            fromfile='main_1006_source.txt',
            tofile='main_13gene_revised.md',
        )
    )
)
assert all(sha(ROOT / k) == v for k, v in hashes.items())
report = {
    'status': 'revised_draft_completed',
    'source_inputs_unchanged': True,
    'input_sha256': hashes,
    'main_docx': 'manuscript/medrxiv_v1_main_13gene_revised_1008.docx',
    'supplement_docx': 'manuscript/medrxiv_v1_supplement_13gene_revised_1008.docx',
    'checks': [
        'DOCX ZIP integrity',
        'all XML parts parse',
        'no obsolete main-panel gene or numeric increments',
        'aggregate model table values from CSV',
        'source hashes unchanged',
    ],
    'limits': [
        'No full Word/LibreOffice pagination render',
        'Author declarations remain placeholders from source',
        'Final variant adjudication remains pending',
        'Family-history coding remains unresolved',
    ],
    'revision_mechanism': 'Rebuilt clean DOCX body using the original package styles; source extraction and unified diff retained. Not Word tracked changes.',
}
(OUT / 'revision_validation.json').write_text(json.dumps(report, indent=2) + '\n')
(OUT / 'revision_log.md').write_text(
    '''# Revision log: 1006 manuscript to thirteen-gene main score

User instruction: use the research skill to revise manuscript 1006 against the thirteen-gene revision outline. Original DOCX and analysis data were preserved; no models were rerun.

| Outline commitment | Completed revision |
|---|---|
| Sole thirteen-gene burden in main | Replaced construction, Results, Abstract and Discussion; removed obsolete main-panel methods and values. |
| Accurate construction | Added all thirteen genes, PTV retained-transcript/LOFTEE rule, contextual eligible missense rules, actual genotype QC and once-per-gene weighted formula. |
| Correct PRS | Distinguished 313 intended, 304 matched and 302 retained loci in PRS311; retained dosage-imputation description. |
| Updated results | Main ORs 1.554/1.269; interaction OR 0.972; mean AUCs 0.7452/0.7682/0.7588/0.7815; mean joint-minus-PRS increment 0.01328. |
| Diagnosis age | PRS −0.43 years/SD, burden −1.82 years/SD, R² 0.02628; conclusions restricted to recorded diagnosis age among cases. |
| Historical comparisons in supplement | All five historical stages and four sequential increments included with descriptive names; historical increment 0.001905 distinguished from the main contrast. |
| Figure consistency | Generated a branching flow, current-only ROC/AUC and observed joint proportions; moved full sensitivity/legacy forest figures to supplement. |
| Negative interaction | Retained exploratory gene interaction table, box plots, robust audit and Barnes 2020 citation; no biological antagonism or positive synergy claim. |
| Literature comparison | Added verified Guare 2025 comparison, with continuous-model efficiency stated only as a possible explanation. |
| Status and limits | Reported chronology, pending adjudication, external-weight mismatch, EHR ascertainment, case conditioning and ancestry scope. |

Independent evidence and draft reviews identified missing historical tables, multipanel figure/caption mismatches and an overstated missingness measurement; these were corrected. The chromosome-21 approximately 1,500-locus proxy is now stated explicitly.

The clean DOCX body was rebuilt using the original document package and styles. `main_revision.diff` is a unified text diff against extracted source text; it is not Word tracked changes or an ARS block-patch artifact. `main_13gene_revised.md` and `supplement_13gene_revised.md` provide readable review copies. Numerical source hashes and package/XML checks are recorded in `revision_validation.json`. The deterministic acronym checker was run; gene symbols, score names, software/filter labels and mathematical subscripts create advisory false positives/coverage limits, recorded separately.

Remaining author/publication work: original correspondence, acknowledgements, ethics, funding, competing interests and contribution placeholders need completion; variant adjudication and family-history phenotype coding remain unresolved. No prospective calibration or clinical-use validation is claimed. Full Word pagination was not rendered because a Word/LibreOffice renderer was unavailable.
'''
)
print(json.dumps(report, indent=2))

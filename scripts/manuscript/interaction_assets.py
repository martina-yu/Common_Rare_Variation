"""Create aggregate figures, supplement tables, and scoped DOCX revision copies."""

from pathlib import Path
import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.functions.paths import PROJECT_ROOT

import hashlib, json, zipfile, re
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import NullLocator, NullFormatter
from lxml import etree as E

ROOT = PROJECT_ROOT
M = ROOT / 'manuscript'
DATA = M / 'datas/interaction'
FIG = M / 'figures'
TAB = M / 'tables'
for p in [DATA, FIG, TAB]:
    p.mkdir(parents=True, exist_ok=True)
source_main = M / 'medrxiv_v1_main_revised_1006.docx'
source_supp = M / 'medrxiv_v1_supplement_1002.docx'
inputs = [
    source_main,
    source_supp,
    ROOT / 'results/table/59_interaction_audit_gene_interactions.csv',
    ROOT / 'results/table/59_interaction_audit_models.csv',
    ROOT / 'data/29_rowlands_extension/section9/frozen_model_inputs.csv.gz',
    ROOT / 'data/29_rowlands_extension/section7/unweighted_gene_carriers.csv.gz',
]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


hashes = {str(p.relative_to(ROOT)): sha(p) for p in inputs}
t = pd.read_csv(inputs[2])
models = pd.read_csv(inputs[3])
t.to_csv(DATA / 'gene_specific_interaction_OR.csv', index=False)
table = pd.DataFrame(
    {
        'Gene category': t.gene.str.replace('Other_pooled_panel_genes', 'Other pooled panel genes'),
        'Interaction beta': t.beta.map(lambda x: f'{x:.3f}'),
        'Interaction OR (95% CI)': [
            f'{r.interaction_OR:.3f} ({r.CI_low:.3f}–{r.CI_high:.3f})' for r in t.itertuples()
        ],
        'Nominal P': t.p_value.map(lambda x: f'{x:.3g}'),
        'BH-adjusted q': t.BH_q_value.map(lambda x: f'{x:.3g}'),
    }
)
table.to_csv(TAB / 'tableS7_gene_PRS_interactions.csv', index=False)
models.to_csv(DATA / 'interaction_sensitivity_models.csv', index=False)
selected = [
    'C_baseline_reproduction',
    'exclude_top100_interaction_influence',
    'exclude_top100_design_leverage',
    'exclude_multigene_carriers',
    'binary_C_carrier',
    'nonlinear_last_observation_age',
    'enrollment_age_replaces_last_observation_age',
    'enrollment_age_50_to74_baseline_covariates',
    'enrollment_age_50_to74_nonlinear_enrollment_age',
    'EHR_log_span_and_date_count_adjustment',
    'EHR_span_ge365days_dates_ge5_age_spline_EHR_adjusted',
]
sens = models.loc[
    models.model.isin(selected), ['model', 'beta', 'OR', 'CI_low', 'CI_high', 'p_value']
].copy()
sens.to_csv(TAB / 'tableS8_interaction_sensitivities.csv', index=False)

d = pd.read_csv(inputs[4]).set_index('person_id')
g = pd.read_csv(inputs[5]).set_index('sample_id').loc[d.index]
assert len(d) == 141204 and d.is_case.sum() == 8528
genes = [
    c[:-5]
    for c in g.columns
    if c.endswith('_main') and not c.endswith('_PTV_main') and not c.endswith('_missense_main')
]
dense = [x for x in t.gene if x != 'Other_pooled_panel_genes']
sparse = [x for x in genes if x not in dense]
groups = {'Noncarrier': g[[x + '_main' for x in genes]].sum(axis=1).eq(0)}
groups.update({x: g[x + '_main'].eq(1) for x in dense})
groups['Other pooled'] = g[[x + '_main' for x in sparse]].max(axis=1).eq(1)
records = []
fig, ax = plt.subplots(figsize=(11, 4.8))
colors = {0: 'tab:orange', 1: 'tab:blue'}
for i, (name, mask) in enumerate(groups.items()):
    for case in [0, 1]:
        x = d.loc[mask & d.is_case.eq(case), 'PRS311'].dropna().to_numpy()
        assert len(x) >= 20, 'Sparse subgroup cannot be plotted individually'
        q1, med, q3 = np.quantile(x, [0.25, 0.5, 0.75])
        iqr = q3 - q1
        within = x[(x >= q1 - 1.5 * iqr) & (x <= q3 + 1.5 * iqr)]
        lo, hi = within.min(), within.max()
        st = {'q1': q1, 'med': med, 'q3': q3, 'whislo': lo, 'whishi': hi, 'fliers': []}
        art = ax.bxp(
            [st],
            positions=[i + (-0.18 if case == 0 else 0.18)],
            widths=0.29,
            showfliers=False,
            patch_artist=True,
            boxprops={'facecolor': colors[case], 'alpha': 0.65},
            medianprops={'color': 'black'},
        )
        records.append(
            dict(
                group=name,
                status='Case' if case else 'Control',
                q25=q1,
                median=med,
                q75=q3,
                whisker_low=lo,
                whisker_high=hi,
            )
        )
ax.set_xticks(range(len(groups)), list(groups))
ax.set_ylabel('PRS311 (raw score)')
ax.set_title('PRS311 distributions by qualifying gene-carrier group')
ax.legend(
    handles=[
        Patch(facecolor=colors[0], alpha=0.65, label='Controls'),
        Patch(facecolor=colors[1], alpha=0.65, label='Cases'),
    ],
    frameon=False,
)
ax.grid(False)
fig.tight_layout()
boxbase = FIG / 'supp_gene_PRS311_boxplots'
fig.savefig(str(boxbase) + '.svg')
fig.savefig(str(boxbase) + '.png', dpi=180)
plt.close(fig)
pd.DataFrame(records).to_csv(DATA / 'gene_PRS311_boxplot_aggregate_stats.csv', index=False)

fig, ax = plt.subplots(figsize=(7.2, 4.5))
pos = np.arange(len(t))
ax.errorbar(
    t.interaction_OR,
    pos,
    xerr=np.vstack([t.interaction_OR - t.CI_low, t.CI_high - t.interaction_OR]),
    fmt='o',
    color='black',
    capsize=3,
)
ax.axvline(1, color='gray', linestyle='--', linewidth=1)
ax.set_yticks(pos, t.gene.str.replace('Other_pooled_panel_genes', 'Other pooled panel genes'))
ax.invert_yaxis()
ax.set_xscale('log')
ax.set_xticks([0.5, 0.75, 1, 1.5], ["0.5", "0.75", "1.0", "1.5"])
ax.xaxis.set_minor_locator(NullLocator())
ax.xaxis.set_minor_formatter(NullFormatter())
ax.set_xlabel('Gene × PRS interaction OR (95% CI)')
ax.set_title('Gene-specific PRS interaction with case status')
ax.grid(False)
fig.tight_layout()
forestbase = FIG / 'supp_gene_PRS_interaction_forest'
fig.savefig(str(forestbase) + '.svg')
fig.savefig(str(forestbase) + '.png', dpi=180)
plt.close(fig)

methods = 'Post hoc interaction audit. We assessed the current 13-gene PTV plus eligible pathogenic/likely pathogenic missense score (version C) using the frozen case-status cohort. Sensitivities included deletion of the 100 observations with greatest absolute one-step interaction influence or design leverage, exclusion of multigene carriers, binary carrier status, cubic-spline PRS and age terms, enrollment-age restriction, and EHR coverage adjustment or restriction. A joint gene-carrier model included gene main effects and PRS-by-gene terms; sparse genes were pooled. HC1 robust Wald tests assessed joint interaction and heterogeneity, and Benjamini–Hochberg correction was applied across seven gene-specific interaction tests. These analyses were exploratory; other sensitivity P values were nominal. Full-cohort PRS and burden scaling was retained in subset fits.'
baseline = 'In the current 13-gene version-C joint model, the OR per standard deviation was 1.554 (95% CI, 1.518–1.591) for PRS311 and 1.269 (1.251–1.286) for rare-variant burden, conditional on the other score at its mean. Their interaction was negative (beta −0.0289; interaction OR, 0.972; 95% CI, 0.958–0.985; P=6.0×10−5). This represents departure from multiplicative odds, rather than an additive absolute-risk interaction.'
results = 'The negative interaction persisted in exploratory analyses using binary carrier status (interaction OR, 0.823; 95% CI, 0.744–0.910; P=1.50×10−4), nonlinear PRS terms (joint interaction P=0.00249), alternative age specifications, and EHR coverage adjustments or restrictions (Supplementary Table S8). Removing the 100 most interaction-influential observations reduced the coefficient magnitude by approximately 27%, while retaining a negative estimate (beta −0.0211; P=0.00229). Interaction coefficients differed between modeled gene categories (heterogeneity P=0.0126). BRCA1 showed the clearest negative interaction (OR, 0.652; 95% CI, 0.514–0.828; BH-adjusted q=0.00301); the BRCA2 nominal P value was 0.0217, but its BH-adjusted q value was 0.0758 (Supplementary Table S7; Supplementary Figure S7B). Excluding BRCA1 carriers attenuated the pooled interaction (beta −0.0206; nominal P=0.0267). Descriptive PRS distributions by gene-carrier group and case status are shown in Supplementary Figure S7A.'
discussion = 'The negative interaction indicates a weaker joint association than expected under a multiplicative-odds model, not a protective effect or biological antagonism. The gene-specific analyses suggest that pooling heterogeneous PRS associations, with a prominent BRCA1 contribution, may partly account for the panel-level result. One possible explanation is the applicability of an overall breast-cancer PRS to the tumor-subtype distribution in BRCA1 carriers. Barnes et al. found that an ER-negative-specific PRS was most strongly associated with breast-cancer risk in BRCA1 carriers, whereas an overall breast-cancer PRS performed best in BRCA2 carriers; PRS associations in carriers were smaller than general-population estimates[5]. These observations support a subtype-related explanation as a hypothesis, but their carrier-cohort hazard ratios are not directly comparable with our EHR case-control odds ratios. We did not evaluate an ER-negative-specific PRS or establish subtype mediation. Influence diagnostics also showed partial attenuation, and ascertainment, prevalent-case selection, and residual model misspecification remain possible contributors. The descriptive box plots do not demonstrate a causal effect of gene-carrier status on PRS; case stratification itself can induce associations between risk factors.'
ref = 'Barnes DR, Rookus MA, McGuffog L, et al. Polygenic risk scores and breast and epithelial ovarian cancer risks for carriers of BRCA1 and BRCA2 pathogenic variants. Genet Med. 2020;22(10):1653–1666. doi:10.1038/s41436-020-0862-x.'
foot = (
    'Interaction OR is exp(beta) for the gene-carrier × standardized PRS311 term: a ratio of PRS odds ratios in carriers versus noncarriers, not the gene’s overall breast-cancer OR. The joint model includes gene-carrier main effects, their PRS interactions, age at last observation, PCs 1–10, sequencing center, and sequencer. Gene categories can overlap for multigene carriers. Other pooled genes: '
    + ', '.join(sparse)
    + '. Each plotted case/control group meets the disclosure threshold. Counts are omitted to preserve joint disclosure protection. Nominal P values and BH q values across seven gene terms are reported; the omnibus heterogeneity P value is nominal. The analysis is post hoc and exploratory.'
)
boxcaption = 'Supplementary Figure S7A. PRS311 distributions by qualifying gene-carrier group and case status in the frozen case-status cohort (N=141,204). Boxes show Q1–Q3 with median lines; whiskers extend to observations within 1.5 times the IQR. Outliers and individual points are not displayed. Noncarriers have no qualifying panel variant. Carriers include current qualifying PTVs or P/LP missense variants; multigene carriers can appear in more than one group. Sparse gene categories are pooled. These unadjusted descriptive distributions neither estimate interaction nor establish an effect of carrier status on PRS; conditioning on case status can induce associations between genetic risk factors.'
forestcaption = 'Supplementary Figure S7B. Gene-specific PRS–case-status interaction ORs and 95% confidence intervals from the joint gene-carrier model. OR below 1 denotes a smaller PRS odds ratio in carriers relative to noncarriers, conditional on modeled covariates; it does not denote a protective gene effect. The dashed line marks no multiplicative-odds interaction. See Supplementary Table S7 for nominal P values, BH q values, pooling and overlap definitions. Exploratory estimates.'
caption_path = M / 'figure_captions/gene_PRS_interaction_supplement.md'
caption_path.write_text(boxcaption + '\n\n' + forestcaption + '\n')
md = (
    '# Scoped additions: negative PRS–burden interaction\n\n## Methods addition\n\n'
    + methods
    + '\n\n## Results replacement and addition\n\n'
    + baseline
    + '\n\n'
    + results
    + '\n\n## Discussion addition\n\n'
    + discussion
    + '\n\n## Reference 5 (existing reference, DOI added)\n\n'
    + ref
    + '\n\n## Supplementary Table S7\n\n'
    + table.to_markdown(index=False)
    + '\n\n'
    + foot
    + '\n\n'
    + boxcaption
    + '\n\n'
    + forestcaption
    + '\n\nScope: these are targeted current-version-C additions. Other legacy sections/tables in the source DOCX files have not undergone a full numeric or methodological refresh. Original DOCX files are preserved.\n'
)
(M / 'interaction_results_discussion_additions.md').write_text(md)

# Preserve all existing DOCX package parts, changing only designated XML parts.
W = 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'
R = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'
REL = 'http://schemas.openxmlformats.org/package/2006/relationships'
CT = 'http://schemas.openxmlformats.org/package/2006/content-types'


def q(n):
    return '{' + W + '}' + n


def para(text, bold=False):
    p = E.Element(q('p'))
    run = E.SubElement(p, q('r'))
    if bold:
        E.SubElement(E.SubElement(run, q('rPr')), q('b'))
    tx = E.SubElement(run, q('t'))
    tx.text = text
    tx.set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')
    return p


def text(p):
    return ''.join(p.itertext())


def xml(b):
    return E.tostring(b, encoding='UTF-8', xml_declaration=True, standalone=True)


def load_doc(p):
    with zipfile.ZipFile(p) as z:
        return {n: z.read(n) for n in z.namelist()}


def save_doc(parts, p):
    with zipfile.ZipFile(p, 'w', zipfile.ZIP_DEFLATED) as z:
        for n, v in parts.items():
            z.writestr(n, v)
    with zipfile.ZipFile(p) as z:
        assert z.testzip() is None
        for n in z.namelist():
            if n.endswith('.xml') or n.endswith('.rels'):
                E.fromstring(z.read(n))


parts = load_doc(source_main)
root = E.fromstring(parts['word/document.xml'])
body = root.find(q('body'))
ps = body.findall(q('p'))
hit = [p for p in ps if text(p).startswith('In controls, PRS311 was not correlated')]
assert len(hit) == 1
idx = body.index(hit[0])
body.remove(hit[0])
body.insert(idx, para(baseline))
body.insert(idx + 1, para(results))
hit = [p for p in ps if text(p).startswith('Sensitivity and post hoc analyses.')]
assert len(hit) == 1
body.insert(body.index(hit[0]) + 1, para(methods))
hit = [p for p in ps if text(p).startswith('Our per-standard-deviation PRS odds ratio')]
assert len(hit) == 1
body.insert(body.index(hit[0]) + 1, para(discussion))
hit = [p for p in ps if text(p).startswith('5. Barnes DR')]
assert len(hit) == 1
idx = body.index(hit[0])
body.remove(hit[0])
body.insert(idx, para('5. ' + ref))
parts['word/document.xml'] = xml(root)
main_out = M / 'medrxiv_v1_main_interaction_revision_1008.docx'
save_doc(parts, main_out)

parts = load_doc(source_supp)
root = E.fromstring(parts['word/document.xml'])
body = root.find(q('body'))
sect = body.find(q('sectPr'))


def append(node):
    body.insert(body.index(sect) if sect is not None else len(body), node)


def page():
    p = E.Element(q('p'))
    E.SubElement(E.SubElement(p, q('r')), q('br')).set(q('type'), 'page')
    append(p)


def addtable(frame):
    tb = E.Element(q('tbl'))
    pr = E.SubElement(tb, q('tblPr'))
    E.SubElement(pr, q('tblW'), {q('w'): '0', q('type'): 'auto'})
    borders = E.SubElement(pr, q('tblBorders'))
    for side in ['top', 'left', 'bottom', 'right', 'insideH', 'insideV']:
        E.SubElement(borders, q(side), {q('val'): 'single', q('sz'): '4', q('color'): 'B0B0B0'})
    for i, row in enumerate([list(frame.columns)] + frame.astype(str).values.tolist()):
        tr = E.SubElement(tb, q('tr'))
        for value in row:
            cell = E.SubElement(tr, q('tc'))
            cell.append(para(value, bold=i == 0))
    append(tb)


rels = E.fromstring(parts['word/_rels/document.xml.rels'])
cts = E.fromstring(parts['[Content_Types].xml'])
if not any(x.get('Extension') == 'png' for x in cts):
    E.SubElement(cts, '{' + CT + '}Default', Extension='png', ContentType='image/png')


def picture(path, num):
    from PIL import Image

    with Image.open(path) as im:
        w, h = im.size
    cx = 5486400
    cy = int(cx * h / w)
    rid = 'rIdInteractionImage' + str(num)
    name = 'interaction_' + str(num) + '.png'
    assert not any(x.get('Id') == rid for x in rels)
    E.SubElement(
        rels, '{' + REL + '}Relationship', Id=rid, Type=R + '/image', Target='media/' + name
    )
    parts['word/media/' + name] = path.read_bytes()
    fragment = f'''<w:p xmlns:w="{W}" xmlns:r="{R}" xmlns:wp="http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing" xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" xmlns:pic="http://schemas.openxmlformats.org/drawingml/2006/picture"><w:r><w:drawing><wp:inline><wp:extent cx="{cx}" cy="{cy}"/><wp:docPr id="{9000+num}" name="Interaction figure {num}"/><a:graphic><a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/picture"><pic:pic><pic:nvPicPr><pic:cNvPr id="0" name="{name}"/><pic:cNvPicPr/></pic:nvPicPr><pic:blipFill><a:blip r:embed="{rid}"/><a:stretch><a:fillRect/></a:stretch></pic:blipFill><pic:spPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="{cx}" cy="{cy}"/></a:xfrm><a:prstGeom prst="rect"><a:avLst/></a:prstGeom></pic:spPr></pic:pic></a:graphicData></a:graphic></wp:inline></w:drawing></w:r></w:p>'''
    append(E.fromstring(fragment.encode()))


page()
append(
    para(
        'Supplementary Table S7. Gene-specific PRS interactions in the current 13-gene analysis.',
        True,
    )
)
addtable(table)
append(para(foot))
page()
append(para('Supplementary Table S8. Exploratory interaction sensitivity analyses.', True))
display = pd.DataFrame(
    {
        'Sensitivity': sens.model.str.replace('_', ' '),
        'Beta': sens.beta.map(lambda x: f'{x:.4f}'),
        'Interaction OR (95% CI)': [
            f'{r.OR:.3f} ({r.CI_low:.3f}–{r.CI_high:.3f})' for r in sens.itertuples()
        ],
        'Nominal P': sens.p_value.map(lambda x: f'{x:.3g}'),
    }
)
addtable(display)
append(
    para(
        'These are current version-C sensitivity analyses with fixed full-cohort scaling. The binary-carrier coefficient uses a 0/1 indicator and is not numerically comparable with continuous per-SD burden coefficients. Subset counts are withheld for joint disclosure protection. Top-influence deletion is an outcome-dependent diagnostic. Other P values are nominal; no causal correction is implied by EHR adjustment.'
    )
)
page()
append(para('Supplementary Figure S7A', True))
picture(Path(str(boxbase) + '.png'), 1)
append(para(boxcaption))
page()
append(para('Supplementary Figure S7B', True))
picture(Path(str(forestbase) + '.png'), 2)
append(para(forestcaption))
parts['word/document.xml'] = xml(root)
parts['word/_rels/document.xml.rels'] = xml(rels)
parts['[Content_Types].xml'] = xml(cts)
supp_out = M / 'medrxiv_v1_supplement_interaction_revision_1008.docx'
save_doc(parts, supp_out)
assert all(sha(ROOT / k) == v for k, v in hashes.items())
meta = {
    'status': 'scoped_manuscript_assets_completed',
    'main_copy': str(main_out.relative_to(ROOT)),
    'supplement_copy': str(supp_out.relative_to(ROOT)),
    'source_inputs_unchanged': True,
    'input_sha256': hashes,
    'boxplot_group_count': len(groups),
    'boxplot_case_control_cells': len(records),
    'minimum_cell_threshold': 20,
    'no_participant_rows_exported': True,
    'reference': {
        'existing_number': 5,
        'doi': '10.1038/s41436-020-0862-x',
        'verified_by': 'primary article search metadata and abstract',
    },
    'limitations': [
        'Boxplots descriptive, not causal or interaction estimates',
        'Other source-manuscript legacy results/methods require broader refresh',
        'DOCX package/XML/image relationships checked; full Word pagination not rendered',
    ],
    'SVG_sizes': {
        x.name: x.stat().st_size
        for x in [Path(str(boxbase) + '.svg'), Path(str(forestbase) + '.svg')]
    },
}
(ROOT / 'results/report/60_interaction_manuscript_assets_validation.json').write_text(
    json.dumps(meta, indent=2) + '\n'
)
print(
    json.dumps(
        {
            'status': meta['status'],
            'main_copy': meta['main_copy'],
            'supplement_copy': meta['supplement_copy'],
            'source_inputs_unchanged': True,
            'boxplot_groups': len(groups),
            'sparse_genes_pooled': sparse,
        },
        indent=2,
    )
)

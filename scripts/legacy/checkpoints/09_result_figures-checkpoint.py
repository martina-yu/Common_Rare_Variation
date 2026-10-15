"""
- Attribute additional carriers.
- Summarize model results.
- Plot age associations and AUC comparisons.
- Plot manuscript score descriptions, discrimination and sensitivity estimates.
- Organize manuscript figures and insert them into Word on explicit task selection.
"""
from __future__ import annotations
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
# ============================================================================
# 1. attribution: attribute_new_carriers
# ============================================================================
from pathlib import Path
import json
import numpy as np
import pandas as pd
import yaml

CARRIER_ATTRIBUTION_ROOT=Path(__file__).resolve().parents[1];CARRIER_ATTRIBUTION_BASE=CARRIER_ATTRIBUTION_ROOT/'data/29_rowlands_extension';CARRIER_ATTRIBUTION_OUT=CARRIER_ATTRIBUTION_BASE/'section11_models'

def attribute_new_carriers():
    cfg=yaml.safe_load((CARRIER_ATTRIBUTION_ROOT/'config/rowlands_13gene_extension.yaml').read_text())
    scores=pd.read_csv(CARRIER_ATTRIBUTION_BASE/'section8/comparable_scores_case_control.csv.gz').set_index('sample_id')
    current=pd.read_csv(CARRIER_ATTRIBUTION_BASE/'section7/unweighted_gene_carriers.csv.gz').set_index('sample_id').loc[scores.index]
    pheno=pd.read_csv(CARRIER_ATTRIBUTION_BASE/'section9/frozen_model_inputs.csv.gz').set_index('person_id').loc[scores.index]
    gene_weights=pd.read_csv(CARRIER_ATTRIBUTION_ROOT/cfg['burden']['external_weight_table']).set_index('gene').ln_OR
    a=scores.carried_gene_count_A.gt(0);q=scores.carried_gene_count_A_prime_QC.gt(0);b=scores.carried_gene_count_B.gt(0);c=scores.carried_gene_count_C.gt(0)
    flags=pd.DataFrame({'A_carrier':a,'current_nine_PTV_carrier':q,'B_carrier':b,'C_carrier':c})
    flags['C_only_vs_A']=c&~a;flags['A_only_vs_C']=a&~c
    flags['QC_PTV_added_vs_A']=q&~a
    flags['missense_added_B_vs_QC_PTV']=b&~q
    flags['added_four_genes_C_vs_B']=c&~b
    # Hierarchical disjoint attribution among actual newly added C carriers.
    flags['new_C_current_nine_PTV']=c&~a&q
    flags['new_C_nine_missense']=c&~a&~q&b
    flags['new_C_added_four_only']=c&~a&~b
    assert np.array_equal(flags.C_only_vs_A.astype(int),flags[['new_C_current_nine_PTV','new_C_nine_missense','new_C_added_four_only']].sum(axis=1))
    flags.index.name='sample_id';flags.to_csv(CARRIER_ATTRIBUTION_OUT/'carrier_transition_flags_internal.csv.gz')
    rows=[];genes=[]
    age=set(pd.read_csv(CARRIER_ATTRIBUTION_BASE/'diagnosis_age_ids.csv.gz').person_id)
    for label,members in [('case_control',set(scores.index)),('diagnosis_age',age)]:
        ix=sorted(members);f=flags.loc[ix];d=current.loc[ix];p=pheno.loc[ix];new=f.C_only_vs_A
        for name in flags:
            x=f[name];rows.append(dict(cohort=label,metric=name,participants=int(x.sum()),cases=int((x&p.is_case.eq(1)).sum()),controls=int((x&p.is_case.eq(0)).sum())))
        for g in cfg['panels']['main_13']:
            x=d[f'{g}_main'].eq(1);px=d[f'{g}_PTV_main'].eq(1);mx=d[f'{g}_missense_main'].eq(1)
            genes.append(dict(cohort=label,gene=g,carrier_n=int(x.sum()),PTV_n=int(px.sum()),missense_n=int(mx.sum()),missense_only_n=int((mx&~px).sum()),
                newly_added_vs_A_n=int((x&new).sum()),new_cases=int((x&new&p.is_case.eq(1)).sum()),
                raw_burden_sum=float(x.sum()*gene_weights[g]),new_raw_burden_sum=float((x&new).sum()*gene_weights[g])))
    pd.DataFrame(rows).to_csv(CARRIER_ATTRIBUTION_OUT/'carrier_increment_attribution_internal.csv',index=False)
    gd=pd.DataFrame(genes);gd['share_raw_burden']=gd.raw_burden_sum/gd.groupby('cohort').raw_burden_sum.transform('sum')
    gd['share_new_raw_burden']=gd.new_raw_burden_sum/gd.groupby('cohort').new_raw_burden_sum.transform('sum')
    gd.to_csv(CARRIER_ATTRIBUTION_OUT/'gene_contribution_dominance_internal.csv',index=False)
    # Kept normalized allele calls; each person/gene weight split across their
    # retained alleles for an additive descriptive contribution (not causal OR).
    columns=['sample_id','gene','normalized_key','PTV/missense','keep_flag','low_penetrance_flag','conflict_flag','VAF','GQ','DP']
    calls=pd.read_csv(CARRIER_ATTRIBUTION_BASE/'section7/sample_variant_audit.csv.gz',usecols=columns,keep_default_na=False)
    calls=calls.loc[calls.keep_flag.eq(1)].drop_duplicates(['sample_id','gene','normalized_key'])
    calls['is_case']=calls.sample_id.map(pheno.is_case).astype(int)
    calls['new_C_vs_A']=calls.sample_id.map(flags.C_only_vs_A).astype(bool)
    calls['weighted_credit']=calls.gene.map(gene_weights)/calls.groupby(['sample_id','gene']).normalized_key.transform('size')
    pos=pd.read_csv(CARRIER_ATTRIBUTION_BASE/'section9/ptv_hc_transcript_position_review.tsv.gz',sep='\t',keep_default_na=False)
    unc=pd.read_csv(CARRIER_ATTRIBUTION_BASE/'section9/LOFTEE_multitranscript_review_internal.tsv.gz',sep='\t',keep_default_na=False)
    def collapse(frame,column):
        return frame.groupby(['gene_symbol','normalized_key'])[column].apply(lambda x:';'.join(sorted(set(v for v in x if v)))).to_dict()
    pf=collapse(pos,'position_review_flags');uf=collapse(unc,'uncertainty_flags')
    calls['position_review_flags']=[pf.get((r.gene,r.normalized_key),'') for r in calls.itertuples()]
    calls['transcript_uncertainty_flags']=[uf.get((r.gene,r.normalized_key),'') for r in calls.itertuples()]
    variant=[]
    for (g,k),d in calls.groupby(['gene','normalized_key']):
        n=d.loc[d.new_C_vs_A]
        variant.append(dict(gene=g,normalized_key=k,variant_class=';'.join(sorted(set(d['PTV/missense']))),carrier_n=d.sample_id.nunique(),
            case_n=d.loc[d.is_case.eq(1),'sample_id'].nunique(),new_C_vs_A_n=n.sample_id.nunique(),new_C_case_n=n.loc[n.is_case.eq(1),'sample_id'].nunique(),
            allocated_raw_burden_credit=d.weighted_credit.sum(),allocated_new_raw_credit=n.weighted_credit.sum(),
            position_review_flags=';'.join(sorted(set(v for v in d.position_review_flags if v))),
            transcript_uncertainty_flags=';'.join(sorted(set(v for v in d.transcript_uncertainty_flags if v))),
            any_low_penetrance_flag=any(str(v).lower()=='true' for v in d.low_penetrance_flag),
            any_conflict_flag=any(str(v).lower()=='true' for v in d.conflict_flag)))
    v=pd.DataFrame(variant).sort_values('allocated_new_raw_credit',ascending=False)
    v.to_csv(CARRIER_ATTRIBUTION_OUT/'variant_contribution_and_suspicion_review_internal.csv.gz',index=False)
    assert np.isclose(calls.weighted_credit.sum(),scores.B_raw_C.sum())
    assert np.isclose(calls.loc[calls.new_C_vs_A,'weighted_credit'].sum(),scores.loc[flags.C_only_vs_A,'B_raw_C'].sum())
    # Disclosure-safe standalone headline metrics, never output complementary
    # gene tables or sparse sensitivity carrier differences publicly.
    r=pd.DataFrame(rows);safe=r.loc[r.cohort.eq('case_control')&r.metric.isin(['C_only_vs_A','QC_PTV_added_vs_A','missense_added_B_vs_QC_PTV','added_four_genes_C_vs_B']),['metric','participants']].copy()
    safe.participants=safe.participants.map(lambda n:str(n) if n>=20 else '<20')
    safe.to_csv(CARRIER_ATTRIBUTION_ROOT/'results/table/46_rowlands_carrier_increment_summary.csv',index=False)
    summary={'status':'carrier_attribution_and_variant_review_completed','variant_calls_reconstruct_C_raw_burden':True,
        'allocation_rule':'split each participant/gene weight equally among retained normalized alleles; descriptive only',
        'manual_biological_adjudication_completed':False,'review_flags_are_not_automatic_exclusion_rules':True}
    (CARRIER_ATTRIBUTION_ROOT/'results/report/46_rowlands_carrier_attribution_validation.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('Internal transition, gene/variant contribution and suspicion-review tables complete.',flush=True)

# ============================================================================
# 2. summary: summarize_model_results
# ============================================================================
import ast
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

RESULT_SUMMARY_ROOT=Path(__file__).resolve().parents[1];RESULT_SUMMARY_BASE=RESULT_SUMMARY_ROOT/'data/29_rowlands_extension';RESULT_SUMMARY_OUT=RESULT_SUMMARY_BASE/'section11_models'

def summarize_model_results():
    status=json.loads((RESULT_SUMMARY_ROOT/'results/report/45_rowlands_models_validation.json').read_text())
    assert status['models_run'] and status['original_PRS311_repeat1_OOF_reproduced']
    cf=pd.read_csv(RESULT_SUMMARY_OUT/'model_coefficients_internal.csv')
    auc=pd.read_csv(RESULT_SUMMARY_OUT/'CV_AUC_summary_internal.csv').set_index('version')
    ds=pd.read_csv(RESULT_SUMMARY_OUT/'CV_delta_summary_internal.csv')
    ar=pd.read_csv(RESULT_SUMMARY_OUT/'age_explained_variance_internal.csv').set_index('version')
    scaling=pd.read_csv(RESULT_SUMMARY_OUT/'actual_model_scaling_internal.csv').set_index(['cohort','version'])
    versions=['A','A_prime','A_prime_QC','B','C','C_without_STK11','C_PTV_only','C_TP53_NF1_strict_VAF','C_ClinVar_min1star','C_without_CDH1']
    labels={'A':'A','A_prime':'A′','A_prime_QC':'A′QC','B':'B','C':'C','C_without_STK11':'C−STK11','C_PTV_only':'C：PTV-only',
        'C_TP53_NF1_strict_VAF':'C: Strict VAF','C_ClinVar_min1star':'C：ClinVar ≥1 star','C_without_CDH1':'C−CDH1'}
    rows=[]
    for v in versions:
        b=cf.loc[cf.version.eq(v)&cf.cohort.eq('case_control')&cf.term.eq('burden_z')].iloc[0]
        a=cf.loc[cf.version.eq(v)&cf.cohort.eq('diagnosis_age')&cf.term.eq('burden_z')].iloc[0]
        dd=ds.loc[ds.version.eq(v)&ds.comparison.eq('joint_minus_PRS311_baseline')].iloc[0]
        rows.append(dict(version=v,OR=b.OR,OR_CI_low=b.CI_low,OR_CI_high=b.CI_high,AUC_mean=auc.loc[v,'AUC_mean'],
            delta_AUC_mean=dd.delta_AUC_mean,delta_AUC_repeat1=dd.delta_AUC_repeat1,delta_AUC_repeat1_CI_low=dd.repeat1_CI_low,delta_AUC_repeat1_CI_high=dd.repeat1_CI_high,
            age_beta=a.beta,age_CI_low=a.CI_beta_low,age_CI_high=a.CI_beta_high,age_R_squared=ar.loc[v,'R_squared'],
            case_model_SD_rounded=scaling.loc[('case_control',v),'model_SD_ddof1'].round(2),
            age_model_SD_rounded=scaling.loc[('diagnosis_age',v),'model_SD_ddof1'].round(2)))
    s=pd.DataFrame(rows);s.to_csv(RESULT_SUMMARY_ROOT/'results/table/47_rowlands_interpretation_summary.csv',index=False)
    # Paired total gain versus legacy A, plus B versus the literal A-prime.
    from scripts.functions.statistics import delong_test
    ns = {'delong_test': delong_test}
    y=pd.read_csv(RESULT_SUMMARY_BASE/'section9/frozen_model_inputs.csv.gz').is_case.to_numpy(dtype=int)
    extra=[]
    for repeat in range(1,11):
        pred=np.load(RESULT_SUMMARY_OUT/f'cv_repeat_{repeat:02d}.npz')['predictions'];names=['PRS311_baseline']+status['versions']
        for a,b in [('C','A'),('B','A_prime')]:
            extra.append(dict(repeat=repeat,comparison=f'{a}_minus_{b}',**ns['delong_test'](y,pred[names.index(a)],pred[names.index(b)])))
    e=pd.DataFrame(extra);e.to_csv(RESULT_SUMMARY_OUT/'additional_total_AUC_comparisons_internal.csv',index=False)
    e.loc[e.repeat.eq(1)].to_csv(RESULT_SUMMARY_ROOT/'results/table/47_rowlands_total_AUC_comparison_repeat1.csv',index=False)
    c=s.set_index('version').loc['C'];a=s.set_index('version').loc['A']
    total=e.loc[e.repeat.eq(1)&e.comparison.eq('C_minus_A')].iloc[0]
    inter=cf.loc[cf.version.eq('C')&cf.cohort.eq('case_control')&cf.term.eq('PRS311_z_x_burden_z')].iloc[0]
    def table(vs):
        lines=['| Version | burden OR/SD (95% CI) | Mean OOF AUC | Mean ΔAUC relative to PRS311 | Within-case age coefficient, years/SD (95% CI) | Age R² |',
            '|---|---|---|---|---|---|']
        for v in vs:
            r=s.set_index('version').loc[v]
            lines.append(f'| {labels[v]} | {r.OR:.3f}（{r.OR_CI_low:.3f}–{r.OR_CI_high:.3f}） | {r.AUC_mean:.4f} | {r.delta_AUC_mean:.5f} | {r.age_beta:.2f}（{r.age_CI_low:.2f}to{r.age_CI_high:.2f}） | {r.age_R_squared*100:.2f}% |')
        return '\n'.join(lines)
    text=f'''# New burden extended analysis: result sorting and judgment

The main model and all scheduled sensitivity analyzes are now complete with N=141,204 case status (8,528 cases) and N=7,996 recorded age at diagnosis within cases. Fixed PRS311, case definition, sample, covariates and original 10×5 split. The case status and age coefficient of A reproduce the old nine-gene results; the original PRS311 first round OOF prediction also reproduces. All analyzes are extension/exploratory and current PTV biological manual adjudication is pending.

## Add carrier and source

Of C\'s 3,209 carriers, 487 were not recognized by old A. Attribution in mutually exclusive order: 148 became PTV carriers of the old nine genes under the current QC/PTV treatment, 239 became carriers only after the addition of qualified missense to the old nine genes, and 100 required the addition of four new genes, CDH1, NF1, PTEN, and STK11, to be identified. This is a breakdown of actual additions beyond the old A, without calling the full 487 missense deltas, or equating stage deltas with net headcount differentials.

The old A is different from the current PTV carrier set, so A′ retains the original carrier and only changes the weight; A′QC is a supplementary bridge. B−A′ contains both QC/PTV processing and missense changes, and B−A′QC only isolates missense.

The contribution of the new weighted score is distributed among multiple genes: BRCA2 is about 29%, BRCA1 is about 18%, and CHEK2 is about 16%. It is not covered by one new gene. The new variant with the largest contribution is CHEK2 missense, involving 60 new carriers (approximately 12%); this variant does not carry the above-mentioned PTV location/transcript review markers. Some other PTVs carry NMD, terminal exon or transcript markers; these are review tips and do not mean that they are determined to be invalid, nor can they claim that all suspected variant effects have been eliminated. The complete mutation-by-mutation sorting and labeling are kept internally, and those with less than 20 cells will not be made public.

## Main comparison on fixed sample{table(['A','A_prime','A_prime_QC','B','C'])}The baseline average OOF AUC of PRS311+original covariates is{auc.loc['PRS311_baseline','AUC_mean']:.4f}. The average AUC increment of C relative to A is{auc.loc['C','AUC_mean']-auc.loc['A','AUC_mean']:.5f},about{(auc.loc['C','AUC_mean']-auc.loc['A','AUC_mean'])*100:.2f}percentage points; the first round matching difference is{total.delta_auc:.5f}，95% DeLong CI {total.ci_low:.5f}–{total.ci_high:.5f}. Relative to the PRS311 baseline, the ΔAUC of the first round of C is{c.delta_AUC_repeat1:.5f}（{c.delta_AUC_repeat1_CI_low:.5f}–{c.delta_AUC_repeat1_CI_high:.5f}).

Decomposed average ΔAUC: changing weight A′−A is about 0.00008; current PTV processing A′QC−A′ is about 0.00034; missense B−A′QC is about 0.00120; new four genes C−B is about 0.00029. The main new discriminative gain comes from missense, but the overall gain magnitude is limited. Age R² is given by A\'s{a.age_R_squared*100:.2f}% raised to C{c.age_R_squared*100:.2f}%, an absolute improvement of approx.{(c.age_R_squared-a.age_R_squared)*100:.2f}percentage points.

The OR of case status is the conditional OR of burden per model SD when PRS is at the mean in the joint model with interactions, not the OR of carrier versus non-carrier. The within-case SD was used for the age coefficient. The actual SD and original unit coefficients of each version have been recorded in the internal comparison table. It cannot be concluded that the original scale effect is enhanced just because the OR/SD increases. CV means are reported separately from the first round CI; DeLong CI is a pairwise comparison of OOF predictions for that round and does not include the full uncertainty of retraining; the ten replicates are not ten independent cohorts. The absolute AUC includes the existing age observation window and technical covariates, focusing on explaining the increment under the same sample and covariates.

## Schedule sensitivity analysis{table(['C_without_STK11','C_PTV_only','C_TP53_NF1_strict_VAF','C_ClinVar_min1star','C_without_CDH1'])}The main estimates for de-STK11, strict TP53/NF1 VAF, one-star ClinVar and de-CDH1 are close to C, with the direction unchanged; PTV-only is still relevant, but the ΔAUC drops to about 0.01198, consistent with missense bringing some additional discriminative information. The above is evidence of stability and does not single out the final protocol based on individual P values, nor does it demonstrate that these rules adequately rule out clonal hematopoiesis.

## Interaction and paper presentation

PRS×burden interaction beta= for C{inter.beta:.5f}（95% CI {inter.CI_beta_low:.5f}to{inter.CI_beta_high:.5f}), interaction OR={inter.OR:.3f}（{inter.CI_low:.3f}–{inter.CI_high:.3f}). A mildly negative interaction on this logistic normalized product scale does not support positive multiplicative synergy; no biological antagonism can be inferred. The age model maintains the original setting without adding interaction. Both main effects are negative, which cannot prove age synergy.

Suggested paper statement: 13-gene PTV+qualified P/LP missense burden Using the approximate external weight of Rowlands all-PV, new carriers not identified by the old nine-gene PTV analysis are added; case status and diagnosis age correlation are recorded within cases in a fixed cohort, and a small additional discriminant gain is brought. The predetermined sensitivity results are generally consistent and can be used as the main result of extended analysis; it is not called a significant improvement in prediction ability, nor is it called a positive collaboration between PRS and burden. Before manual adjudication of variations, a description of the tentative definition must be retained.

The age results are always associated with documented age at diagnosis in breast cancer cases and cannot be interpreted as time-to-onset risk in the general population. External tables, pictures and texts continue to jointly implement small cell and deducible difference limits; accurate carrier/SD/variation tables remain internal.'''
    (RESULT_SUMMARY_ROOT/'results/report/47_rowlands_results_interpretation.md').write_text(text)
    # Exportable forest plot of whole-cohort estimates, no participant points.
    plt.style.use('default');fig,axes=plt.subplots(1,2,figsize=(10,5.5))
    english=['A','A-prime','A-prime QC','B','C','C minus STK11','C PTV-only','C strict VAF','C ClinVar >=1 star','C minus CDH1']
    yy=np.arange(len(s));axes[0].errorbar(s.OR,yy,xerr=np.vstack([s.OR-s.OR_CI_low,s.OR_CI_high-s.OR]),fmt='o',color='tab:blue',capsize=3)
    axes[1].errorbar(s.age_beta,yy,xerr=np.vstack([s.age_beta-s.age_CI_low,s.age_CI_high-s.age_beta]),fmt='o',color='tab:orange',capsize=3)
    axes[0].set_yticks(yy,english);axes[1].set_yticks(yy,[])
    for ax in axes:ax.invert_yaxis();ax.grid(False)
    axes[0].set_xlabel('Case-status burden OR / model SD\nConditional on mean PRS; 95% HC1 CI')
    axes[1].set_xlabel('Recorded diagnosis-age beta (years / case SD)\nWithin breast cancer cases; 95% HC1 CI')
    fig.suptitle('Fixed-cohort extension and prespecified sensitivities\nProvisional variant definition; approximate external weights',fontsize=11)
    fig.tight_layout();fig.savefig(RESULT_SUMMARY_ROOT/'results/figures/47_rowlands_burden_model_comparison.svg');plt.close(fig)
    validation={'status':'results_organized_and_interpreted','models_run':True,'carrier_attribution_completed':True,
        'A_original_case_and_age_coefficients_reproduced':True,'no_significance_based_scheme_selection':True,
        'main_interaction':'negative_on_logistic_product_scale_not_positive_synergy','final_variant_definition_frozen':False,
        'public_small_cells_and_exact_comparative_SD_not_disclosed':True}
    (RESULT_SUMMARY_ROOT/'results/report/47_rowlands_result_interpretation_validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    print('Results report, paired total AUC comparisons, disclosure-safe summary and forest plot complete.',flush=True)

# ============================================================================
# 3. age-figures: plot_age_and_auc_results
# ============================================================================
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from threadpoolctl import threadpool_limits

AGE_FIGURES_ROOT=Path(__file__).resolve().parents[1];AGE_FIGURES_BASE=AGE_FIGURES_ROOT/'data/29_rowlands_extension';AGE_FIGURES_OUT=AGE_FIGURES_BASE/'section12_figures'
FIG=AGE_FIGURES_ROOT/'results/figures';TABLE=AGE_FIGURES_ROOT/'results/table';Z=1.959963984540054
PC=[f'PC{i}' for i in range(1,11)]
V=['A','A_prime','B','C','C_without_STK11']
LABEL=['9 genes: PTV\nOriginal weights','9 genes: PTV\nRowlands weights',
       '9 genes: PTV + qualifying P/LP missense\nRowlands weights',
       '13 genes: PTV + qualifying P/LP missense\nRowlands weights',
       '12 genes (without STK11): PTV + qualifying P/LP missense\nRowlands weights']

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def plot_age_and_auc_results():
    AGE_FIGURES_OUT.mkdir(parents=True,exist_ok=True)
    import subprocess
    assert subprocess.run(['git','check-ignore','-q',str(AGE_FIGURES_OUT)],cwd=AGE_FIGURES_ROOT).returncode==0
    inputs=[AGE_FIGURES_BASE/'section9/frozen_model_inputs.csv.gz',AGE_FIGURES_BASE/'section8/comparable_scores_case_control.csv.gz',
        AGE_FIGURES_BASE/'section11_models/model_coefficients_internal.csv',AGE_FIGURES_BASE/'section11_models/CV_delta_summary_internal.csv',
        AGE_FIGURES_BASE/'section11_models/carrier_transition_flags_internal.csv.gz']
    hashes={str(p.relative_to(AGE_FIGURES_ROOT)):sha(p) for p in inputs}
    all_df=pd.read_csv(inputs[0]).set_index('person_id')
    ids=pd.read_csv(AGE_FIGURES_BASE/'diagnosis_age_ids.csv.gz').person_id
    df=all_df.loc[ids].copy()
    scores=pd.read_csv(inputs[1]).set_index('sample_id')
    carrier=scores.loc[df.index,'carried_gene_count_C'].gt(0).astype(float)
    assert len(df)==7996 and df.is_case.eq(1).all() and int(carrier.sum())==543
    assert int(carrier.sum())>=20 and int((1-carrier).sum())>=20
    prs=df.PRS311;mu=prs.mean();sd=prs.std(ddof=1);p=(prs-mu)/sd
    X=pd.DataFrame({'const':1.,'PRS311_z':p,'C_carrier':carrier,'PRS311_z_x_C_carrier':p*carrier},index=df.index).join(df[PC])
    with threadpool_limits(limits=1):fit=sm.OLS(df.age_at_diagnosis,X).fit(cov_type='HC1')
    cov=fit.cov_params().to_numpy();beta=fit.params.to_numpy()
    assert np.isfinite(beta).all() and np.isfinite(cov).all()
    coeff=pd.DataFrame({'term':X.columns,'beta':beta,'robust_se_HC1':fit.bse.to_numpy(),
        'CI_low':fit.conf_int()[0].to_numpy(),'CI_high':fit.conf_int()[1].to_numpy(),'p_value':fit.pvalues.to_numpy()})
    coeff.to_csv(AGE_FIGURES_OUT/'age_binary_carrier_interaction_all_coefficients_internal.csv',index=False)
    contrasts=[]
    for name,terms in [('Non-carrier PRS slope',{'PRS311_z':1}),('Carrier PRS slope',{'PRS311_z':1,'PRS311_z_x_C_carrier':1}),
        ('Carrier minus non-carrier slope (interaction)',{'PRS311_z_x_C_carrier':1})]:
        c=np.zeros(len(beta))
        for term,value in terms.items():c[X.columns.get_loc(term)]=value
        estimate=c@beta;se=np.sqrt(c@cov@c)
        contrasts.append(dict(contrast=name,beta_years_per_PRS_SD=estimate,robust_se_HC1=se,CI_low=estimate-Z*se,CI_high=estimate+Z*se,
            beta_years_per_raw_PRS_unit=estimate/sd,raw_CI_low=(estimate-Z*se)/sd,raw_CI_high=(estimate+Z*se)/sd))
    ct=pd.DataFrame(contrasts)
    ct.to_csv(AGE_FIGURES_OUT/'age_binary_carrier_PRS_slopes_internal.csv',index=False)
    ct.drop(columns=['beta_years_per_raw_PRS_unit','raw_CI_low','raw_CI_high']).to_csv(TABLE/'48_age_binary_carrier_PRS_slopes.csv',index=False)
    # Common support: central 95% PRS range of each group, intersected.
    ranges=[p.loc[carrier.eq(g)].quantile([.025,.975]).to_numpy() for g in [0,1]]
    lo=max(r[0] for r in ranges);hi=min(r[1] for r in ranges);assert lo<hi
    grid=np.linspace(lo,hi,160);pred=[]
    for group in [0,1]:
        xp=pd.DataFrame({'const':np.ones(len(grid)),'PRS311_z':grid,'C_carrier':float(group),
            'PRS311_z_x_C_carrier':grid*group})
        for pc in PC:xp[pc]=df[pc].mean()
        xp=xp[X.columns].to_numpy();expected=xp@beta
        se=np.sqrt(np.maximum(np.einsum('ij,jk,ik->i',xp,cov,xp),0))
        for j in range(len(grid)):pred.append(dict(C_carrier=group,PRS311_case_z=grid[j],predicted_age=expected[j],CI_low=expected[j]-Z*se[j],CI_high=expected[j]+Z*se[j]))
    pd.DataFrame(pred).to_csv(AGE_FIGURES_OUT/'adjusted_age_prediction_grid_internal.csv',index=False)
    # Checks of group slopes and of uncertainty using full robust covariance.
    for g in [0,1]:
        q=pd.DataFrame(pred).loc[lambda d:d.C_carrier.eq(g)]
        slope=(q.predicted_age.iloc[-1]-q.predicted_age.iloc[0])/(grid[-1]-grid[0])
        assert np.isclose(slope,ct.beta_years_per_PRS_SD.iloc[g])
        assert (q.CI_low<=q.predicted_age).all() and (q.CI_high>=q.predicted_age).all()
    plt.style.use('default')
    colors=['tab:blue','tab:orange']
    fig,axes=plt.subplots(1,2,figsize=(16.8,5.8),gridspec_kw={'width_ratios':[1,1.35]})
    age=pd.read_csv(inputs[2]);age=age.loc[age.cohort.eq('diagnosis_age')&age.term.eq('burden_z')].set_index('version').loc[V]
    yy=np.arange(len(V));axes[0].errorbar(age.beta,yy,xerr=np.vstack([age.beta-age.CI_beta_low,age.CI_beta_high-age.beta]),fmt='o',color='tab:blue',capsize=4)
    axes[0].axvline(0,color='0.55',linestyle='--',linewidth=1);axes[0].set_yticks(yy,LABEL);axes[0].invert_yaxis()
    axes[0].tick_params(axis='y',labelsize=9)
    axes[0].set_xlabel('Difference in recorded diagnosis age\n(years per SD burden)');axes[0].set_title('Burden definitions and weights',loc='left',fontweight='bold')
    for g in [0,1]:
        q=pd.DataFrame(pred).loc[lambda d:d.C_carrier.eq(g)]
        axes[1].plot(q.PRS311_case_z,q.predicted_age,color=colors[g],label='Non-carrier' if g==0 else '13-gene carrier')
        axes[1].fill_between(q.PRS311_case_z,q.CI_low,q.CI_high,color=colors[g],alpha=.18)
    axes[1].set_xlabel('Continuous PRS311 (SD within fixed cases)');axes[1].set_ylabel('PC-adjusted predicted recorded diagnosis age (years)')
    axes[1].set_title('PRS × 13-gene carrier status (exploratory)',loc='left',fontweight='bold')
    annotation=[]
    for r in ct.itertuples():
        name={'Non-carrier PRS slope':'Non-carrier slope','Carrier PRS slope':'Carrier slope','Carrier minus non-carrier slope (interaction)':'Interaction (slope difference)'}[r.contrast]
        annotation.append(f'{name}: {r.beta_years_per_PRS_SD:.2f} [{r.CI_low:.2f}, {r.CI_high:.2f}]')
    axes[1].text(.02,.98,'Years per PRS SD; 95% HC1 CI\n'+'\n'.join(annotation),transform=axes[1].transAxes,va='top',fontsize=9,
        bbox={'facecolor':'white','alpha':.88,'edgecolor':'0.85'})
    prediction_frame=pd.DataFrame(pred)
    axes[1].set_ylim(prediction_frame.CI_low.min()-1.1,prediction_frame.CI_high.max()+2.5)
    axes[1].legend(loc='lower left',frameon=False)
    for ax in axes:ax.grid(False)
    fig.text(.5,.01,'Fixed breast cancer cases (N = 7,996). Burden SD differs by definition. Right panel: additional binary-carrier interaction model.',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.05,1,1]);fig.savefig(FIG/'48_main_recorded_diagnosis_age.svg');plt.close(fig)
    # Exact disjoint counts can be displayed only if all cells (including
    # potentially inferable lost carriers) meet the actual disclosure rule.
    f=pd.read_csv(inputs[4]).set_index('sample_id');shared=int((f.A_carrier&f.C_carrier).sum());lost=int((f.A_carrier&~f.C_carrier).sum())
    names=['Original nine-gene PTV','PTV eligibility / call handling changes','Added qualifying missense (nine genes)','Added CDH1, NF1, PTEN, STK11']
    count=[shared]+[int(f[k].sum()) for k in ['new_C_current_nine_PTV','new_C_nine_missense','new_C_added_four_only']]
    assert lost==0,'A-only carriers need a separate suppression-safe diagram'
    assert all(n>=20 for n in count) and sum(count)==int(f.C_carrier.sum()) and sum(count[1:])==487
    fig,ax=plt.subplots(figsize=(11,3.4));palette=['0.7','tab:purple','tab:blue','tab:orange'];left=0
    for name,n,col in zip(names,count,palette):
        ax.barh(0,n,left=left,color=col,label=f'{name}: {n:,}',height=.5)
        if n>=200:ax.text(left+n/2,0,f'{n:,}',ha='center',va='center',color='black' if col=='0.7' else 'white')
        left+=n
    ax.set_yticks([0],['13-gene carriers']);ax.set_xlabel('Participants, mutually exclusive attribution')
    ax.set_title('Carrier expansion: 487 newly identified carriers',loc='left',fontweight='bold');ax.grid(False)
    ax.legend(loc='upper center',bbox_to_anchor=(.5,-.26),ncol=2,frameon=False)
    fig.tight_layout();fig.savefig(FIG/'48_supp_carrier_increment.svg',bbox_inches='tight');plt.close(fig)
    # Keep point and confidence interval from the same repeat, not the mean
    # of ten repeated AUC values paired with a first-repeat confidence interval.
    dd=pd.read_csv(inputs[3]);dd=dd.loc[dd.comparison.eq('joint_minus_PRS311_baseline')].set_index('version').loc[V[:4]]
    fig,ax=plt.subplots(figsize=(10.4,4.8));yy=np.arange(4)
    ax.errorbar(dd.delta_AUC_repeat1,yy,xerr=np.vstack([dd.delta_AUC_repeat1-dd.repeat1_CI_low,dd.repeat1_CI_high-dd.delta_AUC_repeat1]),fmt='o',capsize=4,color='tab:blue')
    ax.axvline(0,color='0.55',linestyle='--',linewidth=1);ax.set_yticks(yy,LABEL[:4]);ax.invert_yaxis();ax.grid(False)
    ax.tick_params(axis='y',labelsize=9)
    ax.set_xlabel('OOF ΔAUC relative to covariates + PRS311');ax.set_title('Fixed cohort and matched CV folds',loc='left')
    fig.text(.5,.01,'Points and paired DeLong 95% CI: repeat 1 of fixed 10 × 5 CV.',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.05,1,1]);fig.savefig(FIG/'48_supp_delta_AUC.svg');plt.close(fig)
    # Exact supplementary distribution table stays internal; public display
    # precision is chosen so sensitivity carrier differences cannot be recovered.
    dist=pd.read_csv(AGE_FIGURES_BASE/'section10/internal_delivery/03_A_to_C_carrier_and_model_comparison.csv')
    dist=dist.loc[dist.version.isin(V),['cohort','version','n','carrier_fraction','zero_fraction','raw_mean','SD_population_ddof0','SD_sample_ddof1',
        'raw_p25','raw_p50','raw_p75','raw_p95','raw_p99']]
    dist.to_csv(AGE_FIGURES_OUT/'burden_distribution_carrier_supplement_exact_internal.csv',index=False)
    pub=dist.copy()
    for ix,r in pub.iterrows():
        precision=1 if r.cohort=='case_control' else 0
        pub.loc[ix,'carrier_fraction']=round(r.carrier_fraction*100,precision)
        pub.loc[ix,'zero_fraction']=round(r.zero_fraction*100,precision)
    pub=pub.rename(columns={'carrier_fraction':'carrier_percent_coarsened','zero_fraction':'zero_percent_coarsened'})
    definitions=dict(zip(V,[x.replace('\n','; ') for x in LABEL]))
    pub['burden_definition']=pub.version.map(definitions)
    pub=pub.drop(columns='version')
    numeric=[c for c in pub.columns if c.startswith('raw_') or c.startswith('SD_')]
    pub[numeric]=pub[numeric].round(2)
    pub.to_csv(TABLE/'48_burden_distribution_carrier_supplement_coarsened.csv',index=False)
    # Add the new binary-carrier age model explicitly to the internal delivery
    # sheet; these columns do not replace the main continuous-burden estimates.
    comparison_path=AGE_FIGURES_BASE/'section10/internal_delivery/03_A_to_C_carrier_and_model_comparison.csv'
    comparison=pd.read_csv(comparison_path)
    target=comparison.cohort.eq('diagnosis_age')&comparison.version.eq('C')
    assert target.sum()==1
    for j,prefix in enumerate(['age_binary_model_PRS_noncarrier','age_binary_model_PRS_carrier','age_binary_model_interaction']):
        for field,col in [('beta','beta_years_per_PRS_SD'),('CI_low','CI_low'),('CI_high','CI_high')]:
            comparison.loc[target,f'{prefix}_{field}']=ct.iloc[j][col]
    comparison.to_csv(comparison_path,index=False)
    delivery_manifest_path=AGE_FIGURES_BASE/'section10/delivery_manifest_internal.json'
    delivery_manifest=json.loads(delivery_manifest_path.read_text())
    delivery_manifest.setdefault('post_main_model_delivery_sha256',delivery_manifest['delivery'].copy())
    delivery_manifest['delivery']={p.name:sha(p) for p in comparison_path.parent.iterdir() if p.is_file()}
    delivery_manifest['additional_exploratory_age_binary_carrier_model']='data/29_rowlands_extension/section12_figures/figure_and_model_manifest_internal.json'
    delivery_manifest_path.write_text(json.dumps(delivery_manifest,indent=2)+'\n')
    ctrow=ct.iloc[2]
    caption=f'''# Diagnostic age priority chart set: legends and explanations

Presentation order: Supplementary carrier delta plot → Age at diagnosis two-panel main plot → Supplementary ΔAUC plot. The text in the figure is in English to facilitate the export of the paper.

**Program definition**: The original nine genes are ATM, BARD1, BRCA1, BRCA2, CHEK2, PALB2, RAD51C, RAD51D, and TP53. The 13-gene program adds CDH1, NF1, PTEN, and STK11 on this basis; the 12-gene sensitivity program only removes STK11. PTV refers to protein truncating variants, and qualified missense refers to missense variants that pass this program\'s germline P/LP, ClinVar at least 2 stars, low penetrance/risk allele/conflict elimination, related diseases and inheritance patterns, frequency of European ancestry, and calling QC. The Rowlands weight comes from the gene OR of PTV+non-PTV-PV, which is the approximate external weight of the PTV+qualified P/LP missense score in this study. The mutation definitions of the two are not exactly the same.

**Carrier incremental chart**: 3,209 carriers in the 13-gene PTV + qualified P/LP missense program, which are mutually exclusive attributions based on those identified by the original nine-gene PTV (2,722 people), those newly added to the current PTV inclusion/calling process (148 people), those newly added with the nine-gene qualified missense (239 people), and those identified after the addition of CDH1/NF1/PTEN/STK11 (100 people). A total of 487 new people were added, and the difference in PTV processing cannot be counted as a missense increase. All the original nine-gene PTV carriers are retained in the 13-gene plan, and each display unit has at least 20 people; no small gene stratification is drawn.

**Left side of main picture of diagnostic age: burden definition and weight comparison**. Among the fixed 7,996 breast cancer cases, the following comparisons were made sequentially: the original nine-gene PTV + original weight; the same original nine-gene PTV carrier + Rowlands weight; the current nine-gene PTV + qualified P/LP missense + Rowlands weight; the 13-gene PTV + qualified P/LP missense + Rowlands weight; the 12-gene scheme without STK11 + Rowlands weight. The OLS models recording age at diagnosis also included PRS311, adjusted PC1–10, and used HC1 robust standard errors. The points are the burden coefficients, the lines are the 95% CIs, and the zero value has a reference line. Each burden is standardized by its own SD within the case cohort, and the SD of each scheme is different; the original scale effect cannot be judged based solely on the change in the coefficient of each SD. The precise distribution, SD, and carrier ratio are left in the internal supplementary table, and a jointly coarsened version is provided externally.

**Right side of main picture of age at diagnosis: interaction between PRS311 and 13 gene carrier status**. The new expanded model is `record diagnosis age ~ PRS311 standard score + 13 gene carrier index + PRS311 standard score × 13 gene carrier index + PC1–10`, using the same cases and HC1. Carrier is defined as a case carrying at least one 13-gene scheme-eligible PTV or P/LP missense. This binary indicator replaces the continuous burden in this extended model; the original age main model maintains the no-interaction setting. The horizontal axis is the standard score of continuous PRS311 in fixed cases. The two lines set the PC to the mean of the same fixed case cohort; under the linear model, it is equivalent to marginal standardization on the same PC distribution. Shaded are point-wise 95% HC1 confidence intervals predicting the mean recorded age at diagnosis, not individual prediction intervals or simultaneous confidence bands. The horizontal axis uses the intersection of the central 95% PRS ranges of each of the two groups. The slopes and slope differences of the two groups of PRS on the figure are calculated using the complete robust covariance matrix; the vertical distance between the two lines mainly reflects the carrier main effect and cannot be used to claim synergy.

The within-case binary carrier interaction (slope difference) is{ctrow.beta_years_per_PRS_SD:.3f}Year/PRS SD (95% CI{ctrow.CI_low:.3f}to{ctrow.CI_high:.3f}). This is a new age-related interaction estimate. It belongs to a different model from the continuous PRS × burden interaction in the previous case status logistic and cannot be mixed.

**ΔAUC Supplementary Figure**: fixed 141,204 people, same PRS311 and covariates, original 10×5 discount. The original nine-gene PTV+original weight, the same PTV carrier+Rowlands weight, the current nine-gene PTV+qualified missense+Rowlands weight, the 13-gene PTV+qualified missense+Rowlands weight, and the joint model OOF ΔAUC relative to covariates+PRS311 baseline are shown in sequence. The joint model adds burden and the existing PRS×burden interaction. Both points and CI come from the first round of OOF predictions, and CI is paired DeLong; the average and range of 10 repetitions are listed separately in the existing supplementary table. The CI condition is for this round of fitting predictions and does not include the uncertainty of complete model retraining; the 10 repetitions are not regarded as independent samples. The difference between the former nine-gene combined scheme and the original nine-gene PTV weighting scheme includes both PTV processing and missense changes; when isolating the missense contribution, it should be compared with the nine-gene PTV-only bridging scheme using the same current QC.

All age results are documented age-at-diagnosis associations among breast cancer cases only and are not interpreted as time-to-onset risks in the general population. The variant biology review is still pending and the panel is an extended analysis under the tentative definition.

Publicly available data only include overall model/linear comparisons and coarsened supplementary statistics, without participant-by-participant data points. Subtotals, percentages, and small cells that can be recovered from differences must be checked simultaneously. Basis: [All of Us Disclosure Statement](https://support.researchallofus.org/hc/en-us/articles/360043016291-How-to-comply-with-the-All-of-Us-Data-and-Statistics-Dissemination-Policy).'''
    (AGE_FIGURES_ROOT/'results/report/48_age_priority_figure_captions.md').write_text(caption)
    for path,h in hashes.items():assert sha(AGE_FIGURES_ROOT/path)==h,path
    manifest={'status':'age_priority_figures_and_exploratory_age_carrier_interaction_complete','fixed_case_n':7996,
        'fixed_cases_C_carriers':543,'fixed_cases_noncarriers':7453,'main_age_model_unchanged':True,
        'additional_age_model':'PRS311_case_z + C_binary_carrier + interaction + PC1-10; OLS HC1',
        'slope_difference_CI_uses_full_robust_covariance':True,'line_slopes_numerically_verified':True,
        'AUC_points_and_CI_same_repeat':True,'carrier_attribution_mutually_exclusive':True,
        'carrier_cells_all_ge20_or_zero':True,'public_exact_sparse_score_differences_not_reported':True,
        'inputs':hashes,'PRS311_age_scaling':{'mean':mu,'SD_ddof1':sd},'x_range_case_PRS_z':[lo,hi],
        'final_variant_definition_frozen':False}
    (AGE_FIGURES_OUT/'figure_and_model_manifest_internal.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (AGE_FIGURES_ROOT/'results/report/48_age_priority_figures_validation.json').write_text(json.dumps({k:v for k,v in manifest.items() if k not in {'inputs','PRS311_age_scaling','x_range_case_PRS_z'}},indent=2)+'\n')
    print('Age-priority main figure, carrier and AUC supplements, contrasts, captions and validation complete.',flush=True)


# ============================================================================
# 4. score-distributions: plot_raw_score_distributions
# ============================================================================
"""Plot internal-review raw score distributions using the frozen CC cohort."""
from pathlib import Path
import json
import hashlib

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
from threadpoolctl import threadpool_limits

MS_4_ROOT = Path(__file__).resolve().parents[1]
MS_4_OUT = MS_4_ROOT / 'manuscript/figures/internal_score_distributions'
MS_4_BLUE = '#2463B4'
MS_4_RED = '#D54848'


def MS_4_density(values, grid):
    """Exact full-sample Gaussian KDE via duplicate-value compression."""
    values = np.asarray(values, dtype=float)
    unique, counts = np.unique(values, return_counts=True)
    bandwidth = values.std(ddof=1) * len(values) ** (-1 / 5)
    assert bandwidth > 0
    z = (grid[:, None] - unique[None, :]) / bandwidth
    return (np.exp(-.5 * z * z) @ (counts / counts.sum())) / (bandwidth * np.sqrt(2 * np.pi))


def MS_4_save(fig, stem):
    fig.savefig(MS_4_OUT / f'{stem}.svg', bbox_inches='tight')
    png = MS_4_OUT / f'{stem}.png'
    jpg = MS_4_OUT / f'{stem}.jpg'
    fig.savefig(png, dpi=170, bbox_inches='tight')
    fig.savefig(jpg, dpi=170, bbox_inches='tight', pil_kwargs={'quality': 85, 'subsampling': 2})
    # PNG is the inline preview; JPEG is only a temporary size comparison.
    smallest = '.png' if png.stat().st_size <= jpg.stat().st_size else '.jpg'
    if smallest == '.png':
        jpg.unlink()
    else:
        png.unlink()
    plt.close(fig)
    return smallest


def plot_raw_score_distributions():
    MS_4_OUT.mkdir(parents=True, exist_ok=True)
    sources = [MS_4_ROOT / 'manuscript/datas/burden_13gene_case_control.csv',
               MS_4_ROOT / 'manuscript/datas/PRS_scores.csv',
               MS_4_ROOT / 'data/29_rowlands_extension/section9/frozen_model_inputs.csv.gz']
    d = pd.read_csv(sources[2]).set_index('person_id')
    b = pd.read_csv(sources[0]).set_index('sample_id')
    p = pd.read_csv(sources[1]).set_index('person_id')
    assert d.index.is_unique and b.index.is_unique and p.index.is_unique
    assert set(d.index) == set(b.index) and set(d.index) <= set(p.index)
    b = b.loc[d.index]; p = p.loc[d.index]
    assert len(d) == 141204 and int(d.is_case.sum()) == 8528
    assert np.allclose(p.PRS311, d.PRS311, rtol=1e-12, atol=1e-12)
    versions = [
        ('12 genes\nWithout STK11', 'B_raw_12gene_without_STK11_main'),
        ('13 genes\nMain', 'B_raw_13gene_main'),
        ('12 genes\nStrict TP53/NF1 VAF', 'B_raw_12gene_without_STK11_tp53_nf1_vaf_sensitive'),
        ('13 genes\nStrict TP53/NF1 VAF', 'B_raw_13gene_tp53_nf1_vaf_sensitive'),
    ]
    assert np.isfinite(b[[col for _, col in versions]].to_numpy()).all()
    labels = [(1, 'Case', MS_4_BLUE, -.19), (0, 'Control', MS_4_RED, .19)]
    plt.style.use('default')
    fig, axes = plt.subplots(2, 1, figsize=(11, 8.2), sharex=True)
    with threadpool_limits(limits=1):
        for ax, positive_only in zip(axes, [False, True]):
            for j, (_, column) in enumerate(versions):
                for case, label, color, offset in labels:
                    values = b.loc[d.is_case.eq(case), column].to_numpy()
                    if positive_only:
                        values = values[values > 0]
                    assert len(values) >= 20
                    grid = np.linspace(values.min(), values.max(), 350)
                    kde = MS_4_density(values, grid)
                    widths = .16 * kde / kde.max()
                    center = j + offset
                    ax.fill_betweenx(grid, center - widths, center + widths, color=color, alpha=.75,
                                     linewidth=.8, edgecolor=color)
                    median = np.median(values)
                    ax.plot([center - .045, center + .045], [median, median], color='black', linewidth=1.5)
            ax.set_ylabel('Raw weighted burden score')
            ax.grid(False)
            ax.spines[['top', 'right']].set_visible(False)
            ax.set_ylim(bottom=-.1)
        axes[0].set_title('A  All participants, including zero scores', loc='left', fontsize=12)
        axes[1].set_title('B  Non-zero scores only (carriers)', loc='left', fontsize=12)
        axes[0].legend(handles=[Patch(facecolor=MS_4_BLUE, alpha=.75, label='Case'),
                                Patch(facecolor=MS_4_RED, alpha=.75, label='Control')], frameon=False, loc='upper right')
        axes[1].set_xticks(range(len(versions)), [label for label, _ in versions])
        axes[1].set_xlabel('Burden definition / sensitivity analysis')
        fig.suptitle('13-gene PTV + eligible P/LP missense burden', fontsize=14)
        fig.text(.5, .015, 'Same frozen case-control cohort. Violin widths normalized separately; black bars denote medians.\n'
                 'Zero-heavy scores concentrate near 0 in panel A; panel B shows the carrier distribution.', ha='center', fontsize=9)
        fig.tight_layout(rect=[0, .065, 1, .95])
        burden_preview = MS_4_save(fig, 'burden_raw_score_violin')
        fig, ax = plt.subplots(figsize=(8, 5))
        values = p.PRS311.to_numpy(dtype=float)
        assert np.isfinite(values).all()
        grid = np.linspace(values.min(), values.max(), 650)
        for case, label, color, _ in labels:
            v = p.loc[d.is_case.eq(case), 'PRS311'].to_numpy(dtype=float)
            kde = gaussian_kde(v)(grid)
            ax.plot(grid, kde, color=color, linewidth=2, label=label)
            ax.fill_between(grid, 0, kde, color=color, alpha=.12)
        ax.set_xlabel('Raw PRS311 score'); ax.set_ylabel('Probability density')
        ax.set_title('PRS311 distribution in cases and controls')
        ax.legend(frameon=False); ax.grid(False); ax.spines[['top', 'right']].set_visible(False)
        ax.set_ylim(bottom=0); fig.tight_layout()
        prs_preview = MS_4_save(fig, 'PRS311_raw_score_KDE')
    manifest = {'cohort': 'frozen case-control', 'n': len(d), 'case_n': int(d.is_case.sum()),
                'scores': 'raw, not z-standardized', 'burden_versions': [c for _, c in versions],
                'burden_KDE': 'Gaussian Scott bandwidth of full selected sample; duplicates compressed exactly; display limited to score support',
                'PRS_KDE': 'Gaussian Scott bandwidth separately by group; each density integrates to one over real line',
                'violin_width': 'same maximum width for every group, not scaled by sample size',
                'burden_preview_extension': burden_preview, 'PRS_preview_extension': prs_preview,
                'internal_only': True, 'models_refitted': False,
                'input_sha256': {str(q.relative_to(MS_4_ROOT)): hashlib.sha256(q.read_bytes()).hexdigest() for q in sources}}
    (MS_4_OUT / 'plot_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print('Burden violin and PRS311 KDE generated from raw scores; participant data not printed.')

# ============================================================================
# 5. discrimination: plot_four_model_discrimination
# ============================================================================
"""Complete missing four-model CV predictions and plot C sensitivities."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import hashlib
import json
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import RepeatedStratifiedKFold
from threadpoolctl import threadpool_limits

MS_5_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MS_5_ROOT))
from scripts.functions.statistics import auc_ci, delong_test

MS_5_BASE = MS_5_ROOT / 'data/29_rowlands_extension'
MS_5_EXISTING = MS_5_BASE / 'section11_models'
MS_5_CACHE = MS_5_ROOT / 'data/manuscript_refresh_20261006/four_model_CV'
MS_5_FIG = MS_5_ROOT / 'manuscript/figures'
MS_5_DATA = MS_5_ROOT / 'manuscript/datas/case_control'
MS_5_VERSIONS = [
    ('C', '13 genes: main'),
    ('C_without_STK11', '12 genes: without STK11'),
    ('C_PTV_only', '13 genes: PTV only'),
    ('C_TP53_NF1_strict_VAF', '13 genes: strict TP53/NF1 VAF'),
    ('C_ClinVar_min1star', '13 genes: ClinVar >=1 star'),
    ('C_without_CDH1', '12 genes: without CDH1'),
]
MS_5_MODELS = ['Covariates only', 'PRS + covariates', 'Burden + covariates', 'Joint']
MS_5_COLORS = ['#777777', '#E58A27', '#2768B4', '#329461']


def MS_5_sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def MS_5_design(cov, burden=None):
    cols = [np.ones(len(cov))]
    if burden is not None:
        cols.append(burden)
    return np.column_stack(cols + [cov])


def MS_5_fit_predict(y, train_X, test_X):
    fit = sm.GLM(y, train_X, family=sm.families.Binomial()).fit()
    assert fit.converged and np.isfinite(fit.params).all()
    return fit.predict(test_X)


def MS_5_save_figure(fig, name):
    fig.savefig(MS_5_FIG / f'{name}.svg', bbox_inches='tight')
    png, jpg = MS_5_FIG / f'{name}.png', MS_5_FIG / f'{name}.jpg'
    fig.savefig(png, dpi=180, bbox_inches='tight')
    fig.savefig(jpg, dpi=180, bbox_inches='tight', pil_kwargs={'quality': 85, 'subsampling': 2})
    smallest = png if png.stat().st_size <= jpg.stat().st_size else jpg
    (jpg if smallest == png else png).unlink()
    plt.close(fig)
    return str(smallest.relative_to(MS_5_ROOT))


def plot_four_model_discrimination():
    for path in [MS_5_CACHE, MS_5_FIG, MS_5_DATA]:
        path.mkdir(parents=True, exist_ok=True)
    paths = [MS_5_BASE / 'section9/frozen_model_inputs.csv.gz',
             MS_5_BASE / 'section8/comparable_scores_case_control.csv.gz',
             MS_5_BASE / 'section10/prespecified_sensitivity_scores_case_control.csv.gz',
             MS_5_ROOT / 'data/16_joint_model_cv_oof_person_level.csv.gz',
             MS_5_EXISTING / 'model_input_lock_internal.json',
             MS_5_ROOT / 'results/report/45_rowlands_models_validation.json',
             MS_5_ROOT / 'results/table/45_rowlands_CV_AUC.csv',
             MS_5_ROOT / 'scripts/functions/statistics.py', Path(__file__)]
    paths += [MS_5_EXISTING / f'cv_repeat_{r:02d}.npz' for r in range(1, 11)]
    hashes = {str(p.relative_to(MS_5_ROOT)): MS_5_sha(p) for p in paths}
    run_hash = hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest()
    historical_lock = json.loads(paths[4].read_text())
    # Existing model data must match their original fitting inputs exactly.
    for p, digest in historical_lock['inputs'].items():
        if p.startswith('data/'):
            assert MS_5_sha(MS_5_ROOT / p) == digest, 'Existing model input changed: ' + p
    status = json.loads(paths[5].read_text())
    existing_names = ['PRS311_baseline'] + status['versions']
    d = pd.read_csv(paths[0]).set_index('person_id')
    old = pd.read_csv(paths[3]).set_index('person_id')
    assert d.index.equals(old.index) and len(d) == 141204
    scores = pd.read_csv(paths[1]).set_index('sample_id').join(pd.read_csv(paths[2]).set_index('sample_id')).loc[d.index]
    raw = scores[[f'B_raw_{v}' for v, _ in MS_5_VERSIONS]].to_numpy(dtype=float)
    y = d.is_case.to_numpy(dtype=int)
    cov = pd.concat([d[['age_at_last_obs'] + [f'PC{i}' for i in range(1, 11)]],
                     pd.get_dummies(d.seq_site_id, prefix='site', drop_first=True).astype(float),
                     pd.get_dummies(d.seq_sequencer, prefix='seq', drop_first=True).astype(float)], axis=1).to_numpy(dtype=float)
    assert y.sum() == 8528 and np.isfinite(cov).all() and np.isfinite(raw).all()
    strata = y.astype(str) + '_' + old.carrier_main_10gene.to_numpy(dtype=int).astype(str)
    splits = list(RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=20260818).split(np.arange(len(y)), strata))
    fold = np.full(len(y), -1)
    for k, (_, te) in enumerate(splits[:5]):
        fold[te] = k
    assert np.array_equal(fold, old.fold_r1)
    (MS_5_CACHE / 'input_lock_internal.json').write_text(json.dumps({'inputs': hashes, 'run_hash': run_hash,
        'new_models': 'covariates-only and burden-only for six versions; PRS and joint reused',
        'splits': 'original outcome plus legacy carrier strata, 10x5, random_state=20260818'}, indent=2) + '\n')
    auc_rows, delta_rows, roc_rows = [], [], []
    repeat1 = {}
    reference_auc = pd.read_csv(paths[6]).set_index('version')
    with threadpool_limits(limits=1), ThreadPoolExecutor(max_workers=3) as pool:
        for repeat in range(1, 11):
            with np.load(MS_5_EXISTING / f'cv_repeat_{repeat:02d}.npz') as saved:
                assert str(saved['run_hash']) == historical_lock['run_hash']
                reused = saved['predictions']
                assert reused.shape == (len(existing_names), len(y))
            checkpoint = MS_5_CACHE / f'new_models_repeat_{repeat:02d}.npz'
            if checkpoint.exists():
                with np.load(checkpoint) as saved:
                    assert str(saved['run_hash']) == run_hash
                    new = saved['predictions']
            else:
                new = np.full((1 + len(MS_5_VERSIONS), len(y)), np.nan)
                for k in range(5):
                    tr, te = splits[(repeat - 1) * 5 + k]
                    tasks = [pool.submit(MS_5_fit_predict, y[tr], MS_5_design(cov[tr]), MS_5_design(cov[te]))]
                    for j in range(len(MS_5_VERSIONS)):
                        mean, sd = raw[tr, j].mean(), raw[tr, j].std(ddof=1)
                        assert sd > 0
                        tasks.append(pool.submit(MS_5_fit_predict, y[tr], MS_5_design(cov[tr], (raw[tr, j] - mean) / sd),
                                                 MS_5_design(cov[te], (raw[te, j] - mean) / sd)))
                    for j, task in enumerate(tasks):
                        new[j, te] = task.result()
                    print(f'Missing-model CV repeat {repeat}/10 fold {k + 1}/5 completed.', flush=True)
                assert np.isfinite(new).all() and ((new >= 0) & (new <= 1)).all()
                np.savez_compressed(checkpoint, run_hash=run_hash, predictions=new)
            assert new.shape == (1 + len(MS_5_VERSIONS), len(y)) and np.isfinite(new).all()
            if repeat == 1:
                assert np.allclose(new[0], old.oof_r1_covariates_only, atol=1e-8, rtol=1e-8)
                assert np.allclose(reused[0], old.oof_r1_prs_only, atol=1e-8, rtol=1e-8)
            for j, (version, label) in enumerate(MS_5_VERSIONS):
                predictions = [new[0], reused[0], new[j + 1], reused[existing_names.index(version)]]
                for model, pred in zip(MS_5_MODELS, predictions):
                    value = roc_auc_score(y, pred)
                    low, high = auc_ci(y, pred)
                    auc_rows.append({'version': version, 'model': model, 'repeat': repeat,
                                     'AUC': value, 'CI_low': low, 'CI_high': high})
                for a, b in [(3, 1), (3, 2), (3, 0), (1, 0), (2, 0)]:
                    delta_rows.append({'version': version, 'repeat': repeat,
                                       'comparison': MS_5_MODELS[a] + ' minus ' + MS_5_MODELS[b],
                                       **delong_test(y, predictions[a], predictions[b])})
                if repeat == 1:
                    repeat1[version] = predictions
                    assert np.isclose(roc_auc_score(y, predictions[3]), reference_auc.loc[version, 'AUC_repeat1'], atol=1e-12, rtol=1e-12)
                    fpr_grid = np.linspace(0, 1, 1001)
                    for model, pred in zip(MS_5_MODELS, predictions):
                        fpr, tpr, _ = roc_curve(y, pred)
                        interp = np.interp(fpr_grid, fpr, tpr)
                        for f, t in zip(fpr_grid, interp):
                            roc_rows.append({'version': version, 'model': model, 'repeat': 1, 'FPR': f, 'TPR_interpolated': t})
            print(f'Four-model AUC and paired DeLong results repeat {repeat}/10 completed.', flush=True)
    auc = pd.DataFrame(auc_rows)
    summary = auc.groupby(['version', 'model'], sort=False).AUC.agg(AUC_mean='mean', AUC_min='min', AUC_max='max').reset_index()
    summary = summary.merge(auc.loc[auc.repeat.eq(1)].drop(columns='repeat').rename(columns={
        'AUC': 'AUC_repeat1', 'CI_low': 'repeat1_CI_low', 'CI_high': 'repeat1_CI_high'}), on=['version', 'model'])
    for version, _ in MS_5_VERSIONS:
        joint = summary.loc[summary.version.eq(version) & summary.model.eq('Joint')].iloc[0]
        assert np.isclose(joint.AUC_mean, reference_auc.loc[version, 'AUC_mean'], rtol=1e-12, atol=1e-12)
    auc.to_csv(MS_5_DATA / 'four_models_AUC_by_repeat.csv', index=False)
    summary.to_csv(MS_5_DATA / 'four_models_AUC_summary.csv', index=False)
    pd.DataFrame(delta_rows).to_csv(MS_5_DATA / 'four_models_paired_delta_AUC_by_repeat.csv', index=False)
    pd.DataFrame(roc_rows).to_csv(MS_5_DATA / 'four_models_ROC_repeat1_interpolated.csv', index=False)
    plt.style.use('default')
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), sharex=True, sharey=True)
    for j, ((version, label), ax) in enumerate(zip(MS_5_VERSIONS, axes.flat)):
        for model, pred, color in zip(MS_5_MODELS, repeat1[version], MS_5_COLORS):
            fpr, tpr, _ = roc_curve(y, pred)
            ax.plot(fpr, tpr, color=color, linewidth=1.6)
        ax.plot([0, 1], [0, 1], color='.7', linestyle='--', linewidth=.8)
        ax.set_title(f'({chr(97 + j)}) {label}', loc='left', fontsize=11)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.grid(False)
        ax.spines[['top', 'right']].set_visible(False)
        if j >= 3: ax.set_xlabel('False positive rate')
        if j % 3 == 0: ax.set_ylabel('True positive rate')
    fig.legend([Line2D([0], [0], color=c, linewidth=2) for c in MS_5_COLORS], MS_5_MODELS,
               loc='upper center', bbox_to_anchor=(.5, .955), ncol=4, frameon=False)
    fig.suptitle('Discrimination of four case-status models: ROC curves', fontsize=15, y=.995)
    fig.text(.5, .01, 'First-repeat out-of-fold predictions; same frozen cohort and paired folds across all six versions.', ha='center', fontsize=10)
    fig.tight_layout(rect=[0, .04, 1, .92])
    roc_preview = MS_5_save_figure(fig, 'four_case_status_models_ROC_sensitivities')
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), sharey=True)
    ymin = .5
    ymax = max(summary.AUC_mean.max(), summary.repeat1_CI_high.max()) + .012
    for j, ((version, label), ax) in enumerate(zip(MS_5_VERSIONS, axes.flat)):
        rows = summary.loc[summary.version.eq(version)].set_index('model').loc[MS_5_MODELS]
        xx = np.arange(4)
        ax.bar(xx, rows.AUC_mean - ymin, bottom=ymin, color=MS_5_COLORS, width=.58, alpha=.85)
        # Exact CI endpoints are centered on repeat1 AUC, not re-centered on mean.
        ax.vlines(xx, rows.repeat1_CI_low, rows.repeat1_CI_high, color='black', linewidth=1.2)
        ax.hlines(rows.repeat1_CI_low, xx - .075, xx + .075, color='black', linewidth=1.2)
        ax.hlines(rows.repeat1_CI_high, xx - .075, xx + .075, color='black', linewidth=1.2)
        for x, mean, high in zip(xx, rows.AUC_mean, rows.repeat1_CI_high):
            ax.text(x, max(mean, high) + .002, f'{mean:.4f}', ha='center', va='bottom', fontsize=9)
        ax.set_xticks(xx, ['Covariates', 'PRS + cov.', 'Burden + cov.', 'Joint'], rotation=25, ha='right', fontsize=9)
        ax.set_ylim(ymin, ymax); ax.set_title(f'({chr(97 + j)}) {label}', loc='left', fontsize=11)
        ax.grid(False); ax.spines[['top', 'right']].set_visible(False)
        if j % 3 == 0: ax.set_ylabel('Mean out-of-fold AUC')
    fig.suptitle('Discrimination of four case-status models: AUC', fontsize=15)
    fig.text(.5, .015, 'Bars: mean AUC over 10 repeats of 5-fold CV. Error bars: first-repeat DeLong 95% CIs.\n'
             'Y-axis starts at 0.50. First-repeat intervals describe that repeat, not uncertainty in the ten-repeat mean.',
             ha='center', fontsize=10)
    fig.tight_layout(rect=[0, .065, 1, .955])
    auc_preview = MS_5_save_figure(fig, 'four_case_status_models_AUC_sensitivities')
    for p, digest in hashes.items():
        assert MS_5_sha(MS_5_ROOT / p) == digest, 'Input changed during fitting: ' + p
    validation = {'status': 'four_model_discrimination_completed', 'versions': [v for v, _ in MS_5_VERSIONS],
        'models': MS_5_MODELS, 'n': len(y), 'case_n': int(y.sum()), 'CV_repeats': 10, 'CV_folds': 5,
        'PRS_and_joint_reused': True, 'new_models_fitted': 'Covariates only and burden + covariates',
        'all_fits_converged': True, 'original_first_repeat_folds_reproduced': True,
        'covariates_first_repeat_OOF_reproduced': True, 'all_joint_AUCs_match_previous_results': True,
        'score_scaling': 'training fold mean and SD ddof=1; no validation fold fitting',
        'joint_definition': 'PRS + burden + PRS*burden + age_at_last_obs + PC1-10 + sequencing site and sequencer',
        'ROC_repeat': 1, 'AUC_bar_statistic': '10-repeat mean', 'AUC_CI': 'Exact first-repeat DeLong CI endpoints',
        'CI_limit': 'conditional on fitted OOF predictions; excludes full retraining uncertainty',
        'final_variant_definition_frozen': False, 'ROC_preview': roc_preview, 'AUC_preview': auc_preview,
        'inputs_unchanged': True, 'participant_predictions_internal_only': True}
    (MS_5_ROOT / 'results/report/51_four_model_discrimination_validation.json').write_text(json.dumps(validation, indent=2) + '\n')
    (MS_5_FIG / 'four_case_status_models_captions.md').write_text(
        '# Four case-status model discrimination\n\n'
        'Each figure contains six panels: main C; without STK11; PTV-only; strict TP53/NF1 VAF; ClinVar ≥1 star; without CDH1.\n\n'
        '(a) ROC curves of first-repeat out-of-fold predictions from the four case-status models. '
        '(b) Mean out-of-fold AUC over 10 repeats of 5-fold cross-validation. Error bars are DeLong 95% CIs for the first repeat. '
        'These intervals are not confidence intervals for the ten-repeat mean. All models include the same clinical/technical covariates; '
        'the joint model also includes PRS×burden. Same frozen cohort (N=141,204; cases=8,528) and original paired folds. '
        'Existing PRS and joint predictions reused; missing covariates-only and burden-only predictions fitted. '
        'Current variant definition remains provisional.\n')
    print('ROC and AUC figures, unrounded source CSVs and validation saved.', flush=True)

# ============================================================================
# 6. joint-distribution: plot_joint_distribution
# ============================================================================
"""Update PRS-decile/carrier heatmaps using main C binary carrier status."""
from pathlib import Path
import hashlib
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd
import statsmodels.api as sm
from threadpoolctl import threadpool_limits

MS_6_ROOT = Path(__file__).resolve().parents[1]
MS_6_FIG = MS_6_ROOT / 'manuscript/figures'
MS_6_DATA = MS_6_ROOT / 'manuscript/datas/case_control'
MS_6_INTERNAL = MS_6_ROOT / 'data/manuscript_refresh_20261006/joint_distribution'


def plot_joint_distribution():
    for p in [MS_6_FIG, MS_6_DATA, MS_6_INTERNAL]:
        p.mkdir(parents=True, exist_ok=True)
    paths = [MS_6_ROOT / 'data/29_rowlands_extension/section9/frozen_model_inputs.csv.gz',
             MS_6_ROOT / 'manuscript/datas/burden_13gene_case_control.csv',
             MS_6_ROOT / 'manuscript/datas/PRS_scores.csv', Path(__file__)]
    hashes = {str(p.relative_to(MS_6_ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    d = pd.read_csv(paths[0]).set_index('person_id')
    burden = pd.read_csv(paths[1]).set_index('sample_id').loc[d.index]
    prs = pd.read_csv(paths[2], usecols=['person_id', 'PRS311']).set_index('person_id').loc[d.index]
    assert d.index.is_unique and burden.index.is_unique and prs.index.is_unique
    assert len(d) == 141204 and d.is_case.sum() == 8528
    assert np.allclose(d.PRS311, prs.PRS311, atol=1e-12, rtol=1e-12)
    carrier = burden.B_raw_13gene_main.gt(0).astype(int)
    assert carrier.equals(burden.carried_gene_count_13gene_main.gt(0).astype(int))
    decile = pd.qcut(d.PRS311, 10, labels=False) + 1
    pc = [f'PC{i}' for i in range(1, 11)]
    X = pd.DataFrame({'const': 1., 'PRS311_z': (d.PRS311 - d.PRS311.mean()) / d.PRS311.std(ddof=1),
                      'C_carrier': carrier}, index=d.index).join(d[['age_at_last_obs'] + pc])
    assert np.isfinite(X.to_numpy()).all() and not any('interaction' in c for c in X)
    with threadpool_limits(limits=1):
        fit = sm.GLM(d.is_case, X, family=sm.families.Binomial()).fit(cov_type='HC1')
    assert fit.converged and np.isfinite(fit.params).all()
    predicted = fit.predict(X)
    assert np.isclose(predicted.sum(), d.is_case.sum(), atol=1e-5, rtol=0)
    coef = pd.DataFrame({'term': X.columns, 'beta': fit.params.to_numpy(), 'SE_HC1': fit.bse.to_numpy(),
                          'OR': np.exp(fit.params.to_numpy()), 'CI_low': np.exp(fit.conf_int()[0].to_numpy()),
                          'CI_high': np.exp(fit.conf_int()[1].to_numpy()), 'p_value': fit.pvalues.to_numpy()})
    coef.to_csv(MS_6_DATA / 'joint_distribution_no_interaction_coefficients.csv', index=False)
    rows = []
    for k in range(1, 11):
        for group in [0, 1]:
            mask = decile.eq(k) & carrier.eq(group)
            n = int(mask.sum()); cases = int(d.loc[mask].is_case.sum()); controls = n - cases
            expected = float(predicted.loc[mask].sum())
            assert n > 0 and expected > 0
            rows.append({'PRS_decile': k, 'C_carrier': group, 'participants': n, 'cases': cases,
                         'controls': controls, 'observed_case_proportion': cases / n,
                         'expected_cases_no_interaction': expected, 'mean_expected_probability': expected / n,
                         'observed_over_expected': cases / expected,
                         'disclosure_safe': min(n, cases, controls) >= 20})
    cells = pd.DataFrame(rows)
    assert cells.participants.sum() == len(d) and cells.cases.sum() == d.is_case.sum()
    assert np.isclose(cells.expected_cases_no_interaction.sum(), d.is_case.sum(), atol=1e-5, rtol=0)
    for group in [0, 1]:
        sub = cells.loc[cells.C_carrier.eq(group)]
        assert np.isclose(sub.cases.sum(), sub.expected_cases_no_interaction.sum(), atol=1e-5, rtol=0)
    cells.to_csv(MS_6_INTERNAL / 'cells_exact_internal.csv', index=False)
    # Suppress the complementary carrier cell too if any outcome cell is small.
    hidden_deciles = set(cells.loc[~cells.disclosure_safe, 'PRS_decile'])
    public = cells.drop(columns='disclosure_safe').copy()
    public['cell_suppressed'] = public.PRS_decile.isin(hidden_deciles)
    columns = ['participants', 'cases', 'controls', 'observed_case_proportion',
               'expected_cases_no_interaction', 'mean_expected_probability', 'observed_over_expected']
    for c in columns:
        public[c] = public[c].astype(object)
        public.loc[public.cell_suppressed, c] = 'suppressed'
    public.to_csv(MS_6_DATA / 'PRS311_decile_C_carrier_joint_distribution.csv', index=False)
    observed = cells.pivot(index='PRS_decile', columns='C_carrier', values='observed_case_proportion').to_numpy() * 100
    ratio = cells.pivot(index='PRS_decile', columns='C_carrier', values='observed_over_expected').to_numpy()
    for k in hidden_deciles:
        observed[k - 1, :] = np.nan; ratio[k - 1, :] = np.nan
    plt.style.use('default')
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 7.5), sharey=True)
    image_a = axes[0].imshow(observed, origin='lower', aspect='auto', cmap='Reds', vmin=0, vmax=np.nanmax(observed))
    low, high = np.nanmin(ratio), np.nanmax(ratio)
    span = max(abs(low - 1), abs(high - 1), .05)
    image_b = axes[1].imshow(ratio, origin='lower', aspect='auto', cmap='RdBu_r', norm=TwoSlopeNorm(vmin=1-span, vcenter=1, vmax=1+span))
    fig.colorbar(image_a, ax=axes[0], label='Observed case proportion (%)', fraction=.05, pad=.04)
    fig.colorbar(image_b, ax=axes[1], label='Observed / expected cases', fraction=.05, pad=.04)
    maximum = np.nanmax(observed)
    for row in cells.itertuples():
        y, x = row.PRS_decile - 1, row.C_carrier
        if row.PRS_decile in hidden_deciles:
            for ax in axes: ax.text(x, y, 'Suppressed', ha='center', va='center', fontsize=9)
            continue
        color = 'white' if observed[y, x] > maximum * .55 else 'black'
        axes[0].text(x, y, f'{observed[y,x]:.2f}%\nn={row.participants:,}; cases={row.cases:,}',
                     ha='center', va='center', fontsize=8.5, color=color)
        color = 'white' if abs(ratio[y, x] - 1) > span * .65 else 'black'
        axes[1].text(x, y, f'{ratio[y,x]:.2f}', ha='center', va='center', fontsize=10, color=color)
    for ax in axes:
        ax.set_xticks([0, 1], ['Non-carrier', 'Carrier'])
        ax.set_yticks(np.arange(10), np.arange(1, 11))
        ax.set_xlabel('13-gene PTV + eligible missense carrier status')
        ax.grid(False)
        for spine in ax.spines.values(): spine.set_visible(False)
    axes[0].set_ylabel('PRS311 decile (full cohort)')
    axes[0].set_title('(a) Observed case proportion', loc='left', fontsize=12)
    axes[1].set_title('(b) Observed / expected cases', loc='left', fontsize=12)
    fig.suptitle('Joint distribution of PRS311 and 13-gene carrier status\nN = 141,204', fontsize=14)
    fig.text(.5, .02, 'Left: cohort-observed proportions, not absolute risks. Right: expected cases from a logistic model\n'
             'with PRS311 + binary carrier + age at last observation + PC1–10, without interaction.', ha='center', fontsize=9)
    fig.tight_layout(rect=[0, .075, 1, .94])
    stem = MS_6_FIG / 'PRS311_decile_13gene_carrier_joint_distribution'
    fig.savefig(stem.with_suffix('.svg'), bbox_inches='tight')
    fig.savefig(stem.with_suffix('.png'), dpi=180, bbox_inches='tight')
    fig.savefig(stem.with_suffix('.jpg'), dpi=180, bbox_inches='tight', pil_kwargs={'quality':85, 'subsampling':2})
    preview = '.png' if stem.with_suffix('.png').stat().st_size <= stem.with_suffix('.jpg').stat().st_size else '.jpg'
    stem.with_suffix('.jpg' if preview == '.png' else '.png').unlink()
    plt.close(fig)
    for p, digest in hashes.items(): assert hashlib.sha256((MS_6_ROOT / p).read_bytes()).hexdigest() == digest
    result = {'status':'C_joint_distribution_heatmap_completed', 'n':len(d), 'case_n':int(d.is_case.sum()),
              'carrier_definition':'At least one main-C qualifying PTV or P/LP missense in 13 genes',
              'PRS_deciles':'pd.qcut of raw PRS311 in full frozen case-control cohort',
              'expected_model':'unweighted GLM Binomial HC1: PRS311_z + C_carrier + age_at_last_obs + PC1-10; no interaction',
              'expected_model_matches_legacy_adjustment':True, 'expected_predictions':'full-cohort fitted, not OOF',
              'model_converged':True, 'observed_expected_totals_and_carrier_margins_agree':True,
              'raw_scores_match_frozen_input_PRS':True, 'small_outcome_cells_and_complements_suppressed':True,
              'final_variant_definition_frozen':False, 'preview':str(stem.with_suffix(preview).relative_to(MS_6_ROOT))}
    (MS_6_ROOT / 'results/report/52_C_joint_distribution_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    (MS_6_INTERNAL / 'input_and_model_manifest.json').write_text(json.dumps({'inputs':hashes, **result},indent=2)+'\n')
    (MS_6_FIG / 'PRS311_decile_13gene_carrier_joint_distribution_caption.md').write_text(
        'Joint distribution of PRS311 and 13-gene PTV plus eligible P/LP missense carrier status (N=141,204). '
        'Left: observed case proportion by full-cohort PRS311 decile and binary carrier status. These proportions depend on '
        'the case-control composition of the cohort and are not absolute risks. Right: observed/expected cases; expected counts '
        'are sums of fitted probabilities from a logistic model without interaction, adjusted for continuous PRS311, binary carrier status, '
        'age at last EHR observation and PC1–10. The adjustment matches the prior figure; sequencing site and sequencer are not included '
        'in this descriptive expectation model. O/E differences do not independently establish biological interaction. Current variant definition is provisional.\n')
    print('Updated C-carrier joint-distribution heatmap and unrounded aggregate CSVs saved; model and cell totals validated.')

# ============================================================================
# 7. score-deciles: plot_score_decile_boxplots
# ============================================================================
"""Internal score boxplots by full-cohort PRS311 decile."""
from pathlib import Path
import json
import hashlib

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

MS_7_ROOT = Path(__file__).resolve().parents[1]
MS_7_FIG = MS_7_ROOT / 'manuscript/figures'
MS_7_INTERNAL = MS_7_ROOT / 'data/manuscript_refresh_20261006/score_decile_boxplots'


def plot_score_decile_boxplots():
    MS_7_FIG.mkdir(exist_ok=True); MS_7_INTERNAL.mkdir(parents=True, exist_ok=True)
    paths = [MS_7_ROOT / 'data/29_rowlands_extension/section9/frozen_model_inputs.csv.gz',
             MS_7_ROOT / 'manuscript/datas/PRS_scores.csv', MS_7_ROOT / 'manuscript/datas/burden_13gene_case_control.csv']
    d = pd.read_csv(paths[0]).set_index('person_id')
    p = pd.read_csv(paths[1], usecols=['person_id', 'PRS311']).set_index('person_id').loc[d.index]
    b = pd.read_csv(paths[2], usecols=['sample_id', 'B_raw_13gene_main']).set_index('sample_id').loc[d.index]
    assert len(d) == 141204 and int(d.is_case.sum()) == 8528
    assert np.allclose(d.PRS311, p.PRS311, atol=1e-12, rtol=1e-12)
    decile = pd.qcut(d.PRS311, 10, labels=False) + 1
    labels = [f'{10 * i}–{10 * (i + 1)}%' for i in range(10)]
    labels[0] += '\nLowest'; labels[-1] += '\nHighest'
    colors = [(1, 'Case', '#2463B4', -.18), (0, 'Control', '#D54848', .18)]
    plt.style.use('default')
    fig, axes = plt.subplots(2, 1, figsize=(12, 7.5), sharex=True)
    specs = [('PRS311', p.PRS311, False),
             ('13-gene burden: non-zero scores only', b.B_raw_13gene_main, True)]
    summary = []
    for ax, (name, values, positive_only) in zip(axes, specs):
        for case, group, color, offset in colors:
            samples = []
            for k in range(1, 11):
                mask = decile.eq(k) & d.is_case.eq(case)
                if positive_only: mask &= values.gt(0)
                v = values.loc[mask].to_numpy()
                assert len(v) >= 20
                samples.append(v)
                summary.append({'panel': name, 'PRS_decile': k, 'group': group, 'n':len(v),
                                'q25':float(np.quantile(v,.25)), 'median':float(np.median(v)), 'q75':float(np.quantile(v,.75))})
            boxes = ax.boxplot(samples, positions=np.arange(1,11)+offset, widths=.3,
                               patch_artist=True, showfliers=False, manage_ticks=False,
                               medianprops={'color':'black','linewidth':1.2},
                               boxprops={'edgecolor':color,'linewidth':1.1},
                               whiskerprops={'color':color}, capprops={'color':color})
            for patch in boxes['boxes']: patch.set_facecolor(color);patch.set_alpha(.65)
        ax.set_title(name, loc='left', fontsize=12)
        ax.set_ylabel('Raw PRS311 score' if name=='PRS311' else 'Raw weighted burden score')
        ax.grid(False);ax.spines[['top','right']].set_visible(False)
        if name!='PRS311': ax.set_ylim(bottom=-.12)
    axes[0].legend(handles=[Patch(facecolor=c, alpha=.65, label=g) for _,g,c,_ in colors], frameon=False, loc='upper left')
    axes[-1].set_xticks(np.arange(1,11),labels)
    axes[-1].set_xlabel('PRS311 percentile group in the full case-control cohort (deciles 1–10)')
    axes[-1].set_xlim(.4,10.6)
    fig.suptitle('Raw PRS and burden scores across PRS311 deciles', fontsize=15)
    fig.text(.5,.01,'Blue: case; red: control. Boxes: median and IQR; whiskers: within 1.5×IQR. Outliers are not shown.\n'
             'Percentile labels denote intervals, not cumulative groups. Bottom panel includes carriers only.',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.055,1,.96])
    stem=MS_7_FIG/'PRS311_decile_PRS_burden_boxplots'
    fig.savefig(stem.with_suffix('.svg'),bbox_inches='tight')
    fig.savefig(stem.with_suffix('.png'),dpi=180,bbox_inches='tight')
    fig.savefig(stem.with_suffix('.jpg'),dpi=180,bbox_inches='tight',pil_kwargs={'quality':85,'subsampling':2})
    ext='.png' if stem.with_suffix('.png').stat().st_size<=stem.with_suffix('.jpg').stat().st_size else '.jpg'
    stem.with_suffix('.jpg' if ext=='.png' else '.png').unlink();plt.close(fig)
    pd.DataFrame(summary).to_csv(MS_7_INTERNAL/'boxplot_statistics_internal.csv',index=False)
    manifest={'cohort_n':len(d),'deciles':'PRS311 qcut in full frozen cohort','scores':'raw PRS311 and B_raw_13gene_main',
              'panels':'PRS; nonzero burden','case_control_colors':'blue/red',
              'outliers_displayed':False,'internal_review_only':True,
              'preview':str(stem.with_suffix(ext).relative_to(MS_7_ROOT)),
              'inputs':{str(q.relative_to(MS_7_ROOT)):hashlib.sha256(q.read_bytes()).hexdigest() for q in paths}}
    (MS_7_INTERNAL/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Raw score decile boxplots saved; original case-proportion/OE heatmap retained.')

# ============================================================================
# 8. burden-histogram: plot_burden_histogram
# ============================================================================
"""Within-group percentage histogram of all-participant raw C burden."""
from pathlib import Path
import hashlib
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

MS_8_ROOT = Path(__file__).resolve().parents[1]
MS_8_FIG = MS_8_ROOT / 'manuscript/figures'
MS_8_DATA = MS_8_ROOT / 'manuscript/datas/burden'


def plot_burden_histogram():
    MS_8_FIG.mkdir(exist_ok=True); MS_8_DATA.mkdir(parents=True, exist_ok=True)
    inputs = [MS_8_ROOT / 'data/29_rowlands_extension/section9/frozen_model_inputs.csv.gz',
              MS_8_ROOT / 'manuscript/datas/burden_13gene_case_control.csv']
    d = pd.read_csv(inputs[0]).set_index('person_id')
    b = pd.read_csv(inputs[1], usecols=['sample_id', 'B_raw_13gene_main']).set_index('sample_id')
    assert d.index.is_unique and b.index.is_unique and set(d.index) == set(b.index)
    score = b.loc[d.index, 'B_raw_13gene_main']
    assert len(d) == 141204 and int(d.is_case.sum()) == 8528
    assert np.isfinite(score).all() and score.ge(0).all()
    # Prespecified broad intervals. Merge sparse upper bins without revealing
    # the data-derived maximum or exact small counts in the exported table.
    edges = [0., .5, 1., 1.5, 2., 2.5, 3., float('inf')]
    groups = [(1, 'Case', '#2463B4'), (0, 'Control', '#D54848')]
    def counts(e):
        return [np.histogram(score.loc[d.is_case.eq(case) & score.gt(0)], bins=e)[0]
                for case, _, _ in groups]
    while True:
        cnt = counts(edges)
        unsafe = [i for i in range(len(edges)-1) if any(0 < int(c[i]) < 20 for c in cnt)]
        if not unsafe:
            break
        i = unsafe[-1]
        if i > 0:
            del edges[i]
        elif len(edges) > 2:
            del edges[1]
        else:
            raise AssertionError('Insufficient nonzero participants for grouped release')
    rows = []
    for case, group, _ in groups:
        v = score.loc[d.is_case.eq(case)]; n = len(v); zero = int(v.eq(0).sum())
        assert zero >= 20 and n-zero >= 20
        rows.append({'group':group, 'bin':'Exactly 0', 'participants':zero, 'denominator_all_group_participants':n,
                     'percentage':100*zero/n})
        hist = np.histogram(v.loc[v.gt(0)], bins=edges)[0]
        for i, count in enumerate(hist):
            label = f'[{edges[i]:g}, {edges[i+1]:g})' if np.isfinite(edges[i+1]) else f'>= {edges[i]:g}'
            if i == 0:
                label = f'(0, {edges[i+1]:g})' if np.isfinite(edges[i+1]) else '> 0'
            rows.append({'group':group, 'bin':label, 'participants':int(count),
                         'denominator_all_group_participants':n, 'percentage':100*count/n})
        assert zero + int(hist.sum()) == n
    table = pd.DataFrame(rows)
    for group in ['Case','Control']:
        assert np.isclose(table.loc[table.group.eq(group), 'percentage'].sum(), 100)
    table.to_csv(MS_8_DATA/'13gene_burden_histogram_percentages.csv', index=False)
    plt.style.use('default')
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.9), gridspec_kw={'width_ratios':[1, 2.8]})
    for j, (_, group, color) in enumerate(groups):
        r = table.loc[table.group.eq(group)].iloc[0]
        axes[0].bar(j, r.percentage, color=color, alpha=.75, width=.55)
        axes[0].text(j, r.percentage+1, f'{r.percentage:.2f}%', ha='center', fontsize=11)
    axes[0].set_xticks([0,1],['Case','Control'])
    axes[0].set_ylim(0,105); axes[0].set_ylabel('Participants within group (%)')
    axes[0].set_title('(a) Zero burden',loc='left')
    ticklabels = []
    for i in range(len(edges)-1):
        ticklabels.append(f'>{edges[i]:g}–<{edges[i+1]:g}' if i==0 else
                          (f'{edges[i]:g}–<{edges[i+1]:g}' if np.isfinite(edges[i+1]) else f'≥{edges[i]:g}'))
    xx = np.arange(len(edges)-1)
    ymax = 0.
    for j, (_,group,color) in enumerate(groups):
        sub = table.loc[table.group.eq(group)].iloc[1:]
        offset = -.17 if j==0 else .17
        axes[1].bar(xx+offset,sub.percentage,width=.32,color=color,alpha=.75,label=group)
        for x,val in zip(xx+offset,sub.percentage):
            axes[1].text(x,val+.025,f'{val:.2f}%',ha='center',va='bottom',fontsize=9)
        ymax = max(ymax,float(sub.percentage.max()))
    axes[1].set_xticks(xx,ticklabels)
    axes[1].set_ylim(0,ymax*1.22)
    axes[1].set_xlabel('Raw weighted burden score interval')
    axes[1].set_ylabel('Participants within group (%)')
    axes[1].set_title('(b) Non-zero score intervals',loc='left')
    axes[1].legend(frameon=False)
    for ax in axes:
        ax.grid(False);ax.spines[['top','right']].set_visible(False)
    fig.suptitle('13-gene PTV + eligible missense burden: participant proportions',fontsize=14)
    fig.text(.5,.01,'Both panels use all cases or all controls as their denominator; non-zero bars are not renormalized among carriers.\n'
             'Zero and positive-score intervals together sum to 100% within each group. Panel y-axis scales differ.',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.08,1,.94])
    stem = MS_8_FIG/'13gene_burden_histogram_percentages'
    fig.savefig(stem.with_suffix('.svg'),bbox_inches='tight')
    fig.savefig(stem.with_suffix('.png'),dpi=180,bbox_inches='tight')
    fig.savefig(stem.with_suffix('.jpg'),dpi=180,bbox_inches='tight',pil_kwargs={'quality':85,'subsampling':2})
    preview = '.png' if stem.with_suffix('.png').stat().st_size<=stem.with_suffix('.jpg').stat().st_size else '.jpg'
    stem.with_suffix('.jpg' if preview=='.png' else '.png').unlink();plt.close(fig)
    report={'status':'all_participant_C_burden_histogram_completed','denominator':'All cases or all controls separately',
            'cohort_n':len(d),'positive_bins_merged_if_nonzero_counts_below20':True,
            'percentage_totals_verified':True,'zero_CI_or_density_not_used':True,
            'preview':str(stem.with_suffix(preview).relative_to(MS_8_ROOT))}
    (MS_8_ROOT/'results/report/53_burden_histogram_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    manifest_path = MS_8_ROOT/'manuscript/datas/manifest.json'
    manifest=json.loads(manifest_path.read_text());name='burden/13gene_burden_histogram_percentages.csv'
    manifest['files']=[r for r in manifest['files'] if r['file']!=name]
    manifest['files'].append({'file':name,'source':'scripts/09_result_figures.py --task burden-histogram','purpose':'Within-case/control percentages; exact-zero separate; positive bins merged for disclosure',
                              'sha256':hashlib.sha256((MS_8_ROOT/'manuscript/datas'/name).read_bytes()).hexdigest()})
    manifest_path.write_text(json.dumps(manifest,indent=2)+'\n')
    print('Burden percentage histogram and aggregate source CSV saved; percentages sum to 100% within each group.')

# ============================================================================
# 9. age-boxplots: plot_age_burden_boxplots
# ============================================================================
"""Recorded diagnosis age by current 13-gene carrier burden group."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

MS_9_ROOT = Path(__file__).resolve().parents[1]
MS_9_BASE = MS_9_ROOT / 'data/29_rowlands_extension'


def plot_age_burden_boxplots():
    ids = pd.read_csv(MS_9_BASE / 'diagnosis_age_ids.csv.gz').person_id
    d = pd.read_csv(MS_9_BASE / 'section9/frozen_model_inputs.csv.gz').set_index('person_id').loc[ids]
    s = pd.read_csv(MS_9_BASE / 'section8/comparable_scores_case_control.csv.gz').set_index('sample_id').loc[ids]
    b = s.B_raw_C
    age = d.age_at_diagnosis
    assert len(d) == 7996 and d.index.is_unique and d.is_case.eq(1).all()
    assert age.notna().all() and np.isfinite(age).all() and np.isfinite(b).all()
    assert b.ge(0).all() and b.gt(0).equals(s.carried_gene_count_C.gt(0))
    figdir = MS_9_ROOT / 'manuscript/figures'
    figdir.mkdir(parents=True, exist_ok=True)
    rows = []
    manifest = {'n': len(d), 'score': 'B_raw_C: 13-gene PTV + qualifying P/LP missense, sum of ln(OR) per carried gene',
                'outliers_shown': False, 'variants': {}}
    for name, cutoff in [('carrier_median', float(b[b.gt(0)].median()))]:
        masks = [b.gt(0) & b.le(cutoff), b.gt(cutoff), b.eq(0)]
        groups = ['Carrier, lower burden', 'Carrier, higher burden', 'Non-carrier']
        samples = [age.loc[m].to_numpy() for m in masks]
        assert sum(map(len, samples)) == 7996 and min(map(len, samples)) >= 20
        labels = [f'Carrier\nburden ≤ median ({cutoff:.2f})\n(n={len(samples[0]):,})',
                  f'Carrier\nburden > median ({cutoff:.2f})\n(n={len(samples[1]):,})',
                  f'Non-carrier\nburden = 0\n(n={len(samples[2]):,})']
        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.boxplot(samples, tick_labels=labels, showfliers=False, widths=.5,
                   medianprops={'color': 'black'})
        ax.set_ylabel('Age at diagnosis (years)')
        ax.set_title('(a) Age at diagnosis by burden group\n13-gene age-at-diagnosis cohort (n=7,996)', fontsize=12)
        ax.grid(False)
        ax.spines[['top', 'right']].set_visible(False)
        fig.text(.5, .015, 'Boxes: median and IQR; whiskers: within 1.5×IQR. Outliers are not shown.',
                 ha='center', fontsize=8)
        fig.tight_layout(rect=[0, .045, 1, 1])
        stem = figdir / f'age_at_diagnosis_13gene_burden_{name}'
        sizes = {}
        for ext in ['svg', 'png', 'jpg']:
            kwargs = {'pil_kwargs': {'quality': 85, 'subsampling': 2}} if ext == 'jpg' else {}
            fig.savefig(stem.with_suffix('.' + ext), dpi=180, bbox_inches='tight', **kwargs)
            sizes[ext] = stem.with_suffix('.' + ext).stat().st_size
        # Keep the smallest format plus a PNG review preview when SVG is smallest.
        smallest = min(sizes, key=sizes.get)
        for ext in ['svg', 'jpg']:
            if ext != smallest:
                stem.with_suffix('.' + ext).unlink()
        plt.close(fig)
        manifest['variants'][name] = {'cutoff_raw': cutoff, 'group_n': list(map(len, samples)),
                                     'smallest_format': smallest, 'file_sizes': sizes}
        for group, v in zip(groups, samples):
            rows.append({'split': name, 'cutoff_raw': cutoff, 'group': group, 'n': len(v),
                         'q25_age': np.quantile(v, .25), 'median_age': np.median(v), 'q75_age': np.quantile(v, .75)})
    pd.DataFrame(rows).to_csv(MS_9_ROOT / 'results/table/54_age_13gene_boxplot_summary.csv', index=False)
    (MS_9_ROOT / 'results/report/54_age_13gene_boxplot_validation.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(manifest, indent=2))

# ============================================================================
# 10. age-carrier-frequency: plot_age_carrier_frequency
# ============================================================================
"""Reuse the unchanged PRS age panel and plot current 13-gene carrier rates."""
from pathlib import Path
import hashlib
import json
import shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from statsmodels.stats.proportion import proportion_confint

MS_10_ROOT = Path(__file__).resolve().parents[1]
MS_10_BASE = MS_10_ROOT / 'data/29_rowlands_extension'
MS_10_FIG = MS_10_ROOT / 'manuscript/figures'


def plot_age_carrier_frequency():
    MS_10_FIG.mkdir(parents=True, exist_ok=True)
    source = MS_10_ROOT / 'results/figures/17_prs311_vs_age_at_diagnosis_hexbin.jpg'
    target = MS_10_FIG / source.name
    shutil.copy2(source, target)
    digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    assert digest(source) == digest(target)
    ids = pd.read_csv(MS_10_BASE / 'diagnosis_age_ids.csv.gz').person_id
    d = pd.read_csv(MS_10_BASE / 'section9/frozen_model_inputs.csv.gz').set_index('person_id').loc[ids]
    s = pd.read_csv(MS_10_BASE / 'section8/comparable_scores_case_control.csv.gz').set_index('sample_id').loc[ids]
    assert len(d) == 7996 and d.index.is_unique and d.is_case.eq(1).all()
    assert d.age_at_diagnosis.notna().all()
    carrier = s.carried_gene_count_C.gt(0)
    assert carrier.equals(s.B_raw_C.gt(0))
    labels = ['<35', '35–<45', '45–<50', '50–<55', '55–<60', '60–<65', '65–<70', '≥70']
    groups = pd.cut(d.age_at_diagnosis, [-np.inf, 35, 45, 50, 55, 60, 65, 70, np.inf],
                    labels=labels, right=False)
    summary = pd.DataFrame({'age_group': groups, 'carrier': carrier}).groupby('age_group', observed=False).carrier.agg(n_cases='size', n_carriers='sum')
    assert summary.n_cases.sum() == 7996 and summary.n_carriers.sum() == carrier.sum()
    assert (summary.n_carriers >= 20).all() and (summary.n_cases - summary.n_carriers >= 20).all()
    lo, hi = proportion_confint(summary.n_carriers.to_numpy(), summary.n_cases.to_numpy(), alpha=.05, method='wilson')
    summary['carrier_rate_pct'] = summary.n_carriers / summary.n_cases * 100
    summary['wilson95_low_pct'] = lo * 100
    summary['wilson95_high_pct'] = hi * 100
    plt.style.use('default')
    fig, ax = plt.subplots(figsize=(9, 5))
    x = np.arange(len(labels))
    rate = summary.carrier_rate_pct.to_numpy()
    ax.bar(x, rate, width=.65, color='white', edgecolor='black',
           yerr=np.vstack([rate-lo*100, hi*100-rate]), capsize=4)
    ax.set_xticks(x, [f'{label}\n(n={n:,})' for label, n in zip(labels, summary.n_cases)], fontsize=9)
    ax.set_xlabel('Age at diagnosis (years)')
    ax.set_ylabel('Carrier rate (%)')
    ax.set_title('(c) 13-gene carrier frequency by age at diagnosis\nAge-at-diagnosis cohort (n=7,996)', fontsize=12)
    ax.set_ylim(bottom=0)
    ax.grid(False)
    ax.spines[['top', 'right']].set_visible(False)
    fig.text(.5, .015, 'Bars: carrier proportion within each age group. Error bars: 95% Wilson intervals.', ha='center', fontsize=9)
    fig.tight_layout(rect=[0, .045, 1, 1])
    stem = MS_10_FIG / '13gene_carrier_frequency_by_age_at_diagnosis'
    sizes = {}
    for ext in ['svg', 'png', 'jpg']:
        kwargs = {'pil_kwargs': {'quality': 85, 'subsampling': 2}} if ext == 'jpg' else {}
        fig.savefig(stem.with_suffix('.' + ext), dpi=180, bbox_inches='tight', **kwargs)
        sizes[ext] = stem.with_suffix('.' + ext).stat().st_size
    smallest = min(sizes, key=sizes.get)
    for ext in ['svg', 'jpg']:
        if ext != smallest:
            stem.with_suffix('.' + ext).unlink()
    plt.close(fig)
    summary.to_csv(MS_10_ROOT / 'results/table/55_13gene_carrier_frequency_by_diagnosis_age.csv')
    caption = '''Rare and common risk and age at diagnosis among cases (age-at-diagnosis cohort, n=7,996).

(a) Age at diagnosis by 13-gene burden group. Carriers are divided at the carrier burden median (1.73695, displayed as 1.74); non-carriers have zero burden. Boxes show medians and interquartile ranges; whiskers extend to observations within 1.5 times the interquartile range. Outliers are not shown.

(b) PRS311 against age at diagnosis. Hexagonal bins with fewer than 20 cases are not shown. The red line is an unadjusted linear fit. This panel is reused unchanged from the original analysis.

(c) Carrier frequency by age at diagnosis, with 95% Wilson intervals. Bars show the proportion of cases carrying a qualifying variant in the 13-gene panel within each age group. Intervals are left-inclusive and right-exclusive; the first group is <35 and the final group is ≥70 years. Counts below each group are total cases, including carriers and non-carriers.

For (a) and (c), the current 13-gene definition includes PTVs plus qualifying P/LP missense variants; raw burden sums ln(OR) once per carried gene. Age refers to the first qualifying recorded EHR diagnosis, rather than necessarily biological onset. The current variant definition remains provisional.
'''
    (MS_10_FIG / 'rare_common_risk_diagnosis_age_captions.md').write_text(caption)
    report = {'n_cases': len(d), 'n_carriers': int(carrier.sum()), 'age_intervals': labels,
              'boundary_rule': 'left-inclusive, right-exclusive', 'plot': 'bar with 95% Wilson intervals',
              'all_carrier_and_noncarrier_cells_at_least_20': True,
              'PRS_panel_reused': str(source.relative_to(MS_10_ROOT)), 'PRS_panel_sha256': digest(target),
              'smallest_format': smallest, 'file_sizes': sizes}
    (MS_10_ROOT / 'results/report/55_age_carrier_frequency_validation.json').write_text(json.dumps(report, indent=2)+'\n')
    print(summary.to_string())

# ============================================================================
# 11. organize-figures: organize_manuscript_figures
# ============================================================================
"""Order existing figures against the 1002 manuscript; retain SVG only."""
from pathlib import Path
import base64
import csv
import shutil
import xml.etree.ElementTree as ET

MS_11_ROOT = Path(__file__).resolve().parents[1]
MS_11_FIG = MS_11_ROOT / 'manuscript/figures'
MS_11_CAP = MS_11_ROOT / 'manuscript/figure_captions'
MS_11_PREVIEW = MS_11_ROOT / 'data/manuscript_refresh_20261006/figure_previews'
MS_11_MAPPING = {
    'four_case_status_models_ROC_sensitivities.svg': 'Figure_03a_case_status_ROC.svg',
    'four_case_status_models_AUC_sensitivities.svg': 'Figure_03b_case_status_AUC.svg',
    'PRS311_decile_13gene_carrier_joint_distribution.svg': 'Figure_04a_PRS311_carrier_joint_distribution.svg',
    'age_at_diagnosis_13gene_burden_carrier_median.svg': 'Figure_05a_age_by_13gene_burden_group.svg',
    '13gene_carrier_frequency_by_age_at_diagnosis.svg': 'Figure_05c_carrier_frequency_by_diagnosis_age.svg',
    '13gene_burden_histogram_percentages.svg': 'Figure_02c_13gene_burden_histogram_percentages.svg',
    'PRS311_decile_PRS_burden_boxplots.svg': 'Figure_02d_PRS311_decile_PRS_burden_boxplots.svg',
    'internal_score_distributions/PRS311_raw_score_KDE.svg': 'Figure_02b_PRS311_raw_score_KDE.svg',
    'internal_score_distributions/burden_raw_score_violin.svg': 'Figure_02a_burden_raw_score_violin.svg',
}


def organize_manuscript_figures():
    for path in [MS_11_CAP, MS_11_PREVIEW]:
        path.mkdir(parents=True, exist_ok=True)
    rows = []
    for old, new in MS_11_MAPPING.items():
        source, target = MS_11_FIG / old, MS_11_FIG / new
        if source.exists():
            source.replace(target)
        assert target.exists(), target
        ET.parse(target)
        rows.append({'original': old, 'organized': new, 'format': 'vector SVG'})
    # Preserve the original panel exactly, without recalculating PRS or bins.
    source = MS_11_FIG / '17_prs311_vs_age_at_diagnosis_hexbin.jpg'
    if not source.exists():
        source = MS_11_ROOT / 'results/figures/17_prs311_vs_age_at_diagnosis_hexbin.jpg'
    target = MS_11_FIG / 'Figure_05b_PRS311_vs_diagnosis_age.svg'
    from PIL import Image
    with Image.open(source) as im:
        w, h = im.size
    encoded = base64.b64encode(source.read_bytes()).decode('ascii')
    target.write_text(f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="{w}" height="{h}" viewBox="0 0 {w} {h}"><title>PRS311 against age at diagnosis</title><desc>Original JPG panel embedded unchanged; raster content, not a vector redraw.</desc><image width="{w}" height="{h}" xlink:href="data:image/jpeg;base64,{encoded}"/></svg>\n')
    ET.parse(target)
    rows.append({'original': source.name, 'organized': target.name, 'format': 'SVG with embedded original JPG'})
    # Move captions out; preserve previews/checkpoints in ignored internal storage.
    for source in list(MS_11_FIG.rglob('*')):
        if not source.is_file():
            continue
        relative = source.relative_to(MS_11_FIG)
        if '.ipynb_checkpoints' not in relative.parts and source.suffix == '.svg':
            continue
        destination = (MS_11_CAP / relative) if source.suffix == '.md' else (MS_11_PREVIEW / relative)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(source), str(destination))
    for directory in sorted((p for p in MS_11_FIG.rglob('*') if p.is_dir()), key=lambda p: len(p.parts), reverse=True):
        if not any(directory.iterdir()):
            directory.rmdir()
    assert all(p.suffix == '.svg' for p in MS_11_FIG.rglob('*') if p.is_file())
    with (MS_11_CAP / 'figure_file_mapping.csv').open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=['original', 'organized', 'format'])
        writer.writeheader()
        writer.writerows(rows)
    (MS_11_CAP / 'figure_order.md').write_text("# Updated main manuscript figure order\n\nFigure 1: study flow; no current SVG available.\nFigure 2a–2d: descriptive score results: raw burden violin; PRS311 raw score KDE; 13-gene burden histogram; PRS311 and burden boxplots across PRS311 deciles.\nFigure 3a/3b: case-status ROC and AUC.\nFigure 4a: PRS311 and carrier joint distribution. Figure 4b has no current new-method SVG here.\nFigure 5a/5b/5c: burden-group age boxplots; unchanged PRS311 hexbin; updated carrier frequency.\nFigure 6: primary and sensitivity estimates, including 13-, 12- and 9-gene comparisons. Caption: `Figure_06_caption.md`.\n\nFigure 5b embeds the unchanged original JPG in SVG; its content remains raster.\nCaptions are stored alongside this index. Original previews and notebook checkpoints are retained in the ignored internal directory `data/manuscript_refresh_20261006/figure_previews/`.\nThe order reflects the author's addition of descriptive score results before model results. The DOCX files have not yet been renumbered. Existing plotting scripts retain their original output names; rerun the organizer after regenerating those figures.\n")
    print('\n'.join(str(p.relative_to(MS_11_FIG)) for p in sorted(MS_11_FIG.rglob('*.svg'))))

# ============================================================================
# 12. sensitivity-figure: plot_primary_sensitivity_estimates
# ============================================================================
"""Figure 6: existing current-panel primary and sensitivity model estimates."""
from pathlib import Path
import hashlib
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

MS_12_ROOT = Path(__file__).resolve().parents[1]
MS_12_VERSIONS = [
    ('C', '13 genes: primary'),
    ('C_without_STK11', '12 genes: exclude STK11'),
    ('B', '9 genes: PTV + qualifying missense'),
    ('C_PTV_only', '13 genes: PTV only'),
    ('C_TP53_NF1_strict_VAF', '13 genes: stricter TP53/NF1 VAF'),
    ('C_ClinVar_min1star', '13 genes: ClinVar ≥1 star'),
    ('C_without_CDH1', '12 genes: exclude CDH1'),
    ('A', '9 genes: original PTV baseline'),
]


def plot_primary_sensitivity_estimates():
    source = MS_12_ROOT / 'results/table/45_rowlands_fixed_cohort_coefficients.csv'
    d = pd.read_csv(source)
    terms = ['burden_z', 'PRS311_z']
    selected = d[d.version.isin([v for v, _ in MS_12_VERSIONS]) & d.term.isin(terms)].copy()
    assert len(selected) == len(MS_12_VERSIONS)*2*2
    assert not selected.duplicated(['cohort', 'version', 'term']).any()
    plt.style.use('default')
    fig, axes = plt.subplots(1, 2, figsize=(13, 6.2), sharey=True)
    y = np.arange(len(MS_12_VERSIONS))
    records = []
    for ax, cohort, title in zip(axes, ['case_control', 'diagnosis_age'],
                                ['(A) Case status (n=141,204)', '(B) Age at diagnosis (n=7,996)']):
        for term, label, color, marker, offset in [('burden_z', 'Rare-variant burden', 'tab:orange', 'o', -.12),
                                                   ('PRS311_z', 'PRS311', 'tab:green', 's', .12)]:
            sub = selected[(selected.cohort == cohort) & selected.term.eq(term)].set_index('version').loc[[v for v, _ in MS_12_VERSIONS]]
            cols = ['OR', 'CI_low', 'CI_high'] if cohort == 'case_control' else ['beta', 'CI_beta_low', 'CI_beta_high']
            est, lo, hi = [sub[c].to_numpy() for c in cols]
            assert np.isfinite([est, lo, hi]).all() and (lo <= est).all() and (est <= hi).all()
            ax.errorbar(est, y+offset, xerr=np.vstack([est-lo, hi-est]), fmt=marker,
                        color=color, capsize=3, markersize=5, linewidth=1.3, label=label)
            for version, e, l, h in zip(sub.index, est, lo, hi):
                records.append({'cohort': cohort, 'version': version, 'term': term,
                                'estimate': e, 'CI_low': l, 'CI_high': h})
        ax.set_title(title, fontsize=12)
        ax.grid(False)
        for spine in ax.spines.values():
            spine.set_visible(True)
        ax.axvline(1 if cohort == 'case_control' else 0, color='.5', linestyle='--', linewidth=.8)
        ax.set_xlabel('Odds ratio per SD' if cohort == 'case_control' else 'Difference in age at diagnosis (years per SD)')
    axes[0].set_xscale('log')
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
    axes[0].xaxis.set_major_locator(FixedLocator([1, 1.2, 1.4, 1.6]))
    axes[0].xaxis.set_major_formatter(FuncFormatter(lambda value, _: f'{value:g}'))
    axes[0].xaxis.set_minor_locator(NullLocator())
    axes[0].set_xlim(.98, 1.65)
    axes[1].set_xlim(-2.25, .12)
    axes[0].set_yticks(y, [label for _, label in MS_12_VERSIONS], fontsize=10)
    axes[0].invert_yaxis()
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, loc='upper center', bbox_to_anchor=(.6, .93), ncol=2, fontsize=10)
    fig.suptitle('Figure 6. Primary and sensitivity estimates', fontsize=14)
    fig.text(.5, .015, 'Horizontal lines: 95% CIs. Each burden is scaled by its own cohort/model SD.\n'
             'Case-status burden OR is conditional on mean PRS311; the age model has no interaction.', ha='center', fontsize=9)
    fig.tight_layout(rect=[0, .07, 1, .87])
    output = MS_12_ROOT / 'manuscript/figures/Figure_06_primary_and_sensitivity_estimates.svg'
    fig.savefig(output, bbox_inches='tight')
    preview = MS_12_ROOT / 'data/manuscript_refresh_20261006/figure_previews/Figure_06_primary_and_sensitivity_estimates.png'
    preview.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(preview, dpi=160, bbox_inches='tight')
    plt.close(fig)
    pd.DataFrame(records).to_csv(MS_12_ROOT / 'results/table/56_figure5_primary_sensitivity_estimates.csv', index=False)
    caption = '''Figure 6. Primary and sensitivity estimates. (A) Case-status odds ratios per SD from the joint logistic model, including PRS311, rare-variant burden, their interaction, age at last EHR observation and PC1–10 (n=141,204). (B) Difference in recorded age at diagnosis per SD among cases from an OLS model including PRS311, burden and PC1–10, without interaction (n=7,996). Points show estimates for burden (orange circles) and PRS311 (green squares); horizontal lines show 95% CIs based on HC1 robust standard errors.

The primary analysis uses the 13-gene PTV plus qualifying P/LP missense burden. Sensitivities exclude STK11 (12 genes), restrict the panel to the nine original genes with the same current PTV/missense rules and Rowlands weights, include PTVs only, apply stricter heterozygous VAF (0.30–0.70) to TP53/NF1, relax the ClinVar review threshold to ≥1 star, or exclude CDH1 (12 genes). The original nine-gene PTV analysis using its original definitions and weights is shown as a historical comparison baseline. It differs from the current nine-gene PTV plus missense analysis in eligibility/QC and weights as well as variant type.

The same fixed participants and PRS311 are used in each model. Burden is scaled by its own SD within each model cohort, so differences in per-SD coefficients do not isolate effects of changing variant definitions. In panel A, the burden OR is conditional on PRS311 at its mean; it is not a marginal carrier OR. The current variant definition remains provisional.
'''
    (MS_12_ROOT / 'manuscript/figure_captions/Figure_06_caption.md').write_text(caption)
    report = {'source': str(source.relative_to(MS_12_ROOT)), 'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
              'versions': [v for v, _ in MS_12_VERSIONS], 'plotted_estimates': len(records),
              'models_refit': False, 'output': str(output.relative_to(MS_12_ROOT))}
    (MS_12_ROOT / 'results/report/56_figure5_validation.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))

# ============================================================================
# 13. insert-docx-figures: insert_manuscript_figures
# ============================================================================
"""Replace the main DOCX figure section with current SVGs and PNG fallbacks."""
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
from copy import deepcopy
import hashlib
import io
import json
import re
import shutil
import cairosvg
from PIL import Image
from lxml import etree as E

MS_13_ROOT = Path(__file__).resolve().parents[1]
MS_13_DOC = MS_13_ROOT / 'manuscript/medrxiv_v1_main_revised_1002.docx'
MS_13_NS = {'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
      'r': 'http://schemas.openxmlformats.org/officeDocument/2006/relationships',
      'a': 'http://schemas.openxmlformats.org/drawingml/2006/main',
      'wp': 'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing',
      'pic': 'http://schemas.openxmlformats.org/drawingml/2006/picture'}
MS_13_REL = 'http://schemas.openxmlformats.org/package/2006/relationships'
MS_13_ASVG = 'http://schemas.microsoft.com/office/drawing/2016/SVG/main'


def MS_13_text(p):
    return ''.join(p.xpath('.//w:t/text()', namespaces=MS_13_NS))


def MS_13_paragraph(value, template=None):
    p = E.Element('{%s}p' % MS_13_NS['w'])
    if template is not None and template.find('w:pPr', MS_13_NS) is not None:
        p.append(deepcopy(template.find('w:pPr', MS_13_NS)))
    run = E.SubElement(p, '{%s}r' % MS_13_NS['w'])
    E.SubElement(run, '{%s}t' % MS_13_NS['w']).text = value
    return p


def MS_13_renumber(p):
    nodes = p.xpath('.//w:t', namespaces=MS_13_NS)
    values = [n.text or '' for n in nodes]
    joined = ''.join(values)
    # Preserve run formatting while replacing only the matched digit.
    for match in reversed(list(re.finditer(r'\bFigure[\s\u00a0]+([2-5])', joined))):
        pos = match.start(1)
        offset = 0
        for node, value in zip(nodes, values):
            if offset <= pos < offset + len(value):
                i = pos-offset
                current = node.text or ''
                node.text = current[:i] + str(int(match.group(1))+1) + current[i+1:]
                break
            offset += len(value)
    # The old joint-distribution OR panel has no replacement; remove its citation.
    values = [n.text or '' for n in nodes]
    joined = ''.join(values)
    for match in reversed(list(re.finditer(r'\s*\(Figure[\s\u00a0]+4b\)', joined))):
        offset = 0
        for node, value in zip(nodes, values):
            left = max(0, match.start()-offset)
            right = min(len(value), match.end()-offset)
            if left < right:
                current = node.text or ''
                node.text = current[:left] + current[right:]
            offset += len(value)


def insert_manuscript_figures():
    backup = MS_13_ROOT / 'manuscript/archive/2026-10-06_before_figure_replacement' / MS_13_DOC.name
    backup.parent.mkdir(parents=True, exist_ok=True)
    if not backup.exists():
        shutil.copy2(MS_13_DOC, backup)
    with ZipFile(backup) as z:
        files = {name: z.read(name) for name in z.namelist()}
    document = E.fromstring(files['word/document.xml'])
    body = document.find('w:body', MS_13_NS)
    children = list(body)
    start = next(i for i, p in enumerate(children) if MS_13_text(p) == 'Figures')
    first_old = next(i for i in range(start, len(children)) if MS_13_text(children[i]) == '(a)')
    picture_template = deepcopy(children[first_old+1])
    caption_template = children[first_old+4]
    for p in children[:start]:
        MS_13_renumber(p)
    # Retain study-flow placeholder and caption; replace all old result figures.
    for p in children[first_old:]:
        if p.tag != '{%s}sectPr' % MS_13_NS['w']:
            body.remove(p)
    relationships = E.fromstring(files['word/_rels/document.xml.rels'])
    max_rid = max(int(x.get('Id')[3:]) for x in relationships if re.fullmatch(r'rId\d+', x.get('Id', '')))
    id_counter = max([int(x.get('id')) for x in document.xpath('.//wp:docPr', namespaces=MS_13_NS)] + [100])
    sect = body.find('w:sectPr', MS_13_NS)
    figdir = MS_13_ROOT / 'manuscript/figures'
    groups = [
        (2, ['02a_burden_raw_score_violin', '02b_PRS311_raw_score_KDE', '02c_13gene_burden_histogram_percentages', '02d_PRS311_decile_PRS_burden_boxplots']),
        (3, ['03a_case_status_ROC', '03b_case_status_AUC']),
        (4, ['04a_PRS311_carrier_joint_distribution']),
        (5, ['05a_age_by_13gene_burden_group', '05b_PRS311_vs_diagnosis_age', '05c_carrier_frequency_by_diagnosis_age']),
        (6, ['06_primary_and_sensitivity_estimates']),
    ]
    capdir = MS_13_ROOT / 'manuscript/figure_captions'
    captions = {
        2: 'Figure 2. Descriptive distributions of genetic scores in the fixed case-control cohort (n=141,204). (a) Raw burden distributions for the 13-gene and 12-gene (without STK11) panels, with main and stricter TP53/NF1 VAF definitions; upper panel includes all participants and lower panel includes carriers only. Violin widths are normalized separately and black bars denote medians. (b) Kernel density estimates of raw PRS311 by case status. (c) Within-group percentages of zero and positive 13-gene burden scores. (d) Raw PRS311 and nonzero 13-gene burden by full-cohort PRS311 deciles; boxes show median and IQR and whiskers observations within 1.5×IQR; outliers are not shown. Case and control distributions are shown in blue and red, respectively. Burden includes PTVs plus qualifying P/LP missense variants.',
        3: 'Figure 3. Discrimination of the four case-status models. ' + (capdir/'four_case_status_models_captions.md').read_text().split('\n\n', 2)[-1].strip(),
        4: 'Figure 4. ' + (capdir/'PRS311_decile_13gene_carrier_joint_distribution_caption.md').read_text().strip(),
        5: 'Figure 5. Rare and common risk and recorded age at diagnosis among cases (n=7,996). (a) Age by 13-gene burden group, with carriers divided at their burden median (1.73695, displayed as 1.74) and non-carriers shown separately. Boxes show median and IQR; whiskers extend to observations within 1.5×IQR; outliers are not shown. (b) PRS311 against age at diagnosis, reused unchanged. Hexagonal bins with fewer than 20 cases are not shown; the red line is an unadjusted linear fit. (c) Carrier proportions within age groups with 95% Wilson intervals. Age groups are left-inclusive and right-exclusive, with <35 and ≥70 as the extreme groups. Carrier status uses the current 13-gene PTV plus qualifying P/LP missense definition.',
        6: (capdir/'Figure_06_caption.md').read_text().strip(),
    }
    inputs = {}
    for number, names in groups:
        for panel, name in enumerate(names):
            path = figdir / f'Figure_{name}.svg'
            svg = path.read_bytes()
            inputs[path.name] = hashlib.sha256(svg).hexdigest()
            png = cairosvg.svg2png(bytestring=svg, output_width=1950)
            with Image.open(io.BytesIO(png)) as image:
                width, height = image.size
            # Fit each panel within the printable page width and height.
            cx = int(6.5*914400)
            cy = int(cx*height/width)
            if cy > 8.0*914400:
                cx = int(cx*(8.0*914400)/cy)
                cy = int(8.0*914400)
            png_rid, svg_rid = f'rId{max_rid+1}', f'rId{max_rid+2}'
            max_rid += 2
            for ext, data, rid in [('png', png, png_rid), ('svg', svg, svg_rid)]:
                media = f'media/current_Figure_{name}.{ext}'
                files['word/'+media] = data
                E.SubElement(relationships, '{%s}Relationship' % MS_13_REL, Id=rid,
                             Type=MS_13_NS['r']+'/image', Target=media)
            if len(names) > 1 or number == 4:
                body.insert(body.index(sect), MS_13_paragraph(f'({chr(97+panel)})', caption_template))
            p = deepcopy(picture_template)
            for t in p.xpath('.//w:t', namespaces=MS_13_NS):
                t.text = ''
            for extent in p.xpath('.//wp:extent | .//a:xfrm/a:ext', namespaces=MS_13_NS):
                extent.set('cx', str(cx)); extent.set('cy', str(cy))
            id_counter += 1
            for prop in p.xpath('.//wp:docPr | .//pic:cNvPr', namespaces=MS_13_NS):
                prop.set('id', str(id_counter));prop.set('name',path.stem);prop.set('descr',path.stem)
            for blip in p.xpath('.//a:blip', namespaces=MS_13_NS):
                blip.set('{%s}embed' % MS_13_NS['r'], png_rid)
                for child in list(blip):
                    blip.remove(child)
                extensions = E.SubElement(blip, '{%s}extLst' % MS_13_NS['a'])
                extension = E.SubElement(extensions, '{%s}ext' % MS_13_NS['a'], uri='{96DAC541-7B7A-43D3-8B79-37D633B846F1}')
                E.SubElement(extension, '{%s}svgBlip' % MS_13_ASVG, nsmap={'asvg':MS_13_ASVG}, attrib={'{%s}embed' % MS_13_NS['r']:svg_rid})
            body.insert(body.index(sect), p)
        for part in captions[number].split('\n\n'):
            body.insert(body.index(sect), MS_13_paragraph(part, caption_template))
    # Drop old figure media/relationships when unused; preserve other document art.
    used = set(document.xpath('.//@r:embed | .//@r:id | .//@r:link', namespaces=MS_13_NS))
    for relation in list(relationships):
        if relation.get('Type') == MS_13_NS['r']+'/image' and relation.get('Id') not in used:
            target = 'word/'+relation.get('Target')
            relationships.remove(relation)
            files.pop(target, None)
    content_types = E.fromstring(files['[Content_Types].xml'])
    ct = 'http://schemas.openxmlformats.org/package/2006/content-types'
    if not any(e.get('Extension') == 'svg' for e in content_types):
        E.SubElement(content_types, '{%s}Default' % ct, Extension='svg', ContentType='image/svg+xml')
    for name, tree in [('word/document.xml', document), ('word/_rels/document.xml.rels', relationships), ('[Content_Types].xml', content_types)]:
        files[name] = E.tostring(tree, xml_declaration=True, encoding='UTF-8', standalone=True)
    temporary = MS_13_DOC.with_suffix('.updated.docx')
    with ZipFile(temporary, 'w', ZIP_DEFLATED) as z:
        for name, data in files.items():
            z.writestr(name, data)
    with ZipFile(temporary) as z:
        assert z.testzip() is None
        d = E.fromstring(z.read('word/document.xml'))
        assert len(d.findall('.//{%s}svgBlip' % MS_13_ASVG)) == 11
        relmap = {x.get('Id'):x for x in relationships}
        for rid in d.xpath('.//@r:embed', namespaces=MS_13_NS):
            assert 'word/'+relmap[rid].get('Target') in z.namelist()
        # All pre-figure paragraphs retain their complete text except figure digits.
        original = E.fromstring(files['word/document.xml'])
        assert len(original.xpath('.//wp:docPr',namespaces=MS_13_NS)) == 12
    temporary.replace(MS_13_DOC)
    report = {'docx':str(MS_13_DOC.relative_to(MS_13_ROOT)), 'backup':str(backup.relative_to(MS_13_ROOT)),
              'inserted_SVG_count':11, 'PNG_fallback_count':11, 'figure_range':'2–6',
              'figure_references_renumbered':True, 'old_joint_OR_panel_removed':True,
              'main_result_numbers_revised':False, 'supplement_docx_modified':False,
              'inputs':inputs, 'docx_sha256':hashlib.sha256(MS_13_DOC.read_bytes()).hexdigest()}
    (MS_13_ROOT/'results/report/57_docx_figure_replacement_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='inputs'},indent=2))
# Command-line selection of the computations defined above.
def main():
    from scripts.functions.cli import run_local_step
    computations = {
        'score-distributions': plot_raw_score_distributions,
        'discrimination': plot_four_model_discrimination,
        'joint-distribution': plot_joint_distribution,
        'score-deciles': plot_score_decile_boxplots,
        'burden-histogram': plot_burden_histogram,
        'age-boxplots': plot_age_burden_boxplots,
        'age-carrier-frequency': plot_age_carrier_frequency,
        'organize-figures': organize_manuscript_figures,
        'sensitivity-figure': plot_primary_sensitivity_estimates,
        'insert-docx-figures': insert_manuscript_figures,
        'attribution': attribute_new_carriers,
        'summary': summarize_model_results,
        'age-figures': plot_age_and_auc_results,
    }
    run_local_step('Final model interpretation, carrier attribution and age/AUC figures', computations, default_tasks=['attribution', 'summary', 'age-figures'])

if __name__ == '__main__':
    main()

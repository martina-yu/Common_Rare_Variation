"""
- Collapse variants to gene carriers.
- Calculate and standardize weighted burden.
- Compare score definitions and sensitivities.
- Prepare model inputs.
"""
from __future__ import annotations
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
# ============================================================================
# 1. score: calculate_burden_scores
# ============================================================================
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

BURDEN_ROOT=Path(__file__).resolve().parents[1]
BURDEN_BASE=BURDEN_ROOT/'data/29_rowlands_extension'
BURDEN_OUT=BURDEN_BASE/'section7'

def read_weights(path, genes):
    table=pd.read_csv(path)
    assert {'gene','OR'}<=set(table), 'Weight table must contain gene and OR'
    assert not table.gene.duplicated().any(), 'Duplicate gene weights'
    assert set(genes)<=set(table.gene), 'Missing panel gene weight'
    values=table.set_index('gene').loc[genes,'OR'].astype(float)
    assert np.isfinite(values).all() and values.gt(0).all(), 'OR must be finite and positive'
    return np.log(values)

def score_panel(frame, genes, logweights, version):
    gene_cols=[f'{g}_{version}' for g in genes]
    ptv_cols=[f'{g}_PTV_{version}' for g in genes]
    miss_cols=[f'{g}_missense_{version}' for g in genes]
    assert frame[gene_cols+ptv_cols+miss_cols].isin([0,1]).all().all()
    for g in genes:
        assert frame[f'{g}_{version}'].equals(frame[[f'{g}_PTV_{version}',f'{g}_missense_{version}']].max(axis=1)), 'Gene union must count once'
    raw=frame[gene_cols].to_numpy(dtype=float) @ logweights.loc[genes].to_numpy()
    mean=float(raw.mean());sd=float(raw.std(ddof=0))
    assert np.isfinite(raw).all() and sd>0, 'Cannot standardize zero-variance burden'
    ptv=frame[ptv_cols].max(axis=1).astype(bool)
    miss=frame[miss_cols].max(axis=1).astype(bool)
    tag=np.select([ptv&miss,ptv&~miss,miss&~ptv],['both','PTV-only','missense-only'],default='non-carrier')
    return pd.DataFrame({'B_raw':raw,'B_standardized':(raw-mean)/sd,
        'carried_gene_count':frame[gene_cols].sum(axis=1),'PTV_gene_count':frame[ptv_cols].sum(axis=1),
        'missense_gene_count':frame[miss_cols].sum(axis=1),'carrier_class':tag}), {'mean':mean,'SD':sd,'ddof':0}

def calculate_burden_scores():
    parser=argparse.ArgumentParser()
    parser.add_argument('--carriers',type=Path,default=BURDEN_OUT/'unweighted_gene_carriers.csv.gz')
    parser.add_argument('--weights',type=Path,default=BURDEN_ROOT/'results/table/29_rowlands_13gene_external_weights.csv')
    parser.add_argument('--output-dir',type=Path,default=BURDEN_OUT)
    args=parser.parse_args();args.output_dir.mkdir(parents=True,exist_ok=True)
    import subprocess
    assert args.output_dir.resolve().is_relative_to(BURDEN_BASE.resolve()), 'Participant scores must stay in internal data directory'
    assert subprocess.run(['git','check-ignore','-q',str(args.output_dir)],cwd=BURDEN_ROOT).returncode==0
    cfg=yaml.safe_load((BURDEN_ROOT/'config/rowlands_13gene_extension.yaml').read_text())
    frame=pd.read_csv(args.carriers)
    ids=set(pd.read_csv(BURDEN_ROOT/cfg['cohorts']['case_control']['frozen_ids']).person_id)
    age=set(pd.read_csv(BURDEN_ROOT/cfg['cohorts']['diagnosis_age']['frozen_ids']).person_id)
    assert set(frame.sample_id)==ids and not frame.sample_id.duplicated().any() and age<=ids
    logweights=read_weights(args.weights,cfg['panels']['main_13'])
    result=frame.copy();params={}
    for version in ['main','tp53_nf1_vaf_sensitive']:
        for label,panel in [('13gene',cfg['panels']['main_13']),('12gene_without_STK11',cfg['panels']['without_stk11_12'])]:
            scored,params[f'{label}_{version}']=score_panel(frame,panel,logweights,version)
            for col in scored:result[f'{col}_{label}_{version}']=scored[col]
            assert abs(scored.B_standardized.mean())<1e-12
            assert np.isclose(scored.B_standardized.std(ddof=0),1,rtol=1e-10,atol=1e-12)
    result.to_csv(args.output_dir/'carrier_burden_case_control.csv.gz',index=False)
    result.loc[result.sample_id.isin(age)].to_csv(args.output_dir/'carrier_burden_diagnosis_age.csv.gz',index=False)
    for g in cfg['panels']['main_13']:
        assert (result[f'{g}_tp53_nf1_vaf_sensitive']<=result[f'{g}_main']).all()
    manifest={'reference_cohort':'entire_frozen_case_control_including_noncarriers','reference_n':len(ids),
        'diagnosis_age_reuses_reference_parameters':True,'parameters':params,
        'carrier_input_sha256':hashlib.sha256(args.carriers.read_bytes()).hexdigest(),
        'weight_input':str(args.weights.relative_to(BURDEN_ROOT)) if args.weights.is_relative_to(BURDEN_ROOT) else str(args.weights),
        'weight_input_sha256':hashlib.sha256(args.weights.read_bytes()).hexdigest(),
        'weights':{g:float(np.exp(w)) for g,w in logweights.items()},
        'class_definition':'both means at least one qualified PTV and at least one qualified missense within that panel, possibly across different genes',
        'noncarrier_standardized_value':'(0-reference_mean)/reference_SD; not forced to zero',
        'final_variant_definition_frozen':False}
    (args.output_dir/'standardization_and_weights_internal.json').write_text(json.dumps(manifest,indent=2)+'\n')
    if args.output_dir.resolve()==BURDEN_OUT.resolve():
        (BURDEN_ROOT/'results/report/38_rowlands_burden_validation.json').write_text(json.dumps({
            'status':'section7_burden_scores_validated','cohort_n':len(ids),'diagnosis_age_n':len(age),
            'cohort_membership_exact_match':True,'same_gene_multiple_variants_count_once':True,
            'weight_reading_independent_of_variant_filtering':True,'standardization_reference':'all_frozen_case_control_including_noncarriers',
            'standardization_ddof':0,'diagnosis_age_reuses_same_parameters':True,
            'raw_standardized_gene_count_four_classes_saved':True,'versions':['13gene_main','12gene_without_STK11_main',
                '13gene_tp53_nf1_vaf_sensitive','12gene_without_STK11_tp53_nf1_vaf_sensitive'],
            'final_variant_definition_frozen':False,'models_run':False},indent=2)+'\n')
        cfg['protocol_sections_received']=list(dict.fromkeys(cfg['protocol_sections_received']+[7]))
        cfg['status']='executed_through_section7_provisional_definition_pending_final_adjudication'
        cfg['burden']['section7_status']='audit_and_unweighted_carriers_and_four_score_versions_generated_validated'
        cfg['burden']['sample_variant_audit_internal']='data/29_rowlands_extension/section7/sample_variant_audit.csv.gz'
        cfg['burden']['scores_internal']='data/29_rowlands_extension/section7/carrier_burden_case_control.csv.gz'
        cfg['burden']['derived_scores_internal']=cfg['burden']['scores_internal']
        cfg['data']['genotype_source']='same_original_v9_VDS_all_retained_candidate_ALT_calls; LGT_LA_to_native_global_ALT_indices'
        cfg['candidate_extraction']['genotype_alt_verification']='validated_original_global_ALT_indices_and_full_normalized_alleles'
        cfg['new_non_ptv']['remaining_details_pending']=['final_manual_adjudication_of_held_out_records_and_PTV_definition']
        cfg['new_non_ptv']['missing_af_behavior']='unknown_or_disagreeing_AF_held_out'
        cfg['genotype_qc']['af_missing_policy']='unknown_or_disagreeing_AF_held_out'
        cfg['execution_gate']['reason']='Through section7 call audit, current clinical eligibility and scores executed; final variant adjudication and remaining model/output protocol pending'
        cfg['burden']['standardization']={'reference':'entire_frozen_case_control_including_noncarriers','ddof':0,'diagnosis_age_reuses_parameters':True}
        (BURDEN_ROOT/'config/rowlands_13gene_extension.yaml').write_text(yaml.safe_dump(cfg,allow_unicode=True,sort_keys=False))
        (BURDEN_ROOT/'results/report/38_rowlands_section7_summary.md').write_text('\n'.join(report)+'\n')
    print('Raw and standardized burden, gene counts and four carrier classes validated for all score versions; participant outputs internal.')

# ============================================================================
# 2. validate-score: validate_burden_scores
# ============================================================================
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

BURDEN_VALIDATION_ROOT=Path(__file__).resolve().parents[1]
BURDEN_VALIDATION_OUT=BURDEN_VALIDATION_ROOT/'data/29_rowlands_extension/section7'

def validate_burden_scores():
    spec=importlib.util.spec_from_file_location('audit',BURDEN_VALIDATION_ROOT/'scripts/functions/sample_variant_audit.py')
    audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)
    cfg=yaml.safe_load((BURDEN_VALIDATION_ROOT/'config/rowlands_13gene_extension.yaml').read_text())
    unweighted=pd.read_csv(BURDEN_VALIDATION_OUT/'unweighted_gene_carriers.csv.gz').set_index('sample_id')
    ids=set(unweighted.index)
    reconstructed={g:{cl:{v:set() for v in ['main','tp53_nf1_vaf_sensitive']} for cl in ['PTV','missense']} for g in cfg['panels']['main_13']}
    rejected_observed=False
    for chunk in pd.read_csv(BURDEN_VALIDATION_OUT/'sample_variant_audit.csv.gz',chunksize=20000,keep_default_na=False):
        assert set(audit.FIELDS)<=set(chunk.columns)
        assert set(chunk.sample_id)<=ids
        assert chunk['PTV/missense'].isin(['PTV','missense']).all()
        assert chunk.keep_flag.isin([0,1]).all() and chunk.sensitivity_keep_flag.isin([0,1]).all()
        assert (chunk.sensitivity_keep_flag<=chunk.keep_flag).all()
        assert (chunk.loc[chunk.keep_flag==1,'rejection_reason']=='').all()
        assert (chunk.loc[chunk.keep_flag==0,'rejection_reason']!='').all()
        rejected_observed=rejected_observed or chunk.keep_flag.eq(0).any()
        kept=chunk.loc[chunk.keep_flag==1]
        assert kept.site_filter.eq('PASS').all()
        assert kept.genotype_filter.isin(['PASS','.']).all()
        gq=pd.to_numeric(kept.GQ);vaf=pd.to_numeric(kept.VAF);af=pd.to_numeric(kept.EUR_AF)
        assert np.isfinite(gq).all() and gq.ge(20).all() and np.isfinite(af).all() and af.lt(.01).all()
        assert kept.allele_dosage.isin([1,2]).all()
        assert ((kept.allele_dosage.eq(1)&vaf.between(.3,.8))|(kept.allele_dosage.eq(2)&vaf.ge(.8))).all()
        ms=kept.loc[kept['PTV/missense']=='missense']
        assert ms.low_penetrance_flag.astype(str).str.lower().eq('false').all()
        assert ms.conflict_flag.astype(str).str.lower().eq('false').all()
        for version,col in [('main','keep_flag'),('tp53_nf1_vaf_sensitive','sensitivity_keep_flag')]:
            for (gene,cl),d in chunk.loc[chunk[col]==1].groupby(['gene','PTV/missense']):
                reconstructed[gene][cl][version].update(d.sample_id)
    for gene,classes in reconstructed.items():
        for cl,d in classes.items():
            for version,members in d.items():
                assert set(unweighted.index[unweighted[f'{gene}_{cl}_{version}']==1])==members
        for version in ['main','tp53_nf1_vaf_sensitive']:
            assert set(unweighted.index[unweighted[f'{gene}_{version}']==1])==classes['PTV'][version]|classes['missense'][version]
    previous=pd.read_csv(BURDEN_VALIDATION_ROOT/'data/29_rowlands_extension/section6/derived_gene_carriers_and_scores.csv.gz').set_index('person_id')
    for gene in cfg['panels']['main_13']:
        assert set(previous.index[previous[f'{gene}_PTV']==1])==reconstructed[gene]['PTV']['main'], 'Clinical missense fix changed PTV carriers'
    scored=pd.read_csv(BURDEN_VALIDATION_OUT/'carrier_burden_case_control.csv.gz').set_index('sample_id')
    age=pd.read_csv(BURDEN_VALIDATION_OUT/'carrier_burden_diagnosis_age.csv.gz').set_index('sample_id')
    assert set(scored.index)==ids and not scored.index.duplicated().any()
    assert np.allclose(scored.loc[age.index].select_dtypes('number'),age.select_dtypes('number'),atol=1e-12)
    assert scored.loc[age.index].select_dtypes('object').equals(age.select_dtypes('object'))
    params=json.loads((BURDEN_VALIDATION_OUT/'standardization_and_weights_internal.json').read_text())
    for version in ['main','tp53_nf1_vaf_sensitive']:
        for label,panel in [('13gene',cfg['panels']['main_13']),('12gene_without_STK11',cfg['panels']['without_stk11_12'])]:
            suffix=f'{label}_{version}'
            ptv=scored[[f'{g}_PTV_{version}' for g in panel]].max(axis=1).astype(bool)
            miss=scored[[f'{g}_missense_{version}' for g in panel]].max(axis=1).astype(bool)
            expected=np.select([ptv&miss,ptv&~miss,miss&~ptv],['both','PTV-only','missense-only'],default='non-carrier')
            assert (scored[f'carrier_class_{suffix}']==expected).all()
            assert (scored[f'carried_gene_count_{suffix}']==scored[[f'{g}_{version}' for g in panel]].sum(axis=1)).all()
            raw=sum(np.log(params['weights'][g])*scored[f'{g}_{version}'] for g in panel)
            assert np.allclose(scored[f'B_raw_{suffix}'],raw,atol=1e-12)
            p=params['parameters'][suffix]
            assert np.allclose(scored[f'B_standardized_{suffix}'],(raw-p['mean'])/p['SD'],atol=1e-12)
    result={'status':'section7_persisted_outputs_validated','all_requested_audit_fields_present':True,
        'rejected_calls_retained':bool(rejected_observed),'persisted_audit_reconstructs_all_carriers':True,
        'keep_flags_and_QC_boundaries_verified':True,'raw_and_standardized_formulas_verified':True,
        'four_carrier_classes_and_gene_counts_verified':True,'diagnosis_age_values_identical_to_reference_subset':True,
        'participant_outputs_internal_only':True,'final_definition_frozen':False,'association_models_run':False}
    result['PTV_carriers_unchanged_by_clinical_missense_fix']=True
    (BURDEN_VALIDATION_ROOT/'results/report/39_rowlands_section7_execution_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    historical=pd.read_csv(BURDEN_VALIDATION_ROOT/'results/table/35_rowlands_main_carrier_summary.csv').set_index('cohort')
    public_rows=[]
    for label,d in [('case_control_frozen',scored),('diagnosis_age_frozen',age)]:
        n=int(d.carried_gene_count_13gene_main.gt(0).sum())
        delta=abs(n-int(historical.loc[label,'main_13gene_carriers']))
        safe=n>=20 and len(d)-n>=20 and (delta==0 or delta>=20)
        public_rows.append({'cohort':label,'cohort_n':len(d),
            'main_13gene_carriers':str(n) if safe else 'suppressed_to_avoid_disclosing_small_change',
            'status':'current_section7_clinical_reconciliation'})
    pd.DataFrame(public_rows).to_csv(BURDEN_VALIDATION_ROOT/'results/table/39_rowlands_current_main_carrier_summary.csv',index=False)
    summary=BURDEN_VALIDATION_ROOT/'results/report/38_rowlands_section7_summary.md'
    text=summary.read_text().split('## Final acceptance',1)[0]
    text+='## Final acceptance\n\nSave the audit table and reconstruct the carrier, row-by-row keep/QC conditions, four types of markers, number of genes carried, original/standardized formulas and age subsets, all verified. PTV carrier is not affected by this clinical screening revision. The total number of current main 13 genes can be found in `results/table/39_rowlands_current_main_carrier_summary.csv`; if the change from the existing public total is less than 20, the new precise total will not be disclosed.'
    summary.write_text(text)
    print('Persisted audit, carrier matrix, four classifications, gene counts, raw/standardized formulas and age subset all validated.')

# ============================================================================
# 3. compare: calculate_comparison_scores
# ============================================================================
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

SCORE_COMPARISON_ROOT=Path(__file__).resolve().parents[1]
SCORE_COMPARISON_BASE=SCORE_COMPARISON_ROOT/'data/29_rowlands_extension'
SCORE_COMPARISON_OUT=SCORE_COMPARISON_BASE/'section8'

def calculate_comparison_scores():
    parser=argparse.ArgumentParser()
    parser.add_argument('--weights',type=Path,default=SCORE_COMPARISON_ROOT/'results/table/29_rowlands_13gene_external_weights.csv')
    parser.add_argument('--output-dir',type=Path,default=SCORE_COMPARISON_OUT)
    args=parser.parse_args()
    assert args.output_dir.resolve().is_relative_to(SCORE_COMPARISON_BASE.resolve())
    args.output_dir.mkdir(parents=True,exist_ok=True)
    import subprocess
    assert subprocess.run(['git','check-ignore','-q',str(args.output_dir)],cwd=SCORE_COMPARISON_ROOT).returncode==0
    spec=importlib.util.spec_from_file_location('score',SCORE_COMPARISON_ROOT/'scripts/functions/burden.py')
    score=importlib.util.module_from_spec(spec);spec.loader.exec_module(score)
    cfg=yaml.safe_load((SCORE_COMPARISON_ROOT/'config/rowlands_13gene_extension.yaml').read_text())
    legacy_cfg=yaml.safe_load((SCORE_COMPARISON_ROOT/'config/aim1_config.yaml').read_text())['step3_qualifying_variant']
    nine=cfg['panels']['baseline_9'];thirteen=cfg['panels']['main_13'];twelve=cfg['panels']['without_stk11_12']
    assert set(nine)==set(legacy_cfg['panel_genes_sensitivity_9']) and len(nine)==9
    assert set(thirteen)-set(nine)=={'CDH1','NF1','PTEN','STK11'}
    assert set(twelve)==set(thirteen)-{'STK11'}
    assert 'MAP3K1' not in set(thirteen)|set(nine)
    ids=set(pd.read_csv(SCORE_COMPARISON_ROOT/cfg['cohorts']['case_control']['frozen_ids']).person_id)
    age=set(pd.read_csv(SCORE_COMPARISON_ROOT/cfg['cohorts']['diagnosis_age']['frozen_ids']).person_id)
    current=pd.read_csv(SCORE_COMPARISON_BASE/'section7/unweighted_gene_carriers.csv.gz').set_index('sample_id').loc[sorted(ids)].reset_index()
    legacy_path=SCORE_COMPARISON_ROOT/'data/13_carrier_status_main_ptv_maf1pct_person_level.csv.gz'
    legacy=pd.read_csv(legacy_path).set_index('person_id')
    assert not legacy.index.duplicated().any() and ids<=set(legacy.index)
    legacy=legacy.loc[sorted(ids)]
    oldweights=pd.Series({g:float(legacy_cfg['chek2_ptv_or']['main']['or'] if g=='CHEK2' else legacy_cfg['gene_ptv_or'][g]['or']) for g in nine})
    assert np.isfinite(oldweights).all() and oldweights.gt(0).all()
    old_log=np.log(oldweights);rowlands_log=score.read_weights(args.weights,thirteen)
    pd.DataFrame({'gene':nine,'OR':[oldweights[g] for g in nine],'ln_OR':[old_log[g] for g in nine]}).to_csv(
        SCORE_COMPARISON_ROOT/'results/table/41_rowlands_legacy_nine_gene_weights.csv',index=False)
    oldframe=pd.DataFrame({'sample_id':sorted(ids)})
    ptv_current=current.copy()
    for g in nine:
        oldframe[f'{g}_main']=legacy[g].to_numpy()
        oldframe[f'{g}_PTV_main']=legacy[g].to_numpy()
        oldframe[f'{g}_missense_main']=0
        ptv_current[f'{g}_main']=current[f'{g}_PTV_main']
        ptv_current[f'{g}_missense_main']=0
    # The user-defined A-prime changes ONLY weights. A-prime_QC is explicit,
    # supplemental and isolates PTV eligibility/call handling from missense.
    definitions={
        'A':(oldframe,nine,old_log,'original nine-gene PTV carrier matrix and original weights'),
        'A_prime':(oldframe,nine,rowlands_log,'same original nine-gene PTV carriers as A; Rowlands all-PV weights; decomposition comparator'),
        'A_prime_QC':(ptv_current,nine,rowlands_log,'supplemental bridge: current uniformly QC-filtered nine-gene PTV-only carriers; Rowlands weights'),
        'B':(current,nine,rowlands_log,'current nine-gene PTV plus qualifying missense; Rowlands weights'),
        'C':(current,thirteen,rowlands_log,'current thirteen-gene PTV plus qualifying missense; Rowlands weights'),
        'C_without_STK11':(current,twelve,rowlands_log,'same C definitions and weights, excluding STK11')}
    result=pd.DataFrame({'sample_id':sorted(ids)});params={}
    for label,(frame,panel,weights,_) in definitions.items():
        scored,params[label]=score.score_panel(frame,panel,weights,'main')
        for col in scored:result[f'{col}_{label}']=scored[col]
        assert abs(scored.B_standardized.mean())<1e-12
        assert np.isclose(scored.B_standardized.std(ddof=0),1,rtol=1e-10,atol=1e-12)
    original=pd.read_csv(SCORE_COMPARISON_ROOT/'data/14_burden_rare_scores_person_level.csv.gz').set_index('person_id').loc[sorted(ids)]
    assert np.allclose(result.B_raw_A,original.burden_rare_main_ptv_9gene,rtol=1e-12,atol=1e-12),'A does not reproduce original nine-gene score'
    assert result.carried_gene_count_A.equals(result.carried_gene_count_A_prime)
    assert result.carrier_class_A.equals(result.carrier_class_A_prime)
    current_scores=pd.read_csv(SCORE_COMPARISON_BASE/'section7/carrier_burden_case_control.csv.gz').set_index('sample_id').loc[sorted(ids)]
    default_weights=args.weights.resolve()==(SCORE_COMPARISON_ROOT/'results/table/29_rowlands_13gene_external_weights.csv').resolve()
    if default_weights:
        for new,old in [('C','13gene_main'),('C_without_STK11','12gene_without_STK11_main')]:
            for col in ['B_raw','B_standardized','carried_gene_count']:
                assert np.allclose(result[f'{col}_{new}'],current_scores[f'{col}_{old}'],rtol=1e-12,atol=1e-12)
            assert (result[f'carrier_class_{new}'].to_numpy()==current_scores[f'carrier_class_{old}'].to_numpy()).all()
    deltas={
        'delta_weights_A_prime_minus_A':result.B_raw_A_prime-result.B_raw_A,
        'delta_PTV_handling_A_prime_QC_minus_A_prime':result.B_raw_A_prime_QC-result.B_raw_A_prime,
        'delta_missense_B_minus_A_prime_QC':result.B_raw_B-result.B_raw_A_prime_QC,
        'delta_added_genes_C_minus_B':result.B_raw_C-result.B_raw_B,
        'delta_STK11_C_minus_C_without_STK11':result.B_raw_C-result.B_raw_C_without_STK11,
        'delta_total_C_minus_A':result.B_raw_C-result.B_raw_A,
        'delta_B_minus_A_prime_combines_missense_and_PTV_handling':result.B_raw_B-result.B_raw_A_prime}
    for label,values in deltas.items():result[label]=values
    ptv=current[[f'{g}_PTV_main' for g in nine]].to_numpy()
    full=current[[f'{g}_main' for g in nine]].to_numpy()
    assert np.allclose(result.delta_missense_B_minus_A_prime_QC,(full-ptv)@rowlands_log.loc[nine].to_numpy(),rtol=1e-12,atol=1e-12)
    added=[g for g in thirteen if g not in nine]
    assert np.allclose(result.delta_added_genes_C_minus_B,current[[f'{g}_main' for g in added]].to_numpy()@rowlands_log.loc[added].to_numpy(),rtol=1e-12,atol=1e-12)
    assert np.allclose(result.delta_STK11_C_minus_C_without_STK11,current.STK11_main*rowlands_log['STK11'],rtol=1e-12,atol=1e-12)
    assert np.allclose(result.delta_total_C_minus_A,sum(result[c] for c in [
        'delta_weights_A_prime_minus_A','delta_PTV_handling_A_prime_QC_minus_A_prime',
        'delta_missense_B_minus_A_prime_QC','delta_added_genes_C_minus_B']),rtol=1e-12,atol=1e-12)
    result.to_csv(args.output_dir/'comparable_scores_case_control.csv.gz',index=False)
    result.loc[result.sample_id.isin(age)].to_csv(args.output_dir/'comparable_scores_diagnosis_age.csv.gz',index=False)
    saved=pd.read_csv(args.output_dir/'comparable_scores_case_control.csv.gz').set_index('sample_id')
    saved_age=pd.read_csv(args.output_dir/'comparable_scores_diagnosis_age.csv.gz').set_index('sample_id')
    assert set(saved.index)==ids and set(saved_age.index)==age and not saved.index.duplicated().any()
    assert np.allclose(saved.loc[saved_age.index].select_dtypes('number'),saved_age.select_dtypes('number'),rtol=1e-12,atol=1e-12)
    assert saved.loc[saved_age.index].select_dtypes('object').equals(saved_age.select_dtypes('object'))
    assert np.allclose(saved.filter(regex=r'^(B_raw_|B_standardized_|delta_)'),
        result.set_index('sample_id').filter(regex=r'^(B_raw_|B_standardized_|delta_)'),rtol=1e-12,atol=1e-12)
    agreement={g:bool(np.array_equal(legacy[g].to_numpy(),current[f'{g}_PTV_main'].to_numpy())) for g in nine}
    inputs=[legacy_path,SCORE_COMPARISON_ROOT/'data/14_burden_rare_scores_person_level.csv.gz',SCORE_COMPARISON_BASE/'section7/unweighted_gene_carriers.csv.gz',
        SCORE_COMPARISON_ROOT/'config/aim1_config.yaml',args.weights]
    lineage={'standardization_reference':'all_frozen_case_control_including_noncarriers','ddof':0,
        'diagnosis_age_reuses_reference_parameters':True,'parameters':params,'legacy_OR':oldweights.to_dict(),
        'Rowlands_OR':{g:float(np.exp(rowlands_log[g])) for g in thirteen},
        'definitions':{k:v[3] for k,v in definitions.items()},'legacy_and_current_PTV_carriers_identical_by_gene':agreement,
        'inputs_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        'differences_decompose_raw_scores_only':'each standardized version has its own reference mean/SD; standardized differences do not isolate additive contributions'}
    (args.output_dir/'score_definitions_and_standardization_internal.json').write_text(json.dumps(lineage,indent=2)+'\n')
    result.filter(regex=r'^B_raw_').corr().to_csv(args.output_dir/'raw_score_correlations_internal.csv')
    if args.output_dir.resolve()==SCORE_COMPARISON_OUT.resolve():
        metadata={'status':'section8_comparable_scores_validated','cohort_n':len(ids),'diagnosis_age_n':len(age),
            'A_reproduces_original_nine_gene_score':True,'A_prime_changes_weights_only':True,
            'legacy_and_current_PTV_carriers_identical':all(agreement.values()),'A_prime_QC_bridge_generated':True,
            'B_minus_A_prime_is_not_pure_missense_comparison':not all(agreement.values()),
            'B_minus_A_prime_QC_is_missense_only_comparison':True,'C_minus_B_adds_only_CDH1_NF1_PTEN_STK11':True,
            'C_minus_C_without_STK11_is_exact_STK11_contribution':True,'MAP3K1_included':False,
            'C_and_C_without_STK11_match_section7':default_weights,'raw_delta_decomposition_verified':True,
            'age_subset_reuses_standardization':True,'versions_saved':list(definitions),
            'persisted_scores_and_age_subset_validated':True,
            'all_PV_weights_for_PTV_only_are_decomposition_comparator_not_definition_matched_risk_score':True,
            'participant_outputs_internal_only':True,'final_variant_definition_frozen':False,'models_run':False}
        (SCORE_COMPARISON_ROOT/'results/report/41_rowlands_comparable_scores_validation.json').write_text(json.dumps(metadata,indent=2)+'\n')
        cfg['protocol_sections_received']=list(dict.fromkeys(cfg['protocol_sections_received']+[8]))
        cfg['status']='executed_through_section8_provisional_definition_pending_final_adjudication'
        cfg['burden']['section8_status']='five_requested_versions_plus_QC_bridge_generated_and_validated'
        cfg['burden']['comparable_scores_internal']='data/29_rowlands_extension/section8/comparable_scores_case_control.csv.gz'
        cfg['burden']['A_prime_definition']='same_legacy_nine_gene_PTV_carriers_as_A_only_weights_changed_decomposition_comparator'
        cfg['burden']['missense_isolation_comparison']='B_minus_A_prime_QC_not_B_minus_A_prime_when_PTV_handling_differs'
        cfg['execution_gate']['reason']='Through section8 audited calls, clinical eligibility and comparable scores executed; final variant adjudication and remaining model/output protocol pending'
        (SCORE_COMPARISON_ROOT/'config/rowlands_13gene_extension.yaml').write_text(yaml.safe_dump(cfg,allow_unicode=True,sort_keys=False))
        report=['# Section 8: Comparable Fractions','',
            'Five specified versions have been calculated; A strictly reproduces the old nine-gene main PTV score, A′ uses the same old carrier as A, and only replaces the weights. B/C adopts the current unified QC and reviewed clinical qualifications. MAP3K1 does not come into these versions.','',
            '| version | carrier / mutation | weight | purpose |','|---|---|---|---|',
            '| A | Original Nine Genes old PTV carrier | Original weight | Old benchmark, reproduce the original No. 14 Nine Genes score person by person |',
            '| A′ | Exactly the same as A | Rowlands all-PV | Comparison tool with weights replaced only |',
            '| A′_QC (Supplementary Bridge) | Nine Genes’ current unified QC’s PTV-only | Rowlands all-PV | Split PTV qualification/call processing changes |',
            '| B | Nine genes current PTV + qualified missense | Rowlands all-PV | Compare with A′_QC to measure missense alone |',
            '| C | 13-gene current PTV + qualified missense | Rowlands all-PV | Main extension score |',
            '| C−STK11 | C to STK11 | The remaining weights remain unchanged | Check the impact of STK11 imprecision OR |','',
            '**Comparison limitation:** The old PTV carrier is not completely consistent with the current unified QC set, so B−A′ contains both missense and PTV qualification/call processing changes, which cannot be attributed to missense as a single factor. Five user versions are retained and the A′_QC bridge is added.',
            'A′ uses all-PV weights for PTV-only and is only a comparison tool to disassemble the sources of variation; it should not be called a risk score that exactly matches the literature variation definition.','',
            'The original score is decomposed as A→A′→A′_QC→B→C; the sum of the four-step differences is strictly equal to C−A. C−B is strictly equal to the contribution of newly added CDH1, NF1, PTEN, and STK11; C−(C−STK11) is strictly equal to the contribution of STK11.','',
            'The six versions all preserve the original/standardized B, the number of genes carried and the four types of markers; the standardization uses the mean value of each version and the overall SD of all 141,204 people (including non-carriers) as a reference, and the age cohort of 7,996 people reuses the parameters. Source disassembly uses raw score differences; differences after each version is standardized separately cannot be interpreted as a linear contribution.',
            'Person-by-person scores, standardized parameters, differences and correlation matrices are only stored in the internal section8 directory. Undisclosed parallel sparse inferable version counts or moments.',
            'A The old results are retained; B/C are still tentative PTV definitions, and the final manual decision has not been completed. Extended correlation model was not run.']
        (SCORE_COMPARISON_ROOT/'results/report/41_rowlands_comparable_scores_summary.md').write_text('\n'.join(report)+'\n')
    print('A/A-prime legacy identity, current-QC bridge, B/C/STK11 contributions, raw decomposition and six score versions validated; participant outputs internal.')

# ============================================================================
# 4. model-inputs: prepare_model_inputs
# ============================================================================
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

MODEL_PREPARATION_ROOT=Path(__file__).resolve().parents[1]
MODEL_PREPARATION_BASE=MODEL_PREPARATION_ROOT/'data/29_rowlands_extension';MODEL_PREPARATION_OUT=MODEL_PREPARATION_BASE/'section9'
PC=[f'PC{i}' for i in range(1,11)]

def model_preparation_sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def prepare_model_inputs():
    MODEL_PREPARATION_OUT.mkdir(parents=True,exist_ok=True)
    import subprocess
    assert subprocess.run(['git','check-ignore','-q',str(MODEL_PREPARATION_OUT)],cwd=MODEL_PREPARATION_ROOT).returncode==0
    cfg=yaml.safe_load((MODEL_PREPARATION_ROOT/'config/rowlands_13gene_extension.yaml').read_text())
    scores=pd.read_csv(MODEL_PREPARATION_BASE/'section8/comparable_scores_case_control.csv.gz').set_index('sample_id')
    age_scores=pd.read_csv(MODEL_PREPARATION_BASE/'section8/comparable_scores_diagnosis_age.csv.gz').set_index('sample_id')
    current=pd.read_csv(MODEL_PREPARATION_BASE/'section7/unweighted_gene_carriers.csv.gz').set_index('sample_id')
    ids=set(pd.read_csv(MODEL_PREPARATION_ROOT/cfg['cohorts']['case_control']['frozen_ids']).person_id)
    age=set(pd.read_csv(MODEL_PREPARATION_ROOT/cfg['cohorts']['diagnosis_age']['frozen_ids']).person_id)
    assert set(scores.index)==ids and set(current.index)==ids and set(age_scores.index)==age and age<=ids
    assert not scores.index.duplicated().any()
    pheno=pd.read_csv(MODEL_PREPARATION_ROOT/'data/04_phenotype_diagnostics_person_level.csv.gz').set_index('person_id')
    qc=pd.read_csv(MODEL_PREPARATION_ROOT/'data/05_cohort_qc_flags.csv').set_index('person_id')
    pcs=pd.read_csv(MODEL_PREPARATION_ROOT/'data/09_ancestry_pc_person_level.csv.gz').set_index('person_id')
    clinical=pd.read_csv(MODEL_PREPARATION_ROOT/'data/09_extended_clinical_person_level.csv.gz').set_index('person_id')
    prs=pd.read_csv(MODEL_PREPARATION_ROOT/'data/prs_scores.csv').set_index('person_id')
    previous=pd.read_csv(MODEL_PREPARATION_ROOT/'data/16_joint_model_cv_oof_person_level.csv.gz').set_index('person_id')
    assert set(previous.index)==ids
    order=previous.index
    cohort=pheno.loc[order,['group','is_case','age_at_last_obs','age_at_diagnosis','bc_is_first_record']].copy()
    cohort=cohort.join(pcs[PC],how='left').join(clinical[['seq_site_id','seq_sequencer']],how='left').join(prs[['PRS311']],how='left')
    assert cohort[['age_at_last_obs','PRS311']+PC+['seq_site_id','seq_sequencer']].notna().all().all()
    assert qc.loc[order,'pass_qc'].eq(1).all()
    assert cohort.group.equals(previous.group)
    assert (cohort.is_case.astype(int)==cohort.group.eq('case').astype(int)).all()
    assert set(cohort.index[cohort.group.eq('case')&cohort.bc_is_first_record.ne(True)])==age
    assert cohort.loc[sorted(age),'age_at_diagnosis'].notna().all()
    legacy_scores=pd.read_csv(MODEL_PREPARATION_ROOT/'data/14_burden_rare_scores_person_level.csv.gz').set_index('person_id')
    assert np.array_equal(previous.carrier_main_10gene,(legacy_scores.loc[order,'burden_rare_main_ptv_10gene']>0).astype(int))
    from sklearn.model_selection import RepeatedStratifiedKFold
    strata=cohort.is_case.astype(int).astype(str)+'_'+previous.carrier_main_10gene.astype(int).astype(str)
    splits=list(RepeatedStratifiedKFold(n_splits=5,n_repeats=10,random_state=20260818).split(np.arange(len(cohort)),strata))
    fold=np.full(len(cohort),-1)
    for k,(_,test) in enumerate(splits[:5]):fold[test]=k
    assert np.array_equal(fold,previous.fold_r1.to_numpy()),'Original fold assignment not reproduced'
    # Derived model inputs only, fixed cohort and original order for paired CV.
    cohort.index.name='person_id';cohort.to_csv(MODEL_PREPARATION_OUT/'frozen_model_inputs.csv.gz')
    gene_rows=[];added_rows=[];multi_rows=[];dist=[]
    nine=cfg['panels']['baseline_9'];added=['CDH1','NF1','PTEN','STK11']
    old_legacy=pd.read_csv(MODEL_PREPARATION_ROOT/'data/13_carrier_status_main_ptv_maf1pct_person_level.csv.gz').set_index('person_id')
    for label,members in [('case_control',ids),('diagnosis_age',age)]:
        d=current.loc[sorted(members)]
        old_current=d[[f'{g}_main' for g in nine]].max(axis=1).astype(bool)
        old_actual=old_legacy.loc[d.index,nine].max(axis=1).astype(bool)
        for gene in cfg['panels']['main_13']:
            p=d[f'{gene}_PTV_main'].astype(bool);m=d[f'{gene}_missense_main'].astype(bool)
            gene_rows.append({'cohort':label,'gene':gene,'PTV_carriers':int(p.sum()),'missense_carriers':int(m.sum()),
                'PTV_and_missense_carriers':int((p&m).sum()),'missense_only_added_carriers':int((m&~p).sum()),
                'union_carriers':int((p|m).sum())})
            if gene in added:
                c=p|m
                added_rows.append({'cohort':label,'gene':gene,'gene_carriers':int(c.sum()),
                    'overlap_with_current_QC_nine':int((c&old_current).sum()),'outside_current_QC_nine':int((c&~old_current).sum()),
                    'overlap_with_legacy_nine':int((c&old_actual).sum()),'outside_legacy_nine':int((c&~old_actual).sum())})
        counts=d[[f'{g}_main' for g in cfg['panels']['main_13']]].sum(axis=1)
        for k in range(14):multi_rows.append({'cohort':label,'carried_gene_count':k,'participants':int(counts.eq(k).sum())})
        multi_rows.append({'cohort':label,'carried_gene_count':'at_least_2','participants':int(counts.ge(2).sum())})
        s=scores.loc[sorted(members)]
        for col in s.filter(regex='^B_raw_'):
            version=col.removeprefix('B_raw_');x=s[col]
            dist.append({'cohort':label,'version':version,'n':len(x),'carriers':int(s[f'carried_gene_count_{version}'].gt(0).sum()),
                'zero_n':int(x.eq(0).sum()),'zero_fraction':float(x.eq(0).mean()),'mean':float(x.mean()),
                'SD_population_ddof0':float(x.std(ddof=0)),'SD_sample_ddof1':float(x.std(ddof=1)),
                **{f'p{int(q*100):02d}':float(x.quantile(q)) for q in [0,.01,.25,.5,.75,.95,.99,1]}})
    pd.DataFrame(gene_rows).to_csv(MODEL_PREPARATION_OUT/'gene_PTV_missense_carriers_internal.csv',index=False)
    pd.DataFrame(added_rows).to_csv(MODEL_PREPARATION_OUT/'added_gene_contribution_overlap_internal.csv',index=False)
    pd.DataFrame(multi_rows).to_csv(MODEL_PREPARATION_OUT/'multigene_counts_internal.csv',index=False)
    pd.DataFrame(dist).to_csv(MODEL_PREPARATION_OUT/'score_distributions_and_SD_internal.csv',index=False)
    # Extreme records are controlled participant/variant outputs, never public.
    extrema=[];extreme_ids=set()
    for label,ss in [('case_control',scores),('diagnosis_age',age_scores)]:
        for col in ss.filter(regex='^B_raw_'):
            version=col.removeprefix('B_raw_');z=ss.nlargest(20,col)
            for person,row in z.iterrows():extrema.append({'cohort':label,'version':version,'sample_id':person,'B_raw':row[col]});extreme_ids.update(z.index)
    pd.DataFrame(extrema).to_csv(MODEL_PREPARATION_OUT/'extreme_score_samples_internal.csv.gz',index=False)
    with __import__('gzip').open(MODEL_PREPARATION_OUT/'extreme_score_variant_calls_internal.csv.gz','wt') as f:
        first=True
        for chunk in pd.read_csv(MODEL_PREPARATION_BASE/'section7/sample_variant_audit.csv.gz',chunksize=20000,keep_default_na=False,low_memory=False):
            part=chunk.loc[chunk.sample_id.isin(extreme_ids)]
            part.to_csv(f,index=False,header=first);first=False
    # Legacy A/A-prime extremes retain exact old qualifying/carrier evidence;
    # their variant linkage remains separate where old handling differs.
    pd.DataFrame({'sample_id':sorted(extreme_ids)}).to_csv(MODEL_PREPARATION_OUT/'extreme_samples_for_legacy_PTV_review_internal.csv',index=False)
    for name in ['ptv_hc_transcript_position_review.tsv.gz','legacy_ptv_manual_review.tsv.gz','ptv_definition_comparison.tsv.gz']:
        shutil.copyfile(MODEL_PREPARATION_BASE/'ptv_audit'/name,MODEL_PREPARATION_OUT/name)
    candidates=pd.read_csv(MODEL_PREPARATION_BASE/'vat_candidates_by_allele_gene.tsv.gz',sep='\t',dtype=str,keep_default_na=False)
    candidates.loc[candidates.uncertainty_flags.ne('')].to_csv(MODEL_PREPARATION_OUT/'LOFTEE_multitranscript_review_internal.tsv.gz',sep='\t',index=False)
    with __import__('gzip').open(MODEL_PREPARATION_OUT/'TP53_NF1_VAF_review_internal.csv.gz','wt') as f:
        first=True
        for chunk in pd.read_csv(MODEL_PREPARATION_BASE/'section7/sample_variant_audit.csv.gz',chunksize=20000,keep_default_na=False,low_memory=False):
            part=chunk.loc[chunk.gene.isin(['TP53','NF1'])]
            part.to_csv(f,index=False,header=first);first=False
    # Reason counts are candidate allele/gene counts, not carrier counts.
    import importlib.util
    spec=importlib.util.spec_from_file_location('cv',MODEL_PREPARATION_ROOT/'scripts/functions/clinvar.py');cv=importlib.util.module_from_spec(spec);spec.loader.exec_module(cv)
    refs=pd.read_csv(MODEL_PREPARATION_BASE/'section6/clinvar_region_records.tsv.gz',sep='\t',dtype=str,keep_default_na=False)
    evidence=pd.read_csv(MODEL_PREPARATION_BASE/'section6/clinvar_assertion_evidence.tsv.gz',sep='\t',dtype=str,keep_default_na=False)
    reject=pd.read_csv(MODEL_PREPARATION_BASE/'section6/missense_exclusion_audit.tsv.gz',sep='\t',dtype=str,keep_default_na=False)
    rg={k:d for k,d in refs.groupby('normalized_key')};eg={k:d for k,d in evidence.groupby(['gene_symbol','normalized_key'])}
    decisions={k:d.decision.tolist() for k,d in reject.groupby(['gene_symbol','normalized_key'])}
    qualified={(r.gene_symbol,r.normalized_key) for r in pd.read_csv(MODEL_PREPARATION_BASE/'section6/eligible_missense.tsv.gz',sep='\t').itertuples()}
    exclusion_rows=[]
    for r in candidates.loc[candidates.canonical_missense.eq('True')].itertuples():
        pair=(r.gene_symbol,r.normalized_key);ref=rg.get(r.normalized_key);ev=eg.get(pair)
        values=[] if ref is None else ref.sig.tolist()+ref.review.tolist()
        if ref is not None and (~ref.conflict.isin(['','.'])).any():values+=['conflicting classifications']
        if ev is not None:values+=ev.classification.tolist()
        block=cv.exclusions(values)
        if pair in qualified:reason='eligible_before_call_QC'
        elif block:reason=block
        elif ref is None or not ref.sig.map(cv.clean).isin(cv.PLP).any():reason='no_allowed_germline_P_LP'
        elif not (ref.sig.map(cv.clean).isin(cv.PLP)&ref.review.map(cv.clean).map(cv.STARS).fillna(0).ge(2)).any():reason='stars_below_2'
        elif any('inheritance' in x or 'no_relevant' in x for x in decisions.get(pair,[])):reason='disease_or_inheritance_not_suitable_or_unresolved'
        elif any('XML' in x for x in decisions.get(pair,[])):reason='historical_XML_VCF_mismatch'
        else:reason='AF_or_other_eligibility_unresolved'
        exclusion_rows.append({'gene':r.gene_symbol,'normalized_key':r.normalized_key,'primary_reason':reason})
    reasons=pd.DataFrame(exclusion_rows);reasons.to_csv(MODEL_PREPARATION_OUT/'ClinVar_exclusion_by_candidate_internal.tsv.gz',sep='\t',index=False)
    reason_counts=reasons.groupby('primary_reason').size().rename('allele_gene_candidates').reset_index()
    reason_counts.to_csv(MODEL_PREPARATION_OUT/'ClinVar_exclusion_reason_counts_internal.csv',index=False)
    public=reason_counts.copy();public['allele_gene_candidates']=public.allele_gene_candidates.map(lambda n:str(n) if n==0 or n>=20 else '<20')
    public.to_csv(MODEL_PREPARATION_ROOT/'results/table/43_rowlands_ClinVar_exclusion_summary.csv',index=False)
    # Freeze every model input + inherited PRS calculation; no result-dependent changes.
    paths=[MODEL_PREPARATION_OUT/'frozen_model_inputs.csv.gz',MODEL_PREPARATION_BASE/'section8/comparable_scores_case_control.csv.gz',MODEL_PREPARATION_BASE/'section8/comparable_scores_diagnosis_age.csv.gz',
        MODEL_PREPARATION_BASE/'section7/carrier_burden_case_control.csv.gz',MODEL_PREPARATION_ROOT/'data/16_joint_model_cv_oof_person_level.csv.gz',
        MODEL_PREPARATION_ROOT/'data/prs_scores.csv',MODEL_PREPARATION_ROOT/'scripts/03_common_variants_prs.py',MODEL_PREPARATION_ROOT/'scripts/90_legacy_ptv_baseline.py',MODEL_PREPARATION_ROOT/'scripts/90_legacy_ptv_baseline.py']
    # Preserve the authorization already recorded in the protocol across reruns.
    authorized = bool(cfg.get('execution_gate', {}).get('may_run_expanded_models', False))
    paths += sorted((MODEL_PREPARATION_ROOT / 'scripts').glob('[0-9]*.py'))
    paths += sorted((MODEL_PREPARATION_ROOT / 'scripts/functions').glob('*.py'))
    locked={'input_sha256':{str(p.relative_to(MODEL_PREPARATION_ROOT)):model_preparation_sha(p) for p in paths},'cohort_n':len(ids),'diagnosis_age_n':len(age),
        'CV_n_splits':5,'CV_n_repeats':10,'CV_random_state':20260818,'original_repeat1_folds_reproduced':True,
        'PTV_definition':'existing_any_retained_same_gene_transcript_PTV_and_HC_known_consistent_EUR_AF_lt1pct',
        'definition_status':'provisional_with_explicit_manual_review_lists_no_new_terminal_exon_exclusion',
        'PRS311':'existing persisted PRS311; existing CHEK2 linked-tag exclusion unchanged; no SNP selection',
        'saved_score_SD_ddof':0,'model_scaling_SD_ddof':1,'age_model_scaling':'within frozen cases as existing script17',
        'accepted_for_requested_comparison_models':authorized,'user_acceptance_status':'authorized_in_protocol' if authorized else 'pending','computational_checks_passed':True,'final_variant_definition_frozen':False}
    (MODEL_PREPARATION_OUT/'pre_model_input_lock_internal.json').write_text(json.dumps(locked,indent=2)+'\n')
    metadata={k:v for k,v in locked.items() if k!='input_sha256'}
    metadata.update(status='section9_pre_model_computational_acceptance_passed',all_requested_audit_tables_generated=True,
        fixed_phenotype_covariates_PRS_verified=True,review_lists_are_not_completed_biological_adjudications=True,models_run=False)
    (MODEL_PREPARATION_ROOT/'results/report/43_rowlands_pre_model_acceptance.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (MODEL_PREPARATION_ROOT/'results/report/43_rowlands_pre_model_audit_summary.md').write_text('# Section 9: Pre-model verification and acceptance\n\nAll specified checklists have been generated and saved in the controlled internal section9 directory: each gene PTV/missense/new carrier only because of missense; overlap of the new four genes with the current QC nine genes and the old nine genes; multi-gene carriers; six version score zero proportions, distributions and two SDs; extreme samples and candidate calls; ClinVar exclusion reasons; LOFTEE/terminal exons/multiple transcripts and TP53/NF1 VAF list.\n\nFixed case definitions, age cohorts of 141,204 and 7,996, covariate completeness, PRS311, and the original 10 × 50% score were calculated and checked. Input validation values \u200b\u200bare locked; user acceptance has not been completed and the model has not been run.\n\nThis is a computational check under the current tentative PTV definition and is not equivalent to biological manual adjudication or final variant set freezing. Terminal exon/NMD markers are position review agents and no automatic elimination rules will be added; the review list is still explicitly retained. These unresolved issues must be acknowledged during acceptance; if subsequent models follow the tentative definition, the results must be marked as exploratory.\n\nIn order to reproduce the existing No. 16/17 model, the modeling SD uses ddof=1, and the age model is standardized within a fixed subset of cases; the saved reference standard score still uses ddof=0 for the full case/control cohort. Cross-validation only estimates PRS and burden mean/SD at the training fold, and all versions reuse the same original fold score. Each actual model SD is recorded separately, and the original scale effect is reported simultaneously.\n\nDetailed comparisons of counts, distribution moments, extreme variations and exact SDs are kept in-house to avoid recovering sparse numbers from side-by-side version differences. Public ClinVar reason table does not list totals and suppresses <20 counts.')
    print('All section9 audit tables and fixed-cohort/PRS/fold checks passed; model inputs locked; existing protocol authorization preserved; provisional definition and manual lists explicit.')

# ============================================================================
# 5. sensitivities: prepare_sensitivity_scores
# ============================================================================
import gzip
import hashlib
import importlib.util
import json
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

SENSITIVITY_PREPARATION_ROOT=Path(__file__).resolve().parents[1]
SENSITIVITY_PREPARATION_BASE=SENSITIVITY_PREPARATION_ROOT/'data/29_rowlands_extension'
SENSITIVITY_PREPARATION_OUT=SENSITIVITY_PREPARATION_BASE/'section10'
DELIVERY=SENSITIVITY_PREPARATION_OUT/'internal_delivery'

def sensitivity_preparation_sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def prepare_sensitivity_scores():
    SENSITIVITY_PREPARATION_OUT.mkdir(parents=True,exist_ok=True);DELIVERY.mkdir(exist_ok=True)
    assert subprocess.run(['git','check-ignore','-q',str(SENSITIVITY_PREPARATION_OUT)],cwd=SENSITIVITY_PREPARATION_ROOT).returncode==0
    cfg=yaml.safe_load((SENSITIVITY_PREPARATION_ROOT/'config/rowlands_13gene_extension.yaml').read_text())
    acceptance=json.loads((SENSITIVITY_PREPARATION_ROOT/'results/report/43_rowlands_pre_model_acceptance.json').read_text())
    assert not cfg['expanded_models_run'] and not acceptance['models_run'], 'Sensitivity plan must precede C association results'
    plan={
        'analysis':'PTV plus qualifying P/LP missense; approximate Rowlands all-PV weights',
        'locked_before_C_association_results':True,'models_run_at_lock':False,
        'sensitivity_versions':{
            'C_without_STK11':{'panel':'C minus STK11','other_rules':'same as C'},
            'C_PTV_only':{'panel':'13 genes','variant_classes':['PTV'],'other_rules':'current C QC and Rowlands weights; not legacy A'},
            'C_TP53_NF1_strict_VAF':{'panel':'13 genes','heterozygous_VAF_TP53_NF1':[0.30,0.70],
                'heterozygous_VAF_other_genes':[0.30,0.80],'homozygous_ALT_VAF_min':0.80,'other_rules':'same as C'},
            'C_ClinVar_min1star':{'panel':'13 genes','missense_min_stars':1,'relevant_P_LP_RCV_min_stars':1,
                'exclusions_first':['low_penetrance','risk_factor_or_risk_allele','conflict_or_opposing_B_LB'],
                'require_exact_allowed_germline_P_LP':True,'disease_inheritance_review':'same as C',
                'ClinVar_release':cfg['clinvar']['release_date'],'PTV_and_AF_and_call_QC':'same as C'},
            'C_without_CDH1':{'panel':'C minus CDH1','conditional_on_sufficient_information':True,
                'availability_rule':'complete gene carrier and allele-call traceability; estimability assessed without a significance cutoff',
                'no_subtype_claim_without_subtype_data':True,'other_rules':'same as C'}},
        'case_control_model':'existing script16 PRS311 + burden + interaction + existing covariates; HC1; same samples',
        'CV':'existing paired 10 repeats x 5 folds; train-fold-only scaling; increment relative to PRS311 baseline',
        'age_model':'existing script17 case-only recorded diagnosis-age OLS with PRS311 + burden + PCs1-10; HC1',
        'report':'all estimable predefined versions and strata, estimates and CI, actual SD, raw burden distribution, carrier proportion',
        'selection':'no selection of final scheme or new strata using individual significant results',
        'PRS311':'existing score and CHEK2 linked-tag handling unchanged; no reselection',
        'disclosure':'exact comparison and participant/variant tables internal; public sparse cells and inferable complementary differences suppressed',
        'model_gate':'section9 user acceptance required; biological review status must remain explicit',
        'final_variant_definition_frozen':False}
    plan_path=SENSITIVITY_PREPARATION_OUT/'prespecified_sensitivity_plan.json'
    if plan_path.exists():
        previous=json.loads(plan_path.read_text())
        assert {k:v for k,v in previous.items() if k!='locked_at_UTC'}==plan, 'Previously locked sensitivity specification differs'
        plan=previous
    else:
        plan['locked_at_UTC']=datetime.now(timezone.utc).isoformat()
        plan_path.write_text(json.dumps(plan,indent=2,ensure_ascii=False)+'\n')
    spec=importlib.util.spec_from_file_location('cv',SENSITIVITY_PREPARATION_ROOT/'scripts/functions/clinvar.py')
    cv=importlib.util.module_from_spec(spec);spec.loader.exec_module(cv)
    # Review only the fixed reference release; no new allele or PRS selection.
    vat=pd.read_csv(SENSITIVITY_PREPARATION_BASE/'vat_candidates_by_allele_gene.tsv.gz',sep='\t',dtype=str,keep_default_na=False)
    tx=pd.read_csv(SENSITIVITY_PREPARATION_BASE/'vat_selected_transcript_rows.tsv.gz',sep='\t',dtype=str,keep_default_na=False)
    af=tx.groupby(['gene_symbol','normalized_key']).gvs_eur_af.apply(lambda x:pd.to_numeric(x,errors='coerce').dropna().unique().tolist()).to_dict()
    refs=pd.read_csv(SENSITIVITY_PREPARATION_BASE/'section6/clinvar_region_records.tsv.gz',sep='\t',dtype=str,keep_default_na=False)
    refgroups={k:d for k,d in refs.groupby('normalized_key')}
    candidates=[];decisions=[]
    for r in vat.itertuples():
        if r.canonical_missense!='True':continue
        pair=(r.gene_symbol,r.normalized_key);d=refgroups.get(r.normalized_key)
        values=[] if d is None else d.sig.tolist()+d.review.tolist()
        if d is not None and (~d.conflict.isin(['','.'])).any():values.append('conflicting classifications')
        reason=cv.exclusions(values)
        if not reason and not (len(af[pair])==1 and af[pair][0]<.01):reason='EUR_AF_missing_disagreeing_or_not_lt_0p01'
        allowed=None if d is None else d.loc[d.sig.map(cv.clean).isin(cv.PLP)&d.review.map(cv.clean).map(cv.STARS).fillna(0).ge(1)]
        if not reason and (allowed is None or allowed.empty):reason='not_allowed_at_least_1star_germline_P_LP'
        if reason:decisions.append(dict(gene=r.gene_symbol,normalized_key=r.normalized_key,keep=False,rejection_reason=reason))
        else:
            a=allowed.iloc[0];candidates.append((r.gene_symbol,r.normalized_key,a.variation_id,a.sig,a.review))
    print('Prespecification locked; reviewing historical ClinVar assertions for min1-star candidates.',flush=True)
    xml=cv.historical_xml(sorted({r[2] for r in candidates}),cfg['clinvar']['release_date'],hold_unavailable=True)
    eligible=set();evidence=[]
    for gene,key,vid,sig,review in candidates:
        if vid not in xml:
            decisions.append(dict(gene=gene,normalized_key=key,keep=False,rejection_reason='historical_VCV_unavailable_at_frozen_release_pending_review'))
            continue
        ok,reason,rcvs,rows=cv.review_variant(xml[vid],gene,sig,review,min_stars=1)
        decisions.append(dict(gene=gene,normalized_key=key,keep=ok,rejection_reason='' if ok else reason))
        if ok:eligible.add((gene,key))
        for row in rows:row.update(gene=gene,normalized_key=key,variation_id=vid);evidence.append(row)
    pd.DataFrame(decisions).to_csv(SENSITIVITY_PREPARATION_OUT/'min1star_variant_eligibility_audit.tsv.gz',sep='\t',index=False)
    pd.DataFrame(evidence).to_csv(SENSITIVITY_PREPARATION_OUT/'min1star_assertion_evidence.tsv.gz',sep='\t',index=False)
    main_eligible={(r.gene_symbol,r.normalized_key) for r in pd.read_csv(SENSITIVITY_PREPARATION_BASE/'section6/eligible_missense.tsv.gz',sep='\t').itertuples()}
    assert main_eligible<=eligible,'Min1-star must retain all main clinically eligible missense'
    main_records=pd.read_csv(SENSITIVITY_PREPARATION_BASE/'section6/eligible_missense.tsv.gz',sep='\t',dtype=str)
    for r in main_records.itertuples():
        assert cv.review_variant(xml[r.variation_id],r.gene_symbol,r.frozen_sig,r.frozen_review)[0], 'Default2-star review changed'
    current=pd.read_csv(SENSITIVITY_PREPARATION_BASE/'section7/unweighted_gene_carriers.csv.gz').set_index('sample_id')
    genes=cfg['panels']['main_13'];weights=pd.read_csv(SENSITIVITY_PREPARATION_ROOT/cfg['burden']['external_weight_table']).set_index('gene').ln_OR.loc[genes]
    min1={g:set() for g in genes}
    qc_names={'site_filter_not_PASS','genotype_filter_failed','not_complete_diploid_target_ALT_genotype',
        'GQ_missing_nonfinite_or_lt20','AD_missing_invalid_or_zero_total','VAF_outside_main_range'}
    decision_lookup={(r['gene'],r['normalized_key']):r['rejection_reason'] for r in decisions}
    dest=DELIVERY/'02_variant_sample_inclusion_exclusion_audit.csv.gz'
    with gzip.open(dest,'wt',newline='') as f:
        first=True
        for d in pd.read_csv(SENSITIVITY_PREPARATION_BASE/'section7/sample_variant_audit.csv.gz',chunksize=20000,keep_default_na=False,low_memory=False):
            keep=[];why=[]
            for r in d.to_dict('records'):
                pair=(r['gene'],r['normalized_key']);ptv=r['PTV/missense']=='PTV'
                if ptv:k=bool(int(r['keep_flag']));reason=r['rejection_reason']
                else:
                    reasons=[x for x in r['rejection_reason'].split(';') if x in qc_names]
                    if pair not in eligible:reasons.insert(0,decision_lookup.get(pair,'not_eligible_canonical_missense'))
                    k=not reasons;reason=';'.join(reasons)
                    if k:min1[r['gene']].add(int(r['sample_id']))
                keep.append(int(k));why.append(reason)
            d['keep_min1star']=keep;d['rejection_reason_min1star']=why
            assert not (d.keep_flag.astype(int).eq(1)&d.keep_min1star.eq(0)).any()
            d.to_csv(f,index=False,header=first);first=False
    one=current.copy()
    for g in genes:
        one[f'{g}_missense_main']=one.index.isin(min1[g]).astype('int8')
        one[f'{g}_main']=np.maximum(one[f'{g}_PTV_main'],one[f'{g}_missense_main'])
    one.to_csv(SENSITIVITY_PREPARATION_OUT/'min1star_unweighted_gene_carriers.csv.gz')
    # Derive scores from carriers, independently of the preceding eligibility stage.
    ptv=current.copy()
    for g in genes:ptv[f'{g}_main']=ptv[f'{g}_PTV_main']
    definitions={'C_PTV_only':(ptv,genes,'main'),'C_TP53_NF1_strict_VAF':(current,genes,'tp53_nf1_vaf_sensitive'),
        'C_ClinVar_min1star':(one,genes,'main'),'C_without_CDH1':(current,[g for g in genes if g!='CDH1'],'main')}
    sensitive=pd.DataFrame(index=current.index);sensitive.index.name='sample_id'
    for version,(frame,panel,rule) in definitions.items():
        matrix=frame[[f'{g}_{rule}' for g in panel]].to_numpy(dtype=int)
        raw=matrix@weights.loc[panel].to_numpy();sensitive[f'B_raw_{version}']=raw
        sensitive[f'carried_gene_count_{version}']=matrix.sum(axis=1)
        sd=raw.std(ddof=0);assert sd>0
        sensitive[f'B_standardized_{version}']=(raw-raw.mean())/sd
    main=pd.read_csv(SENSITIVITY_PREPARATION_BASE/'section8/comparable_scores_case_control.csv.gz').set_index('sample_id').loc[current.index]
    assert (sensitive.B_raw_C_ClinVar_min1star>=main.B_raw_C-1e-12).all()
    assert np.allclose(main.B_raw_C-sensitive.B_raw_C_without_CDH1,current.CDH1_main*weights.CDH1)
    existing=pd.read_csv(SENSITIVITY_PREPARATION_BASE/'section7/carrier_burden_case_control.csv.gz').set_index('sample_id').loc[current.index]
    assert np.allclose(sensitive.B_raw_C_TP53_NF1_strict_VAF,existing.B_raw_13gene_tp53_nf1_vaf_sensitive)
    sensitive.to_csv(SENSITIVITY_PREPARATION_OUT/'prespecified_sensitivity_scores_case_control.csv.gz')
    age=set(pd.read_csv(SENSITIVITY_PREPARATION_ROOT/cfg['cohorts']['diagnosis_age']['frozen_ids']).person_id)
    sensitive.loc[sorted(age)].to_csv(SENSITIVITY_PREPARATION_OUT/'prespecified_sensitivity_scores_diagnosis_age.csv.gz')
    shutil.copyfile(SENSITIVITY_PREPARATION_ROOT/cfg['burden']['external_weight_table'],DELIVERY/'01_13gene_OR_weights.csv')
    all_scores=main.join(sensitive);rows=[]
    for cohort,members in [('case_control',set(main.index)),('diagnosis_age',age)]:
        s=all_scores.loc[sorted(members)]
        for col in s.filter(regex='^B_raw_'):
            version=col.removeprefix('B_raw_');x=s[col];carrier=s[f'carried_gene_count_{version}']>0
            row=dict(cohort=cohort,version=version,n=len(x),carriers=int(carrier.sum()),carrier_fraction=float(carrier.mean()),
                zero_fraction=float(x.eq(0).mean()),raw_mean=float(x.mean()),SD_population_ddof0=float(x.std(ddof=0)),
                SD_sample_ddof1=float(x.std(ddof=1)),model_status='not_run_pending_section9_user_acceptance',
                OR_per_model_SD=None,OR_CI_low=None,OR_CI_high=None,CV_AUC=None,CV_AUC_increment_over_PRS311=None,
                age_beta_per_model_SD=None,age_CI_low=None,age_CI_high=None,age_R_squared=None,age_R_squared_increment=None)
            row.update({f'raw_p{int(q*100):02d}':float(x.quantile(q)) for q in [0,.01,.25,.5,.75,.95,.99,1]});rows.append(row)
    pd.DataFrame(rows).to_csv(DELIVERY/'03_A_to_C_carrier_and_model_comparison.csv',index=False)
    # Verify the persisted audit, not only the in-memory carrier calculations.
    persisted=pd.read_csv(dest,usecols=['sample_id','gene','PTV/missense','keep_flag','keep_min1star'])
    assert not (persisted.keep_flag.eq(1)&persisted.keep_min1star.eq(0)).any()
    for g in genes:
        observed=set(persisted.loc[persisted.gene.eq(g)&persisted['PTV/missense'].eq('missense')&persisted.keep_min1star.eq(1),'sample_id'])
        assert observed==set(one.index[one[f'{g}_missense_main'].eq(1)])
    saved=pd.read_csv(SENSITIVITY_PREPARATION_OUT/'prespecified_sensitivity_scores_case_control.csv.gz').set_index('sample_id')
    saved_age=pd.read_csv(SENSITIVITY_PREPARATION_OUT/'prespecified_sensitivity_scores_diagnosis_age.csv.gz').set_index('sample_id')
    pd.testing.assert_frame_equal(saved.loc[saved_age.index],saved_age)
    sources=[SENSITIVITY_PREPARATION_ROOT/'config/rowlands_13gene_extension.yaml',SENSITIVITY_PREPARATION_ROOT/'scripts/functions/clinvar.py',Path(__file__),
        SENSITIVITY_PREPARATION_ROOT/cfg['burden']['external_weight_table'],SENSITIVITY_PREPARATION_BASE/'section7/sample_variant_audit.csv.gz',
        SENSITIVITY_PREPARATION_BASE/'section7/unweighted_gene_carriers.csv.gz',SENSITIVITY_PREPARATION_BASE/'section8/comparable_scores_case_control.csv.gz',plan_path]
    manifest={'status':'sensitivity_plan_locked_and_scores_generated_models_pending_acceptance','models_run':False,
        'final_variant_definition_frozen':False,'plan_sha256':sensitivity_preparation_sha(plan_path),'main_min_stars':2,'sensitivity_min_stars':1,
        'all_main_eligible_missense_retained_in_min1star':True,'main_audit_and_scores_not_overwritten':True,
        'all_main_eligible_default2star_reviews_unchanged':True,'persisted_audit_reconstructs_min1star_carriers':True,
        'persisted_age_scores_exact_fixed_subset':True,
        'CDH1_gene_and_call_information_available':True,'CDH1_model_estimability':'not_yet_assessed_no_association_results_viewed',
        'inputs':{str(p.relative_to(SENSITIVITY_PREPARATION_ROOT)):sensitivity_preparation_sha(p) for p in sources},
        'delivery':{p.name:sensitivity_preparation_sha(p) for p in DELIVERY.iterdir() if p.is_file()}}
    (SENSITIVITY_PREPARATION_OUT/'delivery_manifest_internal.json').write_text(json.dumps(manifest,indent=2)+'\n')
    public={k:v for k,v in manifest.items() if k not in {'inputs','delivery'}}
    (SENSITIVITY_PREPARATION_ROOT/'results/report/44_rowlands_section10_validation.json').write_text(json.dumps(public,indent=2)+'\n')
    print('Prespecified scores, three internal delivery tables and source checksums complete; no models fitted.',flush=True)

# Command-line selection of the computations defined above.
def main():
    from scripts.functions.cli import run_local_step
    computations = {
        'score': calculate_burden_scores,
        'validate-score': validate_burden_scores,
        'compare': calculate_comparison_scores,
        'model-inputs': prepare_model_inputs,
        'sensitivities': prepare_sensitivity_scores,
    }
    run_local_step('Gene-collapsed log-OR burden, normalization and comparison scores', computations, default_tasks=['score', 'validate-score', 'compare', 'model-inputs', 'sensitivities'])

if __name__ == '__main__':
    main()

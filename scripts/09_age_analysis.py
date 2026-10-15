"""Fit diagnosis-age OLS in the frozen case cohort and finalize existing model outputs.

Unweighted binomial GLM with HC1 covariance; age analysis uses HC1 OLS.
The existing score versions, cohort assertions and output filenames are retained.
"""

from __future__ import annotations
from pathlib import Path
import sys

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.functions.paths import PROJECT_ROOT, implementation_sources
from scripts.functions.cli import parse_tasks
from scripts.functions.io import sha256 as sha
from scripts.functions.stats import fit_binomial as logistic

TASKS = ['age-models']
if __name__ == '__main__':
    selected_tasks = parse_tasks(__doc__, TASKS, ['age-models'])


PROJECT_ROOT = PROJECT_ROOT

import ast

from concurrent.futures import ThreadPoolExecutor

import hashlib

import json


import time

import numpy as np

import pandas as pd

import statsmodels.api as sm

from scipy import stats

from sklearn.metrics import roc_auc_score

from sklearn.model_selection import RepeatedStratifiedKFold

from threadpoolctl import threadpool_limits

import yaml

ASSOCIATION_MODELS_ROOT = PROJECT_ROOT

BASE = ASSOCIATION_MODELS_ROOT / 'data/29_rowlands_extension'

OUT = BASE / 'section11_models'

PC = [f'PC{i}' for i in range(1, 11)]

Z = 1.959963984540054

# Inputs and outputs used by both model stages. Data locations are unchanged.
INPUT_FROZEN_COHORT = BASE / 'section9/frozen_model_inputs.csv.gz'
INPUT_COMPARISON_SCORES = BASE / 'section8/comparable_scores_case_control.csv.gz'
INPUT_SENSITIVITY_SCORES = BASE / 'section10/prespecified_sensitivity_scores_case_control.csv.gz'
INPUT_BASELINE_OOF = PROJECT_ROOT / 'data/16_joint_model_cv_oof_person_level.csv.gz'
OUTPUT_COEFFICIENTS = OUT / 'model_coefficients_internal.csv'
OUTPUT_CASE_STATE = OUT / 'case_model_stage_internal.json'
CV_FOLDS = 5
CV_REPEATS = 10
CV_RANDOM_SEED = 20260818
CASE_CONTROL_TERMS = [
    'intercept',
    'PRS311_z',
    'burden_z',
    'PRS311_z_x_burden_z',
    'age_at_last_obs',
    *PC,
    'sequencing-site dummies',
    'sequencer dummies',
]
AGE_TERMS = ['intercept', 'PRS311_z', 'burden_z', *PC]


def design(prs_z, burden_z, c, interaction=True):
    cols = [np.ones(len(prs_z)), prs_z]
    if burden_z is not None:
        cols.append(burden_z)
        if interaction:
            cols.append(prs_z * burden_z)
    return np.column_stack(cols + [c])


if __name__ == "__main__":

    # ========================================================================
    # 1. Load the same locked model inputs and verify stage 08 provenance
    # ========================================================================
    started = time.time()

    OUT.mkdir(parents=True, exist_ok=True)

    import subprocess

    assert (
        subprocess.run(
            ['git', 'check-ignore', '-q', str(OUT)], cwd=ASSOCIATION_MODELS_ROOT
        ).returncode
        == 0
    )

    acceptance = json.loads(
        (
            ASSOCIATION_MODELS_ROOT / 'results/report/43_rowlands_pre_model_acceptance.json'
        ).read_text()
    )

    assert acceptance[
        'accepted_for_requested_comparison_models'
    ], 'User model authorization required'

    locked = json.loads((BASE / 'section9/pre_model_input_lock_internal.json').read_text())

    from scripts.functions.provenance import verify_frozen_inputs

    verify_frozen_inputs(locked['input_sha256'])

    plan = BASE / 'section10/prespecified_sensitivity_plan.json'

    pre = json.loads(
        (
            ASSOCIATION_MODELS_ROOT / 'results/report/44_rowlands_section10_validation.json'
        ).read_text()
    )

    assert (
        sha(plan) == pre['plan_sha256']
        and json.loads(plan.read_text())['locked_before_C_association_results']
    )

    inputs = [
        INPUT_FROZEN_COHORT,
        INPUT_COMPARISON_SCORES,
        INPUT_SENSITIVITY_SCORES,
        BASE / 'section7/unweighted_gene_carriers.csv.gz',
        ASSOCIATION_MODELS_ROOT / 'data/16_joint_model_cv_oof_person_level.csv.gz',
        ASSOCIATION_MODELS_ROOT / 'data/prs_scores.csv',
        ASSOCIATION_MODELS_ROOT / 'scripts/03_common_variants_prs.py',
        ASSOCIATION_MODELS_ROOT / 'scripts/legacy/old_9gene_ptv_analysis.py',
        ASSOCIATION_MODELS_ROOT / 'scripts/legacy/old_9gene_ptv_analysis.py',
        ASSOCIATION_MODELS_ROOT / 'scripts/08_case_control_models.py',
        ASSOCIATION_MODELS_ROOT / 'scripts/09_age_analysis.py',
        plan,
    ]

    inputs += implementation_sources(ASSOCIATION_MODELS_ROOT)

    inputs += []

    hashes = {str(prs_z.relative_to(ASSOCIATION_MODELS_ROOT)): sha(prs_z) for prs_z in inputs}

    run_hash = hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest()

    (OUT / 'model_input_lock_internal.json').write_text(
        json.dumps(
            {'inputs': hashes, 'run_hash': run_hash, 'PRS_and_samples_unchanged': True}, indent=2
        )
        + '\n'
    )

    from scripts.functions.stats import delong_auc_variance, delong_test, auc_ci

    ns = {'delong_auc_variance': delong_auc_variance, 'delong_test': delong_test, 'auc_ci': auc_ci}

    model_cohort = pd.read_csv(INPUT_FROZEN_COHORT).set_index('person_id')

    baseline_oof = pd.read_csv(
        ASSOCIATION_MODELS_ROOT / 'data/16_joint_model_cv_oof_person_level.csv.gz'
    ).set_index('person_id')

    assert model_cohort.index.equals(baseline_oof.index)

    scores = (
        pd.read_csv(INPUT_COMPARISON_SCORES)
        .set_index('sample_id')
        .join(pd.read_csv(INPUT_SENSITIVITY_SCORES).set_index('sample_id'))
        .loc[model_cohort.index]
    )

    versions = [c.removeprefix('B_raw_') for c in scores.filter(regex='^B_raw_')]

    burden_matrix = scores[[f'B_raw_{v}' for v in versions]].to_numpy(dtype=float)

    case_status = model_cohort.is_case.to_numpy(dtype=int)

    prs_values = model_cohort.PRS311.to_numpy(dtype=float)

    covariate_matrix = pd.concat(
        [
            model_cohort[['age_at_last_obs'] + PC],
            pd.get_dummies(model_cohort.seq_site_id, prefix='site', drop_first=True).astype(float),
            pd.get_dummies(model_cohort.seq_sequencer, prefix='seq', drop_first=True).astype(float),
        ],
        axis=1,
    ).to_numpy(dtype=float)

    assert (
        len(model_cohort) == 141204
        and case_status.sum() == 8528
        and np.isfinite(covariate_matrix).all()
    )

    age_ids = pd.read_csv(BASE / 'diagnosis_age_ids.csv.gz').person_id

    age_case_mask = model_cohort.index.isin(age_ids)

    assert age_case_mask.sum() == 7996 and case_status[age_case_mask].sum() == 7996

    diagnosis_age = model_cohort.age_at_diagnosis.to_numpy(dtype=float)[age_case_mask]

    age_covariates = model_cohort[PC].to_numpy(dtype=float)[age_case_mask]

    scaling = []

    coefficient_rows = []

    age_stats = []

    diagnostics = []

    case_state = json.loads(OUTPUT_CASE_STATE.read_text())

    assert case_state['versions'] == versions, 'Score versions changed between model stages'

    assert case_state['hashes'] == hashes, 'Model inputs changed between model stages'

    case_diagnostics = {row['version']: row for row in case_state['diagnostics']}

    # ========================================================================
    # 2. Fit within-case age models and assemble coefficients in original order
    # ========================================================================
    with threadpool_limits(limits=1):
        prs_z = (prs_values - prs_values.mean()) / prs_values.std(ddof=1)
        age_prs_z = (prs_values[age_case_mask] - prs_values[age_case_mask].mean()) / prs_values[
            age_case_mask
        ].std(ddof=1)
        age_base = sm.OLS(diagnosis_age, design(age_prs_z, None, age_covariates, False)).fit(
            cov_type='HC1'
        )
        for j, v in enumerate(versions):
            raw_burden = burden_matrix[:, j]
            sd = raw_burden.std(ddof=1)
            burden_z = (raw_burden - raw_burden.mean()) / sd
            case_burden = raw_burden[age_case_mask]
            case_burden_sd = case_burden.std(ddof=1)
            age_burden_z = (case_burden - case_burden.mean()) / case_burden_sd
            age_model = sm.OLS(
                diagnosis_age, design(age_prs_z, age_burden_z, age_covariates, False)
            ).fit(cov_type='HC1')
            for term, k in [('PRS311_z', 1), ('burden_z', 2)]:
                beta, se = (age_model.params[k], age_model.bse[k])
                coefficient_rows.append(
                    dict(
                        cohort='diagnosis_age',
                        version=v,
                        term=term,
                        beta=beta,
                        se=se,
                        CI_beta_low=beta - Z * se,
                        CI_beta_high=beta + Z * se,
                        OR=np.nan,
                        CI_low=np.nan,
                        CI_high=np.nan,
                        p_value=age_model.pvalues[k],
                        beta_raw_burden_unit=beta / case_burden_sd if k == 2 else np.nan,
                    )
                )
            for label, xx, mm, ss in [
                ('case_control', raw_burden, raw_burden.mean(), sd),
                ('diagnosis_age', case_burden, case_burden.mean(), case_burden_sd),
            ]:
                scaling.append(
                    dict(
                        cohort=label,
                        version=v,
                        model_mean=mm,
                        model_SD_ddof1=ss,
                        reference_mean=raw_burden.mean(),
                        reference_SD_ddof0=raw_burden.std(ddof=0),
                    )
                )
            age_stats.append(
                dict(
                    version=v,
                    R_squared=age_model.rsquared,
                    adjusted_R_squared=age_model.rsquared_adj,
                    PRS_PC_R_squared=age_base.rsquared,
                    unique_burden_R_squared=age_model.rsquared - age_base.rsquared,
                )
            )
            diagnostics.append(dict(case_diagnostics[v], age_finite=True))
        age_coefficients = coefficient_rows
        coefficient_rows = []
        for version in versions:
            coefficient_rows.extend(
                (row for row in case_state['coefficients'] if row['version'] == version)
            )
            coefficient_rows.extend((row for row in age_coefficients if row['version'] == version))
        pd.DataFrame(coefficient_rows).to_csv(OUT / 'model_coefficients_internal.csv', index=False)
        pd.DataFrame(scaling).to_csv(OUT / 'actual_model_scaling_internal.csv', index=False)
        pd.DataFrame(age_stats).to_csv(OUT / 'age_explained_variance_internal.csv', index=False)

    # ========================================================================
    # 3. Complete comparison tables, validation and delivery metadata
    # ========================================================================
    summary = pd.read_csv(OUT / 'CV_AUC_summary_internal.csv', float_precision='round_trip')

    ds = pd.read_csv(OUT / 'CV_delta_summary_internal.csv', float_precision='round_trip')

    cf = pd.DataFrame(coefficient_rows)

    sc = pd.DataFrame(scaling).set_index(['cohort', 'version'])

    ar = pd.DataFrame(age_stats).set_index('version')

    table_path = BASE / 'section10/internal_delivery/03_A_to_C_carrier_and_model_comparison.csv'

    table = pd.read_csv(table_path)

    for ix, row in table.iterrows():
        v = row.version
        cc = cf.loc[cf.version.eq(v) & cf.cohort.eq(row.cohort)]
        burden_z = cc.loc[cc.term.eq('burden_z')].iloc[0]
        table.loc[ix, 'model_status'] = 'fitted_provisional_variant_definition_exploratory'
        table.loc[ix, 'actual_model_SD_ddof1'] = sc.loc[(row.cohort, v), 'model_SD_ddof1']
        table.loc[ix, 'beta_raw_burden_unit'] = burden_z.beta_raw_burden_unit
        if row.cohort == 'case_control':
            i = cc.loc[cc.term.eq('PRS311_z_x_burden_z')].iloc[0]
            for key, val in {
                'OR_per_model_SD': burden_z.OR,
                'OR_CI_low': burden_z.CI_low,
                'OR_CI_high': burden_z.CI_high,
                'interaction_OR': i.OR,
                'interaction_CI_low': i.CI_low,
                'interaction_CI_high': i.CI_high,
            }.items():
                table.loc[ix, key] = val
            a = summary.set_index('version').loc[v]
            dd = ds.loc[ds.version.eq(v) & ds.comparison.eq('joint_minus_PRS311_baseline')].iloc[0]
            table.loc[ix, 'CV_AUC'] = a.AUC_mean
            table.loc[ix, 'CV_AUC_increment_over_PRS311'] = dd.delta_AUC_mean
            table.loc[ix, 'CV_delta_repeat1_CI_low'] = dd.repeat1_CI_low
            table.loc[ix, 'CV_delta_repeat1_CI_high'] = dd.repeat1_CI_high
        else:
            for key, val in {
                'age_beta_per_model_SD': burden_z.beta,
                'age_CI_low': burden_z.CI_beta_low,
                'age_CI_high': burden_z.CI_beta_high,
                'age_R_squared': ar.loc[v, 'R_squared'],
                'age_R_squared_increment': ar.loc[v, 'unique_burden_R_squared'],
            }.items():
                table.loc[ix, key] = val

    table.to_csv(table_path, index=False)

    public = cf.drop(columns=['beta_raw_burden_unit']).copy()

    public.to_csv(
        ASSOCIATION_MODELS_ROOT / 'results/table/45_rowlands_fixed_cohort_coefficients.csv',
        index=False,
    )

    summary.to_csv(ASSOCIATION_MODELS_ROOT / 'results/table/45_rowlands_CV_AUC.csv', index=False)

    ds.to_csv(ASSOCIATION_MODELS_ROOT / 'results/table/45_rowlands_CV_delta_AUC.csv', index=False)

    pd.DataFrame(age_stats).to_csv(
        ASSOCIATION_MODELS_ROOT / 'results/table/45_rowlands_age_explained_variance.csv',
        index=False,
    )

    for prs_z, h in hashes.items():
        assert sha(ASSOCIATION_MODELS_ROOT / prs_z) == h, 'Model input changed during fit: ' + prs_z

    manifest_path = BASE / 'section10/delivery_manifest_internal.json'

    manifest = json.loads(manifest_path.read_text())

    manifest['pre_model_delivery_sha256'] = manifest['delivery'].copy()

    manifest['delivery'] = {
        prs_z.name: sha(prs_z) for prs_z in table_path.parent.iterdir() if prs_z.is_file()
    }

    manifest['models_run'] = True

    manifest['model_input_lock'] = (
        'data/29_rowlands_extension/section11_models/model_input_lock_internal.json'
    )

    manifest['status'] = 'comparison_models_completed_provisional_variant_definition'

    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')

    validation = {
        'status': 'fixed_cohort_models_and_all_prespecified_sensitivities_completed',
        'models_run': True,
        'versions': versions,
        'cohort_n': len(case_status),
        'case_n': int(case_status.sum()),
        'diagnosis_age_n': int(age_case_mask.sum()),
        'CV_repeats': 10,
        'CV_folds': 5,
        'original_PRS311_repeat1_OOF_reproduced': True,
        'input_hashes_unchanged': True,
        'training_fold_scaling_only': True,
        'case_model': 'GLM Binomial unweighted HC1; PRS + burden + interaction + existing covariates',
        'burden_OR_interpretation': 'per model SD conditional on PRS at its mean; not a marginal carrier OR',
        'age_model': 'within fixed cases recorded diagnosis-age OLS HC1; PRS + burden + PCs1-10; no interaction',
        'CI_note': 'repeat1 paired DeLong CI conditional on fitted cross-validation predictions; repeat ranges and mean CI bounds descriptive, not independent repetitions or training-uncertainty CI',
        'variant_definition': 'provisional; manual biological adjudication still pending',
        'final_variant_definition_frozen': False,
        'elapsed_seconds': time.time() - started,
        'fit_diagnostics': diagnostics,
    }

    (ASSOCIATION_MODELS_ROOT / 'results/report/45_rowlands_models_validation.json').write_text(
        json.dumps(validation, indent=2) + '\n'
    )

    cfgpath = ASSOCIATION_MODELS_ROOT / 'config/rowlands_13gene_extension.yaml'

    cfg = yaml.safe_load(cfgpath.read_text())

    cfg['expanded_models_run'] = True

    cfg['status'] = (
        'models_and_prespecified_sensitivities_completed_provisional_definition_pending_final_adjudication'
    )

    cfg['internal_delivery']['model_fields_status'] = 'completed_provisional_definition'

    cfgpath.write_text(yaml.safe_dump(cfg, sort_keys=False, allow_unicode=True))

    print(
        'Full models, paired 10x5 CV for all 10 versions, sensitivity estimates and internal model worksheet complete.',
        flush=True,
    )

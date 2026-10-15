"""Fit fixed-cohort case-control models and paired 10x5 cross-validation.

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

TASKS = ['models']
if __name__ == '__main__':
    selected_tasks = parse_tasks(__doc__, TASKS, ['models'])


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


def fold_job(tr, te, j, pt, pe):
    if j < 0:
        bt = be = None
    else:
        mu = burden_matrix[tr, j].mean()
        sd = burden_matrix[tr, j].std(ddof=1)
        bt = (burden_matrix[tr, j] - mu) / sd
        be = (burden_matrix[te, j] - mu) / sd
    case_model = logistic(case_status[tr], design(pt, bt, covariate_matrix[tr]), robust=False)
    return case_model.predict(design(pe, be, covariate_matrix[te]))


if __name__ == "__main__":

    # ========================================================================
    # 1. Verify frozen inputs and the prespecified sensitivity plan
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

    # ========================================================================
    # 2. Fit the joint case-control models and evaluate paired CV
    # ========================================================================
    with threadpool_limits(limits=1):
        prs_z = (prs_values - prs_values.mean()) / prs_values.std(ddof=1)
        for j, v in enumerate(versions):
            raw_burden = burden_matrix[:, j]
            sd = raw_burden.std(ddof=1)
            burden_z = (raw_burden - raw_burden.mean()) / sd
            case_design = design(prs_z, burden_z, covariate_matrix)
            case_model = sm.GLM(case_status, case_design, family=sm.families.Binomial()).fit(
                cov_type='HC1'
            )
            assert (
                case_model.converged
                and np.isfinite(case_model.params).all()
                and np.isfinite(case_model.bse).all()
            )
            for term, k in [('PRS311_z', 1), ('burden_z', 2), ('PRS311_z_x_burden_z', 3)]:
                beta, se = (case_model.params[k], case_model.bse[k])
                coefficient_rows.append(
                    dict(
                        cohort='case_control',
                        version=v,
                        term=term,
                        beta=beta,
                        se=se,
                        CI_beta_low=beta - Z * se,
                        CI_beta_high=beta + Z * se,
                        OR=np.exp(beta),
                        CI_low=np.exp(beta - Z * se),
                        CI_high=np.exp(beta + Z * se),
                        p_value=case_model.pvalues[k],
                        beta_raw_burden_unit=beta / sd if k in [2, 3] else np.nan,
                    )
                )
            diagnostics.append(
                dict(
                    version=v,
                    case_converged=bool(case_model.converged),
                    case_iterations=case_model.fit_history['iteration'],
                )
            )
        print(
            'All fixed-cohort full case and age models fitted; starting paired repeated CV.',
            flush=True,
        )
        strata = (
            case_status.astype(str)
            + np.array(['_'] * len(case_status))
            + baseline_oof.carrier_main_10gene.to_numpy(dtype=int).astype(str)
        )
        splits = list(
            RepeatedStratifiedKFold(
                n_splits=CV_FOLDS, n_repeats=CV_REPEATS, random_state=CV_RANDOM_SEED
            ).split(np.arange(len(case_status)), strata)
        )
        f = np.full(len(case_status), -1)
        for k, (_, te) in enumerate(splits[:5]):
            f[te] = k
        assert np.array_equal(f, baseline_oof.fold_r1)
        auc_rows = []
        delta_rows = []
        fold_scaling = []
        model_names = ['PRS311_baseline'] + versions
        with ThreadPoolExecutor(max_workers=3) as pool:
            for r in range(10):
                cp = OUT / f'cv_repeat_{r + 1:02d}.npz'
                if cp.exists():
                    saved = np.load(cp)
                    assert str(saved['run_hash']) == run_hash
                    pred = saved['predictions']
                else:
                    pred = np.full((len(model_names), len(case_status)), np.nan)
                    for k in range(5):
                        tr, te = splits[r * 5 + k]
                        pm = prs_values[tr].mean()
                        ps = prs_values[tr].std(ddof=1)
                        pt = (prs_values[tr] - pm) / ps
                        pe = (prs_values[te] - pm) / ps
                        jobs = [
                            pool.submit(fold_job, tr, te, j, pt, pe)
                            for j in range(-1, len(versions))
                        ]
                        for j, job in enumerate(jobs):
                            pred[j, te] = job.result()
                        print(f'CV repeat {r + 1}/10 fold {k + 1}/5 complete.', flush=True)
                    assert np.isfinite(pred).all() and ((pred >= 0) & (pred <= 1)).all()
                    np.savez_compressed(cp, run_hash=run_hash, predictions=pred)
                for k in range(5):
                    tr, _ = splits[r * 5 + k]
                    for j, v in enumerate(versions):
                        fold_scaling.append(
                            dict(
                                repeat=r + 1,
                                fold=k + 1,
                                version=v,
                                mean=burden_matrix[tr, j].mean(),
                                SD_ddof1=burden_matrix[tr, j].std(ddof=1),
                                PRS_mean=prs_values[tr].mean(),
                                PRS_SD_ddof1=prs_values[tr].std(ddof=1),
                            )
                        )
                if r == 0:
                    assert np.allclose(
                        pred[0], baseline_oof.oof_r1_prs_only.to_numpy(), atol=1e-08, rtol=1e-08
                    ), 'Existing PRS OOF baseline not reproduced'
                for j, name in enumerate(model_names):
                    lo, hi = auc_ci(case_status, pred[j])
                    auc_rows.append(
                        dict(
                            repeat=r + 1,
                            version=name,
                            AUC=roc_auc_score(case_status, pred[j]),
                            CI_low=lo,
                            CI_high=hi,
                        )
                    )
                    if j:
                        res = delong_test(case_status, pred[j], pred[0])
                        delta_rows.append(
                            dict(
                                repeat=r + 1,
                                version=name,
                                comparison='joint_minus_PRS311_baseline',
                                **res,
                            )
                        )
                pairs = [
                    ('A_prime', 'A'),
                    ('A_prime_QC', 'A_prime'),
                    ('B', 'A_prime_QC'),
                    ('C', 'B'),
                    ('C_without_STK11', 'C'),
                    ('C_PTV_only', 'C'),
                    ('C_TP53_NF1_strict_VAF', 'C'),
                    ('C_ClinVar_min1star', 'C'),
                    ('C_without_CDH1', 'C'),
                ]
                for a, bn in pairs:
                    res = delong_test(
                        case_status, pred[model_names.index(a)], pred[model_names.index(bn)]
                    )
                    delta_rows.append(
                        dict(repeat=r + 1, version=a, comparison=f'{a}_minus_{bn}', **res)
                    )
        pd.DataFrame(fold_scaling).to_csv(
            OUT / 'CV_training_fold_scaling_internal.csv', index=False
        )

    # ========================================================================
    # 3. Save CV summaries and lossless aggregate state for diagnosis-age analysis
    # ========================================================================
    auc = pd.DataFrame(auc_rows)

    delta = pd.DataFrame(delta_rows)

    auc.to_csv(OUT / 'CV_AUC_by_repeat_internal.csv', index=False)

    delta.to_csv(OUT / 'CV_paired_deltas_by_repeat_internal.csv', index=False)

    summary = (
        auc.groupby('version', sort=False)
        .AUC.agg(AUC_mean='mean', AUC_min='min', AUC_max='max')
        .reset_index()
    )

    summary = summary.merge(
        auc.loc[auc.repeat.eq(1), ['version', 'AUC', 'CI_low', 'CI_high']].rename(
            columns={
                'AUC': 'AUC_repeat1',
                'CI_low': 'AUC_repeat1_CI_low',
                'CI_high': 'AUC_repeat1_CI_high',
            }
        ),
        on='version',
    )

    ds = (
        delta.groupby(['version', 'comparison'], sort=False)
        .agg(
            delta_AUC_mean=('delta_auc', 'mean'),
            delta_AUC_min=('delta_auc', 'min'),
            delta_AUC_max=('delta_auc', 'max'),
            CI_low_mean=('ci_low', 'mean'),
            CI_high_mean=('ci_high', 'mean'),
        )
        .reset_index()
    )

    ds = ds.merge(
        delta.loc[
            delta.repeat.eq(1),
            ['version', 'comparison', 'delta_auc', 'ci_low', 'ci_high', 'p_value'],
        ].rename(
            columns={
                'delta_auc': 'delta_AUC_repeat1',
                'ci_low': 'repeat1_CI_low',
                'ci_high': 'repeat1_CI_high',
                'p_value': 'repeat1_p',
            }
        ),
        on=['version', 'comparison'],
    )

    summary.to_csv(OUT / 'CV_AUC_summary_internal.csv', index=False)

    ds.to_csv(OUT / 'CV_delta_summary_internal.csv', index=False)

    case_state = {
        'coefficients': coefficient_rows,
        'diagnostics': diagnostics,
        'versions': versions,
        'hashes': hashes,
        'run_hash': run_hash,
    }

    OUTPUT_CASE_STATE.write_text(
        json.dumps(
            case_state,
            indent=2,
            default=lambda value: value.item() if hasattr(value, 'item') else str(value),
        )
        + '\n'
    )

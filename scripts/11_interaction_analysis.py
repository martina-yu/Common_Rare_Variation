"""Post-hoc interaction robustness and gene-specific interaction analyses.

These analyses remain exploratory. Numbering does not make them prespecified.
Inputs: frozen case-control cohort, C scores, gene carriers and baseline coefficients.
Outputs retain the existing results/table/59_interaction_audit_* filenames.
"""

from pathlib import Path
import sys

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.functions.paths import PROJECT_ROOT
from scripts.functions.cli import parse_tasks

if __name__ == '__main__':
    selected_tasks = parse_tasks(__doc__, ['audit'], [])

import os

os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

import hashlib, json, time, warnings


import numpy as np

import pandas as pd

import statsmodels.api as sm

from patsy import dmatrix

from scipy.stats import chi2

from threadpoolctl import threadpool_limits


def fit(name, x, mask=None, j=3, save=True):
    mask = np.ones(len(y), bool) if mask is None else mask
    print('Fitting ' + name, flush=True)
    f = sm.GLM(y[mask], x[mask], family=sm.families.Binomial()).fit(cov_type='HC1', maxiter=100)
    assert f.converged and np.isfinite(f.params).all() and np.isfinite(f.bse).all()
    if save:
        beta = float(f.params[j])
        se = float(f.bse[j])
        rows.append(
            dict(
                model=name,
                term='PRS_interaction',
                beta=beta,
                se_HC1=se,
                OR=np.exp(beta),
                CI_low=np.exp(beta - 1.95996398454 * se),
                CI_high=np.exp(beta + 1.95996398454 * se),
                p_value=float(f.pvalues[j]),
                converged=True,
                scaling='full_frozen_cohort_ddof1',
                n='withheld_for_joint_suppression',
                case_n='withheld_for_joint_suppression',
            )
        )
        pd.DataFrame(rows).to_csv(OUT / '59_interaction_audit_models.csv', index=False)
    return f


def wald(name, f, indices, R=None):
    if R is None:
        R = np.eye(len(f.params))[indices]
    v = R @ f.params
    V = R @ f.cov_params() @ R.T
    stat = float(v @ np.linalg.pinv(V) @ v)
    df = int(np.linalg.matrix_rank(V))
    tests.append(
        dict(
            test=name,
            method='HC1 robust Wald',
            statistic=stat,
            df=df,
            p_value=float(chi2.sf(stat, df)),
        )
    )
    pd.DataFrame(tests).to_csv(OUT / '59_interaction_audit_tests.csv', index=False)


def spline(a, df=4):
    return np.asarray(dmatrix(f'bs(x, df={df}, degree=3, include_intercept=False)-1', {'x': a}))


# ============================================================================
# 1. Load frozen inputs, reproduce baseline, and fit robustness models
# ============================================================================
if __name__ == "__main__":

    ROOT = PROJECT_ROOT

    BASE = ROOT / 'data/29_rowlands_extension'

    OUT = ROOT / 'results/table'

    REPORT = ROOT / 'results/report'

    INPUTS = [
        BASE / 'section9/frozen_model_inputs.csv.gz',
        BASE / 'section8/comparable_scores_case_control.csv.gz',
        BASE / 'section7/unweighted_gene_carriers.csv.gz',
        ROOT / 'data/04_phenotype_diagnostics_person_level.csv.gz',
        BASE / 'section11_models/model_coefficients_internal.csv',
    ]

    hashes = {str(x.relative_to(ROOT)): hashlib.sha256(x.read_bytes()).hexdigest() for x in INPUTS}

    d = pd.read_csv(INPUTS[0]).set_index('person_id')

    s = pd.read_csv(INPUTS[1]).set_index('sample_id').loc[d.index]

    g = pd.read_csv(INPUTS[2]).set_index('sample_id').loc[d.index]

    e = pd.read_csv(INPUTS[3]).set_index('person_id').loc[d.index]

    y = d.is_case.to_numpy(float)

    raw = s.B_raw_C.to_numpy(float)

    p = (d.PRS311.to_numpy() - d.PRS311.mean()) / d.PRS311.std(ddof=1)

    b = (raw - raw.mean()) / raw.std(ddof=1)

    carrier = (raw > 0).astype(float)

    pc = d[[f'PC{i}' for i in range(1, 11)]].to_numpy(float)

    technical = pd.concat(
        [
            pd.get_dummies(d.seq_site_id, drop_first=True),
            pd.get_dummies(d.seq_sequencer, drop_first=True),
        ],
        axis=1,
    ).to_numpy(float)

    age = d.age_at_last_obs.to_numpy(float)

    cov = np.column_stack([age, pc, technical])

    X = np.column_stack([np.ones(len(y)), p, b, p * b, cov])

    assert len(y) == 141204 and y.sum() == 8528 and np.isfinite(X).all()

    rows = []

    tests = []

    notes = []

    with threadpool_limits(limits=1):
        f = fit('C_baseline_reproduction', X)
        ref = pd.read_csv(INPUTS[4])
        r = ref[
            (ref.cohort == 'case_control')
            & (ref.version == 'C')
            & (ref.term == 'PRS311_z_x_burden_z')
        ].iloc[0]
        assert abs(f.params[3] - r.beta) < 1e-10 and abs(f.pvalues[3] - r.p_value) < 1e-10
        mu = f.fittedvalues
        w = mu * (1 - mu)
        bread = np.linalg.inv(X.T * w @ X)
        h = w * np.einsum('ij,jk,ik->i', X, bread, X)
        influence = (y - mu) * (X @ bread[:, 3]) / (1 - h) / f.bse[3]
        keep = np.ones(len(y), bool)
        keep[np.argsort(np.abs(influence))[-100:]] = False
        fit('exclude_top100_interaction_influence', X, keep)
        highlever = np.ones(len(y), bool)
        highlever[np.argsort(h)[-100:]] = False
        fit('exclude_top100_design_leverage', X, highlever)
        mask = s.carried_gene_count_C.to_numpy() <= 1
        if (~mask).sum() >= 20:
            fit('exclude_multigene_carriers', X, mask)
        threshold = float(np.quantile(raw[carrier == 1], 0.95))
        mask = raw <= threshold
        if (~mask).sum() >= 20:
            fit('exclude_above_carrier_burden_95th_percentile', X, mask)
        binary = np.column_stack([np.ones(len(y)), p, carrier, p * carrier, cov])
        fit('binary_C_carrier', binary)
        sp = spline(p)
        x = np.column_stack([np.ones(len(y)), sp, b, sp * b[:, None], cov])
        sf = fit('PRS_spline_burden_interaction', x, save=False)
        wald('PRS_spline_all_interactions_zero', sf, list(range(6, 10)))
        from scipy.linalg import qr

        resid = (
            sp
            - np.column_stack([np.ones(len(p)), p])
            @ np.linalg.lstsq(np.column_stack([np.ones(len(p)), p]), sp, rcond=None)[0]
        )
        q, rr, piv = qr(resid, mode='economic', pivoting=True)
        nl = resid[:, piv[:3]]
        xn = np.column_stack([X, nl, nl * b[:, None]])
        nf = fit('linear_interaction_with_nonlinear_PRS_terms', xn)
        wald(
            'nonlinear_PRS_main_and_interaction_terms_zero',
            nf,
            list(range(X.shape[1], X.shape[1] + 6)),
        )
        wald(
            'nonlinear_PRS_interaction_terms_zero', nf, list(range(X.shape[1] + 3, X.shape[1] + 6))
        )
        xa = np.column_stack([np.ones(len(y)), p, b, p * b, spline(age), pc, technical])
        fit('nonlinear_last_observation_age', xa)
        enroll = e.age_at_enrollment.to_numpy(float)
        valid = np.isfinite(enroll)
        xe = np.column_stack([np.ones(len(y)), p, b, p * b, np.nan_to_num(enroll), pc, technical])
        fit('enrollment_age_replaces_last_observation_age', xe, valid)
        window = valid & (enroll >= 50) & (enroll < 75)
        fit('enrollment_age_50_to74_baseline_covariates', X, window)
        xe2 = np.column_stack(
            [np.ones(len(y)), p, b, p * b, spline(np.nan_to_num(enroll)), pc, technical]
        )
        fit('enrollment_age_50_to74_nonlinear_enrollment_age', xe2, window)
        span = e.ehr_span_days.to_numpy(float)
        dates = e.n_distinct_ehr_dates.to_numpy(float)
        ev = np.isfinite(span) & np.isfinite(dates) & (span >= 0) & (dates >= 0)
        z = np.column_stack(
            [
                np.log1p(np.maximum(np.nan_to_num(span), 0)),
                np.log1p(np.maximum(np.nan_to_num(dates), 0)),
            ]
        )
        fit('EHR_complete_metadata_baseline', X, ev)
        fit('EHR_log_span_and_date_count_adjustment', np.column_stack([xa, z]), ev)
        coverage = ev & (span >= 365) & (dates >= 5)
        fit('EHR_span_ge365days_dates_ge5_baseline', X, coverage)
        fit(
            'EHR_span_ge365days_dates_ge5_age_spline_EHR_adjusted',
            np.column_stack([xa, z]),
            coverage,
        )
        genes = [
            c[:-5]
            for c in g.columns
            if c.endswith('_main')
            and (not c.endswith('_PTV_main'))
            and (not c.endswith('_missense_main'))
        ]
        gc = g[[gene + '_main' for gene in genes]].to_numpy(float)
        assert np.array_equal(gc.sum(axis=1), s.carried_gene_count_C.to_numpy())
        dense = [
            i for i in range(len(genes)) if gc[y == 1, i].sum() >= 20 and gc[y == 0, i].sum() >= 20
        ]
        sparse = [i for i in range(len(genes)) if i not in dense]
        G = gc[:, dense]
        names = [genes[i] for i in dense]
        if sparse:
            pooled = (gc[:, sparse].sum(axis=1) > 0).astype(float)
            if pooled[y == 1].sum() >= 20 and pooled[y == 0].sum() >= 20:
                G = np.column_stack([G, pooled])
                names.append('Other_pooled_panel_genes')
            else:
                notes.append(
                    'Sparse pooled category lacks 20 cases/controls; its interaction not estimated.'
                )
        kg = G.shape[1]
        xg = np.column_stack([np.ones(len(y)), p, G, p[:, None] * G, cov])
        gf = fit('joint_gene_specific_interactions', xg, save=False)
        inter = list(range(2 + kg, 2 + 2 * kg))
        wald('joint_gene_PRS_interactions_all_zero', gf, inter)
        R = np.zeros((kg - 1, len(gf.params)))
        for i in range(1, kg):
            R[i - 1, inter[i]] = 1
            R[i - 1, inter[0]] = -1
        wald('equality_of_gene_specific_PRS_interactions', gf, inter, R)
        grows = []
        for gene, j in zip(names, inter):
            be = float(gf.params[j])
            se = float(gf.bse[j])
            grows.append(
                dict(
                    gene=gene,
                    beta=be,
                    se_HC1=se,
                    interaction_OR=np.exp(be),
                    CI_low=np.exp(be - 1.95996398454 * se),
                    CI_high=np.exp(be + 1.95996398454 * se),
                    p_value=float(gf.pvalues[j]),
                    counts='not_exported_joint_suppression',
                )
            )
        from statsmodels.stats.multitest import multipletests

        adj = multipletests([a['p_value'] for a in grows], method='fdr_bh')[1]
        for a, qv in zip(grows, adj):
            a['BH_q_value'] = float(qv)
        pd.DataFrame(grows).to_csv(OUT / '59_interaction_audit_gene_interactions.csv', index=False)
        for i in dense:
            fit('exclude_' + genes[i] + '_carriers', X, gc[:, i] == 0)

    meta = dict(
        mode='validate',
        verification_status='VERIFIED_baseline_reproduction_only',
        sensitivity_status='new_exploratory_fits',
        input_sha256=hashes,
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        cohort_n=141204,
        case_n=8528,
        control_n=132676,
        burden_version='C',
        PRS_mean=float(d.PRS311.mean()),
        PRS_SD_ddof1=float(d.PRS311.std(ddof=1)),
        burden_mean=float(raw.mean()),
        burden_SD_ddof1=float(raw.std(ddof=1)),
        notes=notes,
        output_policy='No participant rows; no stratified counts; sparse gene estimates pooled; sample sizes of deletion models withheld to avoid reconstruction',
        influence_method='One-step interaction DFBETA using MLE information, normalized by HC1 SE; delete top100 absolute values and refit; design leverage top100 separately',
        exploratory=True,
    )

    (REPORT / '59_interaction_audit_metadata.json').write_text(json.dumps(meta, indent=2) + '\n')

    print('Completed aggregate audit; baseline reproduced exactly.', flush=True)

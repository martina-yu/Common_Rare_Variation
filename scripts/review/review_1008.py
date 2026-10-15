"""Reproduce review analyses from frozen study inputs; aggregate outputs only.

Run --task interaction, dates, cox, family, plots, or all. Unspecified reviewer
details are recorded as assumptions, never silently treated as exact matches.
"""

from pathlib import Path
import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.functions.paths import PROJECT_ROOT

import os

os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
import hashlib
import json
import warnings
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.duration.hazard_regression import PHReg
from statsmodels.stats.contingency_tables import StratifiedTable
from scipy.stats import norm, chi2
from patsy import dmatrix
from threadpoolctl import threadpool_limits

ROOT = PROJECT_ROOT
sys.path.insert(0, str(ROOT))
BASE = ROOT / 'data/29_rowlands_extension'
INTERNAL = ROOT / 'data/62_review_1008'
TABLE = ROOT / 'results/table'
REPORT = ROOT / 'results/report'
FIGURE = ROOT / 'results/figures'
PC = [f'PC{i}' for i in range(1, 11)]
Z95 = norm.ppf(0.975)
INPUTS = {
    'frozen': BASE / 'section9/frozen_model_inputs.csv.gz',
    'scores': BASE / 'section8/comparable_scores_case_control.csv.gz',
    'genes': BASE / 'section7/unweighted_gene_carriers.csv.gz',
    'phenotype': ROOT / 'data/04_phenotype_diagnostics_person_level.csv.gz',
    'sensitivity': BASE / 'section10/prespecified_sensitivity_scores_case_control.csv.gz',
    'clinical': ROOT / 'data/09_extended_clinical_person_level.csv.gz',
    'age_ids': BASE / 'diagnosis_age_ids.csv.gz',
    'weights': TABLE / '29_rowlands_13gene_external_weights.csv',
    'config': ROOT / 'config/rowlands_13gene_extension.yaml',
}


def write_table(name, rows):
    frame = rows if isinstance(rows, pd.DataFrame) else pd.DataFrame(rows)
    frame.to_csv(TABLE / f'62_review_1008_{name}.csv', index=False)


def save_json(name, obj):
    (REPORT / f'62_review_1008_{name}.json').write_text(
        json.dumps(obj, indent=2, default=lambda x: x.item() if hasattr(x, 'item') else str(x))
        + '\n'
    )


def estimate(beta, se, exponentiate=True):
    trans = np.exp if exponentiate else lambda x: x
    return dict(
        beta=float(beta),
        se=float(se),
        estimate=float(trans(beta)),
        CI_low=float(trans(beta - Z95 * se)),
        CI_high=float(trans(beta + Z95 * se)),
        p_value=float(2 * norm.sf(abs(beta / se))),
    )


def ratio_estimate(beta_c, se_c, beta_n, se_n):
    return estimate(beta_c - beta_n, np.hypot(se_c, se_n))


def disclosed_size(count, total=None):
    if 0 < count <= 20:
        return '≤20'
    if total is not None and 0 < total - count <= 20:
        return 'withheld_complementary_disclosure'
    return int(count)


def zscore(x):
    x = np.asarray(x, float)
    return (x - x.mean()) / x.std(ddof=1)


def load():
    d = pd.read_csv(INPUTS['frozen']).set_index('person_id')
    assert d.index.is_unique
    s = pd.read_csv(INPUTS['scores']).set_index('sample_id').loc[d.index]
    g = pd.read_csv(INPUTS['genes']).set_index('sample_id').loc[d.index]
    e = pd.read_csv(INPUTS['phenotype']).set_index('person_id').loc[d.index]
    sens = pd.read_csv(INPUTS['sensitivity']).set_index('sample_id').loc[d.index]
    technical = pd.concat(
        [
            pd.get_dummies(d.seq_site_id, prefix='site', drop_first=True),
            pd.get_dummies(d.seq_sequencer, prefix='seq', drop_first=True),
        ],
        axis=1,
    ).to_numpy(float)
    raw = s.B_raw_C.to_numpy(float)
    y = d.is_case.to_numpy(int)
    p, b = zscore(d.PRS311), zscore(raw)
    carrier = raw > 0
    brca = (g.BRCA1_main | g.BRCA2_main).to_numpy(bool)
    cov = np.column_stack([d.age_at_last_obs.to_numpy(float), d[PC].to_numpy(float), technical])
    assert np.isfinite(cov).all() and np.isfinite(p).all()
    provenance = dict(
        input_sha256={k: hashlib.sha256(v.read_bytes()).hexdigest() for k, v in INPUTS.items()},
        cohort_n=len(d),
        case_n=int(y.sum()),
        control_n=int((y == 0).sum()),
        PRS_mean=float(d.PRS311.mean()),
        PRS_SD_ddof1=float(d.PRS311.std(ddof=1)),
        burden_mean=float(raw.mean()),
        burden_SD_ddof1=float(raw.std(ddof=1)),
        definition='C: thirteen-gene PTV plus eligible >=2-star ClinVar P/LP germline missense; manual adjudication pending',
        source_reference='ca6bfba; reviewer repository not accessible at implementation time',
        release=os.environ.get('WORKSPACE_CDR', 'not_available').split('.')[-1],
        uncertainty='HC1 logistic and OLS; model-based Breslow Cox; source Cox tie/SE choices unavailable',
        disclosure='No participant rows or stratified counts in published outputs; counts 1-20 suppressed',
    )
    save_json('provenance', provenance)
    return d, s, g, e, sens, y, p, b, raw, carrier, brca, cov, technical


def fit_logistic(y, X, mask):
    x = X[mask]
    # Drop columns constant within a stratum, retaining the intercept.
    keep = np.r_[True, np.ptp(x[:, 1:], axis=0) > 1e-12]
    x = x[:, keep]
    f = sm.GLM(y[mask], x, family=sm.families.Binomial()).fit(cov_type='HC1', maxiter=100)
    if not f.converged or not np.isfinite(f.bse).all():
        raise RuntimeError('Logistic model failed convergence/finite covariance check')
    return f


def interaction(data):
    d, s, g, e, sens, y, p, b, raw, carrier, brca, cov, technical = data
    n = len(d)
    forest, ratios, comparisons, status = [], [], [], []

    def strata(label, covariates, mask=None, c=None, groups=False):
        mask = np.ones(n, bool) if mask is None else mask
        c = carrier if c is None else c
        X = np.column_stack([np.ones(n), p, covariates])
        fits = []
        for name, m in [('Non-carriers', ~c), ('All carriers', c)]:
            mm = mask & m
            if min(y[mm].sum(), (1 - y[mm]).sum()) <= 20:
                status.append(
                    dict(analysis=label + ' ' + name, status='suppressed_insufficient_stratum')
                )
                return
            f = fit_logistic(y, X, mm)
            fits.append(f)
            if groups:
                forest.append(
                    dict(
                        design=label, stratum=name, measure='OR', **estimate(f.params[1], f.bse[1])
                    )
                )
        ratios.append(
            dict(
                analysis=label,
                measure='ratio_OR',
                **ratio_estimate(
                    fits[1].params[1], fits[1].bse[1], fits[0].params[1], fits[0].bse[1]
                ),
            )
        )
        if groups:
            for name, m in [('BRCA1/2 carriers', brca), ('Other panel carriers', carrier & ~brca)]:
                mm = mask & m
                if min(y[mm].sum(), (1 - y[mm]).sum()) > 20:
                    f = fit_logistic(y, X, mm)
                    forest.append(
                        dict(
                            design=label,
                            stratum=name,
                            measure='OR',
                            **estimate(f.params[1], f.bse[1]),
                        )
                    )
        print('Completed interaction: ' + label, flush=True)
        write_table('stratum_effects', forest)
        write_table('sensitivity_ratios', ratios)

    strata('Paper covariates', cov, groups=True)
    spline = np.asarray(
        dmatrix(
            'bs(age, df=4, degree=3, include_intercept=False)-1',
            {'age': d.age_at_last_obs.to_numpy()},
        )
    )
    strata('Age at last EHR as spline (df=4 assumed)', np.column_stack([spline, cov[:, 1:]]))
    strata('No covariates', np.empty((n, 0)))
    age_case = np.where(y == 1, d.age_at_diagnosis, d.age_at_last_obs)
    valid = np.isfinite(age_case)
    strata(
        'Diagnosis age used for cases',
        np.column_stack([np.nan_to_num(age_case), cov[:, 1:]]),
        valid,
    )
    strata('Multigene carriers excluded', cov, s.carried_gene_count_C.to_numpy() <= 1)
    strata('PTV carriers only', cov, c=sens.B_raw_C_PTV_only.to_numpy() > 0)
    strata('Legacy ten-gene panel', cov, c=s.B_raw_A.to_numpy() > 0)
    strata(
        'Rule of two (single-date cases dropped)',
        cov,
        (y == 0) | e.is_case_rule_of_two.to_numpy(bool),
    )
    strata('First-record cases dropped', cov, (y == 0) | ~e.bc_is_first_record.to_numpy(bool))
    enroll = e.age_at_enrollment.to_numpy(float)
    finite = np.isfinite(enroll) & np.isfinite(d.age_at_diagnosis.to_numpy(float))
    early = finite & (d.age_at_diagnosis.to_numpy() < enroll + 1)
    late = finite & (d.age_at_diagnosis.to_numpy() >= enroll + 1)
    strata('Cases before enrolment age + 1 year', cov, (y == 0) | early)
    strata('Cases at/after enrolment age + 1 year', cov, (y == 0) | late)
    exact = INTERNAL / 'dates.csv.gz'
    if exact.exists():
        dates = pd.read_csv(exact).set_index('person_id').loc[d.index]
        dates = dates.apply(lambda x: pd.to_datetime(x) if x.name.endswith('_date') else x)
        strata(
            'Self-reported BC controls excluded',
            cov,
            (y == 1) | ~dates.self_reported_bc.to_numpy(bool),
        )
        mas = dates.first_mastectomy_date
        exclusion = mas.notna() & ((y == 0) | (mas <= dates.first_bc_date - pd.Timedelta(days=183)))
        strata('Mastectomy-related records excluded', cov, ~exclusion.to_numpy())
        valid = dates.enrollment_date.notna() & dates.first_bc_date.notna()
        if valid.any():
            strata(
                'Cases before exact enrolment date',
                cov,
                (y == 0) | (valid & (dates.first_bc_date < dates.enrollment_date)).to_numpy(),
            )
            strata(
                'Cases >=1 year after exact enrolment date',
                cov,
                (y == 0)
                | (
                    valid & (dates.first_bc_date >= dates.enrollment_date + pd.Timedelta(days=365))
                ).to_numpy(),
            )
        else:
            status.append(
                dict(analysis='Exact enrolment sensitivities', status='unavailable_enrollment_date')
            )
    else:
        status.append(dict(analysis='Exact-date sensitivities', status='dates_task_required'))
    status.append(
        dict(analysis='Within-European PCs', status='unavailable_source_PC_specification')
    )
    status.append(dict(analysis='Spline age', status='assumption_df4_not_verified_against_source'))
    status.append(
        dict(
            analysis='PTV-only and multigene carrier grouping',
            status='PDF_interpretation_other_carriers_remain_controls_or_excluded_as_labelled',
        )
    )
    # All interaction forms use the same observations. Include both carrier-class
    # indicators and continuous weight in an encompassing model for valid LRTs.
    other = (carrier & ~brca).astype(float)
    main = np.column_stack([np.ones(n), p, b, cov])
    forms = {
        'PRS_x_burden': np.column_stack([np.ones(n), p, b, p * b, cov]),
        'PRS_x_any_carrier': np.column_stack([np.ones(n), p, b, carrier, p * carrier, cov]),
        'PRS_x_gene_class': np.column_stack([np.ones(n), p, b, carrier, p * brca, p * other, cov]),
        'PRS_x_carried_weight': np.column_stack([np.ones(n), p, b, carrier, p * raw, cov]),
    }
    fitted = {}
    for name, X in forms.items():
        f = fit_logistic(y, X, np.ones(n, bool))
        fitted[name] = f
        comparisons.append(
            dict(
                model=name,
                AIC=float(f.aic),
                log_likelihood=float(f.llf),
                parameters=len(f.params),
                n=len(y),
            )
        )
    # raw and standardized burden interactions span identical models when PRS
    # main effect is retained. This is an identity, not a new model preference.
    encompassing = np.column_stack([forms['PRS_x_gene_class'], p * raw])
    full = fit_logistic(y, encompassing, np.ones(n, bool))
    lr = []
    # Class and weight forms nested in a common model with both main classes.
    class_weight_reduced = forms['PRS_x_carried_weight']
    wf = fit_logistic(y, class_weight_reduced, np.ones(n, bool))
    for name, rf in [
        ('class_plus_weight_vs_class', fitted['PRS_x_gene_class']),
        ('class_plus_weight_vs_weight_with_class_main', wf),
    ]:
        df = len(full.params) - len(rf.params)
        stat = 2 * (full.llf - rf.llf)
        lr.append(dict(comparison=name, LR=float(stat), df=df, p_value=float(chi2.sf(stat, df))))
    write_table('nested_LRT', lr)
    save_json(
        'interaction_checks',
        dict(
            note='Alternative forms retain weighted burden and add any-carrier main effect; carried weight = raw total weight. Main-effect choice is an explicit independent specification until source verification.',
            primary_ratio=ratios[0],
            status=status,
        ),
    )
    write_table('interaction_forms', comparisons)
    write_table('interaction_status', status)


def query(sql, ids=None, cap_gib=100):
    from google.cloud import bigquery

    client = bigquery.Client(project=os.environ['GOOGLE_PROJECT'])
    params = (
        []
        if ids is None
        else [bigquery.ArrayQueryParameter('person_ids', 'INT64', [int(i) for i in ids])]
    )
    cfg = bigquery.QueryJobConfig(dry_run=True, use_query_cache=False, query_parameters=params)
    dry = client.query(sql, job_config=cfg)
    if dry.total_bytes_processed > cap_gib * 1024**3:
        raise RuntimeError('Read-only query exceeds configured scan cap')
    print(f'BigQuery read-only scan: {dry.total_bytes_processed/1024**3:.2f} GiB', flush=True)
    cfg = bigquery.QueryJobConfig(
        query_parameters=params, maximum_bytes_billed=int(cap_gib * 1024**3)
    )
    return client.query(sql, job_config=cfg).result(timeout=600).to_dataframe()


def dates_task(data):
    from scripts.functions.ehr_age_audit import ehr_query

    d = data[0]
    cdr = os.environ['WORKSPACE_CDR']
    sql = ehr_query(cdr)
    dates = query(sql, d.index, cap_gib=2000)
    # exact enrolment comes from person schema; adjusted after schema probe.
    probe = json.loads((REPORT / '62_review_1008_schema_probe.json').read_text())
    person_fields = {x['name'] for x in probe['person']}
    if 'research_id' in person_fields:
        pass
    enrollment = 'enrollment_date' if 'enrollment_date' in person_fields else None
    if enrollment is None:
        dates['enrollment_date'] = pd.NaT
    else:
        extra = query(
            f'''WITH cohort AS (SELECT person_id FROM UNNEST(@person_ids) person_id)
        SELECT person_id, DATE({enrollment}) AS enrollment_date FROM `{cdr}.person` JOIN cohort USING(person_id)''',
            d.index,
        )
        dates = dates.merge(extra, on='person_id', validate='one_to_one')
    sql = f'''WITH cohort AS (SELECT person_id FROM UNNEST(@person_ids) person_id),
    sr AS (SELECT person_id, LOGICAL_OR(question_concept_id=836772 AND answer_concept_id=43528500) AS self_reported_bc
        FROM `{cdr}.ds_survey` JOIN cohort USING(person_id) GROUP BY person_id),
    mas AS (SELECT p.person_id, MIN(p.procedure_date) AS first_mastectomy_date
        FROM `{cdr}.procedure_occurrence` p JOIN cohort USING(person_id)
        JOIN `{cdr}.concept` c ON p.procedure_concept_id=c.concept_id
        WHERE REGEXP_CONTAINS(LOWER(c.concept_name),r'mastectomy')
        AND NOT REGEXP_CONTAINS(LOWER(c.concept_name),r'prosthesis|bra|implant|reconstruction|radiation|radiotherapy|catheter|specimen|pathology|spectroscopy')
        GROUP BY person_id)
    SELECT cohort.person_id, COALESCE(sr.self_reported_bc,FALSE) AS self_reported_bc,mas.first_mastectomy_date
    FROM cohort LEFT JOIN sr USING(person_id) LEFT JOIN mas USING(person_id)'''
    extra = query(sql, d.index, cap_gib=100)
    dates = dates.merge(extra, on='person_id', validate='one_to_one')
    INTERNAL.mkdir(exist_ok=True, parents=True)
    dates.to_csv(INTERNAL / 'dates.csv.gz', index=False)
    save_json(
        'dates_methods',
        dict(
            enrollment='not_available' if enrollment is None else enrollment,
            mastectomy='name-based mapped procedure concepts; reviewer exact code set unavailable',
            self_reported='harmonized question 836772 answer Self; source version-specific definition unavailable',
            observation_domains='condition, procedure, drug, measurement, visit',
            status='Exact enrolment sensitivities unavailable; EHR-only age-scale designs available',
        ),
    )


def riskset_moments(entry, exit, x, eta, times):
    """Breslow risk sets entry < t <= exit, via sorted cumulative moments."""
    x = np.asarray(x, float)
    times = np.asarray(times, float)
    eta = np.asarray(eta, float)
    shift = float(np.max(eta))
    w = np.exp(eta - shift)
    ai = np.argsort(entry)
    bi = np.argsort(exit)
    # entry is strict; exits at the event time stay in the risk set.
    a = np.searchsorted(np.asarray(entry)[ai], times, side='left')
    b = np.searchsorted(np.asarray(exit)[bi], times, side='left')

    def moment(v):
        ca = np.r_[0.0, np.cumsum(v[ai])]
        cb = np.r_[0.0, np.cumsum(v[bi])]
        return ca[a] - cb[b]

    den = moment(w)
    if (den <= 0).any():
        raise RuntimeError('Empty hazard-weighted risk set')
    m = np.column_stack([moment(w * x[:, j]) / den for j in range(x.shape[1])])
    V = np.empty((len(times), x.shape[1], x.shape[1]))
    for j in range(x.shape[1]):
        for k in range(j + 1):
            v = moment(w * x[:, j] * x[:, k]) / den - m[:, j] * m[:, k]
            V[:, j, k] = V[:, k, j] = v
    return m, V, den, shift


def prepare_survival(data, design='first_ehr_plus1'):
    d, s, g, e, sens, y, p, b, raw, carrier, brca, cov, technical = data
    dates = pd.read_csv(INTERNAL / 'dates.csv.gz').set_index('person_id').loc[d.index]
    for c in dates.columns:
        if c.endswith('_date'):
            dates[c] = pd.to_datetime(dates[c])

    def age(c):
        return ((dates[c] - dates.birth_date).dt.total_seconds() / 86400 / 365.25).to_numpy(float)

    first, last, bc, death = (
        age('ehr_first_date'),
        age('ehr_last_date'),
        age('first_bc_date'),
        age('death_date'),
    )
    entry = np.maximum(18, first + 1) if design == 'first_ehr_plus1' else np.full(len(d), 18.0)
    if design == 'first_ehr':
        entry = np.maximum(18, first)
    end = np.fmin(last, death)
    event = (y == 1) & np.isfinite(bc) & (bc <= end)
    exit = np.where(event, bc, end)
    valid = (
        np.isfinite(entry)
        & np.isfinite(exit)
        & (exit > entry)
        & ((y == 0) | (np.isfinite(bc) & (bc > entry)))
    )
    return entry, exit, event, valid, dates


def fit_cox(entry, exit, event, x, mask, strata=None):
    from types import SimpleNamespace

    xx = x[mask]
    keep = np.ptp(xx, axis=0) > 1e-12
    scale = xx[:, keep].std(axis=0)
    model = PHReg(
        exit[mask],
        xx[:, keep] / scale,
        status=event[mask].astype(int),
        entry=np.nextafter(entry[mask], np.inf),
        strata=None if strata is None else strata[mask],
        ties='breslow',
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        try:
            f = model.fit(disp=0, maxiter=100)
        except (np.linalg.LinAlgError, ValueError):
            f = None
        if (
            f is None
            or not np.isfinite(f.params).all()
            or not np.isfinite(f.bse).all()
            or np.max(np.abs(model.score(f.params))) > 1e-3
        ):
            f = model.fit(method='bfgs', disp=0, maxiter=1000, gtol=1e-9)
    if not np.isfinite(f.params).all() or not np.isfinite(f.bse).all():
        raise RuntimeError('Nonfinite Cox fit')
    scoremax = float(np.max(np.abs(model.score(f.params))))
    if scoremax > max(1e-3, 1e-7 * event[mask].sum()):
        raise RuntimeError(f'Cox score convergence check failed (max score {scoremax:.6g})')
    return (
        SimpleNamespace(
            params=f.params / scale,
            bse=f.bse / scale,
            covariance=f.cov_params() / np.outer(scale, scale),
        ),
        keep,
    )


def decompose(case_x, ages, means, V):
    centered = case_x - case_x.mean(axis=0)
    tc = ages - ages.mean()
    inv = np.linalg.inv(centered.T @ centered)
    expected = inv @ (means.T @ tc)
    remainder = inv @ ((case_x - means).T @ tc)
    observed = np.linalg.lstsq(centered, tc, rcond=None)[0]
    variance = inv @ np.einsum('i,ijk->jk', tc**2, V) @ inv
    error = float(np.max(np.abs(observed - expected - remainder)))
    assert error < 1e-8
    return observed, expected, remainder, np.sqrt(np.maximum(0, np.diag(variance))), error


def frailty_hr(beta, theta, H):
    return np.exp(beta) * (1 + theta * np.asarray(H)) / (1 + theta * np.asarray(H) * np.exp(beta))


def split_intervals(entry, exit, event, edges):
    rows = []
    for k, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        a = np.maximum(entry, lo)
        b = np.minimum(exit, hi)
        ok = b > a
        ix = np.flatnonzero(ok)
        rows.extend(
            zip(ix, a[ok], b[ok], (event[ok] & (exit[ok] <= hi)).astype(int), np.full(len(ix), k))
        )
    return np.asarray(rows, float)


def simulate_case_coefficients(
    entry, exit, x, eta, times, H, case_mean, case_sd, simulations, theta=0.0, seed=20261014
):
    """Conditional survival after entry; fixed score scaling from observed cases."""
    rng = np.random.default_rng(seed)
    eta = eta - np.max(eta)
    rr = np.exp(eta)
    Hentry = np.interp(entry, times, H, left=0, right=H[-1])
    rows = []
    for rep in range(simulations):
        u = np.maximum(rng.random(len(entry)), np.finfo(float).tiny)
        if theta == 0:
            target = Hentry - np.log(u) / rr
        else:
            # Gamma mixing conditional on surviving to entry. Inverse of
            # S(t|entry)=[(1+theta rr H(t))/(1+theta rr H(entry))]^(-1/theta).
            target = ((1 + theta * rr * Hentry) * u ** (-theta) - 1) / (theta * rr)
        j = np.searchsorted(H, target, side='left')
        st = np.full(len(entry), np.inf)
        ok = j < len(times)
        st[ok] = times[j[ok]]
        event = (st > entry) & (st <= exit)
        if event.sum() <= x.shape[1] + 2:
            continue
        xx = x[event].copy()
        xx[:, :2] = (xx[:, :2] - case_mean) / case_sd
        X = np.column_stack([np.ones(event.sum()), xx])
        fit = np.linalg.lstsq(X, st[event], rcond=None)[0]
        rows.append([rep, fit[1], fit[2]])
    if not rows:
        raise RuntimeError('No simulated cohorts yielded enough events')
    return np.asarray(rows)


def cox_task(data, simulations=100):
    d, s, g, e, sens, y, p, b, raw, carrier, brca, cov, technical = data
    if not (INTERNAL / 'dates.csv.gz').exists():
        dates_task(data)
    outputs = []
    decomp_rows = []
    age_rows = []
    simulation_rows = []
    checks = []
    age_ids = pd.read_csv(INPUTS['age_ids']).person_id
    case = d.index.isin(age_ids)
    age = d.age_at_diagnosis.to_numpy(float)[case]
    # Table 3 standardises PRS and burden within the diagnosis-age cohort.
    cmean = np.array([d.PRS311.to_numpy()[case].mean(), raw[case].mean()])
    csd = np.array([d.PRS311.to_numpy()[case].std(ddof=1), raw[case].std(ddof=1)])
    age_x = np.column_stack(
        [
            (d.PRS311.to_numpy() - cmean[0]) / csd[0],
            (raw - cmean[1]) / csd[1],
            d[PC].to_numpy(float),
        ]
    )
    ols = sm.OLS(age, np.column_stack([np.ones(case.sum()), age_x[case]])).fit(cov_type='HC1')
    for j, term in enumerate(['PRS311', 'Rare burden']):
        decomp_rows.append(
            dict(
                design='Table3_original',
                term=term,
                expected=np.nan,
                remainder=np.nan,
                remainder_se=np.nan,
                p_against_expected=np.nan,
                **estimate(ols.params[j + 1], ols.bse[j + 1], False),
            )
        )
    x = np.column_stack([p, b, d[PC].to_numpy(float), technical])
    for design in ['first_ehr_plus1', 'adult18']:
        entry, exit, event, valid, dates = prepare_survival(data, design)
        f, keep = fit_cox(entry, exit, event, x, valid)
        eta = x[:, keep] @ f.params
        checks.append(
            dict(
                design=design,
                n=disclosed_size(int(valid.sum()), len(d)),
                events=disclosed_size(int(event[valid].sum()), int(y.sum())),
                ties='breslow',
                score_converged=True,
            )
        )
        print('Completed Cox design: ' + design, flush=True)
        # Fit PRS separately by stratum, with PC/technical nuisance adjustments.
        if design == 'first_ehr_plus1':
            sx = np.column_stack([p, d[PC].to_numpy(float), technical])
            for label, mask in [
                ('Non-carriers', ~carrier),
                ('All carriers', carrier),
                ('BRCA1/2 carriers', brca),
                ('Other panel carriers', carrier & ~brca),
            ]:
                mm = valid & mask
                if event[mm].sum() <= 20:
                    continue
                sf, sk = fit_cox(entry, exit, event, sx, mm)
                outputs.append(
                    dict(
                        design='Age-scale EHR + 1 year',
                        stratum=label,
                        measure='HR',
                        **estimate(sf.params[0], sf.bse[0]),
                    )
                )
            birth_decade = (dates.birth_date.dt.year.to_numpy() // 10) * 10
            for label, mask in [('Non-carriers', ~carrier), ('All carriers', carrier)]:
                sf, sk = fit_cox(entry, exit, event, sx, valid & mask, strata=birth_decade)
                outputs.append(
                    dict(
                        design='Exact dates / birth-decade strata',
                        stratum=label,
                        measure='HR',
                        **estimate(sf.params[0], sf.bse[0]),
                    )
                )
            write_table('cox_stratum_effects', outputs)
        # Means for every Table 3 case age, including cases excluded from
        # incident follow-up. This preserves the case design in the PDF formula.
        means, V, den, shift = riskset_moments(
            entry[valid], exit[valid], age_x[valid], eta[valid], age
        )
        ob, ex, rem, se, error = decompose(age_x[case], age, means, V)
        for j, term in enumerate(['PRS311', 'Rare burden']):
            decomp_rows.append(
                dict(
                    design=design,
                    term=term,
                    expected=ex[j],
                    remainder=rem[j],
                    remainder_se=se[j],
                    p_against_expected=2 * norm.sf(abs(rem[j] / se[j])),
                    **estimate(ob[j], ols.bse[j + 1], False),
                )
            )
        checks[-1]['decomposition_identity_max_error'] = error
        write_table('age_decomposition', decomp_rows)
        # Reconstruct Breslow jumps explicitly: PHReg baseline returns
        # pre-jump cumulative hazards, unsuitable for direct inverse sampling.
        et, counts = np.unique(exit[valid & event], return_counts=True)
        _, _, risk_den, simshift = riskset_moments(
            entry[valid], exit[valid], x[valid, :2], eta[valid], et
        )
        H = np.cumsum(counts / risk_den)
        simx = np.column_stack([d.PRS311.to_numpy(), raw, d[PC].to_numpy(float)])[valid]
        potential_exit = np.fmin(
            ((dates.ehr_last_date - dates.birth_date).dt.days / 365.25).to_numpy(float),
            ((dates.death_date - dates.birth_date).dt.days / 365.25).to_numpy(float),
        )
        sim = simulate_case_coefficients(
            entry[valid], potential_exit[valid], simx, eta[valid], et, H, cmean, csd, simulations
        )
        for j, term in enumerate(['PRS311', 'Rare burden']):
            vals = sim[:, j + 1]
            simulation_rows.append(
                dict(
                    design=design,
                    theta=0,
                    term=term,
                    simulations=len(vals),
                    mean=float(vals.mean()),
                    SD=float(vals.std(ddof=1)),
                    MCSE=float(vals.std(ddof=1) / np.sqrt(len(vals))),
                    CI_low=float(np.quantile(vals, 0.025)),
                    CI_high=float(np.quantile(vals, 0.975)),
                    expected_riskset=float(ex[j]),
                    status='conditional_simulation_assumption_source_unavailable',
                )
            )
        if design == 'first_ehr_plus1':
            # Keep full data only in memory; serialise aggregate coefficients.
            main = (entry, exit, event, valid, dates, eta, et, H, potential_exit)
    write_table('simulation_summary', simulation_rows)
    entry, exit, event, valid, dates, eta, et, H, potential_exit = main
    bandedges = [18, 45, 55, 65, np.inf]
    intervals = split_intervals(entry[valid], exit[valid], event[valid], bandedges)
    ix = intervals[:, 0].astype(int)
    band = intervals[:, 4].astype(int)
    indicator = np.eye(4)[band]
    common = np.column_stack([d[PC].to_numpy(float), technical])[valid][ix]
    for label, variables, names in [
        ('any_carrier', np.column_stack([p, carrier]), ['PRS311', 'Any carrier']),
        (
            'gene_classes',
            np.column_stack([p, brca, carrier & ~brca]),
            ['PRS311_class_adjusted', 'BRCA1/2', 'Other carriers'],
        ),
    ]:
        vx = variables[valid][ix]
        time_x = np.column_stack(
            [indicator * vx[:, j, None] for j in range(vx.shape[1])] + [common]
        )
        bf, bk = fit_cox(
            intervals[:, 1],
            intervals[:, 2],
            intervals[:, 3].astype(bool),
            time_x,
            np.ones(len(ix), bool),
        )
        original = np.flatnonzero(bk)
        for j, term in enumerate(names):
            for k in range(4):
                coefficient = int(np.flatnonzero(original == j * 4 + k)[0])
                age_rows.append(
                    dict(
                        band=['<45', '45–54', '55–64', '≥65'][k],
                        term=term,
                        **estimate(bf.params[coefficient], bf.bse[coefficient]),
                    )
                )
        write_table('age_specific_HR', age_rows)
        print('Completed shared-nuisance interval-split Cox: ' + label, flush=True)
    # Frailty theta is a scenario until the reviewer's estimation code is
    # available. Do not claim reproduction of theta=.28 or its confidence bound.
    for theta in [0.28, 0.73]:
        sim = simulate_case_coefficients(
            entry[valid],
            potential_exit[valid],
            np.column_stack([d.PRS311.to_numpy(), raw, d[PC].to_numpy(float)])[valid],
            eta[valid],
            et,
            H,
            cmean,
            csd,
            simulations,
            theta=theta,
            seed=20261014,
        )
        for j, term in enumerate(['PRS311', 'Rare burden']):
            vals = sim[:, j + 1]
            simulation_rows.append(
                dict(
                    design='first_ehr_plus1',
                    theta=theta,
                    term=term,
                    simulations=len(vals),
                    mean=float(vals.mean()),
                    SD=float(vals.std(ddof=1)),
                    MCSE=float(vals.std(ddof=1) / np.sqrt(len(vals))),
                    CI_low=float(np.quantile(vals, 0.025)),
                    CI_high=float(np.quantile(vals, 0.975)),
                    expected_riskset=np.nan,
                    status='scenario_PDF_theta_not_reestimated_baseline_not_refitted',
                )
            )
    write_table('simulation_summary', simulation_rows)
    curve = []
    # A constant individual effect inferred from a PH carrier model is only a
    # scenario, not an identified individual frailty effect.
    cf, ck = fit_cox(
        entry,
        exit,
        event,
        np.column_stack([brca.astype(float), p, d[PC].to_numpy(float), technical]),
        valid,
    )
    # Cox baseline is normalised by exp(max eta); undo this normalization.
    actual_H = H / np.exp(np.max(eta[valid]))
    for t in np.linspace(25, 85, 121):
        ht = float(np.interp(t, et, actual_H, left=0, right=actual_H[-1]))
        curve.append(
            dict(
                age=t,
                theta=0.73,
                HR=float(frailty_hr(cf.params[0], 0.73, ht)),
                status='illustrative_PH_baseline_not_source_calibrated',
            )
        )
    write_table('frailty_curve', curve)
    save_json(
        'cox_checks',
        dict(
            designs=checks,
            simulations=simulations,
            seed=20261014,
            limitations=[
                'Reviewer second entry design assumed age18',
                'Model-based Cox errors and Breslow ties assumed',
                'Age-band model uses shared PC and technical effects; PRS and carrier class effects vary by age band',
                'Exact-date birth-decade strata use EHR+1year; source EHR-only start specification unavailable',
                'Gamma variance is a PDF scenario, not reestimated; baseline not recalibrated for frailty',
                'Simulation censoring uses clinical end/death rather than observed case diagnosis; source simulation mechanism unverified',
            ],
        ),
    )


def mh_estimate(y, positive, age, strata_extra=None):
    bands = np.digitize(age, [50, 60, 70])
    strata = (
        bands.astype(str)
        if strata_extra is None
        else np.char.add(np.asarray(strata_extra).astype(str), bands.astype(str))
    )
    tables = []
    for level in np.unique(strata):
        mask = strata == level
        tables.append(
            [
                [np.sum(mask & (y == 1) & positive), np.sum(mask & (y == 1) & ~positive)],
                [np.sum(mask & (y == 0) & positive), np.sum(mask & (y == 0) & ~positive)],
            ]
        )
    st = StratifiedTable(np.stack(tables, axis=2))
    lo, hi = st.oddsratio_pooled_confint()
    return dict(
        beta=float(st.logodds_pooled),
        se=float(st.logodds_pooled_se),
        estimate=float(st.oddsratio_pooled),
        CI_low=float(lo),
        CI_high=float(hi),
        p_value=float(st.test_null_odds().pvalue),
    )


def classify_family_flags(relative, self_named, list_answered, bc_ticked, item_answered):
    return np.select(
        [relative, self_named, bc_ticked & ~relative & ~self_named, list_answered],
        [
            'Relative named',
            'Only Self named',
            'BC ticked / relative skipped',
            'List answered / BC not ticked',
        ],
        default='List unanswered',
    )


def family_task(data):
    d, s, g, e, sens, y, p, b, raw, carrier, brca, cov, technical = data
    clinical = pd.read_csv(INPUTS['clinical']).set_index('person_id').loc[d.index]
    original = clinical.family_history_breast_cancer
    # Saved original coding; never use zeros for its missing respondents.
    print(
        'Original family-history labels: ' + str(sorted(original.dropna().unique().tolist())),
        flush=True,
    )
    pos = original.astype(str).str.lower().isin(['positive', 'positive_nonself', '1', '1.0'])
    neg = original.astype(str).str.lower().isin(['negative', 'negative_or_self_only', '0', '0.0'])
    if not pos.any() or not neg.any():
        raise RuntimeError('Unknown or empty original family-history coding')
    answered = (pos | neg).to_numpy()
    rows = []
    table = sm.stats.Table2x2(
        np.array(
            [
                [np.sum((y == 1) & pos), np.sum((y == 1) & neg)],
                [np.sum((y == 0) & pos), np.sum((y == 0) & neg)],
            ]
        )
    )
    lo, hi = table.oddsratio_confint()
    rows.append(
        dict(
            analysis='Table 1 coding / all respondents',
            kind='Original unadjusted OR',
            estimate=float(table.oddsratio),
            CI_low=lo,
            CI_high=hi,
            beta=float(table.log_oddsratio),
            se=float(table.log_oddsratio_se),
            p_value=float(table.oddsratio_pvalue()),
        )
    )
    cdr = os.environ['WORKSPACE_CDR']
    # Standardised ds_survey version columns are NULL in this release. Retain
    # first submission by harmonised survey, and label version inference absent.
    sql = f'''WITH cohort AS (SELECT person_id FROM UNNEST(@person_ids) person_id),
    first_submission AS (
      SELECT s.person_id,MIN(DATE(s.survey_datetime)) AS survey_date
      FROM `{cdr}.ds_survey` s JOIN cohort USING(person_id)
      WHERE s.survey='Personal and Family Health History' GROUP BY person_id),
    flags AS (
      SELECT f.person_id,f.survey_date,
      LOGICAL_OR(s.question_concept_id=836772 AND s.answer_concept_id IN (43528501,43528502,43528503,43528504,43528505,43528506)) AS relative_named,
      LOGICAL_OR(s.question_concept_id=836772 AND s.answer_concept_id=43528500) AS self_named,
      LOGICAL_OR(s.question_concept_id=43528515 AND s.answer_concept_id NOT IN (903079,903096)) AS list_answered,
      LOGICAL_OR(s.question_concept_id=43528515 AND s.answer_concept_id=43528499) AS bc_ticked,
      LOGICAL_OR(s.question_concept_id=836772) AS item_answered
      FROM first_submission f JOIN `{cdr}.ds_survey` s ON s.person_id=f.person_id AND DATE(s.survey_datetime)=f.survey_date
      AND s.survey='Personal and Family Health History'
      GROUP BY f.person_id,f.survey_date)
    SELECT * FROM flags'''
    flags = query(sql, d.index, cap_gib=100).set_index('person_id')
    # No participant-level survey results are written to disk.
    ff = d[['is_case', 'age_at_last_obs']].join(e[['age_at_enrollment']]).join(flags, how='inner')
    for name in ['relative_named', 'self_named', 'list_answered', 'bc_ticked', 'item_answered']:
        ff[name] = ff[name].fillna(False).astype(bool)
    corrected_positive = ff.relative_named.to_numpy()
    # Breast cancer ticked without self/relative is unclassifiable; no negative
    # assignment from a skipped relative item.
    corrected_negative = (
        ff.list_answered & (~ff.bc_ticked | ff.self_named) & ~ff.relative_named
    ).to_numpy()
    mask = (corrected_positive | corrected_negative) & ff.age_at_enrollment.notna().to_numpy()
    est = mh_estimate(
        ff.is_case.to_numpy()[mask], corrected_positive[mask], ff.age_at_enrollment.to_numpy()[mask]
    )
    rows.append(
        dict(
            analysis='Corrected harmonised survey / first submission',
            kind='Age-adjusted OR / version unavailable',
            **est,
        )
    )
    om = ff.item_answered.to_numpy() & ff.age_at_enrollment.notna().to_numpy()
    rows.append(
        dict(
            analysis='Original item / harmonised first submission',
            kind='Age-adjusted OR / version unavailable',
            **mh_estimate(
                ff.is_case.to_numpy()[om],
                corrected_positive[om],
                ff.age_at_enrollment.to_numpy()[om],
            ),
        )
    )
    routes = classify_family_flags(
        ff.relative_named.to_numpy(),
        ff.self_named.to_numpy(),
        ff.list_answered.to_numpy(),
        ff.bc_ticked.to_numpy(),
        ff.item_answered.to_numpy(),
    )
    route_rows = []
    for group, yy in [('Cases', 1), ('Controls', 0)]:
        gm = ff.is_case.to_numpy() == yy
        counts = [
            int(np.sum(gm & (routes == r)))
            for r in [
                'Relative named',
                'Only Self named',
                'BC ticked / relative skipped',
                'List answered / BC not ticked',
                'List unanswered',
            ]
        ]
        # Withhold all route percentages if any cell <=20, to avoid recovering
        # a small complementary cell from the remaining stack.
        disclose = all(c == 0 or c > 20 for c in counts)
        for route, count in zip(
            [
                'Relative named',
                'Only Self named',
                'BC ticked / relative skipped',
                'List answered / BC not ticked',
                'List unanswered',
            ],
            counts,
        ):
            route_rows.append(
                dict(
                    group=group,
                    route=route,
                    percent=100 * count / gm.sum() if disclose else np.nan,
                    count='withheld_joint_disclosure',
                    status='available' if disclose else 'suppressed_entire_stack',
                )
            )
    write_table('family_history_estimates', rows)
    write_table('family_history_routes', route_rows)
    save_json(
        'family_history_methods',
        dict(
            status='partial_version_specific_reproduction_unavailable',
            age_adjustment='Mantel-Haenszel age at enrollment bands <50,50-59,60-69,>=70; source age reference unavailable',
            first_submission='First calendar date in harmonised Personal and Family Health History survey',
            limitations=[
                'Survey version concept IDs and names are NULL; all source questions harmonised',
                'Cannot identify three original survey versions or reproduce their separate/pooled-first-version estimates',
                'Figure 5 includes harmonised first submissions, not a verified 2021-version subset',
                'Self-only answers negative; BC ticked without relative or self unclassifiable',
            ],
            corrected_estimate=est,
            reference='Lancet 2001: RR 1.80, 99% CI 1.69–1.91 as reported in supplied PDF',
        ),
    )


def plots_task():
    import matplotlib

    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter

    plt.rcParams['svg.fonttype'] = 'none'

    def read(name):
        return pd.read_csv(TABLE / f'62_review_1008_{name}.csv')

    def save(fig, number):
        fig.savefig(FIGURE / f'62_review_1008_figure{number}.svg', bbox_inches='tight')
        fig.savefig(INTERNAL / f'figure{number}_preview.png', dpi=110, bbox_inches='tight')
        plt.close(fig)

    def forest(frame, labels, title, xlabel, number):
        fig, ax = plt.subplots(figsize=(10, max(3, 0.42 * len(frame) + 1.2)))
        yy = np.arange(len(frame))[::-1]
        ax.errorbar(
            frame.estimate,
            yy,
            xerr=[frame.estimate - frame.CI_low, frame.CI_high - frame.estimate],
            fmt='o',
            capsize=3,
        )
        ax.set_yticks(yy, labels)
        ax.set_xscale('log')
        ax.axvline(1, color='gray', linestyle='--', linewidth=1)
        ax.xaxis.set_major_formatter(FuncFormatter(lambda value, position: f'{value:g}'))
        ax.xaxis.set_minor_formatter(FuncFormatter(lambda value, position: f'{value:g}'))
        ax.set_xlabel(xlabel)
        ax.set_title(title)
        ax.grid(False)
        for z, row in zip(yy, frame.itertuples()):
            ax.text(
                1.02,
                z,
                f'{row.estimate:.2f} ({row.CI_low:.2f}–{row.CI_high:.2f})',
                transform=ax.get_yaxis_transform(),
                va='center',
            )
        save(fig, number)

    f1 = pd.concat([read('stratum_effects'), read('cox_stratum_effects')], ignore_index=True)
    forest(
        f1,
        [f'{r.design}: {r.stratum}' for r in f1.itertuples()],
        'PRS311 effects within carrier strata',
        'Effect per full-cohort SD (OR / HR)',
        1,
    )
    f2 = read('sensitivity_ratios')
    forest(
        f2,
        f2.analysis.tolist(),
        'Carrier / non-carrier PRS OR ratio\nSource-matched where specified; assumptions documented',
        'Ratio of odds ratios',
        2,
    )
    dec = read('age_decomposition')
    sim = read('frailty_refitted_simulations')
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.set_axisbelow(True)
    for j, term in [(1, 'PRS311'), (0, 'Rare burden')]:
        observed = dec[(dec.design == 'Table3_original') & (dec.term == term)].iloc[0]
        expected = dec[(dec.design != 'Table3_original') & (dec.term == term)].expected
        frail = sim[sim.term == term]
        if len(frail):
            envelope_low = min(expected.min(), frail['mean'].min())
            envelope_high = max(expected.max(), frail['mean'].max())
            ax.barh(
                j,
                envelope_high - envelope_low,
                left=envelope_low,
                height=0.22,
                color='#f6bd87',
                align='center',
                zorder=2,
                label='Gamma scenario envelope / includes no frailty' if j == 1 else None,
            )
        ax.barh(
            j,
            expected.max() - expected.min(),
            left=expected.min(),
            height=0.22,
            color='#86b6d8',
            align='center',
            zorder=3,
            label='Risk-set expectation / two entry designs' if j == 1 else None,
        )
        ax.errorbar(
            observed.estimate,
            j,
            xerr=[[observed.estimate - observed.CI_low], [observed.CI_high - observed.estimate]],
            fmt='D',
            color='black',
            markersize=12,
            capsize=6,
            zorder=4,
        )
        ax.text(
            observed.estimate,
            j + 0.20,
            f'Observed {observed.estimate:.2f} ({observed.CI_low:.2f}–{observed.CI_high:.2f})',
            ha='center',
            va='bottom',
        )
    ax.axvline(0, ls='--', color='gray', linewidth=1)
    ax.set_yticks([0, 1], ['Rare burden', 'PRS311'])
    ax.set_xlabel('Years of recorded diagnosis age per case-cohort SD')
    ax.grid(axis='x', alpha=0.3)
    ax.set_ylim(-0.4, 1.4)
    ax.set_title('Observed and expected diagnosis-age coefficients')
    ax.legend(fontsize=9, loc='lower right')
    save(fig, 3)
    hr = read('age_specific_HR')
    curve = read('frailty_refitted_curve')
    curve = curve[np.isclose(curve.theta, 0.73)]
    fig, ax = plt.subplots(figsize=(8, 5))
    for term in ['BRCA1/2', 'Any carrier', 'Other carriers', 'PRS311']:
        sub = hr[hr.term == term]
        xx = np.arange(len(sub))
        ax.errorbar(
            xx,
            sub.estimate,
            yerr=[sub.estimate - sub.CI_low, sub.CI_high - sub.estimate],
            fmt='o-',
            capsize=3,
            label=term,
        )
    # Band-average x locations for the illustrative depletion scenario.
    at = [35, 50, 60, 72]
    depletion = np.interp(at, curve.age, curve.HR)
    ax.plot(
        np.arange(4), depletion, '--', label='Gamma .73 scenario; independently refitted baseline'
    )
    ax.set_xticks(range(4), ['<45', '45–54', '55–64', '≥65'])
    ax.set_yscale('log')
    ax.set_xlabel('Attained age')
    ax.set_ylabel('Hazard ratio')
    ax.set_title('Age-specific effects; interval-split Cox')
    ax.grid(False)
    ax.legend(fontsize=8)
    save(fig, 4)
    routes = read('family_history_routes')
    fig, ax = plt.subplots(figsize=(8, 4))
    bottom = np.zeros(2)
    if routes.percent.notna().all():
        for route in routes.route.drop_duplicates():
            sub = routes[routes.route == route].set_index('group').loc[['Cases', 'Controls']]
            vals = sub.percent.to_numpy()
            ax.bar(['Cases', 'Controls'], vals, bottom=bottom, label=route)
            bottom += vals
        ax.legend(fontsize=8, bbox_to_anchor=(1.02, 1), loc='upper left')
        ax.set_ylim(0, 100)
    else:
        ax.text(
            0.5,
            0.5,
            'Route stack withheld for complementary disclosure',
            ha='center',
            transform=ax.transAxes,
        )
    ax.set_ylabel('Percent of first survey submissions')
    ax.set_title(
        'Harmonised family-history answer routes\nOriginal survey-version subset unavailable'
    )
    ax.grid(False)
    save(fig, 5)
    fh = read('family_history_estimates')
    reference = pd.DataFrame(
        [
            dict(
                analysis='Published first-degree-relative reference (RR, 99% CI)',
                estimate=1.80,
                CI_low=1.69,
                CI_high=1.91,
            )
        ]
    )
    fh = pd.concat([fh, reference], ignore_index=True)
    forest(
        fh,
        fh.analysis.tolist(),
        'Family-history estimates\nVersion-specific reproduction incomplete',
        'OR (95% CI); published reference RR (99% CI)',
        6,
    )
    write_table('figure1_plot_data', f1)
    write_table('figure6_plot_data', fh)
    print('Six editable SVG figures saved.', flush=True)


def frailty_task(data, simulations):
    from scripts.functions.frailty import profile_theta, fit_gamma, exposure_matrix

    d, s, g, e, sens, y, p, b, raw, carrier, brca, cov, technical = data
    entry, exit, event, valid, dates = prepare_survival(data)
    mask = valid & ~carrier
    nuisance = np.column_stack([d[PC].to_numpy(float), technical])
    x = np.column_stack([p, nuisance])
    scale = x[mask].std(axis=0)
    x = x / scale
    theta, low, high, best, profile = profile_theta(x[mask], entry[mask], exit[mask], event[mask])
    write_table(
        'frailty_theta_profile', [dict(theta=t, negative_log_likelihood=f) for t, f in profile]
    )
    save_json(
        'frailty_theta_estimation',
        dict(
            theta=theta,
            CI_low=low,
            CI_high=high,
            method='Gamma marginal survival likelihood conditional on delayed entry; profile LR 95% interval; piecewise constant baseline 18,35,45,55,65,75,100,130',
            status='independent_baseline_and_profile_specification_not_verified_against_reviewer',
            review_target_theta=0.28,
            review_target_upper=0.73,
        ),
    )
    print('Completed independent gamma-frailty variance profile.', flush=True)
    # Refit the all-cohort baseline for each theta rather than applying frailty
    # to an unchanged Cox baseline. Individual risk effects include gene class.
    gx = np.column_stack([p, brca.astype(float), (carrier & ~brca).astype(float), nuisance])
    gs = gx[valid].std(axis=0)
    xx = gx[valid] / gs
    potential = np.fmin(
        ((dates.ehr_last_date - dates.birth_date).dt.days / 365.25).to_numpy(float),
        ((dates.death_date - dates.birth_date).dt.days / 365.25).to_numpy(float),
    )[valid]
    ids = pd.read_csv(INPUTS['age_ids']).person_id
    case = d.index.isin(ids)
    cm = np.array([d.PRS311.to_numpy()[case].mean(), raw[case].mean()])
    cs = np.array([d.PRS311.to_numpy()[case].std(ddof=1), raw[case].std(ddof=1)])
    simx = np.column_stack([d.PRS311.to_numpy(), raw, d[PC].to_numpy(float)])[valid]
    sim_rows = []
    curves = []
    effects = []
    edges = (18, 35, 45, 55, 65, 75, 100, 130)
    # These theta values are explicitly labelled PDF scenarios. The independent
    # theta estimate is additionally simulated, avoiding substitution of targets.
    scenarios = [0.0, 0.28, 0.73, round(theta, 6)]
    age_pattern_path = REPORT / '62_review_1008_age_pattern_frailty_estimation.json'
    if age_pattern_path.exists():
        scenarios.append(round(json.loads(age_pattern_path.read_text())['theta'], 6))
    scenarios = sorted(set(scenarios))
    for th in scenarios:
        result = fit_gamma(xx, entry[valid], exit[valid], event[valid], th, edges=edges)
        eta = xx @ result.x[: xx.shape[1]]
        times = np.arange(18, 105, 0.05)
        H = exposure_matrix(times, edges) @ np.exp(result.x[xx.shape[1] :])
        # simulation helper internally centres eta, so H must be rescaled too.
        sim = simulate_case_coefficients(
            entry[valid],
            potential,
            simx,
            eta,
            times,
            H * np.exp(np.max(eta)),
            cm,
            cs,
            simulations,
            theta=th,
        )
        for j, term in enumerate(['PRS311', 'Rare burden']):
            vals = sim[:, j + 1]
            sim_rows.append(
                dict(
                    theta=th,
                    term=term,
                    simulations=len(vals),
                    mean=vals.mean(),
                    SD=vals.std(ddof=1),
                    MCSE=vals.std(ddof=1) / np.sqrt(len(vals)),
                    CI_low=np.quantile(vals, 0.025),
                    CI_high=np.quantile(vals, 0.975),
                    status='gamma_refitted_baseline_independent_implementation',
                )
            )
        effects.append(
            dict(
                theta=th,
                PRS_individual_beta=result.x[0] / gs[0],
                BRCA_individual_beta=result.x[1] / gs[1],
                Other_individual_beta=result.x[2] / gs[2],
                negative_log_likelihood=result.fun,
            )
        )
        beta = result.x[1] / gs[1]
        for age, h in zip(times, H):
            curves.append(
                dict(
                    theta=th,
                    age=age,
                    HR=float(frailty_hr(beta, th, h)),
                    status='baseline_refitted_independent_model',
                )
            )
        print('Completed refitted gamma scenario: ' + str(th), flush=True)
    write_table('frailty_refitted_simulations', sim_rows)
    write_table('frailty_refitted_curve', curves)
    write_table('frailty_individual_effects', effects)


def age_pattern_task(data):
    """Calibrate gamma depletion from non-carrier PRS age effects as in PDF.

    Band midpoints and gamma-to-marginal baseline conversion are explicit
    independent choices until the reviewer fitting implementation is available.
    """
    from scipy.optimize import minimize_scalar, brentq
    from scripts.functions.frailty import fit_gamma, exposure_matrix

    d, s, g, e, sens, y, p, b, raw, carrier, brca, cov, technical = data
    entry, exit, event, valid, dates = prepare_survival(data)
    mask = valid & ~carrier
    intervals = split_intervals(entry[mask], exit[mask], event[mask], [18, 45, 55, 65, np.inf])
    ix = intervals[:, 0].astype(int)
    band = intervals[:, 4].astype(int)
    nuisance = np.column_stack([d[PC].to_numpy(float), technical])
    xx = np.column_stack([np.eye(4)[band] * p[mask][ix, None], nuisance[mask][ix]])
    fit, keep = fit_cox(
        intervals[:, 1], intervals[:, 2], intervals[:, 3].astype(bool), xx, np.ones(len(ix), bool)
    )
    if not keep[:4].all():
        raise RuntimeError('Missing non-carrier age-band PRS information')
    beta = fit.params[:4]
    V = fit.covariance[:4, :4]
    inv = np.linalg.inv(V)
    write_table(
        'noncarrier_age_pattern',
        [
            dict(band=['<45', '45–54', '55–64', '≥65'][j], **estimate(beta[j], np.sqrt(V[j, j])))
            for j in range(4)
        ],
    )
    baseline_x = np.column_stack([p, nuisance])[mask]
    scale = baseline_x.std(axis=0)
    baseline_x = baseline_x / scale
    edges = (18, 35, 45, 55, 65, 75, 100, 130)
    constant = fit_gamma(baseline_x, entry[mask], exit[mask], event[mask], 0, edges=edges)
    # Nuisance-marginal reference cumulative hazard; midpoint approximation.
    reference_risk = np.mean(np.exp(baseline_x[:, 1:] @ constant.x[1 : baseline_x.shape[1]]))
    mids = np.array([35.0, 50.0, 60.0, 72.0])
    Hm = exposure_matrix(mids, edges) @ np.exp(constant.x[baseline_x.shape[1] :]) * reference_risk
    profiles = {}

    def objective(theta):
        Hi = Hm if theta == 0 else np.expm1(theta * Hm) / theta

        def loss(individual_beta):
            prediction = np.log(frailty_hr(individual_beta, theta, Hi))
            residual = beta - prediction
            return float(residual @ inv @ residual)

        result = minimize_scalar(loss, bounds=(0, 2), method='bounded')
        profiles[float(theta)] = (float(result.fun), float(result.x))
        return float(result.fun)

    result = minimize_scalar(objective, bounds=(0, 3), method='bounded', options={'xatol': 1e-5})
    theta = float(result.x)
    if objective(0) < objective(theta):
        theta = 0.0
    best = objective(theta)
    cut = chi2.ppf(0.95, 1)
    root = lambda t: objective(t) - best - cut
    low = 0.0 if theta == 0 or root(0) < 0 else float(brentq(root, 0, theta))
    high = None if root(3) < 0 else float(brentq(root, theta, 3))
    write_table(
        'age_pattern_theta_profile',
        [
            dict(theta=t, chi_square=f, individual_PRS_beta=v)
            for t, (f, v) in sorted(profiles.items())
        ],
    )
    trend_x = np.column_stack([np.ones(4), (mids - 55) / 10])
    trend_cov = np.linalg.inv(trend_x.T @ inv @ trend_x)
    trend = trend_cov @ trend_x.T @ inv @ beta
    save_json(
        'age_pattern_frailty_estimation',
        dict(
            theta=theta,
            CI_low=low,
            CI_high=high,
            PRS_HR_per_decade=estimate(trend[1], np.sqrt(trend_cov[1, 1])),
            method='GLS fit of exact gamma population HR formula to non-carrier age-band PRS log HRs; covariance retained; profile chi-square interval',
            baseline='Nuisance-marginal zero-frailty piecewise baseline converted to individual baseline at each theta',
            representative_ages=mids.tolist(),
            status='age_pattern_method_implemented_midpoints_and_baseline_independently_specified',
            review_target_theta=0.28,
            review_target_upper=0.73,
        ),
    )
    print('Completed non-carrier PRS age-pattern gamma calibration.', flush=True)


def validate_task():
    """Compare rounded published targets without claiming method equivalence."""

    def read(name):
        return pd.read_csv(TABLE / f'62_review_1008_{name}.csv')

    test_path = REPORT / '62_review_1008_synthetic_tests.json'
    if not test_path.exists():
        raise RuntimeError('Run tests/test_62_review_1008.py before validation')
    tests = json.loads(test_path.read_text())
    assert (
        tests['passed']
        and tests['analysis_file_sha256'] == hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    )
    # Also repair existing output from development runs before final release.
    checks_path = REPORT / '62_review_1008_cox_checks.json'
    checks = json.loads(checks_path.read_text())
    provenance = json.loads((REPORT / '62_review_1008_provenance.json').read_text())
    from importlib.metadata import version

    provenance['software_versions'] = {
        name: version(name)
        for name in [
            'numpy',
            'pandas',
            'scipy',
            'statsmodels',
            'patsy',
            'matplotlib',
            'google-cloud-bigquery',
        ]
    }
    provenance['internal_dates_sha256'] = hashlib.sha256(
        (INTERNAL / 'dates.csv.gz').read_bytes()
    ).hexdigest()
    provenance['method_helper_sha256'] = hashlib.sha256(
        (ROOT / 'scripts/functions/frailty.py').read_bytes()
    ).hexdigest()
    save_json('provenance', provenance)
    for design in checks['designs']:
        if isinstance(design['n'], int):
            design['n'] = disclosed_size(design['n'], provenance['cohort_n'])
        if isinstance(design['events'], int):
            design['events'] = disclosed_size(design['events'], provenance['case_n'])
    checks_path.write_text(json.dumps(checks, indent=2) + '\n')
    rows = []

    def add(analysis, metric, value, target, decimals=2, method='specified_in_PDF'):
        rows.append(
            dict(
                analysis=analysis,
                metric=metric,
                computed=float(value),
                PDF_target=float(target),
                difference=float(value - target),
                matches_report_rounding=round(float(value), decimals) == round(target, decimals),
                method_status=method,
            )
        )

    f = read('stratum_effects')
    for name, target in [
        ('Non-carriers', 1.56),
        ('All carriers', 1.25),
        ('BRCA1/2 carriers', 1.11),
        ('Other panel carriers', 1.50),
    ]:
        add('Logistic strata', name, f[f.stratum == name].iloc[0].estimate, target)
    ratio = read('sensitivity_ratios').iloc[0]
    add('Interaction', 'OR ratio', ratio.estimate, 0.80)
    forms = read('interaction_forms')
    for name, target in [
        ('PRS_x_burden', 55823),
        ('PRS_x_any_carrier', 55831),
        ('PRS_x_gene_class', 55826),
        ('PRS_x_carried_weight', 55825),
    ]:
        add(
            'Interaction forms',
            name,
            forms[forms.model == name].iloc[0].AIC,
            target,
            0,
            'alternative_main_effect_specification_inferred',
        )
    dec = read('age_decomposition')
    for term, target in [('PRS311', -0.43), ('Rare burden', -1.82)]:
        add(
            'Table 3 observed',
            term,
            dec[(dec.design == 'Table3_original') & (dec.term == term)].iloc[0].estimate,
            target,
        )
    main = dec[dec.design == 'first_ehr_plus1']
    add(
        'Expected from Cox risk sets',
        'PRS311',
        main[main.term == 'PRS311'].iloc[0].expected,
        -0.30,
        2,
        'Cox_covariates_ties_independently_specified',
    )
    add(
        'Expected from Cox risk sets',
        'Rare burden',
        main[main.term == 'Rare burden'].iloc[0].expected,
        -0.72,
        2,
        'Cox_covariates_ties_independently_specified',
    )
    fh = read('family_history_estimates')
    add('Original family history', 'OR', fh.iloc[0].estimate, 0.71)
    add(
        'Corrected family history',
        'OR',
        fh.iloc[1].estimate,
        1.76,
        2,
        'version-specific_source_data_unavailable',
    )
    frail = read('frailty_refitted_simulations')
    add(
        'Gamma theta .73 scenario',
        'Expected burden',
        frail[(np.isclose(frail.theta, 0.73)) & (frail.term == 'Rare burden')].iloc[0]['mean'],
        -1.45,
        2,
        'independent_baseline_and_simulation',
    )
    theta = json.loads((REPORT / '62_review_1008_frailty_theta_estimation.json').read_text())
    add(
        'Frailty theta',
        'Estimate',
        theta['theta'],
        0.28,
        2,
        'independent_likelihood_not_review_age_pattern_estimator',
    )
    add(
        'Frailty theta',
        'Upper 95% bound',
        theta['CI_high'],
        0.73,
        2,
        'independent_likelihood_not_review_age_pattern_estimator',
    )
    pattern = json.loads(
        (REPORT / '62_review_1008_age_pattern_frailty_estimation.json').read_text()
    )
    add(
        'Age-pattern frailty theta',
        'Estimate',
        pattern['theta'],
        0.28,
        2,
        'GLS_age_pattern_midpoint_and_baseline_specification_independent',
    )
    if pattern['CI_high'] is not None:
        add(
            'Age-pattern frailty theta',
            'Upper 95% bound',
            pattern['CI_high'],
            0.73,
            2,
            'GLS_age_pattern_midpoint_and_baseline_specification_independent',
        )
    age_hr = read('age_specific_HR')
    add(
        'Age-specific carrier HR',
        '<45',
        age_hr[(age_hr.term == 'Any carrier') & (age_hr.band == '<45')].iloc[0].estimate,
        8.4,
        1,
    )
    add(
        'Age-specific carrier HR',
        '>=65',
        age_hr[(age_hr.term == 'Any carrier') & (age_hr.band == '≥65')].iloc[0].estimate,
        2.7,
        1,
    )
    write_table('target_comparison', rows)
    outputs = list(FIGURE.glob('62_review_1008_figure*.svg'))
    assert len(outputs) == 6 and all(x.stat().st_size > 1000 for x in outputs)
    for x in TABLE.glob('62_review_1008_*.csv'):
        cols = pd.read_csv(x, nrows=0).columns
        assert not any(c in cols for c in ['person_id', 'sample_id', 'questionnaire_response_id'])
    save_json(
        'validation',
        dict(
            status='IMPLEMENTED_AND_EXECUTED_PARTIAL_SOURCE_REPRODUCTION',
            comparison_rows=len(rows),
            rounded_targets_matched=sum(r['matches_report_rounding'] for r in rows),
            six_svg_figures_verified=True,
            participant_identifier_columns_absent=True,
            source_access='SSH public-key authentication failed; HTTPS repository authentication unavailable',
            pending=[
                'Exact reviewer source commit and parameter verification',
                'Within-European PCs',
                'Exact enrollment dates',
                'Three original family-history survey versions and pooled-first-version selection',
                'Exact reviewer PRS age-pattern calibration choices',
            ],
            figure_sizes_bytes={str(x.relative_to(ROOT)): x.stat().st_size for x in outputs},
            analysis_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            synthetic_test_file='tests/test_62_review_1008.py',
            synthetic_tests=tests,
        ),
    )
    lines = [
        '# Reproduction of the 1008 review analyses',
        '',
        'Implemented and run on the frozen local study cohort. Exact reproduction is incomplete because the reviewer repository is inaccessible and several report specifications are missing. The decision cards and manuscript were not changed.',
        '',
        'Primary logistic stratum estimates, the carrier/non-carrier ratio, Table 3 observed diagnosis-age coefficients, four rounded AICs, age-specific any-carrier HRs and the original family-history OR match the report. Numerical matches alone do not verify undocumented model choices.',
        f"Synthetic verification: {tests['tests_run']} tests passed; test and implementation SHA-256 values are recorded in the validation JSON.",
        '',
        '## Computed estimates against report targets',
        '',
        '| Analysis | Metric | Computed | PDF target | Matches rounding | Method status |',
        '|---|---|---:|---:|---|---|',
    ]
    for r in rows:
        lines.append(
            f"| {r['analysis']} | {r['metric']} | {r['computed']:.6g} | {r['PDF_target']:.6g} | {r['matches_report_rounding']} | {r['method_status']} |"
        )
    lines += [
        '',
        '## What remains unresolved',
        '',
        '- The original source `review_c3_c4_c9/wf/` at reference commit `ca6bfba` could not be read. The reproducible scripts record every inferred model choice.',
        '- The exact enrollment date is absent from inspected person tables. Age-at-enrollment sensitivities were run; exact enrollment-date sensitivities were not substituted with approximate dates.',
        '- Within-European PC definitions are unavailable. The main study PCs remain fixed.',
        '- Survey-version IDs/names are NULL in ds_survey and the queried original observation question IDs map to the same harmonised questions. The current run cannot establish the three original versions. Corrected routing and Figures 5–6 use first harmonised survey submissions; they are not the original 2021-only and three-version results.',
        '- Cox ties use Breslow and model-based covariance; source tie and uncertainty choices are not given. The second entry design is assumed age 18. Exact-date birth-decade fits use the same first-EHR-plus-one-year entry pending source confirmation.',
        '- Gamma variance is estimated in two independent ways: conditional survival likelihood with piecewise baseline, and GLS calibration of non-carrier age-band PRS effects to the gamma depletion formula. The latter follows the report’s age-pattern approach, but representative ages and baseline calibration are independently specified; its estimate/CI differ from the reviewer. The theta=.28 and .73 values are labelled PDF scenarios; baselines and individual effects are refitted, and the independently calibrated theta is additionally simulated.',
        '- Risk-set decomposition is evaluated on the original Table 3 case design, including cases outside incident follow-up; its expected-plus-remainder identity passed. Monte Carlo case coefficients use fixed original case scaling and follow-up ending at clinical last record/death. They need not equal the decomposition under the realised case covariance.',
        '- Mastectomy and self-report sensitivities use documented independently reconstructed code definitions. One nested interaction-form comparison differs from the reviewer; no non-nested likelihood-ratio test is used.',
        '',
        '## Reproduce',
        '',
        'Run from the project root. The dates task uses read-only BigQuery queries and writes an ignored local dates cache. Family-history flags stay in memory. All published outputs are aggregate tables/plots.',
        '',
        '```bash',
        'python scripts/review/probe_1008.py',
        'python scripts/review/review_1008.py --task dates',
        'python scripts/review/review_1008.py --task interaction',
        'python scripts/review/review_1008.py --task cox --simulations 100',
        'python scripts/review/review_1008.py --task age-pattern',
        'python scripts/review/review_1008.py --task frailty --simulations 100',
        'python scripts/review/review_1008.py --task family',
        'python scripts/review/review_1008.py --task plots',
        'python tests/test_62_review_1008.py',
        'python scripts/review/review_1008.py --task validate',
        '```',
        '',
        '## Methods and software',
        '',
        'The supplied PDF pages 1–6 are the primary method source. The implementation uses full-cohort score standardisation for logistic/Cox risk models, diagnosis-case standardisation for Table 3, HC1 logistic/OLS errors, attained-age delayed entry, explicit Breslow risk-set moments, interval-split age-specific coefficients with shared PC/technical nuisance effects, gamma conditional-survival likelihood and MH age adjustment.',
        'The [statsmodels PHReg documentation](https://www.statsmodels.org/stable/generated/statsmodels.formula.api.phreg.html) documents delayed entry, strata and tie selection. [StratifiedTable](https://www.statsmodels.org/stable/generated/statsmodels.stats.contingency_tables.StratifiedTable.html) provides pooled MH odds ratios and confidence intervals. Cumulative baseline jumps are constructed explicitly for simulation, rather than using PHReg’s pre-jump baseline values directly.',
        '',
        'Counts 1–20 are not released. Stratified counts are withheld, and the entire route percentage stack is withheld if a small complementary category occurs. Numerical estimate/CI values are not participant counts. No commit, push or outgoing messages were performed.',
    ]
    (REPORT / '62_review_1008_reproduction.md').write_text('\n'.join(lines) + '\n')
    print('Validation and reproduction report saved.', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--task',
        choices=[
            'interaction',
            'dates',
            'cox',
            'frailty',
            'age-pattern',
            'family',
            'plots',
            'validate',
            'all',
        ],
        required=True,
    )
    parser.add_argument('--simulations', type=int, default=100)
    args = parser.parse_args()
    for path in [INTERNAL, TABLE, REPORT, FIGURE]:
        path.mkdir(parents=True, exist_ok=True)
    import subprocess

    if subprocess.run(['git', 'check-ignore', '-q', str(INTERNAL)], cwd=ROOT).returncode:
        raise RuntimeError('Internal participant cache must be ignored')
    with threadpool_limits(limits=1):
        data = None if args.task in ['plots', 'validate'] else load()
        if args.task in ['dates', 'all']:
            dates_task(data)
        if args.task in ['interaction', 'all']:
            interaction(data)
        if args.task in ['cox', 'all']:
            cox_task(data, args.simulations)
        if args.task in ['age-pattern', 'all']:
            age_pattern_task(data)
        if args.task in ['frailty', 'all']:
            frailty_task(data, args.simulations)
        if args.task in ['family', 'all']:
            family_task(data)
        if args.task in ['plots', 'all']:
            plots_task()
        if args.task == 'all':
            import subprocess

            subprocess.run([sys.executable, str(ROOT / 'tests/test_62_review_1008.py')], check=True)
        if args.task in ['validate', 'all']:
            validate_task()


if __name__ == '__main__':
    main()

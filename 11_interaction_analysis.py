"""Post hoc/exploratory interaction robustness, without plots or result constants."""
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import chi2
from scipy.linalg import qr
from patsy import dmatrix
from functions.models import covariates, logistic
from functions.validation import finite, standardize


def robust_wald(fit, indices):
    r = np.eye(len(fit.params))[indices]
    value = r @ fit.params
    covariance = r @ fit.cov_params() @ r.T
    degrees = np.linalg.matrix_rank(covariance)
    if not degrees:
        raise ValueError('Wald contrast is not estimable')
    statistic = float(value @ np.linalg.pinv(covariance) @ value)
    return {'statistic': statistic, 'df': int(degrees), 'p_value': float(chi2.sf(statistic, degrees))}


def spline(values, df=4):
    return np.asarray(dmatrix(f'bs(x, df={df}, degree=3, include_intercept=False) - 1', {'x': finite(values)}))


def interaction_robustness(cohort, raw_burden, gene_count, top_n=100):
    """Full-cohort scales stay fixed through subset sensitivity fits.

    Default top_n is the original exploratory influence/leverage exclusion size;
    it is a method parameter, never an observed sample count. Returns controlled
    models and Wald tests; no residual, influence or per-person file is written.
    """
    y = finite(cohort.is_case)
    raw = finite(raw_burden)
    p, _ = standardize(cohort.PRS311)
    b, _ = standardize(raw)
    c = covariates(cohort).to_numpy()
    base = np.column_stack([np.ones(len(y)), p, b, p * b, c])
    fits = {'continuous': logistic(y, base)}
    fitted = fits['continuous']
    mu = fitted.fittedvalues
    w = mu * (1 - mu)
    bread = np.linalg.inv(base.T * w @ base)
    leverage = w * np.einsum('ij,jk,ik->i', base, bread, base)
    influence = (y - mu) * (base @ bread[:, 3]) / (1 - leverage) / fitted.bse[3]
    if not 0 < top_n < len(y):
        raise ValueError('Invalid influence exclusion size')
    for name, values in [('influence', np.abs(influence)), ('leverage', leverage)]:
        mask = np.ones(len(y), dtype=bool)
        mask[np.argsort(values)[-top_n:]] = False
        fits['exclude_top_' + name] = logistic(y[mask], base[mask])
    carrier = (raw > 0).astype(float)
    fits['binary_carrier'] = logistic(y, np.column_stack([np.ones(len(y)), p, carrier, p * carrier, c]))
    if not carrier.any():
        raise ValueError('No carriers')
    for name, mask in [('single_gene', np.asarray(gene_count) <= 1),
                       ('below_carrier_p95', raw <= np.quantile(raw[carrier == 1], 0.95))]:
        # Conservative disclosure floor for running an optional subgroup fit.
        if (~mask).sum() > 20:
            fits[name] = logistic(y[mask], base[mask])
    sp = spline(p)
    sx = np.column_stack([np.ones(len(y)), sp, b, sp * b[:, None], c])
    fits['PRS_spline'] = logistic(y, sx)
    residual = sp - np.column_stack([np.ones(len(y)), p]) @ np.linalg.lstsq(np.column_stack([np.ones(len(y)), p]), sp, rcond=None)[0]
    _, _, pivot = qr(residual, mode='economic', pivoting=True)
    nl = residual[:, pivot[:3]]
    fits['nonlinear_PRS'] = logistic(y, np.column_stack([base, nl, nl * b[:, None]]))
    tests = {'spline_interactions': robust_wald(fits['PRS_spline'], list(range(6, 10))),
             'nonlinear_interactions': robust_wald(fits['nonlinear_PRS'], list(range(base.shape[1] + 3, base.shape[1] + 6)))}
    return fits, tests


def gene_interactions(cohort, gene_carriers):
    """Pool sparse panel genes; eligibility needs >20 carriers in each outcome.

    The inclusive-20 source boundary is intentionally tightened for disclosure.
    No per-gene observed carrier counts or result tables are included here.
    """
    if not gene_carriers.index.equals(cohort.index):
        raise ValueError('Gene-carrier order differs')
    y = finite(cohort.is_case)
    p, _ = standardize(cohort.PRS311)
    c = covariates(cohort).to_numpy()
    dense, sparse = [], []
    for gene in gene_carriers:
        x = finite(gene_carriers[gene])
        if not np.isin(x, [0, 1]).all():
            raise ValueError('Expected binary gene carriers')
        (dense if x[y == 1].sum() > 20 and x[y == 0].sum() > 20 else sparse).append(gene)
    indicators = gene_carriers[dense].copy()
    if sparse:
        pooled = gene_carriers[sparse].max(axis=1)
        if pooled[y == 1].sum() > 20 and pooled[y == 0].sum() > 20:
            indicators['other_pooled'] = pooled
    fits = {}
    for gene in indicators:
        indicator = indicators[gene].to_numpy()
        fits[gene] = logistic(y, np.column_stack([np.ones(len(y)), p, indicator, p * indicator, c]))
    if len(indicators.columns):
        g = indicators.to_numpy(dtype=float)
        fit = logistic(y, np.column_stack([np.ones(len(y)), p, g, g * p[:, None], c]))
        fits['joint'] = fit
        fits['omnibus'] = robust_wald(fit, list(range(2 + g.shape[1], 2 + 2 * g.shape[1])))
    return fits


def covariate_sensitivities(cohort, raw_burden, metadata):
    """Post hoc age specification and EHR coverage sensitivity fits.

    Full-cohort PRS/burden scales stay fixed. These are association sensitivities,
    not incident risk or corrected causal estimates. No observed results are
    exported and no specification is selected by its P value.
    """
    if not metadata.index.equals(cohort.index):
        raise ValueError('EHR metadata order differs')
    y = finite(cohort.is_case)
    p, _ = standardize(cohort.PRS311)
    b, _ = standardize(raw_burden)
    c = covariates(cohort).to_numpy()
    technical_and_pc = c[:, 1:]
    age = finite(cohort.age_at_last_obs)
    prefix = np.column_stack([np.ones(len(y)), p, b, p * b])
    baseline = np.column_stack([prefix, c])
    age_spline = np.column_stack([prefix, spline(age), technical_and_pc])
    results = {'nonlinear_observation_age': logistic(y, age_spline)}
    enrollment = metadata.age_at_enrollment.to_numpy(dtype=float)
    valid = np.isfinite(enrollment)
    linear = np.column_stack([prefix, np.nan_to_num(enrollment), technical_and_pc])
    results['enrollment_age'] = logistic(y[valid], linear[valid])
    window = valid & (enrollment >= 50) & (enrollment < 75)
    results['enrollment_50_to74_baseline'] = logistic(y[window], baseline[window])
    # Preserve the source whole-cohort spline construction followed by masking.
    nonlinear = np.column_stack([prefix, spline(np.nan_to_num(enrollment)), technical_and_pc])
    results['enrollment_50_to74_spline'] = logistic(y[window], nonlinear[window])
    span = metadata.ehr_span_days.to_numpy(dtype=float)
    dates = metadata.n_distinct_ehr_dates.to_numpy(dtype=float)
    complete = np.isfinite(span) & np.isfinite(dates) & (span >= 0) & (dates >= 0)
    log_ehr = np.column_stack([np.log1p(np.maximum(np.nan_to_num(span), 0)),
                                np.log1p(np.maximum(np.nan_to_num(dates), 0))])
    adjusted = np.column_stack([age_spline, log_ehr])
    coverage = complete & (span >= 365) & (dates >= 5)
    for name, x, mask in [('EHR_complete_baseline', baseline, complete),
                         ('EHR_complete_adjusted', adjusted, complete),
                         ('EHR_coverage_baseline', baseline, coverage),
                         ('EHR_coverage_adjusted', adjusted, coverage)]:
        results[name] = logistic(y[mask], x[mask])
    return results


def gene_interaction_summary(cohort, gene_carriers):
    """Joint adjusted gene interaction terms, omnibus/equality Wald and BH FDR.

    Coefficients are from the simultaneous gene-main-effect and gene-product
    model, not the optional one-gene marginal models. Output stays controlled.
    """
    from statsmodels.stats.multitest import multipletests
    all_fits = gene_interactions(cohort, gene_carriers)
    if 'joint' not in all_fits:
        raise ValueError('No estimable gene group')
    y = finite(cohort.is_case)
    names = [gene for gene in gene_carriers
             if gene_carriers.loc[y == 1, gene].sum() > 20
             and gene_carriers.loc[y == 0, gene].sum() > 20]
    sparse = [gene for gene in gene_carriers if gene not in names]
    if sparse:
        pooled = gene_carriers[sparse].max(axis=1)
        if pooled[y == 1].sum() > 20 and pooled[y == 0].sum() > 20:
            names.append('other_pooled')
    fit = all_fits['joint']
    indices = list(range(2 + len(names), 2 + 2 * len(names)))
    beta = np.asarray(fit.params)[indices]
    se = np.asarray(fit.bse)[indices]
    p = np.asarray(fit.pvalues)[indices]
    table = pd.DataFrame({'gene': names, 'beta': beta, 'se_HC1': se,
                          'interaction_OR': np.exp(beta), 'p_value': p,
                          'BH_q_value': multipletests(p, method='fdr_bh')[1]})
    table['ci_low'] = np.exp(beta - 1.959963984540054 * se)
    table['ci_high'] = np.exp(beta + 1.959963984540054 * se)
    equality = None
    if len(names) > 1:
        contrast = np.zeros((len(names) - 1, len(fit.params)))
        for row in range(1, len(names)):
            contrast[row - 1, indices[row]] = 1
            contrast[row - 1, indices[0]] = -1
        value = contrast @ fit.params
        covariance = contrast @ fit.cov_params() @ contrast.T
        degrees = np.linalg.matrix_rank(covariance)
        statistic = float(value @ np.linalg.pinv(covariance) @ value)
        equality = {'statistic': statistic, 'df': int(degrees),
                    'p_value': float(chi2.sf(statistic, degrees)) if degrees else np.nan}
    return table, {'omnibus': all_fits['omnibus'], 'equality': equality}

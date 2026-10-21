"""Paired five-model discrimination and case-only age association in memory.

No individual prediction, score, identifier or residual is printed or saved.
Returned fitted models and predictions remain controlled data.
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.metrics import roc_auc_score
from functions.validation import finite, standardize

PC = ['PC' + str(i) for i in range(1, 11)]
MODEL_NAMES = ('covariates', 'prs', 'burden', 'additive', 'interaction')


def covariates(frame):
    if frame[['age_at_last_obs', *PC, 'seq_site_id', 'seq_sequencer']].isna().any().any():
        raise ValueError('Incomplete model covariates; freeze membership before fitting')
    return pd.concat([frame[['age_at_last_obs', *PC]],
        pd.get_dummies(frame.seq_site_id, prefix='site', drop_first=True, dtype=float),
        pd.get_dummies(frame.seq_sequencer, prefix='seq', drop_first=True, dtype=float)], axis=1).astype(float)


def design(prs_z, burden_z, covariate_matrix, model='interaction'):
    if model not in MODEL_NAMES:
        raise ValueError('Unknown model')
    p, b, c = finite(prs_z), finite(burden_z), finite(covariate_matrix)
    columns = [np.ones(len(p))]
    if model in ('prs', 'additive', 'interaction'):
        columns.append(p)
    if model in ('burden', 'additive', 'interaction'):
        columns.append(b)
    if model == 'interaction':
        columns.append(p * b)
    return np.column_stack([*columns, c])


def logistic(y, x, robust=True):
    y, x = finite(y), finite(x)
    if not np.isin(y, [0, 1]).all() or len(np.unique(y)) != 2:
        raise ValueError('Expected both binary outcome groups')
    fitted = sm.GLM(y, x, family=sm.families.Binomial()).fit(**({'cov_type': 'HC1'} if robust else {}))
    if not fitted.converged or not np.isfinite(fitted.params).all() or not np.isfinite(fitted.bse).all():
        raise ValueError('Binomial model failed')
    return fitted


def coefficients(fitted, term_names, odds_ratio=False):
    if len(term_names) != len(fitted.params):
        raise ValueError('Coefficient names do not match model dimensions')
    beta = np.asarray(fitted.params)
    se = np.asarray(fitted.bse)
    z = stats.norm.ppf(0.975)
    table = pd.DataFrame({'term': term_names, 'beta': beta, 'se': se,
                          'ci_low': beta - z * se, 'ci_high': beta + z * se,
                          'p_value': np.asarray(fitted.pvalues)})
    if odds_ratio:
        table['OR'] = np.exp(beta)
        table['OR_ci_low'] = np.exp(table.ci_low)
        table['OR_ci_high'] = np.exp(table.ci_high)
    return table


def repeated_cv(y, prs, burdens, covariate_matrix, legacy_carrier,
                n_splits=5, n_repeats=10, seed=20260818, splits=None):
    """Shared case x legacy-carrier strata, with training-only sample-SD scaling.

    burdens is a dict of raw score versions on identical frozen members. Optional
    splits allow exact reuse of approved fold membership supplied privately.
    Technical dummy encoding follows the source full frozen-cohort encoding.
    """
    y = finite(y).astype(int)
    p, c = finite(prs), finite(covariate_matrix)
    legacy = np.asarray(legacy_carrier)
    if not np.isin(y, [0, 1]).all() or not np.isin(legacy, [0, 1]).all():
        raise ValueError('Invalid outcome or legacy carrier stratum')
    strata = np.asarray([str(a) + '_' + str(b) for a, b in zip(y, legacy)])
    if len(y) != len(p) or len(y) != len(c) or len(y) != len(legacy):
        raise ValueError('CV input dimensions differ')
    if splits is None:
        if pd.Series(strata).value_counts().min() < n_splits:
            raise ValueError('Insufficient stratum size for the fixed fold design')
        splits = list(RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats,
                                              random_state=seed).split(np.arange(len(y)), strata))
    if len(splits) != n_splits * n_repeats:
        raise ValueError('Unexpected number of fixed splits')
    predictions = {(version, model): np.full((n_repeats, len(y)), np.nan)
                   for version in burdens for model in MODEL_NAMES}
    seen = np.zeros((n_repeats, len(y)), dtype=int)
    for k, (train, test) in enumerate(splits):
        train, test = np.asarray(train, dtype=int), np.asarray(test, dtype=int)
        if len(set(train)) != len(train) or len(set(test)) != len(test) or set(train) & set(test):
            raise ValueError('Overlapping or duplicate fold indices')
        if set(train) | set(test) != set(range(len(y))):
            raise ValueError('Fold does not partition the frozen cohort')
        repeat = k // n_splits
        seen[repeat, test] += 1
        pt, _ = standardize(p[train], reference=p[train])
        pe, _ = standardize(p[test], reference=p[train])
        for version, raw in burdens.items():
            raw = finite(raw)
            if len(raw) != len(y):
                raise ValueError('Burden membership length differs')
            bt, _ = standardize(raw[train], reference=raw[train])
            be, _ = standardize(raw[test], reference=raw[train])
            for model in MODEL_NAMES:
                fit = logistic(y[train], design(pt, bt, c[train], model), robust=False)
                predictions[version, model][repeat, test] = fit.predict(design(pe, be, c[test], model))
    if not (seen == 1).all():
        raise ValueError('Each member must be held out once per repeat')
    rows = []
    for (version, model), probability in predictions.items():
        for repeat in range(n_repeats):
            auc = roc_auc_score(y, probability[repeat])
            baseline = roc_auc_score(y, predictions[version, 'prs'][repeat])
            rows.append({'version': version, 'model': model, 'repeat': repeat + 1,
                         'auc': auc, 'paired_delta_vs_prs': auc - baseline})
    # Predictions are returned to the caller in memory only; never auto-exported.
    return pd.DataFrame(rows), predictions


def age_models(age, prs, burden, pcs, covariance='HC1', interaction=False):
    if covariance not in {'HC1', 'HC3'}:
        raise ValueError('Unsupported covariance')
    p, pscale = standardize(prs)
    b, bscale = standardize(burden)
    x = np.column_stack([np.ones(len(p)), p, b, *([p * b] if interaction else []), finite(pcs)])
    fit = sm.OLS(finite(age), x).fit(cov_type=covariance)
    if not np.isfinite(fit.params).all() or not np.isfinite(fit.bse).all():
        raise ValueError('Age model failed')
    return fit, {'prs': pscale, 'burden': bscale}

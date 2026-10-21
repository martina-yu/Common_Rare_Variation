"""Binomial HC1 fits and DeLong placements; all returned arrays are controlled."""
import numpy as np
from scipy import stats
import statsmodels.api as sm

def delong_auc_variance(y_true, y_score):
    """Delong auc variance."""
    y_true = np.asarray(y_true)
    y_score = np.asarray(y_score)
    pos = y_score[y_true == 1]
    neg = y_score[y_true == 0]
    n_pos, n_neg = (len(pos), len(neg))

    def midrank(x):
        return stats.rankdata(x, method='average')
    combined = np.concatenate([pos, neg])
    all_ranks = midrank(combined)
    pos_ranks_in_combined = all_ranks[:n_pos]
    tx = midrank(pos)
    ty = midrank(neg)
    v10 = (pos_ranks_in_combined - tx) / n_neg
    neg_ranks_in_combined = all_ranks[n_pos:]
    v01 = 1.0 - (neg_ranks_in_combined - ty) / n_pos
    auc = v10.mean()
    return (auc, v10, v01)

def delong_test(y_true, score_a, score_b):
    """Delong test."""
    auc_a, v10_a, v01_a = delong_auc_variance(y_true, score_a)
    auc_b, v10_b, v01_b = delong_auc_variance(y_true, score_b)
    n_pos, n_neg = (len(v10_a), len(v01_a))
    v10 = np.vstack([v10_a, v10_b])
    v01 = np.vstack([v01_a, v01_b])
    s10 = np.cov(v10, bias=False) / n_pos
    s01 = np.cov(v01, bias=False) / n_neg
    cov_mat = s10 + s01
    delta = auc_a - auc_b
    var_delta = cov_mat[0, 0] + cov_mat[1, 1] - 2 * cov_mat[0, 1]
    se_delta = np.sqrt(max(var_delta, 0))
    z = delta / se_delta if se_delta > 0 else np.nan
    p = 2 * stats.norm.sf(abs(z)) if se_delta > 0 else np.nan
    ci_low = delta - 1.959963984540054 * se_delta
    ci_high = delta + 1.959963984540054 * se_delta
    return {'auc_a': auc_a, 'auc_b': auc_b, 'delta_auc': delta, 'ci_low': ci_low, 'ci_high': ci_high, 'p_value': p}

def auc_ci(y, score):
    auc, v10, v01 = delong_auc_variance(y, score)
    se = np.sqrt(np.var(v10, ddof=1) / len(v10) + np.var(v01, ddof=1) / len(v01))
    return (auc - 1.959963984540054 * se, auc + 1.959963984540054 * se)

def fit_binomial(yv, X, robust=True):
    fit = sm.GLM(yv, X, family=sm.families.Binomial()).fit(**{'cov_type': 'HC1'} if robust else {})
    assert fit.converged and np.isfinite(fit.params).all() and np.isfinite(fit.bse).all()
    return fit

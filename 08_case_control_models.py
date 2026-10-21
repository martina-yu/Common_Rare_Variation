"""Full-cohort HC1 logistic fits and paired five-model repeated CV.

Comparison/sensitivity scores must be specified before fitting. All versions
use the same frozen members. No outcome-based scheme selection occurs here.
"""
import pandas as pd
from functions.models import covariates, design, logistic, repeated_cv, PC
from functions.validation import aligned, unique_ids, standardize


def fit_models(cohort, burden_scores):
    """cohort and burden_scores have identical unique person_id indices."""
    burden_scores = aligned(burden_scores, cohort.index)
    c = covariates(cohort)
    p, _ = standardize(cohort.PRS311)
    fits = {}
    for version in burden_scores:
        b, _ = standardize(burden_scores[version])
        for model in ('covariates', 'prs', 'burden', 'additive', 'interaction'):
            fits[version, model] = logistic(cohort.is_case, design(p, b, c, model))
    return fits


def evaluate_models(cohort, burden_scores, legacy_carrier, **cv_parameters):
    scores = aligned(burden_scores, cohort.index)
    if not legacy_carrier.index.equals(cohort.index):
        raise ValueError('Legacy stratum order differs from the fixed cohort')
    return repeated_cv(cohort.is_case, cohort.PRS311, {name: scores[name].to_numpy() for name in scores},
                       covariates(cohort), legacy_carrier.to_numpy(), **cv_parameters)

"""Within-case recorded diagnosis-age association, not incident absolute risk.

Additive models use HC1; continuous burden x PRS interaction uses HC3 primarily
and HC1 as a covariance comparison. Main effects in the interaction model are
conditional on the other standardized score being zero.
"""
from functions.models import age_models, PC
from functions.validation import aligned


def fit_age_models(age_cases, burden_scores):
    scores = aligned(burden_scores, age_cases.index)
    if not age_cases.is_case.eq(1).all() or age_cases.bc_is_first_record.eq(True).any():
        raise ValueError('Invalid frozen age-model membership')
    if age_cases.age_at_diagnosis.lt(18).any():
        raise ValueError('Invalid diagnosis age')
    results = {}
    for version in scores:
        for interaction, covariance in ((False, 'HC1'), (True, 'HC3'), (True, 'HC1')):
            results[version, interaction, covariance] = age_models(age_cases.age_at_diagnosis,
                age_cases.PRS311, scores[version], age_cases[PC], covariance, interaction)
    return results

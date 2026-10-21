"""Score definitions to freeze before fitting; no significance-driven selection.

Inherited sensitivity definitions are not a claim of prospective preregistration.
Min1star carriers require a separate full ClinVar/disease/inheritance review at
that threshold, not relabeling the main >=2-star scores.
"""
import numpy as np
import pandas as pd
from functions.validation import aligned


def sensitivity_scores(main_union, ptv_only, strict_vaf_union, min1star_union,
                       odds_ratios, include_without_cdh1=False):
    index = main_union.index
    strict = aligned(strict_vaf_union, index)
    relaxed = aligned(min1star_union, index)
    ptv = aligned(ptv_only, index)
    definitions = {'without_STK11': (main_union, [g for g in odds_ratios if g != 'STK11']),
        'PTV_only': (ptv, list(odds_ratios)),
        'TP53_NF1_strict_VAF': (strict, list(odds_ratios)),
        'ClinVar_min1star': (relaxed, list(odds_ratios))}
    if include_without_cdh1:
        definitions['without_CDH1'] = (main_union, [g for g in odds_ratios if g != 'CDH1'])
    scores = pd.DataFrame(index=index)
    for name, (frame, genes) in definitions.items():
        matrix = frame[genes].to_numpy(dtype=float)
        weights = np.asarray([odds_ratios[g] for g in genes], dtype=float)
        if not np.isin(matrix, [0, 1]).all() or not np.isfinite(weights).all() or (weights <= 0).any():
            raise ValueError('Invalid carriers or external weights')
        scores[name] = matrix @ np.log(weights)
    return scores


def ehr_sensitivity_masks(metadata):
    span, dates = metadata.ehr_span_days, metadata.n_distinct_ehr_dates
    complete = span.notna() & dates.notna() & span.ge(0) & dates.ge(0)
    return {'metadata_complete': complete, 'span_ge365_dates_ge5': complete & span.ge(365) & dates.ge(5)}

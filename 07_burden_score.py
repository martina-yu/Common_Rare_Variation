"""Burden scoring API. Stored-score reference ddof=0 differs from model ddof=1.

Main weights approximate the published PTV+non-PTV PV combined ORs; the study
missense definition and reference variant classes/ClinVar dates are not identical.
"""
import numpy as np
import pandas as pd
from functions.burden import read_weights, score_panel
from functions.validation import aligned, same_membership


def weighted_gene_score(carriers, odds_ratios, reference_parameters=None):
    genes = list(odds_ratios)
    x = carriers[genes].to_numpy(dtype=float)
    weights = np.asarray([odds_ratios[gene] for gene in genes], dtype=float)
    if not np.isin(x, [0, 1]).all() or not np.isfinite(weights).all() or (weights <= 0).any():
        raise ValueError('Expected binary gene indicators and positive finite ORs')
    raw = x @ np.log(weights)
    reference = ({'mean': float(raw.mean()), 'SD': float(raw.std(ddof=0)), 'ddof': 0}
                 if reference_parameters is None else dict(reference_parameters))
    if reference['ddof'] != 0 or not np.isfinite(reference['SD']) or reference['SD'] <= 0:
        raise ValueError('Invalid saved-score reference scale')
    return pd.DataFrame({'B_raw': raw, 'B_standardized': (raw - reference['mean']) / reference['SD'],
                         'carried_gene_count': x.sum(axis=1)}, index=carriers.index), reference


def comparison_scores(legacy_carriers, current_ptv, current_union,
                      legacy_or, combined_or, main_genes):
    """Separate reweighting, updated PTV QC, missense and panel extension.

    A: legacy nine-gene PTVs + original weights.
    A_prime: same exact legacy carriers + combined weights.
    A_prime_QC: current nine-gene PTVs + combined weights.
    B: nine-gene PTV/missense union + combined weights.
    C: main panel union + combined weights.
    Missing AF legacy behavior is retained only through supplied legacy carriers;
    this API does not recreate their historical membership from a new definition.
    """
    nine = list(legacy_or)
    current_ptv = aligned(current_ptv, legacy_carriers.index)
    current_union = aligned(current_union, legacy_carriers.index)
    definitions = {'A': (legacy_carriers, legacy_or),
        'A_prime': (legacy_carriers, {g: combined_or[g] for g in nine}),
        'A_prime_QC': (current_ptv, {g: combined_or[g] for g in nine}),
        'B': (current_union, {g: combined_or[g] for g in nine}),
        'C': (current_union, {g: combined_or[g] for g in main_genes})}
    return {name: weighted_gene_score(frame, weights) for name, (frame, weights) in definitions.items()}

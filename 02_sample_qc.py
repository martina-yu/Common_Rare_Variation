"""EUR, greedy kinship pruning, missingness proxy, control-only HWE and covariates.

The missingness/HWE proxy is three chr21 windows, not whole-genome QC.
Greedy maximum-degree deletion produces an unrelated set, not a proven maximum.
"""
import numpy as np
import pandas as pd
import networkx as nx
from functions.validation import unique_ids


def kinship_exclusions(cohort_ids, pairs, threshold=0.1):
    """Pair input order controls ties, as in the source greedy implementation.

    Available relatedness edges above 0.1 cannot identify weak pairs in the
    0.0884--0.1 interval; 0.1 is a less stringent approximation to 0.0884.
    """
    keep = set(cohort_ids)
    selected = pairs.loc[pairs.id1.isin(keep) & pairs.id2.isin(keep)
                         & pairs.kinship.gt(threshold)]
    graph = nx.Graph()
    graph.add_edges_from(zip(selected.id1, selected.id2))
    excluded = set()
    while graph.number_of_edges():
        degrees = dict(graph.degree())
        node = max(degrees, key=degrees.get)
        excluded.add(node)
        graph.remove_node(node)
    return excluded


def proxy_variant_indices(n_variants):
    windows = (5000, n_variants // 2, n_variants - 5500)
    indices = np.concatenate([np.arange(start, start + 500) for start in windows])
    if len(set(indices)) != len(indices) or indices.min() < 0 or indices.max() >= n_variants:
        raise ValueError('Proxy windows overlap or exceed the source')
    return indices


def missingness_flags(dosage):
    """Input is variants x cohort samples, diploid allele dosage with NaN missing.

    Exclude proxy variants with missing rate >0.2 before applying sample >0.03.
    The same variant exclusion list is used for proxy HWE.
    """
    x = np.asarray(dosage, dtype=float)
    if x.ndim != 2 or not np.isin(x[~np.isnan(x)], [0, 1, 2]).all():
        raise ValueError('Expected diploid dosage matrix')
    retained = np.isnan(x).mean(axis=1) <= 0.2
    if not retained.any():
        raise ValueError('No usable proxy markers')
    rate = np.isnan(x[retained]).mean(axis=0)
    return rate, rate > 0.03, retained


def hwe_command(plink2, bed_prefix, control_keep, excluded_variants, output_prefix):
    """Run only inside authorized storage; PLINK outputs are controlled intermediates.

    Use --hardy P<1e-6 as the proxy marker failure criterion. The source reports
    this proxy diagnostic and does not establish genome-wide HWE filtering.
    """
    return [plink2, '--bfile', str(bed_prefix), '--keep', str(control_keep),
            '--exclude', str(excluded_variants), '--hardy', '--out', str(output_prefix)]


def analysis_membership(cohort, excluded_ids, missingness):
    unique_ids(cohort)
    unique_ids(missingness)
    d = cohort.merge(missingness[['person_id', 'f_miss']], on='person_id', how='left', validate='one_to_one')
    valid_case = d.group.eq('case') & d.age_at_diagnosis.ge(18)
    eligible = (valid_case | d.group.eq('control')) & d.is_eur.eq(True)
    eligible &= ~d.person_id.isin(excluded_ids) & d.f_miss.notna() & d.f_miss.le(0.03)
    qc = d.loc[eligible].copy()
    model_columns = ['PRS311', 'age_at_last_obs', *['PC' + str(i) for i in range(1, 11)], 'seq_site_id', 'seq_sequencer']
    # This is an explicit model-completeness step, separate from QC membership.
    model = qc.dropna(subset=model_columns).copy()
    age = model.loc[model.group.eq('case') & model.bc_is_first_record.ne(True)
                    & model.age_at_diagnosis.notna()].copy()
    return {'qc': qc, 'case_control': model, 'diagnosis_age': age}

def map_menopause(ans):
    if pd.isna(ans):
        return np.nan
    if 'Yes None' in ans or 'Yes, no menstrual periods' in ans:
        return 'postmenopausal'
    if 'Periods Havent Stopped' in ans or ans.strip().endswith(': No'):
        return 'premenopausal'
    if 'Yes But Hormone' in ans or 'hormone' in ans.lower():
        return 'uncertain_hormone_influenced'
    if 'Not Sure' in ans:
        return 'unknown_not_sure'
    if 'Skip' in ans or 'Prefer Not' in ans or 'Dont Know' in ans:
        return np.nan
    return np.nan

def derive_smoking(row):
    a1 = row['ans_100cigs']
    a2 = row['ans_freq']
    if pd.isna(a1) or 'Skip' in str(a1) or 'Prefer Not' in str(a1) or ('Dont Know' in str(a1)):
        return np.nan
    if a1.endswith(': No'):
        return 'never_smoker'
    if a1.endswith(': Yes'):
        if pd.isna(a2):
            return 'ever_smoker_currency_unknown'
        if 'Not At All' in a2:
            return 'former_smoker'
        if 'Every Day' in a2 or 'Some Days' in a2:
            return 'current_smoker'
        return 'ever_smoker_currency_unknown'
    return np.nan


def family_history_label(person_id, positive_nonself_ids, answered_ids):
    if person_id in positive_nonself_ids:
        return 'positive_nonself'
    return 'negative_or_self_only' if person_id in answered_ids else np.nan

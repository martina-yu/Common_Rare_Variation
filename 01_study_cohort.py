"""Female WGS cohort, genetic ancestry and recorded EHR phenotypes.

Inputs are in-memory tables from an authorized environment. No rows are printed
or exported. Queries are supplied separately in functions/queries.py.
"""
import ast
import numpy as np
import pandas as pd
from functions.validation import unique_ids

FEMALE_CONCEPT = 45878463
BREAST_CANCER_CONCEPT = 4112853
ANCESTRY_CLASSES = ('afr', 'amr', 'eas', 'eur', 'mid', 'sas')


def ancestry_table(records, keep_ids):
    """Parse supplied ancestry rows; validate class order and posterior semantics.

    research_id must use the same identifier representation as keep_ids.
    The six-class order and 16 PC fields are release-specific assumptions to
    validate on the chosen release; probability rounding occurs after labeling.
    """
    rows = []
    keep = {str(x) for x in keep_ids}
    for row in records.itertuples(index=False):
        probabilities = ast.literal_eval(row.probabilities) if isinstance(row.probabilities, str) else row.probabilities
        pcs = ast.literal_eval(row.pca_features) if isinstance(row.pca_features, str) else row.pca_features
        if len(probabilities) != 6 or len(pcs) != 16:
            raise ValueError('Unexpected ancestry vector dimensions')
        p = np.asarray(probabilities, dtype=float)
        if not np.isfinite(p).all() or (p < 0).any() or not np.isclose(p.sum(), 1, atol=1e-5):
            raise ValueError('Invalid ancestry probability vector')
        label = ANCESTRY_CLASSES[int(np.argmax(p))]
        other = label if max(p) > 0.75 else 'oth'
        if row.ancestry_pred != label or row.ancestry_pred_other != other:
            raise ValueError('Ancestry class order or threshold semantics changed')
        if str(row.research_id) in keep:
            rows.append([str(row.research_id), other, *np.round(p, 6), *pcs])
    result = pd.DataFrame(rows, columns=['person_id', 'ancestry_pred_other',
        *['p_' + c for c in ANCESTRY_CLASSES], *['PC' + str(i) for i in range(1, 17)]])
    unique_ids(result)
    if set(result.person_id) != keep:
        raise ValueError('Incomplete ancestry coverage')
    result['is_eur'] = result.ancestry_pred_other.eq('eur')
    return result


def select_female_wgs(persons, wgs_ids, ancestry):
    unique_ids(persons)
    unique_ids(ancestry)
    selected = persons.loc[persons.sex_at_birth_concept_id.eq(FEMALE_CONCEPT)
                           & persons.person_id.isin(wgs_ids)].copy()
    result = selected.merge(ancestry, on='person_id', how='left', validate='one_to_one')
    if result.is_eur.isna().any():
        raise ValueError('Incomplete ancestry coverage')
    return result


def derive_phenotypes(universe, diagnosis_dates, ehr_dates, death_dates,
                      consent, diagnosis_counts, any_cancer_ids):
    """One BC date defines a case; >=2 distinct dates is a sensitivity flag.

    diagnosis_dates: person_id, dx_date. ehr_dates: person_id, event_date,
    the union of condition/procedure/drug/measurement/visit dates. death_dates
    has one resolved earliest COALESCE(primary AoU death, OMOP death) per person.
    Controls have no malignant-neoplasm ancestor diagnosis; other cancers are
    excluded from controls. Diagnosis age below 18 is invalidated.
    """
    unique_ids(universe)
    dx = diagnosis_dates.copy()
    dx['dx_date'] = pd.to_datetime(dx.dx_date, utc=True).dt.tz_convert(None).dt.normalize()
    dx = dx.groupby('person_id').dx_date.agg(first_bc_dx_date='min', n_distinct_bc_dx_dates='nunique').reset_index()
    ehr = ehr_dates.copy()
    ehr['event_date'] = pd.to_datetime(ehr.event_date, utc=True).dt.tz_convert(None).dt.normalize()
    ehr = ehr.groupby('person_id').event_date.agg(ehr_first_date='min', ehr_last_date='max', n_distinct_ehr_dates='nunique').reset_index()
    result = universe.copy()
    for table in (dx, ehr, death_dates, consent, diagnosis_counts):
        unique_ids(table)
        result = result.merge(table, on='person_id', how='left', validate='one_to_one')
    result['is_case'] = result.first_bc_dx_date.notna()
    result['is_case_rule_of_two'] = result.n_distinct_bc_dx_dates.ge(2)
    result['is_control'] = ~result.is_case & ~result.person_id.isin(any_cancer_ids)
    result['group'] = np.select([result.is_case, result.is_control], ['case', 'control'], default='other_cancer_excluded')
    for column in ('birth_datetime', 'death_date', 'ehr_last_date', 'ehr_first_date', 'first_bc_dx_date'):
        result[column] = pd.to_datetime(result[column], utc=True).dt.tz_convert(None).dt.normalize()
    result['last_obs_date'] = result[['ehr_last_date', 'death_date']].max(axis=1)
    for output, date in [('age_at_diagnosis', 'first_bc_dx_date'), ('age_at_last_obs', 'last_obs_date')]:
        result[output] = (result[date] - result.birth_datetime).dt.days / 365.25
    result['flag_age_at_diagnosis_invalid_under18'] = result.age_at_diagnosis.lt(18)
    result.loc[~result.is_case | result.flag_age_at_diagnosis_invalid_under18, 'age_at_diagnosis'] = np.nan
    result['age_at_enrollment'] = result.age_at_consent.astype(float)
    result['bc_is_first_record'] = result.first_bc_dx_date.eq(result.ehr_first_date) & result.is_case
    result['ehr_span_days'] = (result.ehr_last_date - result.ehr_first_date).dt.days
    result['dx_code_occurrence_count'] = result.dx_code_occurrence_count.fillna(0).astype(int)
    return result

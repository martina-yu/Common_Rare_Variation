"""Read-only Table 1 audit; export aggregates, never participant records."""

from pathlib import Path
import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.functions.paths import PROJECT_ROOT

import csv
import re
import json
import pandas as pd

BASE = PROJECT_ROOT
OUT = BASE / 'results/table'
table = pd.read_csv(OUT / '09b_supp_table1_cohort.csv', keep_default_na=False)
cols = list(table.columns[1:])
denoms = {c: int(table.iloc[0][c].replace(',', '')) for c in cols}
audit = []
section = ''
for _, row in table.iloc[1:].iterrows():
    label = row['variable'].strip()
    if all(row[c] == '' for c in cols):
        section = label
    for c in cols:
        match = re.fullmatch(r'([\d,]+) \(([\d.]+)\)', row[c])
        if not match or ('n (%)' not in label and not row['variable'].startswith('  ')):
            continue
        n, printed = int(match[1].replace(',', '')), float(match[2])
        calculated = n / denoms[c] * 100
        audit.append(
            dict(
                variable=label,
                section=section if row['variable'].startswith('  ') else '',
                cohort=c,
                numerator=n,
                denominator=denoms[c],
                formula=f'100 * {n} / {denoms[c]}',
                calculated_percent=round(calculated, 6),
                printed_percent=printed,
                rounded_percent=float(f'{calculated:.1f}'),
                arithmetic_matches=float(f'{calculated:.1f}') == printed,
            )
        )
pd.DataFrame(audit).to_csv(OUT / '09b_table1_percentage_audit.csv', index=False)

# Participant records remain in memory; only disclosure-safe summaries leave it.
qc = pd.read_csv(BASE / 'data/05_cohort_qc_flags.csv')
qc = qc.loc[qc.pass_qc.eq(1), ['person_id', 'analysis_group']]
ph = pd.read_csv(BASE / 'data/04_phenotype_diagnostics_person_level.csv.gz')
cl = pd.read_csv(BASE / 'data/09_extended_clinical_person_level.csv.gz').drop(columns='group')
prs = pd.read_csv(BASE / 'data/prs_scores.csv', usecols=['person_id', 'PRS311'])
d = (
    qc.merge(ph.drop(columns='group'), on='person_id', validate='one_to_one')
    .merge(cl, on='person_id', validate='one_to_one')
    .merge(prs, on='person_id', validate='one_to_one')
)
case = d.analysis_group.eq('case')
masks = {
    cols[0]: case,
    cols[1]: d.analysis_group.eq('control'),
    cols[2]: case & d.bc_is_first_record.ne(True),
}
d['menopause'] = d.menopause_status.map(
    {
        'premenopausal': 'Premenopausal',
        'postmenopausal': 'Postmenopausal',
        'unknown_not_sure': 'Uncertain or unknown',
        'uncertain_hormone_influenced': 'Uncertain or unknown',
    }
)
d['family_history'] = d.family_history_breast_cancer.map(
    {
        'positive_nonself': 'Positive (first- or second-degree relative)',
        'negative_or_self_only': 'Negative',
    }
)
d['smoking'] = d.smoking_status.map(
    {
        'never_smoker': 'Never',
        'former_smoker': 'Former or unknown currency',
        'ever_smoker_currency_unknown': 'Former or unknown currency',
        'current_smoker': 'Current',
    }
)
d['seq_site'] = d.seq_site_id.map(
    {'bi': 'Broad Institute', 'uw': 'University of Washington', 'bcm': 'Baylor College of Medicine'}
)
mapping = {
    'Menopausal status, n (%)': 'menopause',
    'Family history of breast cancer, n (%)': 'family_history',
    'Smoking status, n (%)': 'smoking',
    'Sequencing center, n (%)': 'seq_site',
    'Sequencer, n (%)': 'seq_sequencer',
}
cont = [
    ('Age at diagnosis, years, mean (SD)', 'age_at_diagnosis'),
    ('Age at enrollment, years, mean (SD)', 'age_at_enrollment'),
    ('Age at last EHR observation, years, mean (SD)', 'age_at_last_obs'),
    ('BMI, kg/m2, mean (SD)', 'bmi_value'),
    ('PRS311 (raw score), mean (SD)', 'PRS311'),
]
checks = []
for c, m in masks.items():
    checks.append({'check': 'cohort N', 'cohort': c, 'match': int(m.sum()) == denoms[c]})
    for label, v in cont:
        if v == 'age_at_diagnosis' and c == cols[1]:
            continue
        x = d.loc[m, v].dropna()
        value = f'{x.mean():.1f} ({x.std():.1f})'
        original = table.loc[table.variable.eq(label), c].iloc[0]
        checks.append({'check': label, 'cohort': c, 'match': value == original})
    if c != cols[1]:
        x = d.loc[m, 'age_at_diagnosis'].dropna()
        value = f'{x.median():.1f} [{x.quantile(.25):.1f}-{x.quantile(.75):.1f}]'
        checks.append(
            {
                'check': 'diagnosis median/IQR',
                'cohort': c,
                'match': value
                == table.loc[table.variable.str.contains('median', regex=False), c].iloc[0],
            }
        )
    sec = ''
    for _, row in table.iterrows():
        if all(row[k] == '' for k in cols):
            sec = row.variable
        match = re.fullmatch(r'([\d,]+) \(([\d.]+)\)', row[c])
        if not match:
            continue
        if row.variable.startswith('  ') and sec in mapping:
            v = mapping[sec]
            label = row.variable.strip()
            hit = d[v].isna() if label == 'Missing' else d[v].eq(label)
        elif row.variable == 'BMI missing, n (%)':
            hit = d.bmi_value.isna()
        elif row.variable == 'Age at last EHR observation missing, n (%)':
            hit = d.age_at_last_obs.isna()
        else:
            continue
        checks.append(
            {
                'check': sec + ': ' + row.variable.strip(),
                'cohort': c,
                'match': int((m & hit).sum()) == int(match[1].replace(',', '')),
            }
        )
pd.DataFrame(checks).to_csv(OUT / '09b_table1_source_reproduction_audit.csv', index=False)

# Broad menopause/age groups are all safely above the reporting threshold.
age_rows = []
for c, m in masks.items():
    for status in ['Premenopausal', 'Postmenopausal']:
        for v in ['age_at_enrollment', 'age_at_last_obs']:
            x = d.loc[m & d.menopause.eq(status), v].dropna()
            if len(x) < 20:
                continue
            age_rows.append(
                {
                    'cohort': c,
                    'menopause_status': status,
                    'age_variable': v,
                    'n_observed': len(x),
                    'mean': round(x.mean(), 3),
                    'sd': round(x.std(), 3),
                    'median': round(x.median(), 3),
                    'q25': round(x.quantile(0.25), 3),
                    'q75': round(x.quantile(0.75), 3),
                }
            )
pd.DataFrame(age_rows).to_csv(OUT / '09b_table1_menopause_age_audit.csv', index=False)
total_rows = []
for section, v in mapping.items():
    for c, m in masks.items():
        total_rows.append(
            {
                'variable': section,
                'cohort': c,
                'category_count_sum': int(d.loc[m, v].value_counts(dropna=False).sum()),
                'cohort_n': denoms[c],
            }
        )
pd.DataFrame(total_rows).to_csv(OUT / '09b_table1_category_total_audit.csv', index=False)
age_flags = []
for c in cols[:2]:
    m = masks[c] & d.menopause.eq('Premenopausal')
    for cutoff in [55, 60]:
        count = int((m & d.age_at_enrollment.ge(cutoff)).sum())
        safe = str(count) if count >= 20 or count == 0 else ('<5' if count < 5 else '5-19')
        age_flags.append(
            {
                'cohort': c,
                'check': f'Premenopausal with enrollment age >= {cutoff}',
                'count': safe,
                'interpretation': 'Review against age at survey; enrollment age is only a proxy.',
            }
        )
pd.DataFrame(age_flags).to_csv(OUT / '09b_table1_age_flags.csv', index=False)
summary = {
    'published_percentages_checked': len(audit),
    'percentage_mismatches': sum(not a['arithmetic_matches'] for a in audit),
    'source_checks': len(checks),
    'source_mismatches': [a for a in checks if not a['match']],
    'age_rows': age_rows,
    'age_flags': age_flags,
    'diagnosis_age_complete_in_age_cohort': bool(
        d.loc[masks[cols[2]], 'age_at_diagnosis'].notna().all()
    ),
}
print(json.dumps(summary, indent=2))

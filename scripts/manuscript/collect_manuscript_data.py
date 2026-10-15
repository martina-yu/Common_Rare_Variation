"""Collect unrounded aggregate C results; index restricted scores in place."""

from pathlib import Path
import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.functions.paths import PROJECT_ROOT

import csv
import hashlib
import json
import shutil

ROOT = PROJECT_ROOT
OUT = ROOT / 'manuscript/datas'


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    records = []

    def collect(source, destination, purpose, main_only=False):
        src = ROOT / source
        dst = OUT / destination
        dst.parent.mkdir(parents=True, exist_ok=True)
        with src.open(newline='') as f:
            reader = csv.DictReader(f)
            fields = reader.fieldnames
            rows = list(reader)
        assert fields and not {'person_id', 'sample_id', 'normalized_key'} & set(fields)
        if main_only:
            rows = [r for r in rows if r.get('version') in {'C', 'PRS311_baseline'}]
            assert rows
            with dst.open('w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=fields)
                writer.writeheader()
                writer.writerows(rows)
        else:
            shutil.copy2(src, dst)
            assert src.read_bytes() == dst.read_bytes()
        with dst.open(newline='') as f:
            assert list(csv.DictReader(f)) == rows
        records.append(
            {
                'file': destination,
                'source': source,
                'purpose': purpose,
                'selection': 'C and PRS311_baseline only' if main_only else 'unchanged source',
                'sha256': hashlib.sha256(dst.read_bytes()).hexdigest(),
            }
        )

    collect(
        'results/table/29_rowlands_13gene_external_weights.csv',
        'burden/gene_OR_weights.csv',
        'External OR and ln(OR), 13-gene panel',
    )
    collect(
        'results/table/45_rowlands_fixed_cohort_coefficients.csv',
        'case_control/joint_model_coefficients.csv',
        'C joint logistic and age coefficient estimates',
        True,
    )
    # Split the main coefficients by outcome, preserving literal numeric cells.
    mixed = OUT / 'case_control/joint_model_coefficients.csv'
    with mixed.open(newline='') as f:
        reader = csv.DictReader(f)
        fields = reader.fieldnames
        rows = list(reader)
    for cohort, filename in [
        ('case_control', mixed),
        ('diagnosis_age', OUT / 'age/age_model_coefficients.csv'),
    ]:
        filename.parent.mkdir(parents=True, exist_ok=True)
        with filename.open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            writer.writerows(r for r in rows if r['cohort'] == cohort)
    records[1]['purpose'] = 'C joint logistic coefficients: beta, SE, OR, CI, p-value'
    records[1]['selection'] = 'version C, cohort case_control'
    records[1]['sha256'] = hashlib.sha256(mixed.read_bytes()).hexdigest()
    records.append(
        {
            'file': 'age/age_model_coefficients.csv',
            'source': 'results/table/45_rowlands_fixed_cohort_coefficients.csv',
            'selection': 'version C, cohort diagnosis_age',
            'purpose': 'Age OLS beta (years/SD), SE, CI, p-value; OR not applicable',
            'sha256': hashlib.sha256(
                (OUT / 'age/age_model_coefficients.csv').read_bytes()
            ).hexdigest(),
        }
    )
    collect(
        'results/table/45_rowlands_CV_AUC.csv',
        'case_control/CV_AUC_summary.csv',
        'PRS baseline and C joint 10x5 OOF AUC; repeat1 CI',
        True,
    )
    collect(
        'results/table/45_rowlands_CV_delta_AUC.csv',
        'case_control/CV_delta_AUC_summary.csv',
        'C joint-minus-PRS comparison',
        True,
    )
    collect(
        'results/table/45_rowlands_age_explained_variance.csv',
        'age/age_R_squared.csv',
        'C age R-squared and incremental R-squared',
        True,
    )
    collect(
        'results/table/48_age_binary_carrier_PRS_slopes.csv',
        'age/carrier_PRS_slope_contrasts.csv',
        'Exploratory binary carrier age interaction and PRS slopes',
    )
    collect(
        'results/table/48_burden_distribution_carrier_supplement_coarsened.csv',
        'burden/score_distribution_disclosure_coarsened.csv',
        'Disclosure-coarsened published score distribution; not unrounded individual scores',
    )
    collect(
        'results/table/45_rowlands_fixed_cohort_coefficients.csv',
        'comparisons/all_versions_coefficients.csv',
        'Legacy, bridging and new sensitivity models, explicit version labels',
    )
    collect(
        'results/table/45_rowlands_CV_AUC.csv',
        'comparisons/all_versions_CV_AUC_summary.csv',
        'Version comparison and sensitivity AUC',
    )
    collect(
        'results/table/45_rowlands_CV_delta_AUC.csv',
        'comparisons/all_versions_CV_delta_AUC_summary.csv',
        'Version comparison and sensitivity delta AUC',
    )
    collect(
        'results/table/45_rowlands_age_explained_variance.csv',
        'comparisons/all_versions_age_R_squared.csv',
        'Version comparison and sensitivity age R-squared',
    )
    for src, dst, purpose in [
        (
            'CV_AUC_by_repeat_internal.csv',
            'case_control/CV_AUC_by_repeat.csv',
            'Unrounded aggregate AUC by repeat, C and PRS baseline',
        ),
        (
            'CV_paired_deltas_by_repeat_internal.csv',
            'case_control/CV_delta_AUC_by_repeat.csv',
            'Unrounded aggregate paired delta AUC by repeat for version C',
        ),
    ]:
        collect('data/29_rowlands_extension/section11_models/' + src, dst, purpose, True)

    restricted = [
        {
            'score': 'PRS',
            'source': 'data/prs_scores.csv',
            'main_column': 'PRS311',
            'other_columns': ['PRS313', 'PRS311_noimp', 'PRS313_noimp'],
        },
        {
            'score': '13-gene burden, case-control',
            'source': 'data/29_rowlands_extension/section7/carrier_burden_case_control.csv.gz',
            'main_column': 'B_raw_13gene_main',
            'other_columns': [
                'B_standardized_13gene_main',
                'carried_gene_count_13gene_main',
                'carrier_class_13gene_main',
            ],
        },
        {
            'score': '13-gene burden, diagnosis age',
            'source': 'data/29_rowlands_extension/section7/carrier_burden_diagnosis_age.csv.gz',
            'main_column': 'B_raw_13gene_main',
            'other_columns': [
                'B_standardized_13gene_main',
                'carried_gene_count_13gene_main',
                'carrier_class_13gene_main',
            ],
        },
        {
            'score': 'C comparison score',
            'source': 'data/29_rowlands_extension/section8/comparable_scores_case_control.csv.gz',
            'main_column': 'B_raw_C',
            'other_columns': ['B_standardized_C'],
        },
    ]
    for item in restricted:
        assert (ROOT / item['source']).is_file()
        item['internal_review_copy'] = {
            'PRS': 'manuscript/datas/PRS_scores.csv',
            '13-gene burden, case-control': 'manuscript/datas/burden_13gene_case_control.csv',
            '13-gene burden, diagnosis age': 'manuscript/datas/burden_13gene_diagnosis_age.csv',
        }.get(item['score'])
        item['copied_to_publication_directory'] = bool(item['internal_review_copy'])
        item['access'] = (
            'Internal workspace only; participant-level data; not for publication or download'
        )
    (OUT / 'restricted_score_sources.json').write_text(json.dumps(restricted, indent=2) + '\n')
    (OUT / 'manifest.json').write_text(
        json.dumps(
            {
                'main_version': 'C: 13-gene PTV + eligible P/LP missense',
                'final_variant_definition_frozen': False,
                'models_rerun': False,
                'manuscript_1002_modified': False,
                'numeric_policy': 'Preserve source numeric strings; disclosure coarsening remains in the explicitly named distribution file',
                'files': records,
            },
            indent=2,
        )
        + '\n'
    )
    print(
        f'Collected {len(records)} aggregate CSV files in Manuscript/datas; restricted score files indexed in place.'
    )


if __name__ == '__main__':
    main()

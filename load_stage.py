"""Load a numbered method module without executing a data pipeline."""
import importlib.util
from pathlib import Path

STAGES = {
 'cohort': '01_study_cohort.py', 'sample_qc': '02_sample_qc.py',
 'prs': '03_common_variants_prs.py', 'genes': '04_gene_selection.py',
 'annotation': '05_variant_annotation.py', 'genotyping': '06_variant_genotype_qc.py',
 'burden': '07_burden_score.py', 'case_control': '08_case_control_models.py',
 'age': '09_age_analysis.py', 'sensitivity': '10_sensitivity_analysis.py',
 'interaction': '11_interaction_analysis.py'}


def load_stage(name):
    spec = importlib.util.spec_from_file_location('public_method_' + name,
        Path(__file__).resolve().parent / STAGES[name])
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--list', action='store_true')
    args = parser.parse_args()
    if args.list:
        for name, filename in STAGES.items():
            print(name + ': ' + filename)

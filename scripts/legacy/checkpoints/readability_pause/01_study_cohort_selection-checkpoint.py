"""Study cohort selection and phenotype definition"""
from pathlib import Path
import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.functions.cli import run_local_step, task

select_study_cohort = task('cohort', 'select_study_cohort')
define_phenotypes_and_ehr_dates = task('phenotypes', 'define_phenotypes_and_ehr_dates')

def main():
    computations = {
        'select': select_study_cohort,
        'phenotypes': define_phenotypes_and_ehr_dates,
    }
    run_local_step('Study cohort selection and phenotype definition', computations, default_tasks=['select', 'phenotypes'])

if __name__ == "__main__":
    main()

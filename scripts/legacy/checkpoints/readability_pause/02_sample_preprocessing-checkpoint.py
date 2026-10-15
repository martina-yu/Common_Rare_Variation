"""Sample QC and clinical/technical covariates"""
from pathlib import Path
import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.functions.cli import run_local_step, task

preprocess_samples = task('sample_qc', 'preprocess_samples')
derive_clinical_covariates = task('clinical_covariates', 'derive_clinical_covariates')

def main():
    computations = {
        'sample-qc': preprocess_samples,
        'covariates': derive_clinical_covariates,
    }
    run_local_step('Sample QC and clinical/technical covariates', computations, default_tasks=['sample-qc'])

if __name__ == "__main__":
    main()

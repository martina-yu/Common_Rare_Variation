"""Common SNP harmonization, selection and PRS calculation"""
from pathlib import Path
import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.functions.cli import run_local_step, task

calculate_common_variant_prs = task('prs', 'calculate_common_variant_prs')
audit_prs_matches = task('prs', 'audit_prs_matches')

def main():
    computations = {
        'prs': calculate_common_variant_prs,
        'match-audit': audit_prs_matches,
    }
    run_local_step('Common SNP harmonization, selection and PRS calculation', computations, default_tasks=['prs', 'match-audit'])

if __name__ == "__main__":
    main()

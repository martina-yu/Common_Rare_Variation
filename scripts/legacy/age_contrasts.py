"""Historical command alias; the analysis remains visible in exploratory/age_ehr_checks.py."""

from pathlib import Path
import subprocess
import sys

if __name__ == '__main__':
    root = Path(__file__).resolve().parents[2]
    arguments = ['--task', 'legacy-group-contrasts'] + sys.argv[1:]
    subprocess.run(
        [sys.executable, str(root / 'scripts/exploratory/age_ehr_checks.py'), *arguments],
        check=True,
    )

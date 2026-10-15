"""Run the optional review frailty sensitivity with the existing specification."""

from pathlib import Path
import subprocess
import sys

if __name__ == '__main__':
    root = Path(__file__).resolve().parents[2]
    subprocess.run(
        [
            sys.executable,
            str(root / 'scripts/review/review_1008.py'),
            '--task',
            'frailty',
            *sys.argv[1:],
        ],
        check=True,
    )

"""Reproduce the existing PRS matching/gate attrition audit."""

from pathlib import Path
import subprocess
import sys

if __name__ == '__main__':
    root = Path(__file__).resolve().parents[2]
    subprocess.run(
        [
            sys.executable,
            str(root / 'scripts/03_common_variants_prs.py'),
            '--task',
            'match-audit',
            *sys.argv[1:],
        ],
        check=True,
    )

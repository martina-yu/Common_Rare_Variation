"""Historical entry for the explicitly post-hoc numbered interaction analysis."""

from pathlib import Path
import subprocess
import sys

if __name__ == '__main__':
    root = Path(__file__).resolve().parents[2]
    subprocess.run(
        [
            sys.executable,
            str(root / 'scripts/11_interaction_analysis.py'),
            '--task',
            'audit',
            *sys.argv[1:],
        ],
        check=True,
    )

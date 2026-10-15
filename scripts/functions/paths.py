"""Repository paths; scientific output locations remain unchanged."""

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def implementation_sources(root=PROJECT_ROOT):
    return sorted(
        p
        for p in (Path(root) / 'scripts').rglob('*.py')
        if not any(part in {'checkpoints', '.ipynb_checkpoints', '__pycache__'} for part in p.parts)
    )

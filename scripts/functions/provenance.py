# Documentation: scripts/README.md.
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / 'scripts/refactor_manifest.json'


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _verify_historical_layout(inputs, root=ROOT, manifest_path=MANIFEST):
    manifest = json.loads(Path(manifest_path).read_text())
    records = {r['original_path']: r for r in manifest['files']}
    checked_migration = False
    for relative, expected in inputs.items():
        current = Path(root) / relative
        if current.is_file() and sha256(current) == expected:
            continue
        record = records.get(relative)
        assert record and record['active_path'], f'Locked input changed or missing: {relative}'
        assert record['original_sha256'] == expected, f'Unrecognized historical hash: {relative}'
        assert (
            sha256(Path(root) / record['archive_path']) == expected
        ), f'Historical snapshot changed: {relative}'
        assert (
            sha256(Path(root) / record['active_path']) == record['active_sha256']
        ), f'Migrated source changed; create a new input lock: {relative}'
        if not checked_migration:
            for source, digest in manifest['implementation_sha256'].items():
                assert (
                    sha256(Path(root) / source) == digest
                ), f'Migrated implementation changed; create a new input lock: {source}'
            checked_migration = True
    return True


def verify_frozen_inputs(inputs, root=ROOT, manifest_path=MANIFEST):
    """Verify exact data hashes and explicitly recorded source migrations."""
    root = Path(root)
    migration_path = root / 'scripts/readability_manifest.json'
    if not migration_path.exists():
        return _verify_historical_layout(inputs, root, manifest_path)
    migration = json.loads(migration_path.read_text())
    previous = json.loads(
        (root / migration['reference_snapshot'] / 'scripts/refactor_manifest.json').read_text()
    )
    records = {record['original_path']: record for record in previous['files']}
    migrated = False
    for relative, expected in inputs.items():
        current = root / relative
        if current.is_file() and sha256(current) == expected:
            continue
        assert relative.startswith('scripts/'), f'Locked data changed or missing: {relative}'
        snapshots = [
            root / migration[key] / relative
            for key in ['reference_snapshot', 'interrupted_snapshot']
        ]
        if any(path.is_file() and sha256(path) == expected for path in snapshots):
            migrated = True
            continue
        record = records.get(relative)
        assert record and record.get('active_path'), f'Locked source changed or missing: {relative}'
        assert (
            record['original_sha256'] == expected
        ), f'Unrecognized historical source hash: {relative}'
        assert (
            sha256(root / record['archive_path']) == expected
        ), f'Historical source snapshot changed: {relative}'
        reference = root / migration['reference_snapshot']
        assert (
            sha256(reference / record['active_path']) == record['active_sha256']
        ), f'Pre-existing source drift: {relative}'
        for source, digest in previous['implementation_sha256'].items():
            assert (
                sha256(reference / source) == digest
            ), f'Pre-existing implementation drift: {source}'
        migrated = True
    if migrated:
        assert migration.get('implementation_sha256'), 'Readability implementation lock is missing'
        for source, digest in migration['implementation_sha256'].items():
            assert sha256(root / source) == digest, f'Readable implementation changed: {source}'
    return True

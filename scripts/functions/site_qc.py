# Documentation: scripts/README.md.
import json
import argparse
import hashlib
import importlib.util
import os
from pathlib import Path

import hail as hl
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'data/29_rowlands_extension/section6'


def source_site_pass(filters):
    """Empty defined Hail filter set is PASS; missing source information fails."""
    return filters is not None and (len(filters) == 0 or set(filters) == {'PASS'})


def canonical_crosswalk():
    raw_path = OUT / 'original_VDS_site_filters.json'
    source = json.loads(raw_path.read_text())
    records = []
    for key, value in source.items():
        chrom, pos, ref, alts = key.split(':', 3)
        records.append((chrom, pos, ref, alts, 'PASS' if value['PASS'] else 'FAIL'))
    mapping = qc.normalized_native_sites(
        records, ROOT / 'data/refdata/Homo_sapiens_assembly38.fasta'
    )
    desired = set()
    for name in ['ptv_for_call_qc.tsv.gz', 'eligible_missense.tsv.gz']:
        desired.update(pd.read_csv(OUT / name, sep='\t', usecols=['normalized_key']).normalized_key)
    canonical = {}
    for native, key in mapping.items():
        if key in desired:
            entry = canonical.setdefault(key, {'source_native_keys': [], 'PASS': True})
            entry['source_native_keys'].append(':'.join(native[:4]))
            entry['PASS'] &= native[5] == 'PASS'
    (OUT / 'original_VDS_normalized_allele_filters.json').write_text(json.dumps(canonical) + '\n')
    manifest_path = OUT / 'original_VDS_site_filter_manifest.json'
    manifest = json.loads(manifest_path.read_text())
    audit = (
        pd.read_csv(OUT / 'native_site_ALT_QC_audit.tsv.gz', sep='\t', dtype=str)
        if (OUT / 'native_site_ALT_QC_audit.tsv.gz').exists()
        else None
    )
    if audit is not None:
        seen = set(audit.normalized_key)
        manifest['every_found_exome_allele_matches_original_VDS_after_full_normalization'] = all(
            k in canonical for k in seen
        )
        manifest['every_found_exome_allele_original_VDS_PASS'] = all(
            canonical.get(k, {}).get('PASS', False) for k in seen
        )
    manifest['original_VDS_rows_sha256'] = hashlib.sha256(raw_path.read_bytes()).hexdigest()
    manifest['canonical_crosswalk'] = (
        'complete chr:pos:ref:alt normalized with fixed GRCh38; all equivalent source sites must PASS'
    )
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    print('Original VDS filters joined through complete normalized allele keys.', flush=True)

# Documentation: scripts/README.md.
import csv
import gzip
import hashlib
import importlib.util
import json
import os
import re
from pathlib import Path

import hail as hl
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / 'data/29_rowlands_extension'
S6 = BASE / 'section6'
OUT = BASE / 'section7'
FIELDS = [
    'sample_id',
    'gene',
    'chr',
    'pos',
    'ref',
    'alt',
    'transcript',
    'consequence',
    'PTV/missense',
    'LOFTEE',
    'ClinVar_release',
    'classification',
    'stars',
    'condition',
    'low_penetrance_flag',
    'conflict_flag',
    'EUR_AF',
    'site_filter',
    'genotype_filter',
    'GQ',
    'DP',
    'VAF',
    'keep_flag',
    'rejection_reason',
    'normalized_key',
    'source_native_key',
    'source_ALT_index',
    'allele_dosage',
    'AD_total',
    'sensitivity_keep_flag',
    'source',
    'supporting_transcript_rows',
    'RCV_SCV_evidence',
]

from scripts.functions.variants import terms


def load_module(name, filename):
    s = importlib.util.spec_from_file_location(name, ROOT / 'scripts/functions' / filename)
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


def positive_float(value):
    try:
        x = float(value)
        return x if __import__('math').isfinite(x) else None
    except (ValueError, TypeError):
        return None

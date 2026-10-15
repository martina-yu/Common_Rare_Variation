# Documentation: scripts/README.md.
import argparse
from collections import Counter, defaultdict
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'data/29_rowlands_extension/section6'
MOUNT = Path('/tmp/codex_expanded_missense_aou_ro')
SETS = ['raw', 'site_PASS', 'FT_PASS', 'GQ20', 'main', 'CH_sensitive']


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()


def call_qc(gt, ad, gq, target_index, n_alts, gene):
    try:
        alleles = [int(x) for x in gt.replace('|', '/').split('/')]
        if len(alleles) != 2 or any(x < 0 or x > n_alts for x in alleles):
            return None
        dose = alleles.count(target_index)
        if not dose:
            return None
        gq = float(gq)
        depths = [int(x) for x in ad.split(',')]
        if len(depths) != n_alts + 1 or any(x < 0 for x in depths) or sum(depths) <= 0:
            return {'gq_pass': math.isfinite(gq) and gq >= 20, 'valid_AD': False}
        vaf = depths[target_index] / sum(depths)
        main = vaf >= 0.8 if dose == 2 else 0.3 <= vaf <= 0.8
        sensitive = (
            vaf >= 0.8 if dose == 2 else 0.3 <= vaf <= (0.7 if gene in {'TP53', 'NF1'} else 0.8)
        )
        return dict(
            gq_pass=math.isfinite(gq) and gq >= 20,
            valid_AD=True,
            AD_total=sum(depths),
            vaf=vaf,
            dose=dose,
            main=main,
            sensitive=sensitive,
        )
    except (ValueError, IndexError):
        # Still an ALT carrier for attrition, but never a QC-passed carrier.
        try:
            alleles = [int(x) for x in gt.replace('|', '/').split('/')]
            if len(alleles) == 2 and target_index in alleles:
                return {'gq_pass': False, 'valid_AD': False}
        except ValueError:
            pass
        return None


def normalized_native_sites(raw, fasta):
    """Normalize single ALT site copies; leave original GT/AD indices untouched."""
    tokens, lines = {}, ['##fileformat=VCFv4.2']
    for l in Path(str(fasta) + '.fai').read_text().splitlines():
        chrom, length, *_ = l.split('\t')
        lines.append(f'##contig=<ID={chrom},length={length}>')
    lines.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO')
    for native in raw:
        chrom, pos, ref, alts, filt = native
        for index, alt in enumerate(alts.split(','), 1):
            if any(x not in 'ACGTN' for x in alt):
                continue
            token = f'n{len(tokens)}'
            tokens[token] = (chrom, pos, ref, alts, index, filt)
            lines.append(f'{chrom}\t{pos}\t{token}\t{ref}\t{alt}\t.\t.\t.')
    if not tokens:
        return {}
    norm = subprocess.run(
        ['bcftools', 'norm', '-f', str(fasta), '-c', 'e', '-m', '-any', '-Ov', '-'],
        input='\n'.join(lines) + '\n',
        capture_output=True,
        text=True,
        check=True,
    )
    mapping = {}
    for l in norm.stdout.splitlines():
        if l.startswith('#'):
            continue
        chrom, pos, token, ref, alt, *_ = l.split('\t')
        mapping[tokens[token]] = f'{chrom}:{pos}:{ref}:{alt}'
    assert len(mapping) == len(tokens)
    return mapping


def safe_count(n):
    return int(n) if n == 0 or n >= 20 else '<20'

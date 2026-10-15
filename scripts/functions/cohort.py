"""Shared ancestry-file extraction, with the original release assertions."""

from scripts.functions.paths import PROJECT_ROOT
import collections
import io
import os
import re
import subprocess
import sys
import time
import numpy as np
import pandas as pd
from google.cloud import bigquery
from openpyxl import load_workbook

PROJECT_ID = os.environ.get('GOOGLE_PROJECT', '')
MIN_N = 21
ANC_GCS = 'gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/aux/ancestry/'
ANC_PREDS = ANC_GCS + 'ancestry_preds.tsv'
ANCESTRY_CLASSES = ['afr', 'amr', 'eas', 'eur', 'mid', 'sas']
N_ANC_PC = 16
OTH_THRESHOLD = 0.75


def disclose_group_counts(counter, label, min_n=MIN_N):
    """Print a category breakdown under the n>=20 rule.

    Categories whose true n is below the floor are merged into one bucket; if
    the bucket is itself below the floor the smallest remaining category joins
    it, repeatedly, until it clears. A bucket that clears the floor may be
    printed exactly. Same rule v1 used (disclosed_ancestry_rows): it exists so
    that total - sum(printed rows) cannot reconstruct a suppressed cell.
    """
    items = {k: v for k, v in counter.items() if v >= min_n}
    bucket = {k: v for k, v in counter.items() if v < min_n}
    while bucket and sum(bucket.values()) < min_n and items:
        k = min(items, key=items.get)
        bucket[k] = items.pop(k)
    log(f'  {label}')
    for k, v in sorted(items.items(), key=lambda kv: -kv[1]):
        log(f'      {k:<24s} {v:>7,d}')
    if bucket:
        name = 'other_combined(' + '+'.join(sorted(bucket)) + ')'
        log(f'      {name:<24s} {sum(bucket.values()):>7,d}   <- merged to clear n>={min_n}')
    log(f"      {'TOTAL':<24s} {sum(counter.values()):>7,d}")


def extract_ancestry_preds(keep_ids):
    """Stream the release-wide ancestry table, keep our cohort, verify meaning.

    Returns (cohort DataFrame, release-wide label Counter).

    The two assertions are the point of this function. Both facts it relies on
    are undocumented in the AoU paper and were measured from the release, so a
    future release that changes either must fail here rather than quietly
    produce a differently-defined cohort.
    """
    import ast

    it = gs_cat_lines(ANC_PREDS)
    hdr = next(it).rstrip('\n').split('\t')
    expected = [
        'research_id',
        'ancestry_pred',
        'probabilities',
        'pca_features',
        'ancestry_pred_other',
    ]
    log(f'  header: {hdr}')
    if hdr != expected:
        raise SystemExit(f'FATAL: unexpected header; expected {expected}')
    n_all = 0
    release_lab = collections.Counter()
    argmax_check = collections.defaultdict(collections.Counter)
    rows = []
    min_called, max_oth = (1.0, 0.0)
    for line in it:
        f = line.rstrip('\n').split('\t')
        if len(f) < 5:
            continue
        rid, lab, probs, feats, oth = (f[0], f[1], f[2], f[3], f[4])
        n_all += 1
        release_lab[lab] += 1
        pr = ast.literal_eval(probs)
        if len(pr) != len(ANCESTRY_CLASSES):
            raise SystemExit(f'FATAL: {len(pr)} probabilities, expected {len(ANCESTRY_CLASSES)}')
        argmax_check[lab][ANCESTRY_CLASSES[max(range(len(pr)), key=lambda i: pr[i])]] += 1
        m = max(pr)
        if oth == 'oth':
            max_oth = max(max_oth, m)
        else:
            min_called = min(min_called, m)
        if rid in keep_ids:
            fe = ast.literal_eval(feats)
            if len(fe) != N_ANC_PC:
                raise SystemExit(f'FATAL: {len(fe)} PCs, expected {N_ANC_PC}')
            rows.append([rid, lab, oth] + [round(x, 6) for x in pr] + fe)
    log(f'  rows streamed (whole release): {n_all:,d}')
    log(f'  cohort rows matched:           {len(rows):,d} of {len(keep_ids):,d}')
    if len(rows) != len(keep_ids):
        raise SystemExit(
            'FATAL: some cohort samples have no ancestry call; that has to be understood before anything downstream'
        )
    log()
    log(
        f"  ASSERT 1 -- probabilities are ordered ({', '.join(ANCESTRY_CLASSES)}), i.e. ancestry_pred is their argmax:"
    )
    for k in sorted(argmax_check):
        got = list(argmax_check[k])
        log(f'      ancestry_pred={k:<4s} argmax is always {got[0]:<4s} ({release_lab[k]:,d} rows)')
        if len(got) != 1 or got[0] != k:
            raise SystemExit(
                f'FATAL: ancestry_pred={k} is not the argmax under this class order (saw {dict(argmax_check[k])})'
            )
    log('      -> PASS: ancestry_pred == argmax(probabilities), no threshold')
    log()
    log(f"  ASSERT 2 -- 'oth' still splits at {OTH_THRESHOLD} in this release:")
    log(f'      min max-prob among CALLED samples: {min_called}')
    log(f"      max max-prob among 'oth' samples:  {max_oth}")
    if not max_oth <= OTH_THRESHOLD < min_called:
        raise SystemExit(
            f'FATAL: the release no longer splits at {OTH_THRESHOLD}; re-derive it before using either label column to define a cohort'
        )
    log(f"      -> PASS: ancestry_pred_other = argmax if max > {OTH_THRESHOLD}, else 'oth'")
    cols = (
        ['person_id', 'ancestry_pred', 'ancestry_pred_other']
        + [f'p_{c}' for c in ANCESTRY_CLASSES]
        + [f'PC{i + 1}' for i in range(N_ANC_PC)]
    )
    return (pd.DataFrame(rows, columns=cols), release_lab)


def gs_cat_lines(url):
    """Stream a GCS object line by line: read-only, requester-pays, never
    written to disk. Used for the 220 MB ancestry table, of which we keep only
    the cohort's rows."""
    p = subprocess.Popen(
        ['gsutil', '-u', PROJECT_ID, 'cat', url],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1 << 20,
    )
    for line in p.stdout:
        yield line
    p.stdout.close()
    if p.wait() != 0:
        raise SystemExit(f'FATAL: gsutil cat {url} failed:\n{p.stderr.read()}')


def log(m=''):
    """Progress goes to stderr, because stdout carries the workbook bytes."""
    print(m, file=sys.stderr, flush=True)


def parse_hail_scores(series):
    """training_pca.tsv keeps its PC vector as a bracketed string of Java
    floats ('1.0E-4'); float() reads those, np.fromstring does not."""
    return np.array([[float(x) for x in s.strip('[]').split(',')] for s in series])

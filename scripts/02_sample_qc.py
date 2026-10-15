"""Apply sample QC and derive clinical and technical covariates.

Inputs: female WGS cohort, phenotype dates and AoU QC resources.
The existing chr21 missingness proxy and available relatedness threshold are preserved.
"""

from __future__ import annotations

from pathlib import Path
import sys

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.functions.paths import PROJECT_ROOT
from scripts.functions.cli import parse_tasks


# Prespecified scientific constants (values unchanged).
KINSHIP_THRESHOLD_AVAILABLE = 0.1

TASKS = ('sample-qc', 'covariates')
if __name__ == '__main__':
    selected_tasks = parse_tasks(__doc__, TASKS, ['sample-qc'])

import math

import os

import subprocess


import networkx as nx

import numpy as np

import pandas as pd

import json

import matplotlib

import matplotlib.pyplot as plt

import seaborn as sns

from functools import partial
from scripts.functions.disclosure import format_count

disclose = partial(format_count, style='sample_qc')


def run(cmd, **kwargs):
    print(f'$ {cmd}')
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, **kwargs)
    if result.stdout:
        print(result.stdout)
    if result.stderr:
        print(result.stderr, file=sys.stderr)
    if result.returncode != 0:
        print(f'!!! command exited with code {result.returncode}', file=sys.stderr)
    return result


def sample_qc_section(title):
    print('\n' + '=' * 90)
    print(title)
    print('=' * 90)


def kinship_pruning(cohort):
    sample_qc_section('2. Kinship pruning (maximum independent set)')
    print(
        f'Plan requirement threshold: KING kinship >{KINSHIP_THRESHOLD_REQUESTED}(Within second degree of intimacy)'
    )
    local_pairs = f'{SCRATCH}/samples_relatedness.tsv'
    if not os.path.exists(local_pairs):
        run(f'gsutil -u {GOOGLE_PROJECT} cp {AUX_REL}/samples_relatedness.tsv {local_pairs}')
    pairs = pd.read_csv(local_pairs, sep='\t')
    pairs = pairs.rename(columns={'i.s': 'id1', 'j.s': 'id2', 'kin': 'kinship'})
    pairs['id1'] = pairs['id1'].astype('int64')
    pairs['id2'] = pairs['id2'].astype('int64')
    print(
        f"Official samples_relatedness.tsv Total number of pairs:{len(pairs)}, Measured minimum kinship value:{pairs['kinship'].min():.4f}"
    )
    print(
        f"⚠️ Data availability issue: AoU official pairwise relatedness file only exports pairs with kinship > 0.1(Actually measured minimum value ={pairs['kinship'].min():.4f}, not 0.0884), so it cannot be between 0.0884-0.1Identify weak second-degree relationship pairs in a narrow interval. 0.1 falls exactly within the KING second degree interval (0.0884-0.177).is a usable approximation that is more conservative (stricter) than the plan requirement (0.0884) - using it will exclude the theoretical thresholdSlightly less (missing weak second-degree pairs with kin between 0.0884-0.1), but does not introduce false positive exclusion.This analysis uses the actual thresholds available ={KINSHIP_THRESHOLD_AVAILABLE}, and indicate this caveat in the report."
    )
    cohort_ids = set(cohort['person_id'])
    in_cohort = pairs['id1'].isin(cohort_ids) & pairs['id2'].isin(cohort_ids)
    cohort_pairs = pairs.loc[in_cohort & (pairs['kinship'] > KINSHIP_THRESHOLD_AVAILABLE)].copy()
    print(
        f'\n candidate group and kinship>{KINSHIP_THRESHOLD_AVAILABLE}Number of pairs:{disclose(len(cohort_pairs))}'
    )
    G = nx.Graph()
    G.add_edges_from(zip(cohort_pairs['id1'], cohort_pairs['id2']))
    n_edges = G.number_of_edges()
    n_involved = G.number_of_nodes()
    print(
        f'Relationship diagram:{disclose(n_involved)}The individual involves at least one edge, a total of{disclose(n_edges)}strip edge'
    )
    Gc = G.copy()
    excluded = set()
    while Gc.number_of_edges() > 0:
        degrees = dict(Gc.degree())
        node = max(degrees, key=lambda n: degrees[n])
        if degrees[node] == 0:
            break
        excluded.add(node)
        Gc.remove_node(node)
    print(f'Greedy maximum independent set algorithm eliminates people:{disclose(len(excluded))}')
    remaining_edges_check = sum(
        (
            1
            for a, b in zip(cohort_pairs['id1'], cohort_pairs['id2'])
            if a not in excluded and b not in excluded
        )
    )
    print(
        f'The number of remaining edges in the queue after elimination (should be 0):{remaining_edges_check}'
    )
    local_flagged = f'{SCRATCH}/relatedness_flagged_samples.tsv'
    if not os.path.exists(local_flagged):
        run(
            f'gsutil -u {GOOGLE_PROJECT} cp {AUX_REL}/relatedness_flagged_samples.tsv {local_flagged}'
        )
    official_flagged = pd.read_csv(local_flagged, sep='\t')['sample_id'].astype('int64')
    official_flagged_set = set(official_flagged) & cohort_ids
    overlap = len(excluded & official_flagged_set)
    print(
        f'\n [Cross-validation, for reference only] AoU official relatedness_flagged_samples.tsv (the whole queue is about 245k-535k people and MIS based on the picture)Number of people in this cohort:{disclose(len(official_flagged_set))}；The number of people who overlap with this script (only recalculating MIS on this cohort subgraph) removes the list:{disclose(overlap)}。The two do not have to be exactly the same, because the official MIS is calculated on the whole AoU WGS population map (may be excluded due to relatives outside the cohort)Those who could have been retained in the subgraph, and vice versa), this project will be recalculated independently within its own candidate group subgraph according to the plan requirements,It is more consistent with the actual sample set used in this analysis.'
    )
    out = cohort.copy()
    out['flag_related_exclude'] = out['person_id'].isin(excluded)
    return (out[['person_id', 'flag_related_exclude']], n_edges, len(excluded))


def build_representative_acaf_slice(cohort):
    sample_qc_section(
        '3/4 Preparation: Missing rate + ACAF threshold representative variant subsampling shared by HWE'
    )
    print(
        'Target chromosome: chr21 (the smallest among the autosomes, still completely covering the real genotype data structure on this chromosome)'
    )
    bed_path = f'{ACAF_BED_DIR}/acaf_threshold.chr21.bed'
    bim_path = f'{ACAF_BED_DIR}/acaf_threshold.chr21.bim'
    fam_path = f'{ACAF_BED_DIR}/acaf_threshold.chr21.fam'
    du = run(f'gsutil -u {GOOGLE_PROJECT} du {bed_path}')
    total_bytes = int(du.stdout.strip().split()[0])
    print(f'chr21.bed total bytes:{total_bytes:,} ({total_bytes / 1024 ** 3:.1f} GiB)')
    local_fam = f'{SCRATCH}/chr21.fam'
    if not os.path.exists(local_fam):
        run(f'gsutil -u {GOOGLE_PROJECT} cp {fam_path} {local_fam}')
    fam = pd.read_csv(
        local_fam,
        sep='\t',
        header=None,
        na_filter=False,
        dtype=str,
        names=['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENO'],
    )
    n_samples = len(fam)
    block_size = math.ceil(n_samples / 4)
    header_bytes = 3
    n_variants_total = (total_bytes - header_bytes) // block_size
    print(f'chr21 .fam sample number (full AoU queue, non-cohort subset):{n_samples:,}')
    print(f'SNP-major .bed Bytes per mutated block:{block_size:,}')
    print(f'The total variation of chr21 deduced from the file size:{n_variants_total:,}')
    print(
        '⚠️ Scale assessment: The total size of the full ACAF threshold callset (22 autosomal chromosomes, plink_bed format) is about 16TiB.chr21 alone has 247GiB/535,662 individuals. This machine (overlay root partition has about 49GB available, 25GB memory, 4 cores,No local Hail/Spark cluster, no GCS FUSE mount) It is physically impossible to download the entire table or use plink2 for direct streaming.Not to mention doing a full sample-level deletion rate or HWE scan of the 22 chromosomes. If you want to do full-scale, genome-widesample_qc(), the correct path is to set up a Hail-enabled Dataproc cluster (GenomicAnalysis environment) run hl.sample_qc()/hl.variant_qc() on hail.vds or ACAF Hail MatrixTable,The order of magnitude of a job of the same scale (full queue, full genome sample_qc) given in the official document is usually tens of minutes to hours.(Depending on the number of workers), this exceeds the reasonable overhead of this interactive session, so it is not executed and only this evaluation is recorded.'
    )
    print(
        'Alternative plan (actual implementation this time): Make 3 windows based on mutation index inside chr21.bed (500 mutations each at the beginning, middle and end of the file,Accurate byte-range slices of a total of 1500 mutations), recover a small but real .bed/.bim/.fam, useplink2 runs real --missing/--hardy. This is not a genome-wide result, but a limited representative subsampling,This limitation is clearly noted in the report below.'
    )
    windows = [('start', 5000), ('mid', n_variants_total // 2), ('end', n_variants_total - 5500)]
    n_per_window = 500
    local_bim = f'{SCRATCH}/chr21.bim'
    if not os.path.exists(local_bim):
        run(f'gsutil -u {GOOGLE_PROJECT} cp {bim_path} {local_bim}')
    bed_chunks = []
    bim_chunks = []
    for label, start_idx in windows:
        byte_start = header_bytes + start_idx * block_size
        byte_end = byte_start + n_per_window * block_size - 1
        chunk_path = f'{SCRATCH}/chr21_{label}.bedchunk'
        run(f'gsutil -u {GOOGLE_PROJECT} cat -r {byte_start}-{byte_end} {bed_path} > {chunk_path}')
        actual_size = os.path.getsize(chunk_path)
        expected_size = n_per_window * block_size
        print(
            f"window[{label}] variant_idx=[{start_idx},{start_idx + n_per_window}) byte_range=[{byte_start},{byte_end}] Download bytes ={actual_size:,} (expect={expected_size:,}, {('OK' if actual_size == expected_size else 'MISMATCH!')})"
        )
        with open(chunk_path, 'rb') as f:
            bed_chunks.append(f.read())
        bim_slice_path = f'{SCRATCH}/chr21_{label}.bim'
        run(f"sed -n '{start_idx + 1},{start_idx + n_per_window}p' {local_bim} > {bim_slice_path}")
        bim_chunks.append(
            pd.read_csv(
                bim_slice_path, sep='\t', header=None, names=['CHR', 'ID', 'CM', 'BP', 'A1', 'A2']
            )
        )
    mini_bed_path = f'{SCRATCH}/mini.bed'
    with open(mini_bed_path, 'wb') as f:
        f.write(bytes([108, 27, 1]))
        for chunk in bed_chunks:
            f.write(chunk)
    mini_bim = pd.concat(bim_chunks, ignore_index=True)
    mini_bim.to_csv(f'{SCRATCH}/mini.bim', sep='\t', header=False, index=False)
    fam.to_csv(f'{SCRATCH}/mini.fam', sep='\t', header=False, index=False)
    n_total_slice_variants = len(mini_bim)
    print(
        f'\n synthesized mini .bed/.bim/.fam:{n_total_slice_variants}variations x{n_samples:,}individual(mini.bed size={os.path.getsize(mini_bed_path):,} bytes)'
    )
    return f'{SCRATCH}/mini'


def sample_missingness(mini_prefix, cohort):
    sample_qc_section(
        '3. Genotype deletion rate filtering (sample-level, see below for threshold determination)'
    )
    print(
        f'Threshold selection: sample-level missing rate >{SAMPLE_MISSING_THRESHOLD:.0%}then eliminate'
    )
    print(
        'Reason: The plan requires a range of 2%-5%; the common variant used in this project is PRS313/311 (313/311 sites)'
    )
    print(
        '+ HWE filtering (itself the second level of quality control), so the strictest 2% of the interval is not taken (it will cause factor sampling to estimate noise'
    )
    print(
        '(Boundary sample caused by accidental injury), nor the loosest 5%; 3% is the customary value for standard GWAS array QC (Anderson et al. 2010'
    )
    print(
        'Recommend --mind 0.02~0.05, the middle value 0.03 is the most commonly cited) as the implementation choice for this analysis.'
    )
    keep_path = f'{SCRATCH}/cohort_keep.txt'
    cohort[['person_id']].assign(FID=0)[['FID', 'person_id']].to_csv(
        keep_path, sep='\t', header=False, index=False
    )
    vmiss_prefix = f'{SCRATCH}/premiss_variant_check'
    run(
        f'{PLINK2} --bfile {mini_prefix} --keep {keep_path} --missing variant-only --out {vmiss_prefix}'
    )
    vmiss = pd.read_csv(f'{vmiss_prefix}.vmiss', sep='\t')
    artifact_variant_missing_threshold = 0.2
    bad_variants = vmiss.loc[vmiss['F_MISS'] > artifact_variant_missing_threshold, 'ID']
    print(
        f'\n subsampling{len(vmiss)}Among mutations, variant-level deletion rate>{artifact_variant_missing_threshold:.0%}(The number of variants determined to be multi-allelic split represents artifact, not true genotype missing):{len(bad_variants)}'
    )
    if len(bad_variants) > 0:
        print('Examples of excluded variants (top 10):')
        print(
            vmiss.loc[vmiss['ID'].isin(bad_variants)]
            .sort_values('F_MISS', ascending=False)
            .head(10)[['ID', 'F_MISS']]
            .to_string()
        )
    exclude_path = f'{SCRATCH}/artifact_variants_exclude.txt'
    bad_variants.to_csv(exclude_path, header=False, index=False)
    out_prefix = f'{SCRATCH}/missing_out'
    exclude_arg = f'--exclude {exclude_path}' if len(bad_variants) > 0 else ''
    run(
        f'{PLINK2} --bfile {mini_prefix} --keep {keep_path} {exclude_arg} --missing sample-only --out {out_prefix}'
    )
    print(
        f'The final number of markers used for sample-level missing rate estimation:{len(vmiss) - len(bad_variants)} (Original 1500 -{len(bad_variants)}artifact)'
    )
    smiss = pd.read_csv(f'{out_prefix}.smiss', sep='\t')
    smiss = smiss.rename(columns={'IID': 'person_id', 'F_MISS': 'f_miss'})
    smiss['person_id'] = smiss['person_id'].astype('int64')
    print(
        f'\nplink2 --missing output sample number:{len(smiss)}(Should be equal to queue size{len(cohort)})'
    )
    print(
        'f_miss distribution (based on 1500 representative mutant samples, not the whole genome):'
    )
    print(smiss['f_miss'].describe().to_string())
    merged = cohort.merge(smiss[['person_id', 'f_miss']], on='person_id', how='left')
    n_no_result = merged['f_miss'].isna().sum()
    print(f'The number of people who failed to match the output of plink2:{disclose(n_no_result)}')
    merged['flag_high_missing'] = merged['f_miss'].isna() | (
        merged['f_miss'] > SAMPLE_MISSING_THRESHOLD
    )
    n_flag = int(merged['flag_high_missing'].sum())
    print(f'Number of people marked for high missing rate culling:{disclose(n_flag)}')
    return (merged[['person_id', 'f_miss', 'flag_high_missing']], exclude_path, len(bad_variants))


def hwe_control_only(mini_prefix, cohort, exclude_path=None):
    sample_qc_section(
        '4. HWE filtering (control group only, see below for P threshold determination)'
    )
    print(f'P value threshold selection:{HWE_P_THRESHOLD:.0e}')
    print(
        'Reason: The plan gives two levels: 1e-6 or 1e-8. The downstream uses of the common variant of this project are PRS313/311 and 313/311'
    )
    print(
        'For established locations (non-genome-wide newly discovered scans) and array/targeted-QC scenarios, the looser of the two levels is used.'
    )
    print(
        'The more conventional single-queue QC default value is 1e-6 (the customary threshold in the Anderson et al. 2010 GWAS QC tutorial);'
    )
    print(
        'The more stringent 1e-8 is reserved for the genome-wide exploratory scanning scenario of Step 3 rare variant burden (this scenario'
    )
    print(
        'It has a higher multiple testing burden and requires a stricter threshold to control false positives).'
    )
    print(
        "Use the multi-allelic split representation artifact exclusion list identified in step 3, and keep the same set of 'clean' markers used in both steps."
    )
    controls = cohort.loc[cohort['is_case'] == 0, ['person_id']].copy()
    keep_path = f'{SCRATCH}/control_keep.txt'
    controls.assign(FID=0)[['FID', 'person_id']].to_csv(
        keep_path, sep='\t', header=False, index=False
    )
    print(f'Number of people in the control group (for HWE test):{disclose(len(controls))}')
    out_prefix = f'{SCRATCH}/hwe_out'
    exclude_arg = f'--exclude {exclude_path}' if exclude_path else ''
    run(
        f'{PLINK2} --bfile {mini_prefix} --keep {keep_path} {exclude_arg} --hardy --out {out_prefix}'
    )
    hardy = pd.read_csv(f'{out_prefix}.hardy', sep='\t')
    print(
        f'\nplink2 --hardy output variation number:{len(hardy)}(Taken from subsampling, excluding artifact variation)'
    )
    print(hardy[['O(HET_A1)', 'E(HET_A1)', 'P']].describe().to_string())
    n_fail = int((hardy['P'] < HWE_P_THRESHOLD).sum())
    frac_fail = n_fail / len(hardy) if len(hardy) else float('nan')
    print(
        f'\n is in control group, subsampling{len(hardy)}Among the variations, HWE P<{HWE_P_THRESHOLD:.0e}Variation number:{n_fail} / {len(hardy)} ({frac_fail:.4%})'
    )
    print(
        f'⚠️ This number is based solely on{len(hardy)}A subsample of variants, not a true count of the whole genome ACAF threshold callset;'
    )
    print(
        'The reason is the same as step 3: the callset of all 22 chromosomes (16TiB) cannot be processed locally, and the Hail cluster is required to complete genome-wide'
    )
    print(
        "HWE scan. To extrapolate this ratio to the 'order of magnitude of expected deleted variants across the genome', you need to know the ACAF threshold"
    )
    print(
        'The total number of callset variants is not extrapolated here. Only the measured counts within the subsample are reported to avoid giving unverified large numbers.'
    )
    return (hardy, n_fail, frac_fail)


def load_candidates():
    sample_qc_section(
        '0. Candidate group: No. 01 female+WGS (including ancestry annotation) × No. 04 phenotype grouping'
    )
    client = sample_qc_bigquery.Client()
    anc = client.query(
        f'\n        SELECT person_id, is_eur\n        FROM `{GOOGLE_PROJECT}.my_analysis_results.female_wgs_ancestry_cohort`\n    '
    ).to_dataframe()
    pheno = pd.read_csv(
        f'{SAMPLE_QC_DATA_DIR}/04_phenotype_diagnostics_person_level.csv.gz',
        usecols=['person_id', 'group', 'age_at_diagnosis', 'bc_is_first_record'],
    )
    df = anc.merge(pheno, on='person_id', how='left', validate='one_to_one')
    if df['group'].isna().any():
        raise SystemExit('FATAL: The crowds of 01 and 04 are inconsistent.')
    valid_age = df['age_at_diagnosis'].notna() & (df['age_at_diagnosis'] >= 18)
    df['analysis_group'] = np.select(
        [(df['group'] == 'case') & valid_age, df['group'] == 'control'],
        ['case', 'control'],
        default='none',
    )
    print(f'female+WGS candidate pool:{disclose(len(df))}')
    for g in ['case', 'control']:
        print(f"  {g}: {disclose((df['analysis_group'] == g).sum())}")
    return df


def count_row(df, mask):
    sub = df.loc[mask & df['analysis_group'].isin(SAMPLE_QC_COLS), 'analysis_group'].value_counts()
    return {c: int(sub.get(c, 0)) for c in SAMPLE_QC_COLS}


def build_table_s2(df):
    """Build table s2."""
    sample_qc_section('Supplementary Table 2: Stepwise sample QC counts')
    steps = [
        (
            'Female (sex at birth) participants with WGS, case or control',
            df['analysis_group'].isin(SAMPLE_QC_COLS),
            True,
        ),
        ('Genetically inferred European ancestry (ancestry_pred_other = eur)', df['is_eur'], True),
        (
            'Unrelated (KING kinship <= 0.1, maximal unrelated set)',
            ~df['flag_related_exclude'],
            True,
        ),
        (
            f'Sample missingness <= {SAMPLE_MISSING_THRESHOLD:.0%} (Aim 1 analysis cohort)',
            ~df['flag_high_missing'],
            True,
        ),
        (
            'Diagnosis date not equal to first EHR record (Aim 4 analysis cohort)',
            df['bc_is_first_record'] != True,
            False,
        ),
    ]
    keep = pd.Series(True, index=df.index)
    raw = []
    for label, m, ctrl in steps:
        keep = keep & m
        row = count_row(df, keep)
        if not ctrl:
            row['control'] = None
        a = df.loc[keep & (df['analysis_group'] == 'case'), 'age_at_diagnosis']
        row['_age'] = (
            f'{a.median():.1f} [{a.quantile(0.25):.1f}-{a.quantile(0.75):.1f}]',
            f'{a.mean():.1f} ({a.std():.1f})',
        )
        raw.append((label, row))
    raw_df = pd.DataFrame(
        [{'step': l, **{k: v for k, v in r.items() if k != '_age'}} for l, r in raw]
    )
    raw_df.to_csv(f'{SAMPLE_QC_INTERNAL_DIR}/05_table_s2_raw_counts.csv', index=False)
    print(
        f'Raw step count write (without push):{SAMPLE_QC_INTERNAL_DIR}/05_table_s2_raw_counts.csv'
    )
    merged = [raw[0]]
    for label, row in raw[1:]:
        if len(merged) > 1:
            prev = merged[-1][1]
            drops = [
                prev[c] - row[c]
                for c in SAMPLE_QC_COLS
                if row[c] is not None and prev[c] is not None
            ]
            if any((0 < d < 20 for d in drops)):
                merged[-1] = (merged[-1][0] + '; ' + label, row)
                continue
        merged.append((label, row))
    n_merge = len(raw) - len(merged)
    print(f'Number of steps to merge by rule:{n_merge}')
    cells, var_of, published, aggregates = ({}, {}, {}, [])
    for i, (label, row) in enumerate(merged):
        lv = f'R{i}'
        cells[lv] = {c: row[c] if row[c] is not None else 0 for c in SAMPLE_QC_COLS}
        cells[lv]['total'] = sum((v for v in cells[lv].values()))
        published[lv] = {c for c in SAMPLE_QC_COLS if row[c] is not None}
        if all((row[c] is not None for c in SAMPLE_QC_COLS)):
            published[lv] |= {'total'}
        var_of[lv] = lv
        if i > 0:
            ex = f'X{i}'
            prev = cells[f'R{i - 1}']
            cells[ex] = {c: prev[c] - cells[lv][c] for c in list(SAMPLE_QC_COLS) + ['total']}
            published[ex] = set(published[lv]) & set(published[f'R{i - 1}'])
            var_of[ex] = ex
            aggregates.append((f'R{i - 1}', [lv, ex]))
    masked = sample_qc_decide(
        cells,
        var_of,
        published,
        column_partitions=[('total', 'case', 'control')],
        aggregates=aggregates,
    )
    n_masked = sum((len(v) for v in masked.values()))
    print(f'decide() determines the number of covered cells after taking the bottom:{n_masked}')
    out = []
    for i, (label, row) in enumerate(merged):
        lv = f'R{i}'
        rec = {'step': label}
        rec['case_age_at_diagnosis_median_iqr'], rec['case_age_at_diagnosis_mean_sd'] = row['_age']
        for c in list(SAMPLE_QC_COLS) + ['total']:
            if c not in published[lv]:
                rec[c] = 'n/a'
            else:
                rec[c] = sample_qc_render(cells[lv][c], c in masked[lv], None)
        if i > 0:
            ex = f'X{i}'
            for c in list(SAMPLE_QC_COLS) + ['total']:
                rec[f'excluded_{c}'] = (
                    'n/a'
                    if c not in published[ex]
                    else sample_qc_render(cells[ex][c], c in masked[ex], None)
                )
        out.append(rec)
    t = pd.DataFrame(out)
    path = f'{SAMPLE_QC_TABLE_DIR}/05_supp_table2_sample_qc.csv'
    t.to_csv(path, index=False)
    print(t.to_string(index=False))
    print(f'Write:{path}')


def clinical_covariates_log(msg):
    print(msg, flush=True)


def clinical_covariates_section(title):
    clinical_covariates_log('\n' + '=' * 80)
    clinical_covariates_log(title)
    clinical_covariates_log('=' * 80)


def map_menopause(ans):
    if pd.isna(ans):
        return np.nan
    if 'Yes None' in ans or 'Yes, no menstrual periods' in ans:
        return 'postmenopausal'
    if 'Periods Havent Stopped' in ans or ans.strip().endswith(': No'):
        return 'premenopausal'
    if 'Yes But Hormone' in ans or 'hormone' in ans.lower():
        return 'uncertain_hormone_influenced'
    if 'Not Sure' in ans:
        return 'unknown_not_sure'
    if 'Skip' in ans or 'Prefer Not' in ans or 'Dont Know' in ans:
        return np.nan
    return np.nan


def relative_label(ans):
    if 'Self' in ans:
        return 'Self'
    if 'Skip' in ans or 'Prefer Not' in ans or 'Dont Know' in ans:
        return 'Skip'
    for rel in ['Daughter', 'Father', 'Grandparent', 'Mother', 'Sibling', 'Son']:
        if rel in ans:
            return rel
    return 'Other'


def _famhist_label(pid):
    if pid in famhist_positive_ids:
        return 'positive_nonself'
    if pid in famhist_answered_ids:
        return 'negative_or_self_only'
    return np.nan


def derive_smoking(row):
    a1 = row['ans_100cigs']
    a2 = row['ans_freq']
    if pd.isna(a1) or 'Skip' in str(a1) or 'Prefer Not' in str(a1) or ('Dont Know' in str(a1)):
        return np.nan
    if a1.endswith(': No'):
        return 'never_smoker'
    if a1.endswith(': Yes'):
        if pd.isna(a2):
            return 'ever_smoker_currency_unknown'
        if 'Not At All' in a2:
            return 'former_smoker'
        if 'Every Day' in a2 or 'Some Days' in a2:
            return 'current_smoker'
        return 'ever_smoker_currency_unknown'
    return np.nan


def summarize_continuous(df, var, group_label, n_total):
    s = df[var].dropna()
    n = len(s)
    return {
        'group': group_label,
        'variable': var,
        'variable_type': 'continuous',
        'category': np.nan,
        'n_nonmissing': n,
        'missing_rate': 1 - n / n_total,
        'mean': s.mean() if n else np.nan,
        'sd': s.std() if n else np.nan,
        'median': s.median() if n else np.nan,
        'iqr_low': s.quantile(0.25) if n else np.nan,
        'iqr_high': s.quantile(0.75) if n else np.nan,
        'n_category': np.nan,
        'pct_within_nonmissing': np.nan,
        'pct_of_group_total': np.nan,
    }


def summarize_categorical(df, var, group_label, n_total):
    s = df[var].dropna()
    n_nonmissing = len(s)
    rows = []
    if n_nonmissing == 0:
        rows.append(
            {
                'group': group_label,
                'variable': var,
                'variable_type': 'categorical',
                'category': '(no data)',
                'n_nonmissing': 0,
                'missing_rate': 1.0,
                'mean': np.nan,
                'sd': np.nan,
                'median': np.nan,
                'iqr_low': np.nan,
                'iqr_high': np.nan,
                'n_category': 0,
                'pct_within_nonmissing': np.nan,
                'pct_of_group_total': np.nan,
            }
        )
        return rows
    counts = s.value_counts()
    for cat, n_cat in counts.items():
        rows.append(
            {
                'group': group_label,
                'variable': var,
                'variable_type': 'categorical',
                'category': cat,
                'n_nonmissing': n_nonmissing,
                'missing_rate': 1 - n_nonmissing / n_total,
                'mean': np.nan,
                'sd': np.nan,
                'median': np.nan,
                'iqr_low': np.nan,
                'iqr_high': np.nan,
                'n_category': n_cat,
                'pct_within_nonmissing': n_cat / n_nonmissing,
                'pct_of_group_total': n_cat / n_total,
            }
        )
    return rows


if __name__ == "__main__":
    for selected_task in selected_tasks:

        # ============================================================================
        # 1. Preprocess samples
        # ============================================================================
        if selected_task == 'sample-qc':
            from google.cloud import bigquery as sample_qc_bigquery

            sys.path.insert(0, str(PROJECT_ROOT / 'scripts/functions'))

            from scripts.functions.disclosure import (
                decide as sample_qc_decide,
                render as sample_qc_render,
            )

            GOOGLE_PROJECT = os.environ['GOOGLE_PROJECT']

            SAMPLE_QC_CDR = os.environ['WORKSPACE_CDR']

            SAMPLE_QC_CDR_STORAGE_PATH = os.environ['CDR_STORAGE_PATH']

            REPO = str(PROJECT_ROOT)

            SAMPLE_QC_DATA_DIR = f'{REPO}/data'

            SAMPLE_QC_INTERNAL_DIR = f'{SAMPLE_QC_DATA_DIR}/05_internal'

            SAMPLE_QC_TABLE_DIR = f'{REPO}/results/table'

            SCRATCH = '/tmp/com_rare_BC/qc_slice'

            for sample_qc_d in (SCRATCH, SAMPLE_QC_INTERNAL_DIR, SAMPLE_QC_TABLE_DIR):
                os.makedirs(sample_qc_d, exist_ok=True)

            AUX_QC = f'{SAMPLE_QC_CDR_STORAGE_PATH}/wgs/short_read/snpindel/aux/qc'

            AUX_REL = f'{SAMPLE_QC_CDR_STORAGE_PATH}/wgs/short_read/snpindel/aux/relatedness'

            ACAF_BED_DIR = (
                f'{SAMPLE_QC_CDR_STORAGE_PATH}/wgs/short_read/snpindel/acaf_threshold/plink_bed'
            )

            PLINK2 = '/opt/workbench-tools/binaries/bin/plink2'

            KINSHIP_THRESHOLD_REQUESTED = 0.0884

            SAMPLE_MISSING_THRESHOLD = 0.03

            HWE_P_THRESHOLD = 1e-06

            SAMPLE_QC_COLS = ['case', 'control']

            stats = sample_qc_calculation()

            print('\n' + '=' * 90)

            print('Final statistical summary:')

            for sample_qc_k, sample_qc_v in stats.items():
                print(f'  {sample_qc_k}: {sample_qc_v}')

        # ============================================================================
        # 2. Derive clinical covariates
        # ============================================================================
        if selected_task == 'covariates':
            from google.cloud import bigquery as clinical_covariates_bigquery

            matplotlib.use('Agg')

            plt.style.use('default')

            CLINICAL_COVARIATES_BASE = str(PROJECT_ROOT)

            CLINICAL_COVARIATES_DATA_DIR = os.path.join(CLINICAL_COVARIATES_BASE, 'data')

            CLINICAL_COVARIATES_INTERNAL_DIR = os.path.join(
                CLINICAL_COVARIATES_BASE, 'data', '09_internal'
            )

            CLINICAL_COVARIATES_FIG_DIR = CLINICAL_COVARIATES_TABLE_DIR = REPORT_DIR = (
                CLINICAL_COVARIATES_INTERNAL_DIR
            )

            for clinical_covariates_d in [
                CLINICAL_COVARIATES_DATA_DIR,
                CLINICAL_COVARIATES_INTERNAL_DIR,
            ]:
                os.makedirs(clinical_covariates_d, exist_ok=True)

            CLINICAL_COVARIATES_CDR = os.environ['WORKSPACE_CDR']

            CLINICAL_COVARIATES_PROJECT_ID = os.environ['GOOGLE_PROJECT']

            CLINICAL_COVARIATES_CDR_STORAGE_PATH = os.environ['CDR_STORAGE_PATH']

            client = clinical_covariates_bigquery.Client()

            COHORT_TABLE = f'{CLINICAL_COVARIATES_PROJECT_ID}.my_analysis_results.aim1_case_control_cohort_scratch'

            CLINICAL_COVARIATES_CASE_COLOR = 'tab:blue'

            CLINICAL_COVARIATES_CONTROL_COLOR = 'tab:orange'

            CLINICAL_COVARIATES_GROUP_COLORS = {
                'case': CLINICAL_COVARIATES_CASE_COLOR,
                'control': CLINICAL_COVARIATES_CONTROL_COLOR,
            }

            CLINICAL_COVARIATES_GROUP_ORDER = ['case', 'control']

            CLINICAL_COVARIATES_DPI = 150

            clinical_covariates_section(
                'Step 0: Read Aim 1 analysis queue (No. 05 pass_qc==1, Genetics EUR; starting from 2026-09-24)'
            )

            qcf = pd.read_csv(os.path.join(CLINICAL_COVARIATES_DATA_DIR, '05_cohort_qc_flags.csv'))

            cohort_df = (
                qcf.loc[qcf['pass_qc'] == 1, ['person_id', 'analysis_group']]
                .rename(columns={'analysis_group': 'group'})
                .copy()
            )

            clinical_covariates_log(
                f"cohort（case+control）: {cohort_df.shape[0]}people, distribution:{cohort_df['group'].value_counts().to_dict()}"
            )

            job_config = clinical_covariates_bigquery.LoadJobConfig(
                write_disposition='WRITE_TRUNCATE'
            )

            client.load_table_from_dataframe(
                cohort_df, COHORT_TABLE, job_config=job_config
            ).result()

            clinical_covariates_log(
                f'cohort has been (re)written to the scratch table:{COHORT_TABLE}'
            )

            result = cohort_df.copy()

            clinical_covariates_section(
                'Step 1: BMI (measurement concept_id=3038553, the latest record for each person)'
            )

            bmi_query = f'\n    SELECT person_id, bmi_value, measurement_date\n    FROM (\n      SELECT cc.person_id, m.value_as_number AS bmi_value, m.measurement_date,\n             ROW_NUMBER() OVER (PARTITION BY cc.person_id ORDER BY m.measurement_date DESC) AS rn\n      FROM `{COHORT_TABLE}` cc\n      JOIN `{CLINICAL_COVARIATES_CDR}.measurement` m ON cc.person_id = m.person_id\n      WHERE m.measurement_concept_id = 3038553\n        AND m.value_as_number IS NOT NULL\n        AND m.value_as_number BETWEEN 10 AND 100\n    )\n    WHERE rn = 1\n    '

            bmi_df = client.query(bmi_query).to_dataframe()

            clinical_covariates_log(
                f'BMI non-missing number of people:{bmi_df.shape[0]} / {cohort_df.shape[0]} ({bmi_df.shape[0] / cohort_df.shape[0]:.1%})'
            )

            result = result.merge(bmi_df[['person_id', 'bmi_value']], on='person_id', how='left')

            clinical_covariates_section(
                'Step 2: Menopausal status (ds_survey question_concept_id=1585784, Overall Health Questionnaire)'
            )

            menopause_query = f'\n    SELECT person_id, answer\n    FROM (\n      SELECT cc.person_id, s.answer, s.survey_datetime,\n             ROW_NUMBER() OVER (PARTITION BY cc.person_id ORDER BY s.survey_datetime DESC) AS rn\n      FROM `{COHORT_TABLE}` cc\n      JOIN `{CLINICAL_COVARIATES_CDR}.ds_survey` s ON cc.person_id = s.person_id\n      WHERE s.question_concept_id = 1585784\n    )\n    WHERE rn = 1\n    '

            menopause_df = client.query(menopause_query).to_dataframe()

            clinical_covariates_log(
                f'Menopausal status question Number of people with answer records (including skips):{menopause_df.shape[0]} / {cohort_df.shape[0]} ({menopause_df.shape[0] / cohort_df.shape[0]:.1%})'
            )

            clinical_covariates_log('Original answer distribution:')

            clinical_covariates_log(str(menopause_df['answer'].value_counts()))

            menopause_df['menopause_status'] = menopause_df['answer'].apply(map_menopause)

            result = result.merge(
                menopause_df[['person_id', 'menopause_status']], on='person_id', how='left'
            )

            clinical_covariates_log('Post-mapping distribution (menopause_status):')

            clinical_covariates_log(str(result['menopause_status'].value_counts(dropna=False)))

            clinical_covariates_section(
                'Step 3: Reproductive history (parity/age at menarche) - Search results: No structured field found'
            )

            clinical_covariates_log(
                "Retrieved ds_survey table survey IN ('Personal and Family Health History','Overall Health') under"
            )

            clinical_covariates_log(
                'All questions whose text contains pregnan/menarche/period/live birth/children/parity/hysterectomy/ovary,'
            )

            clinical_covariates_log(
                "Only the binary 'Have you ever been pregnant' (OverallHealth_Pregnancy/1PregnancyStatus) was found, with no parity count field and no menarche age field."
            )

            clinical_covariates_log(
                "Conclusion: The current questionnaire of AoU does not provide structured answers to parity or age at menarche. This script does not use 'have you ever been pregnant' instead? Skip and do not force it."
            )

            pregnancy_query = f'\n    SELECT person_id, answer\n    FROM (\n      SELECT cc.person_id, s.answer, s.survey_datetime,\n             ROW_NUMBER() OVER (PARTITION BY cc.person_id ORDER BY s.survey_datetime DESC) AS rn\n      FROM `{COHORT_TABLE}` cc\n      JOIN `{CLINICAL_COVARIATES_CDR}.ds_survey` s ON cc.person_id = s.person_id\n      WHERE s.question_concept_id = 1585811\n    )\n    WHERE rn = 1\n    '

            pregnancy_df = client.query(pregnancy_query).to_dataframe()

            clinical_covariates_log(
                f"\n (circumstantial evidence, non-main analysis variable) Number of answers to the question 'Have you ever been pregnant':{pregnancy_df.shape[0]} / {cohort_df.shape[0]}"
            )

            if pregnancy_df.shape[0] > 0:
                clinical_covariates_log(str(pregnancy_df['answer'].value_counts()))

            clinical_covariates_section(
                'Step 4: Family history of breast cancer (ds_survey question_concept_id=836772, Personal and Family Health History Questionnaire)'
            )

            famhist_query = f'\n    SELECT cc.person_id, s.answer\n    FROM `{COHORT_TABLE}` cc\n    JOIN `{CLINICAL_COVARIATES_CDR}.ds_survey` s ON cc.person_id = s.person_id\n    WHERE s.question_concept_id = 836772\n    '

            famhist_raw = client.query(famhist_query).to_dataframe()

            clinical_covariates_log(
                f"The total number of lines in this question (including multiple-select multiple lines, including skip):{famhist_raw.shape[0]}, Number of people covered:{famhist_raw['person_id'].nunique()}"
            )

            famhist_raw['relative'] = famhist_raw['answer'].apply(relative_label)

            clinical_covariates_log(
                'Relative type distribution (number of people, duplicate removal):'
            )

            clinical_covariates_log(str(famhist_raw.groupby('relative')['person_id'].nunique()))

            non_self_relatives = {'Daughter', 'Father', 'Grandparent', 'Mother', 'Sibling', 'Son'}

            famhist_positive_ids = set(
                famhist_raw.loc[famhist_raw['relative'].isin(non_self_relatives), 'person_id']
            )

            famhist_answered_ids = set(famhist_raw['person_id'])

            famhist_self_ids = set(famhist_raw.loc[famhist_raw['relative'] == 'Self', 'person_id'])

            result['family_history_breast_cancer'] = result['person_id'].apply(_famhist_label)

            clinical_covariates_log(
                'Distribution after mapping (family_history_breast_cancer; NaN=issue not reached/module not completed):'
            )

            clinical_covariates_log(
                str(result['family_history_breast_cancer'].value_counts(dropna=False))
            )

            case_ids = set(result.loc[result['group'] == 'case', 'person_id'])

            control_ids = set(result.loc[result['group'] == 'control', 'person_id'])

            clinical_covariates_log(
                f'\n (circumstantial evidence) Self-reported breast cancer (Self) and it is an EHR case:{len(famhist_self_ids & case_ids)} / EHR case {len(case_ids)}'
            )

            clinical_covariates_log(
                f'(Circumstantial evidence) Self-reported breast cancer (Self) but under EHR control:{len(famhist_self_ids & control_ids)} / EHR control {len(control_ids)}(The latter should be close to 0. If it is obviously >0, please note that the ascertainment is inconsistent. It only prompts that it does not affect the main analysis definition)'
            )

            clinical_covariates_section(
                'Step 5: Smoking status (ds_survey question 1585857+1585860, Lifestyle questionnaire, two-step method derived never/former/current)'
            )

            smk1_query = f'\n    SELECT person_id, answer AS ans_100cigs\n    FROM (\n      SELECT cc.person_id, s.answer, s.survey_datetime,\n             ROW_NUMBER() OVER (PARTITION BY cc.person_id ORDER BY s.survey_datetime DESC) AS rn\n      FROM `{COHORT_TABLE}` cc\n      JOIN `{CLINICAL_COVARIATES_CDR}.ds_survey` s ON cc.person_id = s.person_id\n      WHERE s.question_concept_id = 1585857\n    )\n    WHERE rn = 1\n    '

            smk2_query = f'\n    SELECT person_id, answer AS ans_freq\n    FROM (\n      SELECT cc.person_id, s.answer, s.survey_datetime,\n             ROW_NUMBER() OVER (PARTITION BY cc.person_id ORDER BY s.survey_datetime DESC) AS rn\n      FROM `{COHORT_TABLE}` cc\n      JOIN `{CLINICAL_COVARIATES_CDR}.ds_survey` s ON cc.person_id = s.person_id\n      WHERE s.question_concept_id = 1585860\n    )\n    WHERE rn = 1\n    '

            smk1_df = client.query(smk1_query).to_dataframe()

            smk2_df = client.query(smk2_query).to_dataframe()

            clinical_covariates_log(
                f"Number of respondents to '100 cigarettes lifetime smoking history':{smk1_df.shape[0]} / {cohort_df.shape[0]}"
            )

            clinical_covariates_log(str(smk1_df['ans_100cigs'].value_counts()))

            clinical_covariates_log(
                f"Number of people who responded to \n'Current Smoking Frequency':{smk2_df.shape[0]} / {cohort_df.shape[0]}"
            )

            clinical_covariates_log(str(smk2_df['ans_freq'].value_counts()))

            smk = smk1_df.merge(smk2_df, on='person_id', how='outer')

            smk['smoking_status'] = smk.apply(derive_smoking, axis=1)

            result = result.merge(smk[['person_id', 'smoking_status']], on='person_id', how='left')

            clinical_covariates_log('Derived distribution (smoking_status):')

            clinical_covariates_log(str(result['smoking_status'].value_counts(dropna=False)))

            clinical_covariates_section(
                'Step 6: ER/PR/HER2 receptor status (condition_occurrence, sparsity statistics)'
            )

            RECEPTOR_CONCEPTS = {
                'ER_positive': [4167696, 44837878, 3526950, 3283803],
                'ER_negative': [4261933, 44834368, 3564533, 3183458],
                'PR_positive': [4219694, 3526951, 3315165, 1102831],
                'PR_negative': [4261934, 3562918, 3254239, 1102832],
                'HER2_positive': [4142116, 42542268, 42542428, 45949645],
            }

            receptor_rows = []

            for clinical_covariates_label, ids in RECEPTOR_CONCEPTS.items():
                idlist = ','.join((str(i) for i in ids))
                q = f'\n        SELECT cc.group, COUNT(DISTINCT co.person_id) n_person\n        FROM `{CLINICAL_COVARIATES_CDR}.condition_occurrence` co\n        JOIN `{COHORT_TABLE}` cc ON co.person_id = cc.person_id\n        WHERE co.condition_concept_id IN ({idlist})\n        GROUP BY 1\n        '
                clinical_covariates_df = (
                    client.query(q).to_dataframe().set_index('group')['n_person'].to_dict()
                )
                n_case = clinical_covariates_df.get('case', 0)
                n_control = clinical_covariates_df.get('control', 0)
                receptor_rows.append(
                    {
                        'receptor_status': clinical_covariates_label,
                        'n_case': n_case,
                        'pct_of_case': n_case / cohort_df[cohort_df.group == 'case'].shape[0],
                        'n_control': n_control,
                        'pct_of_control': n_control
                        / cohort_df[cohort_df.group == 'control'].shape[0],
                    }
                )
                clinical_covariates_log(
                    f"{clinical_covariates_label}: case n={n_case} ({n_case / cohort_df[cohort_df.group == 'case'].shape[0]:.2%}), control n={n_control} ({n_control / cohort_df[cohort_df.group == 'control'].shape[0]:.4%})"
                )

            receptor_summary = pd.DataFrame(receptor_rows)

            receptor_path = os.path.join(
                CLINICAL_COVARIATES_TABLE_DIR, '09_receptor_status_sparsity.csv'
            )

            receptor_summary.to_csv(receptor_path, index=False)

            clinical_covariates_log(f'\n has written:{receptor_path}')

            clinical_covariates_log(
                'Conclusion: ER status has structured records in about 15% of the case group, and PR/HER2 coverage is <1% - structured pathological receptor dataIt is indeed extremely sparse in EHR sources and is not recommended as a covariate in the main analysis. It is only for reference in future sensitivity stratified analyzes (if data availability improves).'
            )

            clinical_covariates_section(
                'Step 7: ancestry principal components PC1-10 (GCS streaming filter ancestry_preds.tsv, 210MB full file without disk)'
            )

            ANCESTRY_PATH = f'{CLINICAL_COVARIATES_CDR_STORAGE_PATH}/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv'

            clinical_covariates_log(f'document:{ANCESTRY_PATH}')

            cohort_id_set = set(cohort_df['person_id'].astype(str))

            pc_rows = []

            proc = subprocess.Popen(
                ['gsutil', '-u', CLINICAL_COVARIATES_PROJECT_ID, 'cat', ANCESTRY_PATH],
                stdout=subprocess.PIPE,
                text=True,
                bufsize=1,
            )

            header = None

            n_lines = 0

            for line in proc.stdout:
                n_lines += 1
                parts = line.rstrip('\n').split('\t')
                if header is None:
                    header = parts
                    continue
                research_id = parts[0]
                if research_id not in cohort_id_set:
                    continue
                try:
                    pca_features = json.loads(parts[3])
                except (json.JSONDecodeError, IndexError):
                    continue
                row = {'person_id': int(research_id), 'ancestry_pred': parts[1]}
                for clinical_covariates_k in range(10):
                    row[f'PC{clinical_covariates_k + 1}'] = (
                        pca_features[clinical_covariates_k]
                        if clinical_covariates_k < len(pca_features)
                        else np.nan
                    )
                pc_rows.append(row)

            proc.wait()

            clinical_covariates_log(
                f'Number of lines scanned for the entire file:{n_lines}, the number of people hitting this queue:{len(pc_rows)} / {cohort_df.shape[0]}'
            )

            pc_df = pd.DataFrame(pc_rows)

            result = result.merge(pc_df, on='person_id', how='left')

            clinical_covariates_log(
                'ancestry_pred (AoU official predicted ancestry label) distribution (in this queue, it should be almost all eur, because the queue definition is self-reported white):'
            )

            clinical_covariates_log(str(result['ancestry_pred'].value_counts(dropna=False)))

            clinical_covariates_section(
                "Step 8: Sequencing batch/platform proxy variables (GCS streaming filter genomics_metrics.tsv, no literal 'batch' field, approximated by site_id)"
            )

            QC_PATH = f'{CLINICAL_COVARIATES_CDR_STORAGE_PATH}/wgs/short_read/snpindel/aux/qc/genomics_metrics.tsv'

            clinical_covariates_log(f'document:{QC_PATH}')

            batch_rows = []

            proc = subprocess.Popen(
                ['gsutil', '-u', CLINICAL_COVARIATES_PROJECT_ID, 'cat', QC_PATH],
                stdout=subprocess.PIPE,
                text=True,
                bufsize=1,
            )

            header = None

            n_lines = 0

            for line in proc.stdout:
                n_lines += 1
                parts = line.rstrip('\n').split('\t')
                if header is None:
                    header = parts
                    col_idx = {
                        clinical_covariates_c: i for i, clinical_covariates_c in enumerate(header)
                    }
                    continue
                research_id = parts[0]
                if research_id not in cohort_id_set:
                    continue
                row = {
                    'person_id': int(research_id),
                    'seq_site_id': parts[col_idx['site_id']],
                    'seq_sequencer': parts[col_idx['sequencer']],
                    'seq_collection_date': parts[col_idx['biosample_collection_date']],
                }
                batch_rows.append(row)

            proc.wait()

            clinical_covariates_log(
                f'Number of lines scanned for the entire file:{n_lines}, the number of people hitting this queue:{len(batch_rows)} / {cohort_df.shape[0]}'
            )

            batch_df = pd.DataFrame(batch_rows)

            batch_df['seq_collection_year'] = pd.to_datetime(
                batch_df['seq_collection_date'], errors='coerce'
            ).dt.year

            result = result.merge(
                batch_df[['person_id', 'seq_site_id', 'seq_sequencer', 'seq_collection_year']],
                on='person_id',
                how='left',
            )

            clinical_covariates_log(
                'seq_site_id (sequencing center, used as batch proxy variable) distribution:'
            )

            clinical_covariates_log(str(result['seq_site_id'].value_counts(dropna=False)))

            clinical_covariates_log('seq_sequencer (sequencer model) distribution:')

            clinical_covariates_log(str(result['seq_sequencer'].value_counts(dropna=False)))

            clinical_covariates_log(
                'seq_collection_year (sample collection year, batch time dimension reference) distribution:'
            )

            clinical_covariates_log(
                str(result['seq_collection_year'].value_counts(dropna=False).sort_index())
            )

            clinical_covariates_section('Step 9: Place person-level intermediate data')

            clinical_cols = [
                'person_id',
                'group',
                'bmi_value',
                'menopause_status',
                'family_history_breast_cancer',
                'smoking_status',
                'seq_site_id',
                'seq_sequencer',
                'seq_collection_year',
            ]

            pc_cols = [
                clinical_covariates_c
                for clinical_covariates_c in result.columns
                if clinical_covariates_c.startswith('PC')
            ]

            pc_cols_out = ['person_id', 'group', 'ancestry_pred'] + pc_cols

            clinical_path = os.path.join(
                CLINICAL_COVARIATES_DATA_DIR, '09_extended_clinical_person_level.csv.gz'
            )

            clinical_out = result[clinical_cols].copy()

            clinical_out['bmi_value'] = clinical_out['bmi_value'].round(2)

            clinical_out.to_csv(clinical_path, index=False, compression='gzip')

            clinical_covariates_size_mb = os.path.getsize(clinical_path) / 1000000.0

            clinical_covariates_log(
                f'Has written:{clinical_path} ({clinical_out.shape[0]}row x{clinical_out.shape[1]}List,{clinical_covariates_size_mb:.3f} MB)'
            )

            if clinical_covariates_size_mb >= 5:
                clinical_covariates_log(
                    '!!! WARNING: Exceeding 5MB, you need to streamline the list and re-place it!!!'
                )

            pc_path = os.path.join(
                CLINICAL_COVARIATES_DATA_DIR, '09_ancestry_pc_person_level.csv.gz'
            )

            pc_out = result[pc_cols_out].copy()

            pc_out[pc_cols] = pc_out[pc_cols].round(4)

            pc_out.to_csv(pc_path, index=False, compression='gzip')

            clinical_covariates_size_mb = os.path.getsize(pc_path) / 1000000.0

            clinical_covariates_log(
                f'Has written:{pc_path} ({pc_out.shape[0]}row x{pc_out.shape[1]}List,{clinical_covariates_size_mb:.3f} MB)'
            )

            if clinical_covariates_size_mb >= 5:
                clinical_covariates_log(
                    '!!! WARNING: Exceeding 5MB, you need to streamline the list and re-place it!!!'
                )

            clinical_covariates_section('Step 10: Extended baseline characteristics table')

            clinical_covariates_case_df = result[result['group'] == 'case']

            clinical_covariates_control_df = result[result['group'] == 'control']

            n_case_total = clinical_covariates_case_df.shape[0]

            n_control_total = clinical_covariates_control_df.shape[0]

            CLINICAL_COVARIATES_CONT_VARS = [
                'bmi_value',
                'PC1',
                'PC2',
                'PC3',
                'PC4',
                'PC5',
                'PC6',
                'PC7',
                'PC8',
                'PC9',
                'PC10',
            ]

            CAT_VARS = [
                'menopause_status',
                'family_history_breast_cancer',
                'smoking_status',
                'ancestry_pred',
                'seq_site_id',
                'seq_sequencer',
                'seq_collection_year',
            ]

            table_rows = []

            for group_label, clinical_covariates_df, n_total in [
                ('case', clinical_covariates_case_df, n_case_total),
                ('control', clinical_covariates_control_df, n_control_total),
            ]:
                table_rows.append(
                    {
                        'group': group_label,
                        'variable': 'n_total',
                        'variable_type': 'meta',
                        'category': np.nan,
                        'n_nonmissing': n_total,
                        'missing_rate': 0.0,
                        'mean': np.nan,
                        'sd': np.nan,
                        'median': np.nan,
                        'iqr_low': np.nan,
                        'iqr_high': np.nan,
                        'n_category': np.nan,
                        'pct_within_nonmissing': np.nan,
                        'pct_of_group_total': np.nan,
                    }
                )
                for clinical_covariates_v in CLINICAL_COVARIATES_CONT_VARS:
                    table_rows.append(
                        summarize_continuous(
                            clinical_covariates_df, clinical_covariates_v, group_label, n_total
                        )
                    )
                for clinical_covariates_v in CAT_VARS:
                    table_rows.extend(
                        summarize_categorical(
                            clinical_covariates_df, clinical_covariates_v, group_label, n_total
                        )
                    )

            ext_table = pd.DataFrame(table_rows)

            for clinical_covariates_col in [
                'mean',
                'sd',
                'median',
                'iqr_low',
                'iqr_high',
                'missing_rate',
                'pct_within_nonmissing',
                'pct_of_group_total',
            ]:
                ext_table[clinical_covariates_col] = ext_table[clinical_covariates_col].round(4)

            ext_table_path = os.path.join(
                CLINICAL_COVARIATES_TABLE_DIR, '09_extended_baseline_characteristics.csv'
            )

            ext_table.to_csv(ext_table_path, index=False)

            clinical_covariates_log(f'Has written:{ext_table_path} ({ext_table.shape[0]}OK)')

            clinical_covariates_log(ext_table.to_string(index=False))

            clinical_covariates_section('Step 11: Chart')

            clinical_covariates_fig, clinical_covariates_ax = plt.subplots(figsize=(5, 4.5))

            clinical_covariates_plot_df = result[result['group'].isin(['case', 'control'])][
                ['group', 'bmi_value']
            ].dropna()

            sns.boxplot(
                x='group',
                y='bmi_value',
                data=clinical_covariates_plot_df,
                order=CLINICAL_COVARIATES_GROUP_ORDER,
                hue='group',
                palette=CLINICAL_COVARIATES_GROUP_COLORS,
                legend=False,
                ax=clinical_covariates_ax,
                showfliers=False,
            )

            clinical_covariates_ax.set_title(
                f"BMI: case vs control\n(n_case={clinical_covariates_plot_df[clinical_covariates_plot_df.group == 'case'].shape[0]}, n_control={clinical_covariates_plot_df[clinical_covariates_plot_df.group == 'control'].shape[0]}, outlier points hidden)"
            )

            clinical_covariates_ax.set_xlabel('')

            clinical_covariates_ax.set_ylabel('BMI (kg/m^2)')

            clinical_covariates_ax.grid(False)

            clinical_covariates_fig.tight_layout()

            bmi_path = os.path.join(
                CLINICAL_COVARIATES_FIG_DIR, '09_bmi_boxplot_case_vs_control.png'
            )

            clinical_covariates_fig.savefig(bmi_path, dpi=CLINICAL_COVARIATES_DPI)

            plt.close(clinical_covariates_fig)

            clinical_covariates_log(f'Has written:{bmi_path}')

            cat_plot_specs = [
                (
                    'menopause_status',
                    'Menopause status (proportion within group, non-missing only)',
                    '09_menopause_status_barplot.png',
                ),
                (
                    'smoking_status',
                    'Smoking status (proportion within group, non-missing only)',
                    '09_smoking_status_barplot.png',
                ),
                (
                    'family_history_breast_cancer',
                    'Family history of breast cancer, non-self (proportion within group, non-missing only)',
                    '09_family_history_barplot.png',
                ),
            ]

            for clinical_covariates_var, title, clinical_covariates_fname in cat_plot_specs:
                clinical_covariates_plot_df = (
                    result[result['group'].isin(['case', 'control'])][
                        ['group', clinical_covariates_var]
                    ]
                    .dropna()
                    .copy()
                )
                if clinical_covariates_plot_df.empty:
                    clinical_covariates_log(
                        f'{clinical_covariates_var}: Nothing more than missing data, skip drawing'
                    )
                    continue
                prop_df = (
                    clinical_covariates_plot_df.groupby(['group', clinical_covariates_var])
                    .size()
                    .reset_index(name='n')
                )
                prop_df['total_in_group'] = prop_df.groupby('group')['n'].transform('sum')
                prop_df['proportion'] = prop_df['n'] / prop_df['total_in_group']
                clinical_covariates_fig, clinical_covariates_ax = plt.subplots(figsize=(7, 4.5))
                sns.barplot(
                    x=clinical_covariates_var,
                    y='proportion',
                    hue='group',
                    data=prop_df,
                    order=sorted(prop_df[clinical_covariates_var].unique()),
                    hue_order=CLINICAL_COVARIATES_GROUP_ORDER,
                    palette=CLINICAL_COVARIATES_GROUP_COLORS,
                    ax=clinical_covariates_ax,
                )
                clinical_covariates_ax.set_title(title)
                clinical_covariates_ax.set_xlabel('')
                clinical_covariates_ax.set_ylabel('proportion within group (non-missing)')
                clinical_covariates_ax.tick_params(axis='x', rotation=30)
                clinical_covariates_ax.grid(False)
                clinical_covariates_fig.tight_layout()
                clinical_covariates_out_path = os.path.join(
                    CLINICAL_COVARIATES_FIG_DIR, clinical_covariates_fname
                )
                clinical_covariates_fig.savefig(
                    clinical_covariates_out_path, dpi=CLINICAL_COVARIATES_DPI
                )
                plt.close(clinical_covariates_fig)
                clinical_covariates_log(f'Has written:{clinical_covariates_out_path}')

            clinical_covariates_plot_df = result[result['group'].isin(['case', 'control'])][
                ['group', 'PC1', 'PC2']
            ].dropna()

            clinical_covariates_fig, clinical_covariates_ax = plt.subplots(figsize=(6.5, 5.5))

            for clinical_covariates_g in CLINICAL_COVARIATES_GROUP_ORDER:
                sub = clinical_covariates_plot_df[
                    clinical_covariates_plot_df['group'] == clinical_covariates_g
                ]
                clinical_covariates_ax.scatter(
                    sub['PC1'],
                    sub['PC2'],
                    s=4,
                    alpha=0.25,
                    color=CLINICAL_COVARIATES_GROUP_COLORS[clinical_covariates_g],
                    label=clinical_covariates_g,
                )

            clinical_covariates_ax.set_title(
                f"Ancestry PC1 vs PC2 (AoU official pre-computed, case vs control)\nn_case={clinical_covariates_plot_df[clinical_covariates_plot_df.group == 'case'].shape[0]}, n_control={clinical_covariates_plot_df[clinical_covariates_plot_df.group == 'control'].shape[0]}"
            )

            clinical_covariates_ax.set_xlabel('PC1')

            clinical_covariates_ax.set_ylabel('PC2')

            leg = clinical_covariates_ax.legend(markerscale=4)

            clinical_covariates_ax.grid(False)

            clinical_covariates_fig.tight_layout()

            pc_scatter_path = os.path.join(
                CLINICAL_COVARIATES_FIG_DIR, '09_pc1_vs_pc2_scatter_case_vs_control.jpg'
            )

            clinical_covariates_fig.savefig(
                pc_scatter_path, dpi=CLINICAL_COVARIATES_DPI, pil_kwargs={'quality': 90}
            )

            plt.close(clinical_covariates_fig)

            clinical_covariates_log(f'Has written:{pc_scatter_path}')

            clinical_covariates_section('Step 12: Output file size check (required < 5MB)')

            for clinical_covariates_d in [
                CLINICAL_COVARIATES_FIG_DIR,
                CLINICAL_COVARIATES_TABLE_DIR,
                CLINICAL_COVARIATES_DATA_DIR,
            ]:
                for clinical_covariates_f in sorted(os.listdir(clinical_covariates_d)):
                    if not clinical_covariates_f.startswith('09_'):
                        continue
                    clinical_covariates_path = os.path.join(
                        clinical_covariates_d, clinical_covariates_f
                    )
                    if os.path.isfile(clinical_covariates_path):
                        clinical_covariates_size_mb = (
                            os.path.getsize(clinical_covariates_path) / 1000000.0
                        )
                        clinical_covariates_flag = (
                            'OK' if clinical_covariates_size_mb < 5 else '!!! More than 5MB !!!'
                        )
                        clinical_covariates_log(
                            f'{clinical_covariates_path}: {clinical_covariates_size_mb:.3f} MB [{clinical_covariates_flag}]'
                        )

            clinical_covariates_log('Finish.')

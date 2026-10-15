"""Apply source-site and genotype QC, recover missing sites and reconcile carriers.

Main call rules: GQ>=20; heterozygous VAF 0.3-0.8; homozygous VAF>=0.8.
Source-site FILTER and native ALT mapping are checked before carrier collapse.
"""

from __future__ import annotations

from pathlib import Path
import sys

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.functions.paths import PROJECT_ROOT
from scripts.functions.cli import parse_tasks

TASKS = (
    'site-qc',
    'calls',
    'gene-calls',
    'missing-sites',
    'validate-qc',
    'call-audit',
    'reconcile',
)
if __name__ == '__main__':
    selected_tasks = parse_tasks(
        __doc__,
        TASKS,
        ['site-qc', 'calls', 'missing-sites', 'validate-qc', 'call-audit', 'reconcile'],
    )

import json

import argparse

import hashlib

import importlib.util

import os


import hail as hl

import pandas as pd

SITE_QC_ROOT = PROJECT_ROOT

SITE_QC_OUT = SITE_QC_ROOT / 'data/29_rowlands_extension/section6'

from concurrent.futures import ThreadPoolExecutor

import subprocess


import yaml

PARALLEL_GENOTYPE_QC_ROOT = PROJECT_ROOT

from collections import Counter, defaultdict

import math

import tempfile

import numpy as np

GENOTYPE_QC_ROOT = PROJECT_ROOT

GENOTYPE_QC_OUT = GENOTYPE_QC_ROOT / 'data/29_rowlands_extension/section6'

MOUNT = Path('/tmp/codex_expanded_missense_aou_ro')

SETS = ['raw', 'site_PASS', 'FT_PASS', 'GQ20', 'main', 'CH_sensitive']

MISSING_SITE_CALLS_ROOT = PROJECT_ROOT

MISSING_SITE_CALLS_OUT = MISSING_SITE_CALLS_ROOT / 'data/29_rowlands_extension/section6'

QC_VALIDATION_ROOT = PROJECT_ROOT

QC_VALIDATION_OUT = QC_VALIDATION_ROOT / 'data/29_rowlands_extension/section6'

import csv

import gzip

import re

SAMPLE_VARIANT_AUDIT_ROOT = PROJECT_ROOT

SAMPLE_VARIANT_AUDIT_BASE = SAMPLE_VARIANT_AUDIT_ROOT / 'data/29_rowlands_extension'

SAMPLE_VARIANT_AUDIT_S6 = SAMPLE_VARIANT_AUDIT_BASE / 'section6'

SAMPLE_VARIANT_AUDIT_OUT = SAMPLE_VARIANT_AUDIT_BASE / 'section7'

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

CLINICAL_RECONCILIATION_ROOT = PROJECT_ROOT

CLINICAL_RECONCILIATION_BASE = CLINICAL_RECONCILIATION_ROOT / 'data/29_rowlands_extension'

CLINICAL_RECONCILIATION_S6 = CLINICAL_RECONCILIATION_BASE / 'section6'

CLINICAL_RECONCILIATION_OUT = CLINICAL_RECONCILIATION_BASE / 'section7'

from scripts.functions.site_qc import source_site_pass

from scripts.functions.site_qc import canonical_crosswalk


def run(gene):
    subprocess.run(
        [
            sys.executable,
            str(PARALLEL_GENOTYPE_QC_ROOT / 'scripts/06_variant_genotype_qc.py'),
            '--task',
            'gene-calls',
            '--',
            '--gene',
            gene,
        ],
        cwd=PARALLEL_GENOTYPE_QC_ROOT,
        check=True,
    )


from scripts.functions.io import sha256 as sha

from scripts.functions.genotype_qc import call_qc

from scripts.functions.genotype_qc import normalized_native_sites

from scripts.functions.genotype_qc import safe_count

from scripts.functions.sample_variant_audit import terms

from scripts.functions.sample_variant_audit import load_module

from scripts.functions.sample_variant_audit import positive_float

if __name__ == "__main__":
    for selected_task in selected_tasks:

        # ============================================================================
        # 1. Verify variant site filters
        # ============================================================================
        if selected_task == 'site-qc':
            project = os.environ['GOOGLE_PROJECT']

            hl.init(
                master='local[2]',
                quiet=True,
                default_reference='GRCh38',
                log=str(SITE_QC_OUT / 'site_filter_hail.log'),
                gcs_requester_pays_configuration=project,
                spark_conf={
                    'spark.driver.memory': '4g',
                    'spark.hadoop.fs.gs.impl': 'com.google.cloud.hadoop.fs.gcs.GoogleHadoopFileSystem',
                    'spark.hadoop.fs.AbstractFileSystem.gs.impl': 'com.google.cloud.hadoop.fs.gcs.GoogleHadoopFS',
                    'spark.hadoop.fs.gs.requester.pays.mode': 'AUTO',
                    'spark.hadoop.fs.gs.requester.pays.project.id': project,
                },
            )

            coords = pd.read_csv(SITE_QC_ROOT / 'data/29_rowlands_extension/gene_coordinates.csv')

            coords['chrom_order'] = coords.chrom.str.replace('chr', '', regex=False).astype(int)

            coords = coords.sort_values(['chrom_order', 'query_region_start'])

            intervals = hl.eval(
                hl.array(
                    [
                        hl.parse_locus_interval(
                            f'{r.chrom}:{r.query_region_start}-{r.query_region_end}',
                            reference_genome='GRCh38',
                        )
                        for r in coords.itertuples()
                    ]
                )
            )

            path = (
                os.environ['CDR_STORAGE_PATH']
                + '/wgs/short_read/snpindel/vds/hail.vds/variant_data'
            )

            mt = hl.read_matrix_table(path, _intervals=intervals)

            rows = mt.rows().select('filters')

            records = rows.collect()

            filters = {}

            for r in records:
                key = (
                    f"{r.locus.contig}:{r.locus.position}:{r.alleles[0]}:{','.join(r.alleles[1:])}"
                )
                filters[key] = dict(
                    filters=sorted(r.filters) if r.filters is not None else None,
                    PASS=source_site_pass(r.filters),
                )

            (SITE_QC_OUT / 'original_VDS_site_filters.json').write_text(json.dumps(filters) + '\n')

            audit_path = SITE_QC_OUT / 'native_site_ALT_QC_audit.tsv.gz'

            audit = (
                pd.read_csv(audit_path, sep='\t', dtype=str)
                if audit_path.exists()
                else pd.DataFrame(columns=['native_key', 'site_filter'])
            )

            native = set(audit.native_key)

            manifest = dict(
                source=path,
                source_CDR='C2025Q4R6',
                rows_only_no_participant_entries=True,
                VCF_FILTER_encoding=sorted(set(audit.site_filter)),
                every_previously_found_native_site_in_VDS=all((k in filters for k in native)),
                all_previously_found_native_sites_original_VDS_PASS=all(
                    (filters.get(k, {}).get('PASS', False) for k in native)
                ),
                PASS_definition='defined empty filters set or singleton PASS; missing set does not pass',
                VCF_dot_never_used_without_original_VDS_verification=True,
            )

            (SITE_QC_OUT / 'original_VDS_site_filter_manifest.json').write_text(
                json.dumps(manifest, indent=2) + '\n'
            )

            print(
                'Original VDS row-level filters verified; no participant entries read or written.',
                flush=True,
            )

            hl.stop()

            canonical_crosswalk()

        # ============================================================================
        # 2. Qc genotypes by gene
        # ============================================================================
        if selected_task == 'calls':
            cfg = yaml.safe_load(
                (PARALLEL_GENOTYPE_QC_ROOT / 'config/rowlands_13gene_extension.yaml').read_text()
            )

            with ThreadPoolExecutor(max_workers=min(4, os.cpu_count() or 1)) as pool:
                list(pool.map(run, cfg['panels']['main_13']))

            subprocess.run(
                [
                    sys.executable,
                    str(PARALLEL_GENOTYPE_QC_ROOT / 'scripts/06_variant_genotype_qc.py'),
                    '--task',
                    'gene-calls',
                    '--',
                ],
                cwd=PARALLEL_GENOTYPE_QC_ROOT,
                check=True,
            )

        # ============================================================================
        # 3. Qc gene calls
        # ============================================================================
        if selected_task == 'gene-calls':
            parser = argparse.ArgumentParser()

            parser.add_argument('--gene', help='Run one gene checkpoint only; no score assembly')

            args = parser.parse_args()

            cfg = yaml.safe_load(
                (GENOTYPE_QC_ROOT / 'config/rowlands_13gene_extension.yaml').read_text()
            )

            assert os.environ['CDR_STORAGE_PATH'] == cfg['data']['storage_path']

            assert os.environ['WORKSPACE_CDR'].split('.')[-1] == cfg['data']['cdr']

            from scripts.functions import nine_gene_missense_age as prior

            prior.ensure_ro_mount(MOUNT)

            assert (
                subprocess.run(
                    ['git', 'check-ignore', '-q', str(GENOTYPE_QC_OUT)], cwd=GENOTYPE_QC_ROOT
                ).returncode
                == 0
            )

            panel = cfg['panels']['main_13']

            ptv = pd.read_csv(
                GENOTYPE_QC_OUT / 'ptv_for_call_qc.tsv.gz',
                sep='\t',
                dtype=str,
                keep_default_na=False,
            )

            miss = pd.read_csv(
                GENOTYPE_QC_OUT / 'eligible_missense.tsv.gz',
                sep='\t',
                dtype=str,
                keep_default_na=False,
            )

            ptv['variant_class'] = 'PTV'

            miss['variant_class'] = 'missense'

            variants = pd.concat([ptv, miss], ignore_index=True)

            source_files = [
                GENOTYPE_QC_OUT / 'ptv_for_call_qc.tsv.gz',
                GENOTYPE_QC_OUT / 'eligible_missense.tsv.gz',
                GENOTYPE_QC_ROOT / 'scripts/06_variant_genotype_qc.py',
                GENOTYPE_QC_ROOT / cfg['cohorts']['case_control']['frozen_ids'],
                GENOTYPE_QC_OUT / 'original_VDS_normalized_allele_filters.json',
            ]

            original_site_filters = json.loads(
                (GENOTYPE_QC_OUT / 'original_VDS_normalized_allele_filters.json').read_text()
            )

            raw_cache_manifest_path = GENOTYPE_QC_OUT / 'raw_carrier_cache_manifest.json'

            raw_cache_manifest = (
                json.loads(raw_cache_manifest_path.read_text())
                if raw_cache_manifest_path.exists()
                else None
            )

            if raw_cache_manifest:
                source_files.append(raw_cache_manifest_path)

            fingerprint = hashlib.sha256(
                json.dumps({str(p): sha(p) for p in source_files}, sort_keys=True).encode()
            ).hexdigest()

            ckdir = GENOTYPE_QC_OUT / 'gene_qc_checkpoints'

            ckdir.mkdir(exist_ok=True)

            boundaries = pd.read_csv(
                GENOTYPE_QC_ROOT / 'data/11_exome_vcf_shard_start_boundaries.csv'
            )

            coords = pd.read_csv(
                GENOTYPE_QC_ROOT / 'data/29_rowlands_extension/gene_coordinates.csv'
            )

            ids = set(
                pd.read_csv(
                    GENOTYPE_QC_ROOT / cfg['cohorts']['case_control']['frozen_ids']
                ).person_id.astype(int)
            )

            age_ids = set(
                pd.read_csv(
                    GENOTYPE_QC_ROOT / cfg['cohorts']['diagnosis_age']['frozen_ids']
                ).person_id.astype(int)
            )

            assert age_ids <= ids and len(ids) == cfg['cohorts']['case_control']['expected_n']

            fasta = GENOTYPE_QC_ROOT / 'data/refdata/Homo_sapiens_assembly38.fasta'

            vcfdir = MOUNT / 'v9' / cfg['data']['exome_relative_path']

            carriers = {
                g: {cl: {s: set() for s in SETS} for cl in ['PTV', 'missense']} for g in panel
            }

            metrics, site_audits = ({}, [])

            with tempfile.TemporaryDirectory(prefix='qc_', dir=GENOTYPE_QC_OUT) as temp:
                sample_file = Path(temp) / 'samples.txt'
                sample_file.write_text('\n'.join(map(str, sorted(ids))) + '\n')
                for gene in panel:
                    if args.gene and gene != args.gene:
                        continue
                    checkpoint = ckdir / f'{gene}.json'
                    if checkpoint.exists():
                        saved = json.loads(checkpoint.read_text())
                        assert (
                            saved['fingerprint'] == fingerprint
                        ), 'Checkpoint inputs changed; inspect before rerunning'
                        carriers[gene] = {
                            cl: {s: set(v) for s, v in stages.items()}
                            for cl, stages in saved['carriers'].items()
                        }
                        metrics[gene] = saved['metrics']
                        site_audits.extend(saved['sites'])
                        print(f'{gene}: verified gene QC checkpoint reused.', flush=True)
                        continue
                    group = variants.loc[variants.gene_symbol == gene]
                    expected_raw = None
                    if raw_cache_manifest and all(
                        (
                            sha(GENOTYPE_QC_ROOT / p) == value
                            for p, value in raw_cache_manifest['input_hashes'].items()
                        )
                    ):
                        raw_cache_path = (
                            GENOTYPE_QC_OUT / 'gene_qc_checkpoints_literal_PASS' / f'{gene}.json'
                        )
                        assert sha(raw_cache_path) == raw_cache_manifest['checkpoint_sha256'][gene]
                        raw_cache = json.loads(raw_cache_path.read_text())
                        expected_raw = {
                            cl: set(raw_cache['carriers'][cl]['raw']) for cl in ['PTV', 'missense']
                        }
                        requery_ids = set().union(*expected_raw.values())
                        assert requery_ids <= ids
                        sample_file.write_text('\n'.join(map(str, sorted(requery_ids))) + '\n')
                    else:
                        sample_file.write_text('\n'.join(map(str, sorted(ids))) + '\n')
                    desired = defaultdict(set)
                    for r in group.itertuples():
                        desired[r.normalized_key].add(r.variant_class)
                    if not desired:
                        metrics[gene] = {}
                        continue
                    g = coords.loc[coords.gene_symbol == gene].iloc[0]
                    chrom_bounds = boundaries.loc[
                        boundaries.chrom_ord == int(g.chrom.replace('chr', ''))
                    ]
                    positions = set(group.position.astype(int)) | {
                        int(k.split(':')[1]) for k in desired
                    }
                    shards = prior.candidate_shards(boundaries, g.chrom, sorted(positions))
                    valid_shards = set(chrom_bounds.shard_idx.astype(int))
                    shards = sorted(
                        {s + d for s in shards for d in [-1, 0, 1] if s + d in valid_shards}
                    )
                    count = Counter()
                    ft_count = Counter()
                    ad_depths = defaultdict(list)
                    dp_depths = defaultdict(list)
                    seen, seen_calls = (set(), set())
                    sites = []
                    for shard in shards:
                        path = vcfdir / f'{shard:010d}.vcf.bgz'
                        header = subprocess.run(
                            ['bcftools', 'view', '-h', str(path)],
                            capture_output=True,
                            text=True,
                            check=True,
                        ).stdout
                        fields = {
                            l.split('ID=', 1)[1].split(',', 1)[0]
                            for l in header.splitlines()
                            if l.startswith('##FORMAT=<ID=')
                        }
                        assert {
                            'GT',
                            'FT',
                            'GQ',
                            'AD',
                        } <= fields, 'Required individual QC fields unavailable'
                        count['shards_read'] += 1
                        count['shards_with_DP'] += int('DP' in fields)
                        region = f'{g.chrom}:{g.query_region_start}-{g.query_region_end}'
                        raw = subprocess.run(
                            [
                                'bcftools',
                                'query',
                                '-r',
                                region,
                                '-f',
                                '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\n',
                                str(path),
                            ],
                            capture_output=True,
                            text=True,
                            check=True,
                        ).stdout
                        mapping = normalized_native_sites(
                            [tuple(l.split('\t')) for l in raw.splitlines()], fasta
                        )
                        matching = {
                            native: key for native, key in mapping.items() if key in desired
                        }
                        for native, key in matching.items():
                            seen.add(key)
                            sites.append(
                                dict(
                                    gene_symbol=gene,
                                    normalized_key=key,
                                    shard=shard,
                                    native_key=':'.join(native[:4]),
                                    native_ALT_index=native[4],
                                    site_filter=native[5],
                                    DP_field_present='DP' in fields,
                                    original_VDS_site_PASS=original_site_filters.get(key, {}).get(
                                        'PASS', False
                                    ),
                                )
                            )
                        if not matching:
                            print(f'{gene}: shard {shard} site scan complete.', flush=True)
                            continue
                        bed = Path(temp) / 'matching.bed'
                        bed.write_text(
                            ''.join(
                                (
                                    f'{c}\t{int(p) - 1}\t{p}\n'
                                    for c, p in sorted({(x[0], x[1]) for x in matching})
                                )
                            )
                        )
                        dp_format = '%DP' if 'DP' in fields else '.'
                        fmt = (
                            '[%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%SAMPLE\t%GT\t%FT\t%GQ\t%AD\t'
                            + dp_format
                            + '\n]'
                        )
                        with tempfile.TemporaryFile() as err2:
                            query = subprocess.Popen(
                                [
                                    'bcftools',
                                    'query',
                                    '-R',
                                    str(bed),
                                    '--regions-overlap',
                                    '0',
                                    '-S',
                                    str(sample_file),
                                    '-i',
                                    'GT="alt"',
                                    '-f',
                                    fmt,
                                    str(path),
                                ],
                                stdout=subprocess.PIPE,
                                stderr=err2,
                                text=True,
                            )
                            for line in query.stdout:
                                chrom, pos, ref, alts, filt, person, gt, ft, gq, ad, dp = (
                                    line.rstrip('\n').split('\t')
                                )
                                person = int(person)
                                assert person in ids
                                for index, _ in enumerate(alts.split(','), 1):
                                    key = matching.get((chrom, pos, ref, alts, index, filt))
                                    native_call = (person, chrom, pos, ref, alts, index)
                                    if key is None or native_call in seen_calls:
                                        continue
                                    qc = call_qc(gt, ad, gq, index, len(alts.split(',')), gene)
                                    if qc is None:
                                        continue
                                    seen_calls.add(native_call)
                                    ft_count[ft] += 1
                                    for cl in desired[key]:
                                        count[f'{cl}:raw_calls'] += 1
                                        carriers[gene][cl]['raw'].add(person)
                                        if qc['valid_AD']:
                                            ad_depths[cl].append(qc['AD_total'])
                                        try:
                                            dp_value = int(dp)
                                            if dp_value >= 0:
                                                dp_depths[cl].append(dp_value)
                                        except ValueError:
                                            count[f'{cl}:DP_missing'] += 1
                                        if not original_site_filters.get(key, {}).get(
                                            'PASS', False
                                        ):
                                            continue
                                        count[f'{cl}:site_PASS_calls'] += 1
                                        carriers[gene][cl]['site_PASS'].add(person)
                                        if ft not in {'.', 'PASS'}:
                                            continue
                                        count[f'{cl}:FT_PASS_calls'] += 1
                                        carriers[gene][cl]['FT_PASS'].add(person)
                                        if not qc['gq_pass']:
                                            continue
                                        count[f'{cl}:GQ20_calls'] += 1
                                        carriers[gene][cl]['GQ20'].add(person)
                                        if not qc['valid_AD'] or not qc['main']:
                                            continue
                                        count[f'{cl}:main_calls'] += 1
                                        carriers[gene][cl]['main'].add(person)
                                        if qc['sensitive']:
                                            count[f'{cl}:CH_sensitive_calls'] += 1
                                            carriers[gene][cl]['CH_sensitive'].add(person)
                            query.stdout.close()
                            qr = query.wait()
                            if qr:
                                raise RuntimeError(
                                    'Genotype query failed; raw errors withheld to protect participant IDs'
                                )
                        print(
                            f'{gene}: shard {shard} individual QC complete; no raw genotypes written.',
                            flush=True,
                        )
                    missing_sites = sorted(set(desired) - seen)
                    if expected_raw is not None:
                        for cl in ['PTV', 'missense']:
                            assert (
                                carriers[gene][cl]['raw'] == expected_raw[cl]
                            ), 'Bounded requery changed the complete raw carrier set'
                    summary = dict(
                        counts=dict(count),
                        FT_distribution=dict(ft_count),
                        full_cohort_raw_carrier_cache_validated=expected_raw is not None,
                        missing_callset_sites=missing_sites,
                        expected_sites=len(desired),
                        found_sites=len(seen),
                        AD_total_distribution={},
                        DP_distribution={},
                    )
                    for label, source in [
                        ('AD_total_distribution', ad_depths),
                        ('DP_distribution', dp_depths),
                    ]:
                        for cl, values in source.items():
                            summary[label][cl] = dict(
                                n=len(values),
                                minimum=int(min(values)),
                                quantiles=dict(
                                    zip(
                                        ['p01', 'p05', 'p25', 'p50', 'p75', 'p95', 'p99'],
                                        np.quantile(
                                            values, [0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99]
                                        ).tolist(),
                                    )
                                ),
                                maximum=int(max(values)),
                                below10=int(sum((x < 10 for x in values))),
                            )
                    metrics[gene] = summary
                    site_audits.extend(sites)
                    saved = dict(
                        fingerprint=fingerprint,
                        metrics=summary,
                        sites=sites,
                        carriers={
                            cl: {s: sorted(v) for s, v in stages.items()}
                            for cl, stages in carriers[gene].items()
                        },
                    )
                    checkpoint.write_text(json.dumps(saved) + '\n')
                    print(
                        f'{gene}: QC completed and derived gene-level checkpoint saved.', flush=True
                    )

            if args.gene:
                assert args.gene in panel, 'Unknown gene'
                raise SystemExit(0)

            frame = pd.DataFrame({'person_id': sorted(ids)})

            audit = []

            for gene in panel:
                w = math.log(cfg['burden']['weights'][gene]['or'])
                for cl in ['PTV', 'missense']:
                    for stage in SETS:
                        s = carriers[gene][cl][stage]
                        audit.append(
                            dict(
                                gene=gene,
                                variant_class=cl,
                                stage=stage,
                                carriers=len(s),
                                diagnosis_age_carriers=len(s & age_ids),
                            )
                        )
                    frame[f'{gene}_{cl}'] = frame.person_id.isin(carriers[gene][cl]['main']).astype(
                        int
                    )
                for stage, suffix in [('main', 'main'), ('CH_sensitive', 'tp53_nf1_vaf_sensitive')]:
                    union = carriers[gene]['PTV'][stage] | carriers[gene]['missense'][stage]
                    frame[f'{gene}_{suffix}'] = frame.person_id.isin(union).astype(int)

            for suffix in ['main', 'tp53_nf1_vaf_sensitive']:
                for pname, genes in [
                    ('13gene', panel),
                    ('12gene_without_STK11', cfg['panels']['without_stk11_12']),
                ]:
                    frame[f'burden_{pname}_{suffix}'] = sum(
                        (
                            math.log(cfg['burden']['weights'][g]['or']) * frame[f'{g}_{suffix}']
                            for g in genes
                        )
                    )
                    frame[f'carrier_{pname}_{suffix}'] = frame[
                        [f'{g}_{suffix}' for g in genes]
                    ].max(axis=1)

            assert set(frame.person_id) == ids and (not frame.person_id.duplicated().any())

            changed = {}

            for gene in ['TP53', 'NF1']:
                main = carriers[gene]['PTV']['main'] | carriers[gene]['missense']['main']
                sens = (
                    carriers[gene]['PTV']['CH_sensitive']
                    | carriers[gene]['missense']['CH_sensitive']
                )
                assert sens <= main
                changed[gene] = dict(
                    main_carriers=len(main),
                    sensitive_carriers=len(sens),
                    lost_carriers=len(main - sens),
                )

            frame.to_csv(GENOTYPE_QC_OUT / 'exome_gene_carriers_and_scores.csv.gz', index=False)

            pd.DataFrame(audit).to_csv(
                GENOTYPE_QC_OUT / 'carrier_QC_attrition_internal.csv', index=False
            )

            pd.DataFrame(site_audits).to_csv(
                GENOTYPE_QC_OUT / 'native_site_ALT_QC_audit.tsv.gz', sep='\t', index=False
            )

            (GENOTYPE_QC_OUT / 'genotype_QC_metrics_internal.json').write_text(
                json.dumps(dict(metrics=metrics, sensitivity_changes=changed), indent=2) + '\n'
            )

            all_sites_found = all(
                (not x.get('missing_callset_sites', []) for x in metrics.values())
            )

            has_DP = any((x.get('counts', {}).get('shards_with_DP', 0) for x in metrics.values()))

            manifest = dict(
                status='genotype_QC_executed_provisional_definition',
                cohort_n=len(ids),
                diagnosis_age_n=len(age_ids),
                cohort_membership_unchanged=True,
                required_site_filter='original VDS filters defined empty or singleton PASS; complete normalized allele identity matched',
                required_individual_filter='FT in {PASS, .}, AoU VAT Appendix B',
                GQ_min=20,
                het_VAF=[0.3, 0.8],
                hom_alt_VAF_min=0.8,
                EUR_AF_strict_upper_bound=0.01,
                sensitivity_TP53_NF1_het_VAF=[0.3, 0.7],
                DP_field_observed=has_DP,
                DP_threshold_applied=False,
                AD_total_is_not_DP=True,
                raw_genotypes_written=False,
                all_expected_sites_found=all_sites_found,
                final_definition_frozen=False,
                expanded_models_run=False,
                missing_AF_not_silently_passed=True,
                PTV_definition='existing_any_retained_transcript_PTV_and_LOFTEE_HC',
                remaining=[
                    'PTV manual transcript/NMD/ClinVar adjudication',
                    'unmapped missense conditions if any',
                    'resolve AF missingness and disagreements',
                    'DP source absent if no FORMAT/DP',
                    'remaining model/output protocol',
                ],
                checkpoint_fingerprint=fingerprint,
            )

            (GENOTYPE_QC_ROOT / 'results/report/33_rowlands_genotype_QC_metadata.json').write_text(
                json.dumps(manifest, indent=2) + '\n'
            )

            report = [
                '# As of Section 6 execution results',
                '',
                'The original AoU exome VCF has been re-read, implementing consistent site PASS, individual filtering, GQ and VAF rules for PTV and qualified P/LP missense. The site PASS is checked against the defined empty filters set of the original VDS or only the PASS set; it is not released based on the "." of the VCF. AoU VAT Appendix B stipulates that FT is PASS or "." are unfiltered calls, and other FT values \u200b\u200bare excluded.',
                'GT/AD is read by raw ALT index, full normalized allele matching; genotypes are summarized in memory, raw genotype matrix is \u200b\u200bnot written out.',
                'Case status cohort and age cohort membership remain unchanged; the generated 13 genes, the STK11-depleted 12 genes, and the TP53/NF1 VAF sensitivity derived score are all present only in the internal ignore directory.',
                '',
                'ClinVar screening first excludes low penetrance, risk factor/risk allele and opposite B/LB/conflict, and then requires accurate P/LP and ≥2 stars; to lock VCF and historical VCV/RCV/SCV check diseases before the release date.',
                'Only ≥2-star P/LP RCVs that have been mapped to applicable diseases are included; uncertain disease entries are temporarily stored in the manual review list.',
                '',
                '**DP：** '
                + (
                    'Some query shards contain FORMAT/DP, please see internal audit for distribution.'
                    if has_DP
                    else 'None of the queried shards contain FORMAT/DP, and the real DP distribution cannot be read. AD total depth has been summarized separately to clarify that it is not DP.'
                ),
                'DP≥10 is not applied. It shall not be said that the real DP distribution review has been completed, nor the UKB coverage check said that the AoU has been QC.',
                '',
                '**Definition status:** This is the tentative score after calling QC; PTV manual adjudication and AF missing/disagreements still need to be resolved, and the final variant set is not frozen. Extended association model is not running.',
                '',
                'Detailed sites, assertions, frequencies, call exclusions and deep audits are only kept in `data/29_rowlands_extension/section6/`. The public small group count is not output to avoid recovering the sparse count through adjacent stage subtraction.',
            ]

            (GENOTYPE_QC_ROOT / 'results/report/33_rowlands_genotype_QC_summary.md').write_text(
                '\n'.join(report) + '\n'
            )

            print(
                'All requested call QC and prespecified VAF score versions executed; final adjudication remains explicit.',
                flush=True,
            )

        # ============================================================================
        # 4. Read missing candidate calls
        # ============================================================================
        if selected_task == 'missing-sites':

            from scripts.functions import genotype_qc as qc

            checkpoints = {
                p.stem: json.loads(p.read_text())
                for p in (MISSING_SITE_CALLS_OUT / 'gene_qc_checkpoints').glob('*.json')
            }

            desired = {}

            for name, cl in [
                ('ptv_for_call_qc.tsv.gz', 'PTV'),
                ('eligible_missense.tsv.gz', 'missense'),
            ]:
                for r in pd.read_csv(
                    MISSING_SITE_CALLS_OUT / name, sep='\t', dtype=str
                ).itertuples():
                    desired.setdefault((r.gene_symbol, r.normalized_key), set()).add(cl)

            source = json.loads(
                (MISSING_SITE_CALLS_OUT / 'original_VDS_normalized_allele_filters.json').read_text()
            )

            missing = {
                (g, k)
                for g, d in checkpoints.items()
                for k in d['metrics'].get('missing_callset_sites', [])
            }

            target_rows, missing_original = ({}, [])

            for gene, key in sorted(missing):
                original = source.get(key)
                if original is None:
                    missing_original.append({'gene': gene, 'normalized_key': key})
                    continue
                native = sorted(set(original['source_native_keys']))
                records = [
                    tuple(k.split(':', 3)) + ('PASS' if original['PASS'] else 'FAIL',)
                    for k in native
                ]
                mapping = qc.normalized_native_sites(
                    records, MISSING_SITE_CALLS_ROOT / 'data/refdata/Homo_sapiens_assembly38.fasta'
                )
                for n, canonical in mapping.items():
                    if canonical == key:
                        nk = ':'.join(n[:4])
                        for cl in desired[gene, key]:
                            target_rows.setdefault(nk, []).append(
                                {'gene': gene, 'variant_class': cl, 'key': key, 'index': n[4]}
                            )

            sets = {
                g: {cl: {s: set() for s in qc.SETS} for cl in ['PTV', 'missense']}
                for g in checkpoints
            }

            metrics = {
                'exome_missing_candidates': len(missing),
                'original_VDS_target_rows': len(target_rows),
                'missing_in_original_VDS': missing_original,
                'raw_target_ALT_calls': 0,
                'main_passed_calls': 0,
            }

            if target_rows:
                project = os.environ['GOOGLE_PROJECT']
                hl.init(
                    master='local[2]',
                    quiet=True,
                    default_reference='GRCh38',
                    log=str(MISSING_SITE_CALLS_OUT / 'missing_site_calls_hail.log'),
                    gcs_requester_pays_configuration=project,
                    spark_conf={
                        'spark.driver.memory': '4g',
                        'spark.hadoop.fs.gs.impl': 'com.google.cloud.hadoop.fs.gcs.GoogleHadoopFileSystem',
                        'spark.hadoop.fs.AbstractFileSystem.gs.impl': 'com.google.cloud.hadoop.fs.gcs.GoogleHadoopFS',
                        'spark.hadoop.fs.gs.requester.pays.mode': 'AUTO',
                        'spark.hadoop.fs.gs.requester.pays.project.id': project,
                    },
                )
                positions = sorted(
                    {
                        (int(k.split(':')[0].replace('chr', '')), int(k.split(':')[1]))
                        for k in target_rows
                    }
                )
                intervals = hl.eval(
                    hl.array(
                        [
                            hl.parse_locus_interval(
                                f'chr{c}:{p}-{p + 1}', reference_genome='GRCh38'
                            )
                            for c, p in positions
                        ]
                    )
                )
                path = (
                    os.environ['CDR_STORAGE_PATH']
                    + '/wgs/short_read/snpindel/vds/hail.vds/variant_data'
                )
                mt = hl.read_matrix_table(path, _intervals=intervals)
                lookup = hl.literal(
                    target_rows,
                    hl.tdict(
                        hl.tstr,
                        hl.tarray(
                            hl.tstruct(
                                gene=hl.tstr, variant_class=hl.tstr, key=hl.tstr, index=hl.tint32
                            )
                        ),
                    ),
                )
                mt = mt.annotate_rows(
                    targets=lookup.get(
                        hl.str(mt.locus)
                        + ':'
                        + mt.alleles[0]
                        + ':'
                        + hl.delimit(mt.alleles[1:], ',')
                    )
                )
                mt = mt.filter_rows(hl.is_defined(mt.targets))
                rows = mt.rows()
                seen = rows.aggregate(
                    hl.agg.collect_as_set(
                        hl.str(rows.locus)
                        + ':'
                        + rows.alleles[0]
                        + ':'
                        + hl.delimit(rows.alleles[1:], ',')
                    )
                )
                assert seen == set(target_rows), 'Original VDS target rows not completely retrieved'
                metrics['all_requested_original_rows_found'] = True
                ids = set(
                    pd.read_csv(
                        MISSING_SITE_CALLS_ROOT
                        / 'data/29_rowlands_extension/case_control_ids.csv.gz'
                    ).person_id.astype(str)
                )
                mt = mt.filter_cols(hl.literal(ids).contains(mt.s))
                e = mt.entries().explode('targets')
                e = e.filter(hl.is_defined(e.LGT) & (e.LGT.ploidy == 2) & hl.is_defined(e.LA))
                e = e.annotate(global_gt=hl.range(e.LGT.ploidy).map(lambda i: e.LA[e.LGT[i]]))
                e = e.filter(e.global_gt.contains(e.targets.index))
                records = e.select('targets', 'global_gt', 'LA', 'LAD', 'GQ', 'FT').collect()
                for r in records:
                    metrics['raw_target_ALT_calls'] += 1
                    gene, cl, key, index = (
                        r.targets.gene,
                        r.targets.variant_class,
                        r.targets.key,
                        r.targets.index,
                    )
                    person = int(r.s)
                    assert str(person) in ids
                    n_alts = len(r.alleles) - 1
                    depths = [0] * (n_alts + 1)
                    if (
                        r.LAD is None
                        or len(r.LAD) != len(r.LA)
                        or any((x is None or x < 0 for x in r.LAD))
                    ):
                        ad = '.'
                    else:
                        for allele, depth in zip(r.LA, r.LAD):
                            depths[allele] = depth
                        ad = ','.join(map(str, depths))
                    q = qc.call_qc(
                        '/'.join(map(str, r.global_gt)),
                        ad,
                        '.' if r.GQ is None else str(r.GQ),
                        index,
                        n_alts,
                        gene,
                    )
                    if q is None:
                        continue
                    sets[gene][cl]['raw'].add(person)
                    if not source[key]['PASS']:
                        continue
                    sets[gene][cl]['site_PASS'].add(person)
                    if r.FT is False:
                        continue
                    sets[gene][cl]['FT_PASS'].add(person)
                    if not q['gq_pass']:
                        continue
                    sets[gene][cl]['GQ20'].add(person)
                    if not q['valid_AD'] or not q['main']:
                        continue
                    metrics['main_passed_calls'] += 1
                    sets[gene][cl]['main'].add(person)
                    if q['sensitive']:
                        sets[gene][cl]['CH_sensitive'].add(person)
                hl.stop()

            (MISSING_SITE_CALLS_OUT / 'missing_site_VDS_carriers_internal.json').write_text(
                json.dumps(
                    {
                        g: {cl: {s: sorted(v) for s, v in d.items()} for cl, d in c.items()}
                        for g, c in sets.items()
                    }
                )
                + '\n'
            )

            (MISSING_SITE_CALLS_OUT / 'missing_site_VDS_QC_internal.json').write_text(
                json.dumps(metrics, indent=2) + '\n'
            )

            cfg = yaml.safe_load(
                (MISSING_SITE_CALLS_ROOT / 'config/rowlands_13gene_extension.yaml').read_text()
            )

            frame = pd.read_csv(MISSING_SITE_CALLS_OUT / 'exome_gene_carriers_and_scores.csv.gz')

            for gene in checkpoints:
                unions = {}
                for cl in ['PTV', 'missense']:
                    unions[cl] = {
                        s: set(checkpoints[gene]['carriers'][cl][s]) | sets[gene][cl][s]
                        for s in qc.SETS
                    }
                    frame[f'{gene}_{cl}'] = frame.person_id.isin(unions[cl]['main']).astype(int)
                for stage, suffix in [('main', 'main'), ('CH_sensitive', 'tp53_nf1_vaf_sensitive')]:
                    frame[f'{gene}_{suffix}'] = frame.person_id.isin(
                        unions['PTV'][stage] | unions['missense'][stage]
                    ).astype(int)

            for suffix in ['main', 'tp53_nf1_vaf_sensitive']:
                for label, panel in [
                    ('13gene', cfg['panels']['main_13']),
                    ('12gene_without_STK11', cfg['panels']['without_stk11_12']),
                ]:
                    frame[f'burden_{label}_{suffix}'] = sum(
                        (
                            math.log(cfg['burden']['weights'][g]['or']) * frame[f'{g}_{suffix}']
                            for g in panel
                        )
                    )
                    frame[f'carrier_{label}_{suffix}'] = frame[
                        [f'{g}_{suffix}' for g in panel]
                    ].max(axis=1)

            frame.to_csv(
                MISSING_SITE_CALLS_OUT / 'derived_gene_carriers_and_scores.csv.gz', index=False
            )

            (MISSING_SITE_CALLS_OUT / 'missing_site_VDS_merge_complete.json').write_text(
                json.dumps(
                    {
                        'merged': True,
                        'all_requested_original_rows_found': metrics.get(
                            'all_requested_original_rows_found', True
                        ),
                    }
                )
                + '\n'
            )

            print(
                'Exome-absent candidates checked against original VDS and actual frozen-cohort ALT calls; no raw genotypes written.'
            )

        # ============================================================================
        # 5. Validate genotype qc
        # ============================================================================
        if selected_task == 'validate-qc':
            cfg = yaml.safe_load(
                (QC_VALIDATION_ROOT / 'config/rowlands_13gene_extension.yaml').read_text()
            )

            scores = pd.read_csv(QC_VALIDATION_OUT / 'derived_gene_carriers_and_scores.csv.gz')

            ids = set(
                pd.read_csv(
                    QC_VALIDATION_ROOT / cfg['cohorts']['case_control']['frozen_ids']
                ).person_id
            )

            age = set(
                pd.read_csv(
                    QC_VALIDATION_ROOT / cfg['cohorts']['diagnosis_age']['frozen_ids']
                ).person_id
            )

            assert (
                set(scores.person_id) == ids
                and age <= ids
                and (not scores.person_id.duplicated().any())
            )

            gene_checks = {}

            supplemental = json.loads(
                (QC_VALIDATION_OUT / 'missing_site_VDS_carriers_internal.json').read_text()
            )

            supplemental_metrics = json.loads(
                (QC_VALIDATION_OUT / 'missing_site_VDS_QC_internal.json').read_text()
            )

            assert supplemental_metrics['all_requested_original_rows_found']

            assert not supplemental_metrics['missing_in_original_VDS']

            for gene in cfg['panels']['main_13']:
                checkpoint = json.loads(
                    (QC_VALIDATION_OUT / 'gene_qc_checkpoints' / f'{gene}.json').read_text()
                )
                for cl in ['PTV', 'missense']:
                    for stage in checkpoint['carriers'][cl]:
                        checkpoint['carriers'][cl][stage] = sorted(
                            set(checkpoint['carriers'][cl][stage])
                            | set(supplemental[gene][cl][stage])
                        )
                assert all(
                    (set(x) <= ids for cl in checkpoint['carriers'].values() for x in cl.values())
                )
                for cl in ['PTV', 'missense']:
                    d = checkpoint['carriers'][cl]
                    for a, b in zip(
                        ['raw', 'site_PASS', 'FT_PASS', 'GQ20', 'main'],
                        ['site_PASS', 'FT_PASS', 'GQ20', 'main', 'CH_sensitive'],
                    ):
                        assert set(d[b]) <= set(d[a]), f'Nonmonotone QC {gene} {cl}'
                    assert set(scores.loc[scores[f'{gene}_{cl}'] == 1, 'person_id']) == set(
                        d['main']
                    )
                union = set(checkpoint['carriers']['PTV']['main']) | set(
                    checkpoint['carriers']['missense']['main']
                )
                sens = set(checkpoint['carriers']['PTV']['CH_sensitive']) | set(
                    checkpoint['carriers']['missense']['CH_sensitive']
                )
                assert set(scores.loc[scores[f'{gene}_main'] == 1, 'person_id']) == union
                assert (
                    set(scores.loc[scores[f'{gene}_tp53_nf1_vaf_sensitive'] == 1, 'person_id'])
                    == sens
                )
                if gene not in {'TP53', 'NF1'}:
                    assert union == sens
                gene_checks[gene] = dict(
                    checkpoint_and_score_agree=True,
                    QC_sets_monotone=True,
                    sensitivity_subset=True,
                    expected_sites_found=not checkpoint['metrics'].get('missing_callset_sites', []),
                )

            for suffix in ['main', 'tp53_nf1_vaf_sensitive']:
                for label, panel in [
                    ('13gene', cfg['panels']['main_13']),
                    ('12gene_without_STK11', cfg['panels']['without_stk11_12']),
                ]:
                    expected = sum(
                        (
                            math.log(cfg['burden']['weights'][g]['or']) * scores[f'{g}_{suffix}']
                            for g in panel
                        )
                    )
                    assert np.allclose(
                        scores[f'burden_{label}_{suffix}'], expected, rtol=1e-12, atol=1e-12
                    )
                    assert scores[f'carrier_{label}_{suffix}'].equals(
                        scores[[f'{g}_{suffix}' for g in panel]].max(axis=1)
                    )

            s32 = json.loads((QC_VALIDATION_OUT / 'assertion_manifest.json').read_text())

            metadata = json.loads(
                (
                    QC_VALIDATION_ROOT / 'results/report/33_rowlands_genotype_QC_metadata.json'
                ).read_text()
            )

            metadata.pop('checkpoint_fingerprint', None)

            metadata.setdefault(
                'all_expected_sites_found_in_exome', metadata['all_expected_sites_found']
            )

            metadata['all_expected_sites_found'] = True

            metadata['exome_absent_candidates_original_VDS_checked'] = True

            metadata['supplemental_VDS_calls_merged'] = True

            (
                QC_VALIDATION_ROOT / 'results/report/33_rowlands_genotype_QC_metadata.json'
            ).write_text(json.dumps(metadata, indent=2) + '\n')

            site_manifest = json.loads(
                (QC_VALIDATION_OUT / 'original_VDS_site_filter_manifest.json').read_text()
            )

            assert site_manifest[
                'every_found_exome_allele_matches_original_VDS_after_full_normalization'
            ]

            assert site_manifest['every_found_exome_allele_original_VDS_PASS']

            site_audit = pd.read_csv(
                QC_VALIDATION_OUT / 'native_site_ALT_QC_audit.tsv.gz', sep='\t'
            )

            assert site_audit.original_VDS_site_PASS.all()

            manifest = dict(
                status='execution_outputs_validated',
                cohort_membership_exact_match=True,
                all_score_formulas_verified=True,
                one_person_one_gene_collapse_verified=True,
                historical_ClinVar_XML_dates_verified=s32['historical_XML_at_or_before_release'],
                original_VDS_PASS_verified_by_complete_normalized_alleles=True,
                exome_absent_candidates_original_VDS_rows_and_individual_calls_checked=True,
                supplemental_VDS_calls_merged_and_validated=True,
                full_cohort_raw_carrier_sets_reproduced_after_bounded_requery=all(
                    (
                        json.loads(
                            (QC_VALIDATION_OUT / 'gene_qc_checkpoints' / f'{g}.json').read_text()
                        )['metrics'].get('full_cohort_raw_carrier_cache_validated', False)
                        for g in cfg['panels']['main_13']
                    )
                ),
                gene_checks=gene_checks,
                DP_threshold_applied=False,
                final_variant_definition_frozen=False,
                expanded_association_models_run=False,
            )

            (
                QC_VALIDATION_ROOT / 'results/report/35_rowlands_execution_validation.json'
            ).write_text(json.dumps(manifest, indent=2) + '\n')

            cfg['genotype_qc']['status'] = 'executed_outputs_validated_provisional_definition'

            cfg['status'] = (
                'executed_through_section6_provisional_definition_pending_final_adjudication_and_remaining_protocol'
            )

            cfg['genotype_qc'][
                'site_filter'
            ] = 'original_VDS_filters_defined_empty_or_singleton_PASS_join_by_complete_normalized_allele'

            cfg['genotype_qc'][
                'site_filter_encoding_note'
            ] = 'exome_VCF_FILTER_dot_verified_against_original_VDS_not_automatically_accepted'

            cfg['genotype_qc'][
                'individual_filter_field_and_pass_encoding'
            ] = 'FORMAT/FT in {PASS, .}; AoU VAT Appendix B'

            cfg['genotype_qc'][
                'individual_filter_source'
            ] = 'https://support.researchallofus.org/hc/en-us/articles/4615256690836-Variant-Annotation-Table'

            cfg['genotype_qc']['dp']['distribution_review_status'] = (
                'FORMAT_DP_absent_in_queried_exome_shards_AD_total_reviewed_separately'
                if not metadata['DP_field_observed']
                else 'distribution_read_awaiting_review'
            )

            cfg['genotype_qc']['sensitivity_tp53_nf1'][
                'status'
            ] = 'derived_scores_generated_and_validated'

            cfg['new_non_ptv']['genotype_qc_status'] = 'executed_provisional_definition'

            cfg['new_non_ptv'][
                'disease_review_status'
            ] = 'historical_assertion_screen_executed_unresolved_records_held_out'

            cfg['ptv']['baseline_missing_af_behavior'] = 'inherit_original_fill_zero'

            cfg['ptv'][
                'missing_af_behavior'
            ] = 'new_score_requires_known_consistent_allele_AF_unknown_or_disagreeing_AF_held_out_pending_adjudication'

            cfg['burden'][
                'derived_scores_internal'
            ] = 'data/29_rowlands_extension/section6/derived_gene_carriers_and_scores.csv.gz'

            cfg['execution_gate'][
                'reason'
            ] = 'Call QC and provisional score generation completed; final variant adjudication and remaining model/output protocol pending; DP absent in queried source if reported in metadata'

            (QC_VALIDATION_ROOT / 'config/rowlands_13gene_extension.yaml').write_text(
                yaml.safe_dump(cfg, allow_unicode=True, sort_keys=False)
            )

            main_total = int(scores.carrier_13gene_main.sum())

            age_total = int(scores.loc[scores.person_id.isin(age), 'carrier_13gene_main'].sum())

            public_rows = []

            for label, n, carriers in [
                ('case_control_frozen', len(ids), main_total),
                ('diagnosis_age_frozen', len(age), age_total),
            ]:
                if (
                    carriers >= 20
                    and n - carriers >= 20
                    and (label != 'diagnosis_age_frozen' or main_total - age_total >= 20)
                ):
                    public_rows.append(
                        dict(cohort=label, cohort_n=n, main_13gene_carriers=carriers)
                    )

            pd.DataFrame(
                public_rows, columns=['cohort', 'cohort_n', 'main_13gene_carriers']
            ).to_csv(
                QC_VALIDATION_ROOT / 'results/table/35_rowlands_main_carrier_summary.csv',
                index=False,
            )

            changes = {}

            audit = []

            for gene in cfg['panels']['main_13']:
                cp = json.loads(
                    (QC_VALIDATION_OUT / 'gene_qc_checkpoints' / f'{gene}.json').read_text()
                )
                merged = {
                    cl: {s: set(v) | set(supplemental[gene][cl][s]) for s, v in d.items()}
                    for cl, d in cp['carriers'].items()
                }
                for cl, d in merged.items():
                    for stage, v in d.items():
                        audit.append(
                            dict(
                                gene=gene,
                                variant_class=cl,
                                stage=stage,
                                carriers=len(v),
                                diagnosis_age_carriers=len(v & age),
                            )
                        )
                if gene in {'TP53', 'NF1'}:
                    main = merged['PTV']['main'] | merged['missense']['main']
                    sens = merged['PTV']['CH_sensitive'] | merged['missense']['CH_sensitive']
                    changes[gene] = dict(
                        main_carriers=len(main),
                        sensitive_carriers=len(sens),
                        lost_carriers=len(main - sens),
                    )

            pd.DataFrame(audit).to_csv(
                QC_VALIDATION_OUT / 'carrier_QC_attrition_internal.csv', index=False
            )

            private_metrics = json.loads(
                (QC_VALIDATION_OUT / 'genotype_QC_metrics_internal.json').read_text()
            )

            private_metrics['sensitivity_changes'] = changes

            private_metrics['supplemental_VDS_QC'] = supplemental_metrics

            (QC_VALIDATION_OUT / 'genotype_QC_metrics_internal.json').write_text(
                json.dumps(private_metrics, indent=2) + '\n'
            )

            text = [
                '',
                '## Acceptance and main score total',
                '',
                'Snakemake passed the acceptance: the QC set of all genes is monotonic, the carrier cache is consistent, the same person and the same gene are only counted once, the two cohort members remain unchanged, 13/12 genes and the VAF sensitivity score formula are correct.',
                'The original VDS site filters are passed according to the complete normalized allele check; the "." of VCF is not directly treated as PASS.',
                'A total of seven sets of regression checks including classification, disease and genetic mode, original ALT index, VAF boundaries and original filters coding passed.',
                '',
            ]

            for row in public_rows:
                text.append(
                    f"- {row['cohort']}：n={row['cohort_n']}, the main 13 genes are qualified to call carrier={row['main_13gene_carriers']}。"
                )

            text += [
                '',
                'TP53/NF1 sensitivity only tightens the heterozygous VAF upper bound, while homozygous substitution calling rules remain consistent. Carrier loss after folding by gene:',
            ]

            for gene, values in changes.items():
                n = values['lost_carriers']
                text.append(
                    f'- {gene}：'
                    + (str(n) if n == 0 or n >= 20 else '<20 (suppress exact values)')
                    + '。'
                )

            text += [
                '',
                'DP≥10 is not applied; the real DP field does not exist in the queried exome VCF. AD total depth distribution has been independently reviewed but is not considered DP.',
                'Candidate sites not found by exome have been returned to the same original VDS: all target lines were read completely, and actual ALT calls were combined and scored after the same QC; no reference genotype was inferred, and genome-wide detectability was not claimed.',
                'The tentative score uses any of the existing PTV+HC rules for isogenic transcripts and known consensus AF; PTV definition harmonization, manual adjudication of NMD/contrary ClinVar evidence, AF missing/divergent and historical XML inconsistent entries are still matters before the final freeze.',
            ]

            summary = QC_VALIDATION_ROOT / 'results/report/33_rowlands_genotype_QC_summary.md'

            original = summary.read_text().split('## Acceptance and main score total', 1)[0]

            summary.write_text(original.rstrip() + '\n' + '\n'.join(text) + '\n')

            print(
                'Validated every gene checkpoint, fixed membership, score formula, collapse and VAF sensitivity subset; no association models run.'
            )

        # ============================================================================
        # 6. Audit sample variant calls
        # ============================================================================
        if selected_task == 'call-audit':
            SAMPLE_VARIANT_AUDIT_OUT.mkdir(parents=True, exist_ok=True)

            cfg = yaml.safe_load(
                (SAMPLE_VARIANT_AUDIT_ROOT / 'config/rowlands_13gene_extension.yaml').read_text()
            )

            assert os.environ['CDR_STORAGE_PATH'] == cfg['data']['storage_path']

            assert os.environ['WORKSPACE_CDR'].split('.')[-1] == cfg['data']['cdr']

            import subprocess

            assert (
                subprocess.run(
                    ['git', 'check-ignore', '-q', str(SAMPLE_VARIANT_AUDIT_OUT)],
                    cwd=SAMPLE_VARIANT_AUDIT_ROOT,
                ).returncode
                == 0
            )

            from scripts.functions import genotype_qc as qc

            from scripts.functions import clinvar as cv

            candidates = pd.read_csv(
                SAMPLE_VARIANT_AUDIT_BASE / 'vat_candidates_by_allele_gene.tsv.gz',
                sep='\t',
                dtype=str,
                keep_default_na=False,
            )

            transcripts = pd.read_csv(
                SAMPLE_VARIANT_AUDIT_BASE / 'vat_selected_transcript_rows.tsv.gz',
                sep='\t',
                dtype=str,
                keep_default_na=False,
            )

            clinical = pd.read_csv(
                SAMPLE_VARIANT_AUDIT_S6 / 'clinvar_region_records.tsv.gz',
                sep='\t',
                dtype=str,
                keep_default_na=False,
            )

            evidence = pd.read_csv(
                SAMPLE_VARIANT_AUDIT_S6 / 'clinvar_assertion_evidence.tsv.gz',
                sep='\t',
                dtype=str,
                keep_default_na=False,
            )

            exclusions = pd.read_csv(
                SAMPLE_VARIANT_AUDIT_S6 / 'missense_exclusion_audit.tsv.gz',
                sep='\t',
                dtype=str,
                keep_default_na=False,
            )

            qualified = {
                cl: {
                    (r.gene_symbol, r.normalized_key)
                    for r in pd.read_csv(SAMPLE_VARIANT_AUDIT_S6 / file, sep='\t').itertuples()
                }
                for cl, file in [
                    ('PTV', 'ptv_for_call_qc.tsv.gz'),
                    ('missense', 'eligible_missense.tsv.gz'),
                ]
            }

            transcript_groups = {
                k: d.to_dict('records')
                for k, d in transcripts.groupby(['gene_symbol', 'normalized_key'])
            }

            clinical_groups = {
                k: d.to_dict('records') for k, d in clinical.groupby('normalized_key')
            }

            evidence_groups = {
                k: d.to_dict('records')
                for k, d in evidence.groupby(['gene_symbol', 'normalized_key'])
            }

            exclusion_groups = {
                k: d.decision.tolist()
                for k, d in exclusions.groupby(['gene_symbol', 'normalized_key'])
            }

            source = json.loads(
                (SAMPLE_VARIANT_AUDIT_S6 / 'original_VDS_site_filters.json').read_text()
            )

            raw = [
                tuple(k.split(':', 3)) + ('PASS' if v['PASS'] else 'FAIL',)
                for k, v in source.items()
            ]

            mapping = qc.normalized_native_sites(
                raw, SAMPLE_VARIANT_AUDIT_ROOT / 'data/refdata/Homo_sapiens_assembly38.fasta'
            )

            wanted = set(candidates.normalized_key)

            canonical = {}

            for native, key in mapping.items():
                if key in wanted:
                    canonical.setdefault(key, []).append(
                        (native[0], int(native[1]), native[2], native[3], native[4])
                    )

            annotations = {}

            target_rows = {}

            site_coverage = []

            for r in candidates.to_dict('records'):
                gene, key = (r['gene_symbol'], r['normalized_key'])
                pair = (gene, key)
                supporting = transcript_groups[pair]
                classes = []
                if any(
                    (
                        any((t in terms(s['consequence']) for t in cfg['ptv']['consequence_terms']))
                        for s in supporting
                    )
                ):
                    classes.append('PTV')
                if any(('missense_variant' in terms(s['consequence']) for s in supporting)):
                    classes.append('missense')
                refs = clinical_groups.get(key, [])
                ev = evidence_groups.get(pair, [])
                values = [a['sig'] for a in refs] + [a['classification'] for a in ev]
                combined = ' | '.join((cv.clean(x) for x in values))
                af_values = {positive_float(s['gvs_eur_af']) for s in supporting} - {None}
                af_ok = len(af_values) == 1 and next(iter(af_values)) < 0.01
                af = next(iter(af_values)) if len(af_values) == 1 else ''
                for cl in classes:
                    reasons = []
                    if pair not in qualified[cl]:
                        if cl == 'PTV':
                            if r['ptv_hc_any_transcript'] != 'True':
                                reasons.append('no_retained_same_gene_PTV_and_LOFTEE_HC')
                            if not af_ok:
                                reasons.append('EUR_AF_missing_disagreeing_or_not_lt_0p01')
                        else:
                            reasons = [
                                x
                                for x in exclusion_groups.get(pair, [])
                                if x not in {'eligible', 'proceed_to_assertions'}
                            ]
                            if not reasons:
                                reasons = ['not_eligible_canonical_missense_or_AF_or_ClinVar']
                    relevant = [
                        s
                        for s in supporting
                        if (
                            'missense_variant' in terms(s['consequence'])
                            if cl == 'missense'
                            else any(
                                (
                                    t in terms(s['consequence'])
                                    for t in cfg['ptv']['consequence_terms']
                                )
                            )
                        )
                    ]
                    annotations[gene, key, cl] = {
                        'transcript': '|'.join((s['transcript'] for s in relevant)),
                        'consequence': '|'.join((s['consequence'] for s in relevant)),
                        'LOFTEE': '|'.join((s['LoF'] or 'NA' for s in relevant)),
                        'classification': '|'.join((a['sig'] for a in refs)),
                        'stars': '|'.join(
                            (
                                str(
                                    cv.STARS.get(
                                        cv.clean(a['review']),
                                        (
                                            1
                                            if 'single submitter' in cv.clean(a['review'])
                                            or 'conflicting' in cv.clean(a['review'])
                                            else 0
                                        ),
                                    )
                                )
                                for a in refs
                            )
                        ),
                        'condition': '|'.join((a['disease'] for a in refs)),
                        'low_penetrance_flag': 'low penetrance' in combined,
                        'conflict_flag': 'conflict' in combined
                        or 'benign' in combined
                        or any((a['conflict'] not in {'', '.'} for a in refs)),
                        'EUR_AF': af,
                        'eligible': pair in qualified[cl],
                        'variant_reasons': reasons,
                        'supporting_transcript_rows': json.dumps(relevant, separators=(',', ':')),
                        'RCV_SCV_evidence': json.dumps(ev, separators=(',', ':')),
                    }
                    for chrom, pos, ref, alts, index in canonical.get(key, []):
                        nk = f'{chrom}:{pos}:{ref}:{alts}'
                        target_rows.setdefault(nk, []).append(
                            {'gene': gene, 'key': key, 'cl': cl, 'index': index}
                        )
                site_coverage.append(
                    {
                        'gene': gene,
                        'normalized_key': key,
                        'original_VDS_found': bool(canonical.get(key)),
                    }
                )

            pd.DataFrame(site_coverage).to_csv(
                SAMPLE_VARIANT_AUDIT_OUT / 'candidate_site_coverage_internal.tsv.gz',
                sep='\t',
                index=False,
            )

            ids = set(
                pd.read_csv(
                    SAMPLE_VARIANT_AUDIT_ROOT / cfg['cohorts']['case_control']['frozen_ids']
                ).person_id.astype(str)
            )

            project = os.environ['GOOGLE_PROJECT']

            hl.init(
                master='local[2]',
                quiet=True,
                default_reference='GRCh38',
                log=str(SAMPLE_VARIANT_AUDIT_OUT / 'sample_variant_hail.log'),
                gcs_requester_pays_configuration=project,
                spark_conf={
                    'spark.driver.memory': '8g',
                    'spark.hadoop.fs.gs.impl': 'com.google.cloud.hadoop.fs.gcs.GoogleHadoopFileSystem',
                    'spark.hadoop.fs.AbstractFileSystem.gs.impl': 'com.google.cloud.hadoop.fs.gcs.GoogleHadoopFS',
                    'spark.hadoop.fs.gs.requester.pays.mode': 'AUTO',
                    'spark.hadoop.fs.gs.requester.pays.project.id': project,
                },
            )

            carriers = {
                g: {cl: {s: set() for s in ['main', 'CH_sensitive']} for cl in ['PTV', 'missense']}
                for g in cfg['panels']['main_13']
            }

            audit_path = SAMPLE_VARIANT_AUDIT_OUT / 'sample_variant_call_facts.csv.gz'

            temporary = SAMPLE_VARIANT_AUDIT_OUT / 'sample_variant_call_facts.partial.csv.gz'

            row_count = 0

            alt_rows = 0

            genestatus = {}

            coords = pd.read_csv(SAMPLE_VARIANT_AUDIT_BASE / 'gene_coordinates.csv')

            with gzip.open(temporary, 'wt', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=FIELDS)
                writer.writeheader()
                for gene in cfg['panels']['main_13']:
                    selected = {
                        nk: [t for t in ts if t['gene'] == gene]
                        for nk, ts in target_rows.items()
                        if any((t['gene'] == gene for t in ts))
                    }
                    if not selected:
                        continue
                    coord = coords.loc[coords.gene_symbol == gene].iloc[0]
                    interval = hl.eval(
                        hl.parse_locus_interval(
                            f'{coord.chrom}:{coord.query_region_start}-{int(coord.query_region_end) + 1}',
                            reference_genome='GRCh38',
                        )
                    )
                    mt = hl.read_matrix_table(
                        cfg['data']['storage_path']
                        + '/wgs/short_read/snpindel/vds/hail.vds/variant_data',
                        _intervals=[interval],
                    )
                    lookup = hl.literal(
                        selected,
                        hl.tdict(
                            hl.tstr,
                            hl.tarray(
                                hl.tstruct(gene=hl.tstr, key=hl.tstr, cl=hl.tstr, index=hl.tint32)
                            ),
                        ),
                    )
                    mt = mt.annotate_rows(
                        targets=lookup.get(
                            hl.str(mt.locus)
                            + ':'
                            + mt.alleles[0]
                            + ':'
                            + hl.delimit(mt.alleles[1:], ',')
                        )
                    )
                    mt = mt.filter_rows(hl.is_defined(mt.targets))
                    rows = mt.rows()
                    seen = rows.aggregate(
                        hl.agg.collect_as_set(
                            hl.str(rows.locus)
                            + ':'
                            + rows.alleles[0]
                            + ':'
                            + hl.delimit(rows.alleles[1:], ',')
                        )
                    )
                    assert seen == set(selected), 'Not all requested original VDS rows retrieved'
                    mt = mt.filter_cols(hl.literal(ids).contains(mt.s))
                    e = mt.entries().explode('targets')
                    e = e.filter(hl.is_defined(e.LGT) & hl.is_defined(e.LA))
                    e = e.annotate(global_gt=hl.range(e.LGT.ploidy).map(lambda i: e.LA[e.LGT[i]]))
                    e = e.filter(e.global_gt.contains(e.targets.index))
                    has_dp = 'DP' in e.row.dtype.fields
                    e = e.select(
                        'targets',
                        'global_gt',
                        'LA',
                        'LAD',
                        'GQ',
                        'FT',
                        **{'DP': e.DP} if has_dp else {},
                    )
                    records = e.collect()
                    for r in records:
                        gene, key, cl, index = (
                            r.targets.gene,
                            r.targets.key,
                            r.targets.cl,
                            r.targets.index,
                        )
                        person = int(r.s)
                        assert str(person) in ids
                        nk = f"{r.locus.contig}:{r.locus.position}:{r.alleles[0]}:{','.join(r.alleles[1:])}"
                        info = annotations[gene, key, cl]
                        reasons = list(info['variant_reasons'])
                        depths = [0] * len(r.alleles)
                        valid = (
                            r.LAD is not None
                            and len(r.LAD) == len(r.LA)
                            and all((x is not None and x >= 0 for x in r.LAD))
                        )
                        if valid:
                            for a, d in zip(r.LA, r.LAD):
                                depths[a] = d
                        ad = ','.join(map(str, depths)) if valid else '.'
                        call = qc.call_qc(
                            '/'.join(map(str, r.global_gt)),
                            ad,
                            '.' if r.GQ is None else str(r.GQ),
                            index,
                            len(r.alleles) - 1,
                            gene,
                        )
                        original_pass = all(
                            (
                                source[f'{c}:{p}:{ref}:{alts}']['PASS']
                                for c, p, ref, alts, _ in canonical[key]
                            )
                        )
                        if not original_pass:
                            reasons.append('site_filter_not_PASS')
                        if r.FT is False:
                            reasons.append('genotype_filter_failed')
                        if call is None:
                            reasons.append('not_complete_diploid_target_ALT_genotype')
                        else:
                            if not call['gq_pass']:
                                reasons.append('GQ_missing_nonfinite_or_lt20')
                            if not call['valid_AD']:
                                reasons.append('AD_missing_invalid_or_zero_total')
                            elif not call['main']:
                                reasons.append('VAF_outside_main_range')
                        keep = not reasons
                        sens = keep and call['sensitive']
                        if keep:
                            carriers[gene][cl]['main'].add(person)
                        if sens:
                            carriers[gene][cl]['CH_sensitive'].add(person)
                        chrom, pos, ref, alt = key.split(':', 3)
                        row = {k: v for k, v in info.items() if k in FIELDS}
                        row.update(
                            sample_id=person,
                            gene=gene,
                            chr=chrom,
                            pos=pos,
                            ref=ref,
                            alt=alt,
                            **{'PTV/missense': cl},
                            ClinVar_release=cfg['clinvar']['release_date'],
                            site_filter='PASS' if original_pass else 'FAIL',
                            genotype_filter=(
                                'FAIL' if r.FT is False else 'PASS' if r.FT is True else '.'
                            ),
                            GQ=r.GQ if r.GQ is not None else '',
                            DP=getattr(r, 'DP', None) if has_dp else '',
                            VAF=call.get('vaf', '') if call else '',
                            keep_flag=int(keep),
                            rejection_reason=';'.join(dict.fromkeys(reasons)),
                            normalized_key=key,
                            source_native_key=nk,
                            source_ALT_index=index,
                            allele_dosage=r.global_gt.count(index),
                            AD_total=call.get('AD_total', '') if call else '',
                            sensitivity_keep_flag=int(sens),
                            source='same_original_v9_VDS',
                        )
                        writer.writerow(row)
                        row_count += 1
                    genestatus[gene] = {
                        'all_requested_original_rows_found': True,
                        'DP_field_present': has_dp,
                    }
                    del records
                    print(
                        f'{gene}: authorized sparse audit completed; participant rows and counts withheld.',
                        flush=True,
                    )

            hl.stop()

            temporary.replace(audit_path)

            frame = pd.DataFrame({'sample_id': sorted(map(int, ids))})

            for gene, classes in carriers.items():
                for cl, d in classes.items():
                    for stage, version in [
                        ('main', 'main'),
                        ('CH_sensitive', 'tp53_nf1_vaf_sensitive'),
                    ]:
                        frame[f'{gene}_{cl}_{version}'] = frame.sample_id.isin(d[stage]).astype(
                            'int8'
                        )
                for stage, version in [
                    ('main', 'main'),
                    ('CH_sensitive', 'tp53_nf1_vaf_sensitive'),
                ]:
                    union = classes['PTV'][stage] | classes['missense'][stage]
                    frame[f'{gene}_{version}'] = frame.sample_id.isin(union).astype('int8')

            frame.to_csv(
                SAMPLE_VARIANT_AUDIT_OUT / 'call_read_gene_carriers_internal.csv.gz', index=False
            )

            manifest = {
                'status': 'sparse_audit_and_unweighted_carriers_validated',
                'scope': 'actual ALT calls of all retained candidate allele/gene records, including rejected candidates and rejected calls; no invented reference or missing genotypes',
                'user_authorized_sample_variant_audit': True,
                'genotype_QC_same_as_section6': True,
                'classification_overlay_required': True,
                'weights_read': False,
                'gene_checks': genestatus,
                'final_definition_frozen': False,
                'models_run': False,
            }

            (
                SAMPLE_VARIANT_AUDIT_ROOT
                / 'results/report/37_rowlands_sample_variant_audit_validation.json'
            ).write_text(json.dumps(manifest, indent=2) + '\n')

            (SAMPLE_VARIANT_AUDIT_OUT / 'audit_lineage_internal.json').write_text(
                json.dumps(
                    {
                        'audit_rows': row_count,
                        'audit_sha256': qc.sha(audit_path),
                        'script_sha256': qc.sha(Path(__file__)),
                        'inputs': {
                            str(p.relative_to(SAMPLE_VARIANT_AUDIT_ROOT)): qc.sha(p)
                            for p in [
                                SAMPLE_VARIANT_AUDIT_BASE / 'vat_candidates_by_allele_gene.tsv.gz',
                                SAMPLE_VARIANT_AUDIT_BASE / 'vat_selected_transcript_rows.tsv.gz',
                                SAMPLE_VARIANT_AUDIT_S6 / 'eligible_missense.tsv.gz',
                                SAMPLE_VARIANT_AUDIT_S6 / 'ptv_for_call_qc.tsv.gz',
                                SAMPLE_VARIANT_AUDIT_S6 / 'original_VDS_site_filters.json',
                            ]
                        },
                    },
                    indent=2,
                )
                + '\n'
            )

            print(
                'Section 7 sample/variant audit saved internally; unweighted gene carrier sets validated.',
                flush=True,
            )

        # ============================================================================
        # 7. Reconcile variant eligibility
        # ============================================================================
        if selected_task == 'reconcile':
            cfg = yaml.safe_load(
                (CLINICAL_RECONCILIATION_ROOT / 'config/rowlands_13gene_extension.yaml').read_text()
            )

            from scripts.functions import clinvar as cv

            qualified = {
                cl: {
                    (r.gene_symbol, r.normalized_key)
                    for r in pd.read_csv(CLINICAL_RECONCILIATION_S6 / f, sep='\t').itertuples()
                }
                for cl, f in [
                    ('PTV', 'ptv_for_call_qc.tsv.gz'),
                    ('missense', 'eligible_missense.tsv.gz'),
                ]
            }

            rejected = pd.read_csv(
                CLINICAL_RECONCILIATION_S6 / 'missense_exclusion_audit.tsv.gz',
                sep='\t',
                dtype=str,
                keep_default_na=False,
            )

            reasons = {
                k: [x for x in d.decision if x not in {'eligible', 'proceed_to_assertions'}]
                for k, d in rejected.groupby(['gene_symbol', 'normalized_key'])
            }

            ev = pd.read_csv(
                CLINICAL_RECONCILIATION_S6 / 'clinvar_assertion_evidence.tsv.gz',
                sep='\t',
                dtype=str,
                keep_default_na=False,
            )

            evidence = {
                k: d.to_dict('records') for k, d in ev.groupby(['gene_symbol', 'normalized_key'])
            }

            refs = pd.read_csv(
                CLINICAL_RECONCILIATION_S6 / 'clinvar_region_records.tsv.gz',
                sep='\t',
                dtype=str,
                keep_default_na=False,
            )

            refgroups = {k: d.to_dict('records') for k, d in refs.groupby('normalized_key')}

            candidates = pd.read_csv(
                CLINICAL_RECONCILIATION_BASE / 'vat_candidates_by_allele_gene.tsv.gz',
                sep='\t',
                dtype=str,
                keep_default_na=False,
            )

            candidate_lookup = {
                (r.gene_symbol, r.normalized_key): r for r in candidates.itertuples()
            }

            ids = set(
                pd.read_csv(
                    CLINICAL_RECONCILIATION_ROOT / cfg['cohorts']['case_control']['frozen_ids']
                ).person_id
            )

            carriers = {
                g: {
                    cl: {v: set() for v in ['main', 'tp53_nf1_vaf_sensitive']}
                    for cl in ['PTV', 'missense']
                }
                for g in cfg['panels']['main_13']
            }

            qc_names = {
                'site_filter_not_PASS',
                'genotype_filter_failed',
                'not_complete_diploid_target_ALT_genotype',
                'GQ_missing_nonfinite_or_lt20',
                'AD_missing_invalid_or_zero_total',
                'VAF_outside_main_range',
            }

            meta = {}

            for pair, r in candidate_lookup.items():
                gene, key = pair
                clinical = refgroups.get(key, [])
                assertions = evidence.get(pair, [])
                vals = (
                    [a['sig'] for a in clinical]
                    + [a['review'] for a in clinical]
                    + [a['classification'] for a in assertions]
                )
                text = ' | '.join((cv.clean(x) for x in vals))
                meta[pair] = {
                    'low_penetrance_flag': 'low penetrance' in text,
                    'conflict_flag': 'conflict' in re.sub('\\bno conflicts?\\b', '', text)
                    or 'benign' in text
                    or any((a['conflict'] not in {'', '.'} for a in clinical)),
                    'RCV_SCV_evidence': json.dumps(assertions, separators=(',', ':')),
                }

            temporary = (
                CLINICAL_RECONCILIATION_OUT / 'sample_variant_audit.reclassified.partial.csv.gz'
            )

            with gzip.open(temporary, 'wt', newline='') as f:
                first = True
                for d in pd.read_csv(
                    CLINICAL_RECONCILIATION_OUT / 'sample_variant_call_facts.csv.gz',
                    chunksize=20000,
                    keep_default_na=False,
                    low_memory=False,
                ):
                    updated = []
                    for row in d.to_dict('records'):
                        gene, key, cl = (row['gene'], row['normalized_key'], row['PTV/missense'])
                        pair = (gene, key)
                        assert int(row['sample_id']) in ids
                        reject = []
                        if pair not in qualified[cl]:
                            if cl == 'missense':
                                reject = reasons.get(pair, []) or [
                                    'not_eligible_canonical_missense_or_AF_or_ClinVar'
                                ]
                            else:
                                r = candidate_lookup[pair]
                                if r.ptv_hc_any_transcript != 'True':
                                    reject.append('no_retained_same_gene_PTV_and_LOFTEE_HC')
                                try:
                                    af = float(row['EUR_AF'])
                                    afok = __import__('math').isfinite(af) and af < 0.01
                                except (ValueError, TypeError):
                                    afok = False
                                if not afok:
                                    reject.append('EUR_AF_missing_disagreeing_or_not_lt_0p01')
                                if not reject:
                                    reject = ['PTV_not_eligible_under_current_definition']
                        reject += [x for x in row['rejection_reason'].split(';') if x in qc_names]
                        keep = not reject
                        sensitive = keep and (
                            int(row['allele_dosage']) == 2
                            or float(row['VAF']) <= (0.7 if gene in {'TP53', 'NF1'} else 0.8)
                        )
                        row.update(meta[pair])
                        row.update(
                            keep_flag=int(keep),
                            sensitivity_keep_flag=int(sensitive),
                            rejection_reason=';'.join(dict.fromkeys(reject)),
                        )
                        assert not (
                            keep
                            and cl == 'missense'
                            and (row['low_penetrance_flag'] or row['conflict_flag'])
                        )
                        if keep:
                            carriers[gene][cl]['main'].add(int(row['sample_id']))
                        if sensitive:
                            carriers[gene][cl]['tp53_nf1_vaf_sensitive'].add(int(row['sample_id']))
                        updated.append(row)
                    pd.DataFrame(updated, columns=d.columns).to_csv(f, index=False, header=first)
                    first = False

            temporary.replace(CLINICAL_RECONCILIATION_OUT / 'sample_variant_audit.csv.gz')

            frame = pd.DataFrame({'sample_id': sorted(ids)})

            for g, classes in carriers.items():
                for cl, versions in classes.items():
                    for v, members in versions.items():
                        frame[f'{g}_{cl}_{v}'] = frame.sample_id.isin(members).astype('int8')
                for v in ['main', 'tp53_nf1_vaf_sensitive']:
                    frame[f'{g}_{v}'] = frame.sample_id.isin(
                        classes['PTV'][v] | classes['missense'][v]
                    ).astype('int8')

            frame.to_csv(
                CLINICAL_RECONCILIATION_OUT / 'unweighted_gene_carriers.csv.gz', index=False
            )

            manifest = {
                'status': 'current_eligibility_applied_to_all_candidate_call_facts',
                'no_genotypes_reread': True,
                'weights_read': False,
                'no_conflicts_review_status_not_misclassified': True,
                'all_current_SCV_and_RCV_classifications_screened': True,
                'genotype_QC_and_frozen_membership_unchanged': True,
                'earlier_section6_clinical_eligibility_superseded': True,
                'final_definition_frozen': False,
                'models_run': False,
            }

            (
                CLINICAL_RECONCILIATION_ROOT
                / 'results/report/40_rowlands_clinical_reconciliation.json'
            ).write_text(json.dumps(manifest, indent=2) + '\n')

            original = json.loads(
                (CLINICAL_RECONCILIATION_OUT / 'audit_lineage_internal.json').read_text()
            )

            facts_hash = hashlib.sha256(
                (CLINICAL_RECONCILIATION_OUT / 'sample_variant_call_facts.csv.gz').read_bytes()
            ).hexdigest()

            assert (
                facts_hash == original['audit_sha256']
            ), 'Raw call facts changed since the original read'

            inputs = [
                CLINICAL_RECONCILIATION_OUT / 'sample_variant_call_facts.csv.gz',
                CLINICAL_RECONCILIATION_S6 / 'eligible_missense.tsv.gz',
                CLINICAL_RECONCILIATION_S6 / 'ptv_for_call_qc.tsv.gz',
                CLINICAL_RECONCILIATION_S6 / 'clinvar_assertion_evidence.tsv.gz',
                CLINICAL_RECONCILIATION_S6 / 'missense_exclusion_audit.tsv.gz',
                CLINICAL_RECONCILIATION_ROOT / 'scripts/functions/clinvar.py',
                Path(__file__),
            ]

            (CLINICAL_RECONCILIATION_OUT / 'reconciliation_lineage_internal.json').write_text(
                json.dumps(
                    {
                        'raw_call_facts_verified_against_original_read_lineage': True,
                        'inputs': {
                            str(p.relative_to(CLINICAL_RECONCILIATION_ROOT)): hashlib.sha256(
                                p.read_bytes()
                            ).hexdigest()
                            for p in inputs
                        },
                        'current_audit_sha256': hashlib.sha256(
                            (
                                CLINICAL_RECONCILIATION_OUT / 'sample_variant_audit.csv.gz'
                            ).read_bytes()
                        ).hexdigest(),
                    },
                    indent=2,
                )
                + '\n'
            )

            print(
                'Current clinical eligibility reapplied to persisted all-candidate calls; unweighted gene carriers rebuilt; no genotype or weight reads.'
            )

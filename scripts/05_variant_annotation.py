"""Audit PTV/LOFTEE annotation and select qualifying ClinVar germline P/LP missense.

PTV consequences and HC LOFTEE requirements are retained.
Main missense eligibility requires >=2 review stars after the existing exclusions and disease review.
"""

from __future__ import annotations

from pathlib import Path
import sys

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.functions.paths import PROJECT_ROOT
from scripts.functions.cli import parse_tasks

TASKS = ('ptv', 'missense')
if __name__ == '__main__':
    selected_tasks = parse_tasks(__doc__, TASKS, ['ptv', 'missense'])

import gzip

import hashlib

import json

import re

import subprocess


import pandas as pd

import yaml

PTV_ANNOTATION_ROOT = PROJECT_ROOT

PTV_ANNOTATION_OUT = PTV_ANNOTATION_ROOT / 'data/29_rowlands_extension/ptv_audit'

PTV = {'stop_gained', 'frameshift_variant', 'splice_donor_variant', 'splice_acceptor_variant'}

PTV_ANNOTATION_STARS = {
    'criteria_provided,_multiple_submitters,_no_conflicts': 2,
    'reviewed_by_expert_panel': 3,
    'practice_guideline': 4,
}

import importlib.util

import time

from urllib.parse import unquote

import xml.etree.ElementTree as ET

import requests

CLINVAR_ROOT = PROJECT_ROOT

CLINVAR_OUT = CLINVAR_ROOT / 'data/29_rowlands_extension/section6'

XML = CLINVAR_OUT / 'clinvar_xml'

PLP = {'pathogenic', 'likely pathogenic', 'pathogenic/likely pathogenic'}

CLINVAR_STARS = {
    'criteria provided, multiple submitters, no conflicts': 2,
    'criteria provided, single submitter': 1,
    'reviewed by expert panel': 3,
    'practice guideline': 4,
}

from scripts.functions.io import sha256 as sha

from scripts.functions.clinvar import write_tsv_if_changed

from scripts.functions.clinvar import clean

from scripts.functions.clinvar import exclusions

from scripts.functions.clinvar import fetch

from scripts.functions.clinvar import historical_xml

from scripts.functions.clinvar import disease_rule

from scripts.functions.clinvar import review_variant

if __name__ == "__main__":
    for selected_task in selected_tasks:

        # ============================================================================
        # 1. Audit ptv annotation
        # ============================================================================
        if selected_task == 'ptv':
            cfg = yaml.safe_load(
                (PTV_ANNOTATION_ROOT / 'config/rowlands_13gene_extension.yaml').read_text()
            )

            PTV_ANNOTATION_OUT.mkdir(parents=True, exist_ok=True)

            assert (
                subprocess.run(
                    ['git', 'check-ignore', '-q', str(PTV_ANNOTATION_OUT)], cwd=PTV_ANNOTATION_ROOT
                ).returncode
                == 0
            )

            source = (
                PTV_ANNOTATION_ROOT
                / 'data/29_rowlands_extension/vat_selected_transcript_rows.tsv.gz'
            )

            f = pd.read_csv(source, sep='\t', dtype=str, keep_default_na=False)

            f['ptv_hc'] = f.consequence.map(lambda x: bool(set(x.split(',')) & PTV)) & f.LoF.eq(
                'HC'
            )

            tx_ids = set(f.transcript.str.split('.').str[0])

            gtf = PTV_ANNOTATION_ROOT / 'data/29_ptv_audit/Homo_sapiens.GRCh38.110.gtf.gz'

            exons = {}

            with gzip.open(gtf, 'rt') as handle:
                for line in handle:
                    if line.startswith('#'):
                        continue
                    v = line.rstrip().split('\t')
                    if v[2] != 'exon':
                        continue
                    attrs = dict(re.findall('(\\w+) "([^";]*)"', v[8]))
                    tx = attrs.get('transcript_id')
                    if tx in tx_ids:
                        exons.setdefault(tx, []).append((int(v[3]), int(v[4]), v[6]))

            def position_flags(row):
                parts = row.exon_number.split('/')
                terminal = (
                    len(parts) == 2
                    and parts[0].isdigit()
                    and parts[1].isdigit()
                    and (parts[0] == parts[1])
                )
                tx = exons.get(row.transcript.split('.')[0], [])
                flags = []
                if terminal:
                    flags.append('terminal_exon_VAT')
                if not tx:
                    flags.append('GTF_transcript_unresolved')
                else:
                    ordered = sorted(set(tx), reverse=tx[0][2] == '-')
                    pos = int(row.position)
                    if len(ordered) == 1:
                        flags.append('single_exon_possible_NMD_escape')
                    elif len(ordered) > 1:
                        final, penultimate = (ordered[-1], ordered[-2])
                        if final[0] <= pos <= final[1]:
                            flags.append('terminal_exon_GTF_possible_NMD_escape')
                        if penultimate[0] <= pos <= penultimate[1]:
                            distance = (
                                penultimate[1] - pos
                                if penultimate[2] == '+'
                                else pos - penultimate[0]
                            )
                            if distance <= 55:
                                flags.append('variant_within_55nt_of_final_junction_position_proxy')
                cons = set(row.consequence.split(','))
                if cons & {'frameshift_variant', 'splice_donor_variant', 'splice_acceptor_variant'}:
                    flags.append('NMD_requires_altered_transcript_or_new_stop_review')
                if row.LoF_flags or row.LoF_info:
                    flags.append('LOFTEE_flags_or_info_present_review')
                return ';'.join(flags)

            support = f.loc[f.ptv_hc].copy()

            support['position_review_flags'] = [position_flags(r) for r in support.itertuples()]

            support.to_csv(
                PTV_ANNOTATION_OUT / 'ptv_hc_transcript_position_review.tsv.gz',
                sep='\t',
                index=False,
            )

            cv_path = PTV_ANNOTATION_ROOT / cfg['clinvar']['normalized_vcf']

            raw_cv = PTV_ANNOTATION_ROOT / cfg['clinvar']['full_vcf']

            assert (
                sha(raw_cv) == cfg['clinvar']['sha256']
            ), 'ClinVar source snapshot checksum mismatch'

            with gzip.open(cv_path, 'rt') as handle:
                header_date = None
                for line in handle:
                    if line.startswith('##fileDate='):
                        header_date = line.strip().split('=', 1)[1]
                    if not line.startswith('#'):
                        break

            assert (
                header_date == cfg['clinvar']['release_date']
            ), 'Normalized ClinVar release mismatch'

            clinvar = {}

            coords = pd.read_csv(
                PTV_ANNOTATION_ROOT / 'data/29_rowlands_extension/gene_coordinates.csv'
            )

            fmt = '%CHROM\t%POS\t%REF\t%ALT\t%ID\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/CLNSIGCONF\t%INFO/CLNDN\t%INFO/ORIGIN\n'

            for r in coords.itertuples():
                raw = subprocess.run(
                    [
                        'bcftools',
                        'query',
                        '-r',
                        f'{r.chrom}:{r.query_region_start}-{r.query_region_end}',
                        '-f',
                        fmt,
                        str(cv_path),
                    ],
                    capture_output=True,
                    text=True,
                    check=True,
                ).stdout
                for line in raw.splitlines():
                    chrom, pos, ref, alt, vid, sig, rev, conflict, disease, origin = line.split(
                        '\t'
                    )
                    key = f'{chrom}:{pos}:{ref}:{alt}'
                    record = dict(
                        variation_id=vid,
                        significance=sig,
                        review_status=rev,
                        stars=PTV_ANNOTATION_STARS.get(rev, 0),
                        conflict=conflict,
                        disease=disease,
                        origin=origin,
                    )
                    if record not in clinvar.setdefault(key, []):
                        clinvar[key].append(record)

            rows = []

            for (gene, key), group in f.groupby(['gene_symbol', 'normalized_key']):
                hc = group.loc[group.ptv_hc]
                if hc.empty:
                    continue
                canonical = group.loc[group.is_canonical_transcript.eq('true')]
                representative = (canonical if len(canonical) else group).iloc[0]
                af = pd.to_numeric(group.gvs_eur_af, errors='coerce')
                legacy_af = (
                    pd.to_numeric(pd.Series([representative.gvs_eur_af]), errors='coerce')
                    .fillna(0)
                    .iloc[0]
                )
                known_af = af.dropna()
                strict_af = (
                    not known_af.empty and known_af.nunique() == 1 and (known_af.iloc[0] < 0.01)
                )
                records = clinvar.get(key, [])
                benign = any(
                    (
                        x['significance'] in {'Benign', 'Likely_benign', 'Benign/Likely_benign'}
                        and x['stars'] >= 2
                        for x in records
                    )
                )
                pos_flags = sorted(
                    set(
                        ';'.join(
                            support.loc[
                                (support.gene_symbol == gene) & (support.normalized_key == key),
                                'position_review_flags',
                            ]
                        ).split(';')
                    )
                    - {''}
                )
                flags = list(pos_flags)
                canonical_hc = bool(canonical.ptv_hc.any())
                if not canonical_hc:
                    flags.append('HC_PTV_supported_only_by_noncanonical')
                if af.isna().any():
                    flags.append('EUR_AF_missing_in_retained_rows')
                if known_af.nunique() > 1:
                    flags.append('EUR_AF_disagreement')
                if canonical.transcript.nunique() != 1:
                    flags.append('canonical_selection_unresolved')
                if benign:
                    flags.append('ClinVar_B_LB_at_least_2stars_requires_disease_germline_review')
                if any((x['conflict'] not in {'', '.'} for x in records)):
                    flags.append('ClinVar_conflicting_assertions_review')
                rows.append(
                    dict(
                        gene_symbol=gene,
                        normalized_key=key,
                        original_keys='|'.join(sorted(set(group.orig_key))),
                        legacy_qualifying=bool(legacy_af < 0.01),
                        canonical_ptv_hc=canonical_hc,
                        strict_known_EUR_AF_lt_1pct=bool(strict_af),
                        canonical_rule_annotation_candidate=bool(
                            canonical_hc and strict_af and (canonical.transcript.nunique() == 1)
                        ),
                        genotype_site_QC_verified=False,
                        supporting_transcripts='|'.join(sorted(set(hc.transcript))),
                        canonical_transcripts='|'.join(sorted(set(canonical.transcript))),
                        clinvar_records=json.dumps(records, ensure_ascii=False),
                        review_flags=';'.join(flags),
                    )
                )

            audit = pd.DataFrame(rows)

            audit.to_csv(
                PTV_ANNOTATION_OUT / 'ptv_definition_comparison.tsv.gz', sep='\t', index=False
            )

            audit.loc[audit.legacy_qualifying & audit.review_flags.ne('')].to_csv(
                PTV_ANNOTATION_OUT / 'legacy_ptv_manual_review.tsv.gz', sep='\t', index=False
            )

            audit.loc[audit.legacy_qualifying != audit.canonical_rule_annotation_candidate].to_csv(
                PTV_ANNOTATION_OUT / 'legacy_vs_canonical_differences.tsv.gz', sep='\t', index=False
            )

            old = pd.read_csv(
                PTV_ANNOTATION_ROOT / 'data/12_internal/12_qualifying_variants_main_ptv.csv',
                dtype=str,
            )

            mapping = f[['orig_key', 'normalized_key']].drop_duplicates()

            assert not mapping.orig_key.duplicated().any()

            old['orig_key'] = (
                old.contig + ':' + old.position + ':' + old.ref_allele + ':' + old.alt_allele
            )

            old = old.loc[old.gene_symbol.isin(cfg['panels']['baseline_9'])].merge(
                mapping, on='orig_key', how='left', validate='many_to_one'
            )

            assert old.normalized_key.notna().all(), 'Original PTV has no audited VAT mapping'

            legacy = audit.loc[
                audit.legacy_qualifying & audit.gene_symbol.isin(cfg['panels']['baseline_9'])
            ]

            assert set(zip(old.gene_symbol, old.normalized_key)) == set(
                zip(legacy.gene_symbol, legacy.normalized_key)
            )

            manifest = dict(
                status='annotation_audit_complete_final_QC_freeze_pending',
                baseline_nine_gene_exact_match=True,
                loftee_HC_actually_applied=True,
                canonical_only_matches_legacy=bool(legacy.canonical_ptv_hc.all()),
                final_PTV_definition_frozen=False,
                original_pipeline_modified=False,
                genotypes_read=False,
                clinvar_release_date_verified=header_date,
                clinvar_source_sha256_verified=True,
                clinvar_opposite_evidence_screen='aggregate B/LB >=2 stars; disease/germline and SCV/RCV manual adjudication still required',
                nmd_screen='position proxies only; not a reconstructed premature stop prediction',
                pending=[
                    'section 6 site/genotype QC specification',
                    'genotype/site QC and native ALT verification',
                    'manual transcript/NMD/ClinVar adjudication',
                    'resolve manuscript canonical versus legacy any-transcript definition',
                ],
                inputs={
                    str(p.relative_to(PTV_ANNOTATION_ROOT)): sha(p)
                    for p in [
                        source,
                        gtf,
                        cv_path,
                        PTV_ANNOTATION_ROOT / 'scripts/legacy/old_9gene_ptv_analysis.py',
                    ]
                },
                internal_outputs={
                    p.name: sha(p) for p in sorted(PTV_ANNOTATION_OUT.glob('*.tsv.gz'))
                },
            )

            (PTV_ANNOTATION_OUT / 'annotation_snapshot_manifest.json').write_text(
                json.dumps(manifest, indent=2) + '\n'
            )

            (
                PTV_ANNOTATION_ROOT / 'results/report/31_ptv_definition_audit_metadata.json'
            ).write_text(json.dumps(manifest, indent=2) + '\n')

            print(
                'PTV annotation audit complete; baseline set matched; internal review tables written; final QC freeze pending.'
            )

        # ============================================================================
        # 2. Select clinvar missense
        # ============================================================================
        if selected_task == 'missense':
            CLINVAR_OUT.mkdir(parents=True, exist_ok=True)

            cfg = yaml.safe_load(
                (CLINVAR_ROOT / 'config/rowlands_13gene_extension.yaml').read_text()
            )

            source_hash = hashlib.sha256()

            with (CLINVAR_ROOT / cfg['clinvar']['full_vcf']).open('rb') as handle:
                for block in iter(lambda: handle.read(1048576), b''):
                    source_hash.update(block)

            assert (
                source_hash.hexdigest() == cfg['clinvar']['sha256']
            ), 'Frozen ClinVar checksum mismatch'

            reference_rows = []

            fmt = '%CHROM\t%POS\t%REF\t%ALT\t%ID\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/CLNSIGCONF\t%INFO/CLNDN\t%INFO/ORIGIN\n'

            for r in pd.read_csv(
                CLINVAR_ROOT / 'data/29_rowlands_extension/gene_coordinates.csv'
            ).itertuples():
                raw = subprocess.run(
                    [
                        'bcftools',
                        'query',
                        '-r',
                        f'{r.chrom}:{r.query_region_start}-{r.query_region_end}',
                        '-f',
                        fmt,
                        str(CLINVAR_ROOT / cfg['clinvar']['normalized_vcf']),
                    ],
                    capture_output=True,
                    text=True,
                    check=True,
                ).stdout
                reference_rows.extend((l.split('\t') for l in raw.splitlines()))

            ref = pd.DataFrame(
                reference_rows,
                columns=[
                    'contig',
                    'position',
                    'ref',
                    'alt',
                    'variation_id',
                    'sig',
                    'review',
                    'conflict',
                    'disease',
                    'origin',
                ],
            ).drop_duplicates()

            ref['normalized_key'] = ref.contig + ':' + ref.position + ':' + ref.ref + ':' + ref.alt

            ref.to_csv(CLINVAR_OUT / 'clinvar_region_records.tsv.gz', sep='\t', index=False)

            ref = pd.read_csv(
                CLINVAR_OUT / 'clinvar_region_records.tsv.gz',
                sep='\t',
                dtype=str,
                keep_default_na=False,
            )

            vat = pd.read_csv(
                CLINVAR_ROOT / 'data/29_rowlands_extension/vat_candidates_by_allele_gene.tsv.gz',
                sep='\t',
                dtype=str,
                keep_default_na=False,
            )

            tx = pd.read_csv(
                CLINVAR_ROOT / 'data/29_rowlands_extension/vat_selected_transcript_rows.tsv.gz',
                sep='\t',
                dtype=str,
                keep_default_na=False,
            )

            af = (
                tx.groupby(['gene_symbol', 'normalized_key'])
                .gvs_eur_af.apply(
                    lambda x: pd.to_numeric(x, errors='coerce').dropna().unique().tolist()
                )
                .to_dict()
            )

            vat['known_consistent_af'] = [
                len(af[r.gene_symbol, r.normalized_key]) == 1
                and af[r.gene_symbol, r.normalized_key][0] < 0.01
                for r in vat.itertuples()
            ]

            miss = vat.loc[vat.canonical_missense.eq('True') & vat.known_consistent_af].copy()

            grouped = {k: g for k, g in ref.groupby('normalized_key')}

            candidates, audit = ([], [])

            for r in miss.itertuples():
                records = grouped.get(r.normalized_key)
                if records is None:
                    continue
                values = records.sig.tolist() + records.review.tolist()
                if (~records.conflict.isin(['', '.'])).any():
                    values.append('conflicting classifications')
                reason = exclusions(values)
                allowed = records.loc[
                    records.sig.map(clean).isin(PLP)
                    & records.review.map(clean).map(CLINVAR_STARS).fillna(0).ge(2)
                ]
                if not reason and allowed.empty:
                    reason = 'not_allowed_at_least_2star_germline_P_LP'
                audit.append(
                    dict(
                        gene_symbol=r.gene_symbol,
                        normalized_key=r.normalized_key,
                        stage='VCF',
                        decision=reason or 'proceed_to_assertions',
                    )
                )
                if not reason:
                    a = allowed.iloc[0]
                    row = r._asdict()
                    row.update(
                        variation_id=a.variation_id, frozen_sig=a.sig, frozen_review=a.review
                    )
                    candidates.append(row)

            candidates = pd.DataFrame(candidates)

            records = historical_xml(
                sorted(set(candidates.variation_id)), cfg['clinvar']['release_date']
            )

            eligible, evidence = ([], [])

            for row in candidates.to_dict('records'):
                v = records[row['variation_id']]
                ok, reason, rcvs, review = review_variant(
                    v, row['gene_symbol'], row['frozen_sig'], row['frozen_review']
                )
                for e in review:
                    e.update(
                        gene_symbol=row['gene_symbol'],
                        normalized_key=row['normalized_key'],
                        variation_id=row['variation_id'],
                        VCV_version=v.attrib['Version'],
                        VCV_date_last_updated=v.attrib['DateLastUpdated'],
                    )
                    evidence.append(e)
                row.update(
                    assertion_decision=reason,
                    supporting_RCVs='|'.join(rcvs),
                    VCV_accession=v.attrib['Accession'],
                    VCV_version=v.attrib['Version'],
                )
                audit.append(
                    dict(
                        gene_symbol=row['gene_symbol'],
                        normalized_key=row['normalized_key'],
                        stage='assertions',
                        decision=reason,
                    )
                )
                if ok:
                    eligible.append(row)

            write_tsv_if_changed(
                pd.DataFrame(
                    eligible,
                    columns=list(candidates.columns)
                    + ['assertion_decision', 'supporting_RCVs', 'VCV_accession', 'VCV_version'],
                ),
                CLINVAR_OUT / 'eligible_missense.tsv.gz',
            )

            pd.DataFrame(audit).to_csv(
                CLINVAR_OUT / 'missense_exclusion_audit.tsv.gz', sep='\t', index=False
            )

            pd.DataFrame(evidence).to_csv(
                CLINVAR_OUT / 'clinvar_assertion_evidence.tsv.gz', sep='\t', index=False
            )

            ptv = vat.loc[vat.ptv_hc_any_transcript.eq('True') & vat.known_consistent_af].copy()

            write_tsv_if_changed(ptv, CLINVAR_OUT / 'ptv_for_call_qc.tsv.gz')

            manifest = dict(
                status='assertion_screen_executed',
                frozen_ClinVar_release=cfg['clinvar']['release_date'],
                frozen_VCF_sha256_verified=True,
                historical_XML_at_or_before_release=True,
                exclusions_before_P_LP=True,
                disease_review='gene_specific_curated_named_syndromes_and_relevant_P_LP_RCV_at_least_2stars',
                uncertain_disease_or_inheritance_held_out=True,
                low_penetrance_hard_exclusion=True,
                explicit_recessive_assertions_checked_against_relevant_disease_mapping=True,
                final_PTV_freeze=False,
                models_run=False,
                disease_inheritance_sources={
                    'BRCA1_BRCA2': 'https://www.ncbi.nlm.nih.gov/books/NBK1247/',
                    'TP53': 'https://medlineplus.gov/genetics/condition/li-fraumeni-syndrome/',
                    'PTEN': 'https://www.ncbi.nlm.nih.gov/books/NBK1488/',
                    'CDH1': 'https://www.cancer.gov/publications/pdq/information-summaries/genetics/hereditary-diffuse-gastric-cancer-hp-pdq',
                    'NF1': 'https://www.ncbi.nlm.nih.gov/books/NBK1109/',
                    'STK11': 'https://www.ncbi.nlm.nih.gov/books/NBK1266/',
                    'breast_cancer_panel': 'https://repository.icr.ac.uk/bitstreams/5b90592e-233d-498e-bb48-207a842b40e9/download',
                },
                XML_sha256={
                    p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in XML.glob('*.xml')
                },
            )

            (CLINVAR_OUT / 'assertion_manifest.json').write_text(
                json.dumps(manifest, indent=2) + '\n'
            )

            print(
                'Historical ClinVar assertion screen completed; internal evidence and exclusion audit saved.',
                flush=True,
            )

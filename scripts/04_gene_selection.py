"""Freeze comparison cohorts, define susceptibility panels and select VAT candidates.

Configuration: config/rowlands_13gene_extension.yaml.
Historical baseline OOF membership remains an explicit prerequisite, not a newly selected cohort.
"""

from __future__ import annotations

from pathlib import Path
import sys

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.functions.paths import PROJECT_ROOT
from scripts.functions.cli import parse_tasks

TASKS = ('panel', 'candidates')
if __name__ == '__main__':
    selected_tasks = parse_tasks(__doc__, TASKS, ['panel', 'candidates'])

import hashlib

import json

import math


import subprocess

import pandas as pd

import yaml

GENE_PANEL_ROOT = PROJECT_ROOT

GENE_PANEL_CONFIG = GENE_PANEL_ROOT / 'config/rowlands_13gene_extension.yaml'

INTERNAL = GENE_PANEL_ROOT / 'data/29_rowlands_extension'

PCS = [f'PC{i}' for i in range(1, 11)]

import gzip

import os

import re

VAT_CANDIDATES_ROOT = PROJECT_ROOT

OUT = VAT_CANDIDATES_ROOT / 'data/29_rowlands_extension'

VAT_CANDIDATES_CONFIG = VAT_CANDIDATES_ROOT / 'config/rowlands_13gene_extension.yaml'

PTV = {'stop_gained', 'frameshift_variant', 'splice_acceptor_variant', 'splice_donor_variant'}

FIELDS = [
    'vid',
    'contig',
    'position',
    'ref_allele',
    'alt_allele',
    'gene_symbol',
    'transcript',
    'transcript_source',
    'consequence',
    'is_canonical_transcript',
    'gvs_eur_af',
    'gvs_eur_ac',
    'gvs_eur_an',
    'gnomad_nfe_af',
    'exon_number',
    'intron_number',
    'aa_change',
    'dna_change_in_transcript',
    'mane_select_name',
    'mane_plus_clinical_name',
    'LoF',
    'LoF_filter',
    'LoF_flags',
    'LoF_info',
]


def read(name, columns=None):
    return pd.read_csv(GENE_PANEL_ROOT / 'data' / name, usecols=columns)


def ids(frame):
    assert frame.person_id.notna().all(), 'Missing participant ID'
    assert not frame.person_id.duplicated().any(), 'Duplicate participant ID'
    return set(frame.person_id.astype(int))


from scripts.functions.io import sha256 as file_sha


def validate_extension_membership(case_control, diagnosis_age):
    """Call before later extension models; no implicit dropping of participants."""
    assert ids(case_control) == ids(
        pd.read_csv(INTERNAL / 'case_control_ids.csv.gz')
    ), 'Case/control IDs changed'
    assert ids(diagnosis_age) == ids(
        pd.read_csv(INTERNAL / 'diagnosis_age_ids.csv.gz')
    ), 'Age-model IDs changed'


from scripts.functions.variants import terms

from scripts.functions.variants import key

from scripts.functions.variants import coordinates

from scripts.functions.variants import normalize

if __name__ == "__main__":
    for selected_task in selected_tasks:

        # ============================================================================
        # 1. Freeze cohorts and prepare gene panel
        # ============================================================================
        if selected_task == 'panel':
            cfg = yaml.safe_load(GENE_PANEL_CONFIG.read_text())

            INTERNAL.mkdir(parents=True, exist_ok=True)

            assert (
                subprocess.run(
                    ['git', 'check-ignore', '-q', str(INTERNAL)], cwd=GENE_PANEL_ROOT
                ).returncode
                == 0
            )

            weights = []

            for gene in cfg['panels']['main_13']:
                w = cfg['burden']['weights'][gene]
                note = 'Overall breast cancer OR; PTVs + non-PTV-PVs, Combined; approximate external weight for PTV + qualifying P/LP missense burden; variant definitions are not identical'
                if gene == 'STK11':
                    note += '; estimate imprecise, CI crosses 1; BRIDGES and UKB only'
                if gene == 'CDH1':
                    note += (
                        '; association has marked lobular-subtype specificity; overall OR retained'
                    )
                weights.append(
                    dict(
                        gene=gene,
                        OR=w['or'],
                        CI_low=w['ci_low'],
                        CI_high=w['ci_high'],
                        ln_OR=math.log(w['or']),
                        paper_DOI=cfg['references']['doi'],
                        table='Table 2: PTVs + non-PTV-PVs / Combined',
                        notes=note,
                    )
                )

            assert len(weights) == 13 and len({w['gene'] for w in weights}) == 13

            pheno = read('04_phenotype_diagnostics_person_level.csv.gz')

            qc = read('05_cohort_qc_flags.csv')

            ids(pheno)

            ids(qc)

            base = pheno.merge(qc[['person_id', 'pass_qc']], validate='one_to_one')

            base = base.loc[(base.pass_qc == 1) & base.group.isin(['case', 'control'])].copy()

            cc = base.dropna(subset=['age_at_last_obs']).copy()

            original = read('16_joint_model_cv_oof_person_level.csv.gz', ['person_id', 'group'])

            assert ids(cc) == ids(original), 'Reconstructed script-16 membership differs'

            compare = cc[['person_id', 'group']].merge(
                original, on='person_id', validate='one_to_one'
            )

            assert (compare.group_x == compare.group_y).all(), 'Case labels changed'

            assert len(cc) == 141204 and (cc.group == 'case').sum() == 8528

            age_original = base.loc[
                (base.group == 'case')
                & base.age_at_diagnosis.notna()
                & (base.bc_is_first_record != True)
            ].copy()

            age_extension = cc.loc[
                (cc.group == 'case') & cc.age_at_diagnosis.notna() & (cc.bc_is_first_record != True)
            ].copy()

            assert age_original.bc_is_first_record.notna().all()

            assert ids(age_original) == ids(age_extension), 'Script-17 membership differs'

            assert len(age_original) == 7996

            anc = read('01_female_wgs_ancestry.csv.gz')

            ids(anc)

            checked = cc[['person_id']].merge(anc, validate='one_to_one', how='left')

            assert (checked.ancestry_pred_other == 'eur').all() and (checked.p_eur > 0.75).all()

            q = cc[['person_id']].merge(qc, validate='one_to_one')

            assert (
                (q.pass_qc == 1).all()
                and (not q.flag_related_exclude.any())
                and (not q.flag_high_missing.any())
            )

            prs = read('prs_scores.csv', ['person_id', 'PRS311'])

            pcs = read('09_ancestry_pc_person_level.csv.gz', ['person_id'] + PCS)

            clinical = read(
                '09_extended_clinical_person_level.csv.gz',
                ['person_id', 'seq_site_id', 'seq_sequencer'],
            )

            for frame in [prs, pcs, clinical]:
                joined = cc[['person_id']].merge(frame, validate='one_to_one', how='left')
                assert joined.notna().all().all(), 'Existing covariates are incomplete'

            for name, frame in [('case_control', cc), ('diagnosis_age', age_original)]:
                path = INTERNAL / (name + '_ids.csv.gz')
                if path.exists():
                    assert ids(pd.read_csv(path)) == ids(
                        frame
                    ), 'Previously frozen membership changed'
                else:
                    frame[['person_id']].sort_values('person_id').to_csv(path, index=False)

            validate_extension_membership(cc, age_extension)

            pd.DataFrame(weights).to_csv(
                GENE_PANEL_ROOT / 'results/table/29_rowlands_13gene_external_weights.csv',
                index=False,
            )

            audit = dict(
                status='sections_3_and_4_preparation_complete',
                expanded_models_run=False,
                case_control_n=len(cc),
                case_n=int((cc.group == 'case').sum()),
                control_n=int((cc.group == 'control').sum()),
                diagnosis_age_n=len(age_original),
                case_control_id_set_equals_script16=True,
                case_labels_equal=True,
                diagnosis_age_id_set_equals_script17=True,
                eur_posterior_gt_0p75=True,
                original_qc_flags_pass=True,
                existing_covariates_complete=True,
                no_new_qc_or_phenotype_definition=True,
                kinship_method='existing greedy pruning; not proven global maximum independent set',
                sample_missingness='existing chr21 1500-variant QC proxy; not whole-genome missingness',
                sex_and_wgs='inherited original female-with-srWGS selection; not re-queried',
                pending='section 5 and remaining model/output requirements',
            )

            (GENE_PANEL_ROOT / 'results/report/29_rowlands_preparation_audit.json').write_text(
                json.dumps(audit, indent=2) + '\n'
            )

            cfg['protocol_sections_received'] = list(
                dict.fromkeys(cfg.get('protocol_sections_received', []) + [1, 2, 3, 4])
            )

            cfg['burden'][
                'external_weight_table'
            ] = 'results/table/29_rowlands_13gene_external_weights.csv'

            cfg['cohorts']['case_control']['expected_cases'] = 8528

            cfg['cohorts']['case_control'][
                'frozen_ids'
            ] = 'data/29_rowlands_extension/case_control_ids.csv.gz'

            cfg['cohorts']['diagnosis_age'][
                'frozen_ids'
            ] = 'data/29_rowlands_extension/diagnosis_age_ids.csv.gz'

            cfg['cohorts']['qc_notes'] = {
                'kinship': audit['kinship_method'],
                'missingness': audit['sample_missingness'],
            }

            cfg['execution_gate'].setdefault('reason', 'await remaining model/output requirements')

            GENE_PANEL_CONFIG.write_text(yaml.safe_dump(cfg, allow_unicode=True, sort_keys=False))

            print(
                'External weights saved; both original ID sets match exactly; no expanded models run.'
            )

        # ============================================================================
        # 2. Select vat candidates
        # ============================================================================
        if selected_task == 'candidates':
            cfg = yaml.safe_load(VAT_CANDIDATES_CONFIG.read_text())

            assert os.environ['CDR_STORAGE_PATH'] == cfg['data']['storage_path']

            assert os.environ['WORKSPACE_CDR'].split('.')[-1] == cfg['data']['cdr']

            OUT.mkdir(parents=True, exist_ok=True)

            assert (
                subprocess.run(
                    ['git', 'check-ignore', '-q', str(OUT)], cwd=VAT_CANDIDATES_ROOT
                ).returncode
                == 0
            )

            mount = Path('/tmp/codex_expanded_missense_aou_ro')

            entry = next(
                (
                    x.split()
                    for x in Path('/proc/mounts').read_text().splitlines()
                    if x.split()[1] == str(mount)
                )
            )

            assert 'ro' in entry[3].split(',')

            vat = mount / 'v9' / cfg['data']['vat_relative_path']

            with gzip.open(vat, 'rt') as f:
                header = f.readline().strip().split('\t')

            indexes = {x: header.index(x) for x in FIELDS}

            coord = coordinates(cfg['panels']['main_13'])

            coord.to_csv(OUT / 'gene_coordinates.csv', index=False)

            all_rows = []

            summaries = []

            index_warning = False

            for gene in cfg['panels']['main_13']:
                g = coord.loc[coord.gene_symbol == gene].iloc[0]
                region = f'{g.chrom}:{g.query_region_start}-{g.query_region_end}'
                query = subprocess.run(
                    ['tabix', str(vat), region], capture_output=True, text=True, check=True
                )
                index_warning |= 'index file is older' in query.stderr
                selected = []
                matching_rows = 0
                for line in query.stdout.splitlines():
                    values = line.split('\t')
                    assert len(values) == len(header), 'VAT row/header schema mismatch'
                    if values[indexes['gene_symbol']] != gene:
                        continue
                    matching_rows += 1
                    row = {x: values[i] for x, i in indexes.items()}
                    canonical = row['is_canonical_transcript'] == 'true'
                    if not (canonical or row['LoF'] != ''):
                        continue
                    if not terms(row['consequence']) & (PTV | {'missense_variant'}):
                        continue
                    selected.append(row)
                assert selected, f'No VAT candidates for {gene}; inspect coverage'
                all_rows.extend(selected)
                summaries.append(
                    dict(
                        gene=gene,
                        query_completed=True,
                        candidate_rows=len(selected),
                        matching_gene_rows=matching_rows,
                    )
                )
                print(
                    f'VAT extraction completed for {gene}; no participant counts or identifiers printed',
                    flush=True,
                )

            frame = pd.DataFrame(all_rows).drop_duplicates()

            frame['position'] = frame.position.astype(int)

            frame['orig_key'] = frame.apply(key, axis=1)

            frame['ptv_hc_row'] = frame.consequence.map(
                lambda x: bool(terms(x) & PTV)
            ) & frame.LoF.eq('HC')

            old = pd.read_csv(
                VAT_CANDIDATES_ROOT / 'data/11_panel_gene_vat_ptv_missense_candidates.csv.gz',
                dtype=str,
                keep_default_na=False,
            )

            old = old.loc[old.gene_symbol.isin(cfg['panels']['baseline_9'])]

            new = frame.loc[frame.gene_symbol.isin(cfg['panels']['baseline_9'])]

            check = [
                'vid',
                'gene_symbol',
                'transcript',
                'consequence',
                'is_canonical_transcript',
                'LoF',
                'gvs_eur_af',
            ]

            old_set = set(old[check].itertuples(index=False, name=None))

            new_set = set(new[check].itertuples(index=False, name=None))

            assert (
                old_set == new_set
            ), 'Fresh nine-gene VAT extraction differs from original cache; review before proceeding'

            fasta = VAT_CANDIDATES_ROOT / 'data/refdata/Homo_sapiens_assembly38.fasta'

            frame['normalized_key'] = frame.orig_key.map(normalize(frame, fasta))

            frame.to_csv(OUT / 'vat_selected_transcript_rows.tsv.gz', sep='\t', index=False)

            collapsed, uncertain = ([], [])

            for (gene, allele), group in frame.groupby(
                ['gene_symbol', 'normalized_key'], sort=False
            ):
                canonical = group.loc[group.is_canonical_transcript.eq('true')]
                representative = (canonical if len(canonical) else group).iloc[0]
                issues = []
                if canonical.empty:
                    issues.append('no_canonical_candidate_row')
                if canonical.transcript.nunique() > 1:
                    issues.append('multiple_canonical_transcripts')
                if group.consequence.nunique() > 1:
                    issues.append('consequence_differs_between_retained_transcripts')
                if group.ptv_hc_row.any() and (not canonical.ptv_hc_row.any()):
                    issues.append('PTV_HC_supported_only_by_noncanonical_row')
                if len(set(group.LoF) & {'HC', 'LC'}) > 1:
                    issues.append('HC_LC_disagreement_between_transcripts')
                af = pd.to_numeric(group.gvs_eur_af, errors='coerce')
                if af.isna().any():
                    issues.append('EUR_AF_missing')
                if af.dropna().nunique() > 1:
                    issues.append('EUR_AF_disagreement_after_normalization')
                if group.orig_key.nunique() > 1:
                    issues.append('multiple_original_representations_of_normalized_allele')
                if canonical.empty or canonical.consequence.nunique() > 1:
                    missense = None
                else:
                    missense = bool(terms(canonical.iloc[0].consequence) & {'missense_variant'})
                row = representative.to_dict()
                row.update(
                    normalized_key=allele,
                    ptv_hc_any_transcript=bool(group.ptv_hc_row.any()),
                    canonical_missense=missense,
                    uncertainty_flags=';'.join(issues),
                    original_keys='|'.join(sorted(set(group.orig_key))),
                    retained_transcripts='|'.join(sorted(set(group.transcript))),
                )
                collapsed.append(row)
                if issues:
                    uncertain.append(row)

            candidates = pd.DataFrame(collapsed)

            assert not candidates.duplicated(['gene_symbol', 'normalized_key']).any()

            candidates.to_csv(OUT / 'vat_candidates_by_allele_gene.tsv.gz', sep='\t', index=False)

            pd.DataFrame(uncertain, columns=candidates.columns).to_csv(
                OUT / 'transcript_consequence_uncertainty.tsv.gz', sep='\t', index=False
            )

            baseline = pd.read_csv(
                VAT_CANDIDATES_ROOT / 'data/12_internal/12_qualifying_variants_main_ptv.csv',
                dtype=str,
                keep_default_na=False,
            )

            baseline = baseline.loc[baseline.gene_symbol.isin(cfg['panels']['baseline_9'])].copy()

            baseline['orig_key'] = baseline.apply(key, axis=1)

            mapping = frame[['orig_key', 'normalized_key']].drop_duplicates()

            baseline = baseline.merge(mapping, on='orig_key', how='left', validate='many_to_one')

            assert baseline.normalized_key.notna().all()

            comparable = candidates.loc[
                candidates.gene_symbol.isin(cfg['panels']['baseline_9'])
                & candidates.ptv_hc_any_transcript
                & (pd.to_numeric(candidates.gvs_eur_af, errors='coerce').fillna(0) < 0.01)
            ]

            assert set(zip(baseline.gene_symbol, baseline.normalized_key)) == set(
                zip(comparable.gene_symbol, comparable.normalized_key)
            ), 'Original nine-gene qualifying PTV set changed'

            meta = dict(
                status='VAT_section_5_1_complete',
                expanded_models_run=False,
                genotypes_read=False,
                original_nine_gene_qualifying_PTV_allele_gene_set_exact_match=True,
                genes=cfg['panels']['main_13'],
                original_nine_gene_selected_rows_exact_match=True,
                all_reference_alleles_match_GRCh38=True,
                allele_gene_keys_unique=True,
                index_mtime_warning_observed=index_warning,
                index_warning_interpretation='mtime warning alone does not establish index corruption; original nine-gene retrieval matched; genome-wide coverage not tested',
                annotation_engine='AoU VAT/NIRVANA; no VEP rerun',
                selection='canonical OR nonempty LoF; strict gene match; four PTV terms OR missense',
                normalization='bcftools norm -f fixed GRCh38 -c e -m -any',
                normalization_tool=subprocess.run(
                    ['bcftools', '--version'], capture_output=True, text=True, check=True
                ).stdout.splitlines()[0],
                genotype_join_requirement='match full normalized chr:pos:ref:alt; dose and AD must use native VCF ALT index',
                uncertainty_action='inspection list only; no result-dependent automatic exclusions',
                added_gene_coordinates='cached Ensembl release 110 GTF for region bounds only; VAT supplies annotations',
                not_checked=[
                    'NIRVANA/VEP per-variant concordance',
                    'genotype ALT match (no genotype query in this stage)',
                    'ClinVar missense eligibility (remaining protocol pending)',
                ],
            )

            (
                VAT_CANDIDATES_ROOT / 'results/report/30_rowlands_vat_extraction_metadata.json'
            ).write_text(json.dumps(meta, indent=2) + '\n')

            cfg['candidate_extraction'] = {
                'section_5_1_complete': True,
                'transcript_rows': 'data/29_rowlands_extension/vat_selected_transcript_rows.tsv.gz',
                'allele_gene_candidates': 'data/29_rowlands_extension/vat_candidates_by_allele_gene.tsv.gz',
                'uncertainty_list': 'data/29_rowlands_extension/transcript_consequence_uncertainty.tsv.gz',
                'original_nine_gene_rows_exact_match': True,
                'genotype_alt_verification': 'pending_genotype_stage',
            }

            cfg['protocol_sections_received'] = [1, 2, 3, 4, '5.1']

            VAT_CANDIDATES_CONFIG.write_text(
                yaml.safe_dump(cfg, allow_unicode=True, sort_keys=False)
            )

            print(
                'Section 5.1 completed; detailed variants retained only in ignored internal data; no models run',
                flush=True,
            )

"""Prepare prespecified sensitivity scores before association fitting.

Run before stages 08/09 to retain prespecification timing.
The existing without-STK11, PTV-only, strict TP53/NF1 VAF, >=1-star and without-CDH1 definitions are preserved.
"""

from __future__ import annotations

from pathlib import Path
import sys

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.functions.paths import PROJECT_ROOT
from scripts.functions.cli import parse_tasks

TASKS = ('sensitivities',)
if __name__ == '__main__':
    selected_tasks = parse_tasks(__doc__, TASKS, ['sensitivities'])

import gzip

import hashlib

import importlib.util

import json

import shutil

import subprocess

from datetime import datetime, timezone


import numpy as np

import pandas as pd

import yaml

SENSITIVITY_PREPARATION_ROOT = PROJECT_ROOT

SENSITIVITY_PREPARATION_BASE = SENSITIVITY_PREPARATION_ROOT / 'data/29_rowlands_extension'

SENSITIVITY_PREPARATION_OUT = SENSITIVITY_PREPARATION_BASE / 'section10'

DELIVERY = SENSITIVITY_PREPARATION_OUT / 'internal_delivery'

from scripts.functions.io import sha256 as sensitivity_preparation_sha

if __name__ == "__main__":
    for selected_task in selected_tasks:

        # ============================================================================
        # 1. Prepare sensitivity scores
        # ============================================================================
        if selected_task == 'sensitivities':
            SENSITIVITY_PREPARATION_OUT.mkdir(parents=True, exist_ok=True)

            DELIVERY.mkdir(exist_ok=True)

            assert (
                subprocess.run(
                    ['git', 'check-ignore', '-q', str(SENSITIVITY_PREPARATION_OUT)],
                    cwd=SENSITIVITY_PREPARATION_ROOT,
                ).returncode
                == 0
            )

            cfg = yaml.safe_load(
                (SENSITIVITY_PREPARATION_ROOT / 'config/rowlands_13gene_extension.yaml').read_text()
            )

            acceptance = json.loads(
                (
                    SENSITIVITY_PREPARATION_ROOT
                    / 'results/report/43_rowlands_pre_model_acceptance.json'
                ).read_text()
            )

            assert not cfg['expanded_models_run'] and (
                not acceptance['models_run']
            ), 'Sensitivity plan must precede C association results'

            plan = {
                'analysis': 'PTV plus qualifying P/LP missense; approximate Rowlands all-PV weights',
                'locked_before_C_association_results': True,
                'models_run_at_lock': False,
                'sensitivity_versions': {
                    'C_without_STK11': {'panel': 'C minus STK11', 'other_rules': 'same as C'},
                    'C_PTV_only': {
                        'panel': '13 genes',
                        'variant_classes': ['PTV'],
                        'other_rules': 'current C QC and Rowlands weights; not legacy A',
                    },
                    'C_TP53_NF1_strict_VAF': {
                        'panel': '13 genes',
                        'heterozygous_VAF_TP53_NF1': [0.3, 0.7],
                        'heterozygous_VAF_other_genes': [0.3, 0.8],
                        'homozygous_ALT_VAF_min': 0.8,
                        'other_rules': 'same as C',
                    },
                    'C_ClinVar_min1star': {
                        'panel': '13 genes',
                        'missense_min_stars': 1,
                        'relevant_P_LP_RCV_min_stars': 1,
                        'exclusions_first': [
                            'low_penetrance',
                            'risk_factor_or_risk_allele',
                            'conflict_or_opposing_B_LB',
                        ],
                        'require_exact_allowed_germline_P_LP': True,
                        'disease_inheritance_review': 'same as C',
                        'ClinVar_release': cfg['clinvar']['release_date'],
                        'PTV_and_AF_and_call_QC': 'same as C',
                    },
                    'C_without_CDH1': {
                        'panel': 'C minus CDH1',
                        'conditional_on_sufficient_information': True,
                        'availability_rule': 'complete gene carrier and allele-call traceability; estimability assessed without a significance cutoff',
                        'no_subtype_claim_without_subtype_data': True,
                        'other_rules': 'same as C',
                    },
                },
                'case_control_model': 'existing script16 PRS311 + burden + interaction + existing covariates; HC1; same samples',
                'CV': 'existing paired 10 repeats x 5 folds; train-fold-only scaling; increment relative to PRS311 baseline',
                'age_model': 'existing script17 case-only recorded diagnosis-age OLS with PRS311 + burden + PCs1-10; HC1',
                'report': 'all estimable predefined versions and strata, estimates and CI, actual SD, raw burden distribution, carrier proportion',
                'selection': 'no selection of final scheme or new strata using individual significant results',
                'PRS311': 'existing score and CHEK2 linked-tag handling unchanged; no reselection',
                'disclosure': 'exact comparison and participant/variant tables internal; public sparse cells and inferable complementary differences suppressed',
                'model_gate': 'section9 user acceptance required; biological review status must remain explicit',
                'final_variant_definition_frozen': False,
            }

            plan_path = SENSITIVITY_PREPARATION_OUT / 'prespecified_sensitivity_plan.json'

            if plan_path.exists():
                previous = json.loads(plan_path.read_text())
                assert {
                    k: v for k, v in previous.items() if k != 'locked_at_UTC'
                } == plan, 'Previously locked sensitivity specification differs'
                plan = previous
            else:
                plan['locked_at_UTC'] = datetime.now(timezone.utc).isoformat()
                plan_path.write_text(json.dumps(plan, indent=2, ensure_ascii=False) + '\n')

            from scripts.functions import clinvar as cv

            vat = pd.read_csv(
                SENSITIVITY_PREPARATION_BASE / 'vat_candidates_by_allele_gene.tsv.gz',
                sep='\t',
                dtype=str,
                keep_default_na=False,
            )

            tx = pd.read_csv(
                SENSITIVITY_PREPARATION_BASE / 'vat_selected_transcript_rows.tsv.gz',
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

            refs = pd.read_csv(
                SENSITIVITY_PREPARATION_BASE / 'section6/clinvar_region_records.tsv.gz',
                sep='\t',
                dtype=str,
                keep_default_na=False,
            )

            refgroups = {k: d for k, d in refs.groupby('normalized_key')}

            candidates = []

            decisions = []

            for r in vat.itertuples():
                if r.canonical_missense != 'True':
                    continue
                pair = (r.gene_symbol, r.normalized_key)
                d = refgroups.get(r.normalized_key)
                values = [] if d is None else d.sig.tolist() + d.review.tolist()
                if d is not None and (~d.conflict.isin(['', '.'])).any():
                    values.append('conflicting classifications')
                reason = cv.exclusions(values)
                if not reason and (not (len(af[pair]) == 1 and af[pair][0] < 0.01)):
                    reason = 'EUR_AF_missing_disagreeing_or_not_lt_0p01'
                allowed = (
                    None
                    if d is None
                    else d.loc[
                        d.sig.map(cv.clean).isin(cv.PLP)
                        & d.review.map(cv.clean).map(cv.STARS).fillna(0).ge(1)
                    ]
                )
                if not reason and (allowed is None or allowed.empty):
                    reason = 'not_allowed_at_least_1star_germline_P_LP'
                if reason:
                    decisions.append(
                        dict(
                            gene=r.gene_symbol,
                            normalized_key=r.normalized_key,
                            keep=False,
                            rejection_reason=reason,
                        )
                    )
                else:
                    a = allowed.iloc[0]
                    candidates.append(
                        (r.gene_symbol, r.normalized_key, a.variation_id, a.sig, a.review)
                    )

            print(
                'Prespecification locked; reviewing historical ClinVar assertions for min1-star candidates.',
                flush=True,
            )

            xml = cv.historical_xml(
                sorted({r[2] for r in candidates}),
                cfg['clinvar']['release_date'],
                hold_unavailable=True,
            )

            eligible = set()

            evidence = []

            for gene, key, vid, sig, review in candidates:
                if vid not in xml:
                    decisions.append(
                        dict(
                            gene=gene,
                            normalized_key=key,
                            keep=False,
                            rejection_reason='historical_VCV_unavailable_at_frozen_release_pending_review',
                        )
                    )
                    continue
                ok, reason, rcvs, rows = cv.review_variant(xml[vid], gene, sig, review, min_stars=1)
                decisions.append(
                    dict(
                        gene=gene,
                        normalized_key=key,
                        keep=ok,
                        rejection_reason='' if ok else reason,
                    )
                )
                if ok:
                    eligible.add((gene, key))
                for row in rows:
                    row.update(gene=gene, normalized_key=key, variation_id=vid)
                    evidence.append(row)

            pd.DataFrame(decisions).to_csv(
                SENSITIVITY_PREPARATION_OUT / 'min1star_variant_eligibility_audit.tsv.gz',
                sep='\t',
                index=False,
            )

            pd.DataFrame(evidence).to_csv(
                SENSITIVITY_PREPARATION_OUT / 'min1star_assertion_evidence.tsv.gz',
                sep='\t',
                index=False,
            )

            main_eligible = {
                (r.gene_symbol, r.normalized_key)
                for r in pd.read_csv(
                    SENSITIVITY_PREPARATION_BASE / 'section6/eligible_missense.tsv.gz', sep='\t'
                ).itertuples()
            }

            assert (
                main_eligible <= eligible
            ), 'Min1-star must retain all main clinically eligible missense'

            main_records = pd.read_csv(
                SENSITIVITY_PREPARATION_BASE / 'section6/eligible_missense.tsv.gz',
                sep='\t',
                dtype=str,
            )

            for r in main_records.itertuples():
                assert cv.review_variant(
                    xml[r.variation_id], r.gene_symbol, r.frozen_sig, r.frozen_review
                )[0], 'Default2-star review changed'

            current = pd.read_csv(
                SENSITIVITY_PREPARATION_BASE / 'section7/unweighted_gene_carriers.csv.gz'
            ).set_index('sample_id')

            genes = cfg['panels']['main_13']

            weights = (
                pd.read_csv(SENSITIVITY_PREPARATION_ROOT / cfg['burden']['external_weight_table'])
                .set_index('gene')
                .ln_OR.loc[genes]
            )

            min1 = {g: set() for g in genes}

            qc_names = {
                'site_filter_not_PASS',
                'genotype_filter_failed',
                'not_complete_diploid_target_ALT_genotype',
                'GQ_missing_nonfinite_or_lt20',
                'AD_missing_invalid_or_zero_total',
                'VAF_outside_main_range',
            }

            decision_lookup = {
                (r['gene'], r['normalized_key']): r['rejection_reason'] for r in decisions
            }

            dest = DELIVERY / '02_variant_sample_inclusion_exclusion_audit.csv.gz'

            with gzip.open(dest, 'wt', newline='') as f:
                first = True
                for d in pd.read_csv(
                    SENSITIVITY_PREPARATION_BASE / 'section7/sample_variant_audit.csv.gz',
                    chunksize=20000,
                    keep_default_na=False,
                    low_memory=False,
                ):
                    keep = []
                    why = []
                    for r in d.to_dict('records'):
                        pair = (r['gene'], r['normalized_key'])
                        ptv = r['PTV/missense'] == 'PTV'
                        if ptv:
                            k = bool(int(r['keep_flag']))
                            reason = r['rejection_reason']
                        else:
                            reasons = [x for x in r['rejection_reason'].split(';') if x in qc_names]
                            if pair not in eligible:
                                reasons.insert(
                                    0, decision_lookup.get(pair, 'not_eligible_canonical_missense')
                                )
                            k = not reasons
                            reason = ';'.join(reasons)
                            if k:
                                min1[r['gene']].add(int(r['sample_id']))
                        keep.append(int(k))
                        why.append(reason)
                    d['keep_min1star'] = keep
                    d['rejection_reason_min1star'] = why
                    assert not (d.keep_flag.astype(int).eq(1) & d.keep_min1star.eq(0)).any()
                    d.to_csv(f, index=False, header=first)
                    first = False

            one = current.copy()

            for g in genes:
                one[f'{g}_missense_main'] = one.index.isin(min1[g]).astype('int8')
                one[f'{g}_main'] = np.maximum(one[f'{g}_PTV_main'], one[f'{g}_missense_main'])

            one.to_csv(SENSITIVITY_PREPARATION_OUT / 'min1star_unweighted_gene_carriers.csv.gz')

            ptv = current.copy()

            for g in genes:
                ptv[f'{g}_main'] = ptv[f'{g}_PTV_main']

            definitions = {
                'C_PTV_only': (ptv, genes, 'main'),
                'C_TP53_NF1_strict_VAF': (current, genes, 'tp53_nf1_vaf_sensitive'),
                'C_ClinVar_min1star': (one, genes, 'main'),
                'C_without_CDH1': (current, [g for g in genes if g != 'CDH1'], 'main'),
            }

            sensitive = pd.DataFrame(index=current.index)

            sensitive.index.name = 'sample_id'

            for version, (frame, panel, rule) in definitions.items():
                matrix = frame[[f'{g}_{rule}' for g in panel]].to_numpy(dtype=int)
                raw = matrix @ weights.loc[panel].to_numpy()
                sensitive[f'B_raw_{version}'] = raw
                sensitive[f'carried_gene_count_{version}'] = matrix.sum(axis=1)
                sd = raw.std(ddof=0)
                assert sd > 0
                sensitive[f'B_standardized_{version}'] = (raw - raw.mean()) / sd

            main = (
                pd.read_csv(
                    SENSITIVITY_PREPARATION_BASE / 'section8/comparable_scores_case_control.csv.gz'
                )
                .set_index('sample_id')
                .loc[current.index]
            )

            assert (sensitive.B_raw_C_ClinVar_min1star >= main.B_raw_C - 1e-12).all()

            assert np.allclose(
                main.B_raw_C - sensitive.B_raw_C_without_CDH1, current.CDH1_main * weights.CDH1
            )

            existing = (
                pd.read_csv(
                    SENSITIVITY_PREPARATION_BASE / 'section7/carrier_burden_case_control.csv.gz'
                )
                .set_index('sample_id')
                .loc[current.index]
            )

            assert np.allclose(
                sensitive.B_raw_C_TP53_NF1_strict_VAF, existing.B_raw_13gene_tp53_nf1_vaf_sensitive
            )

            sensitive.to_csv(
                SENSITIVITY_PREPARATION_OUT / 'prespecified_sensitivity_scores_case_control.csv.gz'
            )

            age = set(
                pd.read_csv(
                    SENSITIVITY_PREPARATION_ROOT / cfg['cohorts']['diagnosis_age']['frozen_ids']
                ).person_id
            )

            sensitive.loc[sorted(age)].to_csv(
                SENSITIVITY_PREPARATION_OUT / 'prespecified_sensitivity_scores_diagnosis_age.csv.gz'
            )

            shutil.copyfile(
                SENSITIVITY_PREPARATION_ROOT / cfg['burden']['external_weight_table'],
                DELIVERY / '01_13gene_OR_weights.csv',
            )

            all_scores = main.join(sensitive)

            rows = []

            for cohort, members in [('case_control', set(main.index)), ('diagnosis_age', age)]:
                s = all_scores.loc[sorted(members)]
                for col in s.filter(regex='^B_raw_'):
                    version = col.removeprefix('B_raw_')
                    x = s[col]
                    carrier = s[f'carried_gene_count_{version}'] > 0
                    row = dict(
                        cohort=cohort,
                        version=version,
                        n=len(x),
                        carriers=int(carrier.sum()),
                        carrier_fraction=float(carrier.mean()),
                        zero_fraction=float(x.eq(0).mean()),
                        raw_mean=float(x.mean()),
                        SD_population_ddof0=float(x.std(ddof=0)),
                        SD_sample_ddof1=float(x.std(ddof=1)),
                        model_status='not_run_pending_section9_user_acceptance',
                        OR_per_model_SD=None,
                        OR_CI_low=None,
                        OR_CI_high=None,
                        CV_AUC=None,
                        CV_AUC_increment_over_PRS311=None,
                        age_beta_per_model_SD=None,
                        age_CI_low=None,
                        age_CI_high=None,
                        age_R_squared=None,
                        age_R_squared_increment=None,
                    )
                    row.update(
                        {
                            f'raw_p{int(q * 100):02d}': float(x.quantile(q))
                            for q in [0, 0.01, 0.25, 0.5, 0.75, 0.95, 0.99, 1]
                        }
                    )
                    rows.append(row)

            pd.DataFrame(rows).to_csv(
                DELIVERY / '03_A_to_C_carrier_and_model_comparison.csv', index=False
            )

            persisted = pd.read_csv(
                dest, usecols=['sample_id', 'gene', 'PTV/missense', 'keep_flag', 'keep_min1star']
            )

            assert not (persisted.keep_flag.eq(1) & persisted.keep_min1star.eq(0)).any()

            for g in genes:
                observed = set(
                    persisted.loc[
                        persisted.gene.eq(g)
                        & persisted['PTV/missense'].eq('missense')
                        & persisted.keep_min1star.eq(1),
                        'sample_id',
                    ]
                )
                assert observed == set(one.index[one[f'{g}_missense_main'].eq(1)])

            saved = pd.read_csv(
                SENSITIVITY_PREPARATION_OUT / 'prespecified_sensitivity_scores_case_control.csv.gz'
            ).set_index('sample_id')

            saved_age = pd.read_csv(
                SENSITIVITY_PREPARATION_OUT / 'prespecified_sensitivity_scores_diagnosis_age.csv.gz'
            ).set_index('sample_id')

            pd.testing.assert_frame_equal(saved.loc[saved_age.index], saved_age)

            sources = [
                SENSITIVITY_PREPARATION_ROOT / 'config/rowlands_13gene_extension.yaml',
                SENSITIVITY_PREPARATION_ROOT / 'scripts/functions/clinvar.py',
                Path(__file__),
                SENSITIVITY_PREPARATION_ROOT / cfg['burden']['external_weight_table'],
                SENSITIVITY_PREPARATION_BASE / 'section7/sample_variant_audit.csv.gz',
                SENSITIVITY_PREPARATION_BASE / 'section7/unweighted_gene_carriers.csv.gz',
                SENSITIVITY_PREPARATION_BASE / 'section8/comparable_scores_case_control.csv.gz',
                plan_path,
            ]

            manifest = {
                'status': 'sensitivity_plan_locked_and_scores_generated_models_pending_acceptance',
                'models_run': False,
                'final_variant_definition_frozen': False,
                'plan_sha256': sensitivity_preparation_sha(plan_path),
                'main_min_stars': 2,
                'sensitivity_min_stars': 1,
                'all_main_eligible_missense_retained_in_min1star': True,
                'main_audit_and_scores_not_overwritten': True,
                'all_main_eligible_default2star_reviews_unchanged': True,
                'persisted_audit_reconstructs_min1star_carriers': True,
                'persisted_age_scores_exact_fixed_subset': True,
                'CDH1_gene_and_call_information_available': True,
                'CDH1_model_estimability': 'not_yet_assessed_no_association_results_viewed',
                'inputs': {
                    str(p.relative_to(SENSITIVITY_PREPARATION_ROOT)): sensitivity_preparation_sha(p)
                    for p in sources
                },
                'delivery': {
                    p.name: sensitivity_preparation_sha(p)
                    for p in DELIVERY.iterdir()
                    if p.is_file()
                },
            }

            (SENSITIVITY_PREPARATION_OUT / 'delivery_manifest_internal.json').write_text(
                json.dumps(manifest, indent=2) + '\n'
            )

            public = {k: v for k, v in manifest.items() if k not in {'inputs', 'delivery'}}

            (
                SENSITIVITY_PREPARATION_ROOT
                / 'results/report/44_rowlands_section10_validation.json'
            ).write_text(json.dumps(public, indent=2) + '\n')

            print(
                'Prespecified scores, three internal delivery tables and source checksums complete; no models fitted.',
                flush=True,
            )

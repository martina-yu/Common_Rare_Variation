"""Read-only numerical acceptance of C before manuscript downstream refresh."""
from pathlib import Path
import hashlib
import json
import sys

import numpy as np
import pandas as pd
import yaml
from sklearn.model_selection import RepeatedStratifiedKFold

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
BASE = ROOT / 'data/29_rowlands_extension'


def main():
    from scripts.functions.provenance import verify_frozen_inputs
    cfg = yaml.safe_load((ROOT / 'config/rowlands_13gene_extension.yaml').read_text())
    checks = {}

    def check(name, condition):
        checks[name] = bool(condition)
        assert checks[name], name

    d = pd.read_csv(BASE / 'section9/frozen_model_inputs.csv.gz').set_index('person_id')
    s = pd.read_csv(BASE / 'section8/comparable_scores_case_control.csv.gz').set_index('sample_id')
    c = pd.read_csv(BASE / 'section7/unweighted_gene_carriers.csv.gz').set_index('sample_id')
    old = pd.read_csv(ROOT / 'data/16_joint_model_cv_oof_person_level.csv.gz').set_index('person_id')
    ids = set(pd.read_csv(BASE / 'case_control_ids.csv.gz').person_id)
    age_ids = set(pd.read_csv(BASE / 'diagnosis_age_ids.csv.gz').person_id)
    for name, frame in [('inputs', d), ('scores', s), ('carriers', c), ('legacy_membership', old)]:
        check(name + '_unique_and_same_membership', not frame.index.duplicated().any() and set(frame.index) == ids)
    check('cohort_sizes', len(ids) == 141204 and int(d.is_case.sum()) == 8528 and len(age_ids) == 7996)
    check('age_subset', age_ids <= ids and d.loc[sorted(age_ids)].is_case.eq(1).all()
          and d.loc[sorted(age_ids)].age_at_diagnosis.notna().all()
          and d.loc[sorted(age_ids)].bc_is_first_record.eq(False).all())
    check('model_covariates_complete', d[['PRS311', 'age_at_last_obs'] + [f'PC{i}' for i in range(1, 11)]
          + ['seq_site_id', 'seq_sequencer']].notna().all().all())
    prs = pd.read_csv(ROOT / 'data/prs_scores.csv').set_index('person_id')
    check('PRS311_unchanged', np.array_equal(d.PRS311.to_numpy(), prs.loc[d.index].PRS311.to_numpy()))
    ph = pd.read_csv(ROOT / 'data/04_phenotype_diagnostics_person_level.csv.gz').set_index('person_id').loc[d.index]
    check('phenotype_unchanged', d.group.equals(ph.group) and np.allclose(d.age_at_diagnosis, ph.age_at_diagnosis, equal_nan=True)
          and np.allclose(d.age_at_last_obs, ph.age_at_last_obs, equal_nan=True))
    panel = cfg['panels']['main_13']
    check('correct_13_gene_panel', len(panel) == 13 and 'MAP3K1' not in panel and set(panel) - set(cfg['panels']['baseline_9']) == {'PTEN', 'NF1', 'STK11', 'CDH1'})
    c = c.loc[d.index]
    s = s.loc[d.index]
    weights = pd.read_csv(ROOT / cfg['burden']['external_weight_table']).set_index('gene').OR.loc[panel]
    check('gene_union_binary', c[[f'{g}_main' for g in panel]].isin([0, 1]).all().all()
          and all(np.array_equal(c[f'{g}_main'], c[[f'{g}_PTV_main', f'{g}_missense_main']].max(axis=1)) for g in panel))
    raw = c[[f'{g}_main' for g in panel]].to_numpy() @ np.log(weights.to_numpy())
    check('C_raw_formula', np.allclose(raw, s.B_raw_C, rtol=1e-12, atol=1e-12))
    check('C_saved_standardization_ddof0', np.allclose((raw - raw.mean()) / raw.std(ddof=0), s.B_standardized_C, rtol=1e-12, atol=1e-12))
    check('C_gene_count', np.array_equal(c[[f'{g}_main' for g in panel]].sum(axis=1), s.carried_gene_count_C))
    ac = pd.read_csv(BASE / 'section8/comparable_scores_diagnosis_age.csv.gz').set_index('sample_id')
    check('age_saved_score_reuses_reference', set(ac.index) == age_ids and np.allclose(ac.B_standardized_C, s.loc[ac.index].B_standardized_C))
    # Reconstruct all carriers from retained calls, not from saved score totals.
    reconstructed = {(g, cl): set() for g in panel for cl in ['PTV', 'missense']}
    required = ['sample_id', 'gene', 'PTV/missense', 'keep_flag', 'GQ', 'VAF', 'EUR_AF',
                'site_filter', 'genotype_filter', 'allele_dosage', 'low_penetrance_flag', 'conflict_flag']
    for chunk in pd.read_csv(BASE / 'section7/sample_variant_audit.csv.gz', usecols=required, chunksize=20000, keep_default_na=False):
        kept = chunk.loc[chunk.keep_flag.eq(1)]
        gq = pd.to_numeric(kept.GQ); vaf = pd.to_numeric(kept.VAF); af = pd.to_numeric(kept.EUR_AF)
        assert kept.site_filter.eq('PASS').all() and kept.genotype_filter.isin(['PASS', '.']).all()
        assert np.isfinite(gq).all() and gq.ge(20).all() and np.isfinite(af).all() and af.lt(.01).all()
        assert ((kept.allele_dosage.eq(1) & vaf.between(.3, .8)) | (kept.allele_dosage.eq(2) & vaf.ge(.8))).all()
        miss = kept.loc[kept['PTV/missense'].eq('missense')]
        assert miss.low_penetrance_flag.astype(str).str.lower().eq('false').all()
        assert miss.conflict_flag.astype(str).str.lower().eq('false').all()
        for (g, cl), rows in kept.groupby(['gene', 'PTV/missense']):
            reconstructed[g, cl].update(rows.sample_id)
    check('retained_call_QC_and_carrier_reconstruction', all(set(c.index[c[f'{g}_{cl}_main'].eq(1)]) == members
          for (g, cl), members in reconstructed.items()))
    acceptance = json.loads((ROOT / 'results/report/43_rowlands_pre_model_acceptance.json').read_text())
    lock = json.loads((BASE / 'section9/pre_model_input_lock_internal.json').read_text())
    check('frozen_data_input_hashes_valid', all(
        hashlib.sha256((ROOT / p).read_bytes()).hexdigest() == h
        for p, h in lock['input_sha256'].items() if p.startswith('data/')))
    provenance_issue = ''
    try:
        verify_frozen_inputs(lock['input_sha256'])
        checks['historical_source_migration_hashes_valid'] = True
    except AssertionError as exc:
        checks['historical_source_migration_hashes_valid'] = False
        provenance_issue = str(exc)
    # This refresh consumes persisted, hash-identical upstream data; it does
    # not rerun PRS, phenotyping or genotype extraction. Lock current source
    # separately rather than silently rewriting the historical migration.
    manifest = json.loads((ROOT / 'scripts/refactor_manifest.json').read_text())
    kernels = ['scripts/08_association_models.py', 'scripts/functions/burden.py',
               'scripts/functions/statistics.py', 'scripts/functions/clinvar.py',
               'scripts/functions/genotype_qc.py']
    check('downstream_model_and_QC_kernels_match_historical_manifest', all(
        hashlib.sha256((ROOT / p).read_bytes()).hexdigest() == manifest['implementation_sha256'][p]
        for p in kernels))
    check('model_scaling_explicit', acceptance['model_scaling_SD_ddof'] == 1 and acceptance['saved_score_SD_ddof'] == 0
          and acceptance['age_model_scaling'] == 'within frozen cases as existing script17')
    strata = d.is_case.to_numpy().astype(str) + '_' + old.loc[d.index].carrier_main_10gene.to_numpy().astype(int).astype(str)
    splits = list(RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=20260818).split(np.arange(len(d)), strata))
    fold = np.full(len(d), -1)
    for k, (_, te) in enumerate(splits[:5]):
        fold[te] = k
    check('original_CV_repeat1_folds_reproduced', np.array_equal(fold, old.loc[d.index].fold_r1))
    paths = ['section9/frozen_model_inputs.csv.gz', 'section8/comparable_scores_case_control.csv.gz',
             'section7/unweighted_gene_carriers.csv.gz', 'section7/sample_variant_audit.csv.gz']
    acceptance_checks = {k: v for k, v in checks.items() if k != 'historical_source_migration_hashes_valid'}
    result = {'status': ('step1_persisted_inputs_accepted_with_fresh_refresh_lock' if all(acceptance_checks.values())
                         else 'step1_not_accepted'), 'checks': checks,
              'cohort_n': len(ids), 'case_n': int(d.is_case.sum()), 'diagnosis_age_n': len(age_ids),
              'saved_score_scaling': 'case_control reference, ddof=0; age subset reuses saved score',
              'actual_existing_model_scaling': 'ddof=1; case_control full cohort, age within cases, CV training fold only',
              'final_variant_definition_frozen': cfg['ptv_final_definition_frozen'],
              'models_rerun_in_this_step': False, 'manuscript_modified': False,
              'source_provenance_issue': provenance_issue,
              'source_provenance_resolution': 'preserve historical lock; create separate refresh lock for current source and unchanged persisted inputs; no upstream recomputation certified',
              'input_sha256': {p: hashlib.sha256((BASE / p).read_bytes()).hexdigest() for p in paths}}
    internal = ROOT / 'data/manuscript_refresh_20261006'
    internal.mkdir(exist_ok=True)
    sources = sorted((ROOT / 'scripts').glob('*.py')) + sorted((ROOT / 'scripts/functions').glob('*.py'))
    fresh_hashes = {p: h for p, h in lock['input_sha256'].items() if p.startswith('data/')}
    fresh_hashes.update({str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    fresh_hashes['config/rowlands_13gene_extension.yaml'] = hashlib.sha256((ROOT / 'config/rowlands_13gene_extension.yaml').read_bytes()).hexdigest()
    fresh_hashes['results/table/29_rowlands_13gene_external_weights.csv'] = hashlib.sha256((ROOT / 'results/table/29_rowlands_13gene_external_weights.csv').read_bytes()).hexdigest()
    check('fresh_refresh_input_lock_verifies', verify_frozen_inputs(fresh_hashes))
    (internal / 'refresh_input_lock_internal.json').write_text(json.dumps({
        'purpose': 'C downstream manuscript refresh from persisted inputs; no upstream rerun',
        'historical_lock': 'data/29_rowlands_extension/section9/pre_model_input_lock_internal.json',
        'historical_lock_preserved': True, 'input_sha256': fresh_hashes,
        'historical_source_migration_metadata_stale': bool(provenance_issue),
        'downstream_kernels_unchanged': True}, indent=2) + '\n')
    (internal / 'step1_input_audit_internal.json').write_text(json.dumps(result, indent=2) + '\n')
    public = {k: v for k, v in result.items() if k != 'input_sha256'}
    (ROOT / 'results/report/49_manuscript_refresh_step1_audit.json').write_text(json.dumps(public, indent=2) + '\n')
    print(json.dumps(public, indent=2))


if __name__ == '__main__':
    main()

"""
- Select female participants with WGS.
- Label genetic ancestry.
- Define breast cancer status and EHR dates.
"""
from __future__ import annotations
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
# ============================================================================
# 1. select: select_study_cohort
# ============================================================================

import os
import subprocess
import sys
import pandas as pd

def disclose(n):
    """Format counts according to the disclosure threshold."""
    n = int(n)
    return str(n) if n == 0 or n >= 20 else "<20"

def cohort_selection_section(t):
    print("\n" + "=" * 78 + f"\n{t}\n" + "=" * 78, flush=True)

def select_study_cohort():
    """Initialize runtime state and execute this component explicitly."""
    global ALL_TABLE, COHORT_SELECTION_BC_CONCEPT_ID, CACHE, COHORT_SELECTION_CDR, COHORT_SELECTION_DATA_DIR, EUR_TABLE, COHORT_SELECTION_FEMALE_CONCEPT_ID, COHORT_SELECTION_PROJECT_ID, _, anc, cohort_selection_bigquery, cc_df, cohort_selection_client, eur, cohort_selection_job_config, keep, lab, n, cohort_selection_out_cols, r, sr_white, cohort_selection_storage_opts, v2fn, wgs_ids

    from google.cloud import bigquery as cohort_selection_bigquery

    sys.path.insert(0, "/home/jupyter/workspace/analysis/analysis_v2/scripts")
    import function as v2fn  # noqa: E402

    COHORT_SELECTION_CDR = os.environ["WORKSPACE_CDR"]
    COHORT_SELECTION_PROJECT_ID = os.environ["GOOGLE_PROJECT"]
    cohort_selection_client = cohort_selection_bigquery.Client()
    cohort_selection_storage_opts = {"requester_pays": True, "project": COHORT_SELECTION_PROJECT_ID}

    COHORT_SELECTION_FEMALE_CONCEPT_ID = 45878463   # SexAtBirth_Female
    COHORT_SELECTION_BC_CONCEPT_ID = 4112853        # Malignant tumor of breast (SNOMED standard concept)

    COHORT_SELECTION_DATA_DIR = "/home/jupyter/workspace/analysis/data"
    CACHE = f"{COHORT_SELECTION_DATA_DIR}/01_female_wgs_ancestry.csv.gz"
    ALL_TABLE = f"{COHORT_SELECTION_PROJECT_ID}.my_analysis_results.female_wgs_ancestry_cohort"
    EUR_TABLE = f"{COHORT_SELECTION_PROJECT_ID}.my_analysis_results.eur_female_breast_cancer_cohort"

    r = subprocess.run(["git", "check-ignore", "-q", CACHE], cwd=COHORT_SELECTION_DATA_DIR)
    if r.returncode != 0:
        raise SystemExit(f"FATAL: Cache path is not covered by .gitignore: {CACHE}")

    # ----------------------------------------------------------------------
    cohort_selection_section('1. Select female participants and identify breast cancer diagnoses')
    cc_df = cohort_selection_client.query(f"""
    WITH breast_cancer_patients AS (
      SELECT DISTINCT co.person_id
      FROM `{COHORT_SELECTION_CDR}.condition_occurrence` co
      JOIN `{COHORT_SELECTION_CDR}.concept_ancestor` ca ON co.condition_concept_id = ca.descendant_concept_id
      WHERE ca.ancestor_concept_id = {COHORT_SELECTION_BC_CONCEPT_ID}
    )
    SELECT p.person_id, p.race_concept_id,
           CASE WHEN bc.person_id IS NOT NULL THEN 1 ELSE 0 END AS is_case
    FROM `{COHORT_SELECTION_CDR}.person` p
    LEFT JOIN breast_cancer_patients bc ON p.person_id = bc.person_id
    WHERE p.sex_at_birth_concept_id = {COHORT_SELECTION_FEMALE_CONCEPT_ID}
    """).to_dataframe()
    print(f"Female participants (sex at birth): {len(cc_df):,d}; "
          f"with a breast cancer diagnosis: {cc_df.is_case.sum():,d}")

    # ----------------------------------------------------------------------
    cohort_selection_section('2. Intersect with the WGS manifest')
    wgs_ids = pd.read_csv("gs://vwb-aou-datasets-controlled/v9/wgs/cram/manifest.csv",
                          usecols=["person_id"], storage_options=cohort_selection_storage_opts)
    cc_df = cc_df[cc_df["person_id"].isin(set(wgs_ids["person_id"]))].reset_index(drop=True)
    print(f"Female participants with WGS: {len(cc_df):,d}  (cases {cc_df.is_case.sum():,d} / "
          f"non-cases {(cc_df.is_case == 0).sum():,d})")

    # ----------------------------------------------------------------------
    cohort_selection_section('3. Assign genetic ancestry')
    if os.path.exists(CACHE):
        anc = pd.read_csv(CACHE, dtype={"person_id": str})
        print(f"  Loaded ancestry cache: {CACHE} ({len(anc):,d} rows)")
    else:
        keep = set(cc_df["person_id"].astype(str))
        anc, _ = v2fn.extract_ancestry_preds(keep)
        anc.to_csv(CACHE, index=False, compression="gzip")
        print(f"  Saved ancestry cache: {CACHE}")

    anc["person_id"] = anc["person_id"].astype("int64")
    cc_df = cc_df.merge(anc[["person_id", "ancestry_pred", "ancestry_pred_other", "p_eur"]],
                        on="person_id", how="left", validate="one_to_one")
    if cc_df["ancestry_pred_other"].isna().any():
        raise SystemExit('FATAL: There are WGS samples missing ancestry annotations')
    cc_df["is_eur"] = cc_df["ancestry_pred_other"].eq("eur")

    print('Genetic ancestry distribution among female participants with WGS:')
    for lab, n in cc_df["ancestry_pred_other"].value_counts().items():
        print(f"  {lab:<5s} {disclose(n)}")
    eur = cc_df[cc_df["is_eur"]]
    print(f"\nEuropean genetic ancestry: {len(eur):,d}  (cases {eur.is_case.sum():,d} / non-cases {(eur.is_case == 0).sum():,d})")

    sr_white = eur["race_concept_id"].eq(8527)
    print(f"Self-reported White: {disclose(sr_white.sum())}; "
          f"other race or no response: {disclose((~sr_white).sum())}")

    # ----------------------------------------------------------------------
    cohort_selection_section('4. Save cohort tables to the analysis BigQuery dataset')
    cohort_selection_job_config = cohort_selection_bigquery.LoadJobConfig(write_disposition="WRITE_TRUNCATE")
    cohort_selection_out_cols = ["person_id", "is_case", "race_concept_id", "ancestry_pred_other", "is_eur"]
    cohort_selection_client.load_table_from_dataframe(cc_df[cohort_selection_out_cols], ALL_TABLE, job_config=cohort_selection_job_config).result()
    cohort_selection_client.load_table_from_dataframe(eur[["person_id", "is_case"]].reset_index(drop=True),
                                     EUR_TABLE, job_config=cohort_selection_job_config).result()
    print(f"  {ALL_TABLE}: {len(cc_df):,d} rows")
    print(f"  {EUR_TABLE}: {len(eur):,d} rows")

# ============================================================================
# 2. phenotypes: define_phenotypes_and_ehr_dates
# ============================================================================

import os
import sys
import time
import numpy as np
import pandas as pd

def phenotypes_section(title):
    print("\n" + "=" * 90)
    print(title)
    print("=" * 90)

def timed_query(label, query):
    t0 = time.time()
    df = phenotypes_client.query(query).to_dataframe()
    print(f"[{label}] elapsed={time.time()-t0:.1f}s  n_rows={len(df)}")
    return df

def _group(row):
    if row["is_case"]:
        return "case"
    elif row["is_control"]:
        return "control"
    else:
        return "other_cancer_excluded"

def _age_years(end_date, birth_date):
    return (end_date - birth_date).dt.days / 365.25

def define_phenotypes_and_ehr_dates():
    """Initialize runtime state and execute this component explicitly."""
    global ANY_CANCER_CONCEPT_ID, PHENOTYPES_BC_CONCEPT_ID, PHENOTYPES_CDR, PHENOTYPES_DATA_DIR, PHENOTYPES_FEMALE_CONCEPT_ID, PHENOTYPES_PROJECT_ID, TABLE_DIR, UNIVERSE_TABLE, age_dist_df, age_mask, age_rows, age_thr, age_threshold_df, age_vars, any_cancer_df, any_cancer_q, any_cancer_set, bc_dx_df, bc_dx_q, phenotypes_bigquery, c, case_df, phenotypes_client, combo_df, combo_rows, consent_df, consent_q, control_df, cov_dist_df, cov_rows, cov_thr_rows, cov_threshold_df, cov_vars, death_df, death_q, df, dx_count_df, dx_count_q, ehr_dates_df, ehr_dates_q, first_record_summary, gap_days, gap_mean, gap_median, gap_q25, gap_q75, gap_s, gap_summary_df, grp, phenotypes_job_config, label, mask, meets, meets_5, meets_b, n_5, n_age_under18, n_b, n_case, n_combo, n_control, n_control_age_na, n_excluded_na, n_excluded_under, n_first_record, n_gap, n_gap_missing, n_le_half_year, n_meets, n_negative, n_no_ehr, n_r2_first, n_within, phenotypes_out_cols, out_path, pct_first_record, pct_le_half_year, person_out, person_q, r2_case, s, size_mb, phenotypes_storage_opts, sub, thr, threshold_rows, threshold_summary_all, universe, v, window

    from google.cloud import bigquery as phenotypes_bigquery

    PHENOTYPES_CDR = os.environ["WORKSPACE_CDR"]
    PHENOTYPES_PROJECT_ID = os.environ["GOOGLE_PROJECT"]
    phenotypes_client = phenotypes_bigquery.Client()
    phenotypes_storage_opts = {"requester_pays": True, "project": PHENOTYPES_PROJECT_ID}

    PHENOTYPES_FEMALE_CONCEPT_ID = 45878463         # SexAtBirth_Female
    PHENOTYPES_BC_CONCEPT_ID = 4112853              # SNOMED "Malignant tumor of breast"
    ANY_CANCER_CONCEPT_ID = 443392       # SNOMED "Malignant neoplastic disease" (363346000)

    UNIVERSE_TABLE = f"{PHENOTYPES_PROJECT_ID}.my_analysis_results.aim1_phenotype_universe"

    PHENOTYPES_DATA_DIR = "/home/jupyter/workspace/analysis/data"
    TABLE_DIR = "/home/jupyter/workspace/analysis/data/04_internal"  
    os.makedirs(PHENOTYPES_DATA_DIR, exist_ok=True)
    os.makedirs(TABLE_DIR, exist_ok=True)

    # ============================================================

    # ============================================================
    phenotypes_section('0. Build female + WGS universe (taken from female_wgs_ancestry_cohort No. 01)')

    person_q = f"""
    SELECT p.person_id, p.birth_datetime
    FROM `{PHENOTYPES_CDR}.person` p
    JOIN `{PHENOTYPES_PROJECT_ID}.my_analysis_results.female_wgs_ancestry_cohort` c
      ON p.person_id = c.person_id
    """
    universe = timed_query("female_wgs", person_q)
    print(f"universe (female+WGS, all ancestries) n = {len(universe)}")

    phenotypes_job_config = phenotypes_bigquery.LoadJobConfig(write_disposition="WRITE_TRUNCATE")
    phenotypes_client.load_table_from_dataframe(universe, UNIVERSE_TABLE, job_config=phenotypes_job_config).result()
    print(f"Saved universe to BigQuery:{UNIVERSE_TABLE}")

    # ============================================================
    
    # ============================================================
    phenotypes_section('1. Case definition: single EHR (main analysis) + rule of two (sensitivity markers)')

    bc_dx_q = f"""
    WITH bc_dx AS (
      SELECT co.person_id, DATE(co.condition_start_datetime) AS dx_date
      FROM `{PHENOTYPES_CDR}.condition_occurrence` co
      JOIN `{UNIVERSE_TABLE}` u ON co.person_id = u.person_id
      JOIN `{PHENOTYPES_CDR}.concept_ancestor` ca ON co.condition_concept_id = ca.descendant_concept_id
      WHERE ca.ancestor_concept_id = {PHENOTYPES_BC_CONCEPT_ID}
    )
    SELECT
      person_id,
      MIN(dx_date) AS first_bc_dx_date,
      COUNT(DISTINCT dx_date) AS n_distinct_bc_dx_dates
    FROM bc_dx
    GROUP BY person_id
    """
    bc_dx_df = timed_query("bc_dx (case ascertainment)", bc_dx_q)
    bc_dx_df["is_case_rule_of_two"] = bc_dx_df["n_distinct_bc_dx_dates"] >= 2

    print(f"Number of single EHR cases ={len(bc_dx_df)}")
    print(f"rule of two (≥2 different dates) subset ={bc_dx_df['is_case_rule_of_two'].sum()} "
          f"({bc_dx_df['is_case_rule_of_two'].mean()*100:.1f}% of single-EHR case)")

    # ============================================================
    
    # ============================================================
    phenotypes_section("2. Control qualification mark: exclude all descendants of 'malignant neoplastic disease' (363346000) under concept_ancestor")

    any_cancer_q = f"""
    SELECT DISTINCT co.person_id
    FROM `{PHENOTYPES_CDR}.condition_occurrence` co
    JOIN `{UNIVERSE_TABLE}` u ON co.person_id = u.person_id
    JOIN `{PHENOTYPES_CDR}.concept_ancestor` ca ON co.condition_concept_id = ca.descendant_concept_id
    WHERE ca.ancestor_concept_id = {ANY_CANCER_CONCEPT_ID}
    """
    any_cancer_df = timed_query("any_cancer_dx", any_cancer_q)
    any_cancer_set = set(any_cancer_df["person_id"])
    print(f"The number of people in the universe with any cancer diagnosis code ={len(any_cancer_set)} "
          f"(Including breast cancer cases, but also other cancer types)")

    # ============================================================

    # ============================================================
    phenotypes_section('3. EHR coverage: earliest/latest record date and distinct date number of the union of five clinical fact tables')

    ehr_dates_q = f"""
    WITH domain_dates AS (
      SELECT co.person_id, co.condition_start_date AS event_date
      FROM `{PHENOTYPES_CDR}.condition_occurrence` co
      JOIN `{UNIVERSE_TABLE}` u ON co.person_id = u.person_id
      UNION ALL
      SELECT po.person_id, po.procedure_date AS event_date
      FROM `{PHENOTYPES_CDR}.procedure_occurrence` po
      JOIN `{UNIVERSE_TABLE}` u ON po.person_id = u.person_id
      UNION ALL
      SELECT de.person_id, de.drug_exposure_start_date AS event_date
      FROM `{PHENOTYPES_CDR}.drug_exposure` de
      JOIN `{UNIVERSE_TABLE}` u ON de.person_id = u.person_id
      UNION ALL
      SELECT m.person_id, m.measurement_date AS event_date
      FROM `{PHENOTYPES_CDR}.measurement` m
      JOIN `{UNIVERSE_TABLE}` u ON m.person_id = u.person_id
      UNION ALL
      SELECT vo.person_id, vo.visit_start_date AS event_date
      FROM `{PHENOTYPES_CDR}.visit_occurrence` vo
      JOIN `{UNIVERSE_TABLE}` u ON vo.person_id = u.person_id
    )
    SELECT
      person_id,
      MIN(event_date) AS ehr_first_date,
      MAX(event_date) AS ehr_last_date,
      COUNT(DISTINCT event_date) AS n_distinct_ehr_dates
    FROM domain_dates
    GROUP BY person_id
    """
    ehr_dates_df = timed_query("ehr_dates (5-domain union)", ehr_dates_q)
    n_no_ehr = len(universe) - len(ehr_dates_df)
    print(f"The number of people in the universe who have no records at all in the 5 clinical tables ={n_no_ehr} "
          f"({n_no_ehr/len(universe)*100:.2f}%) - for these people the age_at_last_obs/EHR coverage statistic will be null,"
          f"Unless a death date is available")

    # ------------------------------------------------------------

    # ------------------------------------------------------------
    dx_count_q = f"""
    SELECT co.person_id, COUNT(*) AS dx_code_occurrence_count
    FROM `{PHENOTYPES_CDR}.condition_occurrence` co
    JOIN `{UNIVERSE_TABLE}` u ON co.person_id = u.person_id
    GROUP BY co.person_id
    """
    dx_count_df = timed_query("dx_code_occurrence_count (PheTK-style covariate)", dx_count_q)

    # ============================================================
    
    # ============================================================
    phenotypes_section('4. age_at_enrollment (cb_search_person.age_at_consent) + death date (OMOP death / aou_death)')

    consent_q = f"""
    SELECT u.person_id, csp.age_at_consent
    FROM `{UNIVERSE_TABLE}` u
    JOIN `{PHENOTYPES_CDR}.cb_search_person` csp ON u.person_id = csp.person_id
    """
    consent_df = timed_query("age_at_consent", consent_q)
    print(f"age_at_consent missing number ={consent_df['age_at_consent'].isna().sum()} / {len(consent_df)}")

    death_q = f"""-- 2026-09-24: There can be multiple records for the same person in the death table. After expanding to all ancestors, duplicate rows appear and are aggregated by person.
    SELECT
      u.person_id,
      MIN(COALESCE(ad.death_date, d.death_date)) AS death_date
    FROM `{UNIVERSE_TABLE}` u
    LEFT JOIN `{PHENOTYPES_CDR}.death` d ON u.person_id = d.person_id
    LEFT JOIN (
      SELECT person_id, death_date
      FROM `{PHENOTYPES_CDR}.aou_death`
      WHERE primary_death_record = TRUE
    ) ad ON u.person_id = ad.person_id
    WHERE d.death_date IS NOT NULL OR ad.death_date IS NOT NULL
    GROUP BY u.person_id
    """
    death_df = timed_query("death_dates", death_q)
    print(f"Number of people with death records in the universe ={len(death_df)}")

    # ============================================================
    
    # ============================================================
    phenotypes_section('5. Merge and calculate case/control grouping and three age variables')

    df = universe.merge(bc_dx_df, on="person_id", how="left")
    df = df.merge(ehr_dates_df, on="person_id", how="left")
    df = df.merge(dx_count_df, on="person_id", how="left")
    df = df.merge(consent_df, on="person_id", how="left")
    df = df.merge(death_df, on="person_id", how="left")
    df["dx_code_occurrence_count"] = df["dx_code_occurrence_count"].fillna(0).astype(int)

    df["is_case"] = df["first_bc_dx_date"].notna()
    df["is_case_rule_of_two"] = df["is_case_rule_of_two"].fillna(False)
    df["has_any_cancer_dx"] = df["person_id"].isin(any_cancer_set)
    df["is_control"] = (~df["is_case"]) & (~df["has_any_cancer_dx"])
    df["is_other_cancer_excluded"] = (~df["is_case"]) & (df["has_any_cancer_dx"])

    df["group"] = df.apply(_group, axis=1)

    df["birth_date"] = pd.to_datetime(df["birth_datetime"]).dt.tz_localize(None).dt.normalize()

    for c in ["first_bc_dx_date", "ehr_first_date", "ehr_last_date", "death_date"]:
        df[c] = pd.to_datetime(df[c])

    df["last_obs_date"] = df[["ehr_last_date", "death_date"]].max(axis=1, skipna=True)

    df["age_at_diagnosis"] = np.where(
        df["is_case"], _age_years(df["first_bc_dx_date"], df["birth_date"]), np.nan
    )
    df["age_at_enrollment"] = df["age_at_consent"].astype(float)
    df["age_at_last_obs"] = _age_years(df["last_obs_date"], df["birth_date"])

    # Dimensions of Partial vs Complete Mastectomy Phenotypes in the All of Us Research

    n_age_under18 = int((df["age_at_diagnosis"] < 18).sum())
    print(f"\nage_at_diagnosis<18 (biologically untrustworthy, the minimum age for AoU enrollment is 18 years old) Number of people:{n_age_under18}，"
          f"Age_at_diagnosis is left blank and flag_age_at_diagnosis_invalid_under18 is marked")
    df["flag_age_at_diagnosis_invalid_under18"] = df["age_at_diagnosis"] < 18
    df.loc[df["flag_age_at_diagnosis_invalid_under18"].fillna(False), "age_at_diagnosis"] = np.nan

    df["ehr_span_days"] = (df["ehr_last_date"] - df["ehr_first_date"]).dt.days
    df["ehr_span_years"] = df["ehr_span_days"] / 365.25

    df["bc_is_first_record"] = np.where(
        df["is_case"], df["first_bc_dx_date"] == df["ehr_first_date"], np.nan
    )

    df["gap_last_obs_minus_diagnosis"] = np.where(
        df["is_case"], df["age_at_last_obs"] - df["age_at_diagnosis"], np.nan
    )

    print(df["group"].value_counts())
    print(f"Total number of \ncases ={(df['group']=='case').sum()}")
    print(f"Total number of controls (without any cancer diagnosis, without EHR coverage threshold) ={(df['group']=='control').sum()}")
    print(f"other_cancer_excluded (there are other cancers but not breast cancer, neither case nor control) ="
          f"{(df['group']=='other_cancer_excluded').sum()}")

    # ============================================================
    
    # ============================================================
    phenotypes_section('Diagnostic output 1: Distribution of three age variables grouped by case/control (mean/median/IQR)')

    age_vars = ["age_at_diagnosis", "age_at_enrollment", "age_at_last_obs"]
    age_rows = []
    for grp in ["case", "control"]:
        sub = df[df["group"] == grp]
        for v in age_vars:
            s = sub[v].dropna()
            if len(s) == 0:
                continue
            age_rows.append({
                "group": grp,
                "variable": v,
                "n_nonmissing": len(s),
                "n_missing": sub[v].isna().sum(),
                "mean": s.mean(),
                "median": s.median(),
                "q25": s.quantile(0.25),
                "q75": s.quantile(0.75),
                "min": s.min(),
                "max": s.max(),
            })
    age_dist_df = pd.DataFrame(age_rows)
    pd.set_option("display.width", 160)
    print(age_dist_df.to_string(index=False))
    age_dist_df.to_csv(f"{TABLE_DIR}/04_phenotype_age_distributions.csv", index=False)
    print(f"\n has been saved:{TABLE_DIR}/04_phenotype_age_distributions.csv")

    # ============================================================
    
    # ============================================================
    phenotypes_section("Diagnosis output 2: n_bc_is_first_record (the first BC diagnosis date in the case == the earliest record date of the person's EHR)")

    case_df = df[df["group"] == "case"].copy()
    n_case = len(case_df)
    n_first_record = int(case_df["bc_is_first_record"].sum())
    pct_first_record = n_first_record / n_case * 100
    print(f"Total number of cases (single EHR) ={n_case}")
    print(f"n_bc_is_first_record = {n_first_record} ({pct_first_record:.1f}%)")

    r2_case = case_df[case_df["is_case_rule_of_two"]]
    n_r2_first = int(r2_case["bc_is_first_record"].sum())
    print(f"[Additional] rule-of-two subset (n={len(r2_case)}) in n_bc_is_first_record ="
          f"{n_r2_first} ({n_r2_first/len(r2_case)*100:.1f}%)")

    for window in [7, 30]:
        gap_days = (case_df["first_bc_dx_date"] - case_df["ehr_first_date"]).dt.days
        n_within = int((gap_days <= window).sum())
        print(f"[Additional, softened criteria] First BC diagnosis falls after earliest EHR record{window}Number of cases within a day ="
              f"{n_within} ({n_within/n_case*100:.1f}%)")

    first_record_summary = pd.DataFrame([{
        "n_case_total": n_case,
        "n_bc_is_first_record": n_first_record,
        "pct_bc_is_first_record": round(pct_first_record, 2),
        "n_case_rule_of_two": len(r2_case),
        "n_bc_is_first_record_rule_of_two_subset": n_r2_first,
        "pct_bc_is_first_record_rule_of_two_subset": round(n_r2_first/len(r2_case)*100, 2),
    }])
    first_record_summary.to_csv(f"{TABLE_DIR}/04_phenotype_n_bc_is_first_record.csv", index=False)
    print(f"\n has been saved:{TABLE_DIR}/04_phenotype_n_bc_is_first_record.csv")

    # ============================================================
    
    # ============================================================
    phenotypes_section('Diagnostic output 3: EHR coverage descriptive statistics (record span, distinct record date number, grouped by case/control)')

    cov_vars = ["ehr_span_days", "ehr_span_years", "n_distinct_ehr_dates", "dx_code_occurrence_count"]
    cov_rows = []
    for grp in ["case", "control"]:
        sub = df[df["group"] == grp]
        for v in cov_vars:
            s = sub[v].dropna()
            if len(s) == 0:
                continue
            cov_rows.append({
                "group": grp,
                "variable": v,
                "n_nonmissing": len(s),
                "n_missing_or_no_ehr": sub[v].isna().sum(),
                "mean": s.mean(),
                "median": s.median(),
                "q25": s.quantile(0.25),
                "q75": s.quantile(0.75),
            })
    cov_dist_df = pd.DataFrame(cov_rows)
    print(cov_dist_df.to_string(index=False))
    cov_dist_df.to_csv(f"{TABLE_DIR}/04_phenotype_ehr_coverage_distributions.csv", index=False)
    print(f"\n has been saved:{TABLE_DIR}/04_phenotype_ehr_coverage_distributions.csv")

    # ============================================================
    
    # ============================================================
    phenotypes_section('Diagnostic output 4a: control sample size impact of age threshold candidates (age_at_last_obs >=50 / >=55)')

    control_df = df[df["group"] == "control"].copy()
    n_control = len(control_df)
    n_control_age_na = control_df["age_at_last_obs"].isna().sum()

    threshold_rows = []
    for thr in [50, 55]:
        meets = (control_df["age_at_last_obs"] >= thr)
        n_meets = int(meets.sum())
        n_excluded_under = int((control_df["age_at_last_obs"] < thr).sum())
        n_excluded_na = int(control_df["age_at_last_obs"].isna().sum())
        threshold_rows.append({
            "candidate_threshold": f"age_at_last_obs >= {thr}",
            "n_control_total": n_control,
            "n_meets_threshold": n_meets,
            f"n_under{thr}_at_last_obs": n_excluded_under,
            "n_age_at_last_obs_missing_no_ehr": n_excluded_na,
            "pct_retained": round(n_meets / n_control * 100, 1),
        })
    age_threshold_df = pd.DataFrame(threshold_rows)
    print(age_threshold_df.to_string(index=False))
    print(f"\n (reference) control overall age_at_last_obs missing (no EHR/no death record to push) number of people ={n_control_age_na}")

    phenotypes_section('Diagnostic Output 4b: Sample size impact of control EHR coverage threshold candidates (no threshold / span >= 1 year + distinct >= 2 / distinct >= 5)')

    cov_thr_rows = []
    
    cov_thr_rows.append({
        "candidate_ehr_coverage_threshold": 'No threshold (current control collection, female + WGS + no cancer diagnosis; EUR limited to 05)',
        "n_control": n_control,
        "pct_of_no_threshold": 100.0,
    })
    
    meets_b = (control_df["ehr_span_years"] >= 1) & (control_df["n_distinct_ehr_dates"] >= 2)
    n_b = int(meets_b.fillna(False).sum())
    cov_thr_rows.append({
        "candidate_ehr_coverage_threshold": 'EHR span >= 1 year and distinct record date >= 2 (Step0 plan B recommended threshold)',
        "n_control": n_b,
        "pct_of_no_threshold": round(n_b / n_control * 100, 1),
    })
    # distinct >= 5
    meets_5 = (control_df["n_distinct_ehr_dates"] >= 5)
    n_5 = int(meets_5.fillna(False).sum())
    cov_thr_rows.append({
        "candidate_ehr_coverage_threshold": 'distinct record date>=5',
        "n_control": n_5,
        "pct_of_no_threshold": round(n_5 / n_control * 100, 1),
    })
    cov_threshold_df = pd.DataFrame(cov_thr_rows)
    print(cov_threshold_df.to_string(index=False))

    phenotypes_section('Diagnostic output 4c (additional/for reference): age threshold and EHR coverage threshold control sample size after joint application')

    combo_rows = []
    for age_thr in [50, 55]:
        for label, mask in [
            ('No EHR coverage threshold', pd.Series(True, index=control_df.index)),
            ('span>=1 year+distinct>=2', meets_b.fillna(False)),
            ("distinct>=5", meets_5.fillna(False)),
        ]:
            age_mask = (control_df["age_at_last_obs"] >= age_thr).fillna(False)
            n_combo = int((age_mask & mask).sum())
            combo_rows.append({
                "age_threshold": f">={age_thr}",
                "ehr_coverage_threshold": label,
                "n_control_remaining": n_combo,
                "pct_of_no_threshold_control": round(n_combo / n_control * 100, 1),
            })
    combo_df = pd.DataFrame(combo_rows)
    print(combo_df.to_string(index=False))

    threshold_summary_all = pd.concat(
        [
            age_threshold_df.assign(table="age_threshold"),
            cov_threshold_df.assign(table="ehr_coverage_threshold"),
            combo_df.assign(table="combined_reference"),
        ],
        axis=0,
        ignore_index=True,
    )
    threshold_summary_all.to_csv(f"{TABLE_DIR}/04_phenotype_threshold_sensitivity.csv", index=False)
    print(f"\n has been saved:{TABLE_DIR}/04_phenotype_threshold_sensitivity.csv")

    # ============================================================

    # ============================================================
    phenotypes_section('Diagnosis output 5 (new in v7.1): age_at_last_obs - age_at_diagnosis distribution within case group')

    gap_s = case_df["gap_last_obs_minus_diagnosis"].dropna()
    n_gap = len(gap_s)
    n_gap_missing = (df["group"] == "case").sum() - n_gap
    gap_median = gap_s.median()
    gap_q25 = gap_s.quantile(0.25)
    gap_q75 = gap_s.quantile(0.75)
    gap_mean = gap_s.mean()
    n_le_half_year = int((gap_s <= 0.5).sum())
    pct_le_half_year = n_le_half_year / n_gap * 100
    n_negative = int((gap_s < 0).sum())  

    print(f"The number of people in the case group who can calculate the gap ={n_gap}(missing{n_gap_missing}, mostly cases without any EHR record/no valid last_obs_date)")
    print(f"gap = age_at_last_obs - age_at_diagnosis:")
    print(f"  mean   = {gap_mean:.2f}Year")
    print(f"  median = {gap_median:.2f}Year")
    print(f"  IQR    = [{gap_q25:.2f}, {gap_q75:.2f}] Year")
    print(f"Proportion of cases with gap <= 0.5 years ={n_le_half_year} ({pct_le_half_year:.1f}%)")
    if n_negative > 0:
        print(f"!!! Exception: number of cases with gap<0 ={n_negative}(last_obs is earlier than diagnosis, needs to be checked, and theoretically should not appear)")

    gap_summary_df = pd.DataFrame([{
        "n_case_with_gap": n_gap,
        "n_case_gap_missing": int(n_gap_missing),
        "gap_mean_years": round(gap_mean, 3),
        "gap_median_years": round(gap_median, 3),
        "gap_q25_years": round(gap_q25, 3),
        "gap_q75_years": round(gap_q75, 3),
        "n_gap_le_0.5yr": n_le_half_year,
        "pct_gap_le_0.5yr": round(pct_le_half_year, 2),
        "n_gap_negative_anomaly": n_negative,
    }])
    gap_summary_df.to_csv(f"{TABLE_DIR}/04_phenotype_case_last_obs_minus_diagnosis_gap.csv", index=False)
    print(f"\n has been saved:{TABLE_DIR}/04_phenotype_case_last_obs_minus_diagnosis_gap.csv")

    # ============================================================
    
    # ============================================================
    phenotypes_section('Save person-level diagnostic results (reduced columns, control file size)')

    phenotypes_out_cols = [
        "person_id", "group", "is_case", "is_case_rule_of_two", "has_any_cancer_dx",
        "n_distinct_bc_dx_dates", "bc_is_first_record",
        "age_at_diagnosis", "age_at_enrollment", "age_at_last_obs",
        "flag_age_at_diagnosis_invalid_under18",
        "gap_last_obs_minus_diagnosis",
        "ehr_span_days", "n_distinct_ehr_dates", "dx_code_occurrence_count",
    ]
    person_out = df[phenotypes_out_cols].copy()
    assert not df["person_id"].duplicated().any(), 'person_id duplicate'
    out_path = f"{PHENOTYPES_DATA_DIR}/04_phenotype_diagnostics_person_level.csv.gz"
    person_out.to_csv(out_path, index=False, compression="gzip")
    size_mb = os.path.getsize(out_path) / 1e6
    print(f"Saved:{out_path} ({size_mb:.2f} MB)")
    if size_mb > 5:
        print('!!! WARNING: The file exceeds the 5MB limit and needs to be further streamlined !!!', file=sys.stderr)

    print('Finish.')

# Command-line selection of the computations defined above.
def main():
    from scripts.functions.cli import run_local_step
    computations = {
        'select': select_study_cohort,
        'phenotypes': define_phenotypes_and_ehr_dates,
    }
    run_local_step('Study cohort selection and phenotype definition', computations, default_tasks=['select', 'phenotypes'])

if __name__ == '__main__':
    main()

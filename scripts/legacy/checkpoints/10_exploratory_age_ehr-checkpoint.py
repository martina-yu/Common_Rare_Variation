"""
- Run optional age and PRS contrasts.
- Audit early EHR diagnoses.
- Reproduce exploratory missense analyses.
"""
from __future__ import annotations
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
# ============================================================================
# 1. age-contrasts: analyze_age_group_contrasts
# ============================================================================
from pathlib import Path
import json
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import norm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

AGE_GROUP_CONTRASTS_ROOT = Path(__file__).resolve().parents[1]
AGE_GROUP_CONTRASTS_PCS = [f"PC{age_group_contrasts_i}" for age_group_contrasts_i in range(1, 11)]
TERMS = ["const", "PRS_z", "carrier", "interaction"] + AGE_GROUP_CONTRASTS_PCS

def age_group_contrasts_read(name):
    return pd.read_csv(AGE_GROUP_CONTRASTS_ROOT / "data" / name)

def contrast(fit, vector):
    v = np.asarray(vector, dtype=float)
    estimate = float(v @ fit.params)
    se = float(np.sqrt(max(0, v @ fit.cov_params() @ v)))
    return {"estimate_years": estimate, "se_HC1": se,
            "ci_low": estimate - norm.ppf(.975) * se,
            "ci_high": estimate + norm.ppf(.975) * se}

def design(p, c, pc):
    return np.array([1, p, c, p*c, *pc], dtype=float)

def analyze(df, panel, label, reference, out):
    mu, sd, cutoff, p_low, p_high, pc_ref = reference
    df = df.copy()
    df["carrier"] = (df[f"burden_rare_main_ptv_{panel}gene"] > 0).astype(int)
    df["high"] = (df.PRS311 >= cutoff).astype(int)
    df["PRS_z"] = (df.PRS311 - mu) / sd
    df["interaction"] = df.PRS_z * df.carrier
    counts = df.groupby(["carrier", "high"]).size().reindex(
        pd.MultiIndex.from_product([[0, 1], [0, 1]]), fill_value=0)
    audit = {"analysis": label, "four_group_release": bool((counts >= 20).all())}
    # Suppress the entire table if a cell is small; do not emit recoverable margins.
    if (df.groupby("carrier").size().reindex([0, 1], fill_value=0) < 20).any():
        audit["continuous_model"] = "suppressed"
        return audit
    X = sm.add_constant(df[["PRS_z", "carrier", "interaction"] + AGE_GROUP_CONTRASTS_PCS], has_constant="add")
    fit = sm.OLS(df.age_at_diagnosis.astype(float), X).fit(cov_type="HC1")
    rows = []
    for term in ["PRS_z", "carrier", "interaction"]:
        v = np.zeros(len(TERMS)); v[TERMS.index(term)] = 1
        rows.append({"contrast": term, **contrast(fit, v)})
    v = np.zeros(len(TERMS)); v[1] = v[3] = 1
    rows.append({"contrast": "PRS_per_SD_carriers", **contrast(fit, v)})
    pd.DataFrame(rows).to_csv(out / f"{label}_continuous.csv", index=False)
    if audit["four_group_release"]:
        vectors, summaries = {}, []
        for c in [0, 1]:
            for h in [0, 1]:
                g = df[(df.carrier == c) & (df.high == h)].age_at_diagnosis
                v = design(p_high if h else p_low, c, pc_ref)
                vectors[c, h] = v
                summaries.append({"carrier": c, "high_PRS": h, "n": len(g),
                    "observed_mean": g.mean(), "observed_sd": g.std(),
                    "median": g.median(), "q25": g.quantile(.25), "q75": g.quantile(.75),
                    **contrast(fit, v)})
        pd.DataFrame(summaries).to_csv(out / f"{label}_four_groups.csv", index=False)
        comparisons = {"high_vs_nonhigh_noncarriers": vectors[0,1]-vectors[0,0],
            "high_vs_nonhigh_carriers": vectors[1,1]-vectors[1,0],
            "difference_in_differences": vectors[1,1]-vectors[1,0]-vectors[0,1]+vectors[0,0],
            "double_high_vs_noncarrier_nonhigh": vectors[1,1]-vectors[0,0]}
        pd.DataFrame([{"contrast": k, **contrast(fit, v)} for k,v in comparisons.items()]).to_csv(
            out / f"{label}_contrasts.csv", index=False)
    # Only aggregate model curves; common support avoids carrier-tail extrapolation.
    lo = max(df.loc[df.carrier == c, "PRS_z"].quantile(.01) for c in [0,1])
    hi = min(df.loc[df.carrier == c, "PRS_z"].quantile(.99) for c in [0,1])
    if lo < hi:
        fig, ax = plt.subplots()
        for c in [0,1]:
            xs = np.linspace(lo, hi, 100)
            predictions = pd.DataFrame([contrast(fit, design(p, c, pc_ref)) for p in xs])
            line, = ax.plot(xs, predictions.estimate_years, label=f"Carrier={c}")
            ax.fill_between(xs, predictions.ci_low, predictions.ci_high, alpha=.2, color=line.get_color())
        ax.set(xlabel="PRS311 (SD in primary cases)", ylabel="Adjusted recorded diagnosis age (years)")
        ax.legend(); fig.tight_layout()
        fig.savefig(out / f"{label}_curves.svg"); plt.close(fig)
    event_counts = df.groupby(["carrier", "high", df.age_at_diagnosis < 45]).size().reindex(
        pd.MultiIndex.from_product([[0,1], [0,1], [False,True]]), fill_value=0)
    audit["case_only_early_cells_release_eligible"] = bool((event_counts >= 20).all())
    audit["population_early_risk"] = "not_estimated_observation_window_requires_review"
    return audit

def analyze_age_group_contrasts():
    pheno = age_group_contrasts_read("04_phenotype_diagnostics_person_level.csv.gz")
    qc = age_group_contrasts_read("05_cohort_qc_flags.csv")
    prs = age_group_contrasts_read("prs_scores.csv")
    burden = age_group_contrasts_read("14_burden_rare_scores_person_level.csv.gz")
    pc = age_group_contrasts_read("09_ancestry_pc_person_level.csv.gz")
    membership = age_group_contrasts_read("16_joint_model_cv_oof_person_level.csv.gz")[["person_id"]].drop_duplicates()
    full = membership.merge(prs[["person_id", "PRS311"]], validate="one_to_one")
    assert len(full) == len(membership) and full.PRS311.notna().all()
    cutoff = full.PRS311.quantile(.9)
    base = pheno.loc[pheno.group == "case"].merge(qc[["person_id", "pass_qc"]], validate="one_to_one")
    base = base.loc[base.pass_qc == 1].dropna(subset=["age_at_diagnosis"])
    for source in [prs[["person_id", "PRS311"]], burden[["person_id", "burden_rare_main_ptv_10gene", "burden_rare_main_ptv_9gene"]], pc[["person_id"] + AGE_GROUP_CONTRASTS_PCS]]:
        base = base.merge(source, how="left", validate="one_to_one")
    assert base[["PRS311", "burden_rare_main_ptv_10gene", "burden_rare_main_ptv_9gene"]+AGE_GROUP_CONTRASTS_PCS].notna().all().all()
    assert base.bc_is_first_record.notna().all()
    primary = base.loc[base.bc_is_first_record != True].copy()
    assert len(primary) == 7996, "Primary cohort differs from fixed script-17 cohort; stop and review."
    mu, sd = primary.PRS311.mean(), primary.PRS311.std(ddof=1)
    z = (primary.PRS311-mu)/sd
    reference = (mu, sd, cutoff, z[primary.PRS311 < cutoff].mean(),
                 z[primary.PRS311 >= cutoff].mean(), primary[AGE_GROUP_CONTRASTS_PCS].mean().to_numpy())
    assert np.isfinite(np.array(reference[:5])).all() and sd > 0
    out = AGE_GROUP_CONTRASTS_ROOT / "results" / "common_rare_age"
    out.mkdir(parents=True, exist_ok=True)
    jobs = [(primary, 10, "primary"), (base, 10, "include_first_record"),
            (primary, 9, "nine_gene")]
    if "is_case_rule_of_two" in primary:
        jobs.append((primary.loc[primary.is_case_rule_of_two == 1], 10, "rule_of_two"))
    audits = [analyze(df, panel, label, reference, out) for df,panel,label in jobs]
    metadata = {"status": "post_hoc_exploratory", "primary_n": 7996,
                "PRS_mean_primary": mu, "PRS_SD_primary": sd, "high_PRS_cutoff_full_cohort": cutoff,
                "contrast_PRS_low": reference[3], "contrast_PRS_high": reference[4],
                "CI": "pointwise_95_percent_HC1_conditional_on_empirical_reference",
                "audits": audits}
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print("Completed exploratory aggregate age outputs; see results/common_rare_age/metadata.json.")

# ============================================================================
# 2. ehr-audit: audit_early_ehr_diagnoses
# ============================================================================
import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import norm
from statsmodels.stats.proportion import proportion_confint

EHR_AGE_AUDIT_ROOT = Path(__file__).resolve().parents[1]
EHR_AGE_AUDIT_PCS = [f"PC{ehr_age_audit_i}" for ehr_age_audit_i in range(1, 11)]
EHR_AGE_AUDIT_MIN_N = 21  # Current official policy explicitly prohibits direct counts of 1–20.
OUT = EHR_AGE_AUDIT_ROOT / "results" / "table"
EHR_AGE_AUDIT_REPORT = EHR_AGE_AUDIT_ROOT / "results" / "report"

def ehr_age_audit_read(name, **kwargs):
    return pd.read_csv(EHR_AGE_AUDIT_ROOT / "data" / name, **kwargs)

def releasable(counts):
    # Deliberately suppress zero cells as well, to keep table handling conservative.
    return bool((np.asarray(counts) >= EHR_AGE_AUDIT_MIN_N).all())

def save(frame, suffix):
    frame.to_csv(OUT / f"27_{suffix}.csv", index=False)

def table_counts(df, keys):
    return df.groupby(keys).size().reindex(
        pd.MultiIndex.from_product([[0, 1]] * len(keys)), fill_value=0)

def published_joint_margin_blocks(variable, published_path=None):
    """A released adjacent joint cell plus a new margin reveals a suppressed cell."""
    path = published_path or OUT / "18_C_mixture_category_composition.csv"
    if not path.exists():
        return False
    old = pd.read_csv(path)
    early = old.loc[(old.scheme == "cutoff_45") & (old.group == "early(<45)")]
    adjacent = "carrier only" if variable == "carrier" else "PRS-high only"
    n = pd.to_numeric(early.loc[early.mixture_category.isin(
        [adjacent, "both (PRS-high AND carrier)"]), "n_disclosed"], errors="coerce")
    return bool(n.notna().any())

def binary_table(df, variable):
    counts = table_counts(df, [variable, "early"])
    if not releasable(counts):
        return None
    early_total = int(df.early.sum())
    rows = []
    for group in [0, 1]:
        late, early = int(counts[group, 0]), int(counts[group, 1])
        n = late + early
        within_ci = proportion_confint(early, n, method="wilson")
        composition_ci = proportion_confint(early, early_total, method="wilson")
        rows.append({"variable": variable, "level": group, "n_cases": n,
                     "n_early": early, "n_later": late,
                     "early_fraction_within_group": early/n,
                     "within_ci_low": within_ci[0], "within_ci_high": within_ci[1],
                     "fraction_of_early_cases": early/early_total,
                     "composition_ci_low": composition_ci[0],
                     "composition_ci_high": composition_ci[1]})
    return pd.DataFrame(rows)

def log_odds_contrast(fit, v, name):
    v = np.asarray(v, dtype=float)
    beta = float(v @ fit.params)
    se = float(np.sqrt(max(0., v @ fit.cov_params() @ v)))
    lo, hi = beta - norm.ppf(.975)*se, beta + norm.ppf(.975)*se
    return {"contrast": name, "log_OR": beta, "se_HC1": se, "OR": np.exp(beta),
            "OR_ci_low": np.exp(lo), "OR_ci_high": np.exp(hi)}

def case_only(df, panel, label, ref):
    mu, sd, cutoff = ref
    df = df.copy()
    df["carrier"] = (df[f"burden_rare_main_ptv_{panel}gene"] > 0).astype(int)
    df["high"] = (df.PRS311 >= cutoff).astype(int)
    df["early"] = (df.age_at_diagnosis < 45).astype(int)
    four = table_counts(df, ["carrier", "high", "early"])
    audit = {"analysis": label, "four_groups_early_release": releasable(four)}
    if audit["four_groups_early_release"]:
        rows = []
        for c in [0, 1]:
            for h in [0, 1]:
                late, early = int(four[c, h, 0]), int(four[c, h, 1])
                ci = proportion_confint(early, early+late, method="wilson")
                comp_ci = proportion_confint(early, int(df.early.sum()), method="wilson")
                rows.append({"carrier": c, "high_PRS": h, "n_early": early,
                             "n_later": late, "early_fraction_within_group": early/(early+late),
                             "within_ci_low": ci[0], "within_ci_high": ci[1],
                             "fraction_of_early_cases": early/df.early.sum(),
                             "composition_ci_low": comp_ci[0], "composition_ci_high": comp_ci[1]})
        save(pd.DataFrame(rows), f"{label}_early_four_groups")
    # Different binary margins do not uniquely identify the suppressed intersection.
    # Release the primary margins only: overlapping sensitivity tables can triangulate.
    if label == "primary":
        for variable in ["carrier", "high"]:
            margins = binary_table(df, variable)
            blocked = published_joint_margin_blocks(variable)
            if blocked:
                margins = None
            audit[f"{variable}_margin_blocked_by_published_joint_cells"] = blocked
            audit[f"{variable}_margins_release"] = margins is not None
            if margins is not None:
                save(margins, f"{label}_early_by_{variable}")
            else:
                stale = OUT / f"27_{label}_early_by_{variable}.csv"
                if stale.exists():
                    stale.unlink()  # Only this script's suppressed aggregate artifact.
    if not releasable(table_counts(df, ["carrier", "early"])):
        audit["continuous_case_only_model"] = "suppressed_insufficient_carrier_outcome_cells"
        return audit
    df["PRS_z"] = (df.PRS311-mu)/sd
    df["interaction"] = df.PRS_z*df.carrier
    # PCs rescaled for numerical conditioning; genetic terms retain original units.
    scaled_pc = df[EHR_AGE_AUDIT_PCS].apply(lambda col: (col-col.mean())/col.std(ddof=1))
    assert np.isfinite(scaled_pc.to_numpy()).all()
    X = sm.add_constant(pd.concat([df[["PRS_z", "carrier", "interaction"]], scaled_pc], axis=1),
                        has_constant="add")
    fit = sm.GLM(df.early, X, family=sm.families.Binomial()).fit(cov_type="HC1")
    assert fit.converged and np.isfinite(fit.params).all()
    vectors = {}
    for i, name in [(1, "PRS_per_SD_noncarriers"), (2, "carrier_at_primary_mean_PRS"),
                    (3, "PRS_carrier_interaction")]:
        v = np.zeros(X.shape[1]); v[i] = 1; vectors[name] = v
    v = np.zeros(X.shape[1]); v[1] = v[3] = 1
    vectors["PRS_per_SD_carriers"] = v
    save(pd.DataFrame([log_odds_contrast(fit, v, k) for k, v in vectors.items()]),
         f"{label}_early_continuous_case_only")
    audit["continuous_case_only_model"] = "released_early_vs_later_case_OR_not_population_risk"
    return audit

def ehr_query(cdr):
    domain_columns = [("condition_occurrence", "condition_start_date"),
                      ("procedure_occurrence", "procedure_date"),
                      ("drug_exposure", "drug_exposure_start_date"),
                      ("measurement", "measurement_date"),
                      ("visit_occurrence", "visit_start_date")]
    union = "\nUNION ALL\n".join(
        f"SELECT t.person_id, t.{col} AS event_date FROM `{cdr}.{table}` t JOIN cohort USING(person_id)"
        for table, col in domain_columns)
    return f"""WITH cohort AS (SELECT person_id FROM UNNEST(@person_ids) person_id),
birth AS (SELECT person_id, DATE(birth_datetime) AS birth_date FROM `{cdr}.person` JOIN cohort USING(person_id)),
domain_dates AS ({union}),
bc AS (
 SELECT co.person_id, MIN(DATE(co.condition_start_datetime)) AS first_bc_date
 FROM `{cdr}.condition_occurrence` co JOIN cohort USING(person_id)
 JOIN `{cdr}.concept_ancestor` ca ON co.condition_concept_id=ca.descendant_concept_id
 WHERE ca.ancestor_concept_id=4112853 GROUP BY co.person_id
),
death AS (
 SELECT c.person_id, MIN(COALESCE(ad.death_date,d.death_date)) AS death_date
 FROM cohort c LEFT JOIN `{cdr}.death` d USING(person_id)
 LEFT JOIN (SELECT person_id, death_date FROM `{cdr}.aou_death` WHERE primary_death_record=TRUE) ad USING(person_id)
 GROUP BY c.person_id
),
dates AS (
 SELECT DISTINCT d.person_id, d.event_date, b.birth_date, bc.first_bc_date,
 DATE_DIFF(d.event_date,b.birth_date,DAY)/365.25 AS age
 FROM domain_dates d JOIN birth b USING(person_id) LEFT JOIN bc USING(person_id)
 WHERE d.event_date IS NOT NULL
),
adult_dates AS (
 SELECT *, DATE_DIFF(event_date,LAG(event_date) OVER(PARTITION BY person_id ORDER BY event_date),DAY) AS previous_gap_days
 FROM dates WHERE age>=18 AND age<45
),
ehr AS (
 SELECT person_id, MIN(event_date) AS ehr_first_date, MAX(event_date) AS ehr_last_date,
 MIN(IF(age>=18,event_date,NULL)) AS first_adult_ehr_date,
 COUNT(DISTINCT IF(age>=18 AND age<45,event_date,NULL)) AS n_adult_pre45_dates,
 COUNT(DISTINCT IF(age>=18 AND age<45,CAST(FLOOR(age) AS INT64),NULL)) AS n_adult_pre45_year_bins,
 COUNT(DISTINCT IF(first_bc_date IS NOT NULL AND event_date<first_bc_date
   AND event_date>=DATE_SUB(first_bc_date,INTERVAL 365 DAY),event_date,NULL)) AS n_dates_year_before_bc
 FROM dates GROUP BY person_id
),
gaps AS (
 SELECT person_id, MAX(previous_gap_days) AS max_adult_pre45_gap_days
 FROM adult_dates GROUP BY person_id
),
op AS (
 SELECT o.person_id, COUNT(*) AS n_observation_periods,
 MIN(observation_period_start_date) AS op_first_date,
 MAX(observation_period_end_date) AS op_last_date,
 COUNTIF(observation_period_start_date<=DATE_ADD(b.birth_date,INTERVAL 44 YEAR)
   AND observation_period_end_date>=DATE_ADD(b.birth_date,INTERVAL 45 YEAR)) AS n_periods_covering_44_to45,
 COUNTIF(observation_period_start_date>observation_period_end_date) AS n_invalid_periods
 FROM `{cdr}.observation_period` o JOIN cohort USING(person_id) JOIN birth b USING(person_id)
 GROUP BY o.person_id
)
SELECT c.person_id,b.birth_date,bc.first_bc_date,death.death_date,ehr.* EXCEPT(person_id),
 gaps.max_adult_pre45_gap_days,op.* EXCEPT(person_id)
FROM cohort c JOIN birth b USING(person_id) LEFT JOIN bc USING(person_id)
LEFT JOIN death USING(person_id) LEFT JOIN ehr USING(person_id) LEFT JOIN gaps USING(person_id)
LEFT JOIN op USING(person_id)
"""

def query_windows(ids, maximum_gib):
    from google.cloud import bigquery
    cdr = os.environ["WORKSPACE_CDR"]
    client = bigquery.Client(project=os.environ["GOOGLE_PROJECT"])
    params = [bigquery.ArrayQueryParameter("person_ids", "INT64", [int(i) for i in ids])]
    sql = ehr_query(cdr)
    dry = client.query(sql, job_config=bigquery.QueryJobConfig(
        dry_run=True, use_query_cache=False, query_parameters=params))
    gib = dry.total_bytes_processed/1024**3
    print(f"EHR query dry-run: {gib:.2f} GiB; no participant records printed.", flush=True)
    if gib > maximum_gib:
        raise RuntimeError("EHR query exceeds fixed scan cap; review query before running.")
    config = bigquery.QueryJobConfig(query_parameters=params,
                                    maximum_bytes_billed=int(maximum_gib*1024**3))
    dates = client.query(sql, job_config=config).result(timeout=600).to_dataframe()
    source_sql = f"""
SELECT o.period_type_concept_id, c.concept_name,
 COUNT(DISTINCT o.person_id) AS n_participants
FROM `{cdr}.observation_period` o
JOIN UNNEST(@person_ids) cohort_id ON o.person_id=cohort_id
LEFT JOIN `{cdr}.concept` c ON o.period_type_concept_id=c.concept_id
GROUP BY o.period_type_concept_id,c.concept_name
"""
    sources = client.query(source_sql, job_config=config).result(timeout=120).to_dataframe()
    assert len(dates) == len(ids) and not dates.person_id.duplicated().any()
    return dates, sources, gib

def numeric_summary(values):
    s = pd.Series(values).dropna()
    if len(s) < EHR_AGE_AUDIT_MIN_N:
        return None
    # Avoid exposing small missingness counts through an available n and known total.
    return {"median": s.median(), "q25": s.quantile(.25),
            "q75": s.quantile(.75)}

def span_proxy_masks(df, lookback):
    """Feasibility only: entry after washout; censor at clinical end/death/45."""
    birthday45 = df.birth_date.map(lambda d: d + pd.DateOffset(years=45))
    end_date = pd.concat([df.ehr_last_date, df.death_date, birthday45], axis=1).min(axis=1)
    entry_date = df.first_adult_ehr_date + pd.to_timedelta(lookback, unit="D")
    recorded_age = (df.first_bc_date-df.birth_date).dt.days/365.25
    eligible = (entry_date < end_date) & (entry_date < birthday45) & (
        df.first_bc_date.isna() | (df.first_bc_date > entry_date))
    event = eligible & (df.first_bc_date <= end_date) & (recorded_age < 45)
    return eligible, event

def window_audit(full, dates, sources, cutoff):
    df = full.merge(dates, validate="one_to_one")
    datecols = ["birth_date", "first_bc_date", "death_date", "ehr_first_date", "ehr_last_date",
                "first_adult_ehr_date", "op_first_date", "op_last_date"]
    for col in datecols:
        df[col] = pd.to_datetime(df[col])
    audit = {"n_source_selected": len(df), "other_cancer_excluded_from_source": True,
             "population_risk_ready": False,
             "risk_gate_reason": "outcome_selected_source_and_unvalidated_EHR_capture",
             "observation_period_presence": bool(df.n_observation_periods.notna().any())}
    # Check event identity without publishing a count of discordant records.
    recorded_age = (df.first_bc_date-df.birth_date).dt.days/365.25
    case = df.group == "case"
    delta = (recorded_age[case]-df.loc[case, "age_at_diagnosis"]).abs()
    audit["first_BC_dates_match_saved_case_ages"] = bool(delta.notna().all() and (delta < 1/365.25).all())
    assert audit["first_BC_dates_match_saved_case_ages"], "BC date/age mismatch: review phenotype definition."
    if len(sources) and releasable(sources.n_participants):
        save(sources, "observation_period_provenance")
        audit["observation_period_provenance_released"] = True
    else:
        audit["observation_period_provenance_released"] = False
    df["carrier"] = (df.burden_rare_main_ptv_10gene > 0).astype(int)
    df["high"] = (df.PRS311 >= cutoff).astype(int)
    df["age_first_ehr"] = (df.ehr_first_date-df.birth_date).dt.days/365.25
    df["age_last_ehr"] = (df.ehr_last_date-df.birth_date).dt.days/365.25
    df["pre_BC_record_span_years"] = (df.first_bc_date-df.ehr_first_date).dt.days/365.25
    # Medians of existing dates never imply clinical coverage between them.
    rows = []
    for group, sub in df.groupby("group"):
        for col in ["age_first_ehr", "age_last_ehr", "n_adult_pre45_dates", "n_adult_pre45_year_bins",
                    "max_adult_pre45_gap_days", "pre_BC_record_span_years", "n_dates_year_before_bc"]:
            summary = numeric_summary(sub[col])
            if summary:
                rows.append({"source_group": group, "measure": col, **summary})
    save(pd.DataFrame(rows), "ehr_window_distributions")
    # Mutually exclusive proxy categories, each with its complement above threshold.
    statuses = np.select([df.age_first_ehr.isna(), df.age_first_ehr >= 45,
                         (df.age_first_ehr < 45) & (df.age_last_ehr < 45)],
                        ["no_clinical_EHR", "first_EHR_at_or_after45", "EHR_ends_before45"],
                        default="records_span45")
    counts = pd.Series(statuses).value_counts().reindex(
        ["no_clinical_EHR", "first_EHR_at_or_after45", "EHR_ends_before45", "records_span45"], fill_value=0)
    audit["window_partition_released"] = bool(((counts == 0) | (counts >= EHR_AGE_AUDIT_MIN_N)).all())
    if audit["window_partition_released"]:
        save(counts.rename_axis("window_status").reset_index(name="n"), "ehr_window_partition")
    # Each binary audit reports both cells; no count is recoverable by subtraction.
    checks = {
        "OP_end_after_last_clinical_record": df.op_last_date > df.ehr_last_date,
        "OP_start_before_first_clinical_record": df.op_first_date < df.ehr_first_date,
        "single_inferred_OP_covers_age44_to45": df.n_periods_covering_44_to45.fillna(0) > 0,
    }
    coverage_rows = []
    for name, indicator in checks.items():
        yes, no = int(indicator.sum()), int((~indicator).sum())
        if (yes == 0 or yes >= EHR_AGE_AUDIT_MIN_N) and (no == 0 or no >= EHR_AGE_AUDIT_MIN_N):
            coverage_rows.append({"check": name, "n_yes": yes, "n_no": no,
                                  "fraction_yes": yes/len(df)})
    save(pd.DataFrame(coverage_rows), "inferred_observation_period_checks")
    # Censor at min(last clinical event, death, 45th birthday), never max(EHR,death).
    proxy_rows = []
    event_arrays = []
    for lookback in [365, 730]:
        eligible, event = span_proxy_masks(df, lookback)
        proxy_rows.append({"lookback_days": lookback, "n_span_proxy_eligible": int(eligible.sum()),
                           "n_recorded_early_events": int(event.sum()),
                           "n_proxy_nonevents": int((eligible & ~event).sum())})
        event_arrays.append(event)
        four = df.loc[eligible].assign(early=event.loc[eligible].astype(int))
        audit[f"lookback_{lookback}_four_group_event_gate"] = releasable(table_counts(four, ["carrier", "high", "early"]))
        if lookback == 365:
            all_group_counts = df.groupby(["carrier", "high"]).size().reindex(
                pd.MultiIndex.from_product([[0, 1], [0, 1]]), fill_value=0)
            eligible_counts = four.groupby(["carrier", "high"]).size().reindex(all_group_counts.index, fill_value=0)
            excluded_counts = all_group_counts-eligible_counts
            safe = releasable(eligible_counts) and bool(((excluded_counts == 0) | (excluded_counts >= EHR_AGE_AUDIT_MIN_N)).all())
            audit["one_year_proxy_four_group_denominators_release"] = safe
            if safe:
                save(eligible_counts.rename_axis(["carrier", "high_PRS"]).reset_index(name="n_span_proxy_eligible"),
                     "one_year_span_proxy_four_group_denominators")
        audit[f"lookback_{lookback}_proxy_eligible_release"] = (
            int(eligible.sum()) >= EHR_AGE_AUDIT_MIN_N and int((~eligible).sum()) >= EHR_AGE_AUDIT_MIN_N)
    # Avoid revealing a small nested difference in early events across washouts.
    event_diff = int((event_arrays[0] & ~event_arrays[1]).sum())
    event_release = all(r["n_recorded_early_events"] >= EHR_AGE_AUDIT_MIN_N and r["n_proxy_nonevents"] >= EHR_AGE_AUDIT_MIN_N for r in proxy_rows)
    event_release &= (event_diff == 0 or event_diff >= EHR_AGE_AUDIT_MIN_N)
    eligible_diff = proxy_rows[0]["n_span_proxy_eligible"]-proxy_rows[1]["n_span_proxy_eligible"]
    eligible_release = all(audit[f"lookback_{x}_proxy_eligible_release"] for x in [365,730])
    eligible_release &= (eligible_diff == 0 or eligible_diff >= EHR_AGE_AUDIT_MIN_N)
    if eligible_release:
        released = pd.DataFrame(proxy_rows)
        if not event_release:
            released = released[["lookback_days", "n_span_proxy_eligible"]]
        save(released, "followup_span_proxy_feasibility")
    audit["nested_proxy_event_totals_released"] = bool(event_release and eligible_release)
    audit["followup_span_proxy_totals_released"] = bool(eligible_release)
    audit["requested_additive_risk_contrast"] = "R11-R10-R01+R00"
    audit["requested_four_group_risks_and_additive_CI"] = "not_estimated_failed_observation_source_and_sparse_event_gates"
    return audit

def audit_early_ehr_diagnoses():
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-bigquery", action="store_true", help="case-only output only; EHR audit remains unfinished")
    parser.add_argument("--maximum-scan-gib", type=float, default=200)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True); EHR_AGE_AUDIT_REPORT.mkdir(parents=True, exist_ok=True)
    pheno = ehr_age_audit_read("04_phenotype_diagnostics_person_level.csv.gz")
    membership = ehr_age_audit_read("16_joint_model_cv_oof_person_level.csv.gz", usecols=["person_id"]).drop_duplicates()
    full = membership.merge(pheno, validate="one_to_one")
    for source in [ehr_age_audit_read("prs_scores.csv", usecols=["person_id", "PRS311"]),
                   ehr_age_audit_read("14_burden_rare_scores_person_level.csv.gz", usecols=["person_id", "burden_rare_main_ptv_10gene", "burden_rare_main_ptv_9gene"]),
                   ehr_age_audit_read("09_ancestry_pc_person_level.csv.gz", usecols=["person_id"]+EHR_AGE_AUDIT_PCS)]:
        full = full.merge(source, how="left", validate="one_to_one")
    assert len(full) == len(membership)
    assert full[["PRS311", "burden_rare_main_ptv_10gene", "burden_rare_main_ptv_9gene"]+EHR_AGE_AUDIT_PCS].notna().all().all()
    cases = full.loc[(full.group == "case") & (full.bc_is_first_record != True)].copy()
    assert len(cases) == 7996 and cases.age_at_diagnosis.notna().all()
    ref = (cases.PRS311.mean(), cases.PRS311.std(ddof=1), full.PRS311.quantile(.9))
    age_metadata = json.loads((EHR_AGE_AUDIT_ROOT / "results/common_rare_age/metadata.json").read_text())
    assert np.allclose(ref, [age_metadata["PRS_mean_primary"], age_metadata["PRS_SD_primary"],
                             age_metadata["high_PRS_cutoff_full_cohort"]])
    metadata = {"analysis_status": "post_hoc_exploratory", "minimum_release_n": EHR_AGE_AUDIT_MIN_N,
                "population_risk_estimated": False, "genetic_reference": {
                    "PRS_mean_primary_cases": ref[0], "PRS_SD_primary_cases": ref[1],
                    "high_cutoff_full_Aim1_cohort": ref[2]}}
    metadata["case_only"] = [case_only(cases, 10, "primary", ref),
                             case_only(cases, 9, "nine_gene", ref),
                             case_only(cases.loc[cases.is_case_rule_of_two == 1], 10, "rule_of_two", ref)]
    print("Case-only models and disclosure-checked margins completed; no small-cell counts printed.", flush=True)
    if not args.skip_bigquery:
        dates, sources, gib = query_windows(full.person_id, args.maximum_scan_gib)
        metadata["ehr_audit"] = window_audit(full, dates, sources, ref[2])
        metadata["ehr_query_dry_run_GiB"] = gib
    else:
        metadata["ehr_audit"] = {"status": "not_run"}
    (EHR_AGE_AUDIT_REPORT / "27_early_ehr_audit_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print("Completed. Aggregate outputs: results/table/27_*; audit: results/report/27_early_ehr_audit_metadata.json.", flush=True)

# ============================================================================
# 3. nine-gene-missense: analyze_nine_gene_missense_age
# ============================================================================
import argparse
import gzip
import json
import os
from pathlib import Path
import subprocess
import tempfile
from urllib.parse import unquote

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import norm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NINE_GENE_MISSENSE_AGE_ROOT = Path(__file__).resolve().parents[1]
V2 = NINE_GENE_MISSENSE_AGE_ROOT / "analysis_v2"
ANNOT = V2 / "data/variant_annotation/full"
DATA = NINE_GENE_MISSENSE_AGE_ROOT / "data/28_expanded_missense"
TABLE = NINE_GENE_MISSENSE_AGE_ROOT / "results/table"
FIGURE = NINE_GENE_MISSENSE_AGE_ROOT / "results/figures"
NINE_GENE_MISSENSE_AGE_REPORT = NINE_GENE_MISSENSE_AGE_ROOT / "results/report"
NINE_GENE_MISSENSE_AGE_PCS = [f"PC{nine_gene_missense_age_i}" for nine_gene_missense_age_i in range(1, 11)]
PANEL = ["ATM", "BARD1", "BRCA1", "BRCA2", "CHEK2", "PALB2", "RAD51C", "RAD51D", "TP53"]
PLP = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}
STARS = {"criteria_provided,_single_submitter": 1,
         "criteria_provided,_multiple_submitters,_no_conflicts": 2,
         "reviewed_by_expert_panel": 3, "practice_guideline": 4}
NINE_GENE_MISSENSE_AGE_MIN_N = 21
I157T = "chr22:28725099:A:G"

def nine_gene_missense_age_read(name, **kwargs):
    return pd.read_csv(NINE_GENE_MISSENSE_AGE_ROOT / "data" / name, **kwargs)

def safe_partition(counts):
    return bool(((np.asarray(counts) == 0) | (np.asarray(counts) >= NINE_GENE_MISSENSE_AGE_MIN_N)).all())

def related_condition(name, gene):
    text = unquote(name).replace("_", " ").lower().strip()
    if "breast" in text or "mammary" in text:
        return True
    if text == "hereditary cancer-predisposing syndrome":
        return True
    if gene in {"TP53", "CHEK2"} and ("li-fraumeni" in text or "li fraumeni" in text):
        return True
    return text == f"{gene.lower()}-related cancer predisposition" and gene != "MAP3K1"

def disease_allowed(disease_text, gene):
    names = [x for x in disease_text.split("|") if x not in {".", "", "not_provided", "not_specified"}]
    # The cache is variant-level; exclude ambiguous unrelated co-listed conditions.
    return bool(names) and all(related_condition(x, gene) for x in names)

def germline_allowed(origin):
    try:
        bits = 0
        for x in origin.split(","):
            bits |= int(x)
    except ValueError:
        return True  # CLNSIG is explicitly the germline classification field.
    germline_bits = 1 | 4 | 8 | 16 | 32 | 64 | 128
    return not (bits & 2 and not bits & germline_bits)

def prepare_candidates():
    vat = nine_gene_missense_age_read("11_panel_gene_vat_ptv_missense_candidates.csv.gz", dtype=str, keep_default_na=False)
    vat["_rank"] = np.where(vat.is_canonical_transcript == "true", 0, 1)
    vat = vat.sort_values(["vid", "gene_symbol", "_rank"]).drop_duplicates(["vid", "gene_symbol"])
    vat["eur_af"] = pd.to_numeric(vat.gvs_eur_af, errors="coerce").fillna(0.)
    vat = vat.loc[vat.gene_symbol.isin(PANEL) & vat.consequence.str.contains("missense_variant") &
                  (vat.eur_af < .001)].copy()
    vat["orig_key"] = vat.contig+":"+vat.position+":"+vat.ref_allele+":"+vat.alt_allele
    vat = vat.loc[vat.orig_key != I157T]
    wanted = set(vat.orig_key)
    mapped = []
    columns = ["orig_key", "key", "gene", "exonic_func", "aachange"]
    for chunk in pd.read_csv(ANNOT / "v2_05_annotated_sites.tsv.gz", sep="\t", dtype=str,
                             usecols=columns, chunksize=100000, keep_default_na=False):
        chunk = chunk.loc[chunk.orig_key.isin(wanted)]
        if len(chunk):
            mapped.append(chunk)
    mapping = pd.concat(mapped, ignore_index=True) if mapped else pd.DataFrame(columns=columns)
    mapping = mapping.rename(columns={"gene": "mane_gene", "exonic_func": "mane_exonic_func",
                                      "aachange": "mane_aachange"})
    assert not mapping.orig_key.duplicated().any()
    vat = vat.merge(mapping, on="orig_key", how="left", validate="many_to_one")
    # A single-base substitution needs no indel realignment; cached mapping is used otherwise.
    snv = (vat.ref_allele.str.len() == 1) & (vat.alt_allele.str.len() == 1)
    vat.loc[vat.key.isna() & snv, "key"] = vat.loc[vat.key.isna() & snv, "orig_key"]
    vat = vat.loc[vat.key.notna()].copy()
    reference = []
    shard_map = nine_gene_missense_age_read("11_panel_gene_exome_vcf_shard_map.csv")
    fmt = "%CHROM\t%POS\t%REF\t%ALT\t%ID\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/CLNSIGCONF\t%INFO/CLNDN\t%INFO/CLNDISDB\t%INFO/ORIGIN\t%INFO/MC\n"
    for row in shard_map.itertuples():
        if row.gene_symbol not in PANEL:
            continue
        region = f"{row.chrom}:{row.region_start}-{row.region_end}"
        raw = subprocess.run(["bcftools", "query", "-r", region, "-f", fmt,
                              str(ANNOT / "clinvar.norm.vcf.gz")], check=True, text=True,
                             capture_output=True).stdout
        for line in raw.splitlines():
            chrom, pos, ref, alt, vid, sig, rev, conflict, disease, disease_ids, origin, mc = line.split("\t")
            key = f"{chrom}:{pos}:{ref}:{alt}"
            reference.append({"key": key, "variation_id": vid, "clnsig": sig,
                              "clnrevstat": rev, "clnsigconf": conflict, "disease": disease,
                              "disease_ids": disease_ids, "origin": origin, "clinvar_mc": mc})
    cv = pd.DataFrame(reference).drop_duplicates()
    blocked = set(cv.loc[cv.clnsig.str.contains("Benign|benign|Conflicting|conflicting") |
                         cv.clnrevstat.str.contains("conflicting") |
                         ~cv.clnsigconf.isin(["", "."]), "key"])
    cv["stars"] = cv.clnrevstat.map(STARS).fillna(0).astype(int)
    cv = cv.loc[cv.clnsig.isin(PLP) & (cv.stars >= 1) & ~cv.key.isin(blocked) &
                cv.origin.map(germline_allowed) & cv.clinvar_mc.str.contains("missense_variant")]
    candidates = vat.merge(cv, on="key", how="inner")
    candidates = candidates.loc[[disease_allowed(d, g) for d,g in zip(candidates.disease, candidates.gene_symbol)]]
    # Keep the highest-reviewed germline assertion; all records already passed disease/conflict gates.
    candidates = candidates.sort_values("stars", ascending=False).drop_duplicates(["orig_key", "gene_symbol"])
    candidates["strict"] = candidates.stars >= 2
    mane_available = candidates.mane_exonic_func.notna()
    # Candidates contain cohort AF and must stay in ignored data/, never public tables.
    candidates.to_csv(DATA / "candidates.tsv.gz", sep="\t", index=False)
    header_date = None
    with gzip.open(ANNOT / "clinvar.norm.vcf.gz", "rt") as handle:
        for line in handle:
            if line.startswith("##fileDate="):
                header_date = line.strip().split("=", 1)[1]
            if not line.startswith("#"):
                break
    assert header_date == "2026-09-13"
    provenance = {"status": "post_hoc_exploratory", "ClinVar_snapshot": header_date,
                  "genome": "GRCh38", "annotation_reused_not_rerun": True,
                  "panel": PANEL, "primary_PTV_definition": "original_nine_gene_qualifying_PTV",
                  "reference_annotation": 'normalized ClinVar VCF and existing canonical VAT; MANE .5 audit',
                  "disease_rule": "at_least_one_relevant_named_condition_and_no_ambiguous_unrelated_named_condition",
                  "disease_specific_RCV_cache_available": False,
                  "main_minimum_stars": 1, "strict_minimum_stars": 2,
                  "eur_AF_threshold": .001, "I157T_excluded": True,
                  "GQ_min": 20, "het_VAF": [.30,.80], "hom_alt_VAF_min": .8,
                  "old_v2_genotype_cache_covers_frozen_primary": False,
                  "annotation_audit": {
                      "cached_MANE_available_for_every_candidate": bool(mane_available.all()),
                      "cached_MANE_consequences_match_where_available": bool(candidates.loc[mane_available, "mane_exonic_func"].eq("nonsynonymous SNV").all()),
                      "cached_MANE_genes_match_where_available": bool((candidates.loc[mane_available, "mane_gene"] == candidates.loc[mane_available, "gene_symbol"]).all()),
                      "independent_ClinGen_classification_field_available": False},
                  "minimum_release_n": NINE_GENE_MISSENSE_AGE_MIN_N}
    (DATA / "pre_model_rules.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    print("Prepared fixed clinical missense mask from existing reference files; no cohort allele counts printed.", flush=True)
    return candidates, provenance

def ensure_ro_mount(mount):
    mount.mkdir(parents=True, exist_ok=True)
    existing = [line.split() for line in Path("/proc/mounts").read_text().splitlines()]
    entry = next((x for x in existing if x[1] == str(mount)), None)
    if entry:
        assert "ro" in entry[3].split(","), "Mount is not read-only."
        return
    subprocess.run(["fusermount", "-u", str(mount)], capture_output=True)
    bucket = os.environ["CDR_STORAGE_PATH"].replace("gs://", "").split("/")[0]
    result = subprocess.run(["gcsfuse", "--billing-project", os.environ["GOOGLE_PROJECT"],
                             "-o", "ro", "--implicit-dirs", bucket, str(mount)],
                            capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError("Read-only gcsfuse mount failed: " + result.stderr[-1500:])
    entries = [line.split() for line in Path("/proc/mounts").read_text().splitlines()]
    entry = next(x for x in entries if x[1] == str(mount))
    assert "ro" in entry[3].split(",")

def qualifying_call(gt, ad, gq, alt_field, target):
    try:
        alleles = [int(x) for x in gt.replace("|", "/").split("/")]
        depths = [int(x) for x in ad.split(",")]
        index = alt_field.split(",").index(target)+1
        dosage = alleles.count(index)
        if dosage == 0 or len(alleles) != 2 or sum(depths) <= 0:
            return False
        vaf = depths[index]/sum(depths)
        return float(gq) >= 20 and (vaf >= .8 if dosage == 2 else .3 <= vaf <= .8)
    except (ValueError, IndexError):
        return False

def minimal_allele_key(chrom, pos, ref, alt):
    """Minimal allele key."""
    pos = int(pos)
    if alt == "*" or alt.startswith("<"):
        return chrom, pos, ref, alt
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos+1
    return chrom, pos, ref, alt

def candidate_shards(boundaries, chrom, positions):
    ordinal = int(chrom.removeprefix("chr"))
    bounds = boundaries.loc[boundaries.chrom_ord == ordinal].sort_values("start_pos")
    starts = bounds.start_pos.to_numpy()
    shards = bounds.shard_idx.to_numpy()
    chosen = set()
    for pos in positions:
        left = int(np.searchsorted(starts, pos, side="left"))
        right = int(np.searchsorted(starts, pos, side="right"))
        # Include the prior shard when an exact boundary can split alleles at a locus.
        for index in range(max(0, left-1), max(left, right)):
            chosen.add(int(shards[index]))
    return sorted(chosen)

def extract_carriers(candidates, ids, mount):
    # Sets in memory only. No genotype/person-by-variant output is written.
    main = {gene: set() for gene in PANEL}
    strict = {gene: set() for gene in PANEL}
    if candidates.empty:
        return main, strict, {"no_qualifying_missense_variants": True}
    ensure_ro_mount(mount)
    prefix = "/".join(os.environ["CDR_STORAGE_PATH"].replace("gs://", "").split("/")[1:])
    vcfdir = mount / prefix / "wgs/short_read/snpindel/exome/vcf"
    maps = nine_gene_missense_age_read("11_panel_gene_exome_vcf_shard_map.csv")
    boundaries = nine_gene_missense_age_read("11_exome_vcf_shard_start_boundaries.csv")
    seen_sites = set()
    fmt = "[%CHROM\t%POS\t%REF\t%ALT\t%SAMPLE\t%GT\t%AD\t%GQ\n]"
    with tempfile.TemporaryDirectory(prefix="expanded_missense_", dir=DATA) as temp:
        samples = Path(temp) / "samples.txt"
        samples.write_text("\n".join(map(str, sorted(ids)))+"\n", encoding="utf-8")
        for gene, group in candidates.groupby("gene_symbol"):
            row = maps.loc[maps.gene_symbol == gene].iloc[0]
            regions = Path(temp) / "regions.bed"
            group[["contig", "position"]].drop_duplicates().assign(
                start=lambda d: d.position.astype(int)-1, end=lambda d: d.position.astype(int))[
                    ["contig", "start", "end"]].to_csv(regions, sep="\t", index=False, header=False)
            wanted = {(r.contig, int(r.position), r.ref_allele, r.alt_allele): bool(r.strict)
                      for r in group.itertuples()}
            shards = candidate_shards(boundaries, row.chrom, group.position.astype(int))
            # Native padded representations can start just outside the minimal-key shard.
            shards = sorted({s+d for s in shards for d in [-1,0,1]
                             if int(row.shard_first) <= s+d <= int(row.shard_last)})
            assert shards and min(shards) >= int(row.shard_first) and max(shards) <= int(row.shard_last)
            for shard in shards:
                path = vcfdir / f"{shard:010d}.vcf.bgz"
                assert path.exists(), "Expected exome VCF shard is missing."
                sites = subprocess.run(["bcftools", "query", "-R", str(regions), "--regions-overlap", "1", "-f",
                                        "%CHROM\t%POS\t%REF\t%ALT\n", str(path)],
                                       capture_output=True, text=True, check=True).stdout
                matching = set()
                for site in sites.splitlines():
                    chrom, pos, ref, alts = site.split("\t")
                    for alt in alts.split(","):
                        key = minimal_allele_key(chrom, pos, ref, alt)
                        if key in wanted:
                            seen_sites.add(key); matching.add(key)
                if not matching:
                    continue
                view = subprocess.Popen(["bcftools", "view", "-R", str(regions), "--regions-overlap", "1", "-S", str(samples),
                                         "-Ou", str(path)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                query = subprocess.run(["bcftools", "query", "-i", 'GT="alt"', "-f", fmt],
                                       stdin=view.stdout, capture_output=True, text=True)
                view.stdout.close()
                error = view.stderr.read().decode(); ret = view.wait()
                if ret or query.returncode:
                    raise RuntimeError("Targeted genotype query failed; raw stderr withheld to protect sample IDs.")
                for line in query.stdout.splitlines():
                    chrom, pos, ref, alt_field, person, gt, ad, gq = line.split("\t")
                    person = int(person)
                    assert person in ids
                    for alt in alt_field.split(","):
                        key = minimal_allele_key(chrom, pos, ref, alt)
                        if key in wanted and qualifying_call(gt, ad, gq, alt_field, alt):
                            main[gene].add(person)
                            if wanted[key]:
                                strict[gene].add(person)
            print(f"Completed targeted exome calls for {gene}; no participant or small-cell counts printed.", flush=True)
    expected = {(r.contig, int(r.position), r.ref_allele, r.alt_allele) for r in candidates.itertuples()}
    assert seen_sites == expected, "Some qualifying missense sites were not found in the original callset."
    return main, strict, {"targeted_sites_all_found": True, "source": "original AoU exome VCF",
                          "sample_set_matches_frozen_7996_cases": True, "no_genotypes_written": True,
                          "native_VCF_padding_removed_for_allele_match": True,
                          "AD_uses_native_ALT_index": True}

def carrier_partitions(full, gene_main, gene_strict):
    main_set = set().union(*gene_main.values())
    strict_set = set().union(*gene_strict.values())
    full = full.copy()
    full["ptv"] = (full.burden_rare_main_ptv_9gene > 0).astype(int)
    full["missense"] = full.person_id.isin(main_set).astype(int)
    full["strict_missense"] = full.person_id.isin(strict_set).astype(int)
    full["expanded"] = ((full.ptv == 1) | (full.missense == 1)).astype(int)
    full["expanded_strict"] = ((full.ptv == 1) | (full.strict_missense == 1)).astype(int)
    full["category"] = np.select([(full.ptv == 0) & (full.missense == 0),
                                  (full.ptv == 1) & (full.missense == 0),
                                  (full.ptv == 0) & (full.missense == 1)],
                                 ["neither", "PTV_only", "missense_only"], default="both")
    primary = (full.group == "case") & (full.bc_is_first_record != True)
    audits = {}
    for label, sub in [("primary_cases", full.loc[primary]),
                       ("primary_age_groups", full.loc[primary].assign(age_group=lambda d: np.where(d.age_at_diagnosis < 45, "early", "later")))]:
        keys = (["age_group", "category"] if label == "primary_age_groups" else ["category"])
        levels = ["early", "later"]
        categories = ["neither", "PTV_only", "missense_only", "both"]
        if len(keys) == 1:
            counts = sub.groupby(keys[0]).size().reindex(categories, fill_value=0)
        else:
            counts = sub.groupby(keys).size().reindex(pd.MultiIndex.from_product([levels, categories], names=keys), fill_value=0)
        safe = safe_partition(counts)
        audits[label] = {"partition_release": safe}
        if safe:
            counts.rename_axis(keys).reset_index(name="n").to_csv(TABLE / f"28_{label}_carrier_categories.csv", index=False)
    # Gene summaries omit all numbers when any count or known complement is small.
    gene_rows = []
    for gene in PANEL:
        for label, sub in [("primary_case", full.loc[primary])]:
            n = int(sub.person_id.isin(gene_main[gene]).sum())
            total = len(sub)
            if (n == 0 or n >= NINE_GENE_MISSENSE_AGE_MIN_N) and (total-n == 0 or total-n >= NINE_GENE_MISSENSE_AGE_MIN_N):
                gene_rows.append({"gene": gene, "stratum": label, "n_missense_carriers": n})
            else:
                gene_rows.append({"gene": gene, "stratum": label, "n_missense_carriers": "suppressed"})
    if any(r["n_missense_carriers"] == "suppressed" for r in gene_rows):
        # Do not leave all-but-one gene counts available beside an overall missense total.
        for row in gene_rows:
            row["n_missense_carriers"] = "suppressed_joint_gene_table"
    pd.DataFrame(gene_rows).to_csv(TABLE / "28_missense_carriers_by_gene.csv", index=False)
    return full, audits

def linear_contrast(fit, v):
    v = np.asarray(v, dtype=float)
    beta = float(v @ fit.params); se = float(np.sqrt(max(0, v @ fit.cov_params() @ v)))
    return {"age_difference_years": beta, "se_HC1": se,
            "ci_low": beta-norm.ppf(.975)*se, "ci_high": beta+norm.ppf(.975)*se}

def fit_age(cases, carrier, label, ref):
    mu, sd, pc_ref = ref
    df = cases.copy()
    df["PRS_z"] = (df.PRS311-mu)/sd
    df["carrier"] = df[carrier]
    pc_mu, pc_sd = df[NINE_GENE_MISSENSE_AGE_PCS].mean(), df[NINE_GENE_MISSENSE_AGE_PCS].std(ddof=1)
    df[NINE_GENE_MISSENSE_AGE_PCS] = (df[NINE_GENE_MISSENSE_AGE_PCS]-pc_mu)/pc_sd
    pc_fixed = (pc_ref-pc_mu)/pc_sd
    rows = []
    for interaction in [False, True]:
        df["interaction"] = df.PRS_z*df.carrier
        terms = ["PRS_z", "carrier"] + (["interaction"] if interaction else []) + NINE_GENE_MISSENSE_AGE_PCS
        X = sm.add_constant(df[terms], has_constant="add")
        fit = sm.OLS(df.age_at_diagnosis.astype(float), X).fit(cov_type="HC1")
        for term in ["PRS_z", "carrier"] + (["interaction"] if interaction else []):
            v = np.zeros(X.shape[1]); v[X.columns.get_loc(term)] = 1
            rows.append({"definition": label, "model": "with_interaction" if interaction else "additive",
                         "contrast": term, **linear_contrast(fit, v)})
        if interaction:
            v = np.zeros(X.shape[1]); v[X.columns.get_loc("PRS_z")] = v[X.columns.get_loc("interaction")] = 1
            rows.append({"definition": label, "model": "with_interaction",
                         "contrast": "PRS_per_SD_carriers", **linear_contrast(fit, v)})
            lo = max(df.loc[df.carrier == c, "PRS_z"].quantile(.01) for c in [0,1])
            hi = min(df.loc[df.carrier == c, "PRS_z"].quantile(.99) for c in [0,1])
            assert lo < hi
            fig, ax = plt.subplots()
            for c in [0,1]:
                xs = np.linspace(lo, hi, 100)
                predictions = pd.DataFrame([linear_contrast(fit, [1,p,c,p*c,*pc_fixed]) for p in xs])
                line, = ax.plot(xs, predictions.age_difference_years, label=f"Carrier={c}")
                ax.fill_between(xs, predictions.ci_low, predictions.ci_high, alpha=.2, color=line.get_color())
            ax.set(xlabel="PRS311 (SD in frozen primary cases)", ylabel="Adjusted recorded diagnosis age (years)")
            ax.set_title({"PTV_only": "Nine-gene qualifying PTV",
                          "expanded_star1": "Nine-gene PTV + P/LP missense (review >=1 star)",
                          "expanded_star2": "Nine-gene PTV + P/LP missense (review >=2 stars)"}[label])
            ax.grid(False)
            ax.legend(); fig.tight_layout()
            fig.savefig(FIGURE / f"28_{label}_predicted_age.svg"); plt.close(fig)
    return rows

def finalize_disclosure_and_report():
    """Check released aggregates against earlier ten-gene totals before reporting."""
    meta_path = NINE_GENE_MISSENSE_AGE_REPORT / "28_expanded_missense_metadata.json"
    metadata = json.loads(meta_path.read_text())
    counts_path = TABLE / "28_primary_carrier_increment.csv"
    counts = pd.read_csv(counts_path)
    old_path = NINE_GENE_MISSENSE_AGE_ROOT / "results/common_rare_age/primary_four_groups.csv"
    baseline_safe = metadata.get("nine_gene_baseline_total_release", True)
    if old_path.exists():
        old = pd.read_csv(old_path)
        n10 = int(old.loc[old.carrier == 1, "n"].sum())
        n9 = pd.to_numeric(counts.loc[counts.definition == "PTV_only", "n_carriers"], errors="coerce").iloc[0]
        if pd.notna(n9):
            difference = n10-int(n9)
            baseline_safe = difference == 0 or difference >= NINE_GENE_MISSENSE_AGE_MIN_N
    if not baseline_safe:
        counts["n_carriers"] = "suppressed_cross_panel_difference"
        counts.to_csv(counts_path, index=False)
        for label in ["primary_cases", "primary_age_groups"]:
            path = TABLE / f"28_{label}_carrier_categories.csv"
            if path.exists():
                path.unlink()  # Only new tables whose sums would reveal the excluded small group.
            metadata["counts_audit"][label]["partition_release"] = False
            metadata["counts_audit"][label]["cross_panel_complement_suppressed"] = True
    metadata["nine_gene_baseline_total_release"] = baseline_safe
    genes_path = TABLE / "28_missense_carriers_by_gene.csv"
    genes = pd.read_csv(genes_path)
    if genes.n_missense_carriers.astype(str).str.contains("suppressed").any():
        genes["n_missense_carriers"] = "suppressed_joint_gene_table"
        genes.to_csv(genes_path, index=False)
    meta_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    models = pd.read_csv(TABLE / "28_age_model_comparison.csv")
    names = {"PTV_only": 'Nine genes PTV-only', "expanded_star1": 'Nine genes PTV+P/LP missense (≥1 star)',
             "expanded_star2": 'Strictly Expanded (≥2 stars)'}
    lines = ['# Ninegenes expanded qualifying-variant carrier: exploratory age at diagnosis analysis', "",
             'Using the same batch of 7,996 eligible cases, the same case PRS mean/sample SD and PC1–10, the original PTV-only analysis is still the main analysis.',
             'Reuse ClinVar 2026-09-13/GRCh38 and existing VAT and MANE .5 annotations without re-running ClinVar.', "",
             '## Add carrier', "",
             '| Definition | Total number of carriers | New number of people compared to nine-gene PTV-only |',
             "|---|---:|---:|"]
    for row in counts.itertuples():
        def display(v):
            return 'inhibition' if "suppressed" in str(v) else str(v)
        lines.append(f"| {names[row.definition]} | {display(row.n_carriers)} | {display(row.new_carriers_vs_PTV)} |")
    lines.extend(["", 'New carriers refer to those who are PTV negative for the original nine genes and carry at least one qualified P/LP missense; carriers of both categories are not counted as new additions.',
                  'New increments of less than or equal to 20 people are not directly released; the total number of ten-gene PTVs has been released, so if the total number of nine genes is reversed, a small-scale MAP3K1 difference set will be derived, and the related total number and the grouping table that can be added will also be suppressed. When the difference between the strict and main extended numbers reaches the threshold, the strict number is suppressed.', "",
                  '## Recorded diagnosis age difference corresponding to each SD of PRS (years, 95% CI)', "",
                  '| Carrier definition | Non-carrier slope | Carrier slope | Difference/interaction of two sets of slopes |',
                  "|---|---:|---:|---:|"])
    for label, name in names.items():
        subset = models.loc[(models.definition == label) & (models.model == "with_interaction")].set_index("contrast")
        formatted = []
        for term in ["PRS_z", "PRS_per_SD_carriers", "interaction"]:
            r = subset.loc[term]
            formatted.append(f"{r.age_difference_years:+.2f}（{r.ci_low:+.2f}to{r.ci_high:+.2f}）")
        lines.append(f"| {name} | " + " | ".join(formatted) + " |")
    main = models.loc[(models.definition == "expanded_star1") & (models.model == "with_interaction") &
                      (models.contrast == "interaction")].iloc[0]
    conclusion = ('The 95% CI of the interaction crosses zero, and no clear slope difference is seen; this does not mean that the slopes of the two groups have been proven to be the same.' if main.ci_low <= 0 <= main.ci_high else
                  'The 95% CI of the interaction does not cross zero, which is an exploratory interaction signal that requires independent verification.')
    lines.extend(["", conclusion, "", '## Prediction graph', "",
                  '[Nine-gene expanded carrier predicted age chart](../figures/28_expanded_star1_predicted_age.svg)', "",
                  'The horizontal axis PRS311 is the SD of the fixed case cohort, and the vertical axis is the PC-adjusted average recorded diagnosis age. The two lines represent the extended definition of carrier=0/1 respectively; the shadow is the HC1 point state 95%CI. Only the intersection of the 1–99th percentile range of each group is drawn, using the same PC mean as a reference; the shading is not the individual prediction interval.', "",
                  '## Comments and interpretation restrictions', "",
                  'The main extension requires at least 1 star, conflict-free germline P/LP, clear related diseases, AoU EUR AF<0.1%; the strict extension requires ≥2 stars. The original GQ≥20, het VAF[0.30,0.80], hom-alt VAF≥0.8 and CHEK2 I157T exclusion rules are retained. The new missense only reads fixed points from the original exome callset and folds it in memory; it does not save the new participant genotype or carrier table.', "",
                  'No disease-level RCV classification is cached. Variants that are also listed with other clear diseases and cannot confirm the classification correspondence are conservatively excluded; only not_provided/not_specified will not make them eligible. Therefore, this is a clearly identifiable subset of P/LP missenses and cannot be claimed to exhaust all relevant pathogenic missenses. The original PTV was qualifying PTV, and clinical P/LP was not required one by one, so it was collectively called expanded qualifying-variant carrier.', "",
                  'The MANE site cache does not cover all candidates; the results and genes recorded by MANE all match the original VAT, and the rest are based on the original canonical VAT missense results and the ClinVar missense annotations of the same version. There is no independent ClinGen classification field, the expert review rating comes from ClinVar CLNREVSTAT.', "",
                  'Recording age associations within cases cannot prove early-onset risk synergy or biological interactions in the population; limitations such as retaining EHR dates later than actual onset, index-event bias, and low VAF/clonal hematopoiesis are not completely eliminated. Annotations or PRS thresholds were not adjusted based on age model results.', "",
                  'Verification: The original case n=7,996 and the PRS reference assertion passed, and the newly added candidate site passed the original callset verification; the synthetic test covers disease exclusion, somatic-only, target ALT dose/VAF and the same position across shard boundaries. The main and extended models adopt the same OLS HC1 covariance rule.', "",
                  'Basis: [ClinVar Review Level](https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/), [AoU Disclosure Policy](https://support.researchallofus.org/hc/en-us/articles/22346276580372-Data-and-Statistics-Dissemination-Policy).'])
    (NINE_GENE_MISSENSE_AGE_REPORT / "28_expanded_missense_age_summary.md").write_text("\n".join(lines)+"\n", encoding="utf-8")

def analyze_nine_gene_missense_age():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["prepare", "run"], default="run")
    ap.add_argument("--mount", type=Path, default=Path("/tmp/codex_expanded_missense_aou_ro"))
    args = ap.parse_args()
    for path in [DATA, TABLE, FIGURE, NINE_GENE_MISSENSE_AGE_REPORT]:
        path.mkdir(parents=True, exist_ok=True)
    ignored = subprocess.run(["git", "check-ignore", "-q", str(DATA)], cwd=NINE_GENE_MISSENSE_AGE_ROOT)
    assert ignored.returncode == 0, "Internal candidate folder must be ignored."
    candidates, provenance = prepare_candidates()
    if args.stage == "prepare":
        return
    membership = nine_gene_missense_age_read("16_joint_model_cv_oof_person_level.csv.gz", usecols=["person_id"]).drop_duplicates()
    full = membership.merge(nine_gene_missense_age_read("04_phenotype_diagnostics_person_level.csv.gz"), validate="one_to_one")
    for source in [nine_gene_missense_age_read("prs_scores.csv", usecols=["person_id", "PRS311"]),
                   nine_gene_missense_age_read("14_burden_rare_scores_person_level.csv.gz", usecols=["person_id", "burden_rare_main_ptv_9gene"]),
                   nine_gene_missense_age_read("09_ancestry_pc_person_level.csv.gz", usecols=["person_id"]+NINE_GENE_MISSENSE_AGE_PCS)]:
        full = full.merge(source, how="left", validate="one_to_one")
    assert len(full) == len(membership)
    assert full[["PRS311", "burden_rare_main_ptv_9gene"]+NINE_GENE_MISSENSE_AGE_PCS].notna().all().all()
    frozen_cases = full.loc[(full.group == "case") & (full.bc_is_first_record != True)].copy()
    assert len(frozen_cases) == 7996
    main_carriers, strict_carriers, extraction = extract_carriers(candidates, set(frozen_cases.person_id), args.mount)
    cases, counts_audit = carrier_partitions(frozen_cases, main_carriers, strict_carriers)
    assert len(cases) == 7996 and cases.age_at_diagnosis.notna().all()
    ref = (cases.PRS311.mean(), cases.PRS311.std(ddof=1), cases[NINE_GENE_MISSENSE_AGE_PCS].mean())
    original = json.loads((NINE_GENE_MISSENSE_AGE_ROOT / "results/common_rare_age/metadata.json").read_text())
    assert np.allclose(ref[:2], [original["PRS_mean_primary"], original["PRS_SD_primary"]])
    results = []
    for variable, label in [("ptv", "PTV_only"), ("expanded", "expanded_star1"), ("expanded_strict", "expanded_star2")]:
        assert (cases[variable].value_counts().reindex([0,1], fill_value=0) >= NINE_GENE_MISSENSE_AGE_MIN_N).all()
        results.extend(fit_age(cases, variable, label, ref))
    pd.DataFrame(results).to_csv(TABLE / "28_age_model_comparison.csv", index=False)
    # Total expanded n would reveal a suppressed missense-only increment when PTV n is known.
    count_rows = []
    added_main = int(((cases.expanded == 1) & (cases.ptv == 0)).sum())
    added_strict = int(((cases.expanded_strict == 1) & (cases.ptv == 0)).sum())
    nested_increment_safe = added_main-added_strict == 0 or added_main-added_strict >= NINE_GENE_MISSENSE_AGE_MIN_N
    for variable, label in [("ptv", "PTV_only"), ("expanded", "expanded_star1"), ("expanded_strict", "expanded_star2")]:
        added = int(((cases[variable] == 1) & (cases.ptv == 0)).sum())
        total = int(cases[variable].sum())
        released = (added == 0 or added >= NINE_GENE_MISSENSE_AGE_MIN_N) and total >= NINE_GENE_MISSENSE_AGE_MIN_N and len(cases)-total >= NINE_GENE_MISSENSE_AGE_MIN_N
        if variable == "expanded_strict":
            released &= nested_increment_safe
        count_rows.append({"definition": label, "primary_n": len(cases),
                           "n_carriers": total if released else "suppressed",
                           "new_carriers_vs_PTV": added if released else "suppressed"})
    pd.DataFrame(count_rows).to_csv(TABLE / "28_primary_carrier_increment.csv", index=False)
    metadata = {**provenance, "primary_case_n": 7996, "extraction": extraction,
                "counts_audit": counts_audit, "case_PRS_mean": ref[0], "case_PRS_SD": ref[1],
                "risk_interaction_estimated": False, "records_not_actual_onset": True}
    (NINE_GENE_MISSENSE_AGE_REPORT / "28_expanded_missense_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    finalize_disclosure_and_report()
    print("Completed expanded-carrier age sensitivity; public outputs contain aggregates only.", flush=True)

# ============================================================================
# 4. legacy-prs-age: describe_prs_and_age
# ============================================================================

import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

def prs_age_descriptives_log(msg):
    print(msg, flush=True)

def describe_prs_and_age():
    """Initialize runtime state and execute this component explicitly."""
    global AGE_BINS, AGE_GROUP_COLOR, AGE_LABELS, PRS_AGE_DESCRIPTIVES_BASE, PRS_AGE_DESCRIPTIVES_DATA_DIR, DECADE_BINS, DECADE_LABELS, PRS_AGE_DESCRIPTIVES_DPI, PRS_AGE_DESCRIPTIVES_FIG_DIR, PRS_AGE_DESCRIPTIVES_TABLE_DIR, prs_age_descriptives__, prs_age_descriptives_ax, prs_age_descriptives_axes, prs_age_descriptives_case_df, col, decade_summary_df, decade_summary_path, decade_summary_rows, decade_test_df, decade_test_path, decade_test_rows, decade_violin_path, prs_age_descriptives_df, early, prs_age_descriptives_fig, groups, prs_age_descriptives_grp, kw_p, kw_stat, late, means, mw_p, mw_stat, n_age_missing, n_missing, prs_age_descriptives_pheno, prs_age_descriptives_prs, prs_col, prs_age_descriptives_qc, prs_age_descriptives_rho, prs_age_descriptives_rho_p, s, scatter_path, prs_age_descriptives_sub, summary_df, summary_path, summary_rows, test_df, test_path, test_rows, violin_path

    matplotlib.use("Agg")
    plt.style.use("default")

    PRS_AGE_DESCRIPTIVES_BASE = "/home/jupyter/workspace/analysis"
    PRS_AGE_DESCRIPTIVES_DATA_DIR = os.path.join(PRS_AGE_DESCRIPTIVES_BASE, "data")
    
    PRS_AGE_DESCRIPTIVES_FIG_DIR = PRS_AGE_DESCRIPTIVES_TABLE_DIR = os.path.join(PRS_AGE_DESCRIPTIVES_DATA_DIR, "10_internal")

    os.makedirs(PRS_AGE_DESCRIPTIVES_FIG_DIR, exist_ok=True)
    os.makedirs(PRS_AGE_DESCRIPTIVES_TABLE_DIR, exist_ok=True)

    PRS_AGE_DESCRIPTIVES_DPI = 150
    AGE_GROUP_COLOR = "tab:blue"

    # ---------------------------------------------------------------------------
    
    # ---------------------------------------------------------------------------
    prs_age_descriptives_log("=" * 80)
    prs_age_descriptives_log('Step A: Read three output tables and join them, filter case + pass_qc==1')
    prs_age_descriptives_log("=" * 80)

    prs_age_descriptives_pheno = pd.read_csv(os.path.join(PRS_AGE_DESCRIPTIVES_DATA_DIR, "04_phenotype_diagnostics_person_level.csv.gz"))
    prs_age_descriptives_qc = pd.read_csv(os.path.join(PRS_AGE_DESCRIPTIVES_DATA_DIR, "05_cohort_qc_flags.csv"))
    prs_age_descriptives_prs = pd.read_csv(os.path.join(PRS_AGE_DESCRIPTIVES_DATA_DIR, "prs_scores.csv"))

    prs_age_descriptives_df = prs_age_descriptives_pheno.merge(prs_age_descriptives_qc[["person_id", "pass_qc"]], on="person_id", how="inner")
    prs_age_descriptives_df = prs_age_descriptives_df.merge(prs_age_descriptives_prs, on="person_id", how="inner")
    prs_age_descriptives_log(f"After joining the three tables:{prs_age_descriptives_df.shape}")

    prs_age_descriptives_case_df = prs_age_descriptives_df[(prs_age_descriptives_df["group"] == "case") & (prs_age_descriptives_df["pass_qc"] == 1)].copy()
    prs_age_descriptives_log(f"case and pass_qc==1:{len(prs_age_descriptives_case_df)}people")

    n_age_missing = prs_age_descriptives_case_df["age_at_diagnosis"].isna().sum()
    prs_age_descriptives_log(f"age_at_diagnosis is missing:{n_age_missing}people({n_age_missing/len(prs_age_descriptives_case_df)*100:.2f}%）")
    prs_age_descriptives_case_df = prs_age_descriptives_case_df.dropna(subset=["age_at_diagnosis"])
    prs_age_descriptives_log(f"After removing the missing age_at_diagnosis:{len(prs_age_descriptives_case_df)}people")

    for prs_col in ["PRS313", "PRS311"]:
        n_missing = prs_age_descriptives_case_df[prs_col].isna().sum()
        prs_age_descriptives_log(f"{prs_col}Missing:{n_missing}people")
    prs_age_descriptives_case_df = prs_age_descriptives_case_df.dropna(subset=["PRS313", "PRS311"])
    prs_age_descriptives_log(f"Final analysis sample with none of the three variables missing:{len(prs_age_descriptives_case_df)}people")

    # ---------------------------------------------------------------------------
    
    # ---------------------------------------------------------------------------
    prs_age_descriptives_log("\n" + "=" * 80)
    prs_age_descriptives_log('Step B: age_at_diagnosis clinical cutoff grouping (early onset <45 / moderate 45-64 / late onset >=65)')
    prs_age_descriptives_log("=" * 80)

    AGE_BINS = [0, 45, 65, 150]
    AGE_LABELS = ["early-onset (<45)", "mid-onset (45-64)", "late-onset (>=65)"]
    prs_age_descriptives_case_df["age_group"] = pd.cut(prs_age_descriptives_case_df["age_at_diagnosis"], bins=AGE_BINS,
                                   labels=AGE_LABELS, right=False)
    prs_age_descriptives_log('Age_at_diagnosis range and number of people in each group:')
    prs_age_descriptives_log(str(prs_age_descriptives_case_df.groupby("age_group", observed=True)["age_at_diagnosis"]
            .agg(["count", "min", "median", "max"])))

    # ---------------------------------------------------------------------------
    
    # ---------------------------------------------------------------------------
    prs_age_descriptives_log("\n" + "=" * 80)
    prs_age_descriptives_log('Step C: PRS313/PRS311 summary statistics for each age group + Kruskal-Wallis + early onset vs late onset post-hoc test')
    prs_age_descriptives_log("=" * 80)

    summary_rows = []
    for prs_col in ["PRS313", "PRS311"]:
        for prs_age_descriptives_grp, prs_age_descriptives_sub in prs_age_descriptives_case_df.groupby("age_group", observed=True):
            s = prs_age_descriptives_sub[prs_col]
            summary_rows.append({
                "prs": prs_col, "age_group": prs_age_descriptives_grp, "n": len(s),
                "age_at_diagnosis_median": prs_age_descriptives_sub["age_at_diagnosis"].median(),
                "mean": s.mean(), "median": s.median(), "sd": s.std(),
                "iqr_low": s.quantile(0.25), "iqr_high": s.quantile(0.75),
            })

    summary_df = pd.DataFrame(summary_rows)
    for col in ["age_at_diagnosis_median", "mean", "median", "sd", "iqr_low", "iqr_high"]:
        summary_df[col] = summary_df[col].round(4)

    test_rows = []
    for prs_col in ["PRS313", "PRS311"]:
        groups = [prs_age_descriptives_sub[prs_col].values for prs_age_descriptives__, prs_age_descriptives_sub in prs_age_descriptives_case_df.groupby("age_group", observed=True)]
        kw_stat, kw_p = stats.kruskal(*groups)
        test_rows.append({"prs": prs_col, "test": "Kruskal-Wallis (3 age groups: early/mid/late)",
                           "statistic": round(kw_stat, 4), "p_value": kw_p, "n": len(prs_age_descriptives_case_df)})
        prs_age_descriptives_log(f"\n{prs_col} ~ age_group  Kruskal-Wallis H={kw_stat:.4f}, p={kw_p:.4g}")

        early = prs_age_descriptives_case_df.loc[prs_age_descriptives_case_df["age_group"] == AGE_LABELS[0], prs_col]
        late = prs_age_descriptives_case_df.loc[prs_age_descriptives_case_df["age_group"] == AGE_LABELS[2], prs_col]
        mw_stat, mw_p = stats.mannwhitneyu(early, late, alternative="two-sided")
        test_rows.append({"prs": prs_col, "test": "Mann-Whitney U (early-onset vs late-onset)",
                           "statistic": round(mw_stat, 4), "p_value": mw_p,
                           "n": f"{len(early)} vs {len(late)}"})
        prs_age_descriptives_log(f"{prs_col}  early-onset(n={len(early)}) vs late-onset(n={len(late)})  "
            f"Mann-Whitney U={mw_stat:.1f}, p={mw_p:.4g}")

        prs_age_descriptives_rho, prs_age_descriptives_rho_p = stats.spearmanr(prs_age_descriptives_case_df["age_at_diagnosis"], prs_age_descriptives_case_df[prs_col])
        test_rows.append({"prs": prs_col, "test": "Spearman (age_at_diagnosis, continuous)",
                           "statistic": round(prs_age_descriptives_rho, 4), "p_value": prs_age_descriptives_rho_p, "n": len(prs_age_descriptives_case_df)})
        prs_age_descriptives_log(f"{prs_col}vs age_at_diagnosis (continuous) Spearman rho={prs_age_descriptives_rho:.4f}, p={prs_age_descriptives_rho_p:.4g}")

    test_df = pd.DataFrame(test_rows)

    summary_path = os.path.join(PRS_AGE_DESCRIPTIVES_TABLE_DIR, "10_prs_by_age_at_diagnosis_summary.csv")
    summary_df.to_csv(summary_path, index=False)
    prs_age_descriptives_log(f"\n has written:{summary_path}")
    prs_age_descriptives_log(summary_df.to_string(index=False))

    test_path = os.path.join(PRS_AGE_DESCRIPTIVES_TABLE_DIR, "10_prs_by_age_at_diagnosis_tests.csv")
    test_df.to_csv(test_path, index=False)
    prs_age_descriptives_log(f"\n has written:{test_path}")
    prs_age_descriptives_log(test_df.to_string(index=False))

    # ---------------------------------------------------------------------------
    
    # ---------------------------------------------------------------------------
    prs_age_descriptives_log("\n" + "=" * 80)
    prs_age_descriptives_log('Step D: Publish the picture')
    prs_age_descriptives_log("=" * 80)

    prs_age_descriptives_fig, prs_age_descriptives_axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for prs_age_descriptives_ax, prs_col in zip(prs_age_descriptives_axes, ["PRS313", "PRS311"]):
        sns.violinplot(x="age_group", y=prs_col, data=prs_age_descriptives_case_df, order=AGE_LABELS,
                        color=AGE_GROUP_COLOR, inner="quartile", cut=0, ax=prs_age_descriptives_ax)
        prs_age_descriptives_ax.set_title(f"{prs_col} by age_at_diagnosis group (case only)")
        prs_age_descriptives_ax.set_xlabel("age_at_diagnosis group")
        prs_age_descriptives_ax.set_ylabel(prs_col)
        prs_age_descriptives_ax.grid(False)
    prs_age_descriptives_fig.tight_layout()
    violin_path = os.path.join(PRS_AGE_DESCRIPTIVES_FIG_DIR, "10_prs_by_age_at_diagnosis_group_violin.png")
    prs_age_descriptives_fig.savefig(violin_path, dpi=PRS_AGE_DESCRIPTIVES_DPI)
    plt.close(prs_age_descriptives_fig)
    prs_age_descriptives_log(f"Has written:{violin_path}")

    prs_age_descriptives_fig, prs_age_descriptives_axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for prs_age_descriptives_ax, prs_col in zip(prs_age_descriptives_axes, ["PRS313", "PRS311"]):
        sns.regplot(x="age_at_diagnosis", y=prs_col, data=prs_age_descriptives_case_df, ax=prs_age_descriptives_ax,
                    scatter_kws={"alpha": 0.15, "s": 8, "color": AGE_GROUP_COLOR},
                    line_kws={"color": "tab:red"})
        prs_age_descriptives_ax.set_title(f"{prs_col} vs age_at_diagnosis (case only)")
        prs_age_descriptives_ax.set_xlabel("age_at_diagnosis")
        prs_age_descriptives_ax.set_ylabel(prs_col)
        prs_age_descriptives_ax.grid(False)
    prs_age_descriptives_fig.tight_layout()
    scatter_path = os.path.join(PRS_AGE_DESCRIPTIVES_FIG_DIR, "10_prs_vs_age_at_diagnosis_scatter.jpg")
    prs_age_descriptives_fig.savefig(scatter_path, dpi=PRS_AGE_DESCRIPTIVES_DPI, pil_kwargs={"quality": 90})
    plt.close(prs_age_descriptives_fig)
    prs_age_descriptives_log(f"Has written:{scatter_path}")

    # ---------------------------------------------------------------------------

    # ---------------------------------------------------------------------------
    prs_age_descriptives_log("\n" + "=" * 80)
    prs_age_descriptives_log('Step E: age_at_diagnosis decade segmentation (<40/40-49/50-59/60-69/70-79/80+)')
    prs_age_descriptives_log("=" * 80)

    DECADE_BINS = [0, 40, 50, 60, 70, 80, 150]
    DECADE_LABELS = ["<40", "40-49", "50-59", "60-69", "70-79", "80+"]
    prs_age_descriptives_case_df["age_decade"] = pd.cut(prs_age_descriptives_case_df["age_at_diagnosis"], bins=DECADE_BINS,
                                    labels=DECADE_LABELS, right=False)
    prs_age_descriptives_log('Age_at_diagnosis range and number of people in each decade:')
    prs_age_descriptives_log(str(prs_age_descriptives_case_df.groupby("age_decade", observed=True)["age_at_diagnosis"]
            .agg(["count", "min", "median", "max"])))

    decade_summary_rows = []
    for prs_col in ["PRS313", "PRS311"]:
        for prs_age_descriptives_grp, prs_age_descriptives_sub in prs_age_descriptives_case_df.groupby("age_decade", observed=True):
            s = prs_age_descriptives_sub[prs_col]
            decade_summary_rows.append({
                "prs": prs_col, "age_decade": prs_age_descriptives_grp, "n": len(s),
                "age_at_diagnosis_median": prs_age_descriptives_sub["age_at_diagnosis"].median(),
                "mean": s.mean(), "median": s.median(), "sd": s.std(),
                "iqr_low": s.quantile(0.25), "iqr_high": s.quantile(0.75),
            })
    decade_summary_df = pd.DataFrame(decade_summary_rows)
    for col in ["age_at_diagnosis_median", "mean", "median", "sd", "iqr_low", "iqr_high"]:
        decade_summary_df[col] = decade_summary_df[col].round(4)

    decade_test_rows = []
    for prs_col in ["PRS313", "PRS311"]:
        groups = [prs_age_descriptives_sub[prs_col].values for prs_age_descriptives__, prs_age_descriptives_sub in prs_age_descriptives_case_df.groupby("age_decade", observed=True)]
        kw_stat, kw_p = stats.kruskal(*groups)
        decade_test_rows.append({"prs": prs_col, "test": "Kruskal-Wallis (6 decade groups)",
                                  "statistic": round(kw_stat, 4), "p_value": kw_p, "n": len(prs_age_descriptives_case_df)})
        prs_age_descriptives_log(f"\n{prs_col} ~ age_decade  Kruskal-Wallis H={kw_stat:.4f}, p={kw_p:.4g}")
    decade_test_df = pd.DataFrame(decade_test_rows)

    decade_summary_path = os.path.join(PRS_AGE_DESCRIPTIVES_TABLE_DIR, "10_prs_by_age_decade_summary.csv")
    decade_summary_df.to_csv(decade_summary_path, index=False)
    prs_age_descriptives_log(f"\n has written:{decade_summary_path}")
    prs_age_descriptives_log(decade_summary_df.to_string(index=False))

    decade_test_path = os.path.join(PRS_AGE_DESCRIPTIVES_TABLE_DIR, "10_prs_by_age_decade_tests.csv")
    decade_test_df.to_csv(decade_test_path, index=False)
    prs_age_descriptives_log(f"\n has written:{decade_test_path}")
    prs_age_descriptives_log(decade_test_df.to_string(index=False))

    prs_age_descriptives_fig, prs_age_descriptives_axes = plt.subplots(1, 2, figsize=(13, 4.5))
    for prs_age_descriptives_ax, prs_col in zip(prs_age_descriptives_axes, ["PRS313", "PRS311"]):
        sns.violinplot(x="age_decade", y=prs_col, data=prs_age_descriptives_case_df, order=DECADE_LABELS,
                        color=AGE_GROUP_COLOR, inner="quartile", cut=0, ax=prs_age_descriptives_ax)
        means = prs_age_descriptives_case_df.groupby("age_decade", observed=True)[prs_col].mean().reindex(DECADE_LABELS)
        prs_age_descriptives_ax.plot(range(len(DECADE_LABELS)), means.values, marker="o", color="tab:red",
                linewidth=1.5, label="group mean")
        prs_age_descriptives_ax.set_title(f"{prs_col} by age_at_diagnosis decade (case only)")
        prs_age_descriptives_ax.set_xlabel("age_at_diagnosis decade")
        prs_age_descriptives_ax.set_ylabel(prs_col)
        prs_age_descriptives_ax.legend()
        prs_age_descriptives_ax.grid(False)
    prs_age_descriptives_fig.tight_layout()
    decade_violin_path = os.path.join(PRS_AGE_DESCRIPTIVES_FIG_DIR, "10_prs_by_age_decade_violin.png")
    prs_age_descriptives_fig.savefig(decade_violin_path, dpi=PRS_AGE_DESCRIPTIVES_DPI)
    plt.close(prs_age_descriptives_fig)
    prs_age_descriptives_log(f"Has written:{decade_violin_path}")

    prs_age_descriptives_log('Finish.')

# ============================================================================
# 5. legacy-burden-distribution: describe_legacy_burden
# ============================================================================

import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

def legacy_burden_descriptives_log(msg):
    print(msg, flush=True)

def legacy_burden_descriptives_section(title):
    legacy_burden_descriptives_log("\n" + "=" * 80)
    legacy_burden_descriptives_log(title)
    legacy_burden_descriptives_log("=" * 80)

def legacy_burden_descriptives_suppress(n, min_cell=20):
    """Legacy burden descriptives suppress."""
    if n is None or (isinstance(n, float) and np.isnan(n)):
        return n
    return n if n >= min_cell else f"<{min_cell}"

def legacy_burden_descriptives_wilson_ci(k, n, z=1.959963984540054):
    """Wilson score 95% CI for a binomial proportion."""
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    denom = 1 + z ** 2 / n
    center = (p + z ** 2 / (2 * n)) / denom
    half = (z * np.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2))) / denom
    return (max(0.0, center - half), min(1.0, center + half))

def gene_breakdown(sub_df, group_label):
    """Gene breakdown."""
    carriers = sub_df[sub_df["carrier"] == 1]
    n_total_carriers = len(carriers)
    n_multi = int((carriers["n_genes_carried"] >= 2).sum())
    n_single = int((carriers["n_genes_carried"] == 1).sum())
    
    legacy_burden_descriptives_log(f"\n[{group_label}]Total number of carriers={legacy_burden_descriptives_suppress(n_total_carriers)}")

    rows = []
    singles = carriers[carriers["n_genes_carried"] == 1]
    for g in gene_order:
        n_g = int(singles[g].sum())
        rows.append({"group": group_label, "category": g, "n": n_g, "kind": "single_gene"})
    rows.append({"group": group_label, "category": "compound (2+ genes)", "n": n_multi,
                 "kind": "compound"})

    breakdown_df = pd.DataFrame(rows)
    # sanity check: sum of singles + multi == n_total_carriers
    check_sum = breakdown_df["n"].sum()
    assert check_sum == n_total_carriers, (
        f"{group_label}: sum of singleton+compound({check_sum}) != total number of carriers({n_total_carriers})"
    )
    return breakdown_df, n_total_carriers

def pool_set(breakdown_df, min_cell=20):
    """Pool set."""
    df = breakdown_df.sort_values("n").reset_index(drop=True)
    pool = list(df.loc[df["n"] < min_cell, "category"])
    if not pool:
        return set()
    rest = [c for c in df["category"] if c not in pool]
    while (df.loc[df["category"].isin(pool), "n"].sum() < min_cell or len(pool) < 2) and rest:
        pool.append(rest.pop(0))
    return set(pool)

def merge_categories(breakdown_df, pool):
    """Merge categories."""
    df = breakdown_df.sort_values("n").reset_index(drop=True)
    if not pool:
        out = breakdown_df.copy()
        out["n_display"] = out["n"]
        return out
    in_pool = df["category"].isin(pool)
    pooled_n = int(df.loc[in_pool, "n"].sum())
    assert pooled_n >= LEGACY_BURDEN_DESCRIPTIVES_DISCLOSURE_MIN_CELL and in_pool.sum() >= 2, 'Still <20 or only contains one category after merging'
    pooled = pd.DataFrame([{"group": df["group"].iloc[0],
                            "category": "pooled: " + " + ".join(df.loc[in_pool, "category"]),
                            "n": pooled_n, "kind": "pooled"}])
    out = pd.concat([df.loc[~in_pool], pooled], ignore_index=True)
    out["n_display"] = out["n"]
    return out

def fmt_p(p):
    return f"{p:.2e}" if p < 0.0001 else f"{p:.4g}"

def describe_legacy_burden():
    """Initialize runtime state and execute this component explicitly."""
    global LEGACY_BURDEN_DESCRIPTIVES_BASE, CASE_COLOR, CONTROL_COLOR, LEGACY_BURDEN_DESCRIPTIVES_DATA_DIR, LEGACY_BURDEN_DESCRIPTIVES_DISCLOSURE_MIN_CELL, LEGACY_BURDEN_DESCRIPTIVES_DPI, LEGACY_BURDEN_DESCRIPTIVES_FIG_DIR, GROUP_COLORS, INTERNAL_DIR, PANEL_GENES_10, LEGACY_BURDEN_DESCRIPTIVES_REPORT_DIR, SCORE_COL, LEGACY_BURDEN_DESCRIPTIVES_TABLE_DIR, legacy_burden_descriptives__, _pool, legacy_burden_descriptives_a, age_by_carrier_path, age_by_carrier_summary, age_carrier, age_carrier_only_df, age_noncarrier, age_plot_df, age_tests, age_tests_path, legacy_burden_descriptives_ax, legacy_burden_descriptives_axes, legacy_burden_descriptives_b, bar_colors, legacy_burden_descriptives_burden, legacy_burden_descriptives_c, carrier, carrier_rate_summary, carrier_rate_summary_path, carrier_rate_test, carrier_rate_test_path, case_age_df, case_breakdown, case_breakdown_display, case_breakdown_path, legacy_burden_descriptives_case_df, case_full, case_nz, case_plot_df, case_single_lines, ci_case, ci_control, cohort, control_breakdown, control_breakdown_display, control_breakdown_path, control_full, control_nz, control_plot_df, control_single_lines, legacy_burden_descriptives_d, legacy_burden_descriptives_df, direction, err_high, err_low, fh, legacy_burden_descriptives_fig, fig1_path, fig2_path, fig3_path, fisher_p, full_dist_test, full_dist_test_path, legacy_burden_descriptives_g, gene_or_main, gene_order, groups_order, legacy_burden_descriptives_i, log_or, mw_age_p, mw_age_stat, mw_full_p, mw_full_stat, mw_nz_p, mw_nz_stat, n_carr, n_case, n_case_age_missing, n_case_carrier, n_case_carriers_chk, n_control, n_control_carrier, n_control_carriers_chk, n_g, nonzero_summary, nonzero_summary_path, nz_plot_df, odds_ratio, or_ci_high, or_ci_low, order_case, order_control, legacy_burden_descriptives_p_case, p_control, p_pool, pcts, legacy_burden_descriptives_pheno, legacy_burden_descriptives_qc, legacy_burden_descriptives_report, legacy_burden_descriptives_rho, legacy_burden_descriptives_rho_p, legacy_burden_descriptives_row, se_log_or, se_pool, table_2x2, tick, z_p, z_stat

    matplotlib.use("Agg")
    plt.style.use("default")

    LEGACY_BURDEN_DESCRIPTIVES_BASE = "/home/jupyter/workspace/analysis"
    LEGACY_BURDEN_DESCRIPTIVES_DATA_DIR = os.path.join(LEGACY_BURDEN_DESCRIPTIVES_BASE, "data")
    LEGACY_BURDEN_DESCRIPTIVES_FIG_DIR = os.path.join(LEGACY_BURDEN_DESCRIPTIVES_BASE, "results", "figures")
    LEGACY_BURDEN_DESCRIPTIVES_TABLE_DIR = os.path.join(LEGACY_BURDEN_DESCRIPTIVES_BASE, "results", "table")
    LEGACY_BURDEN_DESCRIPTIVES_REPORT_DIR = os.path.join(LEGACY_BURDEN_DESCRIPTIVES_BASE, "results", "report")

    os.makedirs(LEGACY_BURDEN_DESCRIPTIVES_FIG_DIR, exist_ok=True)
    os.makedirs(LEGACY_BURDEN_DESCRIPTIVES_TABLE_DIR, exist_ok=True)
    os.makedirs(LEGACY_BURDEN_DESCRIPTIVES_REPORT_DIR, exist_ok=True)

    LEGACY_BURDEN_DESCRIPTIVES_DPI = 150
    CASE_COLOR = "tab:blue"
    CONTROL_COLOR = "tab:orange"
    GROUP_COLORS = {"case": CASE_COLOR, "control": CONTROL_COLOR}

    SCORE_COL = "burden_rare_main_ptv_10gene"
    PANEL_GENES_10 = ["ATM", "BARD1", "BRCA1", "BRCA2", "CHEK2", "PALB2",
                       "RAD51C", "RAD51D", "TP53", "MAP3K1"]
    LEGACY_BURDEN_DESCRIPTIVES_DISCLOSURE_MIN_CELL = 20

    # ---------------------------------------------------------------------------
    
    # ---------------------------------------------------------------------------
    legacy_burden_descriptives_section('Step A: Read four output tables and join them, filter case/control + pass_qc==1')

    legacy_burden_descriptives_pheno = pd.read_csv(os.path.join(LEGACY_BURDEN_DESCRIPTIVES_DATA_DIR, "04_phenotype_diagnostics_person_level.csv.gz"))
    legacy_burden_descriptives_qc = pd.read_csv(os.path.join(LEGACY_BURDEN_DESCRIPTIVES_DATA_DIR, "05_cohort_qc_flags.csv"))
    legacy_burden_descriptives_burden = pd.read_csv(os.path.join(LEGACY_BURDEN_DESCRIPTIVES_DATA_DIR, "14_burden_rare_scores_person_level.csv.gz"))
    carrier = pd.read_csv(os.path.join(LEGACY_BURDEN_DESCRIPTIVES_DATA_DIR, "13_carrier_status_main_ptv_maf1pct_person_level.csv.gz"))

    legacy_burden_descriptives_df = (legacy_burden_descriptives_pheno.merge(legacy_burden_descriptives_qc[["person_id", "pass_qc"]], on="person_id", how="inner")
               .merge(legacy_burden_descriptives_burden[["person_id", SCORE_COL]], on="person_id", how="inner")
               .merge(carrier, on="person_id", how="inner"))
    legacy_burden_descriptives_log(f"After joining the four tables:{legacy_burden_descriptives_df.shape}")

    cohort = legacy_burden_descriptives_df[(legacy_burden_descriptives_df["pass_qc"] == 1) & (legacy_burden_descriptives_df["group"].isin(["case", "control"]))].copy()
    legacy_burden_descriptives_log(f"cohort (pass_qc==1 & group in case/control): {cohort.shape}")
    legacy_burden_descriptives_log(str(cohort["group"].value_counts()))

    cohort["carrier"] = (cohort[SCORE_COL] > 0).astype(int)
    n_case = int((cohort["group"] == "case").sum())
    n_control = int((cohort["group"] == "control").sum())
    n_case_carrier = int(cohort.loc[cohort["group"] == "case", "carrier"].sum())
    n_control_carrier = int(cohort.loc[cohort["group"] == "control", "carrier"].sum())
    legacy_burden_descriptives_log(f"\ncase: n={n_case}, nonzero(carrier)={n_case_carrier} ({n_case_carrier/n_case*100:.3f}%)")
    legacy_burden_descriptives_log(f"control: n={n_control}, nonzero(carrier)={n_control_carrier} ({n_control_carrier/n_control*100:.3f}%)")

    # ===========================================================================
    
    # ===========================================================================
    legacy_burden_descriptives_section('Question 1a: carrier rate (whether >0) case vs control — Proportional test + Fisher exact OR')

    legacy_burden_descriptives_p_case = n_case_carrier / n_case
    p_control = n_control_carrier / n_control
    ci_case = legacy_burden_descriptives_wilson_ci(n_case_carrier, n_case)
    ci_control = legacy_burden_descriptives_wilson_ci(n_control_carrier, n_control)
    legacy_burden_descriptives_log(f"case carrier rate ={legacy_burden_descriptives_p_case*100:.4f}%  (95% Wilson CI {ci_case[0]*100:.4f}%-{ci_case[1]*100:.4f}%)")
    legacy_burden_descriptives_log(f"control carrier rate ={p_control*100:.4f}%  (95% Wilson CI {ci_control[0]*100:.4f}%-{ci_control[1]*100:.4f}%)")

    # 2x2 table: rows=carrier(1)/non-carrier(0), cols=case/control
    table_2x2 = np.array([
        [n_case_carrier, n_control_carrier],
        [n_case - n_case_carrier, n_control - n_control_carrier],
    ])
    odds_ratio, fisher_p = stats.fisher_exact(table_2x2)
    legacy_burden_descriptives_log(f"\nFisher exact test: OR(case carrier vs control carrier) = {odds_ratio:.4f}, p = {fisher_p:.4g}")

    # log-OR SE (Woolf/Haldane, standard 2x2 formula) for 95% CI on OR
    legacy_burden_descriptives_a, legacy_burden_descriptives_b, legacy_burden_descriptives_c, legacy_burden_descriptives_d = table_2x2[0, 0], table_2x2[0, 1], table_2x2[1, 0], table_2x2[1, 1]
    log_or = np.log(odds_ratio)
    se_log_or = np.sqrt(1 / legacy_burden_descriptives_a + 1 / legacy_burden_descriptives_b + 1 / legacy_burden_descriptives_c + 1 / legacy_burden_descriptives_d)
    or_ci_low = np.exp(log_or - 1.959963984540054 * se_log_or)
    or_ci_high = np.exp(log_or + 1.959963984540054 * se_log_or)
    legacy_burden_descriptives_log(f"OR 95% CI (Woolf log-OR method): {or_ci_low:.4f} - {or_ci_high:.4f}")

    # Two-proportion z-test (pooled variance under H0)
    p_pool = (n_case_carrier + n_control_carrier) / (n_case + n_control)
    se_pool = np.sqrt(p_pool * (1 - p_pool) * (1 / n_case + 1 / n_control))
    z_stat = (legacy_burden_descriptives_p_case - p_control) / se_pool
    z_p = 2 * stats.norm.sf(abs(z_stat))  
    legacy_burden_descriptives_log(f"\nTwo-proportion z-test: z = {z_stat:.4f}, p = {z_p:.4g}")

    carrier_rate_summary = pd.DataFrame([
        {"group": "case", "n": n_case, "n_carrier": n_case_carrier,
         "pct_carrier": round(legacy_burden_descriptives_p_case * 100, 4),
         "wilson_ci_low_pct": round(ci_case[0] * 100, 4),
         "wilson_ci_high_pct": round(ci_case[1] * 100, 4)},
        {"group": "control", "n": n_control, "n_carrier": n_control_carrier,
         "pct_carrier": round(p_control * 100, 4),
         "wilson_ci_low_pct": round(ci_control[0] * 100, 4),
         "wilson_ci_high_pct": round(ci_control[1] * 100, 4)},
    ])
    carrier_rate_summary_path = os.path.join(LEGACY_BURDEN_DESCRIPTIVES_TABLE_DIR, "15_burden_rare_carrier_rate_by_group.csv")
    carrier_rate_summary.to_csv(carrier_rate_summary_path, index=False)
    legacy_burden_descriptives_log(f"\n has written:{carrier_rate_summary_path}")
    legacy_burden_descriptives_log(carrier_rate_summary.to_string(index=False))

    carrier_rate_test = pd.DataFrame([
        {"test": "Fisher exact (2x2: carrier x case/control)", "statistic_name": "odds_ratio",
         "statistic": round(odds_ratio, 4), "ci_low": round(or_ci_low, 4), "ci_high": round(or_ci_high, 4),
         "p_value": fisher_p},
        {"test": "Two-proportion z-test (pooled)", "statistic_name": "z",
         "statistic": round(z_stat, 4), "ci_low": np.nan, "ci_high": np.nan, "p_value": z_p},
    ])
    carrier_rate_test_path = os.path.join(LEGACY_BURDEN_DESCRIPTIVES_TABLE_DIR, "15_burden_rare_carrier_rate_test.csv")
    carrier_rate_test.to_csv(carrier_rate_test_path, index=False)
    legacy_burden_descriptives_log(f"\n has written:{carrier_rate_test_path}")
    legacy_burden_descriptives_log(carrier_rate_test.to_string(index=False))

    # ===========================================================================
    
    # ===========================================================================
    legacy_burden_descriptives_section('Question 1b: nonzero subset (score>0) Burden_rare distribution case vs control — Mann-Whitney U')

    case_nz = cohort.loc[(cohort["group"] == "case") & (cohort["carrier"] == 1), SCORE_COL]
    control_nz = cohort.loc[(cohort["group"] == "control") & (cohort["carrier"] == 1), SCORE_COL]

    nonzero_summary = pd.DataFrame([
        {"group": "case", "n": len(case_nz), "mean": round(case_nz.mean(), 4),
         "median": round(case_nz.median(), 4), "sd": round(case_nz.std(), 4),
         "q25": round(case_nz.quantile(0.25), 4), "q75": round(case_nz.quantile(0.75), 4),
         "min": round(case_nz.min(), 4), "max": round(case_nz.max(), 4)},
        {"group": "control", "n": len(control_nz), "mean": round(control_nz.mean(), 4),
         "median": round(control_nz.median(), 4), "sd": round(control_nz.std(), 4),
         "q25": round(control_nz.quantile(0.25), 4), "q75": round(control_nz.quantile(0.75), 4),
         "min": round(control_nz.min(), 4), "max": round(control_nz.max(), 4)},
    ])
    legacy_burden_descriptives_log(nonzero_summary.to_string(index=False))
    nonzero_summary_path = os.path.join(LEGACY_BURDEN_DESCRIPTIVES_TABLE_DIR, "15_burden_rare_nonzero_subset_summary.csv")
    nonzero_summary.to_csv(nonzero_summary_path, index=False)
    legacy_burden_descriptives_log(f"\n has written:{nonzero_summary_path}")

    mw_nz_stat, mw_nz_p = stats.mannwhitneyu(case_nz, control_nz, alternative="two-sided")
    legacy_burden_descriptives_log(f"\nMann-Whitney U (nonzero subset, case vs control): U = {mw_nz_stat:.1f}, p = {mw_nz_p:.4g}")

    # ===========================================================================
    
    # ===========================================================================
    legacy_burden_descriptives_section('Question 1c: Full distribution (including 0) Burden_rare case vs control — Mann-Whitney U')

    case_full = cohort.loc[cohort["group"] == "case", SCORE_COL]
    control_full = cohort.loc[cohort["group"] == "control", SCORE_COL]
    mw_full_stat, mw_full_p = stats.mannwhitneyu(case_full, control_full, alternative="two-sided")
    legacy_burden_descriptives_log(f"Mann-Whitney U (full distribution incl. zeros, case vs control): "
        f"U = {mw_full_stat:.1f}, p = {mw_full_p:.4g}")

    full_dist_test = pd.DataFrame([
        {"test": "Mann-Whitney U (nonzero subset only, case vs control)",
         "n_case": len(case_nz), "n_control": len(control_nz),
         "statistic": round(mw_nz_stat, 1), "p_value": mw_nz_p},
        {"test": "Mann-Whitney U (full distribution incl. zeros, case vs control)",
         "n_case": n_case, "n_control": n_control,
         "statistic": round(mw_full_stat, 1), "p_value": mw_full_p},
    ])
    full_dist_test_path = os.path.join(LEGACY_BURDEN_DESCRIPTIVES_TABLE_DIR, "15_burden_rare_distribution_tests.csv")
    full_dist_test.to_csv(full_dist_test_path, index=False)
    legacy_burden_descriptives_log(f"\n has written:{full_dist_test_path}")
    legacy_burden_descriptives_log(full_dist_test.to_string(index=False))

    # ===========================================================================
    
    # ===========================================================================
    legacy_burden_descriptives_section('Figure 1: results/figures/15_burden_rare_case_control_comparison.png')

    legacy_burden_descriptives_fig, legacy_burden_descriptives_axes = plt.subplots(1, 2, figsize=(11, 4.5))

    legacy_burden_descriptives_ax = legacy_burden_descriptives_axes[0]
    groups_order = ["control", "case"]
    pcts = [p_control * 100, legacy_burden_descriptives_p_case * 100]
    err_low = [p_control * 100 - ci_control[0] * 100, legacy_burden_descriptives_p_case * 100 - ci_case[0] * 100]
    err_high = [ci_control[1] * 100 - p_control * 100, ci_case[1] * 100 - legacy_burden_descriptives_p_case * 100]
    bar_colors = [GROUP_COLORS[legacy_burden_descriptives_g] for legacy_burden_descriptives_g in groups_order]
    legacy_burden_descriptives_ax.bar(groups_order, pcts, yerr=[err_low, err_high], capsize=5, color=bar_colors)
    legacy_burden_descriptives_ax.set_ylabel("% carrier (Burden_rare > 0)")
    legacy_burden_descriptives_ax.set_title("Carrier rate by case/control\n(error bars = 95% Wilson CI)")
    legacy_burden_descriptives_ax.grid(False)
    for legacy_burden_descriptives_i, legacy_burden_descriptives_g in enumerate(groups_order):
        n_g = n_control if legacy_burden_descriptives_g == "control" else n_case
        n_carr = n_control_carrier if legacy_burden_descriptives_g == "control" else n_case_carrier
        legacy_burden_descriptives_ax.text(legacy_burden_descriptives_i, pcts[legacy_burden_descriptives_i] + err_high[legacy_burden_descriptives_i] + 0.15, f"{n_carr}/{n_g}", ha="center", fontsize=8)

    legacy_burden_descriptives_ax = legacy_burden_descriptives_axes[1]
    nz_plot_df = pd.concat([
        pd.DataFrame({"group": "case", SCORE_COL: case_nz}),
        pd.DataFrame({"group": "control", SCORE_COL: control_nz}),
    ])
    sns.violinplot(x="group", y=SCORE_COL, hue="group", data=nz_plot_df, order=groups_order[::-1],
                   hue_order=groups_order[::-1], palette=GROUP_COLORS, legend=False,
                   inner="quartile", cut=0, ax=legacy_burden_descriptives_ax)
    legacy_burden_descriptives_ax.set_title(f"Burden_rare among carriers only\n(score > 0, n_case={len(case_nz)}, n_control={len(control_nz)})")
    legacy_burden_descriptives_ax.set_xlabel("group")
    legacy_burden_descriptives_ax.set_ylabel(SCORE_COL)
    legacy_burden_descriptives_ax.grid(False)

    legacy_burden_descriptives_fig.tight_layout()
    fig1_path = os.path.join(LEGACY_BURDEN_DESCRIPTIVES_FIG_DIR, "15_burden_rare_case_control_comparison.png")
    legacy_burden_descriptives_fig.savefig(fig1_path, dpi=LEGACY_BURDEN_DESCRIPTIVES_DPI)
    plt.close(legacy_burden_descriptives_fig)
    legacy_burden_descriptives_log(f"Has written:{fig1_path}")

    # ===========================================================================
    
    # ===========================================================================
    legacy_burden_descriptives_section('Question 2 & 4: Gene attribution distribution of carriers within case/control (10 panel gene 0/1 columns)')

    cohort["n_genes_carried"] = cohort[PANEL_GENES_10].sum(axis=1)

    gene_or_main = {
        "ATM": 2.1, "BARD1": 2.09, "BRCA1": 10.57, "BRCA2": 5.85, "CHEK2": 2.54,
        "PALB2": 5.02, "RAD51C": 1.93, "RAD51D": 1.8, "TP53": 3.06, "MAP3K1": 6.25,
    }
    gene_order = sorted(PANEL_GENES_10, key=lambda g: np.log(gene_or_main[g]))
    legacy_burden_descriptives_log(f"Genes are arranged in ascending order of log(OR):{gene_order}")

    case_breakdown, n_case_carriers_chk = gene_breakdown(cohort[cohort["group"] == "case"], "case")
    control_breakdown, n_control_carriers_chk = gene_breakdown(cohort[cohort["group"] == "control"], "control")

    _pool = pool_set(case_breakdown) | pool_set(control_breakdown)
    case_breakdown_display = merge_categories(case_breakdown, _pool)
    control_breakdown_display = merge_categories(control_breakdown, _pool)
    legacy_burden_descriptives_log('Gene-by-gene attribution tables/figures are only available in data/15_internal/ (gene-by-gene numbers are only disclosed in the specification table No. 19) and will not be printed here.')

    INTERNAL_DIR = os.path.join(LEGACY_BURDEN_DESCRIPTIVES_DATA_DIR, "15_internal")
    os.makedirs(INTERNAL_DIR, exist_ok=True)
    case_breakdown_path = os.path.join(INTERNAL_DIR, "15_burden_rare_case_internal_gene_breakdown.csv")
    case_breakdown_display[["group", "category", "kind", "n_display"]].rename(
        columns={"n_display": "n"}
    ).to_csv(case_breakdown_path, index=False)
    legacy_burden_descriptives_log(f"\n has written:{case_breakdown_path}")

    control_breakdown_path = os.path.join(INTERNAL_DIR, "15_burden_rare_control_internal_gene_breakdown.csv")
    control_breakdown_display[["group", "category", "kind", "n_display"]].rename(
        columns={"n_display": "n"}
    ).to_csv(control_breakdown_path, index=False)
    legacy_burden_descriptives_log(f"Has written:{control_breakdown_path}")

    case_plot_df = case_breakdown_display
    control_plot_df = control_breakdown_display

    # ===========================================================================
    
    # ===========================================================================
    legacy_burden_descriptives_section('Figure 2 (internal): data/15_internal/15_burden_rare_gene_breakdown_case_vs_control.png')

    legacy_burden_descriptives_fig, legacy_burden_descriptives_axes = plt.subplots(1, 2, figsize=(13, 5))

    legacy_burden_descriptives_ax = legacy_burden_descriptives_axes[0]
    order_case = case_plot_df.sort_values("n", ascending=False)["category"].tolist()
    sns.barplot(x="category", y="n", data=case_plot_df, order=order_case, color=CASE_COLOR, ax=legacy_burden_descriptives_ax)
    legacy_burden_descriptives_ax.set_title(f"Case-internal carrier breakdown\n(n={n_case_carriers_chk} carriers)")
    legacy_burden_descriptives_ax.set_xlabel("gene of origin (or compound)")
    legacy_burden_descriptives_ax.set_ylabel("n carriers")
    legacy_burden_descriptives_ax.tick_params(axis="x", rotation=45)
    for tick in legacy_burden_descriptives_ax.get_xticklabels():
        tick.set_ha("right")
    legacy_burden_descriptives_ax.grid(False)

    legacy_burden_descriptives_ax = legacy_burden_descriptives_axes[1]
    order_control = control_plot_df.sort_values("n", ascending=False)["category"].tolist()
    sns.barplot(x="category", y="n", data=control_plot_df, order=order_control, color=CONTROL_COLOR, ax=legacy_burden_descriptives_ax)
    legacy_burden_descriptives_ax.set_title(f"Control-internal carrier breakdown\n(n={n_control_carriers_chk} carriers)")
    legacy_burden_descriptives_ax.set_xlabel("gene of origin (or compound)")
    legacy_burden_descriptives_ax.set_ylabel("n carriers")
    legacy_burden_descriptives_ax.tick_params(axis="x", rotation=45)
    for tick in legacy_burden_descriptives_ax.get_xticklabels():
        tick.set_ha("right")
    legacy_burden_descriptives_ax.grid(False)

    legacy_burden_descriptives_fig.tight_layout()
    fig2_path = os.path.join(INTERNAL_DIR, "15_burden_rare_gene_breakdown_case_vs_control.png")
    legacy_burden_descriptives_fig.savefig(fig2_path, dpi=LEGACY_BURDEN_DESCRIPTIVES_DPI)
    plt.close(legacy_burden_descriptives_fig)
    legacy_burden_descriptives_log(f"Has written:{fig2_path}")

    # ===========================================================================
    
    # ===========================================================================
    legacy_burden_descriptives_section('Question 3: Burden_rare is related to age_at_diagnosis inside the case')

    legacy_burden_descriptives_case_df = cohort[cohort["group"] == "case"].copy()
    n_case_age_missing = legacy_burden_descriptives_case_df["age_at_diagnosis"].isna().sum()
    legacy_burden_descriptives_log(f"age_at_diagnosis is missing in case:{n_case_age_missing}people ({n_case_age_missing/len(legacy_burden_descriptives_case_df)*100:.3f}%)")
    case_age_df = legacy_burden_descriptives_case_df.dropna(subset=["age_at_diagnosis"]).copy()
    legacy_burden_descriptives_log(f"Case analysis sample after excluding age_at_diagnosis missing:{len(case_age_df)}people")

    age_carrier = case_age_df.loc[case_age_df["carrier"] == 1, "age_at_diagnosis"]
    age_noncarrier = case_age_df.loc[case_age_df["carrier"] == 0, "age_at_diagnosis"]
    legacy_burden_descriptives_log(f"\ncarrier (n={len(age_carrier)}): age_at_diagnosis median={age_carrier.median():.2f}, "
        f"mean={age_carrier.mean():.2f}")
    legacy_burden_descriptives_log(f"non-carrier (n={len(age_noncarrier)}): age_at_diagnosis median={age_noncarrier.median():.2f}, "
        f"mean={age_noncarrier.mean():.2f}")

    age_by_carrier_summary = pd.DataFrame([
        {"carrier_status": "carrier (Burden_rare>0)", "n": len(age_carrier),
         "mean": round(age_carrier.mean(), 4), "median": round(age_carrier.median(), 4),
         "sd": round(age_carrier.std(), 4),
         "q25": round(age_carrier.quantile(0.25), 4), "q75": round(age_carrier.quantile(0.75), 4)},
        {"carrier_status": "non-carrier (Burden_rare==0)", "n": len(age_noncarrier),
         "mean": round(age_noncarrier.mean(), 4), "median": round(age_noncarrier.median(), 4),
         "sd": round(age_noncarrier.std(), 4),
         "q25": round(age_noncarrier.quantile(0.25), 4), "q75": round(age_noncarrier.quantile(0.75), 4)},
    ])
    age_by_carrier_path = os.path.join(LEGACY_BURDEN_DESCRIPTIVES_TABLE_DIR, "15_burden_rare_age_at_diagnosis_by_carrier_status.csv")
    age_by_carrier_summary.to_csv(age_by_carrier_path, index=False)
    legacy_burden_descriptives_log(f"\n has written:{age_by_carrier_path}")
    legacy_burden_descriptives_log(age_by_carrier_summary.to_string(index=False))

    mw_age_stat, mw_age_p = stats.mannwhitneyu(age_carrier, age_noncarrier, alternative="two-sided")
    legacy_burden_descriptives_log(f"\nMann-Whitney U (age_at_diagnosis, carrier vs non-carrier): U = {mw_age_stat:.1f}, p = {mw_age_p:.4g}")
    direction = 'earlier' if age_carrier.median() < age_noncarrier.median() else 'Later (or similar)'
    legacy_burden_descriptives_log(f"Direction: carrier group age_at_diagnosis median{direction}In the non-carrier group"
        f"({age_carrier.median():.2f} vs {age_noncarrier.median():.2f})")

    age_carrier_only_df = case_age_df[case_age_df["carrier"] == 1]
    legacy_burden_descriptives_rho, legacy_burden_descriptives_rho_p = stats.spearmanr(age_carrier_only_df[SCORE_COL], age_carrier_only_df["age_at_diagnosis"])
    legacy_burden_descriptives_log(f"\nSpearman (inside the carrier, n={len(age_carrier_only_df)}): "
        f"Burden_rare score vs age_at_diagnosis: rho = {legacy_burden_descriptives_rho:.4f}, p = {legacy_burden_descriptives_rho_p:.4g}")

    age_tests = pd.DataFrame([
        {"test": "Mann-Whitney U (age_at_diagnosis, carrier vs non-carrier)",
         "n_carrier": len(age_carrier), "n_noncarrier": len(age_noncarrier),
         "statistic": round(mw_age_stat, 1), "p_value": mw_age_p},
        {"test": "Spearman (Burden_rare score vs age_at_diagnosis, carriers only)",
         "n_carrier": len(age_carrier_only_df), "n_noncarrier": np.nan,
         "statistic": round(legacy_burden_descriptives_rho, 4), "p_value": legacy_burden_descriptives_rho_p},
    ])
    age_tests_path = os.path.join(LEGACY_BURDEN_DESCRIPTIVES_TABLE_DIR, "15_burden_rare_age_at_diagnosis_tests.csv")
    age_tests.to_csv(age_tests_path, index=False)
    legacy_burden_descriptives_log(f"\n has written:{age_tests_path}")
    legacy_burden_descriptives_log(age_tests.to_string(index=False))

    # ===========================================================================
    
    # ===========================================================================
    legacy_burden_descriptives_section('Figure 3: results/figures/15_burden_rare_age_at_diagnosis_by_carrier_status.png')

    legacy_burden_descriptives_fig, legacy_burden_descriptives_ax = plt.subplots(figsize=(6, 5))
    age_plot_df = case_age_df.copy()
    age_plot_df["carrier_label"] = age_plot_df["carrier"].map(
        {1: "carrier\n(Burden_rare>0)", 0: "non-carrier\n(Burden_rare==0)"}
    )
    sns.violinplot(x="carrier_label", y="age_at_diagnosis", data=age_plot_df,
                   order=["non-carrier\n(Burden_rare==0)", "carrier\n(Burden_rare>0)"],
                   color=CASE_COLOR, inner="quartile", cut=0, ax=legacy_burden_descriptives_ax)
    legacy_burden_descriptives_ax.set_title(f"age_at_diagnosis by carrier status (case only, n={len(case_age_df)})")
    legacy_burden_descriptives_ax.set_xlabel("carrier status")
    legacy_burden_descriptives_ax.set_ylabel("age_at_diagnosis")
    legacy_burden_descriptives_ax.grid(False)
    legacy_burden_descriptives_fig.tight_layout()
    fig3_path = os.path.join(LEGACY_BURDEN_DESCRIPTIVES_FIG_DIR, "15_burden_rare_age_at_diagnosis_by_carrier_status.png")
    legacy_burden_descriptives_fig.savefig(fig3_path, dpi=LEGACY_BURDEN_DESCRIPTIVES_DPI)
    plt.close(legacy_burden_descriptives_fig)
    legacy_burden_descriptives_log(f"Has written:{fig3_path}")

    # ===========================================================================
    
    # ===========================================================================
    legacy_burden_descriptives_section('Writing results/report/15_burden_rare_distributional_analysis_summary.md')

    case_single_lines = []
    for legacy_burden_descriptives__, legacy_burden_descriptives_row in case_breakdown_display.iterrows():
        case_single_lines.append(f"  - {legacy_burden_descriptives_row['category']}: n={legacy_burden_descriptives_row['n_display']}")
    control_single_lines = []
    for legacy_burden_descriptives__, legacy_burden_descriptives_row in control_breakdown_display.iterrows():
        control_single_lines.append(f"  - {legacy_burden_descriptives_row['category']}: n={legacy_burden_descriptives_row['n_display']}")

    legacy_burden_descriptives_report = f"""# Step 15 | Burden_rare (10-gene main analysis score) Exploratory distribution analysis summary report

    Corresponding script: `scripts/10_exploratory_age_ehr.py`. **This analysis is exploratory/descriptive analysis,
    Not an official Aim 4 Cox model** (same property as `scripts/10_exploratory_age_ehr.py`, for intuition checking).

    Main score: `burden_rare_main_ptv_10gene` (main_ptv_maf1pct mask, 10 gene panel, binary carrier ×
    Sum of log(OR), Σ_g ln(OR_g) × carrier_ig). Cohort caliber: group in {{case,control}}& pass_qc==1，
    total{n_case + n_control}people(case={n_case}, control={n_control}), consistent with script 14.

    ---

    ## 1. Case vs Control distribution comparison

    About 98% of people have Burden_rare exactly 0 (do not carry any qualifying PTV in the 10 panel genes), so the full distribution is directly
    There is no information in making a boxplot. This section is divided into three steps to check:

    ### 1a. Carrier rate (Burden_rare > 0) whether > 0 itself

    | group | n | n_carrier | carrier rate | 95% Wilson CI |
    |---|---|---|---|---|
    | case |{n_case} | {n_case_carrier} | {legacy_burden_descriptives_p_case*100:.3f}% | {ci_case[0]*100:.3f}%-{ci_case[1]*100:.3f}% |
    | control | {n_control} | {n_control_carrier} | {p_control*100:.3f}% | {ci_control[0]*100:.3f}%-{ci_control[1]*100:.3f}% |

    - Fisher exact test: OR(case carrier vs control carrier) = **{odds_ratio:.3f}** (95% CI {or_ci_low:.3f}-{or_ci_high:.3f}),
      p = **{fmt_p(fisher_p)}**
    - Two-proportions z-test: z ={z_stat:.3f}, p = {fmt_p(z_p)}- Conclusion: The probability of case carrying qualifying PTV of any panel gene is approximately that of control**{odds_ratio:.2f}Times**, the difference is extremely significant.

    ### 1b. Comparison of score distribution within the Nonzero subset (comparison of weights between carriers)

    | group | n(carriers) | mean | median | sd | q25 | q75 |
    |---|---|---|---|---|---|---|
    | case |{len(case_nz)} | {case_nz.mean():.4f} | {case_nz.median():.4f} | {case_nz.std():.4f} | {case_nz.quantile(0.25):.4f} | {case_nz.quantile(0.75):.4f} |
    | control | {len(control_nz)} | {control_nz.mean():.4f} | {control_nz.median():.4f} | {control_nz.std():.4f} | {control_nz.quantile(0.25):.4f} | {control_nz.quantile(0.75):.4f} |

    Mann-Whitney U (nonzero subset only, case vs control): U = {mw_nz_stat:.1f}, p = **{fmt_p(mw_nz_p)}**

    Conclusion: Even if the comparison is limited to those who are \"already carriers\", the score of case is still significantly higher than that of control - indicating that case carrier
    Not only is the \"carrying rate higher\", the carried genes are also more likely to have high log(OR) genes (such as BRCA1/BRCA2/PALB2), rather than just low-weight genes.
    (such as RAD51C/RAD51D) the carriage rate happens to be increased.

    ### 1c. Mann-Whitney U test for the full distribution (including 0) (a more technically rigorous overall test under zero-inflated distribution)

    Mann-Whitney U (full distribution incl. zeros, case vs control): U ={mw_full_stat:.1f}, p = **{fmt_p(mw_full_p)}**

    The three test directions are completely consistent, and they all point to the case\'s Burden_rare being significantly higher than the control.

    ### 1d. Whether it meets biological expectations

    **Conforms. ** Carriers of rare pathogenic PTV should be enriched in the case (carrier rate OR ={odds_ratio:.2f}), and high penetrance genes
    (BRCA1 log(OR)={np.log(gene_or_main['BRCA1']):.3f}, BRCA2 log(OR)={np.log(gene_or_main['BRCA2']):.3f},
    PALB2 log(OR)={np.log(gene_or_main['PALB2']):.3f}) should be compared with low-penetrance genes
    (RAD51C log(OR)={np.log(gene_or_main['RAD51C']):.3f}, RAD51D log(OR)={np.log(gene_or_main['RAD51D']):.3f})
    More concentrated in the case - this is the actual meaning of the superposition of the two results of 1a \"higher carrier rate\" and 1b \"higher score within the carrier\".
    The gene attribution distributions in Section 2/4 below can further demonstrate this directly.

    Figure: `results/figures/15_burden_rare_case_control_comparison.png`
    Table: `results/table/15_burden_rare_carrier_rate_by_group.csv`,
    `results/table/15_burden_rare_carrier_rate_test.csv`,
    `results/table/15_burden_rare_nonzero_subset_summary.csv`,
    `results/table/15_burden_rare_distribution_tests.csv`

    ---

    ## 2. Inside Case: carrier’s genetic attribution distribution (n={n_case_carriers_chk})

    Attribution is based on the 0/1 carrier column of 10 panel genes: the single-gene carrier is directly attributed to the gene; more than 2 are carried at the same time
    People whose panel genes qualify PTV (compound carriers) are classified into one category separately without further gene subdivision (the sample is too small, and
    \"Which specific genes are these?\" There are many combinations and each grid is far less than 20, so it cannot be disclosed safely).

    Gene-by-gene attribution counts are not disclosed in this report: gene-by-gene carrier counts are only in `results/table/19_gene_tier_carrier_counts.csv`
    (Aim 1 modeling cohort) This specification table is disclosed. Subtracting the gene-by-gene number of this step on the 143,255-person QC cohort will result in small grids.
    See `data/15_internal/` for internal versions.

    Conclusion (qualitative): The genes with the largest number of case carriers are BRCA2, BRCA1, CHEK2, and ATM, among which BRCA1/BRCA2/PALB2
    The share of high-penetrance genes in case is significantly higher than its expected share under the basic prevalence of the entire population (case:control≈1:15.6).
    It is mutually corroborated with the conclusion of Section 1.

    ---

    ## 3. Inside Case: Association between Burden_rare and age_at_diagnosis

    Limit the case and age_at_diagnosis is not missing ({n_case_age_missing}Missing people have been eliminated, and the analysis sample n={len(case_age_df)}).

    ### 3a. Age_at_diagnosis comparison of Carrier vs non-carrier

    | carrier_status | n | mean | median | sd | q25 | q75 |
    |---|---|---|---|---|---|---|
    | carrier (Burden_rare>0) |{len(age_carrier)} | {age_carrier.mean():.2f} | {age_carrier.median():.2f} | {age_carrier.std():.2f} | {age_carrier.quantile(0.25):.2f} | {age_carrier.quantile(0.75):.2f} |
    | non-carrier (Burden_rare==0) | {len(age_noncarrier)} | {age_noncarrier.mean():.2f} | {age_noncarrier.median():.2f} | {age_noncarrier.std():.2f} | {age_noncarrier.quantile(0.25):.2f} | {age_noncarrier.quantile(0.75):.2f} |

    Mann-Whitney U: U = {mw_age_stat:.1f}, p = **{fmt_p(mw_age_p)}**

    **Direction**: Carrier group age_at_diagnosis median ({age_carrier.median():.2f}age){direction}In the non-carrier group
    ({age_noncarrier.median():.2f}age).

    ### 3b. Carrier internal: score (continuous) is related to Spearman of age_at_diagnosis

    n={len(age_carrier_only_df)}(Within the carrier, the sample size is small and no further subdivision and stratification such as ten-year periods will be performed. Please refer to the coordinator instructions).

    Spearman rho = **{legacy_burden_descriptives_rho:.4f}**, p = **{fmt_p(legacy_burden_descriptives_rho_p)}**

    ### 3c. Whether it meets the core hypothesis (README.md Section II: the higher the load, the earlier the onset){'In line with direction expectations: carriers develop disease earlier than non-carriers, and the p value' + ('Significantly' if mw_age_p < 0.05 else 'Not significant') + "。" if age_carrier.median() < age_noncarrier.median() else 'The direction is opposite or not obvious from the core assumption:'}Age difference test p= for carrier vs non-carrier{fmt_p(mw_age_p)}
    {'(Statistically significant, supporting the direction of the H3 early onset hypothesis: patients carrying rare PTV with panel genes have earlier onset)' if mw_age_p < 0.05 else f"(Not statistically significant, with this sample size (n_carrier={len(age_carrier)}) is not enough to confirm the age difference within the carrier, please note that the conclusion may be different from the Aim 4 formal model)"};
    Carrier internal score-Spearman correlation of age rho={legacy_burden_descriptives_rho:.4f}, p={fmt_p(legacy_burden_descriptives_rho_p)}
    {('also points to' + ('The higher the load, the earlier the onset of symptoms' if legacy_burden_descriptives_rho < 0 else '') + (', and reach statistical significance' if legacy_burden_descriptives_rho_p < 0.05 else ', but not statistically significant') + "。") if not np.isnan(legacy_burden_descriptives_rho) else ""}**Note**: This section is only a descriptive exploration and does not control any covariates (such as ancestral principal components, PRS_common, etc.). The conclusion cannot replace Aim 4
    Formal Cox model (`logit/Cox(cancer/onset) ~ PRS_common + Burden_rare + interaction term + covariates`).

    Figure: `results/figures/15_burden_rare_age_at_diagnosis_by_carrier_status.png`
    Table: `results/table/15_burden_rare_age_at_diagnosis_by_carrier_status.csv`,
    `results/table/15_burden_rare_age_at_diagnosis_tests.csv`

    ---

    ## 4. Internal Control: Gene attribution distribution of carrier (n={n_control_carriers_chk})

    Same caliber as section 2 (gene-by-gene counting only in internal version `data/15_internal/`).

    Conclusion (qualitative): Compared with the case internal attribution distribution (Section 2), the genes with the largest share in the control carrier are relatively low/medium genes such as CHEK2 and ATM.
    Penetrance genes, and the relative share of BRCA1 in the control carrier is significantly lower than its share in the case carrier - that is, control
    The genetic makeup of carriers is generally more biased toward low-penetrance genes, which is consistent with Section 1 \"case not only has a higher carrier rate, but also has a higher internal score of carriers.\"
    \"Also higher\" conclusion is consistent and explainable at the genetic level.

    ---

    ## Output file list

    - Table: `results/table/15_burden_rare_carrier_rate_by_group.csv`,
      `results/table/15_burden_rare_carrier_rate_test.csv`,
      `results/table/15_burden_rare_nonzero_subset_summary.csv`,
      `results/table/15_burden_rare_distribution_tests.csv`,
      `results/table/15_burden_rare_age_at_diagnosis_by_carrier_status.csv`,
      `results/table/15_burden_rare_age_at_diagnosis_tests.csv`
    - Figure: `results/figures/15_burden_rare_case_control_comparison.png`,
      `results/figures/15_burden_rare_age_at_diagnosis_by_carrier_status.png`"""

    with open(os.path.join(LEGACY_BURDEN_DESCRIPTIVES_REPORT_DIR, "15_burden_rare_distributional_analysis_summary.md"), "w") as fh:
        fh.write(legacy_burden_descriptives_report)
    legacy_burden_descriptives_log(f"Has written:{os.path.join(LEGACY_BURDEN_DESCRIPTIVES_REPORT_DIR, '15_burden_rare_distributional_analysis_summary.md')}")

    legacy_burden_descriptives_log('Script 15 is executed.')

# ============================================================================
# 6. legacy-group-contrasts: analyze_legacy_age_contrasts
# ============================================================================

import os
import warnings
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from matplotlib.patches import Patch

def legacy_age_contrasts_log(msg):
    print(msg, flush=True)

def legacy_age_contrasts_section(title):
    legacy_age_contrasts_log("\n" + "=" * 84)
    legacy_age_contrasts_log(title)
    legacy_age_contrasts_log("=" * 84)

def legacy_age_contrasts_suppress(n, min_cell=20):
    if n is None or (isinstance(n, float) and np.isnan(n)):
        return n
    return int(n) if n >= min_cell else f"<{min_cell}"

def bucket(n, min_cell=20):
    """Bucket."""
    if n is None or (isinstance(n, float) and np.isnan(n)):
        return n
    n = int(n)
    if n >= min_cell:
        return str(n)
    return "<5" if n < 5 else "5-19"

def sanitize_for_export(df, count_col, derived_cols=(), min_cell=20):
    """Sanitize for export."""
    out = df.copy()
    small = pd.to_numeric(out[count_col], errors="coerce") < min_cell
    for c in derived_cols:
        if c in out.columns:
            out.loc[small, c] = np.nan
    return out.drop(columns=[count_col])

def legacy_age_contrasts_wilson_ci(k, n, z=1.959963984540054):
    """Legacy age contrasts wilson ci."""
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    denom = 1 + z ** 2 / n
    centre = (p + z ** 2 / (2 * n)) / denom
    half = z * np.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / denom
    return (centre - half, centre + half)

def zscore(s):
    mu, sd = s.mean(), s.std()
    return (s - mu) / sd, mu, sd

def smd(x, y):
    """Smd."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    nx, ny = len(x), len(y)
    sp = np.sqrt(((nx - 1) * x.var(ddof=1) + (ny - 1) * y.var(ddof=1)) / (nx + ny - 2))
    return (x.mean() - y.mean()) / sp if sp > 0 else np.nan

def fit_ols_hc1(y, X):
    return sm.OLS(np.asarray(y, float), X).fit(cov_type="HC1")

def fit_logit_hc1(y, X):
    return sm.GLM(np.asarray(y, float), X, family=sm.families.Binomial()).fit(cov_type="HC1")

def coef_rows(fit, label, n, terms=None, exponentiate=False):
    rows = []
    idx = fit.params.index if terms is None else [t for t in terms if t in fit.params.index]
    for t in idx:
        b, se, p = fit.params[t], fit.bse[t], fit.pvalues[t]
        r = {"model": label, "term": t, "beta": b, "robust_se_HC1": se,
             "ci_low": b - Z95 * se, "ci_high": b + Z95 * se, "p_value": p, "n": n}
        if exponentiate:
            r["OR"] = np.exp(b)
            r["OR_ci_low"] = np.exp(b - Z95 * se)
            r["OR_ci_high"] = np.exp(b + Z95 * se)
        rows.append(r)
    return pd.DataFrame(rows)

def build_X(df_in, cols):
    X = df_in[cols].astype(float).copy()
    return sm.add_constant(X, has_constant="add")

def chi2_2x2(k1, n1, k2, n2):
    tbl = np.array([[k1, n1 - k1], [k2, n2 - k2]])
    chi2, p, dof, _ = stats.chi2_contingency(tbl, correction=False)
    return chi2, p, dof

def mean_p(a):
    return float(np.mean(1.0 / (1.0 + np.exp(-(a + lp0)))))

def mix_cat(r):
    if r["prs_high"] == 1 and r["carrier_main_10gene"] == 1:
        return "both (PRS-high AND carrier)"
    if r["prs_high"] == 1:
        return "PRS-high only"
    if r["carrier_main_10gene"] == 1:
        return "carrier only"
    return "neither"

def make_env_X(df_in, cat_cols, extra_terms):
    parts = [df_in[extra_terms + PC_COLS].astype(float)]
    if "bmi_value" in cat_cols:
        parts.append(df_in[["bmi_value"]].astype(float))
    for c in [c for c in cat_cols if c != "bmi_value"]:
        parts.append(pd.get_dummies(df_in[c], prefix=c, drop_first=True).astype(float))
    X = pd.concat(parts, axis=1)
    return sm.add_constant(X, has_constant="add")

def _lab(a, b):
    return f"<{b}" if a == 0 else (f">={a}" if b == 200 else f"{a}-{b - 1}")

def fmt(x, nd=3):
    return "—" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:.{nd}f}"

def pfmt(p):
    return "<1e-300" if p == 0 else f"{p:.3g}"

def _pct_cell(r):
    """ pct cell."""
    if pd.isna(r["pct_of_group"]):
        return "suppressed | suppressed"
    return (f"{r['pct_of_group']:.2f} | "
            f"[{r['pct_ci_low']:.2f}, {r['pct_ci_high']:.2f}]")

def analyze_legacy_age_contrasts():
    """Initialize runtime state and execute this component explicitly."""
    global A2, A_gene, A_gene_rows, A_rows, A_tbl, A_test_rows, A_tests, LEGACY_AGE_CONTRASTS_BASE, B_rows, B_tbl, B_test_tbl, B_tests, CUTOFFS, C_export, C_rows, C_tbl, C_test_tbl, C_tests, LEGACY_AGE_CONTRASTS_DATA_DIR, LEGACY_AGE_CONTRASTS_DISCLOSURE_MIN_CELL, D_rows, D_tbl, EXPECTED_AIM4_N, EXPECTED_FULL_N, E_rows, E_tbl, LEGACY_AGE_CONTRASTS_FIG_DIR, GENE_COLS, MIX_ORDER, N_SIM, PC_COLS, LEGACY_AGE_CONTRASTS_REPORT_DIR, SIM_SEED, LEGACY_AGE_CONTRASTS_TABLE_DIR, X, X2, X3, Xe, Xn, Xs, Z95, legacy_age_contrasts__, _PARTNER, _small, _small2, legacy_age_contrasts_a, a2, a_md, agg, agg_out, agg_p, alpha, atest_md, legacy_age_contrasts_ax, legacy_age_contrasts_b, b_bur, b_inter_d1, b_md, b_prs, blab, bottom, bottoms, btest_md, legacy_age_contrasts_burden, legacy_age_contrasts_c, c2, c_md, carr, cases, cases_e, cat, cc, cell_stats, ch, chi2, cis, cl, clinical, coef16, colors, cols, ct, ctest_md, cut, legacy_age_contrasts_d, d1, d2, d_key, d_md, dec, dof, e, e_key, e_md, edges, en, env_cols, er, errs, f, f3, fe, legacy_age_contrasts_fig, fig_files, fit_d1, fit_d2, fitres, fl, flow_md, flow_rows, flow_tbl, fn, fs, full, full_g, legacy_age_contrasts_g, garr, gene_carrier, gene_map, group_defs, legacy_age_contrasts_grp, h, handles, hi_a, hi_ci, hi_i, hi_q, hi_thr, hi_thr_d, legacy_age_contrasts_i, idx, is_hi, is_lo, is_merged, j, k, k1, k2, kb_e, kb_l, ke, keys, kg, kl, l, lab, lab_e, lab_l, labels, labs, lo_a, lo_ci, lo_i, lo_q, lo_thr, lo_thr_d, lp0, m, mat, mde, merged, merged_edges, mid, miss_md, miss_rows, miss_tbl, mu_b_c, mu_b_f, mu_p_c, mu_p_f, n, n1, n2, n_c_small, n_gene_suppressed, n_na, n_tot, nb_e, nb_l, nch, ncl, ne, need, nh, nl, obs_diff, obs_prev, obs_row_h, obs_row_l, obs_rr, ok, ok_all, p, p2, p_b, legacy_age_contrasts_p_case, p_chi, p_inter_d1, p_p, partner, pc, pct_zero_burden, pg, legacy_age_contrasts_pheno, plot_groups, pow_md, pow_rows, pow_tbl, prelim, legacy_age_contrasts_prs, prs_grp, pt_b, pt_p, q_md, legacy_age_contrasts_qc, quin, r, r1, r2, r3, rates, raw, legacy_age_contrasts_report, report_path, rh, rl, rn, rng, legacy_age_contrasts_row, s1, s2, s3, s4, sch, scheme, sd_b_c, sd_b_f, sd_p_c, sd_p_f, se, se_d, se_inter_d1, sel, sim_ctrl_diff, sim_diff, sim_hi_rate, sim_lo_rate, sim_ncase, sim_p_two, sim_rr, smd_bur, smd_prs, specs, split, strat_rows, legacy_age_contrasts_sub, sub_c, t, t45, t_b, t_p, term, u_b, u_p, v, vals, verdict_null, verdict_sentence, which, xlabels, y, y_bin, z_a, z_b

    matplotlib.use("Agg")
    plt.style.use("default")
    warnings.filterwarnings("ignore", category=FutureWarning)

    LEGACY_AGE_CONTRASTS_BASE = "/home/jupyter/workspace/analysis"
    LEGACY_AGE_CONTRASTS_DATA_DIR = os.path.join(LEGACY_AGE_CONTRASTS_BASE, "data")
    LEGACY_AGE_CONTRASTS_FIG_DIR = os.path.join(LEGACY_AGE_CONTRASTS_BASE, "results", "figures")
    LEGACY_AGE_CONTRASTS_TABLE_DIR = os.path.join(LEGACY_AGE_CONTRASTS_BASE, "results", "table")
    LEGACY_AGE_CONTRASTS_REPORT_DIR = os.path.join(LEGACY_AGE_CONTRASTS_BASE, "results", "report")
    for legacy_age_contrasts_d in (LEGACY_AGE_CONTRASTS_FIG_DIR, LEGACY_AGE_CONTRASTS_TABLE_DIR, LEGACY_AGE_CONTRASTS_REPORT_DIR):
        os.makedirs(legacy_age_contrasts_d, exist_ok=True)

    PC_COLS = [f"PC{legacy_age_contrasts_i}" for legacy_age_contrasts_i in range(1, 11)]
    GENE_COLS = ["ATM", "BARD1", "BRCA1", "BRCA2", "CHEK2", "PALB2",
                 "RAD51C", "RAD51D", "TP53", "MAP3K1"]
    Z95 = 1.959963984540054
    LEGACY_AGE_CONTRASTS_DISCLOSURE_MIN_CELL = 20
    N_SIM = 2000
    SIM_SEED = 20260819
    CUTOFFS = [45]  
    EXPECTED_AIM4_N = 7996   
    EXPECTED_FULL_N = 141204  

    # ===========================================================================
    
    # ===========================================================================
    legacy_age_contrasts_section('Step 0: PRE-STATED, before any results')
    legacy_age_contrasts_log(f"""Aim 4 cases queue n={EXPECTED_AIM4_N:,}, of which the 10-gene panel PTV carrier is expected to be only about 6% (the self-reported white cohort measured 6.2%).
      The SE magnitude of the PRS slope in the carrier layer is about 0.5 years/SD, and the lower detectable limit of the PRS x carrier interaction term is therefore very high.
      This script will use the actual fitted interaction term SE to calculate the smallest detectable effect at 80% power and two-sided alpha=0.05
      MDE = (z_0.975 + z_0.80) * SE = 2.8016 * SE.
      **The interaction term is not significant == \"Insufficient power to rule out interaction\" does not mean \"no interaction\". **
      Also remind: All comparisons within the case are affected by index-event/collider bias (see script header).""")

    # ===========================================================================
    
    # ===========================================================================
    legacy_age_contrasts_section(f"Step A0: Read data and build the full queue (n={EXPECTED_FULL_N:,}) with Aim4 cases queue(n={EXPECTED_AIM4_N:,})")

    legacy_age_contrasts_pheno = pd.read_csv(os.path.join(LEGACY_AGE_CONTRASTS_DATA_DIR, "04_phenotype_diagnostics_person_level.csv.gz"))
    legacy_age_contrasts_qc = pd.read_csv(os.path.join(LEGACY_AGE_CONTRASTS_DATA_DIR, "05_cohort_qc_flags.csv"))
    pc = pd.read_csv(os.path.join(LEGACY_AGE_CONTRASTS_DATA_DIR, "09_ancestry_pc_person_level.csv.gz"))
    legacy_age_contrasts_burden = pd.read_csv(os.path.join(LEGACY_AGE_CONTRASTS_DATA_DIR, "14_burden_rare_scores_person_level.csv.gz"))
    legacy_age_contrasts_prs = pd.read_csv(os.path.join(LEGACY_AGE_CONTRASTS_DATA_DIR, "prs_scores.csv"))
    split = pd.read_csv(os.path.join(LEGACY_AGE_CONTRASTS_DATA_DIR, "16_joint_model_cv_oof_person_level.csv.gz"))  
    clinical = pd.read_csv(os.path.join(LEGACY_AGE_CONTRASTS_DATA_DIR, "09_extended_clinical_person_level.csv.gz"))
    gene_carrier = pd.read_csv(os.path.join(LEGACY_AGE_CONTRASTS_DATA_DIR, "13_carrier_status_main_ptv_maf1pct_person_level.csv.gz"))

    full = split[["person_id", "group", "carrier_main_10gene"]].copy()
    full = full.merge(legacy_age_contrasts_prs[["person_id", "PRS311"]], on="person_id", how="left")
    full = full.merge(legacy_age_contrasts_burden[["person_id", "burden_rare_main_ptv_10gene"]], on="person_id", how="left")
    full = full.merge(pc[["person_id"] + PC_COLS], on="person_id", how="left")
    legacy_age_contrasts_log(f"Full queue n={len(full)} (case={int((full['group']=='case').sum())}, "
        f"control={int((full['group']=='control').sum())})")
    assert len(full) == EXPECTED_FULL_N, f"Full queue n={len(full)} != {EXPECTED_FULL_N}"
    for legacy_age_contrasts_c in ["PRS311", "burden_rare_main_ptv_10gene", "carrier_main_10gene"]:
        assert full[legacy_age_contrasts_c].isna().sum() == 0, f"{legacy_age_contrasts_c}There are NAs in the whole queue"

    full["PRS311_z_full"], mu_p_f, sd_p_f = zscore(full["PRS311"])
    full["burden_z_full"], mu_b_f, sd_b_f = zscore(full["burden_rare_main_ptv_10gene"])
    legacy_age_contrasts_log(f"[z-normalized@full cohort] PRS311 mean={mu_p_f:.5f} sd={sd_p_f:.5f}; "
        f"burden mean={mu_b_f:.6f} sd={sd_b_f:.6f}")

    flow_rows = [{"step": '0. Overall phenotype table', "n": len(legacy_age_contrasts_pheno)}]
    s1 = legacy_age_contrasts_pheno[legacy_age_contrasts_pheno["group"] == "case"].copy()
    flow_rows.append({"step": "1. group=='case'", "n": len(s1)})
    s2 = s1.merge(legacy_age_contrasts_qc[["person_id", "pass_qc"]], on="person_id", how="inner")
    s2 = s2[s2["pass_qc"] == 1].copy()
    flow_rows.append({"step": "2. & pass_qc==1", "n": len(s2)})
    s3 = s2.dropna(subset=["age_at_diagnosis"]).copy()
    flow_rows.append({"step": '3. & age_at_diagnosis is not missing', "n": len(s3)})
    s4 = s3[s3["bc_is_first_record"] != True].copy()  # noqa: E712
    flow_rows.append({"step": "4. & bc_is_first_record != True", "n": len(s4)})

    cases = s4.merge(legacy_age_contrasts_prs[["person_id", "PRS311"]], on="person_id", how="left")
    cases = cases.merge(legacy_age_contrasts_burden, on="person_id", how="left")
    cases = cases.merge(pc[["person_id"] + PC_COLS], on="person_id", how="left")
    need = ["PRS311", "burden_rare_main_ptv_10gene", "burden_rare_main_ptv_9gene"] + PC_COLS
    cases = cases.dropna(subset=need).reset_index(drop=True)
    flow_rows.append({"step": '5. & PRS311/Burden_rare/PC1-10 is not missing (complete-case)', "n": len(cases)})

    flow_tbl = pd.DataFrame(flow_rows)
    flow_tbl["n_disclosed"] = flow_tbl["n"].apply(legacy_age_contrasts_suppress)
    flow_tbl.to_csv(os.path.join(LEGACY_AGE_CONTRASTS_TABLE_DIR, "18_cohort_flow.csv"), index=False)
    legacy_age_contrasts_log('Aim 4 cases queue construction process:')
    legacy_age_contrasts_log(flow_tbl.to_string(index=False))
    assert len(cases) == EXPECTED_AIM4_N, (
        f"Aim4 cases queue n={len(cases)}, with n= from script 17{EXPECTED_AIM4_N}inconsistent,"
        'The queue building logic has deviated and is stopped.')
    legacy_age_contrasts_log(f"ASSERT PASS: Aim4 cases queue n={len(cases)}== script 17 of{EXPECTED_AIM4_N}")

    cases["carrier_main_10gene"] = (cases["burden_rare_main_ptv_10gene"] > 0).astype(int)
    
    cases["PRS311_z"], mu_p_c, sd_p_c = zscore(cases["PRS311"])
    cases["burden_z"], mu_b_c, sd_b_c = zscore(cases["burden_rare_main_ptv_10gene"])
    legacy_age_contrasts_log(f"[z-normalized @Aim4 cases queue n={EXPECTED_AIM4_N:,}] PRS311 mean={mu_p_c:.5f} sd={sd_p_c:.5f}; "
        f"burden mean={mu_b_c:.6f} sd={sd_b_c:.6f}")

    cases["PRS311_z_fullref"] = (cases["PRS311"] - mu_p_f) / sd_p_f
    cases["burden_z_fullref"] = (cases["burden_rare_main_ptv_10gene"] - mu_b_f) / sd_b_f
    legacy_age_contrasts_log(f"[z-normalized @full queue reference *_z_fullref] cases queue overall PRS311_z_fullref mean="
        f"{cases['PRS311_z_fullref'].mean():.4f}, burden_z_fullrefmean={cases['burden_z_fullref'].mean():.4f}")
    pct_zero_burden = float((cases["burden_rare_main_ptv_10gene"] == 0).mean() * 100)
    legacy_age_contrasts_log(f"The proportion of burden==0 in the cases queue:{pct_zero_burden:.2f}% "
        f"(carrier rate{100-pct_zero_burden:.2f}%, n_carrier={int(cases['carrier_main_10gene'].sum())})")

    # ===========================================================================
    
    # ===========================================================================
    legacy_age_contrasts_section(f"Step A: 2x2 discriminant group comparison (PRS311 extreme group x case/control, full queue n={EXPECTED_FULL_N:,})")

    gene_map = gene_carrier.set_index("person_id")
    full_g = full.merge(gene_carrier, on="person_id", how="left")

    A_rows, A_test_rows, A_gene_rows = [], [], []
    for scheme, lo_q, hi_q in [('decile (10/90, main analysis)', 0.10, 0.90),
                               ('quintile (20/80, sensitivity)', 0.20, 0.80)]:
        lo_thr = full["PRS311"].quantile(lo_q)
        hi_thr = full["PRS311"].quantile(hi_q)
        legacy_age_contrasts_log(f"\n  --- {scheme}: PRS311 cut point low<={lo_thr:.4f}, high>={hi_thr:.4f}(The cut point is calculated on the entire queue) ---")
        prs_grp = np.where(full_g["PRS311"] >= hi_thr, "high",
                           np.where(full_g["PRS311"] <= lo_thr, "low", "mid"))
        full_g["_prs_grp"] = prs_grp
        cell_stats = {}
        for legacy_age_contrasts_grp in ["case", "control"]:
            for pg in ["high", "low"]:
                legacy_age_contrasts_sub = full_g[(full_g["group"] == legacy_age_contrasts_grp) & (full_g["_prs_grp"] == pg)]
                n = len(legacy_age_contrasts_sub)
                k = int(legacy_age_contrasts_sub["carrier_main_10gene"].sum())
                lo_ci, hi_ci = legacy_age_contrasts_wilson_ci(k, n)
                cell_stats[(legacy_age_contrasts_grp, pg)] = (k, n)
                A_rows.append({
                    "scheme": scheme, "group": legacy_age_contrasts_grp, "prs_group": pg,
                    "n": n, "n_disclosed": legacy_age_contrasts_suppress(n),
                    "n_carrier": k, "n_carrier_disclosed": legacy_age_contrasts_suppress(k),
                    "carrier_rate_pct": 100 * k / n if n else np.nan,
                    "carrier_rate_ci_low_pct": 100 * lo_ci, "carrier_rate_ci_high_pct": 100 * hi_ci,
                    "burden_mean": legacy_age_contrasts_sub["burden_rare_main_ptv_10gene"].mean(),
                    "burden_median": legacy_age_contrasts_sub["burden_rare_main_ptv_10gene"].median(),
                    "burden_mean_among_carriers": legacy_age_contrasts_sub.loc[legacy_age_contrasts_sub["carrier_main_10gene"] == 1,
                                                          "burden_rare_main_ptv_10gene"].mean(),
                    "prs311_mean_z_full": legacy_age_contrasts_sub["PRS311_z_full"].mean(),
                })
                
                if scheme.startswith("decile"):
                    for legacy_age_contrasts_g in GENE_COLS:
                        kg = int(legacy_age_contrasts_sub[legacy_age_contrasts_g].sum()) if legacy_age_contrasts_g in legacy_age_contrasts_sub.columns else np.nan
                        A_gene_rows.append({
                            "scheme": scheme, "group": legacy_age_contrasts_grp, "prs_group": pg, "gene": legacy_age_contrasts_g,
                            "n_cell": n, "n_carrier_gene": kg,
                            "n_carrier_gene_disclosed": legacy_age_contrasts_suppress(kg),
                            "carrier_rate_pct": 100 * kg / n if n else np.nan,
                        })
        for legacy_age_contrasts_grp in ["case", "control"]:
            k1, n1 = cell_stats[(legacy_age_contrasts_grp, "high")]
            k2, n2 = cell_stats[(legacy_age_contrasts_grp, "low")]
            chi2, p, dof = chi2_2x2(k1, n1, k2, n2)
            r1, r2 = k1 / n1, k2 / n2
            
            se_d = np.sqrt(r1 * (1 - r1) / n1 + r2 * (1 - r2) / n2)
            legacy_age_contrasts_d = r1 - r2
            A_test_rows.append({
                "scheme": scheme, "contrast": f"{legacy_age_contrasts_grp}: PRS-high vs PRS-low (carrier rate)",
                "n_high": n1, "carrier_rate_high_pct": 100 * r1,
                "n_low": n2, "carrier_rate_low_pct": 100 * r2,
                "diff_pct_points": 100 * legacy_age_contrasts_d,
                "diff_ci_low_pct": 100 * (legacy_age_contrasts_d - Z95 * se_d), "diff_ci_high_pct": 100 * (legacy_age_contrasts_d + Z95 * se_d),
                "rate_ratio_high_over_low": r1 / r2 if r2 > 0 else np.nan,
                "chi2": chi2, "df": dof, "p_value": p,
            })
            legacy_age_contrasts_log(f"    [{legacy_age_contrasts_grp}] high {k1}/{n1}={100*r1:.2f}% vs low {k2}/{n2}={100*r2:.2f}%  "
                f"diff={100*legacy_age_contrasts_d:+.2f}pp  RR={r1/r2:.3f}  chi2={chi2:.2f}, p={p:.3g}")

    A_tbl = pd.DataFrame(A_rows)
    A_tests = pd.DataFrame(A_test_rows)
    A_gene = pd.DataFrame(A_gene_rows)
    A_tbl.to_csv(os.path.join(LEGACY_AGE_CONTRASTS_TABLE_DIR, "18_A_2x2_prs_extreme_cells.csv"), index=False)
    A_tests.to_csv(os.path.join(LEGACY_AGE_CONTRASTS_TABLE_DIR, "18_A_2x2_carrier_rate_tests.csv"), index=False)
    sanitize_for_export(A_gene, "n_carrier_gene", ["carrier_rate_pct"]).to_csv(
        os.path.join(LEGACY_AGE_CONTRASTS_TABLE_DIR, "18_A_2x2_gene_level_breakdown.csv"), index=False)
    legacy_age_contrasts_log('Summary of each grid of 2x2:')
    legacy_age_contrasts_log(A_tbl[["scheme", "group", "prs_group", "n", "n_carrier", "carrier_rate_pct",
               "carrier_rate_ci_low_pct", "carrier_rate_ci_high_pct",
               "burden_mean", "burden_median"]].to_string(index=False))
    n_gene_suppressed = int((A_gene["n_carrier_gene"] < LEGACY_AGE_CONTRASTS_DISCLOSURE_MIN_CELL).sum())
    legacy_age_contrasts_log(f"\n sub-gene breakdown:{len(A_gene)}grid{n_gene_suppressed}Cell n<20 has been marked as \'<20\' under the disclosure rules"
        f"(The original exact count column is not written to the external table; the derived carrier_rate_pct is left blank in these cells to prevent the pct x denominator from being inverted)")

    prelim = {("case", "high"): (1606, 4.42), ("case", "low"): (417, 10.55),
              ("control", "high"): (12835, 1.74), ("control", "low"): (14024, 1.87)}
    legacy_age_contrasts_log('Consistency check with expected value (decile):')
    ok_all = True
    dec = A_tbl[A_tbl["scheme"].str.startswith("decile")]
    for (legacy_age_contrasts_g, pg), (en, er) in prelim.items():
        legacy_age_contrasts_row = dec[(dec["group"] == legacy_age_contrasts_g) & (dec["prs_group"] == pg)].iloc[0]
        ok = (legacy_age_contrasts_row["n"] == en) and (abs(legacy_age_contrasts_row["carrier_rate_pct"] - er) < 0.02)
        ok_all &= ok
        legacy_age_contrasts_log(f"    {legacy_age_contrasts_g}-{pg}: n={legacy_age_contrasts_row['n']}(expect{en}), carrier={legacy_age_contrasts_row['carrier_rate_pct']:.2f}% "
            f"(expect{er}%) -> {'MATCH' if ok else 'MISMATCH'}")
    legacy_age_contrasts_log(f"Expected numerical reproduction:{'All consistent' if ok_all else 'There are inconsistencies that require manual verification.'}")

    # ===========================================================================
    
    # ===========================================================================
    legacy_age_contrasts_section('Step A2: Additive (no interaction) zero model simulation benchmark - determine whether the within-case comparison of A is just a selection bias')

    coef16 = pd.read_csv(os.path.join(LEGACY_AGE_CONTRASTS_TABLE_DIR, "16_joint_model_coefficients_primary.csv"))
    b_prs = float(coef16.loc[coef16["term"] == "PRS311_z", "beta"].iloc[0])
    b_bur = float(coef16.loc[coef16["term"] == "burden_z", "beta"].iloc[0])
    legacy_age_contrasts_log(f"Get the Aim 1 main model main effects (results/table/16_joint_model_coefficients_primary.csv):")
    legacy_age_contrasts_log(f"    beta_PRS311_z = {b_prs:.4f},  beta_burden_z = {b_bur:.4f}(Interaction terms are always set to 0 = zero model)")

    obs_prev = float((full["group"] == "case").mean())
    lp0 = b_prs * full["PRS311_z_full"].values + b_bur * full["burden_z_full"].values

    lo_a, hi_a = -20.0, 20.0
    for legacy_age_contrasts__ in range(200):
        mid = 0.5 * (lo_a + hi_a)
        if mean_p(mid) < obs_prev:
            lo_a = mid
        else:
            hi_a = mid
    alpha = 0.5 * (lo_a + hi_a)
    legacy_age_contrasts_p_case = 1.0 / (1.0 + np.exp(-(alpha + lp0)))
    legacy_age_contrasts_log(f"Simulation settings: logistic (liability on logit scale),"
        f"logit P(case) = alpha + {b_prs:.4f}*PRS_z + {b_bur:.4f}*burden_z, **No interaction terms**;"
        f"PRS_z/burden_z directly takes the individual value of the real observation (retaining the real marginal distribution and the real joint structure of the two),"
        f"The case state is redrawn according to Bernoulli(p_i).")
    legacy_age_contrasts_log(f"  alpha={alpha:.5f}, simulated average prevalence ={legacy_age_contrasts_p_case.mean()*100:.4f}% (observation{obs_prev*100:.4f}%)")

    hi_thr_d = full["PRS311"].quantile(0.90)
    lo_thr_d = full["PRS311"].quantile(0.10)
    is_hi = (full["PRS311"].values >= hi_thr_d)
    is_lo = (full["PRS311"].values <= lo_thr_d)
    carr = full["carrier_main_10gene"].values.astype(bool)

    rng = np.random.default_rng(SIM_SEED)
    sim_diff, sim_rr, sim_hi_rate, sim_lo_rate, sim_ctrl_diff, sim_ncase = [], [], [], [], [], []
    for legacy_age_contrasts__ in range(N_SIM):
        y = rng.random(len(legacy_age_contrasts_p_case)) < legacy_age_contrasts_p_case
        ch = carr & y & is_hi
        cl = carr & y & is_lo
        nh, nl = (y & is_hi).sum(), (y & is_lo).sum()
        rh = ch.sum() / nh if nh else np.nan
        rl = cl.sum() / nl if nl else np.nan
        sim_hi_rate.append(rh); sim_lo_rate.append(rl)
        sim_diff.append(100 * (rh - rl))
        sim_rr.append(rh / rl if rl > 0 else np.nan)
        nch, ncl = (~y & is_hi).sum(), (~y & is_lo).sum()
        sim_ctrl_diff.append(100 * ((carr & ~y & is_hi).sum() / nch - (carr & ~y & is_lo).sum() / ncl))
        sim_ncase.append(int(y.sum()))

    sim_diff = np.array(sim_diff); sim_rr = np.array(sim_rr)
    obs_row_h = dec[(dec["group"] == "case") & (dec["prs_group"] == "high")].iloc[0]
    obs_row_l = dec[(dec["group"] == "case") & (dec["prs_group"] == "low")].iloc[0]
    obs_diff = obs_row_h["carrier_rate_pct"] - obs_row_l["carrier_rate_pct"]
    obs_rr = obs_row_h["carrier_rate_pct"] / obs_row_l["carrier_rate_pct"]

    sim_p_two = 2 * min((sim_diff <= obs_diff).mean(), (sim_diff >= obs_diff).mean())
    sim_p_two = min(1.0, max(sim_p_two, 1.0 / N_SIM))
    A2 = pd.DataFrame([{
        "n_sim": N_SIM, "seed": SIM_SEED, "sim_model": "logistic, no interaction",
        "beta_prs_used": b_prs, "beta_burden_used": b_bur, "alpha_calibrated": alpha,
        "target_prevalence_pct": 100 * obs_prev, "sim_mean_prevalence_pct": 100 * legacy_age_contrasts_p_case.mean(),
        "sim_mean_n_case": float(np.mean(sim_ncase)),
        "null_mean_carrier_rate_high_pct": 100 * np.nanmean(sim_hi_rate),
        "null_mean_carrier_rate_low_pct": 100 * np.nanmean(sim_lo_rate),
        "null_mean_diff_pct_points": sim_diff.mean(),
        "null_ci95_low_diff_pct_points": np.percentile(sim_diff, 2.5),
        "null_ci95_high_diff_pct_points": np.percentile(sim_diff, 97.5),
        "null_mean_rate_ratio": np.nanmean(sim_rr),
        "null_ci95_low_rate_ratio": np.nanpercentile(sim_rr, 2.5),
        "null_ci95_high_rate_ratio": np.nanpercentile(sim_rr, 97.5),
        "observed_diff_pct_points": obs_diff, "observed_rate_ratio": obs_rr,
        "empirical_two_sided_p_vs_null": sim_p_two,
        "observed_more_extreme_than_null": bool(obs_diff < np.percentile(sim_diff, 2.5)
                                                or obs_diff > np.percentile(sim_diff, 97.5)),
        "null_mean_control_diff_pct_points": float(np.mean(sim_ctrl_diff)),
    }])
    A2.to_csv(os.path.join(LEGACY_AGE_CONTRASTS_TABLE_DIR, "18_A2_additive_null_simulation.csv"), index=False)
    legacy_age_contrasts_log(f"\n case-high carrier rate under zero model ={100*np.nanmean(sim_hi_rate):.3f}%, "
        f"case-low = {100*np.nanmean(sim_lo_rate):.3f}%")
    legacy_age_contrasts_log(f"Null model diff(high-low) mean ={sim_diff.mean():+.3f} pp, "
        f"95% interval [{np.percentile(sim_diff,2.5):+.3f}, {np.percentile(sim_diff,97.5):+.3f}] pp")
    legacy_age_contrasts_log(f"Observation diff ={obs_diff:+.3f}pp(RRobservation={obs_rr:.3f}, zero model RR mean={np.nanmean(sim_rr):.3f})")
    legacy_age_contrasts_log(f"Whether the observation falls outside the 95% interval of the null model:"
        f"{'Yes -> beyond what can be explained by pure selection effects' if A2['observed_more_extreme_than_null'].iloc[0] else 'No -> indistinguishable from pure selection effect'}")
    legacy_age_contrasts_log(f"Empirical two-sided p (observed vs null distribution) ={sim_p_two:.4g}(Resolution lower limit 1/{N_SIM})")

    # ===========================================================================
    
    # ===========================================================================
    legacy_age_contrasts_section('Step B: Comparison between early onset vs late onset groups (cut point <45)')
    legacy_age_contrasts_log('⚠️ collider/index-event bias: The following are all case internal comparisons and are descriptive.')

    B_rows, B_tests = [], []
    group_defs = []
    for cut in CUTOFFS:
        lab_e, lab_l = f"early(<{cut})", f"late(>={cut})"
        legacy_age_contrasts_g = np.where(cases["age_at_diagnosis"] < cut, lab_e, lab_l)
        group_defs.append((f"cutoff_{cut}", legacy_age_contrasts_g, [lab_e, lab_l]))

    for scheme, garr, labs in group_defs:
        cases["_g"] = garr
        for lab in labs:
            legacy_age_contrasts_sub = cases[cases["_g"] == lab]
            n, k = len(legacy_age_contrasts_sub), int(legacy_age_contrasts_sub["carrier_main_10gene"].sum())
            lo_ci, hi_ci = legacy_age_contrasts_wilson_ci(k, n)
            B_rows.append({
                "scheme": scheme, "group": lab, "n": n, "n_disclosed": legacy_age_contrasts_suppress(n),
                "n_carrier": k, "n_carrier_disclosed": legacy_age_contrasts_suppress(k),
                "carrier_rate_pct": 100 * k / n,
                "carrier_rate_ci_low_pct": 100 * lo_ci, "carrier_rate_ci_high_pct": 100 * hi_ci,
                "burden_mean": legacy_age_contrasts_sub["burden_rare_main_ptv_10gene"].mean(),
                "burden_median": legacy_age_contrasts_sub["burden_rare_main_ptv_10gene"].median(),
                "burden_z_mean": legacy_age_contrasts_sub["burden_z"].mean(),
                "prs311_z_mean": legacy_age_contrasts_sub["PRS311_z"].mean(),
                "prs311_z_sd": legacy_age_contrasts_sub["PRS311_z"].std(),
                "prs311_z_fullref_mean": legacy_age_contrasts_sub["PRS311_z_fullref"].mean(),
                "burden_z_fullref_mean": legacy_age_contrasts_sub["burden_z_fullref"].mean(),
                "age_mean": legacy_age_contrasts_sub["age_at_diagnosis"].mean(),
            })
        if len(labs) == 2:
            e = cases[cases["_g"] == labs[0]]
            l = cases[cases["_g"] == labs[1]]
            ke, ne = int(e["carrier_main_10gene"].sum()), len(e)
            kl, nl = int(l["carrier_main_10gene"].sum()), len(l)
            chi2, p_chi, dof = chi2_2x2(ke, ne, kl, nl)
            u_b, p_b = stats.mannwhitneyu(e["burden_rare_main_ptv_10gene"],
                                          l["burden_rare_main_ptv_10gene"], alternative="two-sided")
            u_p, p_p = stats.mannwhitneyu(e["PRS311_z"], l["PRS311_z"], alternative="two-sided")
            t_p, pt_p = stats.ttest_ind(e["PRS311_z"], l["PRS311_z"], equal_var=False)
            t_b, pt_b = stats.ttest_ind(e["burden_z"], l["burden_z"], equal_var=False)
            smd_prs = smd(e["PRS311_z"], l["PRS311_z"])
            smd_bur = smd(e["burden_z"], l["burden_z"])
            B_tests.append({
                "scheme": scheme, "n_early": ne, "n_late": nl,
                "carrier_rate_early_pct": 100 * ke / ne, "carrier_rate_late_pct": 100 * kl / nl,
                "carrier_rate_diff_pp": 100 * (ke / ne - kl / nl),
                "carrier_rate_ratio": (ke / ne) / (kl / nl),
                "chi2_carrier": chi2, "p_chi2_carrier": p_chi,
                "burden_mean_early": e["burden_rare_main_ptv_10gene"].mean(),
                "burden_mean_late": l["burden_rare_main_ptv_10gene"].mean(),
                "mannwhitney_U_burden": u_b, "p_mannwhitney_burden": p_b,
                "prs311_z_mean_early": e["PRS311_z"].mean(), "prs311_z_mean_late": l["PRS311_z"].mean(),
                "prs311_z_fullref_mean_early": e["PRS311_z_fullref"].mean(),
                "prs311_z_fullref_mean_late": l["PRS311_z_fullref"].mean(),
                "prs311_z_diff": e["PRS311_z"].mean() - l["PRS311_z"].mean(),
                "mannwhitney_U_prs": u_p, "p_mannwhitney_prs": p_p,
                "t_stat_prs_welch": t_p, "p_ttest_prs_welch": pt_p,
                "t_stat_burden_welch": t_b, "p_ttest_burden_welch": pt_b,
                "smd_prs311_early_vs_late": smd_prs,
                "smd_burden_early_vs_late": smd_bur,
                "smd_ratio_burden_over_prs": abs(smd_bur) / abs(smd_prs) if smd_prs else np.nan,
            })
            legacy_age_contrasts_log(f"\n  [{scheme}] early n={ne} carrier={100*ke/ne:.2f}% PRS_z={e['PRS311_z'].mean():+.3f} | "
                f"late n={nl} carrier={100*kl/nl:.2f}% PRS_z={l['PRS311_z'].mean():+.3f}")
            legacy_age_contrasts_log(f"      carrier chi2={chi2:.2f} p={p_chi:.3g} (RR={((ke/ne)/(kl/nl)):.2f}); "
                f"burden MWU p={p_b:.3g}; PRS MWU p={p_p:.3g}; PRS Welch-t p={pt_p:.3g}")
            legacy_age_contrasts_log(f"Effect size comparison SMD(burden)={smd_bur:+.3f} vs SMD(PRS311)={smd_prs:+.3f} "
                f"-> |SMD_burden|/|SMD_PRS| = {abs(smd_bur)/abs(smd_prs):.1f}x")

    B_tbl = pd.DataFrame(B_rows)
    B_test_tbl = pd.DataFrame(B_tests)
    B_tbl.to_csv(os.path.join(LEGACY_AGE_CONTRASTS_TABLE_DIR, "18_B_early_late_group_summary.csv"), index=False)
    B_test_tbl.to_csv(os.path.join(LEGACY_AGE_CONTRASTS_TABLE_DIR, "18_B_early_late_tests.csv"), index=False)
    legacy_age_contrasts_log('Group summary table:')
    legacy_age_contrasts_log(B_tbl[["scheme", "group", "n", "n_carrier", "carrier_rate_pct", "burden_mean",
               "burden_median", "prs311_z_mean", "prs311_z_fullref_mean"]].to_string(index=False))

    for legacy_age_contrasts__, r in B_tbl.iterrows():
        for v in (r["n"], r["n_carrier"], r["n"] - r["n_carrier"]):
            assert v == 0 or v >= LEGACY_AGE_CONTRASTS_DISCLOSURE_MIN_CELL, f"{r['scheme']} {r['group']}There are cells 1-19"

    # ===========================================================================
    
    # ===========================================================================
    legacy_age_contrasts_section('Step C: Operationalization of the hybrid mechanism (PRS-high [full queue top decile] x carrier composed of four categories)')
    legacy_age_contrasts_log(f"The PRS-high threshold is taken from the **whole queue** PRS311 90th percentile ={hi_thr_d:.4f}(Not the quantile inside cases)")

    cases["prs_high"] = (cases["PRS311"] >= hi_thr_d).astype(int)

    cases["mix_cat"] = cases.apply(mix_cat, axis=1)
    MIX_ORDER = ["both (PRS-high AND carrier)", "PRS-high only", "carrier only", "neither"]

    C_rows, C_tests = [], []
    for cut in CUTOFFS:
        lab_e, lab_l = f"early(<{cut})", f"late(>={cut})"
        cases["_g"] = np.where(cases["age_at_diagnosis"] < cut, lab_e, lab_l)
        ct = pd.crosstab(cases["_g"], cases["mix_cat"]).reindex(columns=MIX_ORDER, fill_value=0)
        chi2, p, dof, legacy_age_contrasts__ = stats.chi2_contingency(ct.values)

        _small = int(ct.values.min()) < LEGACY_AGE_CONTRASTS_DISCLOSURE_MIN_CELL
        if _small:
            chi2 = p = np.nan
        C_tests.append({"scheme": f"cutoff_{cut}", "table": "2x4 (early/late x 4 mixture cats)",
                        "chi2": chi2, "df": dof, "p_value": p,
                        "n_total": int(ct.values.sum()),
                        "min_cell_n_bucket": bucket(int(ct.values.min())),
                        "note": "not reported (input cell <20)" if _small else ""})
        legacy_age_contrasts_log(f"\n  [cutoff_{cut}] 2x4 " + ('Inspection is not reported (input contains <20 spaces)' if _small else f"chi2={chi2:.2f}, df={dof}, p={p:.3g}")
            + f"; Minimum grid n={bucket(int(ct.values.min()))}")
        for lab in [lab_e, lab_l]:
            n_tot = int(ct.loc[lab].sum())
            for cat in MIX_ORDER:
                k = int(ct.loc[lab, cat])
                lo_ci, hi_ci = legacy_age_contrasts_wilson_ci(k, n_tot)
                C_rows.append({"scheme": f"cutoff_{cut}", "group": lab, "mixture_category": cat,
                               "n": k, "n_disclosed": legacy_age_contrasts_suppress(k), "n_group_total": n_tot,
                               "pct_of_group": 100 * k / n_tot,
                               "pct_ci_low": 100 * lo_ci, "pct_ci_high": 100 * hi_ci})
        
        kb_e = int(ct.loc[lab_e, MIX_ORDER[0]]); nb_e = int(ct.loc[lab_e].sum())
        kb_l = int(ct.loc[lab_l, MIX_ORDER[0]]); nb_l = int(ct.loc[lab_l].sum())
        c2, p2, d2 = chi2_2x2(kb_e, nb_e, kb_l, nb_l)
        _small2 = min(kb_e, kb_l) < LEGACY_AGE_CONTRASTS_DISCLOSURE_MIN_CELL
        if _small2:
            c2 = p2 = np.nan
        C_tests.append({"scheme": f"cutoff_{cut}", "table": "2x2 both-vs-rest (early vs late)",
                        "chi2": c2, "df": d2, "p_value": p2, "n_total": nb_e + nb_l,
                        "min_cell_n_bucket": bucket(min(kb_e, kb_l)),
                        "note": "not reported (input cell <20)" if _small2 else ""})
        if _small2:
            legacy_age_contrasts_log(f"\'both\' category early vs late: neither proportion nor test is reported (with grid <20) (minimum grid n={bucket(min(kb_e,kb_l))})")
        else:
            legacy_age_contrasts_log(f"\'both\' class early vs late: early{100*kb_e/nb_e:.2f}% vs late {100*kb_l/nb_l:.2f}%, "
                f"chi2={c2:.2f}, p={p2:.3g}(minimum grid n={bucket(min(kb_e,kb_l))})")

    C_tbl = pd.DataFrame(C_rows)
    C_test_tbl = pd.DataFrame(C_tests)

    _PARTNER = {MIX_ORDER[0]: MIX_ORDER[2], MIX_ORDER[2]: MIX_ORDER[0],
                MIX_ORDER[1]: MIX_ORDER[3], MIX_ORDER[3]: MIX_ORDER[1]}
    C_tbl["_secondary"] = False
    for (sch, legacy_age_contrasts_grp), legacy_age_contrasts_g in C_tbl.groupby(["scheme", "group"]):
        for idx, r in legacy_age_contrasts_g.iterrows():
            if r["n"] < LEGACY_AGE_CONTRASTS_DISCLOSURE_MIN_CELL:
                partner = legacy_age_contrasts_g.index[legacy_age_contrasts_g["mixture_category"] == _PARTNER[r["mixture_category"]]][0]
                if C_tbl.at[partner, "n"] >= LEGACY_AGE_CONTRASTS_DISCLOSURE_MIN_CELL:
                    C_tbl.at[partner, "_secondary"] = True
    C_export = sanitize_for_export(C_tbl, "n", ["pct_of_group", "pct_ci_low", "pct_ci_high"])
    C_export.loc[C_tbl["_secondary"], "n_disclosed"] = "suppressed (secondary)"
    C_export.loc[C_tbl["_secondary"], ["pct_of_group", "pct_ci_low", "pct_ci_high"]] = np.nan
    C_export = C_export.drop(columns=["_secondary"])
    C_export.to_csv(os.path.join(LEGACY_AGE_CONTRASTS_TABLE_DIR, "18_C_mixture_category_composition.csv"), index=False)
    legacy_age_contrasts_log('Four categories (after disclosure processing):')
    legacy_age_contrasts_log(C_export[["scheme", "group", "mixture_category", "n_disclosed", "n_group_total"]].to_string(index=False))
    C_test_tbl.to_csv(os.path.join(LEGACY_AGE_CONTRASTS_TABLE_DIR, "18_C_mixture_category_tests.csv"), index=False)
    n_c_small = int((C_tbl["n"] < LEGACY_AGE_CONTRASTS_DISCLOSURE_MIN_CELL).sum())
    legacy_age_contrasts_log(f"\n Four categories of composition table{len(C_tbl)}lattice, where n<20{n_c_small}Format: The external table only retains the n_disclosed position, and the original n and pct/CI are not included in the disk.")

    # ===========================================================================
    
    # ===========================================================================
    legacy_age_contrasts_section('Step D: Interaction test (the part not done in script 17) + power calculation')
    legacy_age_contrasts_log("⚠️ A prior statement has been given in Step 0: Any non-significant results in this section can only be expressed as 'insufficient power to rule out interaction'.")
    legacy_age_contrasts_log('⚠️ collider bias: The interaction term of cases-only regression is also descriptive.')

    D_rows = []
    # --- D1: age ~ PRS_z * carrier + PC1-10 (OLS HC1) ---
    cases["PRSxCarrier"] = cases["PRS311_z"] * cases["carrier_main_10gene"]
    X = build_X(cases, ["PRS311_z", "carrier_main_10gene", "PRSxCarrier"] + PC_COLS)
    fit_d1 = fit_ols_hc1(cases["age_at_diagnosis"], X)
    d1 = coef_rows(fit_d1, "D1_age ~ PRS311_z * carrier + PC1-10 (OLS HC1)", len(cases))
    d1["r_squared"] = fit_d1.rsquared
    D_rows.append(d1)
    legacy_age_contrasts_log("\n  [D1] age_at_diagnosis ~ PRS311_z * carrier_main_10gene + PC1-10 (n=%d)" % len(cases))
    legacy_age_contrasts_log(d1[d1["term"].isin(["PRS311_z", "carrier_main_10gene", "PRSxCarrier"])]
        [["term", "beta", "robust_se_HC1", "ci_low", "ci_high", "p_value"]].to_string(index=False))
    se_inter_d1 = float(d1.loc[d1["term"] == "PRSxCarrier", "robust_se_HC1"].iloc[0])
    b_inter_d1 = float(d1.loc[d1["term"] == "PRSxCarrier", "beta"].iloc[0])
    p_inter_d1 = float(d1.loc[d1["term"] == "PRSxCarrier", "p_value"].iloc[0])

    strat_rows = []
    for lab, legacy_age_contrasts_sub in [("carrier (burden>0)", cases[cases["carrier_main_10gene"] == 1]),
                     ("non-carrier", cases[cases["carrier_main_10gene"] == 0])]:
        Xs = build_X(legacy_age_contrasts_sub, ["PRS311_z"] + PC_COLS)
        fs = fit_ols_hc1(legacy_age_contrasts_sub["age_at_diagnosis"], Xs)
        r = coef_rows(fs, f"D1_stratified: {lab}", len(legacy_age_contrasts_sub), terms=["PRS311_z"])
        strat_rows.append(r)
        legacy_age_contrasts_log(f"Hierarchical PRS slope [{lab}] n={len(legacy_age_contrasts_sub)}: beta={fs.params['PRS311_z']:+.4f} "
            f"Year/SD, SE={fs.bse['PRS311_z']:.4f}, p={fs.pvalues['PRS311_z']:.3g}")
    D_rows.extend(strat_rows)

    # --- D2: age ~ PRS_z * burden_z + PC1-10 ---
    legacy_age_contrasts_log(f"\n [D2] age ~ PRS311_z * burden_z + PC1-10 —— ⚠️Note: The proportion of burden==0 inside cases is"
        f"{pct_zero_burden:.1f}%, the continuous x-continuous interaction term actually consists of only a few (n={int(cases['carrier_main_10gene'].sum())})"
        f"Non-zero individual drive, limited explanatory power.")
    cases["PRSxBurden"] = cases["PRS311_z"] * cases["burden_z"]
    X2 = build_X(cases, ["PRS311_z", "burden_z", "PRSxBurden"] + PC_COLS)
    fit_d2 = fit_ols_hc1(cases["age_at_diagnosis"], X2)
    d2 = coef_rows(fit_d2, "D2_age ~ PRS311_z * burden_z + PC1-10 (OLS HC1)", len(cases))
    d2["r_squared"] = fit_d2.rsquared
    D_rows.append(d2)
    legacy_age_contrasts_log(d2[d2["term"].isin(["PRS311_z", "burden_z", "PRSxBurden"])]
        [["term", "beta", "robust_se_HC1", "ci_low", "ci_high", "p_value"]].to_string(index=False))

    # --- D3: early(0/1) ~ PRS_z * burden_z + PC1-10 (logistic) ---
    for cut in CUTOFFS:
        y_bin = (cases["age_at_diagnosis"] < cut).astype(int)
        X3 = build_X(cases, ["PRS311_z", "burden_z", "PRSxBurden"] + PC_COLS)
        f3 = fit_logit_hc1(y_bin, X3)
        r3 = coef_rows(f3, f"D3_logistic early(<{cut}) ~ PRS311_z * burden_z + PC1-10",
                       len(cases), exponentiate=True)
        D_rows.append(r3)
        sel = r3[r3["term"].isin(["PRS311_z", "burden_z", "PRSxBurden"])]
        legacy_age_contrasts_log(f"\n  [D3 cut<{cut}] n_event={int(y_bin.sum())}/{len(cases)}")
        legacy_age_contrasts_log(sel[["term", "beta", "robust_se_HC1", "OR", "OR_ci_low", "OR_ci_high",
                 "p_value"]].to_string(index=False))

    D_tbl = pd.concat(D_rows, ignore_index=True)
    D_tbl.to_csv(os.path.join(LEGACY_AGE_CONTRASTS_TABLE_DIR, "18_D_interaction_models.csv"), index=False)

    legacy_age_contrasts_section('Step D4: Interaction term power - minimum detectable effect (MDE) at 80% power')
    z_a, z_b = stats.norm.ppf(0.975), stats.norm.ppf(0.80)
    pow_rows = []
    for lab, fitres, term in [('D1 PRS311_z x carrier (year/SD)', fit_d1, "PRSxCarrier"),
                              ('D2 PRS311_z x burden_z (year/SD^2)', fit_d2, "PRSxBurden")]:
        se = float(fitres.bse[term]); legacy_age_contrasts_b = float(fitres.params[term])
        mde = (z_a + z_b) * se
        pow_rows.append({"model": lab, "observed_beta": legacy_age_contrasts_b, "robust_se_HC1": se,
                         "p_value": float(fitres.pvalues[term]),
                         "mde_80pct_power_alpha0.05": mde,
                         "observed_abs_beta_over_mde": abs(legacy_age_contrasts_b) / mde,
                         "n": len(cases), "n_carrier": int(cases["carrier_main_10gene"].sum()),
                         "carrier_pct": 100 * cases["carrier_main_10gene"].mean()})
        legacy_age_contrasts_log(f"  [{lab}] beta={legacy_age_contrasts_b:+.4f}, SE={se:.4f}, p={float(fitres.pvalues[term]):.3g}  "
            f"-> MDE (80% efficacy) ={mde:.4f};observation|beta|/MDE ={abs(legacy_age_contrasts_b)/mde:.2f}")
    pow_tbl = pd.DataFrame(pow_rows)
    pow_tbl.to_csv(os.path.join(LEGACY_AGE_CONTRASTS_TABLE_DIR, "18_D_interaction_power.csv"), index=False)
    legacy_age_contrasts_log('Interpretation: Only when |true effect| >= MDE, this sample size can be detected with 80% confidence. Interactions below the MDE cannot be excluded from this study.')

    # ===========================================================================
    
    # ===========================================================================
    legacy_age_contrasts_section('Step E: Total-effect model adding environmental/clinical covariates (BMI/menopause/family history/smoking)')
    legacy_age_contrasts_log("⚠️ This is the only alternative pathway available to this project for the 'environmental exposure' half of the Aim 4 narrative - these are the baseline clinical/questionnaire variables,"
        'Not a true measure of environmental exposure, labeled total-effect covariate adjustment rather than causal mediation analysis.')

    env_cols = ["bmi_value", "menopause_status", "family_history_breast_cancer", "smoking_status"]
    cases_e = cases.merge(clinical[["person_id"] + env_cols], on="person_id", how="left")
    miss_rows = []
    for legacy_age_contrasts_c in env_cols:
        n_na = int(cases_e[legacy_age_contrasts_c].isna().sum())
        miss_rows.append({"covariate": legacy_age_contrasts_c, "n_missing": n_na, "n_total": len(cases_e),
                          "pct_missing": 100 * n_na / len(cases_e)})
        legacy_age_contrasts_log(f"    {legacy_age_contrasts_c}: missing{n_na}/{len(cases_e)} ({100*n_na/len(cases_e):.2f}%)")
    miss_tbl = pd.DataFrame(miss_rows)
    miss_tbl.to_csv(os.path.join(LEGACY_AGE_CONTRASTS_TABLE_DIR, "18_E_env_covariate_missingness.csv"), index=False)

    E_rows = []
    specs = [
        ('E1_all 4 covariates (including family history, complete cases)', env_cols),
        ('E2_No family history (BMI/menopause/smoking, complete case)', ["bmi_value", "menopause_status", "smoking_status"]),
    ]
    for lab, cc in specs:
        legacy_age_contrasts_sub = cases_e.dropna(subset=cc).reset_index(drop=True)
        legacy_age_contrasts_log(f"\n  [{lab}] Complete cases n={len(legacy_age_contrasts_sub)}(accounting for cases queue{100*len(legacy_age_contrasts_sub)/len(cases_e):.1f}%)")
        
        Xe = make_env_X(legacy_age_contrasts_sub, cc, ["burden_z", "PRS311_z"])
        fe = fit_ols_hc1(legacy_age_contrasts_sub["age_at_diagnosis"], Xe)
        r = coef_rows(fe, f"{lab} | age ~ burden_z + PRS311_z + env + PC1-10", len(legacy_age_contrasts_sub))
        r["r_squared"] = fe.rsquared
        E_rows.append(r)
        legacy_age_contrasts_log(r[r["term"].isin(["burden_z", "PRS311_z"])]
            [["term", "beta", "robust_se_HC1", "ci_low", "ci_high", "p_value"]].to_string(index=False))
        
        Xn = build_X(legacy_age_contrasts_sub, ["burden_z", "PRS311_z"] + PC_COLS)
        fn = fit_ols_hc1(legacy_age_contrasts_sub["age_at_diagnosis"], Xn)
        rn = coef_rows(fn, f"{lab}| Control with the same subsample (excluding env)", len(legacy_age_contrasts_sub))
        rn["r_squared"] = fn.rsquared
        E_rows.append(rn)
        legacy_age_contrasts_log('Same sub-sample without env control:' + ", ".join(
            f"{t}={fn.params[t]:+.4f}(p={fn.pvalues[t]:.3g})" for t in ["burden_z", "PRS311_z"]))
        
        for cut in CUTOFFS:
            y_bin = (legacy_age_contrasts_sub["age_at_diagnosis"] < cut).astype(int)
            fl = fit_logit_hc1(y_bin, Xe)
            rl = coef_rows(fl, f"{lab} | logistic early(<{cut}) ~ burden_z + PRS311_z + env + PC1-10",
                           len(legacy_age_contrasts_sub), exponentiate=True)
            E_rows.append(rl)
            sel = rl[rl["term"].isin(["burden_z", "PRS311_z"])]
            legacy_age_contrasts_log(f"    early(<{cut}) logistic (n_event={int(y_bin.sum())}): " + ", ".join(
                f"{row.term} OR={row.OR:.3f} [{row.OR_ci_low:.3f},{row.OR_ci_high:.3f}] p={row.p_value:.3g}"
                for row in sel.itertuples()))

    E_tbl = pd.concat(E_rows, ignore_index=True)
    E_tbl.to_csv(os.path.join(LEGACY_AGE_CONTRASTS_TABLE_DIR, "18_E_env_covariate_models.csv"), index=False)

    # ===========================================================================
    
    # ===========================================================================
    legacy_age_contrasts_section('Step F: Plot (matplotlib default style, no grid, SVG priority)')
    fig_files = []

    legacy_age_contrasts_fig, legacy_age_contrasts_ax = plt.subplots(figsize=(6, 4.2))
    legacy_age_contrasts_d = dec.set_index(["group", "prs_group"])
    labels = ["case\nPRS-high", "case\nPRS-low", "control\nPRS-high", "control\nPRS-low"]
    keys = [("case", "high"), ("case", "low"), ("control", "high"), ("control", "low")]
    vals = [legacy_age_contrasts_d.loc[k, "carrier_rate_pct"] for k in keys]
    errs = [[legacy_age_contrasts_d.loc[k, "carrier_rate_pct"] - legacy_age_contrasts_d.loc[k, "carrier_rate_ci_low_pct"] for k in keys],
            [legacy_age_contrasts_d.loc[k, "carrier_rate_ci_high_pct"] - legacy_age_contrasts_d.loc[k, "carrier_rate_pct"] for k in keys]]
    colors = ["tab:blue", "tab:blue", "tab:orange", "tab:orange"]
    legacy_age_contrasts_ax.bar(labels, vals, yerr=errs, capsize=4, color=colors)
    for legacy_age_contrasts_i, k in enumerate(keys):
        legacy_age_contrasts_ax.text(legacy_age_contrasts_i, vals[legacy_age_contrasts_i] + errs[1][legacy_age_contrasts_i] + 0.3, f"n={int(legacy_age_contrasts_d.loc[k,'n'])}", ha="center", fontsize=8)
    legacy_age_contrasts_ax.set_ylabel("PTV carrier rate (%)")
    legacy_age_contrasts_ax.set_title(f"PTV carrier rate by PRS311 decile group and case status\n(full cohort, n={len(full):,}; error bars = 95% Wilson CI)",
                 fontsize=10)
    legacy_age_contrasts_ax.grid(False)
    legacy_age_contrasts_fig.tight_layout()
    p = os.path.join(LEGACY_AGE_CONTRASTS_FIG_DIR, "18_A_carrier_rate_2x2_bar.svg"); legacy_age_contrasts_fig.savefig(p); plt.close(legacy_age_contrasts_fig)
    fig_files.append(p)

    legacy_age_contrasts_fig, legacy_age_contrasts_ax = plt.subplots(figsize=(6, 4.2))
    legacy_age_contrasts_ax.hist(sim_diff, bins=40, color="tab:gray")
    legacy_age_contrasts_ax.axvline(obs_diff, color="tab:red", linewidth=2)
    legacy_age_contrasts_ax.text(obs_diff, legacy_age_contrasts_ax.get_ylim()[1] * 0.9, f" observed = {obs_diff:+.2f} pp",
            color="tab:red", fontsize=9)
    legacy_age_contrasts_ax.set_xlabel("carrier rate difference, case PRS-high minus case PRS-low (percentage points)")
    legacy_age_contrasts_ax.set_ylabel("simulation count")
    legacy_age_contrasts_ax.set_title(f"Additive (no-interaction) null benchmark, {N_SIM} simulations", fontsize=10)
    legacy_age_contrasts_ax.grid(False)
    legacy_age_contrasts_fig.tight_layout()
    p = os.path.join(LEGACY_AGE_CONTRASTS_FIG_DIR, "18_A2_additive_null_simulation_hist.svg"); legacy_age_contrasts_fig.savefig(p); plt.close(legacy_age_contrasts_fig)
    fig_files.append(p)

    edges = [0, 35, 40, 45, 50, 55, 60, 65, 70, 200]
    raw = pd.cut(cases["age_at_diagnosis"], bins=edges, right=False, labels=False)
    merged_edges, lo_i = [edges[0]], 0
    for hi_i in range(1, len(edges)):
        m = (raw >= lo_i) & (raw < hi_i)
        k = int(cases.loc[m, "carrier_main_10gene"].sum()); n = int(m.sum())
        if k >= LEGACY_AGE_CONTRASTS_DISCLOSURE_MIN_CELL and n - k >= LEGACY_AGE_CONTRASTS_DISCLOSURE_MIN_CELL:
            merged_edges.append(edges[hi_i]); lo_i = hi_i
    if merged_edges[-1] != edges[-1]:
        merged_edges[-1] = edges[-1]  

    blab = [_lab(legacy_age_contrasts_a, legacy_age_contrasts_b) for legacy_age_contrasts_a, legacy_age_contrasts_b in zip(merged_edges[:-1], merged_edges[1:])]
    legacy_age_contrasts_log(f"Age binning (after merge):{blab}")
    cases["_agebin"] = pd.cut(cases["age_at_diagnosis"], bins=merged_edges, labels=blab, right=False)
    agg = cases.groupby("_agebin", observed=False)["carrier_main_10gene"].agg(["size", "sum"])
    assert ((agg["sum"] >= LEGACY_AGE_CONTRASTS_DISCLOSURE_MIN_CELL) & (agg["size"] - agg["sum"] >= LEGACY_AGE_CONTRASTS_DISCLOSURE_MIN_CELL)).all()
    rates = 100 * agg["sum"] / agg["size"]
    cis = [legacy_age_contrasts_wilson_ci(k, n) for k, n in zip(agg["sum"], agg["size"])]
    legacy_age_contrasts_fig, legacy_age_contrasts_ax = plt.subplots(figsize=(6.5, 4.2))
    legacy_age_contrasts_ax.errorbar(range(len(agg)), rates.values,
                yerr=[rates.values - np.array([100 * legacy_age_contrasts_c[0] for legacy_age_contrasts_c in cis]),
                      np.array([100 * legacy_age_contrasts_c[1] for legacy_age_contrasts_c in cis]) - rates.values],
                fmt="o-", capsize=4, color="tab:blue")
    legacy_age_contrasts_ax.set_xticks(range(len(agg))); legacy_age_contrasts_ax.set_xticklabels(agg.index.astype(str))
    legacy_age_contrasts_ax.set_xlabel("Age at diagnosis (years)")
    legacy_age_contrasts_ax.set_ylabel("PTV carrier rate (%)")
    legacy_age_contrasts_fig.tight_layout()
    p = os.path.join(LEGACY_AGE_CONTRASTS_FIG_DIR, "18_B_carrier_rate_by_age_bin.svg"); legacy_age_contrasts_fig.savefig(p); plt.close(legacy_age_contrasts_fig)
    fig_files.append(p)
    agg_out = pd.DataFrame({"n": agg["size"], "n_carrier": agg["sum"], "carrier_rate_pct": rates,
                            "ci_low_pct": [100 * legacy_age_contrasts_c[0] for legacy_age_contrasts_c in cis], "ci_high_pct": [100 * legacy_age_contrasts_c[1] for legacy_age_contrasts_c in cis]})
    agg_out.reset_index().rename(columns={"_agebin": "age_bin"}).to_csv(
        os.path.join(LEGACY_AGE_CONTRASTS_TABLE_DIR, "18_B_carrier_rate_by_age_bin.csv"), index=False)

    agg_p = cases.groupby("_agebin", observed=False)["PRS311_z"].agg(["size", "mean", "sem"])
    legacy_age_contrasts_fig, legacy_age_contrasts_ax = plt.subplots(figsize=(6.5, 4.2))
    legacy_age_contrasts_ax.errorbar(range(len(agg_p)), agg_p["mean"].values, yerr=Z95 * agg_p["sem"].values,
                fmt="s-", capsize=4, color="tab:orange")
    legacy_age_contrasts_ax.axhline(0, color="black", linewidth=0.8)
    legacy_age_contrasts_ax.set_xticks(range(len(agg_p))); legacy_age_contrasts_ax.set_xticklabels(agg_p.index.astype(str))
    legacy_age_contrasts_ax.set_xlabel("Age at diagnosis (years)")
    legacy_age_contrasts_ax.set_ylabel("Mean PRS311 (z, within cases)")
    legacy_age_contrasts_fig.tight_layout()
    p = os.path.join(LEGACY_AGE_CONTRASTS_FIG_DIR, "18_B_prs311_by_age_bin.svg"); legacy_age_contrasts_fig.savefig(p); plt.close(legacy_age_contrasts_fig)
    fig_files.append(p)

    legacy_age_contrasts_fig, legacy_age_contrasts_ax = plt.subplots(figsize=(7, 4.2))
    plot_groups, bottoms = [], None
    sub_c = C_tbl[C_tbl["scheme"].isin([f"cutoff_{legacy_age_contrasts_c}" for legacy_age_contrasts_c in CUTOFFS])]

    xlabels, mat, merged = [], [], []
    for cut in CUTOFFS:
        for which in ["early", "late"]:
            lab = f"early(<{cut})" if which == "early" else f"late(>={cut})"
            legacy_age_contrasts_row = sub_c[(sub_c["scheme"] == f"cutoff_{cut}") & (sub_c["group"] == lab)].set_index("mixture_category")
            is_merged = bool((legacy_age_contrasts_row["n"] < LEGACY_AGE_CONTRASTS_DISCLOSURE_MIN_CELL).any())
            mat.append([float(legacy_age_contrasts_row.loc[cat, "pct_of_group"]) for cat in MIX_ORDER])
            merged.append(is_merged)
            xlabels.append(f"<{cut}\n{which}")
    mat = np.array(mat)
    bottom = np.zeros(len(mat))
    cols = ["tab:red", "tab:orange", "tab:blue", "lightgrey"]
    for j, cat in enumerate(MIX_ORDER):
        for legacy_age_contrasts_i in range(len(mat)):
            if merged[legacy_age_contrasts_i] and j == 0:
                continue  
            if merged[legacy_age_contrasts_i] and j == 2:
                h = mat[legacy_age_contrasts_i, 0] + mat[legacy_age_contrasts_i, 2]
                legacy_age_contrasts_ax.bar(legacy_age_contrasts_i, h, bottom=bottom[legacy_age_contrasts_i], color="white", edgecolor="tab:blue", hatch="//")
                bottom[legacy_age_contrasts_i] += h
                continue
            legacy_age_contrasts_ax.bar(legacy_age_contrasts_i, mat[legacy_age_contrasts_i, j], bottom=bottom[legacy_age_contrasts_i], color=cols[j])
            bottom[legacy_age_contrasts_i] += mat[legacy_age_contrasts_i, j]

    handles = [Patch(color=legacy_age_contrasts_c, label=l) for legacy_age_contrasts_c, l in zip(cols, MIX_ORDER)]
    if any(merged):
        handles.append(Patch(facecolor="white", edgecolor="tab:blue", hatch="//",
                             label="both + carrier only (not split: a cell <20)"))
    legacy_age_contrasts_ax.set_xticks(range(len(mat))); legacy_age_contrasts_ax.set_xticklabels(xlabels, fontsize=8)
    legacy_age_contrasts_ax.set_ylabel("% of group")
    legacy_age_contrasts_ax.set_title("Mixture-mechanism category composition, early vs late onset (cases)", fontsize=10)
    legacy_age_contrasts_ax.legend(handles=handles, fontsize=7, loc="lower right")
    legacy_age_contrasts_ax.grid(False)
    legacy_age_contrasts_fig.tight_layout()
    p = os.path.join(LEGACY_AGE_CONTRASTS_FIG_DIR, "18_C_mixture_composition_stacked_bar.svg"); legacy_age_contrasts_fig.savefig(p); plt.close(legacy_age_contrasts_fig)
    fig_files.append(p)

    for f in fig_files:
        legacy_age_contrasts_log(f"  saved {os.path.basename(f)}  ({os.path.getsize(f)/1024:.1f} KB)")
    legacy_age_contrasts_log(f"Total image size:{sum(os.path.getsize(f) for f in fig_files)/1024:.1f} KB")

    # ===========================================================================
    
    # ===========================================================================
    legacy_age_contrasts_section('Step G: Generate report results/report/18_aim4_group_contrasts_summary.md')

    t45 = B_test_tbl[B_test_tbl["scheme"] == "cutoff_45"].iloc[0]
    a2 = A2.iloc[0]

    flow_md = "\n".join(f"| {r['step']} | {r['n']:,} |" for _, r in flow_tbl.iterrows())

    a_md = "\n".join(
        f"| {r['group']} | {r['prs_group']} | {r['n']:,} | {r['n_carrier']:,} | "
        f"{r['carrier_rate_pct']:.2f} | [{r['carrier_rate_ci_low_pct']:.2f}, {r['carrier_rate_ci_high_pct']:.2f}] | "
        f"{r['burden_mean']:.5f} | {r['burden_median']:.3f} |"
        for _, r in dec.iterrows())
    quin = A_tbl[A_tbl["scheme"].str.startswith("quintile")]
    q_md = "\n".join(
        f"| {r['group']} | {r['prs_group']} | {r['n']:,} | {r['n_carrier']:,} | {r['carrier_rate_pct']:.2f} |"
        for _, r in quin.iterrows())
    atest_md = "\n".join(
        f"| {r['scheme']} | {r['contrast']} | {r['carrier_rate_high_pct']:.2f} | {r['carrier_rate_low_pct']:.2f} | "
        f"{r['diff_pct_points']:+.2f} | {r['rate_ratio_high_over_low']:.3f} | {r['chi2']:.2f} | {pfmt(r['p_value'])} |"
        for _, r in A_tests.iterrows())

    b_md = "\n".join(
        f"| {r['scheme']} | {r['group']} | {r['n']:,} | {r['n_carrier']:,} | {r['carrier_rate_pct']:.2f} | "
        f"{r['burden_mean']:.5f} | {r['burden_median']:.3f} | {r['prs311_z_mean']:+.3f} | "
        f"{r['prs311_z_fullref_mean']:+.3f} |"
        for _, r in B_tbl.iterrows())
    btest_md = "\n".join(
        f"| {r['scheme']} | {r['n_early']:,} / {r['n_late']:,} | {r['carrier_rate_early_pct']:.2f} / {r['carrier_rate_late_pct']:.2f} | "
        f"{r['carrier_rate_ratio']:.2f} | {pfmt(r['p_chi2_carrier'])} | {r['prs311_z_mean_early']:+.3f} / {r['prs311_z_mean_late']:+.3f} | "
        f"{pfmt(r['p_ttest_prs_welch'])} | {r['smd_burden_early_vs_late']:+.3f} | {r['smd_prs311_early_vs_late']:+.3f} | "
        f"{r['smd_ratio_burden_over_prs']:.1f}x |"
        for _, r in B_test_tbl.iterrows())

    c_md = "\n".join(
        f"| {r['scheme']} | {r['group']} | {r['mixture_category']} | {r['n_disclosed']} | {_pct_cell(r)} |"
        for _, r in C_export.iterrows())
    ctest_md = "\n".join(
        (f"| {r['scheme']} | {r['table']} | not reported | {r['df']} | not reported | {r['min_cell_n_bucket']} |"
         if pd.isna(r['chi2']) else
         f"| {r['scheme']} | {r['table']} | {r['chi2']:.2f} | {r['df']} | {pfmt(r['p_value'])} | {r['min_cell_n_bucket']} |")
        for _, r in C_test_tbl.iterrows())

    d_key = D_tbl[D_tbl["term"].isin(["PRS311_z", "carrier_main_10gene", "PRSxCarrier",
                                      "burden_z", "PRSxBurden"])]
    d_md = "\n".join(
        f"| {r['model']} | {r['term']} | {r['beta']:+.4f} | {r['robust_se_HC1']:.4f} | "
        f"[{r['ci_low']:+.4f}, {r['ci_high']:+.4f}] | {pfmt(r['p_value'])} |"
        for _, r in d_key.iterrows())
    pow_md = "\n".join(
        f"| {r['model']} | {r['observed_beta']:+.4f} | {r['robust_se_HC1']:.4f} | {pfmt(r['p_value'])} | "
        f"{r['mde_80pct_power_alpha0.05']:.4f} | {r['observed_abs_beta_over_mde']:.2f} |"
        for _, r in pow_tbl.iterrows())
    miss_md = "\n".join(
        f"| {r['covariate']} | {r['n_missing']:,} | {r['pct_missing']:.2f} |" for _, r in miss_tbl.iterrows())
    e_key = E_tbl[E_tbl["term"].isin(["burden_z", "PRS311_z"])]
    e_md = "\n".join(
        f"| {r['model']} | {r['term']} | {r['beta']:+.4f} | {r['robust_se_HC1']:.4f} | "
        f"[{r['ci_low']:+.4f}, {r['ci_high']:+.4f}] | {pfmt(r['p_value'])} | {r['n']:,} |"
        for _, r in e_key.iterrows())

    verdict_null = ('**Exceeds** the 95% interval of the null model' if bool(a2["observed_more_extreme_than_null"])
                    else '**Falls within the 95% interval of** the null model')
    verdict_sentence = (
        'Comparison of observed PRS-high vs PRS-low carrier rates within cases' +
        (''
         if bool(a2["observed_more_extreme_than_null"]) else
         'Not beyond what would be induced by a purely additive model + case selection itself, so **cannot** be interpreted as evidence of interaction/mechanism')
    )

    legacy_age_contrasts_report = f"""# Aim 4 Supplementary analysis｜Comparison between groups, operationalization of mixed mechanism and interaction test summary report

    Corresponding script: `scripts/10_exploratory_age_ehr.py`. Running log:
    `scripts/script_output_py/output_18_aim4_group_contrasts.txt`.

    This report supplements `results/report/17_step5_aim4_age_of_onset_model_summary.md` and does not replace it.
    `scripts/90_legacy_ptv_baseline.py` (pre-registered additive confirmatory model
    `age_at_diagnosis ~ Burden_rare + PRS311 + PC1-10`, no interaction terms) was left as is without any modification.
    This script fills in the parts not covered by script 17 in the original Aim 4 specification: comparison between groups (A/B), the literal meaning of \"mixing mechanism\"
    Operationalization (C), interaction test (D), environmental covariate total-effect model (E), and an interpretation of A
    Additive null model simulation benchmark (A2).

    ---

    ## 0. Three qualifications that must be read first

    ### (a) somatic mutation burden not analyzed - Aim 4 spec only partially implemented

    The original Aim 4 specification required a comparison of three exposures: rare variant burden, common PRS, and **somatic mutation burden**.
    This environment only has **germline whole blood WGS** from All of Us, but no tumor tissue sequencing data;
    `results/report/02_pusztai_gerstein_wes_data_check.md` has confirmed that the external WES queue is not available in this environment.
    Therefore, **the somatic mutation load branch has not been done at all, and no surrogate indicators have been used to approximate it** (no unfounded proxies are made).
    The specification for Aim 4 is therefore only **partially implemented**, and any statement involving \"comparison of three mechanisms\" must be downgraded to \"comparison of two germline mechanisms\".

    ### (b) collider / index-event bias (throughout all cases-only comparisons of A/B/C/D/E)

    Once conditioned on \"disease\", multiple independent risk factors will be automatically induced to have negative correlations, even if they are mutually exclusive in the population.
    Independent, biologically uninteractive. Therefore, all case-only between-group comparisons in this report (B/C/D/E) and A
    Within-case comparisons are **descriptive** and cannot be read directly as evidence of mechanism. Aim 1 (`scripts/16`, case-control design)
    Not affected by this bias - **the two cannot be interpreted at the same level**. Section A2 gives an idea of \"how much contrast a pure selection effect can cause\"
    The purpose of quantitative benchmarks is to make this warning measurable rather than just empty talk.

    ### (c) Advance statement of efficacy (written before looking at the results)

    Aim 4 cases queue n={len(cases):,}, of which the PTV carrier of the 10-gene panel only{100*cases['carrier_main_10gene'].mean():.2f}%（n={int(cases['carrier_main_10gene'].sum()):,}).
    The PRS slope SE within the carrier layer is of the order of 0.5 years/SD, so the PRS x carrier interaction term only occurs when the effect is large
    Can be detected. **If the conclusion in this section is not significant, the correct statement is \"insufficient efficacy and interaction cannot be ruled out\", not \"no interaction\". **
    The actual calculated minimum detectable effect (MDE) with 80% power is given in Section 5.

    ---

    ## 1. Queue construction

    ### Aim 4 cases queue (completely reuse `scripts/17` logic, hard assertion n= in script{EXPECTED_AIM4_N:,})

    | steps | n |
    |---|---|{flow_md}Assertion passed: n={len(cases):,}, exactly the same as script 17.

    ### Full queue (for A/A2)

    Directly adopt the modeling queue of `scripts/90_legacy_ptv_baseline.py` (`data/16_joint_model_cv_oof_person_level.csv.gz`
    person_id collection), n={len(full):,}（case={int((full['group']=='case').sum()):,}，
    control={int((full['group']=='control').sum()):,}).

    ### z Standardized caliber (two sets, do not mix with each other)

    | z variable | queue for calculation | mean | sd |
    |---|---|---|---|
    | `PRS311_z_full` / `burden_z_full` | Full queue n={len(full):,} | {mu_p_f:.5f} / {mu_b_f:.6f} | {sd_p_f:.5f} / {sd_b_f:.6f}|
    | `PRS311_z` / `burden_z` | Aim 4 cases queue n={len(cases):,} | {mu_p_c:.5f} / {mu_b_c:.6f} | {sd_p_c:.5f} / {sd_b_c:.6f}|

    A/A2 uses the whole queue z and the whole queue quantile cut point; the internal z of the cases of B/C/D/E is consistent with `scripts/17`.
    **PRS-high threshold exception** for Section C: Taken from **Full Cohort** PRS311 90th percentile ({hi_thr_d:.4f}),
    Not the quantile within cases to ensure comparability with the grouping definition of A .

    ---

    ## 2. A. 2x2 discriminant group (PRS311 decile extreme group x case/control, full queue)

    The main analysis uses deciles (top/bottom 10%) and quintiles (20/80) for sensitivity analysis.

    ### Tenth place (main analysis) each grid

    | group | PRS group | n | n_carrier | carrier rate (%) | 95% CI (Wilson) | Burden mean | Burden median |
    |---|---|---|---|---|---|---|---|{a_md}### Quintile (sensitivity) each grid

    | group | PRS group | n | n_carrier | carrier rate (%) |
    |---|---|---|---|---|{q_md}### Between-group test

    | Plan | Comparison | high rate (%) | low rate (%) | difference (pp) | rate ratio | chi2 | p |
    |---|---|---|---|---|---|---|---|{atest_md}**Reading**: Within **cases**, the PTV carrier rate in the lowest decile of PRS ({obs_row_l['carrier_rate_pct']:.2f}%)
    Significantly higher than the top PRS decile ({obs_row_h['carrier_rate_pct']:.2f}%), that is, the so-called \"patients with low PRS are more likely to rely on
    Rare variant disease.\" But within **Control** this contrast almost disappears
    ({A_tests[(A_tests['scheme'].str.startswith('decile'))&(A_tests['contrast'].str.startswith('control'))].iloc[0]['carrier_rate_high_pct']:.2f}%
    vs {A_tests[(A_tests['scheme'].str.startswith('decile'))&(A_tests['contrast'].str.startswith('control'))].iloc[0]['carrier_rate_low_pct']:.2f}%，
    p={pfmt(A_tests[(A_tests['scheme'].str.startswith('decile'))&(A_tests['contrast'].str.startswith('control'))].iloc[0]['p_value'])}).
    **This pattern of \"present in cases and absent in controls\" is exactly the textbook feature of collider bias** and cannot be directly regarded as mechanism evidence——
    So you must look at the zero model benchmark in the next section.

    For gene carrier breakdown, see `results/table/18_A_2x2_gene_level_breakdown.csv`;
    Among them{n_gene_suppressed}/{len(A_gene)}Cell n<20, has been disclosed according to CLAUDE.md rules in `n_carrier_gene_disclosed`
    The column is marked `<20` (please use this column only for external sharing).

    ---

    ## 3. A2. Additive (no interaction) zero-model simulation benchmark

    **Settings**: Take the main effects of the Aim 1 primary model (`results/table/16_joint_model_coefficients_primary.csv`)
    beta_PRS311_z={b_prs:.4f}、beta_burden_z={b_bur:.4f}, **force the interaction term to 0**, constructed on the logit scale
    `logit P(case) = alpha + beta_PRS*PRS_z + beta_burden*burden_z`. Individual PRS_z and burden_z are taken directly
    True observations (preserving the true marginal distribution and the true joint structure of the two), alpha is dichotomously calibrated to the observed prevalence{100*obs_prev:.4f}% (simulated average prevalence{100*legacy_age_contrasts_p_case.mean():.4f}%). Then press Bernoulli(p_i) to re-extract the case state,
    Compute the exact same case-high vs case-low carrier rate comparison as A in simulated data, repeat{N_SIM:,}times
    (seed={SIM_SEED}).

    | Quantity | Null model mean | Null model 95% interval | Observations |
    |---|---|---|---|
    | case-high carrier rate (%) |{a2['null_mean_carrier_rate_high_pct']:.3f} | — | {obs_row_h['carrier_rate_pct']:.3f}|
    | case-low carrier rate (%) |{a2['null_mean_carrier_rate_low_pct']:.3f} | — | {obs_row_l['carrier_rate_pct']:.3f}|
    | Difference high−low (pp) |{a2['null_mean_diff_pct_points']:+.3f} | [{a2['null_ci95_low_diff_pct_points']:+.3f}, {a2['null_ci95_high_diff_pct_points']:+.3f}] | **{obs_diff:+.3f}** |
    | Rate ratio high/low |{a2['null_mean_rate_ratio']:.3f} | [{a2['null_ci95_low_rate_ratio']:.3f}, {a2['null_ci95_high_rate_ratio']:.3f}] | **{obs_rr:.3f}** |

    Empirical two-sided p (observed vs null distribution) ={pfmt(a2['empirical_two_sided_p_vs_null'])}(Resolution lower limit 1/{N_SIM}).

    **Conclusion**: Observation comparison{verdict_null}. That is to say,{verdict_sentence}. **Purely additive effects + case selection** under the null model
    itself will produce{a2['null_mean_diff_pct_points']:+.3f}pp\'s negative within-case contrast - this is why \"Section A should only be used as a descriptive
    Quantitative basis for \"results\".

    **Both points must be stated together**:

    1. The **absolute** carrier rate of the null model does not coincide with the observation (null model case-high{a2['null_mean_carrier_rate_high_pct']:.2f}%
       vs observation{obs_row_h['carrier_rate_pct']:.2f}%). The reason why the observations are \"more extreme\" mainly comes from high PRS cases.
       The carrier rate is lower than expected from an additive model, i.e., sub-additive in direction.
    2. This is the interaction term actually fitted by the Aim 1 main model.
       (`inter_z` beta={float(coef16.loc[coef16['term']=='inter_z','beta'].iloc[0]):.4f}，
       p={pfmt(float(coef16.loc[coef16['term']=='inter_z','p_value'].iloc[0]))}, is **negative**) in the same direction——
       This simulation deliberately sets this interaction term to 0. Therefore the correct conclusion of this section is: \"Observational contrasts exceed the *no interaction* additive model
       The range that can be explained by additive selection bias, and the direction beyond is consistent with the negative (subadditive) interaction estimated by Aim 1\",
       And **not** \"discovered new mechanism evidence\". The zero-model benchmark excludes the idea that pure collider bias can explain everything.
       The weakest hypothesis, nothing more.

    ---

    ## 4. B. Early onset vs late onset (cases queue n={len(cases):,})

    Cut-off point <45 years old (determined by user on 2026-09-24; the original three cut-off points of <40/<45/<50 have been deleted).

    | Scheme | Group | n | n_carrier | carrier rate (%) | Burden mean | Burden median | PRS311_z mean (cases internal reference) | PRS311_z mean (full cohort reference) |
    |---|---|---|---|---|---|---|---|---|{b_md}Note: PRS311_z has two reference systems. The overall mean value of z within cases is always 0 by definition; the entire queue refers to z
    (`PRS311_z_fullref`, transformed using the mean/sd of the whole cohort) reflects \"the overall PRS of the cases is higher than that of the general population by about{cases['PRS311_z_fullref'].mean():.2f}SD\" fact. The only difference between the two is a linear transformation,
    Standardized effect sizes such as between-group differences and SMD were not affected at all.

    ### Test and effect size comparison

    | Scheme | n early/late | carrier rate (%) early/late | rate ratio | p(chi2) | PRS_z mean early/late | p(Welch t, PRS) | SMD(burden) | SMD(PRS) | \\|SMD_b\\|/\\|SMD_p\\| |
    |---|---|---|---|---|---|---|---|---|---|{btest_md}### **Core discovery: Asymmetry between rare and common**

    Take <45 cut-off point as an example: early-onset group carrier rate{t45['carrier_rate_early_pct']:.2f}%, late onset group{t45['carrier_rate_late_pct']:.2f}%, rate ratio **{t45['carrier_rate_ratio']:.2f}x**
    （chi2 p={pfmt(t45['p_chi2_carrier'])}); while under the same comparison, the mean value of PRS311_z is only from{t45['prs311_z_mean_late']:+.3f}change to{t45['prs311_z_mean_early']:+.3f},
    Difference{t45['prs311_z_diff']:+.3f} SD（Welch t p={pfmt(t45['p_ttest_prs_welch'])}).
    Standardized effect size comparison: **SMD(burden)={t45['smd_burden_early_vs_late']:+.3f} vs
    SMD(PRS311)={t45['smd_prs311_early_vs_late']:+.3f}, the former belongs to the latter{t45['smd_ratio_burden_over_prs']:.1f}Times**.

    **This means**: the earlier the age of onset, the patient will specifically carry more **rare components**, while
    **Not a \"common + rare mixture of two ingredients at the same time\"**. The original Aim 4 hypothesis stated that \"younger patients are more dependent on mixed mechanisms\"
    **not supported** in this data: the other half of the mix (common PRS) shows little difference between early versus late onset.
    This is an **asymmetric** conclusion, not a \"mixed\" conclusion - it should not be reported vaguely as \"younger patients have a higher genetic load\".

    ---

    ## 5. C. Literal operationalization of \"hybrid mechanism\" (four categories)

    Classification rule: PRS-high = PRS311 >= 90th percentile of the entire cohort ({hi_thr_d:.4f}); carrier = 10 gene panel
    PTV burden > 0. Category four: both/PRS-high only/carrier only/neither.

    | Plan | Group | Category | n (disclosure standard) | Accounting for this group (%) | 95% CI |
    |---|---|---|---|---|---|{c_md}| Scheme | Table | chi2 | df | p | Minimum grid n (gear) |
    |---|---|---|---|---|---|{ctest_md}**Reading**: If \"young patients are more dependent on mixed mechanisms\" is true, the proportion of **both** classes in the early-onset group should be significantly higher than that in the late-onset group.
    The actual results are shown in the table above - the increase in the proportion of the early-onset group is mainly the **carrier only** category, not the both category;
    Moreover, the absolute number of people in the both category is n<20 at most cut-off points, which is itself on the edge of the disclosure threshold and the estimate is extremely unstable.
    **Not enough to support any positive assertion about the \"hybrid mechanism\"**. This is consistent with the asymmetry findings of Section 4 .

    ---

    ## 6. D. Interaction test (the omitted part of script 17)

    | model | term | beta | robust SE (HC1) | 95% CI | p |
    |---|---|---|---|---|---|{d_md}Carrier layered PRS slopes (fitted separately for each `age ~ PRS311_z + PC1-10`) see
    The two lines starting with `D1_stratified:` in `results/table/18_D_interaction_models.csv`.

    D2 (continuous x continuous) Note: The proportion of burden==0 within cases is up to **{pct_zero_burden:.1f}%**, the interaction term is actually
    only by{int(cases['carrier_main_10gene'].sum()):,}The name is non-zero individual drive and has limited explanatory power.

    The complete coefficients for D3 (binary outcome logistic, three cut points) are shown in the same table, row starting with `D3_logistic`.

    ### Power (MDE, alpha=0.05 bilateral, 80% power)

    | model | observed beta | SE | p | MDE | observed\\|beta\\|/MDE |
    |---|---|---|---|---|---|{pow_md}**Interpretation**: Only if the absolute value of the true interaction effect >= the MDE in the table can this sample size be 80% confident of detecting it.
    Interactions below the MDE cannot be excluded from this study. Prior statement in connection with Section 0(c): non-significant results in this section always read
    \"Insufficient power\" does not read \"no interaction\".

    ---

    ## 7. E. Total-effect model with environmental/clinical covariates

    **Positioning Notes**: This is the only available alternative path** to the \"environmental exposure\" half of the Aim 4 narrative for this project. BMI, menopausal status,
    Family history of breast cancer and smoking status are baseline clinical/questionnaire variables, **not true environmental exposure measures**; the model is total-effect
    Covariate adjustment, not causal mediation analysis. Family history is particularly important to note: it is also a proxy for genetic susceptibility and is put into the model
    It may be an over-adjustment in terms of cause and effect, so a version without family history was set up for comparison.

    ### Missing covariates (cases queue n={len(cases_e):,})

    | Covariates | Missing n | Missing (%) |
    |---|---|---|{miss_md}The missing family history rate is as high as{miss_tbl[miss_tbl['covariate']=='family_history_breast_cancer']['pct_missing'].iloc[0]:.1f}%,
    Full case analysis would lose more than half the sample, so E1 (with family history) is reported side by side with E2 (without family history);
    Each specification is also accompanied by **Control fitting on the same subsample without environmental covariates** to distinguish whether the \"coefficient change\" comes from the covariate adjustment.
    Or just from subsample variation.

    | model | term | beta | robust SE | 95% CI | p | n |
    |---|---|---|---|---|---|---|{e_md}---

    ## 8. Conclusion

    1. **Aim 4 spec only partially implemented**: somatic mutation burden not analyzed at all (AoU no tumor sequencing, external WES
       not available), the statement referring to \"three mechanisms\" must be downgraded to \"two germline mechanisms\".
    2. **The main empirical finding is rare/common asymmetry, not \"mixing\"**: the earlier the age of onset, the higher the PTV carrier rate.
       Steep (<45 cut point rate ratio{t45['carrier_rate_ratio']:.2f}x), whereas PRS311 showed little difference between early and late onset
       (SMD{t45['smd_prs311_early_vs_late']:+.3f}). Effect size ratio approx.{t45['smd_ratio_burden_over_prs']:.1f}times.
       Therefore, the hypothesis that \"younger patients are more dependent on common+rare mixed mechanisms\" is **not supported**; the data supports that \"early-onset cases
       Specific enrichment of rare high-penetrance variants\". The four-category composition analysis in Section C gives the same conclusion: the early-onset group has more
       **carrier only**, not **both**.
    3. **The within-case PRS-burden negative correlation in section A cannot be read as mechanistic evidence**: Additive null model simulations of A2 show a purely additive effect
       In addition, case selection itself will induce{a2['null_mean_diff_pct_points']:+.3f}Within-case comparison of pp; observed value{verdict_null}.
    4. **Interaction test is limited by power**: carrier only accounts for cases{100*cases['carrier_main_10gene'].mean():.2f}%,
       MDE See Section 6; insignificant = cannot be excluded, does not mean not present.
    5. All cases-only comparisons are affected by collider/index-event bias, are descriptive results, and are case-control with Aim 1
       The evidence is not at the same level of evidence.

    ---

    ## 9. Output file list

    Table (`results/table/`):

    - `18_cohort_flow.csv` - Aim 4 cases queue construction process
    - `18_A_2x2_prs_extreme_cells.csv` - 2x2 cells (tenths + fifths)
    - `18_A_2x2_carrier_rate_tests.csv` —— Section A chi2 test
    - `18_A_2x2_gene_level_breakdown.csv` —— Gene carrier breakdown (including disclosure column)
    - `18_A2_additive_null_simulation.csv` - Additive null model simulation benchmark
    - `18_B_early_late_group_summary.csv` - early/late group summary (cut point <45)
    - `18_B_early_late_tests.csv` – early/late onset test versus SMD effect size
    - `18_B_carrier_rate_by_age_bin.csv` - Carrier rate binned by age of onset (bins are automatically merged by carrier and non-carrier >=20)
    - `18_C_mixture_category_composition.csv` —— Four categories of composition
    - `18_C_mixture_category_tests.csv` - 2x4 with both-vs-rest test
    - `18_D_interaction_models.csv` —— D1/D1 layered/D2/D3 full coefficient
    - `18_D_interaction_power.csv` - Interaction MDE
    - `18_E_env_covariate_missingness.csv` - Missing environmental covariates
    - `18_E_env_covariate_models.csv` - E1/E2 and same subsample control, logistic version

    Figures (`results/figures/`, all SVG, matplotlib default style, no grid):

    - `18_A_carrier_rate_2x2_bar.svg`
    - `18_A2_additive_null_simulation_hist.svg`
    - `18_B_carrier_rate_by_age_bin.svg`
    - `18_B_prs311_by_age_bin.svg`
    - `18_C_mixture_composition_stacked_bar.svg`

    Run log: `scripts/script_output_py/output_18_aim4_group_contrasts.txt`"""

    report_path = os.path.join(LEGACY_AGE_CONTRASTS_REPORT_DIR, "18_aim4_group_contrasts_summary.md")
    with open(report_path, "w") as f:
        f.write(legacy_age_contrasts_report)
    legacy_age_contrasts_log(f"Write a report:{report_path} ({os.path.getsize(report_path)/1024:.1f} KB)")

    legacy_age_contrasts_section("DONE")
    legacy_age_contrasts_log('All outputs are aggregated results, and no participant-level intermediate files are placed on disk.')

# Command-line selection of the computations defined above.
def main():
    from scripts.functions.cli import run_local_step
    computations = {
        'age-contrasts': analyze_age_group_contrasts,
        'ehr-audit': audit_early_ehr_diagnoses,
        'nine-gene-missense': analyze_nine_gene_missense_age,
        'legacy-prs-age': describe_prs_and_age,
        'legacy-burden-distribution': describe_legacy_burden,
        'legacy-group-contrasts': analyze_legacy_age_contrasts,
    }
    run_local_step('Optional age/PRS contrasts and EHR feasibility analyses', computations, default_tasks=[])

if __name__ == '__main__':
    main()

# Documentation: scripts/README.md.
import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import norm
from statsmodels.stats.proportion import proportion_confint

ROOT = Path(__file__).resolve().parents[2]
PCS = [f"PC{i}" for i in range(1, 11)]
MIN_N = 21  # Current official policy explicitly prohibits direct counts of 1–20.
OUT = ROOT / "results" / "table"
REPORT = ROOT / "results" / "report"


def read(name, **kwargs):
    return pd.read_csv(ROOT / "data" / name, **kwargs)


def releasable(counts):
    # Deliberately suppress zero cells as well, to keep table handling conservative.
    return bool((np.asarray(counts) >= MIN_N).all())


def save(frame, suffix):
    frame.to_csv(OUT / f"27_{suffix}.csv", index=False)


def table_counts(df, keys):
    return (
        df.groupby(keys)
        .size()
        .reindex(pd.MultiIndex.from_product([[0, 1]] * len(keys)), fill_value=0)
    )


def published_joint_margin_blocks(variable, published_path=None):
    """A released adjacent joint cell plus a new margin reveals a suppressed cell."""
    path = published_path or OUT / "18_C_mixture_category_composition.csv"
    if not path.exists():
        return False
    old = pd.read_csv(path)
    early = old.loc[(old.scheme == "cutoff_45") & (old.group == "early(<45)")]
    adjacent = "carrier only" if variable == "carrier" else "PRS-high only"
    n = pd.to_numeric(
        early.loc[
            early.mixture_category.isin([adjacent, "both (PRS-high AND carrier)"]), "n_disclosed"
        ],
        errors="coerce",
    )
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
        rows.append(
            {
                "variable": variable,
                "level": group,
                "n_cases": n,
                "n_early": early,
                "n_later": late,
                "early_fraction_within_group": early / n,
                "within_ci_low": within_ci[0],
                "within_ci_high": within_ci[1],
                "fraction_of_early_cases": early / early_total,
                "composition_ci_low": composition_ci[0],
                "composition_ci_high": composition_ci[1],
            }
        )
    return pd.DataFrame(rows)


def log_odds_contrast(fit, v, name):
    v = np.asarray(v, dtype=float)
    beta = float(v @ fit.params)
    se = float(np.sqrt(max(0.0, v @ fit.cov_params() @ v)))
    lo, hi = beta - norm.ppf(0.975) * se, beta + norm.ppf(0.975) * se
    return {
        "contrast": name,
        "log_OR": beta,
        "se_HC1": se,
        "OR": np.exp(beta),
        "OR_ci_low": np.exp(lo),
        "OR_ci_high": np.exp(hi),
    }


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
                ci = proportion_confint(early, early + late, method="wilson")
                comp_ci = proportion_confint(early, int(df.early.sum()), method="wilson")
                rows.append(
                    {
                        "carrier": c,
                        "high_PRS": h,
                        "n_early": early,
                        "n_later": late,
                        "early_fraction_within_group": early / (early + late),
                        "within_ci_low": ci[0],
                        "within_ci_high": ci[1],
                        "fraction_of_early_cases": early / df.early.sum(),
                        "composition_ci_low": comp_ci[0],
                        "composition_ci_high": comp_ci[1],
                    }
                )
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
    df["PRS_z"] = (df.PRS311 - mu) / sd
    df["interaction"] = df.PRS_z * df.carrier
    # PCs rescaled for numerical conditioning; genetic terms retain original units.
    scaled_pc = df[PCS].apply(lambda col: (col - col.mean()) / col.std(ddof=1))
    assert np.isfinite(scaled_pc.to_numpy()).all()
    X = sm.add_constant(
        pd.concat([df[["PRS_z", "carrier", "interaction"]], scaled_pc], axis=1), has_constant="add"
    )
    fit = sm.GLM(df.early, X, family=sm.families.Binomial()).fit(cov_type="HC1")
    assert fit.converged and np.isfinite(fit.params).all()
    vectors = {}
    for i, name in [
        (1, "PRS_per_SD_noncarriers"),
        (2, "carrier_at_primary_mean_PRS"),
        (3, "PRS_carrier_interaction"),
    ]:
        v = np.zeros(X.shape[1])
        v[i] = 1
        vectors[name] = v
    v = np.zeros(X.shape[1])
    v[1] = v[3] = 1
    vectors["PRS_per_SD_carriers"] = v
    save(
        pd.DataFrame([log_odds_contrast(fit, v, k) for k, v in vectors.items()]),
        f"{label}_early_continuous_case_only",
    )
    audit["continuous_case_only_model"] = "released_early_vs_later_case_OR_not_population_risk"
    return audit


def ehr_query(cdr):
    domain_columns = [
        ("condition_occurrence", "condition_start_date"),
        ("procedure_occurrence", "procedure_date"),
        ("drug_exposure", "drug_exposure_start_date"),
        ("measurement", "measurement_date"),
        ("visit_occurrence", "visit_start_date"),
    ]
    union = "\nUNION ALL\n".join(
        f"SELECT t.person_id, t.{col} AS event_date FROM `{cdr}.{table}` t JOIN cohort USING(person_id)"
        for table, col in domain_columns
    )
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
    dry = client.query(
        sql,
        job_config=bigquery.QueryJobConfig(
            dry_run=True, use_query_cache=False, query_parameters=params
        ),
    )
    gib = dry.total_bytes_processed / 1024**3
    print(f"EHR query dry-run: {gib:.2f} GiB; no participant records printed.", flush=True)
    if gib > maximum_gib:
        raise RuntimeError("EHR query exceeds fixed scan cap; review query before running.")
    config = bigquery.QueryJobConfig(
        query_parameters=params, maximum_bytes_billed=int(maximum_gib * 1024**3)
    )
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
    if len(s) < MIN_N:
        return None
    # Avoid exposing small missingness counts through an available n and known total.
    return {"median": s.median(), "q25": s.quantile(0.25), "q75": s.quantile(0.75)}


def span_proxy_masks(df, lookback):
    """Feasibility only: entry after washout; censor at clinical end/death/45."""
    birthday45 = df.birth_date.map(lambda d: d + pd.DateOffset(years=45))
    end_date = pd.concat([df.ehr_last_date, df.death_date, birthday45], axis=1).min(axis=1)
    entry_date = df.first_adult_ehr_date + pd.to_timedelta(lookback, unit="D")
    recorded_age = (df.first_bc_date - df.birth_date).dt.days / 365.25
    eligible = (
        (entry_date < end_date)
        & (entry_date < birthday45)
        & (df.first_bc_date.isna() | (df.first_bc_date > entry_date))
    )
    event = eligible & (df.first_bc_date <= end_date) & (recorded_age < 45)
    return eligible, event


def window_audit(full, dates, sources, cutoff):
    df = full.merge(dates, validate="one_to_one")
    datecols = [
        "birth_date",
        "first_bc_date",
        "death_date",
        "ehr_first_date",
        "ehr_last_date",
        "first_adult_ehr_date",
        "op_first_date",
        "op_last_date",
    ]
    for col in datecols:
        df[col] = pd.to_datetime(df[col])
    audit = {
        "n_source_selected": len(df),
        "other_cancer_excluded_from_source": True,
        "population_risk_ready": False,
        "risk_gate_reason": "outcome_selected_source_and_unvalidated_EHR_capture",
        "observation_period_presence": bool(df.n_observation_periods.notna().any()),
    }
    # Check event identity without publishing a count of discordant records.
    recorded_age = (df.first_bc_date - df.birth_date).dt.days / 365.25
    case = df.group == "case"
    delta = (recorded_age[case] - df.loc[case, "age_at_diagnosis"]).abs()
    audit["first_BC_dates_match_saved_case_ages"] = bool(
        delta.notna().all() and (delta < 1 / 365.25).all()
    )
    assert audit[
        "first_BC_dates_match_saved_case_ages"
    ], "BC date/age mismatch: review phenotype definition."
    if len(sources) and releasable(sources.n_participants):
        save(sources, "observation_period_provenance")
        audit["observation_period_provenance_released"] = True
    else:
        audit["observation_period_provenance_released"] = False
    df["carrier"] = (df.burden_rare_main_ptv_10gene > 0).astype(int)
    df["high"] = (df.PRS311 >= cutoff).astype(int)
    df["age_first_ehr"] = (df.ehr_first_date - df.birth_date).dt.days / 365.25
    df["age_last_ehr"] = (df.ehr_last_date - df.birth_date).dt.days / 365.25
    df["pre_BC_record_span_years"] = (df.first_bc_date - df.ehr_first_date).dt.days / 365.25
    # Medians of existing dates never imply clinical coverage between them.
    rows = []
    for group, sub in df.groupby("group"):
        for col in [
            "age_first_ehr",
            "age_last_ehr",
            "n_adult_pre45_dates",
            "n_adult_pre45_year_bins",
            "max_adult_pre45_gap_days",
            "pre_BC_record_span_years",
            "n_dates_year_before_bc",
        ]:
            summary = numeric_summary(sub[col])
            if summary:
                rows.append({"source_group": group, "measure": col, **summary})
    save(pd.DataFrame(rows), "ehr_window_distributions")
    # Mutually exclusive proxy categories, each with its complement above threshold.
    statuses = np.select(
        [
            df.age_first_ehr.isna(),
            df.age_first_ehr >= 45,
            (df.age_first_ehr < 45) & (df.age_last_ehr < 45),
        ],
        ["no_clinical_EHR", "first_EHR_at_or_after45", "EHR_ends_before45"],
        default="records_span45",
    )
    counts = (
        pd.Series(statuses)
        .value_counts()
        .reindex(
            ["no_clinical_EHR", "first_EHR_at_or_after45", "EHR_ends_before45", "records_span45"],
            fill_value=0,
        )
    )
    audit["window_partition_released"] = bool(((counts == 0) | (counts >= MIN_N)).all())
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
        if (yes == 0 or yes >= MIN_N) and (no == 0 or no >= MIN_N):
            coverage_rows.append(
                {"check": name, "n_yes": yes, "n_no": no, "fraction_yes": yes / len(df)}
            )
    save(pd.DataFrame(coverage_rows), "inferred_observation_period_checks")
    # Censor at min(last clinical event, death, 45th birthday), never max(EHR,death).
    proxy_rows = []
    event_arrays = []
    for lookback in [365, 730]:
        eligible, event = span_proxy_masks(df, lookback)
        proxy_rows.append(
            {
                "lookback_days": lookback,
                "n_span_proxy_eligible": int(eligible.sum()),
                "n_recorded_early_events": int(event.sum()),
                "n_proxy_nonevents": int((eligible & ~event).sum()),
            }
        )
        event_arrays.append(event)
        four = df.loc[eligible].assign(early=event.loc[eligible].astype(int))
        audit[f"lookback_{lookback}_four_group_event_gate"] = releasable(
            table_counts(four, ["carrier", "high", "early"])
        )
        if lookback == 365:
            all_group_counts = (
                df.groupby(["carrier", "high"])
                .size()
                .reindex(pd.MultiIndex.from_product([[0, 1], [0, 1]]), fill_value=0)
            )
            eligible_counts = (
                four.groupby(["carrier", "high"])
                .size()
                .reindex(all_group_counts.index, fill_value=0)
            )
            excluded_counts = all_group_counts - eligible_counts
            safe = releasable(eligible_counts) and bool(
                ((excluded_counts == 0) | (excluded_counts >= MIN_N)).all()
            )
            audit["one_year_proxy_four_group_denominators_release"] = safe
            if safe:
                save(
                    eligible_counts.rename_axis(["carrier", "high_PRS"]).reset_index(
                        name="n_span_proxy_eligible"
                    ),
                    "one_year_span_proxy_four_group_denominators",
                )
        audit[f"lookback_{lookback}_proxy_eligible_release"] = (
            int(eligible.sum()) >= MIN_N and int((~eligible).sum()) >= MIN_N
        )
    # Avoid revealing a small nested difference in early events across washouts.
    event_diff = int((event_arrays[0] & ~event_arrays[1]).sum())
    event_release = all(
        r["n_recorded_early_events"] >= MIN_N and r["n_proxy_nonevents"] >= MIN_N
        for r in proxy_rows
    )
    event_release &= event_diff == 0 or event_diff >= MIN_N
    eligible_diff = proxy_rows[0]["n_span_proxy_eligible"] - proxy_rows[1]["n_span_proxy_eligible"]
    eligible_release = all(audit[f"lookback_{x}_proxy_eligible_release"] for x in [365, 730])
    eligible_release &= eligible_diff == 0 or eligible_diff >= MIN_N
    if eligible_release:
        released = pd.DataFrame(proxy_rows)
        if not event_release:
            released = released[["lookback_days", "n_span_proxy_eligible"]]
        save(released, "followup_span_proxy_feasibility")
    audit["nested_proxy_event_totals_released"] = bool(event_release and eligible_release)
    audit["followup_span_proxy_totals_released"] = bool(eligible_release)
    audit["requested_additive_risk_contrast"] = "R11-R10-R01+R00"
    audit["requested_four_group_risks_and_additive_CI"] = (
        "not_estimated_failed_observation_source_and_sparse_event_gates"
    )
    return audit

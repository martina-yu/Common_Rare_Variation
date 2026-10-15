"""Functions behind 01_get_clinical_summary_table.py and 02_sample_QC.py.

The driver script tells the story; this file holds the machinery. Groups 1-4
serve 01 (the clinical summary table); group 5 serves 02 (sample-level QC):

  1. LOAD      -- BigQuery / GCS access, with a local cache so that changing the
                  table never re-runs the queries (CLAUDE.md 8: extraction and
                  analysis stay decoupled).
  2. DERIVE    -- one clinical variable per function, each returning a column
                  that is defined for every patient, so no row is ever dropped.
  3. COUNT     -- the cross-tabulation, plus a hard check that every block is
                  exhaustive.
  4. DISCLOSE  -- the small-count protection, and the rendering / writing of the
                  finished table.
  5. SAMPLE QC -- sex/age filters, the WGS join, PLINK missingness, and KING
                  relatedness, for the pipeline in 02_sample_QC.py. Every group
                  count this section prints is checked against the numbers
                  already published in results/tables/clinical_summary_table.csv
                  before it is printed (CLAUDE.md analysis_v2 #13: published
                  tables and this script's own log must be judged jointly).

Every derivation here is copied unchanged from the scripts this folder replaces
(v2_19 / v2_21), so the numbers are identical; only the packaging is new.
"""

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

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
# The project root, derived from this file's own location, so a clone works
# wherever it is checked out.
V2 = str(PROJECT_ROOT)
TABLE_DIR = os.path.join(V2, "results", "tables")
CACHE_DIR = os.path.join(V2, "data", "clinical_table_cache")

CDR = os.environ["WORKSPACE_CDR"]  # the release we are reading
PROJECT_ID = os.environ["GOOGLE_PROJECT"]
WGS_MANIFEST = "gs://vwb-aou-datasets-controlled/v9/wgs/cram/manifest.csv"
STORAGE_OPTS = {"requester_pays": True, "project": PROJECT_ID}

BC_CONCEPT_ID = 4112853  # SNOMED "Malignant tumor of breast"
BMI_CONCEPT_ID = 3038553
Q_MENOPAUSE, Q_SMK100, Q_SMKFREQ, Q_FAMHIST = (
    1585784,
    1585857,
    1585860,
    836772,
)  # menopause SNOMED code = 1585784
EARLY_ONSET_CUTOFF = 45  # years; locked before any testing (D-table)

MIN_N = 21  # the floor for an exact published count (AoU bars 1..20)
BOUND = 20  # what we admit to when we cannot be exact
SECONDARY = "suppressed (secondary)"
UNKNOWN = "Unknown / not reported"
COLS = ["early", "late", "all"]

# Receptor status is read from diagnosis codes, the only place this release
# carries it. Concept lists are v1's, unchanged.
RECEPTOR_CONCEPTS = {
    "ER_positive": [4167696, 44837878, 3526950, 3283803],
    "ER_negative": [4261933, 44834368, 3564533, 3183458],
    "PR_positive": [4219694, 3526951, 3315165, 1102831],
    "PR_negative": [4261934, 3562918, 3254239, 1102832],
    "HER2_positive": [4142116, 42542268, 42542428, 45949645],
    # No HER2-negative diagnosis concept exists beyond this one highly specific
    # term, which nobody here carries. Kept so the row is a measured zero rather
    # than an omission.
    "HER2_negative": [35624616],
}

# Category collapsing, applied identically in every column.
SEX_OTHER = "Male, intersex, or unreported"
RACE_MAP = {
    "White": "White",
    "Black or African American": "Black or African American",
    "Asian": "Asian",
}
RACE_UNKNOWN = {
    "None Indicated",
    "PMI: Skip",
    "I prefer not to answer",
    "unmapped",
    "No matching concept",
}

# The table's rows: (block label, dataframe column, level order, parent level).
# Every block is EXHAUSTIVE -- its levels sum to the column total -- so an
# unrecorded value is carried by an explicit level instead of vanishing. A
# parent makes the block nested inside one level of another block.
BLOCKS = [
    ("Sex at birth", "sex_at_birth", ["Female", SEX_OTHER], None),
    (
        "Race (self-reported)",
        "race",
        ["White", "Black or African American", "Asian", "Other", UNKNOWN],
        None,
    ),
    (
        "Menopausal status",
        "menopause",
        ["Premenopausal", "Postmenopausal", "Unknown / uncertain", UNKNOWN],
        None,
    ),
    ("Smoking status", "smoking2", ["Ever", "Never", UNKNOWN], None),
    (
        "Smoking, detail among ever smokers",
        "smoking4",
        ["Current smoker", "Former smoker", "Ever, currency unknown"],
        ("Smoking status", "Ever"),
    ),
    ("Body mass index", "bmi_recorded", ["Recorded", "Not recorded"], None),
    (
        "Family history of breast cancer",
        "famhist",
        ["Positive (non-self)", "Negative or self only", UNKNOWN],
        None,
    ),
    ("Tumour stage", "stage", ["Any staging code recorded", "No staging code recorded"], None),
    ("ER status", "ER_status", ["Positive", "Negative", UNKNOWN], None),
    ("PR status", "PR_status", ["Positive", "Negative", UNKNOWN], None),
    ("HER2 status", "HER2_status", ["Positive", "Negative", UNKNOWN], None),
]

# Footnote markers, kept because two rows are actively misleading without one:
# HER2 negative is a true measured zero about coding, not about tumours, and the
# staging row is an upper bound on availability rather than a stage count. The
# note text itself lives with the report, not here.
FOOTNOTES = {("HER2 status", "Negative"): "a", ("Tumour stage", "Any staging code recorded"): "b"}

# Levels meaning "everyone not accounted for above". When a block must give up a
# cell to protect a small count, one of these is the cheapest to lose.
RESIDUAL_LEVELS = {
    "Unknown / not reported",
    "Unknown / uncertain",
    "Not recorded",
    "No staging code recorded",
    "Ever, currency unknown",
}

_client = None
USE_CACHE = True  # --no-cache turns this off


def log(m=""):
    """Progress goes to stderr, because stdout carries the workbook bytes."""
    print(m, file=sys.stderr, flush=True)


def section(t):
    log("\n" + "=" * 78)
    log(t)
    log("=" * 78)


# ==========================================================================
# 1. LOAD
# ==========================================================================
def bq(sql):
    """Run one query against the CDR and return a DataFrame."""
    global _client
    if _client is None:
        _client = bigquery.Client(project=PROJECT_ID)
    return _client.query(sql).to_dataframe()


BC_CTE = f"""
WITH bc AS (
  SELECT co.person_id, MIN(DATE(co.condition_start_datetime)) AS first_bc_dx_date
  FROM `{CDR}.condition_occurrence` co
  JOIN `{CDR}.concept_ancestor` ca
    ON co.condition_concept_id = ca.descendant_concept_id
  WHERE ca.ancestor_concept_id = {BC_CONCEPT_ID}
  GROUP BY co.person_id
)
"""


def cached(name, build):
    """Read `name` from the cache, or build it once and store it.

    The cache lives under data/, which .gitignore covers, so nothing
    person-level can reach the repository. Delete a file to force a refresh.
    """
    os.makedirs(CACHE_DIR, exist_ok=True)
    path = os.path.join(CACHE_DIR, name + ".csv.gz")
    if USE_CACHE and os.path.exists(path):
        log(f"  cache hit: {name}")
        return pd.read_csv(path)
    log(f"  querying: {name}")
    df = build()
    df.to_csv(path, index=False)
    return df


def survey_latest(
    question_id, colname
):  # extract values from survey, logic: only take most recent(survey_datetime) answer
    """Most recent answer to one survey question, for every patient in `bc`."""
    return bq(
        BC_CTE
        + f"""
    SELECT person_id, answer AS {colname} FROM (
      SELECT b.person_id, s.answer,
             ROW_NUMBER() OVER (PARTITION BY b.person_id
                                ORDER BY s.survey_datetime DESC) AS rn
      FROM bc b JOIN `{CDR}.ds_survey` s ON b.person_id = s.person_id
      WHERE s.question_concept_id = {question_id}
    ) WHERE rn = 1
    """
    )


def load_patients():
    """Every person with a breast cancer diagnosis, with the three fields that
    are recorded for all of them: date of first diagnosis, sex at birth, race."""
    return cached(
        "patients",
        lambda: bq(
            BC_CTE
            + f"""
    SELECT b.person_id, b.first_bc_dx_date, p.birth_datetime,
           COALESCE(cs.concept_name, 'unmapped') AS sex_at_birth,
           COALESCE(cr.concept_name, 'unmapped') AS race
    FROM bc b JOIN `{CDR}.person` p USING (person_id)
    LEFT JOIN `{CDR}.concept` cs ON p.sex_at_birth_concept_id = cs.concept_id
    LEFT JOIN `{CDR}.concept` cr ON p.race_concept_id = cr.concept_id
    """
        ),
    )


def load_wgs_ids():
    """The person_ids with whole-genome sequencing, from the CRAM manifest."""
    df = cached(
        "wgs_manifest",
        lambda: pd.read_csv(WGS_MANIFEST, usecols=["person_id"], storage_options=STORAGE_OPTS),
    )
    return set(df["person_id"])


# ==========================================================================
# 2. DERIVE
# ==========================================================================
def add_wgs(df, wgs_ids):
    """Sequencing availability as a COLUMN, not a filter: this table's
    denominator is every breast cancer patient, sequenced or not."""
    df = df.copy()
    df["has_wgs"] = df["person_id"].isin(wgs_ids)
    return df


def add_age(df):
    """Age at FIRST breast cancer diagnosis = that date minus date of birth.

    Not age at enrolment and not current age. The early-onset split is made on
    this variable at the cutoff fixed before any testing. Patients younger than
    18 are kept, so that early + late = all holds exactly.
    """
    df = df.copy()
    birth = pd.to_datetime(df["birth_datetime"], utc=True).dt.tz_localize(None)
    df["age_at_diagnosis"] = (pd.to_datetime(df["first_bc_dx_date"]) - birth).dt.days / 365.25
    if df["age_at_diagnosis"].isna().any():
        raise SystemExit(
            "FATAL: an age is not computable; early + late = all no "
            "longer holds -- fix this before publishing"
        )
    df["age_group"] = np.where(df["age_at_diagnosis"] < EARLY_ONSET_CUTOFF, "early", "late")
    return df


def add_menopause(
    df,
):  # function: extract menopause status from survey (function lies in line 163: survey_latest())
    """Menopausal status, from the latest answer to the menopause survey
    question. Never answered -> "Unknown / not reported"."""
    m = cached("menopause", lambda: survey_latest(Q_MENOPAUSE, "ans"))

    def level(a):
        if pd.isna(a):
            return UNKNOWN
        if "Yes None" in a or "Yes, no menstrual periods" in a:
            return "Postmenopausal"
        if "Periods Havent Stopped" in a or a.strip().endswith(": No"):
            return "Premenopausal"
        if "Yes But Hormone" in a or "hormone" in a.lower() or "Not Sure" in a:
            return "Unknown / uncertain"
        return UNKNOWN

    m["menopause"] = m["ans"].apply(level)
    df = df.merge(m[["person_id", "menopause"]], on="person_id", how="left")
    df["menopause"] = df["menopause"].fillna(UNKNOWN)
    return df


def add_smoking(df):
    """Smoking, from two survey questions: "100 cigarettes in your life?" gives
    ever/never, and "how often now?" splits ever into current and former."""
    s1 = cached("smoking_100cigs", lambda: survey_latest(Q_SMK100, "ans_100cigs"))
    s2 = cached("smoking_freq", lambda: survey_latest(Q_SMKFREQ, "ans_freq"))
    smk = s1.merge(s2, on="person_id", how="outer")

    def level(row):
        a1, a2 = row["ans_100cigs"], row["ans_freq"]
        if pd.isna(a1) or any(k in str(a1) for k in ["Skip", "Prefer Not", "Dont Know"]):
            return UNKNOWN
        if a1.endswith(": No"):
            return "Never smoker"
        if a1.endswith(": Yes"):
            if pd.isna(a2):
                return "Ever, currency unknown"
            if "Not At All" in a2:
                return "Former smoker"
            if "Every Day" in a2 or "Some Days" in a2:
                return "Current smoker"
            return "Ever, currency unknown"
        return UNKNOWN

    smk["smoking4"] = smk.apply(level, axis=1)
    df = df.merge(smk[["person_id", "smoking4"]], on="person_id", how="left")
    df["smoking4"] = df["smoking4"].fillna(UNKNOWN)
    df["smoking2"] = np.where(
        df["smoking4"] == UNKNOWN,
        UNKNOWN,
        np.where(df["smoking4"] == "Never smoker", "Never", "Ever"),
    )
    return df


def add_famhist(df):
    """Family history of breast cancer, from the multi-select relatives
    question. Positive means a non-self relative was named; a patient who
    answered but named only themselves is negative; no answer is unknown."""
    fh = cached(
        "famhist",
        lambda: bq(
            BC_CTE
            + f"""
    SELECT b.person_id, s.answer FROM bc b
    JOIN `{CDR}.ds_survey` s ON b.person_id = s.person_id
    WHERE s.question_concept_id = {Q_FAMHIST}
    """
        ),
    )
    nonself = ["Daughter", "Father", "Grandparent", "Mother", "Sibling", "Son"]
    pos = set(fh.loc[fh["answer"].apply(lambda a: any(r in a for r in nonself)), "person_id"])
    answered = set(fh["person_id"])
    df = df.copy()
    df["famhist"] = np.where(
        df["person_id"].isin(pos),
        "Positive (non-self)",
        np.where(df["person_id"].isin(answered), "Negative or self only", UNKNOWN),
    )
    return df


def add_bmi(df):
    """BMI: the most recent recorded measurement per person, keeping only
    physiologically possible values (10-100)."""
    bmi = cached(
        "bmi",
        lambda: bq(
            BC_CTE
            + f"""
    SELECT person_id, bmi_value FROM (
      SELECT b.person_id, m.value_as_number AS bmi_value,
             ROW_NUMBER() OVER (PARTITION BY b.person_id
                                ORDER BY m.measurement_date DESC) AS rn
      FROM bc b JOIN `{CDR}.measurement` m ON b.person_id = m.person_id
      WHERE m.measurement_concept_id = {BMI_CONCEPT_ID}
        AND m.value_as_number IS NOT NULL
        AND m.value_as_number BETWEEN 10 AND 100
    ) WHERE rn = 1
    """
        ),
    )
    df = df.merge(bmi, on="person_id", how="left")
    df["bmi_recorded"] = np.where(df["bmi_value"].notna(), "Recorded", "Not recorded")
    return df


def add_receptors(df):
    """ER / PR / HER2 status from diagnosis codes.

    This measures CODING BEHAVIOUR, not tumour biology. A positive receptor
    status is a billable entity tied to a targeted therapy and does get coded; a
    negative one is not a disease entity and is essentially never entered, which
    is why HER2 negative comes out as zero. A patient carrying both a positive
    and a negative code has no determinable status and joins the unknowns.
    """
    df = df.copy()
    for label, cids in RECEPTOR_CONCEPTS.items():
        hits = cached(
            "receptor_" + label,
            lambda cids=cids: bq(
                BC_CTE
                + f"""
        SELECT DISTINCT co.person_id FROM `{CDR}.condition_occurrence` co
        JOIN bc USING (person_id)
        WHERE co.condition_concept_id IN ({','.join(str(i) for i in cids)})
        """
            ),
        )
        df["rec_" + label] = df["person_id"].isin(set(hits["person_id"]))
    for rec in ["ER", "PR", "HER2"]:
        p, n = df["rec_" + rec + "_positive"], df["rec_" + rec + "_negative"]
        df[rec + "_status"] = np.where(p & ~n, "Positive", np.where(n & ~p, "Negative", UNKNOWN))
    return df


def add_stage(df):
    """Whether ANY staging-like code is recorded -- not a stage distribution.

    Stage is not a diagnosis in this release and the registry fields that would
    carry it resolve to no matching concept. The search over concept names is
    deliberately permissive, so the result is an upper bound on availability and
    is enough to establish that stage is unusable here.
    """

    def build():
        cand = bq(
            f"""
        SELECT concept_id FROM `{CDR}.concept`
        WHERE standard_concept = 'S'
          AND (LOWER(concept_name) LIKE '%ajcc%' OR LOWER(concept_name) LIKE '%tnm%'
               OR REGEXP_CONTAINS(LOWER(concept_name),
                                  r'(^|[^a-z])stage (0|i|ii|iii|iv)([^a-z]|$)'))
        """
        )
        return bq(
            BC_CTE
            + f"""
        SELECT DISTINCT co.person_id FROM `{CDR}.condition_occurrence` co
        JOIN bc USING (person_id)
        WHERE co.condition_concept_id IN
              ({','.join(str(i) for i in cand['concept_id'].tolist()[:10000])})
        """
        )

    ids = set(cached("stage_any", build)["person_id"])
    df = df.copy()
    df["stage"] = np.where(
        df["person_id"].isin(ids), "Any staging code recorded", "No staging code recorded"
    )
    return df


def collapse_categories(df):
    """Merge the categories that would otherwise sit below the floor. Applied
    identically in every column, so no differencing route is opened."""
    df = df.copy()
    df["sex_at_birth"] = df["sex_at_birth"].where(df["sex_at_birth"] == "Female", SEX_OTHER)
    df["race"] = df["race"].apply(
        lambda r: RACE_MAP.get(r, UNKNOWN if r in RACE_UNKNOWN else "Other")
    )
    return df


# ==========================================================================
# 3. COUNT
# ==========================================================================
def split_columns(df):
    """The table's three columns: early onset, later onset, everyone."""
    return {c: (df if c == "all" else df[df.age_group == c]) for c in COLS}


def count_blocks(subs, denom):
    """counts[(block, level)][column], with exhaustiveness enforced.

    If a block's levels stop summing to their declared total, every disclosure
    argument below is void, so this fails hard rather than warning.
    """
    counts, totals = {}, {}
    for lab, col, order, parent in BLOCKS:
        for v in order:
            counts[(lab, v)] = {c: int((subs[c][col] == v).sum()) for c in COLS}
        totals[lab] = {c: sum(counts[(lab, v)][c] for v in order) for c in COLS}
        want = denom if parent is None else counts[parent]
        for c in COLS:
            if totals[lab][c] != want[c]:
                raise SystemExit(
                    f"FATAL: block '{lab}' is not exhaustive in "
                    f"column '{c}' ({totals[lab][c]} vs {want[c]})"
                )
    return counts, totals


# ==========================================================================
# 4. DISCLOSE
# ==========================================================================
def read_published(path, log_path):
    """Everything an earlier, already-public run states EXACTLY.

    Optional. Pass None, or point at files that do not exist, and this returns
    nothing: a first run of this table has published no earlier version, so
    there is no difference for a reader to take. Only a repository that has
    already released a version of this table over a different population needs
    the constraint.

    Those quantities are a constraint on this table: the two describe the same
    people under two definitions, so a value stated exactly in both states their
    difference. A cell published there as a bound pins nothing and is ignored.

    The earlier run's LOG counts too, because scripts/script_output_py/ is
    tracked by git and pushed. It states one quantity its own table withholds --
    the number of patients carrying any staging code -- and that, against the
    table's exact later-onset cell, also fixes the early-onset one. Both are
    read out of the log rather than written here as literals, so this file never
    restates a protected count in its own source.
    """
    if not path or not os.path.exists(path):
        log(
            "  no earlier published version given; the closure protects this "
            "table's own counts only"
        )
        return {}
    pub = pd.read_csv(path).fillna("")
    pub["Phenotype"] = pub["Phenotype"].replace("", pd.NA).ffill()
    non_cell = {"Number of patients", "Age at diagnosis, years"}
    colmap = {"early": "Early-onset (<45)", "late": "Late-onset (>=45)", "all": "All patients"}
    exact = {}
    for _, r in pub.iterrows():
        lab, cat = str(r["Phenotype"]).strip(), str(r["Category"])
        if lab in non_cell or cat.startswith("Mean (SD)"):
            continue
        level = cat.strip().split(" [")[0]
        for c in COLS:
            v = str(r[colmap[c]]).strip()
            m = v.split(" (")[0].replace(",", "")
            if m.isdigit():
                exact[(lab, level, c)] = int(m)

    if not log_path or not os.path.exists(log_path):
        log("  the earlier run's log was not given; only its table constrains " "this one")
        return exact
    hit = re.search(r"patients with any: ([\d,]+)", open(log_path).read())
    if not hit:
        raise SystemExit(
            "FATAL: cannot find the staging total in the published "
            "log; if its format changed, re-derive what it makes "
            "public before trusting the closure below"
        )
    stage, level = "Tumour stage", "Any staging code recorded"
    exact[(stage, level, "all")] = int(hit.group(1).replace(",", ""))
    exact[(stage, level, "early")] = exact[(stage, level, "all")] - exact[(stage, level, "late")]
    return exact


def decide_masking(counts, totals, counts_ref, totals_ref, denom, denom_ref, published):
    """Which cells cannot be stated exactly. Returns the masked set.

    Two layers, both closed to a fixed point:

      layer 0  this table's own counts. A count of 1..20 cannot be stated, and
               once a cell is hidden, at least one partner in every exact
               relation must be hidden too, or subtraction recovers it.
      layer 1  the difference against the already-published table over the
               sequenced subset. A cell can be far above the floor in BOTH and
               still be unpublishable, because the two exact values state a
               small difference -- the count of unsequenced patients at that
               level. That has to be seeded, since nothing about it is hidden.

    The exact relations are the same two in both layers: the levels of a block
    sum to its total, and early + late = all.
    """

    def n_of(g):
        return counts[(g[0], g[1])][g[2]]

    def resid(g):
        return n_of(g) - counts_ref[(g[0], g[1])][g[2]]

    rels = []
    for lab, col, order, parent in BLOCKS:
        for c in COLS:
            rels.append(([(lab, v, c) for v in order], (totals[lab][c], totals_ref[lab][c])))
        for v in order:
            rels.append(([(lab, v, c) for c in COLS], None))

    masked = {(lab, v, c) for (lab, v), cc in counts.items() for c in COLS if 0 < cc[c] < MIN_N}
    own = len(masked)
    if not published:
        log(f"  masked because their own count is under the floor: {own}")
    masked |= {g for g in published if 0 < resid(g) < MIN_N}
    if published:
        log(f"  masked because their own count is under the floor: {own}")
    if published:
        log(
            f"  masked because their difference from the published table is "
            f"under the floor: {len(masked) - own}"
        )

    def pub_exact(g):
        # A zero is printed as 0 whether masked or not: zero is not a
        # disclosure. So masking a zero protects nothing, in either layer.
        return g not in masked or n_of(g) == 0

    def resid_known(g):
        return g in published and pub_exact(g)

    def hide(candidates):
        """Give up one more cell. Residual levels ("unknown", "not recorded")
        go first: a reader already reads them as "the rest", so bounding one
        costs least. Otherwise the smallest substantive level."""
        pool = [g for g in candidates if g not in masked and n_of(g) > 0]
        if not pool:
            raise SystemExit(
                "FATAL: a relation needs another protected cell "
                "and has none available -- collapse the block"
            )
        masked.add(min(pool, key=lambda g: (g[1] not in RESIDUAL_LEVELS, n_of(g))))

    changed = True
    while changed:
        changed = False
        for cells, tot in rels:
            hidden = [g for g in cells if not pub_exact(g)]
            need = len(hidden) == 1
            if not need and tot is not None:
                left = tot[0] - sum(n_of(g) for g in cells if pub_exact(g))
                need = 0 < left < MIN_N
            if need:
                hide(cells)
                changed = True
                continue
            if not published:
                continue  # no earlier version, so no difference to hide
            hidden_r = [g for g in cells if not resid_known(g)]
            need = len(hidden_r) == 1
            if not need:
                if tot is not None:
                    left = (tot[0] - tot[1]) - sum(resid(g) for g in cells if resid_known(g))
                    need = 0 < left < MIN_N
                elif resid_known(cells[2]):  # early + late = all
                    left = resid(cells[2]) - sum(resid(g) for g in cells[:2] if resid_known(g))
                    need = 0 < left < MIN_N
            if need:
                hide([g for g in cells if resid_known(g)])
                changed = True

    # Asserted, not trusted: three separate routes have slipped past a hand
    # review in this table's history.
    for g in published:
        if resid_known(g) and 0 < resid(g) < MIN_N:
            raise SystemExit(
                f"FATAL: {g} is exact in both tables and the two "
                f"state a difference of {resid(g)}"
            )
    for cells, tot in rels:
        for label, is_known, value, total in (
            ("count", pub_exact, n_of, None if tot is None else tot[0]),
            ("residual", resid_known, resid, None if tot is None else tot[0] - tot[1]),
        ):
            hidden = [g for g in cells if not is_known(g)]
            if len(hidden) == 1:
                raise SystemExit(
                    f"FATAL: the {label} for {hidden[0]} is the only "
                    "hidden quantity in an exact relation"
                )
            if total is not None:
                left = total - sum(value(g) for g in cells if is_known(g))
                if 0 < left < MIN_N:
                    raise SystemExit(
                        f"FATAL: the hidden {label}s in {cells[0][0]}/"
                        f"{cells[0][2]} sum to {left}, itself under the floor"
                    )
    log(
        "  self-check passed: no hidden count or difference, and no group of "
        "them, is recoverable"
    )
    return masked, rels


def render_cells(counts, denom, masked, rels, counts_ref, published):
    """Turn every cell into the text the table prints.

    Exact where allowed; "<=20" for a genuinely small count; ">=n" for a large
    count that had to be protected to shield a small one; and
    "suppressed (secondary)" when even a lower bound would read as a claim that
    the count is small. A bound that is false is worse than no number at all.
    """

    def n_of(g):
        return counts[(g[0], g[1])][g[2]]

    small = {
        g
        for g in masked
        if n_of(g) > 0
        and (
            n_of(g) < MIN_N
            or (g in published and 0 < n_of(g) - counts_ref[(g[0], g[1])][g[2]] < MIN_N)
        )
    }
    cover = {
        g: max(1, len({h for cells, _ in rels if g in cells for h in cells if h in small}))
        for g in masked
    }

    rendered, n_up, n_lo, n_sup = {}, 0, 0, 0
    for (lab, v), cc in counts.items():
        for c in COLS:
            g, n = (lab, v, c), cc[c]
            if g not in masked:
                rendered[g] = f"{n:,} ({100 * n / denom[c]:.1f}%)"
            elif n == 0:
                rendered[g] = "0 (0.0%)"
            elif n < MIN_N:
                rendered[g] = f"<={BOUND} (<={100 * BOUND / denom[c]:.2f}%)"
                n_up += 1
            else:
                lo = n - BOUND * cover[g]
                if lo < MIN_N:
                    rendered[g] = SECONDARY
                    n_sup += 1
                else:
                    rendered[g] = f">={lo:,} (>={100 * lo / denom[c]:.2f}%)"
                    n_lo += 1
    log(
        f"  bounded cells: {n_up} upper, {n_lo} lower, {n_sup} suppressed; "
        f"{3 * len(counts) - n_up - n_lo - n_sup} exact"
    )
    return rendered


def build_table(subs, denom, denom_wgs, rendered):
    """Assemble the finished table: header rows, then one row per block level."""
    na = "N/A"
    rows = [
        {"Phenotype": "Number of patients", "Category": "n", **{c: f"{denom[c]:,}" for c in COLS}}
    ]
    for i, (cat, val) in enumerate([("WES", None), ("WGS", denom_wgs), ("WES overlap WGS", None)]):
        rows.append(
            {
                "Phenotype": "Available genotype data" if i == 0 else "",
                "Category": cat,
                **{
                    c: (na if val is None else f"{val[c]:,} ({100 * val[c] / denom[c]:.1f}%)")
                    for c in COLS
                },
            }
        )
    rows.append(
        {
            "Phenotype": "Age at diagnosis, years",
            "Category": "Mean (SD)",
            **{
                c: f"{pd.to_numeric(subs[c]['age_at_diagnosis']).mean():.1f} "
                f"({pd.to_numeric(subs[c]['age_at_diagnosis']).std():.1f})"
                for c in COLS
            },
        }
    )
    for lab, col, order, parent in BLOCKS:
        first = True
        for v in order:
            mark = FOOTNOTES.get((lab, v), "")
            rows.append(
                {
                    "Phenotype": lab if first else "",
                    "Category": (("   " + v) if parent else v) + (f" [{mark}]" if mark else ""),
                    **{c: rendered[(lab, v, c)] for c in COLS},
                }
            )
            first = False
        if col == "bmi_recorded":
            rows.append(
                {
                    "Phenotype": "",
                    "Category": "Mean (SD), where recorded",
                    **{
                        c: (lambda x: f"{x.mean():.1f} ({x.std():.1f})")(
                            pd.to_numeric(subs[c]["bmi_value"]).dropna()
                        )
                        for c in COLS
                    },
                }
            )
    return pd.DataFrame(rows).rename(
        columns={"early": "Early-onset (<45)", "late": "Late-onset (>=45)", "all": "All patients"}
    )


def check_against_published(out, path):
    """Compare the finished table with the version already in the repository.

    Drift is not a formatting matter: the earlier table is public, so a changed
    cell states the difference between the two. This stops rather than
    overwriting, so the change is made deliberately and its disclosure
    consequences are re-checked.
    """
    if not path or not os.path.exists(path):
        log("  no baseline given; the table is not compared against anything")
        return
    # keep_default_na=False: the WES rows hold the literal text "N/A", which
    # pandas would otherwise read back as a missing value and call a difference.
    old = pd.read_csv(path, keep_default_na=False).astype(str)
    new = out.fillna("").astype(str)
    if old.shape != new.shape or not (old.values == new.values).all():
        raise SystemExit(
            f"FATAL: the table differs from the published {os.path.basename(path)}. "
            "That is a real change, not a formatting one -- re-check the "
            "disclosure consequences, then delete or move the old file to accept it"
        )
    log(f"  identical to the published {os.path.basename(path)}")


# ---- the cross-cohort workbook -------------------------------------------
# The template carries the row labels in columns A and B and the All of Us
# header in C..E, with every data cell blank. This fills those blanks; no other
# cohort's column is touched, so a template that already holds TCGA, CPTAC or
# UKBB numbers keeps them.
STAGE_NA = "Not codeable in All of Us"
STAGE_LEVELS = ["Stage I", "Stage II", "Stage III", "Stage IV"]
WES_NOTE = "Not available, but exome smaller callset is provided for exome"
WES_ROWS = {("Available genotype data", "WES"), ("Available genotype data", "WES overlap WGS")}
BLOCK_ALIAS = {"avaliable genotype data": "Available genotype data"}
AOU_COLS = ("C", "D", "E")  # early onset, later onset, all patients
FIRST_DATA_ROW = 2  # row 1 is the header, row 2 its second line


def fill_workbook(out, template):
    """Fill the All of Us columns of the template and return the .xlsx bytes.

    Nothing is written to disk here: the bytes are handed back so the caller can
    send them to stdout or to a file. The template is opened read-only in
    effect, since it is never saved over.

    Three kinds of row do not come straight from the table:

      * WES and WES overlap WGS. All of Us has run no whole-exome sequencing --
        the release's exome product is the whole-genome callset restricted to
        the exome calling regions -- so these carry a sentence, not a number.
      * Stage I..IV. Stage is not codeable in this release, so writing the count
        of patients with ANY staging code into the Stage I row would read as a
        stage I count and is not one. They get an explicit "not codeable", and
        the unknown row carries everyone.
      * Merged cells. The template merges the four stage rows in each All of Us
        column; the value is written once, into the anchor of the merge.
    """
    vals = {
        (str(r["Phenotype"]).strip(), str(r["Category"]).split(" [")[0].strip()): [
            r["Early-onset (<45)"],
            r["Late-onset (>=45)"],
            r["All patients"],
        ]
        for _, r in out.assign(Phenotype=out["Phenotype"].replace("", pd.NA).ffill()).iterrows()
    }

    sheet = pd.read_excel(template).fillna("")
    for col in ("Phenotype", "Category"):
        if col not in sheet.columns:
            raise SystemExit(
                f"FATAL: the template has no '{col}' column; its "
                "first row must be the header of the label columns"
            )
    sheet["Phenotype"] = sheet["Phenotype"].replace("", pd.NA).ffill()
    wb = load_workbook(template)
    ws = wb[wb.sheetnames[0]]
    anchor = {
        c.coordinate: rng.coord.split(":")[0]
        for rng in ws.merged_cells.ranges
        for row in ws[rng.coord]
        for c in row
    }

    written, missing = 0, []
    for i, r in sheet.iterrows():
        block = BLOCK_ALIAS.get(str(r["Phenotype"]).strip(), str(r["Phenotype"]).strip())
        cat, row = str(r["Category"]).strip(), i + FIRST_DATA_ROW
        if not cat or cat == "nan":
            continue
        if (block, cat) in WES_ROWS:
            src = [WES_NOTE] * 3
        elif block == "Tumour stage":
            src = (
                [STAGE_NA] * 3
                if cat in STAGE_LEVELS
                else vals[("Tumour stage", "No staging code recorded")]
            )
        elif (block, cat) in vals:
            src = vals[(block, cat)]
        else:
            missing.append((block, cat))
            continue
        skipped = 0
        for col, v in zip(AOU_COLS, src):
            coord = f"{col}{row}"
            if anchor.get(coord, coord) != coord:
                skipped += 1
                continue
            ws[coord] = v
        if skipped < len(AOU_COLS):
            written += 1
    if missing:
        raise SystemExit(
            f"FATAL: {len(missing)} template rows have no counterpart "
            f"in the computed table, so they would stay blank: {missing}"
        )
    buf = io.BytesIO()
    wb.save(buf)
    log(f"  template rows filled: {written}")
    return buf.getvalue()


# ==========================================================================
# 5. SAMPLE QC (02_sample_QC.py)
# ==========================================================================
# Paths. Everything person-level lives under data/, which .gitignore covers
# (checked explicitly with `git check-ignore -v` at every write site in the
# driver, per the task's compliance requirement).
QC_DATA_DIR = os.path.join(V2, "data", "sample_qc")
QC_SCRATCH_DIR = os.path.join(V2, "data", "sample_qc_scratch")
PUBLISHED_TABLE1 = os.path.join(TABLE_DIR, "clinical_summary_table.csv")

# gs://.../exome/pgen/exome.chr{N}.{pgen,pvar,psam} -- AoU v9 srWGS, exome
# calling regions, PLINK2 binary format. Path verified 2026-09-10 (gsutil ls).
PGEN_BUCKET = "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/" "snpindel/exome/pgen"


# ---- disclosure helpers for this script's own log/tables -----------------
def disclose_n(n):
    """A single count with no exact partner to protect (own-floor rule only):
    exact if 0 or >=MIN_N, else a bound. Same floor as render_cells's own-count
    branch, just for a lone number rather than a table cell."""
    n = int(n)
    if n == 0 or n >= MIN_N:
        return f"{n:,}"
    return f"<={BOUND}"


def load_published_anchors(path=PUBLISHED_TABLE1):
    """The handful of cells in the already-published clinical table that are
    stated EXACTLY (no '>=' / '<=' / 'suppressed'), keyed the way this script
    needs them. Used to (a) sanity-check this script's own recomputation of
    the same denominator, and (b) seed the cross-table difference check
    (analysis_v2 CLAUDE.md #13): if a count computed here stands in an exact
    arithmetic relation to one of these, the difference must be 0 or >=21.
    """
    if not os.path.exists(path):
        raise SystemExit(
            f"FATAL: published anchor table not found: {path}; "
            "run 01_get_clinical_summary_table.py first so this "
            "script has something to cross-check against"
        )
    t = pd.read_csv(path).fillna("")
    t["Phenotype"] = t["Phenotype"].replace("", pd.NA).ffill()

    def exact(phenotype, category, col):
        row = t[(t["Phenotype"] == phenotype) & (t["Category"] == category)]
        if row.empty:
            raise SystemExit(f"FATAL: anchor row not found: {phenotype}/{category}")
        v = str(row.iloc[0][col]).strip()
        m = v.split(" (")[0].replace(",", "")
        if not m.isdigit():
            raise SystemExit(
                f"FATAL: anchor {phenotype}/{category}/{col} is not "
                f"published exactly ('{v}'); this script must not "
                "treat it as a known value"
            )
        return int(m)

    cols = {"early": "Early-onset (<45)", "late": "Late-onset (>=45)", "all": "All patients"}
    return {
        "denom": {c: exact("Number of patients", "n", cols[c]) for c in COLS},
        "wgs": {c: exact("Available genotype data", "WGS", cols[c]) for c in COLS},
        # Only the "all" column of the sex block is published exactly; early
        # and late are themselves bounded there (Male, intersex... is small),
        # so this script never recomputes or prints an early/late split of
        # sex-at-birth -- see the driver's note at step 1.
        "female_all": exact("Sex at birth", "Female", cols["all"]),
    }


def check_diff(new_exact, published_exact, label):
    """One count computed fresh here, and one already public: their difference
    states a quantity nobody agreed to publish (analysis_v2 CLAUDE.md #13,
    point 7). 0 is fine (nothing removed); >=21 is fine (safely large); 1..20
    is not, and is fatal rather than silently printed.
    """
    diff = published_exact - new_exact
    if 0 < diff < MIN_N:
        raise SystemExit(
            f"FATAL: {label} -- new count {new_exact:,} vs "
            f"published {published_exact:,} differ by {diff}, "
            "which is in 1..20; this would disclose the removed "
            "group. Bound the new count instead of printing it."
        )
    log(
        f"  cross-check ok: {label} differs from published by "
        f"{'0' if diff == 0 else f'{diff:,} (>=21)'}"
    )
    return diff


# ---- steps 1-2 and the WGS join -------------------------------------------
def compute_age_soft(df):
    """Age at first BC diagnosis, same formula as add_age, but WITHOUT the
    hard failure on a missing value: this script's whole point is to find out
    whether any sample lacks it, not to assume the answer.
    """
    df = df.copy()
    birth = pd.to_datetime(df["birth_datetime"], utc=True).dt.tz_localize(None)
    df["age_at_diagnosis"] = (pd.to_datetime(df["first_bc_dx_date"]) - birth).dt.days / 365.25
    return df


def qc_step1_female(df, anchors):
    """Step 1: keep sex_at_birth == 'Female' only.

    Early/late are NOT reported for this step. The published clinical table
    already states the exact 'all' total for Female (15,415) but gives only
    lower bounds for the early and late cells of that same row, because the
    complementary "Male, intersex, or unreported" cells are themselves small
    enough to need protection in at least one column. Recomputing the exact
    early/late split here, from a fresh query, would state the very numbers
    that table deliberately declined to state -- same underlying people, same
    exact quantity, different route to it. So this step is reported as one
    total only.
    """
    out = df[df["sex_at_birth"] == "Female"].copy()
    check_diff(
        len(out),
        anchors["female_all"],
        "step1 Female-only vs published " "Sex at birth / Female / All patients",
    )
    log(f"  step 1 (sex_at_birth == Female): {len(out):,} remain " f"(of {len(df):,})")
    return out


def qc_step2_age(df):
    """Step 2: keep samples with a computable age at diagnosis.

    01_get_clinical_summary_table.py's add_age() already asserts, over the
    FULL 15,812-patient denominator, that every patient has a computable age
    (it raises SystemExit otherwise, and that script has run successfully).
    A subset of a set with zero missing values also has zero missing values,
    so this step is expected to remove nobody; it is still computed and
    checked explicitly rather than assumed, since the task asks for the
    removal count at this step by name.
    """
    out = compute_age_soft(df)
    n_before = len(out)
    out = out[out["age_at_diagnosis"].notna()].copy()
    n_removed = n_before - len(out)
    log(
        f"  step 2 (age at diagnosis available): {len(out):,} remain "
        f"(of {n_before:,}; removed {disclose_n(n_removed)})"
    )
    return out


def qc_wgs_join(df, wgs_ids, anchors):
    """The structurally-required step between steps 1-2 and genotype QC:
    samples must actually have WGS data before missingness/relatedness can be
    computed on them at all. Reported as its own explicit flow line.
    """
    out = df.copy()
    out["has_wgs"] = out["person_id"].isin(wgs_ids)
    n_before = len(out)
    out = out[out["has_wgs"]].copy()
    log(f"  WGS join: {len(out):,} of {n_before:,} have whole-genome data")
    # Sanity bound only (not an exact-relation check): this count must be far
    # from the published all-patient WGS total in either direction, since the
    # non-Female group is itself tiny (397 of 15,812) -- if it ever came out
    # implausibly close, that would mean the join is wrong, not that a real
    # disclosure risk appeared.
    if not (anchors["wgs"]["all"] - 500 <= len(out) <= anchors["wgs"]["all"]):
        raise SystemExit(
            f"FATAL: Female&WGS count {len(out):,} is outside the "
            f"plausible range given published WGS-all "
            f"{anchors['wgs']['all']:,}; check the join before "
            "trusting this number"
        )
    return out


def write_keep_list(person_ids, path):
    """A plink2 --keep file: one #IID per line, matching the psam's ID scheme
    (checked against the actual psam header before this is used -- the psam
    for this release carries only '#IID', no separate FID).

    Person-level (one row per person_id), so this goes under data/ and is
    checked with `git check-ignore -v` at the call site, never under results/.
    """
    os.makedirs(os.path.dirname(path), exist_ok=True)
    pd.DataFrame({"#IID": sorted(int(i) for i in person_ids)}).to_csv(path, index=False, sep="\t")
    return path


# ---- step 3: PLINK missingness on AoU v9 exome pgen -----------------------
def pgen_gcs_paths(chrom):
    base = f"{PGEN_BUCKET}/exome.chr{chrom}"
    return {"pgen": f"{base}.pgen", "pvar": f"{base}.pvar", "psam": f"{base}.psam"}


def gcs_size(path):
    """Byte size of one GCS object, via `gsutil -u PROJECT ls -l` (read-only)."""
    r = subprocess.run(
        ["gsutil", "-u", PROJECT_ID, "ls", "-l", path], capture_output=True, text=True, check=True
    )
    line = [l for l in r.stdout.splitlines() if path in l][0]
    return int(line.split()[0])


def download_pgen(chrom, dest_dir=QC_SCRATCH_DIR):
    """Download one chromosome's pgen/pvar/psam triple to a local scratch dir.

    Read-only `gsutil cp` from AoU's own controlled-tier bucket to local disk
    inside this analysis environment -- not an egress of participant-level
    data out of the workspace, and required because plink2 cannot operate on
    gs:// paths directly. Deleted by cleanup_pgen() right after use, never
    committed (git check-ignore -v is run at the driver's write sites).
    """
    os.makedirs(dest_dir, exist_ok=True)
    paths = pgen_gcs_paths(chrom)
    local = {}
    total_bytes = 0
    for kind, gpath in paths.items():
        size = gcs_size(gpath)
        total_bytes += size
        lpath = os.path.join(dest_dir, os.path.basename(gpath))
        log(f"  downloading {gpath} ({size / 1e9:.2f} GB) ...")
        # -q suppresses gsutil's progress bar: this script's stdout/stderr is
        # redirected into a git-tracked log, and the spinner writes tens of KB
        # of carriage-return noise per chromosome.
        subprocess.run(["gsutil", "-q", "-u", PROJECT_ID, "-m", "cp", gpath, lpath], check=True)
        local[kind] = lpath
    log(f"  chr{chrom} download total: {total_bytes / 1e9:.2f} GB")
    return os.path.join(dest_dir, f"exome.chr{chrom}"), total_bytes


def cleanup_pgen(prefix):
    """Delete one chromosome's local pgen/pvar/psam (and any plink2 outputs
    sharing the prefix) once its missingness has been extracted. Disk here is
    the scarce resource (~41 GB free vs. 124 GB of pgen across all chroms), so
    this MUST run before the next chromosome downloads."""
    removed = []
    d = os.path.dirname(prefix)
    b = os.path.basename(prefix)
    for f in os.listdir(d):
        if f.startswith(b):
            fp = os.path.join(d, f)
            os.remove(fp)
            removed.append(f)
    log(f"  cleaned up {len(removed)} local files for {b}")
    return removed


def plink_missing(pfile_prefix, keep_path, out_prefix, threads=4):
    """Run `plink2 --pfile ... --keep ... --missing --out ...` and return the
    exact command, its captured stdout/stderr, and the parsed .imiss table.
    """
    cmd = [
        "plink2",
        "--pfile",
        pfile_prefix,
        "--keep",
        keep_path,
        "--missing",
        "--threads",
        str(threads),
        "--out",
        out_prefix,
    ]
    log("  $ " + " ".join(cmd))
    r = subprocess.run(cmd, capture_output=True, text=True)
    log(r.stdout)
    if r.returncode != 0:
        log(r.stderr)
        raise SystemExit(f"FATAL: plink2 exited {r.returncode}")
    imiss_path = out_prefix + ".smiss"
    if not os.path.exists(imiss_path):
        imiss_path = out_prefix + ".imiss"  # older plink2 naming
    imiss = pd.read_csv(imiss_path, sep=r"\s+")
    return " ".join(cmd), r.stdout + r.stderr, imiss


# ---- step 3, scaled: all 22 autosomes, resumable -------------------------
# chrX/chrY/chrM are deliberately excluded (locked with the coordinator
# 2026-09-10): the cohort is all-female so chrY is near-universally no-call
# and would inflate F_MISS for reasons unrelated to sample quality; chrX
# ploidy/hardcall handling differs and makes OBS_CT vary by sample, breaking
# the "OBS_CT constant per chromosome" assertion below; chrM is a separate
# callset entirely.
AUTOSOMES = [str(i) for i in range(1, 23)]
SMISS_DIR = os.path.join(QC_DATA_DIR, "smiss")

# gsutil -u $GOOGLE_PROJECT ls -l .../pgen/*.pgen, measured 2026-09-10,
# autosomes only (chr1-chr22). Reused for the pre-run size estimate; the
# per-chromosome size actually downloaded is re-measured every run.
AUTOSOME_PGEN_BYTES = {
    1: 15025433749,
    2: 8150068204,
    3: 5790328100,
    4: 4307759687,
    5: 5266710641,
    6: 5848308912,
    7: 7988828410,
    8: 5143176157,
    9: 5077496311,
    10: 5135107311,
    11: 6135240523,
    12: 5903215893,
    13: 2112728101,
    14: 3969181112,
    15: 7503604765,
    16: 7094402826,
    17: 7615803723,
    18: 1944650240,
    19: 9629059502,
    20: 2573317249,
    21: 2310777959,
    22: 3323118317,
}


def load_cached_smiss(chrom):
    """A previously-cached .smiss for one chromosome, if present and it
    actually parses -- the resumability contract for the full-autosome run.
    Returns None (not an exception) on any problem, so the caller just
    re-downloads and re-runs plink2 for that chromosome.
    """
    path = os.path.join(SMISS_DIR, f"chr{chrom}.smiss")
    if not os.path.exists(path):
        return None
    try:
        df = pd.read_csv(path, sep="\t")
        if not {"#IID", "MISSING_CT", "OBS_CT"}.issubset(df.columns) or len(df) == 0:
            log(
                f"  chr{chrom}: cache at {path} exists but does not parse "
                "as expected -- re-running"
            )
            return None
        return df
    except Exception as e:
        log(f"  chr{chrom}: cache at {path} failed to read ({e}) -- re-running")
        return None


def run_chrom_missing(chrom, keep_path):
    """Full step-3 pipeline for ONE chromosome, resumable:
      1. if data/sample_qc/smiss/chr{N}.smiss already exists and parses, skip
         straight to returning it (no download, no plink2).
      2. else: download this chromosome's pgen/pvar/psam, run plink2
         --missing, cache the .smiss under data/sample_qc/smiss/, then DELETE
         the pgen/pvar/psam immediately (disk is the scarce resource: at most
         one chromosome's pgen on disk at a time).
    Asserts MISSING_CT/OBS_CT reproduces F_MISS, and that OBS_CT (the site
    count) is constant across samples for this chromosome, as it must be on
    an autosome with no per-sample ploidy variation.
    """
    os.makedirs(SMISS_DIR, exist_ok=True)
    cache_path = os.path.join(SMISS_DIR, f"chr{chrom}.smiss")

    cached_df = load_cached_smiss(chrom)
    if cached_df is not None:
        n_sites = int(cached_df["OBS_CT"].iloc[0])
        assert (
            cached_df["OBS_CT"] == n_sites
        ).all(), f"chr{chrom}: cached OBS_CT not constant across samples"
        log(
            f"  chr{chrom}: RESUMED from cache ({cache_path}), "
            f"{len(cached_df):,} samples, {n_sites:,} sites -- "
            "no download, no plink2 this run"
        )
        return {
            "chrom": chrom,
            "smiss": cached_df,
            "n_sites": n_sites,
            "t_download": None,
            "t_plink": None,
            "cached": True,
        }

    t0 = time.time()
    prefix, total_bytes = download_pgen(chrom)
    t_download = time.time() - t0

    t1 = time.time()
    out_prefix = os.path.join(QC_SCRATCH_DIR, f"chr{chrom}_full")
    cmd, raw_out, smiss = plink_missing(prefix, keep_path, out_prefix)
    t_plink = time.time() - t1

    if "F_MISS" not in smiss.columns:
        smiss["F_MISS"] = smiss["MISSING_CT"] / smiss["OBS_CT"]
    assert (
        abs(smiss["MISSING_CT"].iloc[0] / smiss["OBS_CT"].iloc[0] - smiss["F_MISS"].iloc[0]) < 1e-6
    ), f"chr{chrom}: F_MISS is not MISSING_CT/OBS_CT; re-check this plink2's .smiss spec"
    n_sites = int(smiss["OBS_CT"].iloc[0])
    assert (
        smiss["OBS_CT"] == n_sites
    ).all(), f"chr{chrom}: OBS_CT is not constant across samples -- unexpected on an autosome"

    smiss.to_csv(cache_path, sep="\t", index=False)
    log(
        f"  chr{chrom}: cached .smiss -> {cache_path} ({len(smiss):,} samples, "
        f"{n_sites:,} sites)"
    )
    log(
        f"  chr{chrom}: download {t_download:.1f}s for {total_bytes / 1e9:.2f} GB, "
        f"plink2 --missing {t_plink:.1f}s"
    )

    # Delete the pgen/pvar/psam (and any stray plink2 outputs) NOW, before the
    # next chromosome downloads -- peak disk is one chromosome at a time.
    cleanup_pgen(prefix)
    for ext in (".smiss", ".imiss", ".vmiss", ".log"):
        p = out_prefix + ext
        if os.path.exists(p):
            os.remove(p)

    return {
        "chrom": chrom,
        "smiss": smiss,
        "cmd": cmd,
        "n_sites": n_sites,
        "t_download": t_download,
        "t_plink": t_plink,
        "total_bytes": total_bytes,
        "cached": False,
    }


def exome_wide_fmiss(chrom_results, n_expected_samples):
    """Exome-wide F_MISS per sample: sum MISSING_CT across all 22
    chromosomes' .smiss, sum OBS_CT across the same 22, THEN divide -- one
    division per sample, not an average of 22 per-chromosome F_MISS values
    (which would silently reweight toward small chromosomes).

    Asserts every one of n_expected_samples appears in all 22 files (an
    outer merge with a post-hoc completeness check, so a sample missing from
    even one chromosome raises loudly instead of silently dropping out of
    the sum).
    """
    merged = None
    for r in chrom_results:
        c = r["chrom"]
        s = r["smiss"][["#IID", "MISSING_CT", "OBS_CT"]].rename(
            columns={"MISSING_CT": f"MISS_{c}", "OBS_CT": f"OBS_{c}"}
        )
        merged = s if merged is None else merged.merge(s, on="#IID", how="outer")

    n_incomplete = int(merged.isna().any(axis=1).sum())
    assert n_incomplete == 0, (
        f"exome_wide_fmiss: {n_incomplete} samples are missing from at least "
        "one of the 22 chromosome .smiss files"
    )
    assert len(merged) == n_expected_samples, (
        f"exome_wide_fmiss: {len(merged)} distinct samples across all "
        f"chromosomes, expected {n_expected_samples}"
    )

    miss_cols = [c for c in merged.columns if c.startswith("MISS_")]
    obs_cols = [c for c in merged.columns if c.startswith("OBS_")]
    assert len(miss_cols) == len(obs_cols) == len(chrom_results) == 22

    out = merged[["#IID"]].copy()
    out["MISSING_CT_total"] = merged[miss_cols].sum(axis=1)  # calculate for missing genotype rate
    out["OBS_CT_total"] = merged[obs_cols].sum(axis=1)
    out["F_MISS"] = out["MISSING_CT_total"] / out["OBS_CT_total"]
    return out


# ---- step 3 figure --------------------------------------------------------
FIG_DIR = os.path.join(V2, "results", "figures")


def plot_fmiss_distribution(fmiss, thresh, path_prefix, source_label):
    """Histogram of exome-wide per-sample F_MISS.

    DISCLOSURE: a histogram bar is a group count, so a bin holding 1..20
    samples would state a small count graphically -- the same thing the n>=20
    rule forbids in a table. The right tail is therefore merged into a single
    ">= x" bin, widened until every plotted bin is either empty or >= MIN_N.
    The merge point and the number of merged samples are returned so the
    driver can log what was done.

    Matplotlib defaults, no grid, no styling beyond what the data needs.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    v = np.sort(np.asarray(fmiss["F_MISS"], dtype=float))
    edges = np.arange(0.014, v.max() + 0.001, 0.001)
    counts, edges = np.histogram(v, bins=edges)

    # Fold each tail inward until that tail bar itself clears the floor. Both
    # tails need it: the sparse bins are at the two ends, the mode holds
    # thousands. Interior bins are then asserted clear rather than assumed.
    lo = 0
    while lo < len(counts) - 1 and counts[: lo + 1].sum() < MIN_N:
        lo += 1
    hi = len(counts)
    while hi > lo + 1 and counts[hi - 1 :].sum() < MIN_N:
        hi -= 1
    interior = counts[lo + 1 : hi - 1]
    assert not (
        (interior > 0) & (interior < MIN_N)
    ).any(), "an interior histogram bin holds 1..20 samples; widen the bins"

    plot_edges = [edges[0]] + list(edges[lo + 1 : hi]) + [edges[-1]]
    plot_counts = (
        [counts[: lo + 1].sum()] + list(counts[lo + 1 : hi - 1]) + [counts[hi - 1 :].sum()]
    )
    merge_from = edges[hi - 1]
    n_merged = int(counts[hi - 1 :].sum())

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.bar(
        plot_edges[:-1],
        plot_counts,
        width=np.diff(plot_edges),
        align="edge",
        edgecolor="white",
        linewidth=0.3,
    )
    ax.set_xlabel("Exome-wide missing genotype rate (F_MISS)")
    ax.set_ylabel("Number of samples")
    ax.set_title(f"Per-sample missingness, {len(v):,} samples\n{source_label}", fontsize=10)
    ax.annotate(
        f"exclusion threshold F_MISS $\\geq$ {thresh:.0%}\n"
        f"maximum observed {v.max():.4f} — no sample excluded",
        xy=(0.97, 0.85),
        xycoords="axes fraction",
        ha="right",
        fontsize=9,
    )
    ax.annotate(
        f"$\\geq${merge_from:.3f} merged",
        xy=(merge_from, 0),
        xytext=(merge_from, max(plot_counts) * 0.12),
        fontsize=8,
        ha="center",
        arrowprops=dict(arrowstyle="-", lw=0.6),
    )
    fig.tight_layout()

    written = {}
    for ext in ("svg", "png", "jpg"):
        p = f"{path_prefix}.{ext}"
        if ext == "jpg":
            fig.savefig(p, pil_kwargs={"quality": 85, "subsampling": 2})
        else:
            fig.savefig(p)
        written[ext] = os.path.getsize(p)
    plt.close(fig)
    return written, merge_from, n_merged


# ==========================================================================
# step 4: relatedness / duplicate filter
# ==========================================================================
# SOURCE DECISION (user, 2026-09-14): use All of Us's own pre-computed
# pairwise kinship rather than re-running KING ourselves.
#   - AoU estimates kinship with Hail's pc_relate over the full srWGS
#     callset and RELEASES ONLY PAIRS WITH kin > 0.1 (file below, 1.8 MB).
#   - The released truncation at 0.1 is therefore the operative cut. It sits
#     INSIDE KING's 2nd-degree band [0.0884, 0.177], near its lower edge, so
#     "kin > 0.1" means "duplicates + 1st + (almost all) 2nd degree".
#     Mirroring the observed 2nd-degree peak about its mode puts the loss at
#     ~4.4% of 2nd-degree pairs AoU-wide, i.e. <1 pair at this cohort's size.
#   - Same convention as the AoU flagship genomic paper (Nature 2024,
#     PMC10937371: "maximally unrelated set ... for kinship scored greater
#     than 0.1") and the AoU polygenic-prediction paper (PMC11326295: "The
#     estimates of kinship coefficient were provided by AoU", threshold 0.1).
#   - 3rd degree is NOT OBSERVABLE in this file; those cells stay blank
#     rather than being reported as zero.
RELATEDNESS_GCS = (
    "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/"
    "snpindel/aux/relatedness/samples_relatedness.tsv"
)

# KING's published bands (kingrelatedness.com/manual.shtml): ">0.354,
# [0.177, 0.354], [0.0884, 0.177] and [0.0442, 0.0884] corresponds to
# duplicate/MZ twin, 1st-degree, 2nd-degree, and 3rd-degree relationships".
# KIN_2ND is 0.1 not 0.0884 purely because that is where AoU truncated.
KIN_DUP, KIN_1ST, KIN_2ND = 0.354, 0.177, 0.10

REL_CLASSES = ["Duplicate/Monozygotic twins", "1st degree", "2nd degree", "3rd degree", "Unrelated"]
REL_EXPECTED_PHI = {
    "Duplicate/Monozygotic twins": "0.5",
    "1st degree": "0.25",
    "2nd degree": "0.125",
    "3rd degree": "0.0625",
    "Unrelated": "0",
}
REL_KING_RANGE = {
    "Duplicate/Monozygotic twins": "> 0.354",
    "1st degree": "0.177 - 0.354",
    "2nd degree": "0.1 - 0.177",
    "3rd degree": "not provided",
    "Unrelated": "not provided",
}


def relationship_class(kin):
    """Map one kinship coefficient onto KING's relationship bands."""
    if kin > KIN_DUP:
        return "Duplicate/Monozygotic twins"
    if kin >= KIN_1ST:
        return "1st degree"
    return "2nd degree"  # the file cannot contain anything below 0.1


FLAGGED_SAMPLES_GCS = (
    "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/" "aux/qc/flagged_samples.tsv"
)


def load_flagged_samples():
    """The person_ids All of Us flags as population outliers.

    SEMANTICS (checked against the file itself, 2026-09-21): this file is a
    LIST, not a pass/fail table -- it contains only the flagged samples
    (1,012 rows in v9), so membership IS the flag and there is no boolean
    column to test. Every row additionally carries eight
    `fail_<metric>_residual` columns (lowercase 'true'/'false' STRINGS, not
    booleans) and a `qc_metrics_filters` summary that lists exactly the
    metrics with 'true'; it is never empty, i.e. every listed sample fails at
    least one metric. The flag is the OR of those eight: a residual at least
    eight MADs from the median, after regressing out sixteen PCs, on SNP /
    insertion / deletion counts, ins:del ratio, Ti/Tv, SNP and indel
    het:hom ratios, and singleton count.

    WE remove them; All of Us does NOT ask us to. Its report flags rather
    than removes "since we cannot validate whether samples (especially
    population outliers) are problematic or are just a part of a
    poorly-sampled ancestry" -- the opposite of its wording for variant
    sites, which it does recommend excluding. The justification for removing
    them here is ours and belongs in the methods: this analysis is restricted
    to the European-ancestry group, for which that caveat is weakest, and the
    most frequently failed metric (indel het:hom ratio, 822/1,012) bears
    directly on the rare-variant burden endpoint.

    Read with `gsutil cat` into memory -- never copied to disk (CLAUDE.md
    allows byte-range/stream reads, not cp of participant-level data).
    """
    import csv as _csv
    import io as _io

    cmd = ["gsutil", "-q", "-u", os.environ["GOOGLE_CLOUD_PROJECT"], "cat", FLAGGED_SAMPLES_GCS]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"FATAL: cannot read {FLAGGED_SAMPLES_GCS}: {r.stderr[:400]}")
    rows = list(_csv.reader(_io.StringIO(r.stdout), delimiter="\t"))
    hdr, body = rows[0], rows[1:]
    if hdr[0] != "s":
        raise SystemExit(f"FATAL: unexpected first column in flagged_samples.tsv: {hdr[0]!r}")
    fi = hdr.index("qc_metrics_filters")
    n_empty = sum(1 for x in body if not x[fi].strip() or x[fi].strip() == "[]")
    if n_empty:
        raise SystemExit(
            f"FATAL: {n_empty} rows in flagged_samples.tsv have an "
            "empty qc_metrics_filters; the membership==flagged "
            "assumption no longer holds"
        )
    ids = {int(x[0]) for x in body}
    log(
        f"  AoU flagged_samples.tsv: {len(body):,} rows, {len(ids):,} distinct "
        "person_ids (membership == flagged; no PASS rows in this file)"
    )
    from collections import Counter as _C

    top = _C(x[fi] for x in body).most_common(3)
    for k, v in top:
        log(f"    {v:>5} x {k[:64]}")
    return ids


def load_aou_relatedness(dest_dir=QC_DATA_DIR):
    """The AoU pairwise kinship table, cached locally (person-level -> data/).

    Columns are i.s, j.s, kin. There is NO IBS0 / HETHET / NSNP column, so the
    usual kinship-vs-IBS0 scatter cannot be drawn and PO cannot be told apart
    from FS; contamination is instead read off AoU's own verify_bam_id2 /
    dragen contamination metrics if it is ever needed.
    """
    os.makedirs(dest_dir, exist_ok=True)
    path = os.path.join(dest_dir, "aou_samples_relatedness.tsv")
    if not os.path.exists(path):
        subprocess.run(["gsutil", "-q", "-u", PROJECT_ID, "cp", RELATEDNESS_GCS, path], check=True)
    rel = pd.read_csv(path, sep="\t")
    assert list(rel.columns) == ["i.s", "j.s", "kin"], rel.columns
    assert rel["kin"].min() > KIN_2ND, (
        "the released table is supposed to be truncated at kin > 0.1; a value "
        f"at or below it ({rel['kin'].min()}) means the file changed and the "
        "0.1-vs-0.0884 reasoning above has to be revisited"
    )
    return rel, path


def split_relatedness(rel, cohort_ids):
    """Partition the AoU pairs against our cohort.

    BOTH ends in the cohort -> a related pair we must break.
    ONE end in the cohort   -> the relative is a control / non-BC participant,
                               outside the analysis set. Such a pair induces
                               no correlation WITHIN our data and is NOT
                               removed. Getting this wrong is the single
                               easiest way to delete thousands of samples for
                               nothing, so it is asserted, not assumed.
    """
    ids = set(int(i) for i in cohort_ids)
    in_i = rel["i.s"].isin(ids)
    in_j = rel["j.s"].isin(ids)
    both = rel[in_i & in_j].copy()
    one = rel[in_i ^ in_j].copy()
    assert len(both) + len(one) <= len(rel)
    return both, one


def build_adjacency(pairs):
    """{sample: {relative: kin}} over the within-cohort pairs."""
    adj = collections.defaultdict(dict)
    for a, b, k in zip(pairs["i.s"], pairs["j.s"], pairs["kin"]):
        adj[int(a)][int(b)] = float(k)
        adj[int(b)][int(a)] = float(k)
    return adj


def greedy_unrelated(adj, fmiss):
    """Greedy maximal independent set: which samples to DROP.

    Removal order, locked before looking at any downstream result:
      1. most remaining relatives  (breaks 3+ clusters with one deletion)
      2. higher F_MISS             (the rule used by the reference papers:
                                    "removed a sample from each such pair
                                    that had the highest fraction of missing
                                    data")
      3. larger person_id          (pure tie-break, for determinism)
    Degrees are recomputed after every removal, so a sample whose only
    relative has already gone is never removed.

    plink2's own --king-cutoff is the same family of greedy algorithm ("this
    maximum independent set problem is NP-hard, so we use a greedy algorithm
    which does not guarantee an optimal result") but gives no way to say
    WHICH member to keep, which is why the rule is implemented here.
    """
    removed, order = set(), []
    while True:
        live = [
            (len(set(adj[v]) - removed), fmiss[v], -v)
            for v in adj
            if v not in removed and (set(adj[v]) - removed)
        ]
        if not live:
            break
        live.sort(reverse=True)
        n_rel, fm, neg_id = live[0]
        victim = -neg_id
        removed.add(victim)
        order.append(
            {
                "person_id": victim,
                "n_relatives_at_removal": n_rel,
                "F_MISS": fm,
                "removal_round": len(order) + 1,
            }
        )
    return removed, pd.DataFrame(order)


def relationship_tables(cohort_ids, adj, pairs_both, is_early):
    """The four internal tables. All of them live under data/ and are never
    published, so no n>=20 disclosure handling is applied here (user
    instruction, 2026-09-14); the pushable log must not restate their cells.
    """
    ids = [int(i) for i in cohort_ids]
    closest = {v: relationship_class(max(adj[v].values())) for v in adj}

    # t1: one row per sample, assigned to its CLOSEST relationship, so the
    # rows partition the cohort exactly and Early + Late = All.
    rows = []
    for cls in REL_CLASSES:
        if cls == "Unrelated":
            members = [v for v in ids if v not in closest]
        elif cls == "3rd degree":
            rows.append(
                {
                    "Relationship": cls,
                    "Expected phi": REL_EXPECTED_PHI[cls],
                    "KING Range": REL_KING_RANGE[cls],
                    "All Patients": "",
                    "Early-Onset": "",
                    "Late-Onset": "",
                }
            )
            continue
        else:
            members = [v for v in ids if closest.get(v) == cls]
        e = sum(1 for v in members if is_early.get(v))
        rows.append(
            {
                "Relationship": cls,
                "Expected phi": REL_EXPECTED_PHI[cls],
                "KING Range": REL_KING_RANGE[cls],
                "All Patients": len(members),
                "Early-Onset": e,
                "Late-Onset": len(members) - e,
            }
        )
    t1 = pd.DataFrame(rows)
    tot = t1[t1["All Patients"] != ""]
    t1.loc[len(t1)] = {
        "Relationship": "Total",
        "Expected phi": "",
        "KING Range": "",
        "All Patients": tot["All Patients"].sum(),
        "Early-Onset": tot["Early-Onset"].sum(),
        "Late-Onset": tot["Late-Onset"].sum(),
    }

    # t1b: the same bands counted as PAIRS, KING's native unit. A pair can
    # straddle the onset cut, hence the explicit discordant column.
    rows = []
    for cls in REL_CLASSES[:3]:
        s = pairs_both[pairs_both["kin"].map(relationship_class) == cls]
        ee = sum(
            1 for a, b in zip(s["i.s"], s["j.s"]) if is_early.get(int(a)) and is_early.get(int(b))
        )
        ll = sum(
            1
            for a, b in zip(s["i.s"], s["j.s"])
            if not is_early.get(int(a)) and not is_early.get(int(b))
        )
        rows.append(
            {
                "Relationship": cls,
                "Pairs": len(s),
                "Both early": ee,
                "Both late": ll,
                "Discordant": len(s) - ee - ll,
            }
        )
    t1b = pd.DataFrame(rows)
    t1b.loc[len(t1b)] = ["Total", *t1b[["Pairs", "Both early", "Both late", "Discordant"]].sum()]

    # t2: "number of relations".
    counts = collections.Counter(len(v) for v in adj.values())
    t2 = pd.DataFrame(sorted(counts.items()), columns=["n_relatives", "n_samples"])
    t2.loc[len(t2)] = ["any (>=1)", int(t2["n_samples"].sum())]
    return t1, t1b, t2


def plot_kinship_density(kin, path_prefix, title, subtitle):
    """Kernel-density curve of kinship, with KING's band edges marked.

    A density curve rather than a histogram (user instruction): bar heights
    would be group counts. Matplotlib defaults, no grid.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.stats import gaussian_kde

    v = np.asarray(kin, dtype=float)
    grid = np.linspace(v.min(), v.max(), 512)
    dens = gaussian_kde(v)(grid)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(grid, dens)
    ax.fill_between(grid, dens, alpha=0.25)
    # label sides in plot order: the band to the LEFT of the edge first
    for x, lab in ((KIN_1ST, "0.177\n2nd | 1st"), (KIN_DUP, "0.354\n1st | dup")):
        if v.min() < x < v.max():
            ax.axvline(x, linestyle="--", linewidth=0.8, color="0.3")
            ax.annotate(lab, xy=(x, ax.get_ylim()[1] * 0.92), fontsize=8, ha="center", va="top")
    ax.set_xlabel("Kinship coefficient (AoU pc_relate)")
    ax.set_ylabel("Density")
    ax.set_title(f"{title}\n{subtitle}", fontsize=10)
    fig.tight_layout()

    written = {}
    for ext in ("svg", "png", "jpg"):
        p = f"{path_prefix}.{ext}"
        if ext == "jpg":
            fig.savefig(p, pil_kwargs={"quality": 85, "subsampling": 2})
        else:
            fig.savefig(p)
        written[ext] = os.path.getsize(p)
    plt.close(fig)
    return written


def keep_smallest_format(path_prefix, written):
    """Project rule: measure PNG/JPG/SVG, keep the smallest, delete the rest."""
    keep = min(written, key=written.get)
    for ext in written:
        if ext != keep:
            os.remove(f"{path_prefix}.{ext}")
    return keep


# ==========================================================================
# 6. ANCESTRY -- AoU's pre-computed genetic-ancestry call, for 03_ancestry_pca.py
# ==========================================================================
# AoU v9 ships a finished ancestry call for every participant in the release.
# We do not recompute it; we read it, verify that it still means what the
# release documentation says it means, and restrict the analysis cohort with it.
#
#   aux/ancestry/ancestry_preds.tsv   one row per research_id, 5 columns
#   aux/ancestry/training_pca.tsv     the HGDP+1000G reference samples
#   aux/ancestry/eigenvalues.txt      the 16 eigenvalues of the reference PCA
#
# Method (AoU, Nature 2024, PMC10937371, Methods):
#   "We trained a random forest classifier on a training set of the HGDP and
#    1000 Genomes samples variants on the autosome"
#   "We generated the first 16 principal components (PCs) of the training
#    sample genotypes (using the hwe_normalized_pca in Hail) at the
#    high-quality variant sites"
#   "The ancestry categories are based on the same labels used in gnomAD:
#    African (AFR); Latino/admixed American (AMR); East Asian (EAS); Middle
#    Eastern (MID); European (EUR), composed of Finnish (FIN) and Non-Finnish
#    European (NFE); Other (OTH), not belonging to one of the other ancestries
#    or is an admixture; South Asian (SAS)."
#
# So `pca_features` holds our samples PROJECTED onto the reference panel's 16
# axes -- not PCs computed inside our cohort. That is why they can be plotted
# against the labelled reference, and why they are the wrong thing to hand a
# within-EUR regression unchanged (see plot_pc_variance_within_eur).

ANC_DIR = os.path.join(V2, "data", "ancestry")
ANC_FIG_DIR = os.path.join(ANC_DIR, "figures")
ANC_GCS = "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/" "aux/ancestry/"
ANC_PREDS = ANC_GCS + "ancestry_preds.tsv"
ANC_TRAIN = ANC_GCS + "training_pca.tsv"
ANC_EIGEN = ANC_GCS + "eigenvalues.txt"

ANC_CACHE = os.path.join(ANC_DIR, "v2_03_ancestry_cohort.csv.gz")
ANC_REF_CACHE = os.path.join(ANC_DIR, "v2_03_reference_training_pca.tsv")
ANC_EIG_CACHE = os.path.join(ANC_DIR, "v2_03_reference_eigenvalues.txt")
EUR_KEEP = os.path.join(ANC_DIR, "step5_eur_keep.tsv")

# The order of the `probabilities` vector. Undocumented in the paper; derived
# from the release (for every row, ancestry_pred is the argmax under this
# order) and re-asserted on every extraction by extract_ancestry_preds.
ANCESTRY_CLASSES = ["afr", "amr", "eas", "eur", "mid", "sas"]
N_ANC_PC = 16
# max-prob strictly above this keeps the label in ancestry_pred_other; at or
# below it the sample becomes 'oth'. Also undocumented, also re-asserted.
OTH_THRESHOLD = 0.75


def git_check_ignored(path, must_be_ignored=True):
    """Assert that `path` is covered by .gitignore, and log the matching rule.

    Person-level files (keep lists, cached calls, unsuppressed figures) must
    never sit anywhere git can reach. Failing loudly here is the whole point:
    a silent miss is how participant data reaches a public repository.
    """
    r = subprocess.run(["git", "check-ignore", "-v", path], cwd=V2, capture_output=True, text=True)
    out = r.stdout.strip()
    log(f"  git check-ignore -v {path} -> {out or '[NOT IGNORED]'}")
    if must_be_ignored and r.returncode != 0:
        raise SystemExit(
            f"FATAL: {path} is NOT covered by .gitignore; refusing "
            "to write person-level data that could be pushed"
        )
    return out


from scripts.functions.cohort import gs_cat_lines


def read_keep_ids(path):
    """Read a keep list (plink2 '#IID' format or a one-column csv) as strings."""
    s = pd.read_csv(path, dtype=str).iloc[:, 0]
    return set(s.str.strip())


from scripts.functions.cohort import extract_ancestry_preds


def fetch_ancestry_reference():
    """Download the reference panel and its eigenvalues (inbound public data:
    HGDP+1000G samples, not AoU participants). Cached; -q so no progress bar
    lands in the git-tracked log."""
    for url, dest in [(ANC_TRAIN, ANC_REF_CACHE), (ANC_EIGEN, ANC_EIG_CACHE)]:
        if os.path.exists(dest):
            log(f"  CACHE HIT: {dest}")
            continue
        git_check_ignored(dest)
        subprocess.run(["gsutil", "-q", "-u", PROJECT_ID, "cp", url, dest], check=True)
        log(f"  downloaded: {dest} ({os.path.getsize(dest)/1e6:.2f} MB)")
    return ANC_REF_CACHE, ANC_EIG_CACHE


from scripts.functions.cohort import parse_hail_scores


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
    log(f"  {label}")
    for k, v in sorted(items.items(), key=lambda kv: -kv[1]):
        log(f"      {k:<24s} {v:>7,d}")
    if bucket:
        name = "other_combined(" + "+".join(sorted(bucket)) + ")"
        log(f"      {name:<24s} {sum(bucket.values()):>7,d}" f"   <- merged to clear n>={min_n}")
    log(f"      {'TOTAL':<24s} {sum(counter.values()):>7,d}")


def demotion_breakdown(df):
    """For the samples argmax-called EUR that the confidence threshold sends to
    'oth', which class ranks second. This is the evidence behind the locked
    cohort definition: it shows the exclusion is directional, not a random
    thinning, which is a limitation that has to reach the methods.
    """
    P = df[[f"p_{c}" for c in ANCESTRY_CLASSES]].to_numpy()
    second = [ANCESTRY_CLASSES[i] for i in np.argsort(-P, axis=1)[:, 1]]
    d = df.assign(second=second, maxp=P.max(1))
    dropped = d[(d.ancestry_pred == "eur") & (d.ancestry_pred_other == "oth")]
    kept = d[d.ancestry_pred_other == "eur"]
    return dropped, kept, collections.Counter(dropped["second"])


# ---- ancestry figures ----------------------------------------------------
# One colour per ancestry, fixed across every panel so the reader can carry the
# legend from figure to figure. matplotlib's default cycle, no restyling, and
# no grid anywhere (CLAUDE.md).
def _anc_colors():
    import matplotlib.pyplot as plt

    c = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    m = {
        "afr": c[0],
        "amr": c[1],
        "eas": c[2],
        "eur": c[3],
        "mid": c[4],
        "sas": c[5],
        "oth": "0.55",
    }
    m["mid+sas"] = m["mid"]
    return m


def anc_legend(ax, **kw):
    """Legend swatches at full opacity: the scatter alpha is for the cloud, not
    for the key, and at alpha=0.3 the swatches are unreadable."""
    lg = ax.legend(fontsize=8, frameon=False, **kw)
    for h in lg.legend_handles:
        h.set_alpha(1.0)
    return lg


def anc_scatter(ax, x, y, labels, order, size, alpha, colors, min_n=MIN_N, show_n=False):
    """Draw groups largest-first so small ones stay visible on top.

    Every label handed in must be drawn or explicitly excluded by the floor. A
    group silently missing from `order` would make the panel show fewer points
    than the caller believes, and two panels of the same cohort would stop
    being comparable with nothing on the figure saying so.
    """
    unlisted = set(labels) - set(order)
    if unlisted:
        raise SystemExit(
            f"FATAL: labels {sorted(unlisted)} are not in the draw "
            "order; merge them into a category or drop them on "
            "purpose, do not let them vanish"
        )
    drawn = 0
    for g in order:
        m = labels == g
        n = int(m.sum())
        if n == 0:
            continue
        if n < min_n:
            log(f"      [not drawn] {g}: below the n>={min_n} floor")
            continue
        drawn += n
        ax.scatter(
            x[m],
            y[m],
            s=size,
            alpha=alpha,
            color=colors[g],
            label=f"{g} (n={n:,d})" if show_n else g,
            linewidths=0,
        )
    log(f"      points drawn: {drawn:,d} of {len(labels):,d}")


def save_fig_smallest(fig, out_dir, stem, internal=False):
    """Write svg/png/jpg, measure all three, keep the smallest, delete the rest.

    Choose the output format by measured file size. Dense scatters are the
    documented case where a raster beats SVG, and these panels carry 13,455 +
    4,150 points. JPG, when it wins, gets quality=85/subsampling=2 because
    matplotlib's default of 75 is visibly worse.
    """
    import matplotlib.pyplot as plt

    prefix = os.path.join(out_dir, stem)
    written = {}
    for ext in ("svg", "png", "jpg"):
        p = f"{prefix}.{ext}"
        if ext == "jpg":
            fig.savefig(p, pil_kwargs={"quality": 85, "subsampling": 2})
        else:
            fig.savefig(p)
        written[ext] = os.path.getsize(p)
    plt.close(fig)
    best = keep_smallest_format(prefix, written)
    detail = "  ".join(f"{e}={written[e]:,d}B" for e in ("svg", "png", "jpg"))
    tag = ""
    if internal:
        tag = "  [INTERNAL, " + git_check_ignored(f"{prefix}.{best}") + "]"
    log(f"  {stem}: {detail}  -> kept .{best}{tag}")
    return written[best]


def plot_reference_vs_cohort(pcs, lab, order, refs, ref_lab, colors, show_n, min_n=MIN_N):
    """The figure that makes the axes interpretable: the labelled HGDP+1000G
    panel that DEFINES the PC space, next to our samples projected into it."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharex=True, sharey=True)
    ref_order = [g for g, _ in sorted(collections.Counter(ref_lab).items(), key=lambda kv: -kv[1])]
    anc_scatter(
        axes[0], refs[:, 0], refs[:, 1], ref_lab, ref_order, 8, 0.7, colors, min_n=0, show_n=show_n
    )
    axes[0].set_title(f"HGDP + 1000 Genomes reference (n={len(ref_lab):,d})")
    anc_scatter(
        axes[1], pcs[:, 0], pcs[:, 1], lab, order, 3, 0.3, colors, min_n=min_n, show_n=show_n
    )
    axes[1].set_title(f"Breast cancer cohort, projected (n={len(lab):,d})")
    for ax in axes:
        ax.set_xlabel("PC1")
        anc_legend(ax, markerscale=3)
    axes[0].set_ylabel("PC2")
    fig.tight_layout()
    return fig


def plot_cohort_pc_pairs(pcs, lab, order, colors, show_n, min_n=MIN_N):
    """PC1-PC2 and PC3-PC4 for the cohort: the standard supplementary panel."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for ax, (a, b) in zip(axes, [(0, 1), (2, 3)]):
        anc_scatter(
            ax, pcs[:, a], pcs[:, b], lab, order, 3, 0.3, colors, min_n=min_n, show_n=show_n
        )
        ax.set_xlabel(f"PC{a+1}")
        ax.set_ylabel(f"PC{b+1}")
        ax.set_title(f"PC{a+1} vs PC{b+1}")
    anc_legend(axes[0], markerscale=3)
    fig.tight_layout()
    return fig


def plot_two_label_rules(pcs, lab, lab_o, order, colors, show_n, suppress, min_n=MIN_N):
    """The same points twice: identical coordinates, identical classifier, only
    the decision rule differs. This is the figure the cohort definition rests
    on -- the 'oth' points sit at the junctions between arms, where admixed
    samples belong, not scattered at random."""
    import matplotlib.pyplot as plt
    import numpy as np

    fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharex=True, sharey=True)
    anc_scatter(
        axes[0], pcs[:, 0], pcs[:, 1], lab, order, 3, 0.3, colors, min_n=min_n, show_n=show_n
    )
    axes[0].set_title("ancestry_pred (argmax, no threshold)")
    if suppress:
        # mid is below the floor under this rule; merging it with the
        # next-smallest category both lifts the bucket over 20 and removes the
        # sas difference against the other panel. Drawn, not dropped, so both
        # panels still hold every sample.
        lab_draw = np.where(np.isin(lab_o, ["mid", "sas"]), "mid+sas", lab_o)
        order_o = ["eur", "afr", "amr", "oth", "eas", "mid+sas"]
    else:
        lab_draw = lab_o
        order_o = ["eur", "afr", "amr", "oth", "eas", "sas", "mid"]
    anc_scatter(
        axes[1], pcs[:, 0], pcs[:, 1], lab_draw, order_o, 3, 0.3, colors, min_n=min_n, show_n=show_n
    )
    axes[1].set_title(f"ancestry_pred_other (max-prob > {OTH_THRESHOLD}, " "else oth)")
    if suppress:
        axes[1].annotate(
            "mid and sas merged:\nmid alone is below the\n" f"n>={MIN_N-1} floor",
            xy=(0.02, 0.02),
            xycoords="axes fraction",
            fontsize=7,
            va="bottom",
        )
    for ax in axes:
        ax.set_xlabel("PC1")
        anc_legend(ax, markerscale=3)
    axes[0].set_ylabel("PC2")
    fig.tight_layout()
    return fig


def plot_reference_scree(eig):
    """Variance explained by each reference PC. Bar heights are eigenvalue
    fractions of the PUBLIC reference panel, not counts of participants, so
    this figure carries no disclosure implication."""
    import matplotlib.pyplot as plt
    import numpy as np

    ve = eig / eig.sum() * 100
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(np.arange(1, len(eig) + 1), ve)
    ax.set_xlabel("Principal component (reference panel)")
    ax.set_ylabel("Variance explained (%)")
    ax.set_xticks(np.arange(1, len(eig) + 1))
    fig.tight_layout()
    return fig, ve


def plot_confidence_density(maxp, lab, order, colors, min_n=MIN_N):
    """Kernel densities, not histograms: a bar height would be a group count.
    The flat curve is the diagnostic -- a class the classifier never calls
    confidently is a class this cohort cannot analyse."""
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.stats import gaussian_kde

    fig, ax = plt.subplots(figsize=(7, 4.5))
    grid = np.linspace(max(0.15, maxp.min() - 0.01), 1.0, 400)
    for g in order:
        m = lab == g
        if m.sum() < min_n or maxp[m].std() == 0:
            continue
        ax.plot(grid, gaussian_kde(maxp[m])(grid), color=colors[g], label=g)
    ax.axvline(OTH_THRESHOLD, color="k", linestyle="--", linewidth=1)
    ax.annotate(
        f"{OTH_THRESHOLD} threshold\n(below -> 'oth')",
        xy=(OTH_THRESHOLD, ax.get_ylim()[1] * 0.75),
        xytext=(0.36, ax.get_ylim()[1] * 0.78),
        fontsize=8,
        arrowprops=dict(arrowstyle="->", linewidth=0.8),
    )
    ax.set_xlabel("Random forest confidence, max over the 6 ancestry classes")
    ax.set_ylabel("Density")
    anc_legend(ax)
    fig.tight_layout()
    return fig


def plot_pc_variance_within_eur(pcs, masks):
    """How much of each projected axis survives the restriction to EUR.

    Ratio = SD within the subset / SD over the whole cohort. Near 0 means the
    axis separates continents and is nearly constant inside EUR, so entering it
    as a covariate adds a noise column; near 1 means it still carries structure.
    This is the quantitative reason PC1-PC10 is the wrong covariate set for a
    EUR-only analysis, and it is a variance argument, not a p-value one, so
    reading it does not compromise the pre-specified tests.
    """
    import matplotlib.pyplot as plt
    import numpy as np

    sd_all = pcs.std(axis=0, ddof=1)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    xs = np.arange(1, pcs.shape[1] + 1)
    w = 0.4
    ratios = {}
    for k, (name, mask) in enumerate(masks.items()):
        r = pcs[mask].std(axis=0, ddof=1) / sd_all
        ratios[name] = (int(mask.sum()), r)
        ax.bar(xs + (k - 0.5) * w, r, width=w, label=name)
    ax.axhline(1.0, color="k", linewidth=0.8)
    ax.set_xlabel("Projected principal component")
    ax.set_ylabel("SD within EUR / SD in the full cohort")
    ax.set_xticks(xs)
    anc_legend(ax)
    fig.tight_layout()
    return fig, ratios


def plot_within_eur_axes(pcs, m_eur, sub_lab, colors, pairs, show_n, min_n=MIN_N):
    """The surviving axes, EUR subset only, coloured by the confidence rule so
    the samples the threshold declines to label can be seen in place."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for ax, (a, b) in zip(axes, pairs):
        anc_scatter(
            ax,
            pcs[m_eur, a],
            pcs[m_eur, b],
            sub_lab,
            ["eur", "oth"],
            4,
            0.35,
            colors,
            min_n=min_n,
            show_n=show_n,
        )
        ax.set_xlabel(f"PC{a+1}")
        ax.set_ylabel(f"PC{b+1}")
        ax.set_title(f"PC{a+1} vs PC{b+1}, EUR subset")
    anc_legend(axes[0], markerscale=3, title="ancestry_pred_other", title_fontsize=8)
    fig.tight_layout()
    return fig

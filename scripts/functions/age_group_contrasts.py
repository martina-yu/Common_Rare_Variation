# Documentation: scripts/README.md.
from pathlib import Path
import json
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import norm
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
PCS = [f"PC{i}" for i in range(1, 11)]
TERMS = ["const", "PRS_z", "carrier", "interaction"] + PCS


def read(name):
    return pd.read_csv(ROOT / "data" / name)


def contrast(fit, vector):
    v = np.asarray(vector, dtype=float)
    estimate = float(v @ fit.params)
    se = float(np.sqrt(max(0, v @ fit.cov_params() @ v)))
    return {
        "estimate_years": estimate,
        "se_HC1": se,
        "ci_low": estimate - norm.ppf(0.975) * se,
        "ci_high": estimate + norm.ppf(0.975) * se,
    }


def design(p, c, pc):
    return np.array([1, p, c, p * c, *pc], dtype=float)


def analyze(df, panel, label, reference, out):
    mu, sd, cutoff, p_low, p_high, pc_ref = reference
    df = df.copy()
    df["carrier"] = (df[f"burden_rare_main_ptv_{panel}gene"] > 0).astype(int)
    df["high"] = (df.PRS311 >= cutoff).astype(int)
    df["PRS_z"] = (df.PRS311 - mu) / sd
    df["interaction"] = df.PRS_z * df.carrier
    counts = (
        df.groupby(["carrier", "high"])
        .size()
        .reindex(pd.MultiIndex.from_product([[0, 1], [0, 1]]), fill_value=0)
    )
    audit = {"analysis": label, "four_group_release": bool((counts >= 20).all())}
    # Suppress the entire table if a cell is small; do not emit recoverable margins.
    if (df.groupby("carrier").size().reindex([0, 1], fill_value=0) < 20).any():
        audit["continuous_model"] = "suppressed"
        return audit
    X = sm.add_constant(df[["PRS_z", "carrier", "interaction"] + PCS], has_constant="add")
    fit = sm.OLS(df.age_at_diagnosis.astype(float), X).fit(cov_type="HC1")
    rows = []
    for term in ["PRS_z", "carrier", "interaction"]:
        v = np.zeros(len(TERMS))
        v[TERMS.index(term)] = 1
        rows.append({"contrast": term, **contrast(fit, v)})
    v = np.zeros(len(TERMS))
    v[1] = v[3] = 1
    rows.append({"contrast": "PRS_per_SD_carriers", **contrast(fit, v)})
    pd.DataFrame(rows).to_csv(out / f"{label}_continuous.csv", index=False)
    if audit["four_group_release"]:
        vectors, summaries = {}, []
        for c in [0, 1]:
            for h in [0, 1]:
                g = df[(df.carrier == c) & (df.high == h)].age_at_diagnosis
                v = design(p_high if h else p_low, c, pc_ref)
                vectors[c, h] = v
                summaries.append(
                    {
                        "carrier": c,
                        "high_PRS": h,
                        "n": len(g),
                        "observed_mean": g.mean(),
                        "observed_sd": g.std(),
                        "median": g.median(),
                        "q25": g.quantile(0.25),
                        "q75": g.quantile(0.75),
                        **contrast(fit, v),
                    }
                )
        pd.DataFrame(summaries).to_csv(out / f"{label}_four_groups.csv", index=False)
        comparisons = {
            "high_vs_nonhigh_noncarriers": vectors[0, 1] - vectors[0, 0],
            "high_vs_nonhigh_carriers": vectors[1, 1] - vectors[1, 0],
            "difference_in_differences": vectors[1, 1]
            - vectors[1, 0]
            - vectors[0, 1]
            + vectors[0, 0],
            "double_high_vs_noncarrier_nonhigh": vectors[1, 1] - vectors[0, 0],
        }
        pd.DataFrame([{"contrast": k, **contrast(fit, v)} for k, v in comparisons.items()]).to_csv(
            out / f"{label}_contrasts.csv", index=False
        )
    # Only aggregate model curves; common support avoids carrier-tail extrapolation.
    lo = max(df.loc[df.carrier == c, "PRS_z"].quantile(0.01) for c in [0, 1])
    hi = min(df.loc[df.carrier == c, "PRS_z"].quantile(0.99) for c in [0, 1])
    if lo < hi:
        fig, ax = plt.subplots()
        for c in [0, 1]:
            xs = np.linspace(lo, hi, 100)
            predictions = pd.DataFrame([contrast(fit, design(p, c, pc_ref)) for p in xs])
            (line,) = ax.plot(xs, predictions.estimate_years, label=f"Carrier={c}")
            ax.fill_between(
                xs, predictions.ci_low, predictions.ci_high, alpha=0.2, color=line.get_color()
            )
        ax.set(
            xlabel="PRS311 (SD in primary cases)", ylabel="Adjusted recorded diagnosis age (years)"
        )
        ax.legend()
        fig.tight_layout()
        fig.savefig(out / f"{label}_curves.svg")
        plt.close(fig)
    event_counts = (
        df.groupby(["carrier", "high", df.age_at_diagnosis < 45])
        .size()
        .reindex(pd.MultiIndex.from_product([[0, 1], [0, 1], [False, True]]), fill_value=0)
    )
    audit["case_only_early_cells_release_eligible"] = bool((event_counts >= 20).all())
    audit["population_early_risk"] = "not_estimated_observation_window_requires_review"
    return audit

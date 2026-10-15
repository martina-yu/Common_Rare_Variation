"""Generate the existing score, discrimination, diagnosis-age and sensitivity figures.

Figure data, bins, model specifications, labels and output paths are preserved.
Discrimination and joint-distribution sections retain their explicitly visible model calculations.
"""

from __future__ import annotations

from pathlib import Path
import sys

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.functions.paths import PROJECT_ROOT
from scripts.functions.cli import parse_tasks

TASKS = (
    'age-figures',
    'score-distributions',
    'discrimination',
    'joint-distribution',
    'score-deciles',
    'burden-histogram',
    'age-boxplots',
    'age-carrier-frequency',
    'sensitivity-figure',
)
if __name__ == '__main__':
    selected_tasks = parse_tasks(__doc__, TASKS, ['age-figures'])

import hashlib

import json


import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt

import numpy as np

import pandas as pd

import statsmodels.api as sm

from threadpoolctl import threadpool_limits

AGE_FIGURES_ROOT = PROJECT_ROOT

AGE_FIGURES_BASE = AGE_FIGURES_ROOT / 'data/29_rowlands_extension'

AGE_FIGURES_OUT = AGE_FIGURES_BASE / 'section12_figures'

FIG = AGE_FIGURES_ROOT / 'results/figures'

TABLE = AGE_FIGURES_ROOT / 'results/table'

Z = 1.959963984540054

PC = [f'PC{i}' for i in range(1, 11)]

V = ['A', 'A_prime', 'B', 'C', 'C_without_STK11']

LABEL = [
    '9 genes: PTV\nOriginal weights',
    '9 genes: PTV\nRowlands weights',
    '9 genes: PTV + qualifying P/LP missense\nRowlands weights',
    '13 genes: PTV + qualifying P/LP missense\nRowlands weights',
    '12 genes (without STK11): PTV + qualifying P/LP missense\nRowlands weights',
]

matplotlib.use('Agg')

from matplotlib.patches import Patch

from scipy.stats import gaussian_kde

MS_4_ROOT = PROJECT_ROOT

SCORE_DISTRIBUTION_FIGURE_DIR = MS_4_ROOT / 'manuscript/figures/internal_score_distributions'

MS_4_BLUE = '#2463B4'

MS_4_RED = '#D54848'

from concurrent.futures import ThreadPoolExecutor


matplotlib.use('Agg')

from matplotlib.lines import Line2D

from sklearn.metrics import roc_auc_score, roc_curve

from sklearn.model_selection import RepeatedStratifiedKFold

MS_5_ROOT = PROJECT_ROOT

sys.path.insert(0, str(MS_5_ROOT))

from scripts.functions.stats import auc_ci, delong_test

MS_5_BASE = MS_5_ROOT / 'data/29_rowlands_extension'

MS_5_EXISTING = MS_5_BASE / 'section11_models'

DISCRIMINATION_CACHE_DIR = MS_5_ROOT / 'data/manuscript_refresh_20261006/four_model_CV'

MS_5_FIG = MS_5_ROOT / 'manuscript/figures'

MS_5_DATA = MS_5_ROOT / 'manuscript/datas/case_control'

DISCRIMINATION_VERSIONS = [
    ('C', '13 genes: main'),
    ('C_without_STK11', '12 genes: without STK11'),
    ('C_PTV_only', '13 genes: PTV only'),
    ('C_TP53_NF1_strict_VAF', '13 genes: strict TP53/NF1 VAF'),
    ('C_ClinVar_min1star', '13 genes: ClinVar >=1 star'),
    ('C_without_CDH1', '12 genes: without CDH1'),
]

DISCRIMINATION_MODELS = ['Covariates only', 'PRS + covariates', 'Burden + covariates', 'Joint']

MS_5_COLORS = ['#777777', '#E58A27', '#2768B4', '#329461']

matplotlib.use('Agg')

from matplotlib.colors import TwoSlopeNorm

MS_6_ROOT = PROJECT_ROOT

MS_6_FIG = MS_6_ROOT / 'manuscript/figures'

MS_6_DATA = MS_6_ROOT / 'manuscript/datas/case_control'

JOINT_DISTRIBUTION_CACHE_DIR = MS_6_ROOT / 'data/manuscript_refresh_20261006/joint_distribution'

matplotlib.use('Agg')

MS_7_ROOT = PROJECT_ROOT

MS_7_FIG = MS_7_ROOT / 'manuscript/figures'

SCORE_DECILE_CACHE_DIR = MS_7_ROOT / 'data/manuscript_refresh_20261006/score_decile_boxplots'

matplotlib.use('Agg')

MS_8_ROOT = PROJECT_ROOT

MS_8_FIG = MS_8_ROOT / 'manuscript/figures'

BURDEN_FIGURE_DATA_DIR = MS_8_ROOT / 'manuscript/datas/burden'

matplotlib.use('Agg')

MS_9_ROOT = PROJECT_ROOT

MS_9_BASE = MS_9_ROOT / 'data/29_rowlands_extension'

import shutil

matplotlib.use('Agg')

from statsmodels.stats.proportion import proportion_confint

MS_10_ROOT = PROJECT_ROOT

MS_10_BASE = MS_10_ROOT / 'data/29_rowlands_extension'

MS_10_FIG = MS_10_ROOT / 'manuscript/figures'

matplotlib.use('Agg')

MS_12_ROOT = PROJECT_ROOT

SENSITIVITY_FIGURE_VERSIONS = [
    ('C', '13 genes: primary'),
    ('C_without_STK11', '12 genes: exclude STK11'),
    ('B', '9 genes: PTV + qualifying missense'),
    ('C_PTV_only', '13 genes: PTV only'),
    ('C_TP53_NF1_strict_VAF', '13 genes: stricter TP53/NF1 VAF'),
    ('C_ClinVar_min1star', '13 genes: ClinVar ≥1 star'),
    ('C_without_CDH1', '12 genes: exclude CDH1'),
    ('A', '9 genes: original PTV baseline'),
]

from scripts.functions.io import sha256 as sha


from scripts.functions.plotting import compressed_kde


def save_score_distribution_figure(fig, stem):
    fig.savefig(SCORE_DISTRIBUTION_FIGURE_DIR / f'{stem}.svg', bbox_inches='tight')
    png = SCORE_DISTRIBUTION_FIGURE_DIR / f'{stem}.png'
    jpg = SCORE_DISTRIBUTION_FIGURE_DIR / f'{stem}.jpg'
    fig.savefig(png, dpi=170, bbox_inches='tight')
    fig.savefig(jpg, dpi=170, bbox_inches='tight', pil_kwargs={'quality': 85, 'subsampling': 2})
    smallest = '.png' if png.stat().st_size <= jpg.stat().st_size else '.jpg'
    if smallest == '.png':
        jpg.unlink()
    else:
        png.unlink()
    plt.close(fig)
    return smallest


from scripts.functions.io import sha256 as MS_5_sha


def discrimination_design_matrix(cov, burden=None):
    cols = [np.ones(len(cov))]
    if burden is not None:
        cols.append(burden)
    return np.column_stack(cols + [cov])


def fit_discrimination_model(y, train_X, test_X):
    fit = sm.GLM(y, train_X, family=sm.families.Binomial()).fit()
    assert fit.converged and np.isfinite(fit.params).all()
    return fit.predict(test_X)


def save_discrimination_figure(fig, name):
    fig.savefig(MS_5_FIG / f'{name}.svg', bbox_inches='tight')
    png, jpg = (MS_5_FIG / f'{name}.png', MS_5_FIG / f'{name}.jpg')
    fig.savefig(png, dpi=180, bbox_inches='tight')
    fig.savefig(jpg, dpi=180, bbox_inches='tight', pil_kwargs={'quality': 85, 'subsampling': 2})
    smallest = png if png.stat().st_size <= jpg.stat().st_size else jpg
    (jpg if smallest == png else png).unlink()
    plt.close(fig)
    return str(smallest.relative_to(MS_5_ROOT))


if __name__ == "__main__":
    for selected_task in selected_tasks:

        # ============================================================================
        # 1. Plot age and auc results
        # ============================================================================
        if selected_task == 'age-figures':
            AGE_FIGURES_OUT.mkdir(parents=True, exist_ok=True)

            import subprocess

            assert (
                subprocess.run(
                    ['git', 'check-ignore', '-q', str(AGE_FIGURES_OUT)], cwd=AGE_FIGURES_ROOT
                ).returncode
                == 0
            )

            inputs = [
                AGE_FIGURES_BASE / 'section9/frozen_model_inputs.csv.gz',
                AGE_FIGURES_BASE / 'section8/comparable_scores_case_control.csv.gz',
                AGE_FIGURES_BASE / 'section11_models/model_coefficients_internal.csv',
                AGE_FIGURES_BASE / 'section11_models/CV_delta_summary_internal.csv',
                AGE_FIGURES_BASE / 'section11_models/carrier_transition_flags_internal.csv.gz',
            ]

            hashes = {str(p.relative_to(AGE_FIGURES_ROOT)): sha(p) for p in inputs}

            all_df = pd.read_csv(inputs[0]).set_index('person_id')

            ids = pd.read_csv(AGE_FIGURES_BASE / 'diagnosis_age_ids.csv.gz').person_id

            df = all_df.loc[ids].copy()

            scores = pd.read_csv(inputs[1]).set_index('sample_id')

            carrier = scores.loc[df.index, 'carried_gene_count_C'].gt(0).astype(float)

            assert len(df) == 7996 and df.is_case.eq(1).all() and (int(carrier.sum()) == 543)

            assert int(carrier.sum()) >= 20 and int((1 - carrier).sum()) >= 20

            prs = df.PRS311

            mu = prs.mean()

            sd = prs.std(ddof=1)

            p = (prs - mu) / sd

            X = pd.DataFrame(
                {
                    'const': 1.0,
                    'PRS311_z': p,
                    'C_carrier': carrier,
                    'PRS311_z_x_C_carrier': p * carrier,
                },
                index=df.index,
            ).join(df[PC])

            with threadpool_limits(limits=1):
                fit = sm.OLS(df.age_at_diagnosis, X).fit(cov_type='HC1')

            cov = fit.cov_params().to_numpy()

            beta = fit.params.to_numpy()

            assert np.isfinite(beta).all() and np.isfinite(cov).all()

            coeff = pd.DataFrame(
                {
                    'term': X.columns,
                    'beta': beta,
                    'robust_se_HC1': fit.bse.to_numpy(),
                    'CI_low': fit.conf_int()[0].to_numpy(),
                    'CI_high': fit.conf_int()[1].to_numpy(),
                    'p_value': fit.pvalues.to_numpy(),
                }
            )

            coeff.to_csv(
                AGE_FIGURES_OUT / 'age_binary_carrier_interaction_all_coefficients_internal.csv',
                index=False,
            )

            contrasts = []

            for name, terms in [
                ('Non-carrier PRS slope', {'PRS311_z': 1}),
                ('Carrier PRS slope', {'PRS311_z': 1, 'PRS311_z_x_C_carrier': 1}),
                ('Carrier minus non-carrier slope (interaction)', {'PRS311_z_x_C_carrier': 1}),
            ]:
                c = np.zeros(len(beta))
                for term, value in terms.items():
                    c[X.columns.get_loc(term)] = value
                estimate = c @ beta
                se = np.sqrt(c @ cov @ c)
                contrasts.append(
                    dict(
                        contrast=name,
                        beta_years_per_PRS_SD=estimate,
                        robust_se_HC1=se,
                        CI_low=estimate - Z * se,
                        CI_high=estimate + Z * se,
                        beta_years_per_raw_PRS_unit=estimate / sd,
                        raw_CI_low=(estimate - Z * se) / sd,
                        raw_CI_high=(estimate + Z * se) / sd,
                    )
                )

            ct = pd.DataFrame(contrasts)

            ct.to_csv(AGE_FIGURES_OUT / 'age_binary_carrier_PRS_slopes_internal.csv', index=False)

            ct.drop(columns=['beta_years_per_raw_PRS_unit', 'raw_CI_low', 'raw_CI_high']).to_csv(
                TABLE / '48_age_binary_carrier_PRS_slopes.csv', index=False
            )

            ranges = [p.loc[carrier.eq(g)].quantile([0.025, 0.975]).to_numpy() for g in [0, 1]]

            lo = max((r[0] for r in ranges))

            hi = min((r[1] for r in ranges))

            assert lo < hi

            grid = np.linspace(lo, hi, 160)

            pred = []

            for group in [0, 1]:
                xp = pd.DataFrame(
                    {
                        'const': np.ones(len(grid)),
                        'PRS311_z': grid,
                        'C_carrier': float(group),
                        'PRS311_z_x_C_carrier': grid * group,
                    }
                )
                for pc in PC:
                    xp[pc] = df[pc].mean()
                xp = xp[X.columns].to_numpy()
                expected = xp @ beta
                se = np.sqrt(np.maximum(np.einsum('ij,jk,ik->i', xp, cov, xp), 0))
                for j in range(len(grid)):
                    pred.append(
                        dict(
                            C_carrier=group,
                            PRS311_case_z=grid[j],
                            predicted_age=expected[j],
                            CI_low=expected[j] - Z * se[j],
                            CI_high=expected[j] + Z * se[j],
                        )
                    )

            pd.DataFrame(pred).to_csv(
                AGE_FIGURES_OUT / 'adjusted_age_prediction_grid_internal.csv', index=False
            )

            for g in [0, 1]:
                q = pd.DataFrame(pred).loc[lambda d: d.C_carrier.eq(g)]
                slope = (q.predicted_age.iloc[-1] - q.predicted_age.iloc[0]) / (grid[-1] - grid[0])
                assert np.isclose(slope, ct.beta_years_per_PRS_SD.iloc[g])
                assert (q.CI_low <= q.predicted_age).all() and (q.CI_high >= q.predicted_age).all()

            plt.style.use('default')

            colors = ['tab:blue', 'tab:orange']

            fig, axes = plt.subplots(
                1, 2, figsize=(16.8, 5.8), gridspec_kw={'width_ratios': [1, 1.35]}
            )

            age = pd.read_csv(inputs[2])

            age = (
                age.loc[age.cohort.eq('diagnosis_age') & age.term.eq('burden_z')]
                .set_index('version')
                .loc[V]
            )

            yy = np.arange(len(V))

            axes[0].errorbar(
                age.beta,
                yy,
                xerr=np.vstack([age.beta - age.CI_beta_low, age.CI_beta_high - age.beta]),
                fmt='o',
                color='tab:blue',
                capsize=4,
            )

            axes[0].axvline(0, color='0.55', linestyle='--', linewidth=1)

            axes[0].set_yticks(yy, LABEL)

            axes[0].invert_yaxis()

            axes[0].tick_params(axis='y', labelsize=9)

            axes[0].set_xlabel('Difference in recorded diagnosis age\n(years per SD burden)')

            axes[0].set_title('Burden definitions and weights', loc='left', fontweight='bold')

            for g in [0, 1]:
                q = pd.DataFrame(pred).loc[lambda d: d.C_carrier.eq(g)]
                axes[1].plot(
                    q.PRS311_case_z,
                    q.predicted_age,
                    color=colors[g],
                    label='Non-carrier' if g == 0 else '13-gene carrier',
                )
                axes[1].fill_between(
                    q.PRS311_case_z, q.CI_low, q.CI_high, color=colors[g], alpha=0.18
                )

            axes[1].set_xlabel('Continuous PRS311 (SD within fixed cases)')

            axes[1].set_ylabel('PC-adjusted predicted recorded diagnosis age (years)')

            axes[1].set_title(
                'PRS × 13-gene carrier status (exploratory)', loc='left', fontweight='bold'
            )

            annotation = []

            for r in ct.itertuples():
                name = {
                    'Non-carrier PRS slope': 'Non-carrier slope',
                    'Carrier PRS slope': 'Carrier slope',
                    'Carrier minus non-carrier slope (interaction)': 'Interaction (slope difference)',
                }[r.contrast]
                annotation.append(
                    f'{name}: {r.beta_years_per_PRS_SD:.2f} [{r.CI_low:.2f}, {r.CI_high:.2f}]'
                )

            axes[1].text(
                0.02,
                0.98,
                'Years per PRS SD; 95% HC1 CI\n' + '\n'.join(annotation),
                transform=axes[1].transAxes,
                va='top',
                fontsize=9,
                bbox={'facecolor': 'white', 'alpha': 0.88, 'edgecolor': '0.85'},
            )

            prediction_frame = pd.DataFrame(pred)

            axes[1].set_ylim(
                prediction_frame.CI_low.min() - 1.1, prediction_frame.CI_high.max() + 2.5
            )

            axes[1].legend(loc='lower left', frameon=False)

            for ax in axes:
                ax.grid(False)

            fig.text(
                0.5,
                0.01,
                'Fixed breast cancer cases (N = 7,996). Burden SD differs by definition. Right panel: additional binary-carrier interaction model.',
                ha='center',
                fontsize=9,
            )

            fig.tight_layout(rect=[0, 0.05, 1, 1])

            fig.savefig(FIG / '48_main_recorded_diagnosis_age.svg')

            plt.close(fig)

            f = pd.read_csv(inputs[4]).set_index('sample_id')

            shared = int((f.A_carrier & f.C_carrier).sum())

            lost = int((f.A_carrier & ~f.C_carrier).sum())

            names = [
                'Original nine-gene PTV',
                'PTV eligibility / call handling changes',
                'Added qualifying missense (nine genes)',
                'Added CDH1, NF1, PTEN, STK11',
            ]

            count = [shared] + [
                int(f[k].sum())
                for k in ['new_C_current_nine_PTV', 'new_C_nine_missense', 'new_C_added_four_only']
            ]

            assert lost == 0, 'A-only carriers need a separate suppression-safe diagram'

            assert (
                all((n >= 20 for n in count))
                and sum(count) == int(f.C_carrier.sum())
                and (sum(count[1:]) == 487)
            )

            fig, ax = plt.subplots(figsize=(11, 3.4))

            palette = ['0.7', 'tab:purple', 'tab:blue', 'tab:orange']

            left = 0

            for name, n, col in zip(names, count, palette):
                ax.barh(0, n, left=left, color=col, label=f'{name}: {n:,}', height=0.5)
                if n >= 200:
                    ax.text(
                        left + n / 2,
                        0,
                        f'{n:,}',
                        ha='center',
                        va='center',
                        color='black' if col == '0.7' else 'white',
                    )
                left += n

            ax.set_yticks([0], ['13-gene carriers'])

            ax.set_xlabel('Participants, mutually exclusive attribution')

            ax.set_title(
                'Carrier expansion: 487 newly identified carriers', loc='left', fontweight='bold'
            )

            ax.grid(False)

            ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.26), ncol=2, frameon=False)

            fig.tight_layout()

            fig.savefig(FIG / '48_supp_carrier_increment.svg', bbox_inches='tight')

            plt.close(fig)

            dd = pd.read_csv(inputs[3])

            dd = (
                dd.loc[dd.comparison.eq('joint_minus_PRS311_baseline')]
                .set_index('version')
                .loc[V[:4]]
            )

            fig, ax = plt.subplots(figsize=(10.4, 4.8))

            yy = np.arange(4)

            ax.errorbar(
                dd.delta_AUC_repeat1,
                yy,
                xerr=np.vstack(
                    [
                        dd.delta_AUC_repeat1 - dd.repeat1_CI_low,
                        dd.repeat1_CI_high - dd.delta_AUC_repeat1,
                    ]
                ),
                fmt='o',
                capsize=4,
                color='tab:blue',
            )

            ax.axvline(0, color='0.55', linestyle='--', linewidth=1)

            ax.set_yticks(yy, LABEL[:4])

            ax.invert_yaxis()

            ax.grid(False)

            ax.tick_params(axis='y', labelsize=9)

            ax.set_xlabel('OOF ΔAUC relative to covariates + PRS311')

            ax.set_title('Fixed cohort and matched CV folds', loc='left')

            fig.text(
                0.5,
                0.01,
                'Points and paired DeLong 95% CI: repeat 1 of fixed 10 × 5 CV.',
                ha='center',
                fontsize=9,
            )

            fig.tight_layout(rect=[0, 0.05, 1, 1])

            fig.savefig(FIG / '48_supp_delta_AUC.svg')

            plt.close(fig)

            dist = pd.read_csv(
                AGE_FIGURES_BASE
                / 'section10/internal_delivery/03_A_to_C_carrier_and_model_comparison.csv'
            )

            dist = dist.loc[
                dist.version.isin(V),
                [
                    'cohort',
                    'version',
                    'n',
                    'carrier_fraction',
                    'zero_fraction',
                    'raw_mean',
                    'SD_population_ddof0',
                    'SD_sample_ddof1',
                    'raw_p25',
                    'raw_p50',
                    'raw_p75',
                    'raw_p95',
                    'raw_p99',
                ],
            ]

            dist.to_csv(
                AGE_FIGURES_OUT / 'burden_distribution_carrier_supplement_exact_internal.csv',
                index=False,
            )

            pub = dist.copy()

            for ix, r in pub.iterrows():
                precision = 1 if r.cohort == 'case_control' else 0
                pub.loc[ix, 'carrier_fraction'] = round(r.carrier_fraction * 100, precision)
                pub.loc[ix, 'zero_fraction'] = round(r.zero_fraction * 100, precision)

            pub = pub.rename(
                columns={
                    'carrier_fraction': 'carrier_percent_coarsened',
                    'zero_fraction': 'zero_percent_coarsened',
                }
            )

            definitions = dict(zip(V, [x.replace('\n', '; ') for x in LABEL]))

            pub['burden_definition'] = pub.version.map(definitions)

            pub = pub.drop(columns='version')

            numeric = [c for c in pub.columns if c.startswith('raw_') or c.startswith('SD_')]

            pub[numeric] = pub[numeric].round(2)

            pub.to_csv(
                TABLE / '48_burden_distribution_carrier_supplement_coarsened.csv', index=False
            )

            comparison_path = (
                AGE_FIGURES_BASE
                / 'section10/internal_delivery/03_A_to_C_carrier_and_model_comparison.csv'
            )

            comparison = pd.read_csv(comparison_path)

            target = comparison.cohort.eq('diagnosis_age') & comparison.version.eq('C')

            assert target.sum() == 1

            for j, prefix in enumerate(
                [
                    'age_binary_model_PRS_noncarrier',
                    'age_binary_model_PRS_carrier',
                    'age_binary_model_interaction',
                ]
            ):
                for field, col in [
                    ('beta', 'beta_years_per_PRS_SD'),
                    ('CI_low', 'CI_low'),
                    ('CI_high', 'CI_high'),
                ]:
                    comparison.loc[target, f'{prefix}_{field}'] = ct.iloc[j][col]

            comparison.to_csv(comparison_path, index=False)

            delivery_manifest_path = AGE_FIGURES_BASE / 'section10/delivery_manifest_internal.json'

            delivery_manifest = json.loads(delivery_manifest_path.read_text())

            delivery_manifest.setdefault(
                'post_main_model_delivery_sha256', delivery_manifest['delivery'].copy()
            )

            delivery_manifest['delivery'] = {
                p.name: sha(p) for p in comparison_path.parent.iterdir() if p.is_file()
            }

            delivery_manifest['additional_exploratory_age_binary_carrier_model'] = (
                'data/29_rowlands_extension/section12_figures/figure_and_model_manifest_internal.json'
            )

            delivery_manifest_path.write_text(json.dumps(delivery_manifest, indent=2) + '\n')

            ctrow = ct.iloc[2]

            caption = f"# Diagnostic age priority chart set: legends and explanations\n\nPresentation order: Supplementary carrier delta plot → Age at diagnosis two-panel main plot → Supplementary ΔAUC plot. The text in the figure is in English to facilitate the export of the paper.\n\n**Program definition**: The original nine genes are ATM, BARD1, BRCA1, BRCA2, CHEK2, PALB2, RAD51C, RAD51D, and TP53. The 13-gene program adds CDH1, NF1, PTEN, and STK11 on this basis; the 12-gene sensitivity program only removes STK11. PTV refers to protein truncating variants, and qualified missense refers to missense variants that pass this program's germline P/LP, ClinVar at least 2 stars, low penetrance/risk allele/conflict elimination, related diseases and inheritance patterns, frequency of European ancestry, and calling QC. The Rowlands weight comes from the gene OR of PTV+non-PTV-PV, which is the approximate external weight of the PTV+qualified P/LP missense score in this study. The mutation definitions of the two are not exactly the same.\n\n**Carrier incremental chart**: 3,209 carriers in the 13-gene PTV + qualified P/LP missense program, which are mutually exclusive attributions based on those identified by the original nine-gene PTV (2,722 people), those newly added to the current PTV inclusion/calling process (148 people), those newly added with the nine-gene qualified missense (239 people), and those identified after the addition of CDH1/NF1/PTEN/STK11 (100 people). A total of 487 new people were added, and the difference in PTV processing cannot be counted as a missense increase. All the original nine-gene PTV carriers are retained in the 13-gene plan, and each display unit has at least 20 people; no small gene stratification is drawn.\n\n**Left side of main picture of diagnostic age: burden definition and weight comparison**. Among the fixed 7,996 breast cancer cases, the following comparisons were made sequentially: the original nine-gene PTV + original weight; the same original nine-gene PTV carrier + Rowlands weight; the current nine-gene PTV + qualified P/LP missense + Rowlands weight; the 13-gene PTV + qualified P/LP missense + Rowlands weight; the 12-gene scheme without STK11 + Rowlands weight. The OLS models recording age at diagnosis also included PRS311, adjusted PC1–10, and used HC1 robust standard errors. The points are the burden coefficients, the lines are the 95% CIs, and the zero value has a reference line. Each burden is standardized by its own SD within the case cohort, and the SD of each scheme is different; the original scale effect cannot be judged based solely on the change in the coefficient of each SD. The precise distribution, SD, and carrier ratio are left in the internal supplementary table, and a jointly coarsened version is provided externally.\n\n**Right side of main picture of age at diagnosis: interaction between PRS311 and 13 gene carrier status**. The new expanded model is `record diagnosis age ~ PRS311 standard score + 13 gene carrier index + PRS311 standard score × 13 gene carrier index + PC1–10`, using the same cases and HC1. Carrier is defined as a case carrying at least one 13-gene scheme-eligible PTV or P/LP missense. This binary indicator replaces the continuous burden in this extended model; the original age main model maintains the no-interaction setting. The horizontal axis is the standard score of continuous PRS311 in fixed cases. The two lines set the PC to the mean of the same fixed case cohort; under the linear model, it is equivalent to marginal standardization on the same PC distribution. Shaded are point-wise 95% HC1 confidence intervals predicting the mean recorded age at diagnosis, not individual prediction intervals or simultaneous confidence bands. The horizontal axis uses the intersection of the central 95% PRS ranges of each of the two groups. The slopes and slope differences of the two groups of PRS on the figure are calculated using the complete robust covariance matrix; the vertical distance between the two lines mainly reflects the carrier main effect and cannot be used to claim synergy.\n\nThe within-case binary carrier interaction (slope difference) is{ctrow.beta_years_per_PRS_SD:.3f}Year/PRS SD (95% CI{ctrow.CI_low:.3f}to{ctrow.CI_high:.3f}). This is a new age-related interaction estimate. It belongs to a different model from the continuous PRS × burden interaction in the previous case status logistic and cannot be mixed.\n\n**ΔAUC Supplementary Figure**: fixed 141,204 people, same PRS311 and covariates, original 10×5 discount. The original nine-gene PTV+original weight, the same PTV carrier+Rowlands weight, the current nine-gene PTV+qualified missense+Rowlands weight, the 13-gene PTV+qualified missense+Rowlands weight, and the joint model OOF ΔAUC relative to covariates+PRS311 baseline are shown in sequence. The joint model adds burden and the existing PRS×burden interaction. Both points and CI come from the first round of OOF predictions, and CI is paired DeLong; the average and range of 10 repetitions are listed separately in the existing supplementary table. The CI condition is for this round of fitting predictions and does not include the uncertainty of complete model retraining; the 10 repetitions are not regarded as independent samples. The difference between the former nine-gene combined scheme and the original nine-gene PTV weighting scheme includes both PTV processing and missense changes; when isolating the missense contribution, it should be compared with the nine-gene PTV-only bridging scheme using the same current QC.\n\nAll age results are documented age-at-diagnosis associations among breast cancer cases only and are not interpreted as time-to-onset risks in the general population. The variant biology review is still pending and the panel is an extended analysis under the tentative definition.\n\nPublicly available data only include overall model/linear comparisons and coarsened supplementary statistics, without participant-by-participant data points. Subtotals, percentages, and small cells that can be recovered from differences must be checked simultaneously. Basis: [All of Us Disclosure Statement](https://support.researchallofus.org/hc/en-us/articles/360043016291-How-to-comply-with-the-All-of-Us-Data-and-Statistics-Dissemination-Policy)."

            (AGE_FIGURES_ROOT / 'results/report/48_age_priority_figure_captions.md').write_text(
                caption
            )

            for path, h in hashes.items():
                assert sha(AGE_FIGURES_ROOT / path) == h, path

            manifest = {
                'status': 'age_priority_figures_and_exploratory_age_carrier_interaction_complete',
                'fixed_case_n': 7996,
                'fixed_cases_C_carriers': 543,
                'fixed_cases_noncarriers': 7453,
                'main_age_model_unchanged': True,
                'additional_age_model': 'PRS311_case_z + C_binary_carrier + interaction + PC1-10; OLS HC1',
                'slope_difference_CI_uses_full_robust_covariance': True,
                'line_slopes_numerically_verified': True,
                'AUC_points_and_CI_same_repeat': True,
                'carrier_attribution_mutually_exclusive': True,
                'carrier_cells_all_ge20_or_zero': True,
                'public_exact_sparse_score_differences_not_reported': True,
                'inputs': hashes,
                'PRS311_age_scaling': {'mean': mu, 'SD_ddof1': sd},
                'x_range_case_PRS_z': [lo, hi],
                'final_variant_definition_frozen': False,
            }

            (AGE_FIGURES_OUT / 'figure_and_model_manifest_internal.json').write_text(
                json.dumps(manifest, indent=2) + '\n'
            )

            (
                AGE_FIGURES_ROOT / 'results/report/48_age_priority_figures_validation.json'
            ).write_text(
                json.dumps(
                    {
                        k: v
                        for k, v in manifest.items()
                        if k not in {'inputs', 'PRS311_age_scaling', 'x_range_case_PRS_z'}
                    },
                    indent=2,
                )
                + '\n'
            )

            print(
                'Age-priority main figure, carrier and AUC supplements, contrasts, captions and validation complete.',
                flush=True,
            )

        # ============================================================================
        # 2. Plot raw score distributions
        # ============================================================================
        if selected_task == 'score-distributions':
            SCORE_DISTRIBUTION_FIGURE_DIR.mkdir(parents=True, exist_ok=True)

            sources = [
                MS_4_ROOT / 'manuscript/datas/burden_13gene_case_control.csv',
                MS_4_ROOT / 'manuscript/datas/PRS_scores.csv',
                MS_4_ROOT / 'data/29_rowlands_extension/section9/frozen_model_inputs.csv.gz',
            ]

            d = pd.read_csv(sources[2]).set_index('person_id')

            b = pd.read_csv(sources[0]).set_index('sample_id')

            p = pd.read_csv(sources[1]).set_index('person_id')

            assert d.index.is_unique and b.index.is_unique and p.index.is_unique

            assert set(d.index) == set(b.index) and set(d.index) <= set(p.index)

            b = b.loc[d.index]

            p = p.loc[d.index]

            assert len(d) == 141204 and int(d.is_case.sum()) == 8528

            assert np.allclose(p.PRS311, d.PRS311, rtol=1e-12, atol=1e-12)

            versions = [
                ('12 genes\nWithout STK11', 'B_raw_12gene_without_STK11_main'),
                ('13 genes\nMain', 'B_raw_13gene_main'),
                (
                    '12 genes\nStrict TP53/NF1 VAF',
                    'B_raw_12gene_without_STK11_tp53_nf1_vaf_sensitive',
                ),
                ('13 genes\nStrict TP53/NF1 VAF', 'B_raw_13gene_tp53_nf1_vaf_sensitive'),
            ]

            assert np.isfinite(b[[col for _, col in versions]].to_numpy()).all()

            labels = [(1, 'Case', MS_4_BLUE, -0.19), (0, 'Control', MS_4_RED, 0.19)]

            plt.style.use('default')

            fig, axes = plt.subplots(2, 1, figsize=(11, 8.2), sharex=True)

            with threadpool_limits(limits=1):
                for ax, positive_only in zip(axes, [False, True]):
                    for j, (_, column) in enumerate(versions):
                        for case, label, color, offset in labels:
                            values = b.loc[d.is_case.eq(case), column].to_numpy()
                            if positive_only:
                                values = values[values > 0]
                            assert len(values) >= 20
                            grid = np.linspace(values.min(), values.max(), 350)
                            kde = compressed_kde(values, grid)
                            widths = 0.16 * kde / kde.max()
                            center = j + offset
                            ax.fill_betweenx(
                                grid,
                                center - widths,
                                center + widths,
                                color=color,
                                alpha=0.75,
                                linewidth=0.8,
                                edgecolor=color,
                            )
                            median = np.median(values)
                            ax.plot(
                                [center - 0.045, center + 0.045],
                                [median, median],
                                color='black',
                                linewidth=1.5,
                            )
                    ax.set_ylabel('Raw weighted burden score')
                    ax.grid(False)
                    ax.spines[['top', 'right']].set_visible(False)
                    ax.set_ylim(bottom=-0.1)
                axes[0].set_title(
                    'A  All participants, including zero scores', loc='left', fontsize=12
                )
                axes[1].set_title('B  Non-zero scores only (carriers)', loc='left', fontsize=12)
                axes[0].legend(
                    handles=[
                        Patch(facecolor=MS_4_BLUE, alpha=0.75, label='Case'),
                        Patch(facecolor=MS_4_RED, alpha=0.75, label='Control'),
                    ],
                    frameon=False,
                    loc='upper right',
                )
                axes[1].set_xticks(range(len(versions)), [label for label, _ in versions])
                axes[1].set_xlabel('Burden definition / sensitivity analysis')
                fig.suptitle('13-gene PTV + eligible P/LP missense burden', fontsize=14)
                fig.text(
                    0.5,
                    0.015,
                    'Same frozen case-control cohort. Violin widths normalized separately; black bars denote medians.\nZero-heavy scores concentrate near 0 in panel A; panel B shows the carrier distribution.',
                    ha='center',
                    fontsize=9,
                )
                fig.tight_layout(rect=[0, 0.065, 1, 0.95])
                burden_preview = save_score_distribution_figure(fig, 'burden_raw_score_violin')
                fig, ax = plt.subplots(figsize=(8, 5))
                values = p.PRS311.to_numpy(dtype=float)
                assert np.isfinite(values).all()
                grid = np.linspace(values.min(), values.max(), 650)
                for case, label, color, _ in labels:
                    v = p.loc[d.is_case.eq(case), 'PRS311'].to_numpy(dtype=float)
                    kde = gaussian_kde(v)(grid)
                    ax.plot(grid, kde, color=color, linewidth=2, label=label)
                    ax.fill_between(grid, 0, kde, color=color, alpha=0.12)
                ax.set_xlabel('Raw PRS311 score')
                ax.set_ylabel('Probability density')
                ax.set_title('PRS311 distribution in cases and controls')
                ax.legend(frameon=False)
                ax.grid(False)
                ax.spines[['top', 'right']].set_visible(False)
                ax.set_ylim(bottom=0)
                fig.tight_layout()
                prs_preview = save_score_distribution_figure(fig, 'PRS311_raw_score_KDE')

            manifest = {
                'cohort': 'frozen case-control',
                'n': len(d),
                'case_n': int(d.is_case.sum()),
                'scores': 'raw, not z-standardized',
                'burden_versions': [c for _, c in versions],
                'burden_KDE': 'Gaussian Scott bandwidth of full selected sample; duplicates compressed exactly; display limited to score support',
                'PRS_KDE': 'Gaussian Scott bandwidth separately by group; each density integrates to one over real line',
                'violin_width': 'same maximum width for every group, not scaled by sample size',
                'burden_preview_extension': burden_preview,
                'PRS_preview_extension': prs_preview,
                'internal_only': True,
                'models_refitted': False,
                'input_sha256': {
                    str(q.relative_to(MS_4_ROOT)): hashlib.sha256(q.read_bytes()).hexdigest()
                    for q in sources
                },
            }

            (SCORE_DISTRIBUTION_FIGURE_DIR / 'plot_manifest.json').write_text(
                json.dumps(manifest, indent=2) + '\n'
            )

            print(
                'Burden violin and PRS311 KDE generated from raw scores; participant data not printed.'
            )

        # ============================================================================
        # 3. Plot four model discrimination
        # ============================================================================
        if selected_task == 'discrimination':
            for path in [DISCRIMINATION_CACHE_DIR, MS_5_FIG, MS_5_DATA]:
                path.mkdir(parents=True, exist_ok=True)

            paths = [
                MS_5_BASE / 'section9/frozen_model_inputs.csv.gz',
                MS_5_BASE / 'section8/comparable_scores_case_control.csv.gz',
                MS_5_BASE / 'section10/prespecified_sensitivity_scores_case_control.csv.gz',
                MS_5_ROOT / 'data/16_joint_model_cv_oof_person_level.csv.gz',
                MS_5_EXISTING / 'model_input_lock_internal.json',
                MS_5_ROOT / 'results/report/45_rowlands_models_validation.json',
                MS_5_ROOT / 'results/table/45_rowlands_CV_AUC.csv',
                MS_5_ROOT / 'scripts/functions/stats.py',
                Path(__file__),
            ]

            paths += [MS_5_EXISTING / f'cv_repeat_{r:02d}.npz' for r in range(1, 11)]

            hashes = {str(p.relative_to(MS_5_ROOT)): MS_5_sha(p) for p in paths}

            run_hash = hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest()

            historical_lock = json.loads(paths[4].read_text())

            for p, digest in historical_lock['inputs'].items():
                if p.startswith('data/'):
                    assert MS_5_sha(MS_5_ROOT / p) == digest, 'Existing model input changed: ' + p

            status = json.loads(paths[5].read_text())

            existing_names = ['PRS311_baseline'] + status['versions']

            d = pd.read_csv(paths[0]).set_index('person_id')

            old = pd.read_csv(paths[3]).set_index('person_id')

            assert d.index.equals(old.index) and len(d) == 141204

            scores = (
                pd.read_csv(paths[1])
                .set_index('sample_id')
                .join(pd.read_csv(paths[2]).set_index('sample_id'))
                .loc[d.index]
            )

            raw = scores[[f'B_raw_{v}' for v, _ in DISCRIMINATION_VERSIONS]].to_numpy(dtype=float)

            y = d.is_case.to_numpy(dtype=int)

            cov = pd.concat(
                [
                    d[['age_at_last_obs'] + [f'PC{i}' for i in range(1, 11)]],
                    pd.get_dummies(d.seq_site_id, prefix='site', drop_first=True).astype(float),
                    pd.get_dummies(d.seq_sequencer, prefix='seq', drop_first=True).astype(float),
                ],
                axis=1,
            ).to_numpy(dtype=float)

            assert y.sum() == 8528 and np.isfinite(cov).all() and np.isfinite(raw).all()

            strata = y.astype(str) + '_' + old.carrier_main_10gene.to_numpy(dtype=int).astype(str)

            splits = list(
                RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=20260818).split(
                    np.arange(len(y)), strata
                )
            )

            fold = np.full(len(y), -1)

            for k, (_, te) in enumerate(splits[:5]):
                fold[te] = k

            assert np.array_equal(fold, old.fold_r1)

            (DISCRIMINATION_CACHE_DIR / 'input_lock_internal.json').write_text(
                json.dumps(
                    {
                        'inputs': hashes,
                        'run_hash': run_hash,
                        'new_models': 'covariates-only and burden-only for six versions; PRS and joint reused',
                        'splits': 'original outcome plus legacy carrier strata, 10x5, random_state=20260818',
                    },
                    indent=2,
                )
                + '\n'
            )

            auc_rows, delta_rows, roc_rows = ([], [], [])

            repeat1 = {}

            reference_auc = pd.read_csv(paths[6]).set_index('version')

            with threadpool_limits(limits=1), ThreadPoolExecutor(max_workers=3) as pool:
                for repeat in range(1, 11):
                    with np.load(MS_5_EXISTING / f'cv_repeat_{repeat:02d}.npz') as saved:
                        assert str(saved['run_hash']) == historical_lock['run_hash']
                        reused = saved['predictions']
                        assert reused.shape == (len(existing_names), len(y))
                    checkpoint = DISCRIMINATION_CACHE_DIR / f'new_models_repeat_{repeat:02d}.npz'
                    if checkpoint.exists():
                        with np.load(checkpoint) as saved:
                            assert str(saved['run_hash']) == run_hash
                            new = saved['predictions']
                    else:
                        new = np.full((1 + len(DISCRIMINATION_VERSIONS), len(y)), np.nan)
                        for k in range(5):
                            tr, te = splits[(repeat - 1) * 5 + k]
                            tasks = [
                                pool.submit(
                                    fit_discrimination_model,
                                    y[tr],
                                    discrimination_design_matrix(cov[tr]),
                                    discrimination_design_matrix(cov[te]),
                                )
                            ]
                            for j in range(len(DISCRIMINATION_VERSIONS)):
                                mean, sd = (raw[tr, j].mean(), raw[tr, j].std(ddof=1))
                                assert sd > 0
                                tasks.append(
                                    pool.submit(
                                        fit_discrimination_model,
                                        y[tr],
                                        discrimination_design_matrix(
                                            cov[tr], (raw[tr, j] - mean) / sd
                                        ),
                                        discrimination_design_matrix(
                                            cov[te], (raw[te, j] - mean) / sd
                                        ),
                                    )
                                )
                            for j, task in enumerate(tasks):
                                new[j, te] = task.result()
                            print(
                                f'Missing-model CV repeat {repeat}/10 fold {k + 1}/5 completed.',
                                flush=True,
                            )
                        assert np.isfinite(new).all() and ((new >= 0) & (new <= 1)).all()
                        np.savez_compressed(checkpoint, run_hash=run_hash, predictions=new)
                    assert (
                        new.shape == (1 + len(DISCRIMINATION_VERSIONS), len(y))
                        and np.isfinite(new).all()
                    )
                    if repeat == 1:
                        assert np.allclose(
                            new[0], old.oof_r1_covariates_only, atol=1e-08, rtol=1e-08
                        )
                        assert np.allclose(reused[0], old.oof_r1_prs_only, atol=1e-08, rtol=1e-08)
                    for j, (version, label) in enumerate(DISCRIMINATION_VERSIONS):
                        predictions = [
                            new[0],
                            reused[0],
                            new[j + 1],
                            reused[existing_names.index(version)],
                        ]
                        for model, pred in zip(DISCRIMINATION_MODELS, predictions):
                            value = roc_auc_score(y, pred)
                            low, high = auc_ci(y, pred)
                            auc_rows.append(
                                {
                                    'version': version,
                                    'model': model,
                                    'repeat': repeat,
                                    'AUC': value,
                                    'CI_low': low,
                                    'CI_high': high,
                                }
                            )
                        for a, b in [(3, 1), (3, 2), (3, 0), (1, 0), (2, 0)]:
                            delta_rows.append(
                                {
                                    'version': version,
                                    'repeat': repeat,
                                    'comparison': DISCRIMINATION_MODELS[a]
                                    + ' minus '
                                    + DISCRIMINATION_MODELS[b],
                                    **delong_test(y, predictions[a], predictions[b]),
                                }
                            )
                        if repeat == 1:
                            repeat1[version] = predictions
                            assert np.isclose(
                                roc_auc_score(y, predictions[3]),
                                reference_auc.loc[version, 'AUC_repeat1'],
                                atol=1e-12,
                                rtol=1e-12,
                            )
                            fpr_grid = np.linspace(0, 1, 1001)
                            for model, pred in zip(DISCRIMINATION_MODELS, predictions):
                                fpr, tpr, _ = roc_curve(y, pred)
                                interp = np.interp(fpr_grid, fpr, tpr)
                                for f, t in zip(fpr_grid, interp):
                                    roc_rows.append(
                                        {
                                            'version': version,
                                            'model': model,
                                            'repeat': 1,
                                            'FPR': f,
                                            'TPR_interpolated': t,
                                        }
                                    )
                    print(
                        f'Four-model AUC and paired DeLong results repeat {repeat}/10 completed.',
                        flush=True,
                    )

            auc = pd.DataFrame(auc_rows)

            summary = (
                auc.groupby(['version', 'model'], sort=False)
                .AUC.agg(AUC_mean='mean', AUC_min='min', AUC_max='max')
                .reset_index()
            )

            summary = summary.merge(
                auc.loc[auc.repeat.eq(1)]
                .drop(columns='repeat')
                .rename(
                    columns={
                        'AUC': 'AUC_repeat1',
                        'CI_low': 'repeat1_CI_low',
                        'CI_high': 'repeat1_CI_high',
                    }
                ),
                on=['version', 'model'],
            )

            for version, _ in DISCRIMINATION_VERSIONS:
                joint = summary.loc[summary.version.eq(version) & summary.model.eq('Joint')].iloc[0]
                assert np.isclose(
                    joint.AUC_mean, reference_auc.loc[version, 'AUC_mean'], rtol=1e-12, atol=1e-12
                )

            auc.to_csv(MS_5_DATA / 'four_models_AUC_by_repeat.csv', index=False)

            summary.to_csv(MS_5_DATA / 'four_models_AUC_summary.csv', index=False)

            pd.DataFrame(delta_rows).to_csv(
                MS_5_DATA / 'four_models_paired_delta_AUC_by_repeat.csv', index=False
            )

            pd.DataFrame(roc_rows).to_csv(
                MS_5_DATA / 'four_models_ROC_repeat1_interpolated.csv', index=False
            )

            plt.style.use('default')

            fig, axes = plt.subplots(2, 3, figsize=(15, 9), sharex=True, sharey=True)

            for j, ((version, label), ax) in enumerate(zip(DISCRIMINATION_VERSIONS, axes.flat)):
                for model, pred, color in zip(DISCRIMINATION_MODELS, repeat1[version], MS_5_COLORS):
                    fpr, tpr, _ = roc_curve(y, pred)
                    ax.plot(fpr, tpr, color=color, linewidth=1.6)
                ax.plot([0, 1], [0, 1], color='.7', linestyle='--', linewidth=0.8)
                ax.set_title(f'({chr(97 + j)}) {label}', loc='left', fontsize=11)
                ax.set_xlim(0, 1)
                ax.set_ylim(0, 1)
                ax.grid(False)
                ax.spines[['top', 'right']].set_visible(False)
                if j >= 3:
                    ax.set_xlabel('False positive rate')
                if j % 3 == 0:
                    ax.set_ylabel('True positive rate')

            fig.legend(
                [Line2D([0], [0], color=c, linewidth=2) for c in MS_5_COLORS],
                DISCRIMINATION_MODELS,
                loc='upper center',
                bbox_to_anchor=(0.5, 0.955),
                ncol=4,
                frameon=False,
            )

            fig.suptitle(
                'Discrimination of four case-status models: ROC curves', fontsize=15, y=0.995
            )

            fig.text(
                0.5,
                0.01,
                'First-repeat out-of-fold predictions; same frozen cohort and paired folds across all six versions.',
                ha='center',
                fontsize=10,
            )

            fig.tight_layout(rect=[0, 0.04, 1, 0.92])

            roc_preview = save_discrimination_figure(
                fig, 'four_case_status_models_ROC_sensitivities'
            )

            fig, axes = plt.subplots(2, 3, figsize=(15, 9), sharey=True)

            ymin = 0.5

            ymax = max(summary.AUC_mean.max(), summary.repeat1_CI_high.max()) + 0.012

            for j, ((version, label), ax) in enumerate(zip(DISCRIMINATION_VERSIONS, axes.flat)):
                rows = (
                    summary.loc[summary.version.eq(version)]
                    .set_index('model')
                    .loc[DISCRIMINATION_MODELS]
                )
                xx = np.arange(4)
                ax.bar(
                    xx, rows.AUC_mean - ymin, bottom=ymin, color=MS_5_COLORS, width=0.58, alpha=0.85
                )
                ax.vlines(
                    xx, rows.repeat1_CI_low, rows.repeat1_CI_high, color='black', linewidth=1.2
                )
                ax.hlines(rows.repeat1_CI_low, xx - 0.075, xx + 0.075, color='black', linewidth=1.2)
                ax.hlines(
                    rows.repeat1_CI_high, xx - 0.075, xx + 0.075, color='black', linewidth=1.2
                )
                for x, mean, high in zip(xx, rows.AUC_mean, rows.repeat1_CI_high):
                    ax.text(
                        x,
                        max(mean, high) + 0.002,
                        f'{mean:.4f}',
                        ha='center',
                        va='bottom',
                        fontsize=9,
                    )
                ax.set_xticks(
                    xx,
                    ['Covariates', 'PRS + cov.', 'Burden + cov.', 'Joint'],
                    rotation=25,
                    ha='right',
                    fontsize=9,
                )
                ax.set_ylim(ymin, ymax)
                ax.set_title(f'({chr(97 + j)}) {label}', loc='left', fontsize=11)
                ax.grid(False)
                ax.spines[['top', 'right']].set_visible(False)
                if j % 3 == 0:
                    ax.set_ylabel('Mean out-of-fold AUC')

            fig.suptitle('Discrimination of four case-status models: AUC', fontsize=15)

            fig.text(
                0.5,
                0.015,
                'Bars: mean AUC over 10 repeats of 5-fold CV. Error bars: first-repeat DeLong 95% CIs.\nY-axis starts at 0.50. First-repeat intervals describe that repeat, not uncertainty in the ten-repeat mean.',
                ha='center',
                fontsize=10,
            )

            fig.tight_layout(rect=[0, 0.065, 1, 0.955])

            auc_preview = save_discrimination_figure(
                fig, 'four_case_status_models_AUC_sensitivities'
            )

            for p, digest in hashes.items():
                assert MS_5_sha(MS_5_ROOT / p) == digest, 'Input changed during fitting: ' + p

            validation = {
                'status': 'four_model_discrimination_completed',
                'versions': [v for v, _ in DISCRIMINATION_VERSIONS],
                'models': DISCRIMINATION_MODELS,
                'n': len(y),
                'case_n': int(y.sum()),
                'CV_repeats': 10,
                'CV_folds': 5,
                'PRS_and_joint_reused': True,
                'new_models_fitted': 'Covariates only and burden + covariates',
                'all_fits_converged': True,
                'original_first_repeat_folds_reproduced': True,
                'covariates_first_repeat_OOF_reproduced': True,
                'all_joint_AUCs_match_previous_results': True,
                'score_scaling': 'training fold mean and SD ddof=1; no validation fold fitting',
                'joint_definition': 'PRS + burden + PRS*burden + age_at_last_obs + PC1-10 + sequencing site and sequencer',
                'ROC_repeat': 1,
                'AUC_bar_statistic': '10-repeat mean',
                'AUC_CI': 'Exact first-repeat DeLong CI endpoints',
                'CI_limit': 'conditional on fitted OOF predictions; excludes full retraining uncertainty',
                'final_variant_definition_frozen': False,
                'ROC_preview': roc_preview,
                'AUC_preview': auc_preview,
                'inputs_unchanged': True,
                'participant_predictions_internal_only': True,
            }

            (MS_5_ROOT / 'results/report/51_four_model_discrimination_validation.json').write_text(
                json.dumps(validation, indent=2) + '\n'
            )

            (MS_5_FIG / 'four_case_status_models_captions.md').write_text(
                '# Four case-status model discrimination\n\nEach figure contains six panels: main C; without STK11; PTV-only; strict TP53/NF1 VAF; ClinVar ≥1 star; without CDH1.\n\n(a) ROC curves of first-repeat out-of-fold predictions from the four case-status models. (b) Mean out-of-fold AUC over 10 repeats of 5-fold cross-validation. Error bars are DeLong 95% CIs for the first repeat. These intervals are not confidence intervals for the ten-repeat mean. All models include the same clinical/technical covariates; the joint model also includes PRS×burden. Same frozen cohort (N=141,204; cases=8,528) and original paired folds. Existing PRS and joint predictions reused; missing covariates-only and burden-only predictions fitted. Current variant definition remains provisional.\n'
            )

            print('ROC and AUC figures, unrounded source CSVs and validation saved.', flush=True)

        # ============================================================================
        # 4. Plot joint distribution
        # ============================================================================
        if selected_task == 'joint-distribution':
            for p in [MS_6_FIG, MS_6_DATA, JOINT_DISTRIBUTION_CACHE_DIR]:
                p.mkdir(parents=True, exist_ok=True)

            paths = [
                MS_6_ROOT / 'data/29_rowlands_extension/section9/frozen_model_inputs.csv.gz',
                MS_6_ROOT / 'manuscript/datas/burden_13gene_case_control.csv',
                MS_6_ROOT / 'manuscript/datas/PRS_scores.csv',
                Path(__file__),
            ]

            hashes = {
                str(p.relative_to(MS_6_ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in paths
            }

            d = pd.read_csv(paths[0]).set_index('person_id')

            burden = pd.read_csv(paths[1]).set_index('sample_id').loc[d.index]

            prs = (
                pd.read_csv(paths[2], usecols=['person_id', 'PRS311'])
                .set_index('person_id')
                .loc[d.index]
            )

            assert d.index.is_unique and burden.index.is_unique and prs.index.is_unique

            assert len(d) == 141204 and d.is_case.sum() == 8528

            assert np.allclose(d.PRS311, prs.PRS311, atol=1e-12, rtol=1e-12)

            carrier = burden.B_raw_13gene_main.gt(0).astype(int)

            assert carrier.equals(burden.carried_gene_count_13gene_main.gt(0).astype(int))

            decile = pd.qcut(d.PRS311, 10, labels=False) + 1

            pc = [f'PC{i}' for i in range(1, 11)]

            X = pd.DataFrame(
                {
                    'const': 1.0,
                    'PRS311_z': (d.PRS311 - d.PRS311.mean()) / d.PRS311.std(ddof=1),
                    'C_carrier': carrier,
                },
                index=d.index,
            ).join(d[['age_at_last_obs'] + pc])

            assert np.isfinite(X.to_numpy()).all() and (not any(('interaction' in c for c in X)))

            with threadpool_limits(limits=1):
                fit = sm.GLM(d.is_case, X, family=sm.families.Binomial()).fit(cov_type='HC1')

            assert fit.converged and np.isfinite(fit.params).all()

            predicted = fit.predict(X)

            assert np.isclose(predicted.sum(), d.is_case.sum(), atol=1e-05, rtol=0)

            coef = pd.DataFrame(
                {
                    'term': X.columns,
                    'beta': fit.params.to_numpy(),
                    'SE_HC1': fit.bse.to_numpy(),
                    'OR': np.exp(fit.params.to_numpy()),
                    'CI_low': np.exp(fit.conf_int()[0].to_numpy()),
                    'CI_high': np.exp(fit.conf_int()[1].to_numpy()),
                    'p_value': fit.pvalues.to_numpy(),
                }
            )

            coef.to_csv(
                MS_6_DATA / 'joint_distribution_no_interaction_coefficients.csv', index=False
            )

            rows = []

            for k in range(1, 11):
                for group in [0, 1]:
                    mask = decile.eq(k) & carrier.eq(group)
                    n = int(mask.sum())
                    cases = int(d.loc[mask].is_case.sum())
                    controls = n - cases
                    expected = float(predicted.loc[mask].sum())
                    assert n > 0 and expected > 0
                    rows.append(
                        {
                            'PRS_decile': k,
                            'C_carrier': group,
                            'participants': n,
                            'cases': cases,
                            'controls': controls,
                            'observed_case_proportion': cases / n,
                            'expected_cases_no_interaction': expected,
                            'mean_expected_probability': expected / n,
                            'observed_over_expected': cases / expected,
                            'disclosure_safe': min(n, cases, controls) >= 20,
                        }
                    )

            cells = pd.DataFrame(rows)

            assert cells.participants.sum() == len(d) and cells.cases.sum() == d.is_case.sum()

            assert np.isclose(
                cells.expected_cases_no_interaction.sum(), d.is_case.sum(), atol=1e-05, rtol=0
            )

            for group in [0, 1]:
                sub = cells.loc[cells.C_carrier.eq(group)]
                assert np.isclose(
                    sub.cases.sum(), sub.expected_cases_no_interaction.sum(), atol=1e-05, rtol=0
                )

            cells.to_csv(JOINT_DISTRIBUTION_CACHE_DIR / 'cells_exact_internal.csv', index=False)

            hidden_deciles = set(cells.loc[~cells.disclosure_safe, 'PRS_decile'])

            public = cells.drop(columns='disclosure_safe').copy()

            public['cell_suppressed'] = public.PRS_decile.isin(hidden_deciles)

            columns = [
                'participants',
                'cases',
                'controls',
                'observed_case_proportion',
                'expected_cases_no_interaction',
                'mean_expected_probability',
                'observed_over_expected',
            ]

            for c in columns:
                public[c] = public[c].astype(object)
                public.loc[public.cell_suppressed, c] = 'suppressed'

            public.to_csv(MS_6_DATA / 'PRS311_decile_C_carrier_joint_distribution.csv', index=False)

            observed = (
                cells.pivot(
                    index='PRS_decile', columns='C_carrier', values='observed_case_proportion'
                ).to_numpy()
                * 100
            )

            ratio = cells.pivot(
                index='PRS_decile', columns='C_carrier', values='observed_over_expected'
            ).to_numpy()

            for k in hidden_deciles:
                observed[k - 1, :] = np.nan
                ratio[k - 1, :] = np.nan

            plt.style.use('default')

            fig, axes = plt.subplots(1, 2, figsize=(11.5, 7.5), sharey=True)

            image_a = axes[0].imshow(
                observed,
                origin='lower',
                aspect='auto',
                cmap='Reds',
                vmin=0,
                vmax=np.nanmax(observed),
            )

            low, high = (np.nanmin(ratio), np.nanmax(ratio))

            span = max(abs(low - 1), abs(high - 1), 0.05)

            image_b = axes[1].imshow(
                ratio,
                origin='lower',
                aspect='auto',
                cmap='RdBu_r',
                norm=TwoSlopeNorm(vmin=1 - span, vcenter=1, vmax=1 + span),
            )

            fig.colorbar(
                image_a, ax=axes[0], label='Observed case proportion (%)', fraction=0.05, pad=0.04
            )

            fig.colorbar(
                image_b, ax=axes[1], label='Observed / expected cases', fraction=0.05, pad=0.04
            )

            maximum = np.nanmax(observed)

            for row in cells.itertuples():
                y, x = (row.PRS_decile - 1, row.C_carrier)
                if row.PRS_decile in hidden_deciles:
                    for ax in axes:
                        ax.text(x, y, 'Suppressed', ha='center', va='center', fontsize=9)
                    continue
                color = 'white' if observed[y, x] > maximum * 0.55 else 'black'
                axes[0].text(
                    x,
                    y,
                    f'{observed[y, x]:.2f}%\nn={row.participants:,}; cases={row.cases:,}',
                    ha='center',
                    va='center',
                    fontsize=8.5,
                    color=color,
                )
                color = 'white' if abs(ratio[y, x] - 1) > span * 0.65 else 'black'
                axes[1].text(
                    x, y, f'{ratio[y, x]:.2f}', ha='center', va='center', fontsize=10, color=color
                )

            for ax in axes:
                ax.set_xticks([0, 1], ['Non-carrier', 'Carrier'])
                ax.set_yticks(np.arange(10), np.arange(1, 11))
                ax.set_xlabel('13-gene PTV + eligible missense carrier status')
                ax.grid(False)
                for spine in ax.spines.values():
                    spine.set_visible(False)

            axes[0].set_ylabel('PRS311 decile (full cohort)')

            axes[0].set_title('(a) Observed case proportion', loc='left', fontsize=12)

            axes[1].set_title('(b) Observed / expected cases', loc='left', fontsize=12)

            fig.suptitle(
                'Joint distribution of PRS311 and 13-gene carrier status\nN = 141,204', fontsize=14
            )

            fig.text(
                0.5,
                0.02,
                'Left: cohort-observed proportions, not absolute risks. Right: expected cases from a logistic model\nwith PRS311 + binary carrier + age at last observation + PC1–10, without interaction.',
                ha='center',
                fontsize=9,
            )

            fig.tight_layout(rect=[0, 0.075, 1, 0.94])

            stem = MS_6_FIG / 'PRS311_decile_13gene_carrier_joint_distribution'

            fig.savefig(stem.with_suffix('.svg'), bbox_inches='tight')

            fig.savefig(stem.with_suffix('.png'), dpi=180, bbox_inches='tight')

            fig.savefig(
                stem.with_suffix('.jpg'),
                dpi=180,
                bbox_inches='tight',
                pil_kwargs={'quality': 85, 'subsampling': 2},
            )

            preview = (
                '.png'
                if stem.with_suffix('.png').stat().st_size
                <= stem.with_suffix('.jpg').stat().st_size
                else '.jpg'
            )

            stem.with_suffix('.jpg' if preview == '.png' else '.png').unlink()

            plt.close(fig)

            for p, digest in hashes.items():
                assert hashlib.sha256((MS_6_ROOT / p).read_bytes()).hexdigest() == digest

            result = {
                'status': 'C_joint_distribution_heatmap_completed',
                'n': len(d),
                'case_n': int(d.is_case.sum()),
                'carrier_definition': 'At least one main-C qualifying PTV or P/LP missense in 13 genes',
                'PRS_deciles': 'pd.qcut of raw PRS311 in full frozen case-control cohort',
                'expected_model': 'unweighted GLM Binomial HC1: PRS311_z + C_carrier + age_at_last_obs + PC1-10; no interaction',
                'expected_model_matches_legacy_adjustment': True,
                'expected_predictions': 'full-cohort fitted, not OOF',
                'model_converged': True,
                'observed_expected_totals_and_carrier_margins_agree': True,
                'raw_scores_match_frozen_input_PRS': True,
                'small_outcome_cells_and_complements_suppressed': True,
                'final_variant_definition_frozen': False,
                'preview': str(stem.with_suffix(preview).relative_to(MS_6_ROOT)),
            }

            (MS_6_ROOT / 'results/report/52_C_joint_distribution_validation.json').write_text(
                json.dumps(result, indent=2) + '\n'
            )

            (JOINT_DISTRIBUTION_CACHE_DIR / 'input_and_model_manifest.json').write_text(
                json.dumps({'inputs': hashes, **result}, indent=2) + '\n'
            )

            (MS_6_FIG / 'PRS311_decile_13gene_carrier_joint_distribution_caption.md').write_text(
                'Joint distribution of PRS311 and 13-gene PTV plus eligible P/LP missense carrier status (N=141,204). Left: observed case proportion by full-cohort PRS311 decile and binary carrier status. These proportions depend on the case-control composition of the cohort and are not absolute risks. Right: observed/expected cases; expected counts are sums of fitted probabilities from a logistic model without interaction, adjusted for continuous PRS311, binary carrier status, age at last EHR observation and PC1–10. The adjustment matches the prior figure; sequencing site and sequencer are not included in this descriptive expectation model. O/E differences do not independently establish biological interaction. Current variant definition is provisional.\n'
            )

            print(
                'Updated C-carrier joint-distribution heatmap and unrounded aggregate CSVs saved; model and cell totals validated.'
            )

        # ============================================================================
        # 5. Plot score decile boxplots
        # ============================================================================
        if selected_task == 'score-deciles':
            MS_7_FIG.mkdir(exist_ok=True)

            SCORE_DECILE_CACHE_DIR.mkdir(parents=True, exist_ok=True)

            paths = [
                MS_7_ROOT / 'data/29_rowlands_extension/section9/frozen_model_inputs.csv.gz',
                MS_7_ROOT / 'manuscript/datas/PRS_scores.csv',
                MS_7_ROOT / 'manuscript/datas/burden_13gene_case_control.csv',
            ]

            d = pd.read_csv(paths[0]).set_index('person_id')

            p = (
                pd.read_csv(paths[1], usecols=['person_id', 'PRS311'])
                .set_index('person_id')
                .loc[d.index]
            )

            b = (
                pd.read_csv(paths[2], usecols=['sample_id', 'B_raw_13gene_main'])
                .set_index('sample_id')
                .loc[d.index]
            )

            assert len(d) == 141204 and int(d.is_case.sum()) == 8528

            assert np.allclose(d.PRS311, p.PRS311, atol=1e-12, rtol=1e-12)

            decile = pd.qcut(d.PRS311, 10, labels=False) + 1

            labels = [f'{10 * i}–{10 * (i + 1)}%' for i in range(10)]

            labels[0] += '\nLowest'

            labels[-1] += '\nHighest'

            colors = [(1, 'Case', '#2463B4', -0.18), (0, 'Control', '#D54848', 0.18)]

            plt.style.use('default')

            fig, axes = plt.subplots(2, 1, figsize=(12, 7.5), sharex=True)

            specs = [
                ('PRS311', p.PRS311, False),
                ('13-gene burden: non-zero scores only', b.B_raw_13gene_main, True),
            ]

            summary = []

            for ax, (name, values, positive_only) in zip(axes, specs):
                for case, group, color, offset in colors:
                    samples = []
                    for k in range(1, 11):
                        mask = decile.eq(k) & d.is_case.eq(case)
                        if positive_only:
                            mask &= values.gt(0)
                        v = values.loc[mask].to_numpy()
                        assert len(v) >= 20
                        samples.append(v)
                        summary.append(
                            {
                                'panel': name,
                                'PRS_decile': k,
                                'group': group,
                                'n': len(v),
                                'q25': float(np.quantile(v, 0.25)),
                                'median': float(np.median(v)),
                                'q75': float(np.quantile(v, 0.75)),
                            }
                        )
                    boxes = ax.boxplot(
                        samples,
                        positions=np.arange(1, 11) + offset,
                        widths=0.3,
                        patch_artist=True,
                        showfliers=False,
                        manage_ticks=False,
                        medianprops={'color': 'black', 'linewidth': 1.2},
                        boxprops={'edgecolor': color, 'linewidth': 1.1},
                        whiskerprops={'color': color},
                        capprops={'color': color},
                    )
                    for patch in boxes['boxes']:
                        patch.set_facecolor(color)
                        patch.set_alpha(0.65)
                ax.set_title(name, loc='left', fontsize=12)
                ax.set_ylabel(
                    'Raw PRS311 score' if name == 'PRS311' else 'Raw weighted burden score'
                )
                ax.grid(False)
                ax.spines[['top', 'right']].set_visible(False)
                if name != 'PRS311':
                    ax.set_ylim(bottom=-0.12)

            axes[0].legend(
                handles=[Patch(facecolor=c, alpha=0.65, label=g) for _, g, c, _ in colors],
                frameon=False,
                loc='upper left',
            )

            axes[-1].set_xticks(np.arange(1, 11), labels)

            axes[-1].set_xlabel(
                'PRS311 percentile group in the full case-control cohort (deciles 1–10)'
            )

            axes[-1].set_xlim(0.4, 10.6)

            fig.suptitle('Raw PRS and burden scores across PRS311 deciles', fontsize=15)

            fig.text(
                0.5,
                0.01,
                'Blue: case; red: control. Boxes: median and IQR; whiskers: within 1.5×IQR. Outliers are not shown.\nPercentile labels denote intervals, not cumulative groups. Bottom panel includes carriers only.',
                ha='center',
                fontsize=9,
            )

            fig.tight_layout(rect=[0, 0.055, 1, 0.96])

            stem = MS_7_FIG / 'PRS311_decile_PRS_burden_boxplots'

            fig.savefig(stem.with_suffix('.svg'), bbox_inches='tight')

            fig.savefig(stem.with_suffix('.png'), dpi=180, bbox_inches='tight')

            fig.savefig(
                stem.with_suffix('.jpg'),
                dpi=180,
                bbox_inches='tight',
                pil_kwargs={'quality': 85, 'subsampling': 2},
            )

            ext = (
                '.png'
                if stem.with_suffix('.png').stat().st_size
                <= stem.with_suffix('.jpg').stat().st_size
                else '.jpg'
            )

            stem.with_suffix('.jpg' if ext == '.png' else '.png').unlink()

            plt.close(fig)

            pd.DataFrame(summary).to_csv(
                SCORE_DECILE_CACHE_DIR / 'boxplot_statistics_internal.csv', index=False
            )

            manifest = {
                'cohort_n': len(d),
                'deciles': 'PRS311 qcut in full frozen cohort',
                'scores': 'raw PRS311 and B_raw_13gene_main',
                'panels': 'PRS; nonzero burden',
                'case_control_colors': 'blue/red',
                'outliers_displayed': False,
                'internal_review_only': True,
                'preview': str(stem.with_suffix(ext).relative_to(MS_7_ROOT)),
                'inputs': {
                    str(q.relative_to(MS_7_ROOT)): hashlib.sha256(q.read_bytes()).hexdigest()
                    for q in paths
                },
            }

            (SCORE_DECILE_CACHE_DIR / 'manifest.json').write_text(
                json.dumps(manifest, indent=2) + '\n'
            )

            print('Raw score decile boxplots saved; original case-proportion/OE heatmap retained.')

        # ============================================================================
        # 6. Plot burden histogram
        # ============================================================================
        if selected_task == 'burden-histogram':
            MS_8_FIG.mkdir(exist_ok=True)

            BURDEN_FIGURE_DATA_DIR.mkdir(parents=True, exist_ok=True)

            inputs = [
                MS_8_ROOT / 'data/29_rowlands_extension/section9/frozen_model_inputs.csv.gz',
                MS_8_ROOT / 'manuscript/datas/burden_13gene_case_control.csv',
            ]

            d = pd.read_csv(inputs[0]).set_index('person_id')

            b = pd.read_csv(inputs[1], usecols=['sample_id', 'B_raw_13gene_main']).set_index(
                'sample_id'
            )

            assert d.index.is_unique and b.index.is_unique and (set(d.index) == set(b.index))

            score = b.loc[d.index, 'B_raw_13gene_main']

            assert len(d) == 141204 and int(d.is_case.sum()) == 8528

            assert np.isfinite(score).all() and score.ge(0).all()

            edges = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, float('inf')]

            groups = [(1, 'Case', '#2463B4'), (0, 'Control', '#D54848')]

            def counts(e):
                return [
                    np.histogram(score.loc[d.is_case.eq(case) & score.gt(0)], bins=e)[0]
                    for case, _, _ in groups
                ]

            while True:
                cnt = counts(edges)
                unsafe = [
                    i for i in range(len(edges) - 1) if any((0 < int(c[i]) < 20 for c in cnt))
                ]
                if not unsafe:
                    break
                i = unsafe[-1]
                if i > 0:
                    del edges[i]
                elif len(edges) > 2:
                    del edges[1]
                else:
                    raise AssertionError('Insufficient nonzero participants for grouped release')

            rows = []

            for case, group, _ in groups:
                v = score.loc[d.is_case.eq(case)]
                n = len(v)
                zero = int(v.eq(0).sum())
                assert zero >= 20 and n - zero >= 20
                rows.append(
                    {
                        'group': group,
                        'bin': 'Exactly 0',
                        'participants': zero,
                        'denominator_all_group_participants': n,
                        'percentage': 100 * zero / n,
                    }
                )
                hist = np.histogram(v.loc[v.gt(0)], bins=edges)[0]
                for i, count in enumerate(hist):
                    label = (
                        f'[{edges[i]:g}, {edges[i + 1]:g})'
                        if np.isfinite(edges[i + 1])
                        else f'>= {edges[i]:g}'
                    )
                    if i == 0:
                        label = f'(0, {edges[i + 1]:g})' if np.isfinite(edges[i + 1]) else '> 0'
                    rows.append(
                        {
                            'group': group,
                            'bin': label,
                            'participants': int(count),
                            'denominator_all_group_participants': n,
                            'percentage': 100 * count / n,
                        }
                    )
                assert zero + int(hist.sum()) == n

            table = pd.DataFrame(rows)

            for group in ['Case', 'Control']:
                assert np.isclose(table.loc[table.group.eq(group), 'percentage'].sum(), 100)

            table.to_csv(
                BURDEN_FIGURE_DATA_DIR / '13gene_burden_histogram_percentages.csv', index=False
            )

            plt.style.use('default')

            fig, axes = plt.subplots(
                1, 2, figsize=(12, 4.9), gridspec_kw={'width_ratios': [1, 2.8]}
            )

            for j, (_, group, color) in enumerate(groups):
                r = table.loc[table.group.eq(group)].iloc[0]
                axes[0].bar(j, r.percentage, color=color, alpha=0.75, width=0.55)
                axes[0].text(j, r.percentage + 1, f'{r.percentage:.2f}%', ha='center', fontsize=11)

            axes[0].set_xticks([0, 1], ['Case', 'Control'])

            axes[0].set_ylim(0, 105)

            axes[0].set_ylabel('Participants within group (%)')

            axes[0].set_title('(a) Zero burden', loc='left')

            ticklabels = []

            for i in range(len(edges) - 1):
                ticklabels.append(
                    f'>{edges[i]:g}–<{edges[i + 1]:g}'
                    if i == 0
                    else (
                        f'{edges[i]:g}–<{edges[i + 1]:g}'
                        if np.isfinite(edges[i + 1])
                        else f'≥{edges[i]:g}'
                    )
                )

            xx = np.arange(len(edges) - 1)

            ymax = 0.0

            for j, (_, group, color) in enumerate(groups):
                sub = table.loc[table.group.eq(group)].iloc[1:]
                offset = -0.17 if j == 0 else 0.17
                axes[1].bar(
                    xx + offset, sub.percentage, width=0.32, color=color, alpha=0.75, label=group
                )
                for x, val in zip(xx + offset, sub.percentage):
                    axes[1].text(
                        x, val + 0.025, f'{val:.2f}%', ha='center', va='bottom', fontsize=9
                    )
                ymax = max(ymax, float(sub.percentage.max()))

            axes[1].set_xticks(xx, ticklabels)

            axes[1].set_ylim(0, ymax * 1.22)

            axes[1].set_xlabel('Raw weighted burden score interval')

            axes[1].set_ylabel('Participants within group (%)')

            axes[1].set_title('(b) Non-zero score intervals', loc='left')

            axes[1].legend(frameon=False)

            for ax in axes:
                ax.grid(False)
                ax.spines[['top', 'right']].set_visible(False)

            fig.suptitle(
                '13-gene PTV + eligible missense burden: participant proportions', fontsize=14
            )

            fig.text(
                0.5,
                0.01,
                'Both panels use all cases or all controls as their denominator; non-zero bars are not renormalized among carriers.\nZero and positive-score intervals together sum to 100% within each group. Panel y-axis scales differ.',
                ha='center',
                fontsize=9,
            )

            fig.tight_layout(rect=[0, 0.08, 1, 0.94])

            stem = MS_8_FIG / '13gene_burden_histogram_percentages'

            fig.savefig(stem.with_suffix('.svg'), bbox_inches='tight')

            fig.savefig(stem.with_suffix('.png'), dpi=180, bbox_inches='tight')

            fig.savefig(
                stem.with_suffix('.jpg'),
                dpi=180,
                bbox_inches='tight',
                pil_kwargs={'quality': 85, 'subsampling': 2},
            )

            preview = (
                '.png'
                if stem.with_suffix('.png').stat().st_size
                <= stem.with_suffix('.jpg').stat().st_size
                else '.jpg'
            )

            stem.with_suffix('.jpg' if preview == '.png' else '.png').unlink()

            plt.close(fig)

            report = {
                'status': 'all_participant_C_burden_histogram_completed',
                'denominator': 'All cases or all controls separately',
                'cohort_n': len(d),
                'positive_bins_merged_if_nonzero_counts_below20': True,
                'percentage_totals_verified': True,
                'zero_CI_or_density_not_used': True,
                'preview': str(stem.with_suffix(preview).relative_to(MS_8_ROOT)),
            }

            (MS_8_ROOT / 'results/report/53_burden_histogram_validation.json').write_text(
                json.dumps(report, indent=2) + '\n'
            )

            manifest_path = MS_8_ROOT / 'manuscript/datas/manifest.json'

            manifest = json.loads(manifest_path.read_text())

            name = 'burden/13gene_burden_histogram_percentages.csv'

            manifest['files'] = [r for r in manifest['files'] if r['file'] != name]

            manifest['files'].append(
                {
                    'file': name,
                    'source': 'scripts/13_figures.py --task burden-histogram',
                    'purpose': 'Within-case/control percentages; exact-zero separate; positive bins merged for disclosure',
                    'sha256': hashlib.sha256(
                        (MS_8_ROOT / 'manuscript/datas' / name).read_bytes()
                    ).hexdigest(),
                }
            )

            manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')

            print(
                'Burden percentage histogram and aggregate source CSV saved; percentages sum to 100% within each group.'
            )

        # ============================================================================
        # 7. Plot age burden boxplots
        # ============================================================================
        if selected_task == 'age-boxplots':
            ids = pd.read_csv(MS_9_BASE / 'diagnosis_age_ids.csv.gz').person_id

            d = (
                pd.read_csv(MS_9_BASE / 'section9/frozen_model_inputs.csv.gz')
                .set_index('person_id')
                .loc[ids]
            )

            s = (
                pd.read_csv(MS_9_BASE / 'section8/comparable_scores_case_control.csv.gz')
                .set_index('sample_id')
                .loc[ids]
            )

            b = s.B_raw_C

            age = d.age_at_diagnosis

            assert len(d) == 7996 and d.index.is_unique and d.is_case.eq(1).all()

            assert age.notna().all() and np.isfinite(age).all() and np.isfinite(b).all()

            assert b.ge(0).all() and b.gt(0).equals(s.carried_gene_count_C.gt(0))

            figdir = MS_9_ROOT / 'manuscript/figures'

            figdir.mkdir(parents=True, exist_ok=True)

            rows = []

            manifest = {
                'n': len(d),
                'score': 'B_raw_C: 13-gene PTV + qualifying P/LP missense, sum of ln(OR) per carried gene',
                'outliers_shown': False,
                'variants': {},
            }

            for name, cutoff in [('carrier_median', float(b[b.gt(0)].median()))]:
                masks = [b.gt(0) & b.le(cutoff), b.gt(cutoff), b.eq(0)]
                groups = ['Carrier, lower burden', 'Carrier, higher burden', 'Non-carrier']
                samples = [age.loc[m].to_numpy() for m in masks]
                assert sum(map(len, samples)) == 7996 and min(map(len, samples)) >= 20
                labels = [
                    f'Carrier\nburden ≤ median ({cutoff:.2f})\n(n={len(samples[0]):,})',
                    f'Carrier\nburden > median ({cutoff:.2f})\n(n={len(samples[1]):,})',
                    f'Non-carrier\nburden = 0\n(n={len(samples[2]):,})',
                ]
                plt.style.use('default')
                fig, ax = plt.subplots(figsize=(7, 5))
                ax.boxplot(
                    samples,
                    tick_labels=labels,
                    showfliers=False,
                    widths=0.5,
                    medianprops={'color': 'black'},
                )
                ax.set_ylabel('Age at diagnosis (years)')
                ax.set_title(
                    '(a) Age at diagnosis by burden group\n13-gene age-at-diagnosis cohort (n=7,996)',
                    fontsize=12,
                )
                ax.grid(False)
                ax.spines[['top', 'right']].set_visible(False)
                fig.text(
                    0.5,
                    0.015,
                    'Boxes: median and IQR; whiskers: within 1.5×IQR. Outliers are not shown.',
                    ha='center',
                    fontsize=8,
                )
                fig.tight_layout(rect=[0, 0.045, 1, 1])
                stem = figdir / f'age_at_diagnosis_13gene_burden_{name}'
                sizes = {}
                for ext in ['svg', 'png', 'jpg']:
                    kwargs = (
                        {'pil_kwargs': {'quality': 85, 'subsampling': 2}} if ext == 'jpg' else {}
                    )
                    fig.savefig(stem.with_suffix('.' + ext), dpi=180, bbox_inches='tight', **kwargs)
                    sizes[ext] = stem.with_suffix('.' + ext).stat().st_size
                smallest = min(sizes, key=sizes.get)
                for ext in ['svg', 'jpg']:
                    if ext != smallest:
                        stem.with_suffix('.' + ext).unlink()
                plt.close(fig)
                manifest['variants'][name] = {
                    'cutoff_raw': cutoff,
                    'group_n': list(map(len, samples)),
                    'smallest_format': smallest,
                    'file_sizes': sizes,
                }
                for group, v in zip(groups, samples):
                    rows.append(
                        {
                            'split': name,
                            'cutoff_raw': cutoff,
                            'group': group,
                            'n': len(v),
                            'q25_age': np.quantile(v, 0.25),
                            'median_age': np.median(v),
                            'q75_age': np.quantile(v, 0.75),
                        }
                    )

            pd.DataFrame(rows).to_csv(
                MS_9_ROOT / 'results/table/54_age_13gene_boxplot_summary.csv', index=False
            )

            (MS_9_ROOT / 'results/report/54_age_13gene_boxplot_validation.json').write_text(
                json.dumps(manifest, indent=2) + '\n'
            )

            print(json.dumps(manifest, indent=2))

        # ============================================================================
        # 8. Plot age carrier frequency
        # ============================================================================
        if selected_task == 'age-carrier-frequency':
            MS_10_FIG.mkdir(parents=True, exist_ok=True)

            source = MS_10_ROOT / 'results/figures/17_prs311_vs_age_at_diagnosis_hexbin.jpg'

            target = MS_10_FIG / source.name

            shutil.copy2(source, target)

            digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()

            assert digest(source) == digest(target)

            ids = pd.read_csv(MS_10_BASE / 'diagnosis_age_ids.csv.gz').person_id

            d = (
                pd.read_csv(MS_10_BASE / 'section9/frozen_model_inputs.csv.gz')
                .set_index('person_id')
                .loc[ids]
            )

            s = (
                pd.read_csv(MS_10_BASE / 'section8/comparable_scores_case_control.csv.gz')
                .set_index('sample_id')
                .loc[ids]
            )

            assert len(d) == 7996 and d.index.is_unique and d.is_case.eq(1).all()

            assert d.age_at_diagnosis.notna().all()

            carrier = s.carried_gene_count_C.gt(0)

            assert carrier.equals(s.B_raw_C.gt(0))

            labels = ['<35', '35–<45', '45–<50', '50–<55', '55–<60', '60–<65', '65–<70', '≥70']

            groups = pd.cut(
                d.age_at_diagnosis,
                [-np.inf, 35, 45, 50, 55, 60, 65, 70, np.inf],
                labels=labels,
                right=False,
            )

            summary = (
                pd.DataFrame({'age_group': groups, 'carrier': carrier})
                .groupby('age_group', observed=False)
                .carrier.agg(n_cases='size', n_carriers='sum')
            )

            assert summary.n_cases.sum() == 7996 and summary.n_carriers.sum() == carrier.sum()

            assert (summary.n_carriers >= 20).all() and (
                summary.n_cases - summary.n_carriers >= 20
            ).all()

            lo, hi = proportion_confint(
                summary.n_carriers.to_numpy(),
                summary.n_cases.to_numpy(),
                alpha=0.05,
                method='wilson',
            )

            summary['carrier_rate_pct'] = summary.n_carriers / summary.n_cases * 100

            summary['wilson95_low_pct'] = lo * 100

            summary['wilson95_high_pct'] = hi * 100

            plt.style.use('default')

            fig, ax = plt.subplots(figsize=(9, 5))

            x = np.arange(len(labels))

            rate = summary.carrier_rate_pct.to_numpy()

            ax.bar(
                x,
                rate,
                width=0.65,
                color='white',
                edgecolor='black',
                yerr=np.vstack([rate - lo * 100, hi * 100 - rate]),
                capsize=4,
            )

            ax.set_xticks(
                x, [f'{label}\n(n={n:,})' for label, n in zip(labels, summary.n_cases)], fontsize=9
            )

            ax.set_xlabel('Age at diagnosis (years)')

            ax.set_ylabel('Carrier rate (%)')

            ax.set_title(
                '(c) 13-gene carrier frequency by age at diagnosis\nAge-at-diagnosis cohort (n=7,996)',
                fontsize=12,
            )

            ax.set_ylim(bottom=0)

            ax.grid(False)

            ax.spines[['top', 'right']].set_visible(False)

            fig.text(
                0.5,
                0.015,
                'Bars: carrier proportion within each age group. Error bars: 95% Wilson intervals.',
                ha='center',
                fontsize=9,
            )

            fig.tight_layout(rect=[0, 0.045, 1, 1])

            stem = MS_10_FIG / '13gene_carrier_frequency_by_age_at_diagnosis'

            sizes = {}

            for ext in ['svg', 'png', 'jpg']:
                kwargs = {'pil_kwargs': {'quality': 85, 'subsampling': 2}} if ext == 'jpg' else {}
                fig.savefig(stem.with_suffix('.' + ext), dpi=180, bbox_inches='tight', **kwargs)
                sizes[ext] = stem.with_suffix('.' + ext).stat().st_size

            smallest = min(sizes, key=sizes.get)

            for ext in ['svg', 'jpg']:
                if ext != smallest:
                    stem.with_suffix('.' + ext).unlink()

            plt.close(fig)

            summary.to_csv(
                MS_10_ROOT / 'results/table/55_13gene_carrier_frequency_by_diagnosis_age.csv'
            )

            caption = 'Rare and common risk and age at diagnosis among cases (age-at-diagnosis cohort, n=7,996).\n\n(a) Age at diagnosis by 13-gene burden group. Carriers are divided at the carrier burden median (1.73695, displayed as 1.74); non-carriers have zero burden. Boxes show medians and interquartile ranges; whiskers extend to observations within 1.5 times the interquartile range. Outliers are not shown.\n\n(b) PRS311 against age at diagnosis. Hexagonal bins with fewer than 20 cases are not shown. The red line is an unadjusted linear fit. This panel is reused unchanged from the original analysis.\n\n(c) Carrier frequency by age at diagnosis, with 95% Wilson intervals. Bars show the proportion of cases carrying a qualifying variant in the 13-gene panel within each age group. Intervals are left-inclusive and right-exclusive; the first group is <35 and the final group is ≥70 years. Counts below each group are total cases, including carriers and non-carriers.\n\nFor (a) and (c), the current 13-gene definition includes PTVs plus qualifying P/LP missense variants; raw burden sums ln(OR) once per carried gene. Age refers to the first qualifying recorded EHR diagnosis, rather than necessarily biological onset. The current variant definition remains provisional.\n'

            (MS_10_FIG / 'rare_common_risk_diagnosis_age_captions.md').write_text(caption)

            report = {
                'n_cases': len(d),
                'n_carriers': int(carrier.sum()),
                'age_intervals': labels,
                'boundary_rule': 'left-inclusive, right-exclusive',
                'plot': 'bar with 95% Wilson intervals',
                'all_carrier_and_noncarrier_cells_at_least_20': True,
                'PRS_panel_reused': str(source.relative_to(MS_10_ROOT)),
                'PRS_panel_sha256': digest(target),
                'smallest_format': smallest,
                'file_sizes': sizes,
            }

            (MS_10_ROOT / 'results/report/55_age_carrier_frequency_validation.json').write_text(
                json.dumps(report, indent=2) + '\n'
            )

            print(summary.to_string())

        # ============================================================================
        # 9. Plot primary sensitivity estimates
        # ============================================================================
        if selected_task == 'sensitivity-figure':
            source = MS_12_ROOT / 'results/table/45_rowlands_fixed_cohort_coefficients.csv'

            d = pd.read_csv(source)

            terms = ['burden_z', 'PRS311_z']

            selected = d[
                d.version.isin([v for v, _ in SENSITIVITY_FIGURE_VERSIONS]) & d.term.isin(terms)
            ].copy()

            assert len(selected) == len(SENSITIVITY_FIGURE_VERSIONS) * 2 * 2

            assert not selected.duplicated(['cohort', 'version', 'term']).any()

            plt.style.use('default')

            fig, axes = plt.subplots(1, 2, figsize=(13, 6.2), sharey=True)

            y = np.arange(len(SENSITIVITY_FIGURE_VERSIONS))

            records = []

            for ax, cohort, title in zip(
                axes,
                ['case_control', 'diagnosis_age'],
                ['(A) Case status (n=141,204)', '(B) Age at diagnosis (n=7,996)'],
            ):
                for term, label, color, marker, offset in [
                    ('burden_z', 'Rare-variant burden', 'tab:orange', 'o', -0.12),
                    ('PRS311_z', 'PRS311', 'tab:green', 's', 0.12),
                ]:
                    sub = (
                        selected[(selected.cohort == cohort) & selected.term.eq(term)]
                        .set_index('version')
                        .loc[[v for v, _ in SENSITIVITY_FIGURE_VERSIONS]]
                    )
                    cols = (
                        ['OR', 'CI_low', 'CI_high']
                        if cohort == 'case_control'
                        else ['beta', 'CI_beta_low', 'CI_beta_high']
                    )
                    est, lo, hi = [sub[c].to_numpy() for c in cols]
                    assert (
                        np.isfinite([est, lo, hi]).all() and (lo <= est).all() and (est <= hi).all()
                    )
                    ax.errorbar(
                        est,
                        y + offset,
                        xerr=np.vstack([est - lo, hi - est]),
                        fmt=marker,
                        color=color,
                        capsize=3,
                        markersize=5,
                        linewidth=1.3,
                        label=label,
                    )
                    for version, e, l, h in zip(sub.index, est, lo, hi):
                        records.append(
                            {
                                'cohort': cohort,
                                'version': version,
                                'term': term,
                                'estimate': e,
                                'CI_low': l,
                                'CI_high': h,
                            }
                        )
                ax.set_title(title, fontsize=12)
                ax.grid(False)
                for spine in ax.spines.values():
                    spine.set_visible(True)
                ax.axvline(
                    1 if cohort == 'case_control' else 0, color='.5', linestyle='--', linewidth=0.8
                )
                ax.set_xlabel(
                    'Odds ratio per SD'
                    if cohort == 'case_control'
                    else 'Difference in age at diagnosis (years per SD)'
                )

            axes[0].set_xscale('log')

            from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

            axes[0].xaxis.set_major_locator(FixedLocator([1, 1.2, 1.4, 1.6]))

            axes[0].xaxis.set_major_formatter(FuncFormatter(lambda value, _: f'{value:g}'))

            axes[0].xaxis.set_minor_locator(NullLocator())

            axes[0].set_xlim(0.98, 1.65)

            axes[1].set_xlim(-2.25, 0.12)

            axes[0].set_yticks(y, [label for _, label in SENSITIVITY_FIGURE_VERSIONS], fontsize=10)

            axes[0].invert_yaxis()

            handles, labels = axes[0].get_legend_handles_labels()

            fig.legend(
                handles,
                labels,
                frameon=False,
                loc='upper center',
                bbox_to_anchor=(0.6, 0.93),
                ncol=2,
                fontsize=10,
            )

            fig.suptitle('Figure 6. Primary and sensitivity estimates', fontsize=14)

            fig.text(
                0.5,
                0.015,
                'Horizontal lines: 95% CIs. Each burden is scaled by its own cohort/model SD.\nCase-status burden OR is conditional on mean PRS311; the age model has no interaction.',
                ha='center',
                fontsize=9,
            )

            fig.tight_layout(rect=[0, 0.07, 1, 0.87])

            output = (
                MS_12_ROOT / 'manuscript/figures/Figure_06_primary_and_sensitivity_estimates.svg'
            )

            fig.savefig(output, bbox_inches='tight')

            preview = (
                MS_12_ROOT
                / 'data/manuscript_refresh_20261006/figure_previews/Figure_06_primary_and_sensitivity_estimates.png'
            )

            preview.parent.mkdir(parents=True, exist_ok=True)

            fig.savefig(preview, dpi=160, bbox_inches='tight')

            plt.close(fig)

            pd.DataFrame(records).to_csv(
                MS_12_ROOT / 'results/table/56_figure5_primary_sensitivity_estimates.csv',
                index=False,
            )

            caption = 'Figure 6. Primary and sensitivity estimates. (A) Case-status odds ratios per SD from the joint logistic model, including PRS311, rare-variant burden, their interaction, age at last EHR observation and PC1–10 (n=141,204). (B) Difference in recorded age at diagnosis per SD among cases from an OLS model including PRS311, burden and PC1–10, without interaction (n=7,996). Points show estimates for burden (orange circles) and PRS311 (green squares); horizontal lines show 95% CIs based on HC1 robust standard errors.\n\nThe primary analysis uses the 13-gene PTV plus qualifying P/LP missense burden. Sensitivities exclude STK11 (12 genes), restrict the panel to the nine original genes with the same current PTV/missense rules and Rowlands weights, include PTVs only, apply stricter heterozygous VAF (0.30–0.70) to TP53/NF1, relax the ClinVar review threshold to ≥1 star, or exclude CDH1 (12 genes). The original nine-gene PTV analysis using its original definitions and weights is shown as a historical comparison baseline. It differs from the current nine-gene PTV plus missense analysis in eligibility/QC and weights as well as variant type.\n\nThe same fixed participants and PRS311 are used in each model. Burden is scaled by its own SD within each model cohort, so differences in per-SD coefficients do not isolate effects of changing variant definitions. In panel A, the burden OR is conditional on PRS311 at its mean; it is not a marginal carrier OR. The current variant definition remains provisional.\n'

            (MS_12_ROOT / 'manuscript/figure_captions/Figure_06_caption.md').write_text(caption)

            report = {
                'source': str(source.relative_to(MS_12_ROOT)),
                'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
                'versions': [v for v, _ in SENSITIVITY_FIGURE_VERSIONS],
                'plotted_estimates': len(records),
                'models_refit': False,
                'output': str(output.relative_to(MS_12_ROOT)),
            }

            (MS_12_ROOT / 'results/report/56_figure5_validation.json').write_text(
                json.dumps(report, indent=2) + '\n'
            )

            print(json.dumps(report, indent=2))

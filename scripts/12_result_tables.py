"""Construct cohort descriptions, carrier attribution and model result tables.

Outputs remain under results/table and results/report.
The existing summary operation also writes its associated model-comparison figure.
"""

from __future__ import annotations

from pathlib import Path
import sys

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.functions.paths import PROJECT_ROOT
from scripts.functions.cli import parse_tasks

TASKS = ('descriptives', 'cohort-table', 'attribution', 'summary')
if __name__ == '__main__':
    selected_tasks = parse_tasks(
        __doc__, TASKS, ['descriptives', 'cohort-table', 'attribution', 'summary']
    )

import os


import matplotlib

import matplotlib.pyplot as plt

import numpy as np

import pandas as pd

import seaborn as sns


import json

import yaml

CARRIER_ATTRIBUTION_ROOT = PROJECT_ROOT

CARRIER_ATTRIBUTION_BASE = CARRIER_ATTRIBUTION_ROOT / 'data/29_rowlands_extension'

CARRIER_ATTRIBUTION_OUT = CARRIER_ATTRIBUTION_BASE / 'section11_models'

import ast

matplotlib.use('Agg')

from scipy import stats

RESULT_SUMMARY_ROOT = PROJECT_ROOT

RESULT_SUMMARY_BASE = RESULT_SUMMARY_ROOT / 'data/29_rowlands_extension'

RESULT_SUMMARY_OUT = RESULT_SUMMARY_BASE / 'section11_models'


def cohort_descriptives_log(msg):
    print(msg, flush=True)


def summarize(sub, label):
    rows = []
    rows.append(
        {
            'group': label,
            'variable': 'n_total',
            'n_nonmissing': len(sub),
            'mean': np.nan,
            'median': np.nan,
            'sd': np.nan,
            'iqr_low': np.nan,
            'iqr_high': np.nan,
            'min': np.nan,
            'max': np.nan,
        }
    )
    for v in COHORT_DESCRIPTIVES_CONT_VARS:
        s = sub[v].dropna()
        if len(s) == 0:
            rows.append(
                {
                    'group': label,
                    'variable': v,
                    'n_nonmissing': 0,
                    'mean': np.nan,
                    'median': np.nan,
                    'sd': np.nan,
                    'iqr_low': np.nan,
                    'iqr_high': np.nan,
                    'min': np.nan,
                    'max': np.nan,
                }
            )
            continue
        rows.append(
            {
                'group': label,
                'variable': v,
                'n_nonmissing': len(s),
                'mean': s.mean(),
                'median': s.median(),
                'sd': s.std(),
                'iqr_low': s.quantile(0.25),
                'iqr_high': s.quantile(0.75),
                'min': s.min(),
                'max': s.max(),
            }
        )
    return rows


def fmt_cat(count, is_masked, denom):
    s = cohort_table_render(count, is_masked, denom)
    if s.endswith('%'):
        n, pct = s.split(' ')
        return f'{int(n):,d} ({pct[:-1]})'
    return s


def fmt_cont(x, kind):
    x = x.dropna()
    if kind == 'mean':
        return f'{x.mean():.1f} ({x.std():.1f})'
    return f'{x.median():.1f} [{x.quantile(0.25):.1f}-{x.quantile(0.75):.1f}]'


if __name__ == "__main__":
    for selected_task in selected_tasks:

        # ============================================================================
        # 1. Describe cohort
        # ============================================================================
        if selected_task == 'descriptives':
            matplotlib.use('Agg')

            plt.style.use('default')

            COHORT_DESCRIPTIVES_BASE = str(PROJECT_ROOT)

            COHORT_DESCRIPTIVES_DATA_DIR = os.path.join(COHORT_DESCRIPTIVES_BASE, 'data')

            COHORT_DESCRIPTIVES_FIG_DIR = COHORT_DESCRIPTIVES_TABLE_DIR = os.path.join(
                COHORT_DESCRIPTIVES_DATA_DIR, '08_internal'
            )

            os.makedirs(COHORT_DESCRIPTIVES_FIG_DIR, exist_ok=True)

            COHORT_DESCRIPTIVES_CASE_COLOR = 'tab:blue'

            COHORT_DESCRIPTIVES_CONTROL_COLOR = 'tab:orange'

            COHORT_DESCRIPTIVES_DPI = 150

            COHORT_DESCRIPTIVES_GROUP_COLORS = {
                'case': COHORT_DESCRIPTIVES_CASE_COLOR,
                'control': COHORT_DESCRIPTIVES_CONTROL_COLOR,
            }

            COHORT_DESCRIPTIVES_GROUP_ORDER = ['case', 'control']

            cohort_descriptives_log('=' * 80)

            cohort_descriptives_log('Step A: Read three output tables and join')

            cohort_descriptives_log('=' * 80)

            pheno = pd.read_csv(
                os.path.join(
                    COHORT_DESCRIPTIVES_DATA_DIR, '04_phenotype_diagnostics_person_level.csv.gz'
                )
            )

            cohort_descriptives_qc = pd.read_csv(
                os.path.join(COHORT_DESCRIPTIVES_DATA_DIR, '05_cohort_qc_flags.csv')
            )

            cohort_descriptives_prs = pd.read_csv(
                os.path.join(COHORT_DESCRIPTIVES_DATA_DIR, 'prs_scores.csv')
            )

            cohort_descriptives_log(
                f'phenotype table:{pheno.shape}, QC form:{cohort_descriptives_qc.shape}, PRS table:{cohort_descriptives_prs.shape}'
            )

            cohort_descriptives_qc = cohort_descriptives_qc[cohort_descriptives_qc['is_eur'] == 1]

            set_qc = set(cohort_descriptives_qc['person_id'])

            set_prs = set(cohort_descriptives_prs['person_id'])

            assert set_qc == set_prs, 'EUR list is inconsistent with PRS table person_id collection'

            assert set_qc <= set(
                pheno['person_id']
            ), 'There is someone on the EUR list who is not in list 04'

            qc_cols = ['person_id', 'flag_related_exclude', 'flag_high_missing', 'pass_qc']

            cohort_descriptives_df = pheno.merge(
                cohort_descriptives_qc[qc_cols], on='person_id', how='inner'
            )

            cohort_descriptives_df = cohort_descriptives_df.merge(
                cohort_descriptives_prs, on='person_id', how='inner'
            )

            cohort_descriptives_log(
                f'Wide table after join (Genetics EUR):{cohort_descriptives_df.shape}'
            )

            assert cohort_descriptives_df.shape[0] == len(
                set_qc
            ), 'The number of rows after join is not as expected, please check'

            cohort_descriptives_df['ehr_span_years'] = (
                cohort_descriptives_df['ehr_span_days'] / 365.25
            )

            cohort_descriptives_log('group distribution:')

            cohort_descriptives_log(str(cohort_descriptives_df['group'].value_counts(dropna=False)))

            cohort_descriptives_log('pass_qc distribution (by group):')

            cohort_descriptives_log(
                str(pd.crosstab(cohort_descriptives_df['group'], cohort_descriptives_df['pass_qc']))
            )

            cohort_descriptives_log('\n' + '=' * 80)

            cohort_descriptives_log('Step B: Queue overview summary statistics table')

            cohort_descriptives_log('=' * 80)

            COHORT_DESCRIPTIVES_CONT_VARS = [
                'age_at_diagnosis',
                'age_at_enrollment',
                'age_at_last_obs',
                'ehr_span_years',
                'n_distinct_ehr_dates',
                'dx_code_occurrence_count',
                'PRS313',
                'PRS311',
            ]

            summary_rows = []

            for cohort_descriptives_g in ['case', 'control', 'other_cancer_excluded']:
                summary_rows.extend(
                    summarize(
                        cohort_descriptives_df[
                            cohort_descriptives_df['group'] == cohort_descriptives_g
                        ],
                        cohort_descriptives_g,
                    )
                )

            summary_rows.extend(summarize(cohort_descriptives_df, 'all_174381'))

            summary_df = pd.DataFrame(summary_rows)

            for cohort_descriptives_col in [
                'mean',
                'median',
                'sd',
                'iqr_low',
                'iqr_high',
                'min',
                'max',
            ]:
                summary_df[cohort_descriptives_col] = summary_df[cohort_descriptives_col].round(3)

            summary_path = os.path.join(
                COHORT_DESCRIPTIVES_TABLE_DIR, '08_cohort_summary_statistics.csv'
            )

            summary_df.to_csv(summary_path, index=False)

            cohort_descriptives_log(f'Has written:{summary_path} ({summary_df.shape[0]}OK)')

            cohort_descriptives_log(summary_df.to_string(index=False))

            qc_overview = (
                cohort_descriptives_df[cohort_descriptives_df['group'].isin(['case', 'control'])]
                .groupby(['group', 'pass_qc'])
                .size()
                .reset_index(name='n')
            )

            qc_overview_path = os.path.join(
                COHORT_DESCRIPTIVES_TABLE_DIR, '08_cohort_qc_before_after_n.csv'
            )

            qc_overview.to_csv(qc_overview_path, index=False)

            cohort_descriptives_log(f'\n has written:{qc_overview_path}')

            cohort_descriptives_log(qc_overview.to_string(index=False))

            cohort_descriptives_log('\n' + '=' * 80)

            cohort_descriptives_log('Step C: Age distribution chart (histogram+KDE, 3-panel)')

            cohort_descriptives_log('=' * 80)

            cohort_descriptives_case_df = cohort_descriptives_df[
                cohort_descriptives_df['group'] == 'case'
            ]

            cohort_descriptives_control_df = cohort_descriptives_df[
                cohort_descriptives_df['group'] == 'control'
            ]

            cohort_descriptives_fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

            cohort_descriptives_ax = axes[0]

            sns.histplot(
                cohort_descriptives_case_df['age_at_diagnosis'].dropna(),
                stat='density',
                color=COHORT_DESCRIPTIVES_CASE_COLOR,
                alpha=0.5,
                label='case',
                ax=cohort_descriptives_ax,
                bins=30,
            )

            sns.kdeplot(
                cohort_descriptives_case_df['age_at_diagnosis'].dropna(),
                color=COHORT_DESCRIPTIVES_CASE_COLOR,
                ax=cohort_descriptives_ax,
            )

            cohort_descriptives_ax.set_title('age_at_diagnosis (case only)')

            cohort_descriptives_ax.set_xlabel('age at diagnosis (years)')

            cohort_descriptives_ax.set_ylabel('density')

            cohort_descriptives_ax.legend()

            cohort_descriptives_ax.grid(False)

            cohort_descriptives_ax = axes[1]

            sns.histplot(
                cohort_descriptives_case_df['age_at_enrollment'].dropna(),
                stat='density',
                color=COHORT_DESCRIPTIVES_CASE_COLOR,
                alpha=0.5,
                label='case',
                ax=cohort_descriptives_ax,
                bins=30,
            )

            sns.kdeplot(
                cohort_descriptives_case_df['age_at_enrollment'].dropna(),
                color=COHORT_DESCRIPTIVES_CASE_COLOR,
                ax=cohort_descriptives_ax,
            )

            sns.histplot(
                cohort_descriptives_control_df['age_at_enrollment'].dropna(),
                stat='density',
                color=COHORT_DESCRIPTIVES_CONTROL_COLOR,
                alpha=0.5,
                label='control',
                ax=cohort_descriptives_ax,
                bins=30,
            )

            sns.kdeplot(
                cohort_descriptives_control_df['age_at_enrollment'].dropna(),
                color=COHORT_DESCRIPTIVES_CONTROL_COLOR,
                ax=cohort_descriptives_ax,
            )

            cohort_descriptives_ax.set_title('age_at_enrollment (case vs control)')

            cohort_descriptives_ax.set_xlabel('age at enrollment (years)')

            cohort_descriptives_ax.set_ylabel('density')

            cohort_descriptives_ax.legend()

            cohort_descriptives_ax.grid(False)

            cohort_descriptives_ax = axes[2]

            sns.histplot(
                cohort_descriptives_case_df['age_at_last_obs'].dropna(),
                stat='density',
                color=COHORT_DESCRIPTIVES_CASE_COLOR,
                alpha=0.5,
                label='case',
                ax=cohort_descriptives_ax,
                bins=30,
            )

            sns.kdeplot(
                cohort_descriptives_case_df['age_at_last_obs'].dropna(),
                color=COHORT_DESCRIPTIVES_CASE_COLOR,
                ax=cohort_descriptives_ax,
            )

            sns.histplot(
                cohort_descriptives_control_df['age_at_last_obs'].dropna(),
                stat='density',
                color=COHORT_DESCRIPTIVES_CONTROL_COLOR,
                alpha=0.5,
                label='control',
                ax=cohort_descriptives_ax,
                bins=30,
            )

            sns.kdeplot(
                cohort_descriptives_control_df['age_at_last_obs'].dropna(),
                color=COHORT_DESCRIPTIVES_CONTROL_COLOR,
                ax=cohort_descriptives_ax,
            )

            cohort_descriptives_ax.set_title('age_at_last_obs (case vs control)')

            cohort_descriptives_ax.set_xlabel('age at last observation (years)')

            cohort_descriptives_ax.set_ylabel('density')

            cohort_descriptives_ax.legend()

            cohort_descriptives_ax.grid(False)

            cohort_descriptives_fig.suptitle(
                'Age distributions: case vs control (All of Us WGS cohort)'
            )

            cohort_descriptives_fig.tight_layout()

            age_hist_path = os.path.join(
                COHORT_DESCRIPTIVES_FIG_DIR, 'age_distribution_histogram_kde.png'
            )

            cohort_descriptives_fig.savefig(age_hist_path, dpi=COHORT_DESCRIPTIVES_DPI)

            plt.close(cohort_descriptives_fig)

            cohort_descriptives_log(f'Has written:{age_hist_path}')

            cohort_descriptives_fig, axes = plt.subplots(1, 3, figsize=(12, 4.5))

            cohort_descriptives_ax = axes[0]

            sns.boxplot(
                y=cohort_descriptives_case_df['age_at_diagnosis'].dropna(),
                ax=cohort_descriptives_ax,
                color=COHORT_DESCRIPTIVES_CASE_COLOR,
            )

            cohort_descriptives_ax.set_xticks([0])

            cohort_descriptives_ax.set_xticklabels(['case'])

            cohort_descriptives_ax.set_title('age_at_diagnosis\n(case only, no control value)')

            cohort_descriptives_ax.set_ylabel('age (years)')

            cohort_descriptives_ax.grid(False)

            cohort_descriptives_ax = axes[1]

            cohort_descriptives_plot_df = cohort_descriptives_df[
                cohort_descriptives_df['group'].isin(['case', 'control'])
            ][['group', 'age_at_enrollment']].dropna()

            sns.boxplot(
                x='group',
                y='age_at_enrollment',
                data=cohort_descriptives_plot_df,
                order=COHORT_DESCRIPTIVES_GROUP_ORDER,
                hue='group',
                palette=COHORT_DESCRIPTIVES_GROUP_COLORS,
                legend=False,
                ax=cohort_descriptives_ax,
            )

            cohort_descriptives_ax.set_title('age_at_enrollment')

            cohort_descriptives_ax.set_xlabel('')

            cohort_descriptives_ax.set_ylabel('age (years)')

            cohort_descriptives_ax.grid(False)

            cohort_descriptives_ax = axes[2]

            cohort_descriptives_plot_df = cohort_descriptives_df[
                cohort_descriptives_df['group'].isin(['case', 'control'])
            ][['group', 'age_at_last_obs']].dropna()

            sns.boxplot(
                x='group',
                y='age_at_last_obs',
                data=cohort_descriptives_plot_df,
                order=COHORT_DESCRIPTIVES_GROUP_ORDER,
                hue='group',
                palette=COHORT_DESCRIPTIVES_GROUP_COLORS,
                legend=False,
                ax=cohort_descriptives_ax,
            )

            cohort_descriptives_ax.set_title('age_at_last_obs')

            cohort_descriptives_ax.set_xlabel('')

            cohort_descriptives_ax.set_ylabel('age (years)')

            cohort_descriptives_ax.grid(False)

            cohort_descriptives_fig.suptitle('Age boxplots: case vs control')

            cohort_descriptives_fig.tight_layout()

            age_box_path = os.path.join(
                COHORT_DESCRIPTIVES_FIG_DIR, 'age_boxplot_case_vs_control.png'
            )

            cohort_descriptives_fig.savefig(age_box_path, dpi=COHORT_DESCRIPTIVES_DPI)

            plt.close(cohort_descriptives_fig)

            cohort_descriptives_log(f'Has written:{age_box_path}')

            cohort_descriptives_log('\n' + '=' * 80)

            cohort_descriptives_log(
                'Step D: EHR coverage intensity boxplot (log1p scale, case vs control)'
            )

            cohort_descriptives_log('=' * 80)

            cohort_descriptives_log(
                'Note: ehr_span_years/n_distinct_ehr_dates/dx_code_occurrence_count are severely right-skewed and contain 0 values.'
            )

            cohort_descriptives_log(
                'Direct use of log scale cannot handle 0, so the original value is transformed by log1p(x)=log(1+x) before drawing the boxplot.'
            )

            cohort_descriptives_log(
                'It not only retains the visual position of the 0 value, but also compresses the long tail to prevent the graph from being flattened by extreme values. The coordinate axes are marked with log1p(x+1) scale.'
            )

            ehr_vars = ['ehr_span_years', 'n_distinct_ehr_dates', 'dx_code_occurrence_count']

            cohort_descriptives_fig, axes = plt.subplots(1, 3, figsize=(12, 4.5))

            for cohort_descriptives_ax, cohort_descriptives_v in zip(axes, ehr_vars):
                cohort_descriptives_plot_df = (
                    cohort_descriptives_df[
                        cohort_descriptives_df['group'].isin(['case', 'control'])
                    ][['group', cohort_descriptives_v]]
                    .dropna()
                    .copy()
                )
                cohort_descriptives_plot_df['log1p_val'] = np.log1p(
                    cohort_descriptives_plot_df[cohort_descriptives_v]
                )
                sns.boxplot(
                    x='group',
                    y='log1p_val',
                    data=cohort_descriptives_plot_df,
                    order=COHORT_DESCRIPTIVES_GROUP_ORDER,
                    hue='group',
                    palette=COHORT_DESCRIPTIVES_GROUP_COLORS,
                    legend=False,
                    ax=cohort_descriptives_ax,
                )
                cohort_descriptives_ax.set_title(cohort_descriptives_v)
                cohort_descriptives_ax.set_xlabel('')
                cohort_descriptives_ax.set_ylabel(f'log1p({cohort_descriptives_v})')
                cohort_descriptives_ax.grid(False)

            cohort_descriptives_fig.suptitle(
                'EHR coverage intensity boxplots: case vs control (log1p-transformed)'
            )

            cohort_descriptives_fig.tight_layout()

            ehr_box_path = os.path.join(
                COHORT_DESCRIPTIVES_FIG_DIR, 'ehr_coverage_boxplot_case_vs_control.png'
            )

            cohort_descriptives_fig.savefig(ehr_box_path, dpi=COHORT_DESCRIPTIVES_DPI)

            plt.close(cohort_descriptives_fig)

            cohort_descriptives_log(f'Has written:{ehr_box_path}')

            cohort_descriptives_fig, axes = plt.subplots(1, 3, figsize=(12, 4.5))

            for cohort_descriptives_ax, cohort_descriptives_v in zip(axes, ehr_vars):
                cohort_descriptives_plot_df = cohort_descriptives_df[
                    cohort_descriptives_df['group'].isin(['case', 'control'])
                ][['group', cohort_descriptives_v]].dropna()
                ylim_top = cohort_descriptives_plot_df[cohort_descriptives_v].quantile(0.99)
                sns.boxplot(
                    x='group',
                    y=cohort_descriptives_v,
                    data=cohort_descriptives_plot_df,
                    order=COHORT_DESCRIPTIVES_GROUP_ORDER,
                    hue='group',
                    palette=COHORT_DESCRIPTIVES_GROUP_COLORS,
                    legend=False,
                    ax=cohort_descriptives_ax,
                    showfliers=False,
                )
                cohort_descriptives_ax.set_ylim(bottom=-0.02 * ylim_top, top=ylim_top)
                cohort_descriptives_ax.set_title(
                    f'{cohort_descriptives_v}\n(y capped at 99th pct, outlier points hidden)'
                )
                cohort_descriptives_ax.set_xlabel('')
                cohort_descriptives_ax.grid(False)

            cohort_descriptives_fig.suptitle(
                'EHR coverage intensity boxplots: case vs control (raw scale, 99th pct cap)'
            )

            cohort_descriptives_fig.tight_layout()

            ehr_box_raw_path = os.path.join(
                COHORT_DESCRIPTIVES_FIG_DIR,
                'ehr_coverage_boxplot_case_vs_control_rawscale_capped.png',
            )

            cohort_descriptives_fig.savefig(ehr_box_raw_path, dpi=COHORT_DESCRIPTIVES_DPI)

            plt.close(cohort_descriptives_fig)

            cohort_descriptives_log(f'Has written:{ehr_box_raw_path}')

            cohort_descriptives_log('\n' + '=' * 80)

            cohort_descriptives_log('Step E: PRS313/PRS311 distribution diagram')

            cohort_descriptives_log('=' * 80)

            for prs_col, cohort_descriptives_fname in [
                ('PRS313', 'prs313_distribution_case_vs_control.png'),
                ('PRS311', 'prs311_distribution_case_vs_control.png'),
            ]:
                cohort_descriptives_fig, cohort_descriptives_ax = plt.subplots(figsize=(6, 4.5))
                sns.histplot(
                    cohort_descriptives_case_df[prs_col].dropna(),
                    stat='density',
                    color=COHORT_DESCRIPTIVES_CASE_COLOR,
                    alpha=0.5,
                    label='case',
                    ax=cohort_descriptives_ax,
                    bins=40,
                )
                sns.kdeplot(
                    cohort_descriptives_case_df[prs_col].dropna(),
                    color=COHORT_DESCRIPTIVES_CASE_COLOR,
                    ax=cohort_descriptives_ax,
                )
                sns.histplot(
                    cohort_descriptives_control_df[prs_col].dropna(),
                    stat='density',
                    color=COHORT_DESCRIPTIVES_CONTROL_COLOR,
                    alpha=0.5,
                    label='control',
                    ax=cohort_descriptives_ax,
                    bins=40,
                )
                sns.kdeplot(
                    cohort_descriptives_control_df[prs_col].dropna(),
                    color=COHORT_DESCRIPTIVES_CONTROL_COLOR,
                    ax=cohort_descriptives_ax,
                )
                cohort_descriptives_ax.axvline(
                    cohort_descriptives_case_df[prs_col].mean(),
                    color=COHORT_DESCRIPTIVES_CASE_COLOR,
                    linestyle='--',
                    linewidth=1,
                )
                cohort_descriptives_ax.axvline(
                    cohort_descriptives_control_df[prs_col].mean(),
                    color=COHORT_DESCRIPTIVES_CONTROL_COLOR,
                    linestyle='--',
                    linewidth=1,
                )
                cohort_descriptives_ax.set_title(
                    f'{prs_col} distribution: case vs control\n(dashed lines = group means)'
                )
                cohort_descriptives_ax.set_xlabel(prs_col)
                cohort_descriptives_ax.set_ylabel('density')
                cohort_descriptives_ax.legend()
                cohort_descriptives_ax.grid(False)
                cohort_descriptives_fig.tight_layout()
                cohort_descriptives_out_path = os.path.join(
                    COHORT_DESCRIPTIVES_FIG_DIR, cohort_descriptives_fname
                )
                cohort_descriptives_fig.savefig(
                    cohort_descriptives_out_path, dpi=COHORT_DESCRIPTIVES_DPI
                )
                plt.close(cohort_descriptives_fig)
                cohort_descriptives_log(f'Has written:{cohort_descriptives_out_path}')
                cohort_descriptives_log(
                    f'{prs_col}: case mean={cohort_descriptives_case_df[prs_col].mean():.3f} median={cohort_descriptives_case_df[prs_col].median():.3f}; control mean={cohort_descriptives_control_df[prs_col].mean():.3f} median={cohort_descriptives_control_df[prs_col].median():.3f}; mean diff (case-control)={cohort_descriptives_case_df[prs_col].mean() - cohort_descriptives_control_df[prs_col].mean():.3f}'
                )

            cohort_descriptives_fig, axes = plt.subplots(
                1, 2, figsize=(11, 4.5), sharex=True, sharey=True
            )

            cohort_descriptives_ax = axes[0]

            sns.kdeplot(
                cohort_descriptives_case_df['PRS313'].dropna(),
                label='PRS313',
                ax=cohort_descriptives_ax,
            )

            sns.kdeplot(
                cohort_descriptives_case_df['PRS311'].dropna(),
                label='PRS311',
                ax=cohort_descriptives_ax,
            )

            cohort_descriptives_ax.set_title('case')

            cohort_descriptives_ax.set_xlabel('PRS score')

            cohort_descriptives_ax.set_ylabel('density')

            cohort_descriptives_ax.legend()

            cohort_descriptives_ax.grid(False)

            cohort_descriptives_ax = axes[1]

            sns.kdeplot(
                cohort_descriptives_control_df['PRS313'].dropna(),
                label='PRS313',
                ax=cohort_descriptives_ax,
            )

            sns.kdeplot(
                cohort_descriptives_control_df['PRS311'].dropna(),
                label='PRS311',
                ax=cohort_descriptives_ax,
            )

            cohort_descriptives_ax.set_title('control')

            cohort_descriptives_ax.set_xlabel('PRS score')

            cohort_descriptives_ax.legend()

            cohort_descriptives_ax.grid(False)

            cohort_descriptives_fig.suptitle(
                'PRS313 vs PRS311: distribution comparison within case / control'
            )

            cohort_descriptives_fig.tight_layout()

            prs_compare_path = os.path.join(
                COHORT_DESCRIPTIVES_FIG_DIR, 'prs313_vs_prs311_comparison.png'
            )

            cohort_descriptives_fig.savefig(prs_compare_path, dpi=COHORT_DESCRIPTIVES_DPI)

            plt.close(cohort_descriptives_fig)

            cohort_descriptives_log(f'Has written:{prs_compare_path}')

            cohort_descriptives_log('\n' + '=' * 80)

            cohort_descriptives_log('Step F: Output file size check (required < 5MB)')

            cohort_descriptives_log('=' * 80)

            for cohort_descriptives_d in [
                COHORT_DESCRIPTIVES_FIG_DIR,
                COHORT_DESCRIPTIVES_TABLE_DIR,
            ]:
                for cohort_descriptives_f in sorted(os.listdir(cohort_descriptives_d)):
                    if (
                        not cohort_descriptives_f.startswith('08_')
                        and 'case_vs_control' not in cohort_descriptives_f
                        and ('distribution' not in cohort_descriptives_f)
                        and ('comparison' not in cohort_descriptives_f)
                    ):
                        continue
                    cohort_descriptives_path = os.path.join(
                        cohort_descriptives_d, cohort_descriptives_f
                    )
                    if os.path.isfile(cohort_descriptives_path):
                        cohort_descriptives_size_mb = (
                            os.path.getsize(cohort_descriptives_path) / 1000000.0
                        )
                        cohort_descriptives_flag = (
                            'OK' if cohort_descriptives_size_mb < 5 else '!!! More than 5MB !!!'
                        )
                        cohort_descriptives_log(
                            f'{cohort_descriptives_path}: {cohort_descriptives_size_mb:.3f} MB [{cohort_descriptives_flag}]'
                        )

            cohort_descriptives_log('Finish.')

        # ============================================================================
        # 2. Write cohort table
        # ============================================================================
        if selected_task == 'cohort-table':
            from google.cloud import bigquery as cohort_table_bigquery

            sys.path.insert(0, str(PROJECT_ROOT / 'scripts/functions'))

            from scripts.functions.disclosure import (
                decide as cohort_table_decide,
                render as cohort_table_render,
            )

            COHORT_TABLE_BASE = str(PROJECT_ROOT)

            COHORT_TABLE_DATA_DIR = f'{COHORT_TABLE_BASE}/data'

            COHORT_TABLE_TABLE_DIR = f'{COHORT_TABLE_BASE}/results/table'

            COHORT_TABLE_PROJECT_ID = os.environ['GOOGLE_PROJECT']

            COHORT_TABLE_COLS = ['aim1_case', 'aim1_control', 'aim4_case']

            HIDDEN = 'aim1_not_aim4'

            COL_LABEL = {
                'aim1_case': 'Aim 1 cases',
                'aim1_control': 'Aim 1 controls',
                'aim4_case': 'Aim 4 cases',
            }

            cohort_table_qc = pd.read_csv(f'{COHORT_TABLE_DATA_DIR}/05_cohort_qc_flags.csv')

            cohort_table_qc = cohort_table_qc[cohort_table_qc['pass_qc'] == 1][
                ['person_id', 'analysis_group']
            ]

            ph = pd.read_csv(
                f'{COHORT_TABLE_DATA_DIR}/04_phenotype_diagnostics_person_level.csv.gz',
                usecols=[
                    'person_id',
                    'age_at_diagnosis',
                    'age_at_enrollment',
                    'age_at_last_obs',
                    'bc_is_first_record',
                ],
            )

            cl = pd.read_csv(
                f'{COHORT_TABLE_DATA_DIR}/09_extended_clinical_person_level.csv.gz'
            ).drop(columns=['group'])

            cohort_table_prs = pd.read_csv(
                f'{COHORT_TABLE_DATA_DIR}/prs_scores.csv', usecols=['person_id', 'PRS311']
            )

            race = (
                cohort_table_bigquery.Client()
                .query(
                    f'\n        SELECT person_id, race_concept_id\n        FROM `{COHORT_TABLE_PROJECT_ID}.my_analysis_results.female_wgs_ancestry_cohort`\n    '
                )
                .to_dataframe()
            )

            cohort_table_df = (
                cohort_table_qc.merge(ph, on='person_id', how='left', validate='one_to_one')
                .merge(cl, on='person_id', how='left', validate='one_to_one')
                .merge(cohort_table_prs, on='person_id', how='left', validate='one_to_one')
                .merge(race, on='person_id', how='left', validate='one_to_one')
            )

            assert len(cohort_table_df) == len(cohort_table_qc)

            is_case = cohort_table_df['analysis_group'] == 'case'

            in_aim4 = is_case & (cohort_table_df['bc_is_first_record'] != True)

            masks = {
                'aim1_case': is_case,
                'aim1_control': cohort_table_df['analysis_group'] == 'control',
                'aim4_case': in_aim4,
                HIDDEN: is_case & ~in_aim4,
            }

            n_col = {cohort_table_c: int(m.sum()) for cohort_table_c, m in masks.items()}

            assert n_col['aim1_case'] == n_col['aim4_case'] + n_col[HIDDEN]

            cohort_table_df['menopause'] = cohort_table_df['menopause_status'].map(
                {
                    'postmenopausal': 'Postmenopausal',
                    'premenopausal': 'Premenopausal',
                    'uncertain_hormone_influenced': 'Uncertain or unknown',
                    'unknown_not_sure': 'Uncertain or unknown',
                }
            )

            cohort_table_df['smoking'] = cohort_table_df['smoking_status'].map(
                {
                    'never_smoker': 'Never',
                    'former_smoker': 'Former or unknown currency',
                    'ever_smoker_currency_unknown': 'Former or unknown currency',
                    'current_smoker': 'Current',
                }
            )

            cohort_table_df['family_history'] = cohort_table_df['family_history_breast_cancer'].map(
                {
                    'positive_nonself': 'Positive (first- or second-degree relative)',
                    'negative_or_self_only': 'Negative',
                }
            )

            cohort_table_df['seq_site'] = cohort_table_df['seq_site_id'].map(
                {
                    'bi': 'Broad Institute',
                    'uw': 'University of Washington',
                    'bcm': 'Baylor College of Medicine',
                }
            )

            cohort_table_df['self_reported_race'] = np.where(
                cohort_table_df['race_concept_id'] == 8527,
                'White',
                'Other, multiple, or not reported',
            )

            CATEGORICAL = [
                (
                    'Menopausal status',
                    'menopause',
                    ['Premenopausal', 'Postmenopausal', 'Uncertain or unknown'],
                ),
                (
                    'Family history of breast cancer',
                    'family_history',
                    ['Positive (first- or second-degree relative)', 'Negative'],
                ),
                ('Smoking status', 'smoking', ['Never', 'Former or unknown currency', 'Current']),
                (
                    'Self-reported race',
                    'self_reported_race',
                    ['White', 'Other, multiple, or not reported'],
                ),
                (
                    'Sequencing center',
                    'seq_site',
                    ['Broad Institute', 'University of Washington', 'Baylor College of Medicine'],
                ),
                ('Sequencer', 'seq_sequencer', ['NovaSeq6000 S4', 'NovaSeq X 25B']),
            ]

            cells, var_of, published = ({}, {}, {})

            for cohort_table_var, cohort_table_col, levels in CATEGORICAL:
                unexpected = set(cohort_table_df[cohort_table_col].dropna()) - set(levels)
                assert not unexpected, f'{cohort_table_col}There are unexpected values:{unexpected}'
                for lv in levels + ['Missing']:
                    key = f'{cohort_table_var}: {lv}'
                    hit = (
                        cohort_table_df[cohort_table_col].isna()
                        if lv == 'Missing'
                        else cohort_table_df[cohort_table_col].eq(lv)
                    )
                    cells[key] = {
                        cohort_table_c: int((m & hit).sum()) for cohort_table_c, m in masks.items()
                    }
                    var_of[key] = cohort_table_var
                    published[key] = set(COHORT_TABLE_COLS) | {HIDDEN}

            masked = cohort_table_decide(
                cells, var_of, published, column_partitions=[('aim1_case', 'aim4_case', HIDDEN)]
            )

            n_masked = sum((len(v - {HIDDEN}) for v in masked.values()))

            print(
                f'decide(): print the number of shaded cells in the column ={n_masked}(Occlusion of hidden columns does not count)'
            )

            rows = [
                {
                    'variable': 'N',
                    **{
                        COL_LABEL[cohort_table_c]: f'{n_col[cohort_table_c]:,d}'
                        for cohort_table_c in COHORT_TABLE_COLS
                    },
                }
            ]

            CONTINUOUS = [
                ('Age at diagnosis, years, mean (SD)', 'age_at_diagnosis', 'mean', True),
                ('Age at diagnosis, years, median [IQR]', 'age_at_diagnosis', 'median', True),
                ('Age at enrollment, years, mean (SD)', 'age_at_enrollment', 'mean', False),
                ('Age at last EHR observation, years, mean (SD)', 'age_at_last_obs', 'mean', False),
                ('BMI, kg/m2, mean (SD)', 'bmi_value', 'mean', False),
                ('PRS311 (raw score), mean (SD)', 'PRS311', 'mean', False),
            ]

            for cohort_table_label, cohort_table_col, kind, cases_only in CONTINUOUS:
                r = {'variable': cohort_table_label}
                for cohort_table_c in COHORT_TABLE_COLS:
                    r[COL_LABEL[cohort_table_c]] = (
                        'n/a'
                        if cases_only and cohort_table_c == 'aim1_control'
                        else fmt_cont(
                            cohort_table_df.loc[masks[cohort_table_c], cohort_table_col], kind
                        )
                    )
                rows.append(r)

            for cohort_table_label, cohort_table_col in [
                ('BMI missing', 'bmi_value'),
                ('Age at last EHR observation missing', 'age_at_last_obs'),
            ]:
                key = f'{cohort_table_label}'
                cells[key] = {
                    cohort_table_c: int(
                        (masks[cohort_table_c] & cohort_table_df[cohort_table_col].isna()).sum()
                    )
                    for cohort_table_c in masks
                }
                miss_masked = cohort_table_decide(
                    {
                        key: cells[key],
                        f'{key} (observed)': {
                            cohort_table_c: n_col[cohort_table_c] - cells[key][cohort_table_c]
                            for cohort_table_c in masks
                        },
                    },
                    {key: key, f'{key} (observed)': key},
                    {
                        key: set(COHORT_TABLE_COLS) | {HIDDEN},
                        f'{key} (observed)': set(COHORT_TABLE_COLS) | {HIDDEN},
                    },
                    column_partitions=[('aim1_case', 'aim4_case', HIDDEN)],
                )
                rows.append(
                    {
                        'variable': f'{cohort_table_label}, n (%)',
                        **{
                            COL_LABEL[cohort_table_c]: fmt_cat(
                                cells[key][cohort_table_c],
                                cohort_table_c in miss_masked[key],
                                n_col[cohort_table_c],
                            )
                            for cohort_table_c in COHORT_TABLE_COLS
                        },
                    }
                )

            for cohort_table_var, cohort_table_col, levels in CATEGORICAL:
                rows.append(
                    {
                        'variable': f'{cohort_table_var}, n (%)',
                        **{COL_LABEL[cohort_table_c]: '' for cohort_table_c in COHORT_TABLE_COLS},
                    }
                )
                for lv in levels + ['Missing']:
                    key = f'{cohort_table_var}: {lv}'
                    if lv == 'Missing' and all((cells[key][c] == 0 for c in COHORT_TABLE_COLS)):
                        continue
                    rows.append(
                        {
                            'variable': f'  {lv}',
                            **{
                                COL_LABEL[cohort_table_c]: fmt_cat(
                                    cells[key][cohort_table_c],
                                    cohort_table_c in masked[key],
                                    n_col[cohort_table_c],
                                )
                                for cohort_table_c in COHORT_TABLE_COLS
                            },
                        }
                    )

            t = pd.DataFrame(rows)

            out = f'{COHORT_TABLE_TABLE_DIR}/09b_supp_table1_cohort.csv'

            t.to_csv(out, index=False)

            print(t.to_string(index=False))

            print(f'\n writes:{out}')

        # ============================================================================
        # 3. Attribute new carriers
        # ============================================================================
        if selected_task == 'attribution':
            cfg = yaml.safe_load(
                (CARRIER_ATTRIBUTION_ROOT / 'config/rowlands_13gene_extension.yaml').read_text()
            )

            scores = pd.read_csv(
                CARRIER_ATTRIBUTION_BASE / 'section8/comparable_scores_case_control.csv.gz'
            ).set_index('sample_id')

            current = (
                pd.read_csv(CARRIER_ATTRIBUTION_BASE / 'section7/unweighted_gene_carriers.csv.gz')
                .set_index('sample_id')
                .loc[scores.index]
            )

            pheno = (
                pd.read_csv(CARRIER_ATTRIBUTION_BASE / 'section9/frozen_model_inputs.csv.gz')
                .set_index('person_id')
                .loc[scores.index]
            )

            gene_weights = (
                pd.read_csv(CARRIER_ATTRIBUTION_ROOT / cfg['burden']['external_weight_table'])
                .set_index('gene')
                .ln_OR
            )

            a = scores.carried_gene_count_A.gt(0)

            q = scores.carried_gene_count_A_prime_QC.gt(0)

            b = scores.carried_gene_count_B.gt(0)

            c = scores.carried_gene_count_C.gt(0)

            flags = pd.DataFrame(
                {'A_carrier': a, 'current_nine_PTV_carrier': q, 'B_carrier': b, 'C_carrier': c}
            )

            flags['C_only_vs_A'] = c & ~a

            flags['A_only_vs_C'] = a & ~c

            flags['QC_PTV_added_vs_A'] = q & ~a

            flags['missense_added_B_vs_QC_PTV'] = b & ~q

            flags['added_four_genes_C_vs_B'] = c & ~b

            flags['new_C_current_nine_PTV'] = c & ~a & q

            flags['new_C_nine_missense'] = c & ~a & ~q & b

            flags['new_C_added_four_only'] = c & ~a & ~b

            assert np.array_equal(
                flags.C_only_vs_A.astype(int),
                flags[
                    ['new_C_current_nine_PTV', 'new_C_nine_missense', 'new_C_added_four_only']
                ].sum(axis=1),
            )

            flags.index.name = 'sample_id'

            flags.to_csv(CARRIER_ATTRIBUTION_OUT / 'carrier_transition_flags_internal.csv.gz')

            rows = []

            genes = []

            age = set(pd.read_csv(CARRIER_ATTRIBUTION_BASE / 'diagnosis_age_ids.csv.gz').person_id)

            for label, members in [('case_control', set(scores.index)), ('diagnosis_age', age)]:
                ix = sorted(members)
                f = flags.loc[ix]
                d = current.loc[ix]
                p = pheno.loc[ix]
                new = f.C_only_vs_A
                for name in flags:
                    x = f[name]
                    rows.append(
                        dict(
                            cohort=label,
                            metric=name,
                            participants=int(x.sum()),
                            cases=int((x & p.is_case.eq(1)).sum()),
                            controls=int((x & p.is_case.eq(0)).sum()),
                        )
                    )
                for g in cfg['panels']['main_13']:
                    x = d[f'{g}_main'].eq(1)
                    px = d[f'{g}_PTV_main'].eq(1)
                    mx = d[f'{g}_missense_main'].eq(1)
                    genes.append(
                        dict(
                            cohort=label,
                            gene=g,
                            carrier_n=int(x.sum()),
                            PTV_n=int(px.sum()),
                            missense_n=int(mx.sum()),
                            missense_only_n=int((mx & ~px).sum()),
                            newly_added_vs_A_n=int((x & new).sum()),
                            new_cases=int((x & new & p.is_case.eq(1)).sum()),
                            raw_burden_sum=float(x.sum() * gene_weights[g]),
                            new_raw_burden_sum=float((x & new).sum() * gene_weights[g]),
                        )
                    )

            pd.DataFrame(rows).to_csv(
                CARRIER_ATTRIBUTION_OUT / 'carrier_increment_attribution_internal.csv', index=False
            )

            gd = pd.DataFrame(genes)

            gd['share_raw_burden'] = gd.raw_burden_sum / gd.groupby(
                'cohort'
            ).raw_burden_sum.transform('sum')

            gd['share_new_raw_burden'] = gd.new_raw_burden_sum / gd.groupby(
                'cohort'
            ).new_raw_burden_sum.transform('sum')

            gd.to_csv(
                CARRIER_ATTRIBUTION_OUT / 'gene_contribution_dominance_internal.csv', index=False
            )

            columns = [
                'sample_id',
                'gene',
                'normalized_key',
                'PTV/missense',
                'keep_flag',
                'low_penetrance_flag',
                'conflict_flag',
                'VAF',
                'GQ',
                'DP',
            ]

            calls = pd.read_csv(
                CARRIER_ATTRIBUTION_BASE / 'section7/sample_variant_audit.csv.gz',
                usecols=columns,
                keep_default_na=False,
            )

            calls = calls.loc[calls.keep_flag.eq(1)].drop_duplicates(
                ['sample_id', 'gene', 'normalized_key']
            )

            calls['is_case'] = calls.sample_id.map(pheno.is_case).astype(int)

            calls['new_C_vs_A'] = calls.sample_id.map(flags.C_only_vs_A).astype(bool)

            calls['weighted_credit'] = calls.gene.map(gene_weights) / calls.groupby(
                ['sample_id', 'gene']
            ).normalized_key.transform('size')

            pos = pd.read_csv(
                CARRIER_ATTRIBUTION_BASE / 'section9/ptv_hc_transcript_position_review.tsv.gz',
                sep='\t',
                keep_default_na=False,
            )

            unc = pd.read_csv(
                CARRIER_ATTRIBUTION_BASE / 'section9/LOFTEE_multitranscript_review_internal.tsv.gz',
                sep='\t',
                keep_default_na=False,
            )

            def collapse(frame, column):
                return (
                    frame.groupby(['gene_symbol', 'normalized_key'])[column]
                    .apply(lambda x: ';'.join(sorted(set((v for v in x if v)))))
                    .to_dict()
                )

            pf = collapse(pos, 'position_review_flags')

            uf = collapse(unc, 'uncertainty_flags')

            calls['position_review_flags'] = [
                pf.get((r.gene, r.normalized_key), '') for r in calls.itertuples()
            ]

            calls['transcript_uncertainty_flags'] = [
                uf.get((r.gene, r.normalized_key), '') for r in calls.itertuples()
            ]

            variant = []

            for (g, k), d in calls.groupby(['gene', 'normalized_key']):
                n = d.loc[d.new_C_vs_A]
                variant.append(
                    dict(
                        gene=g,
                        normalized_key=k,
                        variant_class=';'.join(sorted(set(d['PTV/missense']))),
                        carrier_n=d.sample_id.nunique(),
                        case_n=d.loc[d.is_case.eq(1), 'sample_id'].nunique(),
                        new_C_vs_A_n=n.sample_id.nunique(),
                        new_C_case_n=n.loc[n.is_case.eq(1), 'sample_id'].nunique(),
                        allocated_raw_burden_credit=d.weighted_credit.sum(),
                        allocated_new_raw_credit=n.weighted_credit.sum(),
                        position_review_flags=';'.join(
                            sorted(set((v for v in d.position_review_flags if v)))
                        ),
                        transcript_uncertainty_flags=';'.join(
                            sorted(set((v for v in d.transcript_uncertainty_flags if v)))
                        ),
                        any_low_penetrance_flag=any(
                            (str(v).lower() == 'true' for v in d.low_penetrance_flag)
                        ),
                        any_conflict_flag=any((str(v).lower() == 'true' for v in d.conflict_flag)),
                    )
                )

            v = pd.DataFrame(variant).sort_values('allocated_new_raw_credit', ascending=False)

            v.to_csv(
                CARRIER_ATTRIBUTION_OUT
                / 'variant_contribution_and_suspicion_review_internal.csv.gz',
                index=False,
            )

            assert np.isclose(calls.weighted_credit.sum(), scores.B_raw_C.sum())

            assert np.isclose(
                calls.loc[calls.new_C_vs_A, 'weighted_credit'].sum(),
                scores.loc[flags.C_only_vs_A, 'B_raw_C'].sum(),
            )

            r = pd.DataFrame(rows)

            safe = r.loc[
                r.cohort.eq('case_control')
                & r.metric.isin(
                    [
                        'C_only_vs_A',
                        'QC_PTV_added_vs_A',
                        'missense_added_B_vs_QC_PTV',
                        'added_four_genes_C_vs_B',
                    ]
                ),
                ['metric', 'participants'],
            ].copy()

            safe.participants = safe.participants.map(lambda n: str(n) if n >= 20 else '<20')

            safe.to_csv(
                CARRIER_ATTRIBUTION_ROOT
                / 'results/table/46_rowlands_carrier_increment_summary.csv',
                index=False,
            )

            summary = {
                'status': 'carrier_attribution_and_variant_review_completed',
                'variant_calls_reconstruct_C_raw_burden': True,
                'allocation_rule': 'split each participant/gene weight equally among retained normalized alleles; descriptive only',
                'manual_biological_adjudication_completed': False,
                'review_flags_are_not_automatic_exclusion_rules': True,
            }

            (
                CARRIER_ATTRIBUTION_ROOT
                / 'results/report/46_rowlands_carrier_attribution_validation.json'
            ).write_text(json.dumps(summary, indent=2) + '\n')

            print(
                'Internal transition, gene/variant contribution and suspicion-review tables complete.',
                flush=True,
            )

        # ============================================================================
        # 4. Summarize model results
        # ============================================================================
        if selected_task == 'summary':
            status = json.loads(
                (
                    RESULT_SUMMARY_ROOT / 'results/report/45_rowlands_models_validation.json'
                ).read_text()
            )

            assert status['models_run'] and status['original_PRS311_repeat1_OOF_reproduced']

            cf = pd.read_csv(RESULT_SUMMARY_OUT / 'model_coefficients_internal.csv')

            auc = pd.read_csv(RESULT_SUMMARY_OUT / 'CV_AUC_summary_internal.csv').set_index(
                'version'
            )

            ds = pd.read_csv(RESULT_SUMMARY_OUT / 'CV_delta_summary_internal.csv')

            ar = pd.read_csv(RESULT_SUMMARY_OUT / 'age_explained_variance_internal.csv').set_index(
                'version'
            )

            scaling = pd.read_csv(
                RESULT_SUMMARY_OUT / 'actual_model_scaling_internal.csv'
            ).set_index(['cohort', 'version'])

            versions = [
                'A',
                'A_prime',
                'A_prime_QC',
                'B',
                'C',
                'C_without_STK11',
                'C_PTV_only',
                'C_TP53_NF1_strict_VAF',
                'C_ClinVar_min1star',
                'C_without_CDH1',
            ]

            labels = {
                'A': 'A',
                'A_prime': 'A′',
                'A_prime_QC': 'A′QC',
                'B': 'B',
                'C': 'C',
                'C_without_STK11': 'C−STK11',
                'C_PTV_only': 'C：PTV-only',
                'C_TP53_NF1_strict_VAF': 'C: Strict VAF',
                'C_ClinVar_min1star': 'C：ClinVar ≥1 star',
                'C_without_CDH1': 'C−CDH1',
            }

            rows = []

            for v in versions:
                b = cf.loc[
                    cf.version.eq(v) & cf.cohort.eq('case_control') & cf.term.eq('burden_z')
                ].iloc[0]
                a = cf.loc[
                    cf.version.eq(v) & cf.cohort.eq('diagnosis_age') & cf.term.eq('burden_z')
                ].iloc[0]
                dd = ds.loc[
                    ds.version.eq(v) & ds.comparison.eq('joint_minus_PRS311_baseline')
                ].iloc[0]
                rows.append(
                    dict(
                        version=v,
                        OR=b.OR,
                        OR_CI_low=b.CI_low,
                        OR_CI_high=b.CI_high,
                        AUC_mean=auc.loc[v, 'AUC_mean'],
                        delta_AUC_mean=dd.delta_AUC_mean,
                        delta_AUC_repeat1=dd.delta_AUC_repeat1,
                        delta_AUC_repeat1_CI_low=dd.repeat1_CI_low,
                        delta_AUC_repeat1_CI_high=dd.repeat1_CI_high,
                        age_beta=a.beta,
                        age_CI_low=a.CI_beta_low,
                        age_CI_high=a.CI_beta_high,
                        age_R_squared=ar.loc[v, 'R_squared'],
                        case_model_SD_rounded=scaling.loc[
                            ('case_control', v), 'model_SD_ddof1'
                        ].round(2),
                        age_model_SD_rounded=scaling.loc[
                            ('diagnosis_age', v), 'model_SD_ddof1'
                        ].round(2),
                    )
                )

            s = pd.DataFrame(rows)

            s.to_csv(
                RESULT_SUMMARY_ROOT / 'results/table/47_rowlands_interpretation_summary.csv',
                index=False,
            )

            from scripts.functions.stats import delong_test

            ns = {'delong_test': delong_test}

            y = pd.read_csv(
                RESULT_SUMMARY_BASE / 'section9/frozen_model_inputs.csv.gz'
            ).is_case.to_numpy(dtype=int)

            extra = []

            for repeat in range(1, 11):
                pred = np.load(RESULT_SUMMARY_OUT / f'cv_repeat_{repeat:02d}.npz')['predictions']
                names = ['PRS311_baseline'] + status['versions']
                for a, b in [('C', 'A'), ('B', 'A_prime')]:
                    extra.append(
                        dict(
                            repeat=repeat,
                            comparison=f'{a}_minus_{b}',
                            **ns['delong_test'](y, pred[names.index(a)], pred[names.index(b)]),
                        )
                    )

            e = pd.DataFrame(extra)

            e.to_csv(
                RESULT_SUMMARY_OUT / 'additional_total_AUC_comparisons_internal.csv', index=False
            )

            e.loc[e.repeat.eq(1)].to_csv(
                RESULT_SUMMARY_ROOT / 'results/table/47_rowlands_total_AUC_comparison_repeat1.csv',
                index=False,
            )

            c = s.set_index('version').loc['C']

            a = s.set_index('version').loc['A']

            total = e.loc[e.repeat.eq(1) & e.comparison.eq('C_minus_A')].iloc[0]

            inter = cf.loc[
                cf.version.eq('C')
                & cf.cohort.eq('case_control')
                & cf.term.eq('PRS311_z_x_burden_z')
            ].iloc[0]

            def table(vs):
                lines = [
                    '| Version | burden OR/SD (95% CI) | Mean OOF AUC | Mean ΔAUC relative to PRS311 | Within-case age coefficient, years/SD (95% CI) | Age R² |',
                    '|---|---|---|---|---|---|',
                ]
                for v in vs:
                    r = s.set_index('version').loc[v]
                    lines.append(
                        f'| {labels[v]} | {r.OR:.3f}（{r.OR_CI_low:.3f}–{r.OR_CI_high:.3f}） | {r.AUC_mean:.4f} | {r.delta_AUC_mean:.5f} | {r.age_beta:.2f}（{r.age_CI_low:.2f}to{r.age_CI_high:.2f}） | {r.age_R_squared * 100:.2f}% |'
                    )
                return '\n'.join(lines)

            text = f"# New burden extended analysis: result sorting and judgment\n\nThe main model and all scheduled sensitivity analyzes are now complete with N=141,204 case status (8,528 cases) and N=7,996 recorded age at diagnosis within cases. Fixed PRS311, case definition, sample, covariates and original 10×5 split. The case status and age coefficient of A reproduce the old nine-gene results; the original PRS311 first round OOF prediction also reproduces. All analyzes are extension/exploratory and current PTV biological manual adjudication is pending.\n\n## Add carrier and source\n\nOf C's 3,209 carriers, 487 were not recognized by old A. Attribution in mutually exclusive order: 148 became PTV carriers of the old nine genes under the current QC/PTV treatment, 239 became carriers only after the addition of qualified missense to the old nine genes, and 100 required the addition of four new genes, CDH1, NF1, PTEN, and STK11, to be identified. This is a breakdown of actual additions beyond the old A, without calling the full 487 missense deltas, or equating stage deltas with net headcount differentials.\n\nThe old A is different from the current PTV carrier set, so A′ retains the original carrier and only changes the weight; A′QC is a supplementary bridge. B−A′ contains both QC/PTV processing and missense changes, and B−A′QC only isolates missense.\n\nThe contribution of the new weighted score is distributed among multiple genes: BRCA2 is about 29%, BRCA1 is about 18%, and CHEK2 is about 16%. It is not covered by one new gene. The new variant with the largest contribution is CHEK2 missense, involving 60 new carriers (approximately 12%); this variant does not carry the above-mentioned PTV location/transcript review markers. Some other PTVs carry NMD, terminal exon or transcript markers; these are review tips and do not mean that they are determined to be invalid, nor can they claim that all suspected variant effects have been eliminated. The complete mutation-by-mutation sorting and labeling are kept internally, and those with less than 20 cells will not be made public.\n\n## Main comparison on fixed sample{table(['A', 'A_prime', 'A_prime_QC', 'B', 'C'])}The baseline average OOF AUC of PRS311+original covariates is{auc.loc['PRS311_baseline', 'AUC_mean']:.4f}. The average AUC increment of C relative to A is{auc.loc['C', 'AUC_mean'] - auc.loc['A', 'AUC_mean']:.5f},about{(auc.loc['C', 'AUC_mean'] - auc.loc['A', 'AUC_mean']) * 100:.2f}percentage points; the first round matching difference is{total.delta_auc:.5f}，95% DeLong CI {total.ci_low:.5f}–{total.ci_high:.5f}. Relative to the PRS311 baseline, the ΔAUC of the first round of C is{c.delta_AUC_repeat1:.5f}（{c.delta_AUC_repeat1_CI_low:.5f}–{c.delta_AUC_repeat1_CI_high:.5f}).\n\nDecomposed average ΔAUC: changing weight A′−A is about 0.00008; current PTV processing A′QC−A′ is about 0.00034; missense B−A′QC is about 0.00120; new four genes C−B is about 0.00029. The main new discriminative gain comes from missense, but the overall gain magnitude is limited. Age R² is given by A's{a.age_R_squared * 100:.2f}% raised to C{c.age_R_squared * 100:.2f}%, an absolute improvement of approx.{(c.age_R_squared - a.age_R_squared) * 100:.2f}percentage points.\n\nThe OR of case status is the conditional OR of burden per model SD when PRS is at the mean in the joint model with interactions, not the OR of carrier versus non-carrier. The within-case SD was used for the age coefficient. The actual SD and original unit coefficients of each version have been recorded in the internal comparison table. It cannot be concluded that the original scale effect is enhanced just because the OR/SD increases. CV means are reported separately from the first round CI; DeLong CI is a pairwise comparison of OOF predictions for that round and does not include the full uncertainty of retraining; the ten replicates are not ten independent cohorts. The absolute AUC includes the existing age observation window and technical covariates, focusing on explaining the increment under the same sample and covariates.\n\n## Schedule sensitivity analysis{table(['C_without_STK11', 'C_PTV_only', 'C_TP53_NF1_strict_VAF', 'C_ClinVar_min1star', 'C_without_CDH1'])}The main estimates for de-STK11, strict TP53/NF1 VAF, one-star ClinVar and de-CDH1 are close to C, with the direction unchanged; PTV-only is still relevant, but the ΔAUC drops to about 0.01198, consistent with missense bringing some additional discriminative information. The above is evidence of stability and does not single out the final protocol based on individual P values, nor does it demonstrate that these rules adequately rule out clonal hematopoiesis.\n\n## Interaction and paper presentation\n\nPRS×burden interaction beta= for C{inter.beta:.5f}（95% CI {inter.CI_beta_low:.5f}to{inter.CI_beta_high:.5f}), interaction OR={inter.OR:.3f}（{inter.CI_low:.3f}–{inter.CI_high:.3f}). A mildly negative interaction on this logistic normalized product scale does not support positive multiplicative synergy; no biological antagonism can be inferred. The age model maintains the original setting without adding interaction. Both main effects are negative, which cannot prove age synergy.\n\nSuggested paper statement: 13-gene PTV+qualified P/LP missense burden Using the approximate external weight of Rowlands all-PV, new carriers not identified by the old nine-gene PTV analysis are added; case status and diagnosis age correlation are recorded within cases in a fixed cohort, and a small additional discriminant gain is brought. The predetermined sensitivity results are generally consistent and can be used as the main result of extended analysis; it is not called a significant improvement in prediction ability, nor is it called a positive collaboration between PRS and burden. Before manual adjudication of variations, a description of the tentative definition must be retained.\n\nThe age results are always associated with documented age at diagnosis in breast cancer cases and cannot be interpreted as time-to-onset risk in the general population. External tables, pictures and texts continue to jointly implement small cell and deducible difference limits; accurate carrier/SD/variation tables remain internal."

            (
                RESULT_SUMMARY_ROOT / 'results/report/47_rowlands_results_interpretation.md'
            ).write_text(text)

            plt.style.use('default')

            fig, axes = plt.subplots(1, 2, figsize=(10, 5.5))

            english = [
                'A',
                'A-prime',
                'A-prime QC',
                'B',
                'C',
                'C minus STK11',
                'C PTV-only',
                'C strict VAF',
                'C ClinVar >=1 star',
                'C minus CDH1',
            ]

            yy = np.arange(len(s))

            axes[0].errorbar(
                s.OR,
                yy,
                xerr=np.vstack([s.OR - s.OR_CI_low, s.OR_CI_high - s.OR]),
                fmt='o',
                color='tab:blue',
                capsize=3,
            )

            axes[1].errorbar(
                s.age_beta,
                yy,
                xerr=np.vstack([s.age_beta - s.age_CI_low, s.age_CI_high - s.age_beta]),
                fmt='o',
                color='tab:orange',
                capsize=3,
            )

            axes[0].set_yticks(yy, english)

            axes[1].set_yticks(yy, [])

            for ax in axes:
                ax.invert_yaxis()
                ax.grid(False)

            axes[0].set_xlabel(
                'Case-status burden OR / model SD\nConditional on mean PRS; 95% HC1 CI'
            )

            axes[1].set_xlabel(
                'Recorded diagnosis-age beta (years / case SD)\nWithin breast cancer cases; 95% HC1 CI'
            )

            fig.suptitle(
                'Fixed-cohort extension and prespecified sensitivities\nProvisional variant definition; approximate external weights',
                fontsize=11,
            )

            fig.tight_layout()

            fig.savefig(
                RESULT_SUMMARY_ROOT / 'results/figures/47_rowlands_burden_model_comparison.svg'
            )

            plt.close(fig)

            validation = {
                'status': 'results_organized_and_interpreted',
                'models_run': True,
                'carrier_attribution_completed': True,
                'A_original_case_and_age_coefficients_reproduced': True,
                'no_significance_based_scheme_selection': True,
                'main_interaction': 'negative_on_logistic_product_scale_not_positive_synergy',
                'final_variant_definition_frozen': False,
                'public_small_cells_and_exact_comparative_SD_not_disclosed': True,
            }

            (
                RESULT_SUMMARY_ROOT
                / 'results/report/47_rowlands_result_interpretation_validation.json'
            ).write_text(json.dumps(validation, indent=2) + '\n')

            print(
                'Results report, paired total AUC comparisons, disclosure-safe summary and forest plot complete.',
                flush=True,
            )

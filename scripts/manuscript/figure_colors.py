"""Apply the requested palette only to manuscript/figures/.

Color-only SVG updates preserve every path, coordinate and numerical label.
Review Figure 2 and manuscript Figure 6 are redrawn from aggregate tables for
a publication forest layout and a narrower two-panel layout, respectively.
Run this styling pass after regenerating figures with the statistical scripts.
"""

from pathlib import Path
import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.functions.paths import PROJECT_ROOT

import os

os.environ['OPENBLAS_NUM_THREADS'] = '1'
import hashlib
import argparse
import json
import re
import xml.etree.ElementTree as ET

ROOT = PROJECT_ROOT
sys.path.insert(0, str(ROOT))
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scripts.functions.plotting import (
    PALETTE,
    MODELS,
    SCORE_ESTIMATES,
    PRIMARY_SENSITIVITY_COLORS,
    REVIEW_GROUPS,
    REVIEW_MARKERS,
    REVIEW_SUBGROUPS,
    grouped_review_forest,
    grouped_review_strata,
)

FIG = ROOT / 'manuscript/figures'
CAP = ROOT / 'manuscript/figure_captions'
PREVIEW = ROOT / 'data/manuscript_figure_colors'
HEX = re.compile(r'#[0-9a-fA-F]{6}\b')


def geometry(text):
    # Normalise just RGB tokens; every other SVG byte must remain unchanged.
    return HEX.sub('#COLOR_', text)


def recolor(path, mapping):
    text = path.read_text()
    ET.fromstring(text)
    mapping = {key.lower(): value for key, value in mapping.items()}
    before = hashlib.sha256(text.encode()).hexdigest()
    updated = HEX.sub(lambda m: mapping.get(m.group().lower(), m.group()), text)
    assert geometry(text) == geometry(updated), 'Non-color SVG content changed'
    ET.fromstring(updated)
    path.write_text(updated)
    return dict(
        file=str(path.relative_to(ROOT)),
        before_sha256=before,
        after_sha256=hashlib.sha256(updated.encode()).hexdigest(),
        geometry_and_numbers_unchanged=True,
        color_mapping=mapping,
    )


def redraw_figure6():
    versions = [
        ('C', '13 genes: primary'),
        ('C_without_STK11', '12 genes: exclude STK11'),
        ('B', '9 genes: PTV + qualifying missense'),
        ('C_PTV_only', '13 genes: PTV only'),
        ('C_TP53_NF1_strict_VAF', '13 genes: stricter TP53/NF1 VAF'),
        ('C_ClinVar_min1star', '13 genes: ClinVar ≥1 star'),
        ('C_without_CDH1', '12 genes: exclude CDH1'),
        ('A', '9 genes: original PTV baseline'),
    ]
    candidates = [
        ROOT / '9_gene_results/table/56_figure5_primary_sensitivity_estimates.csv',
        ROOT / 'results/table/56_figure5_primary_sensitivity_estimates.csv',
    ]
    source = next((path for path in candidates if path.exists()), None)
    if source is None:
        raise FileNotFoundError('Figure 6 aggregate estimates not found')
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    frame = pd.read_csv(source)
    assert len(frame) == 32 and not frame.duplicated(['cohort', 'version', 'term']).any()
    plt.rcParams['svg.fonttype'] = 'none'
    fig, axes = plt.subplots(1, 2, figsize=(10, 6.2), sharey=True)
    y = np.arange(len(versions))
    for ax, cohort, title in zip(
        axes,
        ['case_control', 'diagnosis_age'],
        ['(A) Case status (n=141,204)', '(B) Age at diagnosis (n=7,996)'],
    ):
        for term, label, marker, offset in [
            ('burden_z', 'Rare-variant burden', 'o', -0.12),
            ('PRS311_z', 'PRS311', 's', 0.12),
        ]:
            sub = (
                frame[(frame.cohort == cohort) & (frame.term == term)]
                .set_index('version')
                .loc[[v for v, label in versions]]
            )
            est, lo, hi = [sub[c].to_numpy() for c in ['estimate', 'CI_low', 'CI_high']]
            assert np.isfinite([est, lo, hi]).all() and (lo <= est).all() and (est <= hi).all()
            ax.errorbar(
                est,
                y + offset,
                xerr=np.vstack([est - lo, hi - est]),
                fmt=marker,
                color=PRIMARY_SENSITIVITY_COLORS[label],
                capsize=3,
                markersize=5,
                linewidth=1.3,
                label=label,
            )
        ax.set_title(title, fontsize=12)
        ax.grid(False)
        ax.axvline(1 if cohort == 'case_control' else 0, color='.5', linestyle='--', linewidth=0.8)
        ax.set_xlabel(
            'Odds ratio per SD'
            if cohort == 'case_control'
            else 'Difference in age at diagnosis (years per SD)'
        )
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

    axes[0].set_xscale('log')
    axes[0].xaxis.set_major_locator(FixedLocator([1, 1.2, 1.4, 1.6]))
    axes[0].xaxis.set_major_formatter(FuncFormatter(lambda v, p: f'{v:g}'))
    axes[0].xaxis.set_minor_locator(NullLocator())
    axes[0].set_xlim(0.98, 1.65)
    axes[1].set_xlim(-2.25, 0.12)
    axes[0].set_yticks(y, [label for version, label in versions], fontsize=10)
    axes[0].invert_yaxis()
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        frameon=False,
        loc='upper center',
        bbox_to_anchor=(0.65, 0.93),
        ncol=2,
        fontsize=10,
    )
    fig.suptitle('Figure 6. Primary and sensitivity estimates', fontsize=14)
    fig.text(
        0.5,
        0.015,
        'Horizontal lines: 95% CIs. Each burden is scaled by its own cohort/model SD.\n'
        'Case-status burden OR is conditional on mean PRS311; the age model has no interaction.',
        ha='center',
        fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.07, 1, 0.87])
    output = FIG / 'Figure_06_primary_and_sensitivity_estimates.svg'
    fig.savefig(output, bbox_inches='tight')
    fig.savefig(PREVIEW / (output.stem + '.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    assert digest == hashlib.sha256(source.read_bytes()).hexdigest()
    caption = CAP / 'Figure_06_caption.md'
    text = caption.read_text().replace(
        'burden (orange circles) and PRS311 (green squares)',
        'burden (blue circles, #4F86AD) and PRS311 (orange squares, #D8893D)',
    )
    text = text.replace(
        'burden (blue 2 circles, #86B6D8) and PRS311 (orange squares, #F6BD87)',
        'burden (darker blue 2 circles, #4F86AD) and PRS311 (darker orange squares, #D8893D)',
    )
    caption.write_text(text)
    record = dict(
        file=str(output.relative_to(ROOT)),
        aggregate_source=str(source.relative_to(ROOT)),
        aggregate_source_sha256=digest,
        after_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
        figure_size_inches=[10, 6.2],
        score_colors=PRIMARY_SENSITIVITY_COLORS,
        estimates_unchanged=True,
        models_refitted=False,
    )
    (CAP / 'Figure_06_style_validation.json').write_text(json.dumps(record, indent=2) + '\n')
    return record


def redraw_review_figure1():
    base = ROOT / '9_gene_results/table'
    if not base.exists():
        base = ROOT / 'results/table'
    paths = [
        base / '62_review_1008_stratum_effects.csv',
        base / '62_review_1008_cox_stratum_effects.csv',
    ]
    hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    frame = pd.concat([pd.read_csv(p) for p in paths], ignore_index=True)
    assert not frame.duplicated(['design', 'stratum']).any()
    plt.rcParams['svg.fonttype'] = 'none'
    fig, plotted = grouped_review_strata(frame)
    output = FIG / '62_review_1008_figure1.svg'
    fig.savefig(output, bbox_inches='tight')
    fig.savefig(PREVIEW / (output.stem + '.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    for p in paths:
        assert hashes[str(p.relative_to(ROOT))] == hashlib.sha256(p.read_bytes()).hexdigest()
    for row in plotted.itertuples():
        original = frame[(frame.design == row.design) & (frame.stratum == row.stratum)].iloc[0]
        assert (
            row.estimate == original.estimate
            and row.CI_low == original.CI_low
            and row.CI_high == original.CI_high
        )
    plotted.to_csv(CAP / '62_review_1008_figure1_grouped_plot_data.csv', index=False)
    record = dict(
        file=str(output.relative_to(ROOT)),
        aggregate_sources_sha256=hashes,
        after_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
        figure_size_inches=[9, 6.4],
        subgroup_encoding={
            key: dict(label=label, color=color, marker=marker)
            for key, (label, color, marker) in REVIEW_SUBGROUPS.items()
        },
        rows=len(plotted),
        sections=3,
        estimates_unchanged=True,
        models_refitted=False,
        exact_date_label_note='Requested section title retained; source EHR-only entry specification remains unverified as documented in the reproduction report',
    )
    (CAP / '62_review_1008_figure1_style_validation.json').write_text(
        json.dumps(record, indent=2) + '\n'
    )
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--figure6-only', action='store_true')
    parser.add_argument('--review1-only', action='store_true')
    args = parser.parse_args()
    CAP.mkdir(exist_ok=True)
    PREVIEW.mkdir(parents=True, exist_ok=True)
    if args.review1_only:
        record = redraw_review_figure1()
        manifest_path = CAP / 'figure_color_manifest.json'
        if manifest_path.exists():
            manifest = json.loads(manifest_path.read_text())
            manifest['review_figure1_subgroups'] = record['subgroup_encoding']
            manifest['files'] = [x for x in manifest['files'] if x['file'] != record['file']] + [
                record
            ]
            manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
        print(
            'Review Figure 1 redrawn with three sections, consistent subgroup colors/shapes, and unchanged estimates.'
        )
        return
    if args.figure6_only:
        record = redraw_figure6()
        manifest_path = CAP / 'figure_color_manifest.json'
        if manifest_path.exists():
            manifest = json.loads(manifest_path.read_text())
            manifest['primary_sensitivity_scores'] = PRIMARY_SENSITIVITY_COLORS
            manifest['files'] = [
                record if x['file'] == record['file'] else x for x in manifest['files']
            ]
            manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
        print(
            'Figure 6 redrawn at 10 × 6.2 inches with darker blue/orange; aggregate estimates unchanged.'
        )
        return
    records = []
    records.append(redraw_review_figure1())
    descriptive = [
        'Figure_02a_burden_raw_score_violin.svg',
        'Figure_02b_PRS311_raw_score_KDE.svg',
        'Figure_02c_13gene_burden_histogram_percentages.svg',
        'Figure_02d_PRS311_decile_PRS_burden_boxplots.svg',
    ]
    for name in descriptive:
        records.append(
            recolor(FIG / name, {'#2463B4': PALETTE['blue1'], '#D54848': PALETTE['red']})
        )
    model_mapping = {
        '#777777': MODELS['Covariates only'],
        '#E58A27': MODELS['PRS + covariates'],
        '#2768B4': MODELS['Burden + covariates'],
        '#329461': MODELS['Joint'],
    }
    for name in ['Figure_03a_case_status_ROC.svg', 'Figure_03b_case_status_AUC.svg']:
        records.append(recolor(FIG / name, model_mapping))
    records.append(redraw_figure6())
    records.append(
        recolor(
            FIG / '62_review_1008_figure3.svg',
            {'#86b6d8': PALETTE['blue2'], '#f6bd87': PALETTE['orange']},
        )
    )
    plt.rcParams['svg.fonttype'] = 'none'
    candidates = [
        ROOT / '9_gene_results/table/62_review_1008_sensitivity_ratios.csv',
        ROOT / 'results/table/62_review_1008_sensitivity_ratios.csv',
    ]
    source = next((path for path in candidates if path.exists()), None)
    if source is None:
        raise FileNotFoundError(
            'Review sensitivity aggregate table not found in either result location'
        )
    source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    frame = pd.read_csv(source)
    assert (frame.CI_low <= frame.estimate).all() and (frame.estimate <= frame.CI_high).all()
    fig, plotted = grouped_review_forest(frame)
    output = FIG / '62_review_1008_figure2.svg'
    fig.savefig(output)
    fig.savefig(output.with_suffix('.pdf'))
    fig.savefig(output.with_suffix('.png'), dpi=600)
    fig.savefig(PREVIEW / '62_review_1008_figure2.png', dpi=130, bbox_inches='tight')
    plt.close(fig)
    plotted.to_csv(CAP / '62_review_1008_figure2_grouped_plot_data.csv', index=False)
    assert source_hash == hashlib.sha256(source.read_bytes()).hexdigest()
    records.append(
        dict(
            file=str(output.relative_to(ROOT)),
            aggregate_source=str(source.relative_to(ROOT)),
            aggregate_source_sha256=source_hash,
            groups=[x[0] for x in REVIEW_GROUPS],
            estimates_unchanged=True,
            after_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
        )
    )
    # Match the Figure 6 caption to its new score-specific colors.
    caption = CAP / 'Figure_06_caption.md'
    text = caption.read_text()
    text = text.replace(
        'burden (orange circles) and PRS311 (green squares)',
        'burden (blue 2 circles, #86B6D8) and PRS311 (orange squares, #F6BD87)',
    )
    caption.write_text(text)
    manifest = dict(
        palette=PALETTE,
        case_control={'Case': PALETTE['blue1'], 'Control': PALETTE['red']},
        discrimination_models=MODELS,
        primary_sensitivity_scores=PRIMARY_SENSITIVITY_COLORS,
        review_groups={group: color for group, wrapped, color in REVIEW_GROUPS},
        review_group_markers=REVIEW_MARKERS,
        files=records,
        models_refitted=False,
        source_numeric_tables_modified=False,
        figure_location='manuscript/figures',
        assumption='Figure 6 burden uses blue 2 to match review Figure 3, following the last color reference in the request',
    )
    (CAP / 'figure_color_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    # Preview SVGs without changing their underlying geometry or numerical data.
    import cairosvg

    for record in records:
        path = ROOT / record['file']
        ET.parse(path)
        if path.name != '62_review_1008_figure2.svg':
            cairosvg.svg2png(
                url=str(path), write_to=str(PREVIEW / (path.stem + '.png')), output_width=1500
            )
    print(
        f'Updated {len(records)} manuscript SVGs; grouped review forests added. Numerical results preserved.'
    )


if __name__ == '__main__':
    main()

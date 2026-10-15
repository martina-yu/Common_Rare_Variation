"""Shared five-color manuscript palette and grouped review forest plot."""

import numpy as np

PALETTE = {
    'red': '#D54848',
    'blue1': '#2463B4',
    'blue2': '#86B6D8',
    'dark_grey': '#4D4D4D',
    'orange': '#F6BD87',
}
CASE_CONTROL = {'Case': PALETTE['blue1'], 'Control': PALETTE['red']}
MODELS = {
    'Covariates only': PALETTE['dark_grey'],
    'PRS + covariates': PALETTE['red'],
    'Burden + covariates': PALETTE['blue1'],
    'Joint': PALETTE['orange'],
}
SCORE_ESTIMATES = {'PRS311': PALETTE['orange'], 'Rare-variant burden': PALETTE['blue2']}
# Darker shades of the same two hues for thin Figure 6 estimate/CI lines.
PRIMARY_SENSITIVITY_COLORS = {'PRS311': '#D8893D', 'Rare-variant burden': '#4F86AD'}
REVIEW_SUBGROUPS = {
    'Non-carriers': ('Non-carriers', '#4D4D4D', 's'),
    'All carriers': ('All carriers', '#2463B4', 'o'),
    'BRCA1/2 carriers': ('BRCA1 or BRCA2 carriers', '#D8893D', 'D'),
    'Other panel carriers': ('Carriers of other panel genes', '#D54848', '^'),
}


def grouped_review_strata(frame):
    """Compact three-section forest plot with fixed subgroup colors/shapes."""
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
    import pandas as pd

    sections = [
        ('Paper covariates', 'Odds ratio, case–control\n(covariates of the paper’s joint model)'),
        (
            'Age-scale EHR + 1 year',
            'Hazard ratio, age-scale Cox\n(follow-up from 1 y after the first record)',
        ),
        (
            'Exact dates / birth-decade strata',
            'Hazard ratio, exact dates\n(EHR-only follow-up, birth-decade strata)',
        ),
    ]
    fig = plt.figure(figsize=(9, 6.4), facecolor='white')
    gs = fig.add_gridspec(1, 3, width_ratios=[3.1, 3.45, 2.05], wspace=0.025)
    labels = fig.add_subplot(gs[0, 0])
    ax = fig.add_subplot(gs[0, 1], sharey=labels)
    numbers = fig.add_subplot(gs[0, 2], sharey=labels)
    records = []
    position = 0.0
    headings = []
    for design, heading in sections:
        sub = frame[frame.design == design].set_index('stratum')
        if sub.empty:
            raise ValueError('Missing Figure 1 analysis section: ' + design)
        headings.append((-position, heading))
        position += 1.2
        for key, (label, color, marker) in REVIEW_SUBGROUPS.items():
            if key not in sub.index:
                continue
            row = sub.loc[key]
            yy = -position
            ax.errorbar(
                float(row.estimate),
                yy,
                xerr=[[row.estimate - row.CI_low], [row.CI_high - row.estimate]],
                fmt=marker,
                color=color,
                markersize=8,
                capsize=0,
                linewidth=1.8,
                zorder=3,
            )
            labels.text(0.03, yy, label, ha='left', va='center', fontsize=10, fontweight='normal')
            numbers.text(
                0.98,
                yy,
                f'{row.estimate:.2f} ({row.CI_low:.2f}–{row.CI_high:.2f})',
                ha='right',
                va='center',
                fontsize=10,
            )
            records.append(
                dict(
                    design=design,
                    stratum=key,
                    label=label,
                    y=yy,
                    color=color,
                    marker=marker,
                    estimate=row.estimate,
                    CI_low=row.CI_low,
                    CI_high=row.CI_high,
                )
            )
            position += 0.85
        position += 0.85
    for side in [labels, numbers]:
        side.set_xlim(0, 1)
        side.set_axis_off()
    for panel in [labels, ax, numbers]:
        panel.set_ylim(-position + 0.65, 0.65)
    ax.set_xscale('log')
    ax.set_xlim(0.9, 1.8)
    ax.xaxis.set_major_locator(FixedLocator([0.9, 1.0, 1.2, 1.5, 1.8]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, p: f'{v:g}'))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_yticks([])
    ax.set_axisbelow(True)
    ax.grid(axis='x', color='#E5E5E5', linewidth=0.6)
    ax.axvline(1, ls='--', color='#777777', linewidth=1, zorder=2)
    ax.spines[['top', 'right', 'left']].set_visible(False)
    ax.spines['bottom'].set_color('#AAAAAA')
    ax.set_xlabel('PRS311 effect per SD (log scale)', fontsize=10)
    numbers.text(
        0.98,
        1.015,
        'Estimate (95% CI)',
        transform=numbers.transAxes,
        ha='right',
        va='bottom',
        fontsize=10,
    )
    fig.suptitle('PRS311 effects within carrier strata', fontsize=13, y=0.975)
    fig.subplots_adjust(left=0.035, right=0.97, top=0.90, bottom=0.10)
    # Figure-level headings stay above adjacent axes patches, so long headings
    # are not clipped at the label-panel boundary.
    for yy, heading in headings:
        xx, fy = fig.transFigure.inverted().transform(labels.transData.transform((0, yy)))
        fig.text(
            xx,
            fy,
            heading,
            ha='left',
            va='center',
            fontsize=10,
            fontweight='bold',
            linespacing=1.25,
            zorder=20,
            bbox=dict(facecolor='white', edgecolor='none', pad=2),
        )
    return fig, pd.DataFrame(records)


REVIEW_GROUPS = [
    ('Main analysis', 'Main analysis', PALETTE['dark_grey']),
    ('Model and panel variants', 'Model and\npanel variants', PALETTE['blue2']),
    ('Case definition and timing', 'Case definition\nand timing', PALETTE['red']),
    ('New, from the exact-date query', 'New, from the\nexact-date query', PALETTE['orange']),
]
REVIEW_MARKERS = {
    'Main analysis': 's',
    'Model and panel variants': 'o',
    'Case definition and timing': 'D',
    'New, from the exact-date query': '^',
}


def review_group(analysis):
    if analysis == 'Paper covariates':
        return REVIEW_GROUPS[0][0]
    if analysis in [
        'Self-reported BC controls excluded',
        'Mastectomy-related records excluded',
        'Cases before exact enrolment date',
        'Cases >=1 year after exact enrolment date',
    ]:
        return REVIEW_GROUPS[3][0]
    if analysis in [
        'Rule of two (single-date cases dropped)',
        'First-record cases dropped',
        'Cases before enrolment age + 1 year',
        'Cases at/after enrolment age + 1 year',
    ]:
        return REVIEW_GROUPS[2][0]
    if analysis in [
        'Age at last EHR as spline (df=4 assumed)',
        'No covariates',
        'Within-European principal components',
        'Within-European PCs',
        'Diagnosis age used for cases',
        'Multigene carriers excluded',
        'PTV carriers only',
        'Legacy ten-gene panel',
    ]:
        return REVIEW_GROUPS[1][0]
    raise ValueError('Unknown review sensitivity group: ' + analysis)


def grouped_review_forest(frame):
    """Publication forest with aligned labels, estimates, and numerical CIs."""
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
    import pandas as pd

    frame = frame.copy()
    frame['group'] = frame.analysis.map(review_group)
    pieces = []
    headings = []
    position = 0.0
    for group, wrapped, color in REVIEW_GROUPS:
        sub = frame[frame.group == group].copy()
        if sub.empty:
            raise ValueError('Missing review forest group: ' + group)
        headings.append((-position, group))
        position += 0.95
        sub['y'] = -(position + np.arange(len(sub)) * 0.82)
        sub['color'] = color
        sub['marker'] = REVIEW_MARKERS[group]
        pieces.append(sub)
        position += len(sub) * 0.82 + 0.6
    plotted = pd.concat(pieces, ignore_index=True)
    with plt.rc_context(
        {
            'font.family': 'sans-serif',
            'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
            'pdf.fonttype': 42,
            'svg.fonttype': 'none',
        }
    ):
        fig = plt.figure(figsize=(7.8, 6.8), facecolor='white')
        gs = fig.add_gridspec(1, 3, width_ratios=[3.7, 2.45, 1.45], wspace=0.045)
        label_ax = fig.add_subplot(gs[0, 0])
        ax = fig.add_subplot(gs[0, 1], sharey=label_ax)
        estimate_ax = fig.add_subplot(gs[0, 2], sharey=label_ax)
        for group, wrapped, color in REVIEW_GROUPS:
            sub = plotted[plotted.group == group]
            ax.errorbar(
                sub.estimate.to_numpy(),
                sub.y.to_numpy(),
                xerr=[
                    (sub.estimate - sub.CI_low).to_numpy(),
                    (sub.CI_high - sub.estimate).to_numpy(),
                ],
                fmt=REVIEW_MARKERS[group],
                color=color,
                markersize=6 if group == 'Main analysis' else 5.5,
                capsize=0,
                linewidth=1.35,
                zorder=3,
            )
        for yy, heading in headings:
            label_ax.text(
                0, yy, heading, ha='left', va='center', fontsize=9, fontweight='bold', color='black'
            )
        for row in plotted.itertuples():
            label = (
                'Age at last EHR record as a spline'
                if row.analysis == 'Age at last EHR as spline (df=4 assumed)'
                else row.analysis
            )
            label_ax.text(0.035, row.y, label, ha='left', va='center', fontsize=8.2, color='black')
            estimate_ax.text(
                0,
                row.y,
                f'{row.estimate:.2f} ({row.CI_low:.2f}–{row.CI_high:.2f})',
                ha='left',
                va='center',
                fontsize=8.2,
                color='black',
            )
        lower, upper = plotted.y.min() - 0.65, 0.7
        for side in [label_ax, estimate_ax]:
            side.set_xlim(0, 1)
            side.set_ylim(lower, upper)
            side.set_axis_off()
        ax.set_ylim(lower, upper)
        ax.set_xscale('log')
        ax.set_xlim(0.6, 1.12)
        ax.xaxis.set_major_locator(FixedLocator([0.6, 0.7, 0.8, 0.9, 1.0, 1.1]))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, p: f'{v:g}'))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.set_yticks([])
        ax.set_axisbelow(True)
        ax.grid(axis='x', color='#E8E8E8', linewidth=0.5)
        ax.axvline(1, ls='--', color='#888888', linewidth=0.75, zorder=2)
        primary = float(plotted.loc[plotted.group == 'Main analysis', 'estimate'].iloc[0])
        ax.axvline(primary, ls=':', color='#999999', linewidth=0.65, zorder=2)
        ax.spines[['top', 'right', 'left']].set_visible(False)
        ax.spines['bottom'].set_color('#AAAAAA')
        ax.spines['bottom'].set_linewidth(0.6)
        ax.tick_params(axis='x', labelsize=8, width=0.6, length=3, color='#AAAAAA')
        ax.set_xlabel('Carrier / non-carrier PRS odds-ratio ratio', fontsize=8.5, labelpad=8)
        estimate_ax.text(
            0,
            1.015,
            'Estimate (95% CI)',
            transform=estimate_ax.transAxes,
            ha='left',
            va='bottom',
            fontsize=8.2,
        )
        fig.subplots_adjust(left=0.025, right=0.985, top=0.965, bottom=0.10)
    return fig, plotted


from scipy.stats import gaussian_kde


def compressed_kde(values, grid):
    """Exact full-sample Gaussian KDE via duplicate-value compression."""
    values = np.asarray(values, dtype=float)
    unique, counts = np.unique(values, return_counts=True)
    bandwidth = values.std(ddof=1) * len(values) ** (-1 / 5)
    assert bandwidth > 0
    z = (grid[:, None] - unique[None, :]) / bandwidth
    return np.exp(-0.5 * z * z) @ (counts / counts.sum()) / (bandwidth * np.sqrt(2 * np.pi))

"""Export the publication layout from the existing Figure 2 aggregate data."""

from pathlib import Path
import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.functions.paths import PROJECT_ROOT

import hashlib
import json
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

ROOT = PROJECT_ROOT
sys.path.insert(0, str(ROOT))
from scripts.functions.plotting import grouped_review_forest, REVIEW_MARKERS


def main():
    source = ROOT / 'manuscript/figure_captions/62_review_1008_figure2_grouped_plot_data.csv'
    source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    frame = pd.read_csv(source)
    fig, plotted = grouped_review_forest(frame)
    cols = ['analysis', 'estimate', 'CI_low', 'CI_high', 'color']
    pd.testing.assert_frame_equal(
        frame[cols].reset_index(drop=True), plotted[cols].reset_index(drop=True)
    )
    assert plotted.marker.tolist() == [REVIEW_MARKERS[g] for g in frame.group]
    # Check text extents before export: row labels must fit their column.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    labels, forest, numbers = fig.axes
    assert all(t.get_window_extent(renderer).x1 < forest.bbox.x0 for t in labels.texts)
    assert all(t.get_window_extent(renderer).x1 <= fig.bbox.x1 for t in numbers.texts)
    stem = ROOT / 'manuscript/figures/62_review_1008_figure2'
    for suffix in ['svg', 'pdf', 'png']:
        fig.savefig(stem.with_suffix('.' + suffix), dpi=600, facecolor='white')
    preview = ROOT / 'data/manuscript_figure_colors/62_review_1008_figure2.png'
    preview.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(preview, dpi=150, facecolor='white')
    assert hashlib.sha256(source.read_bytes()).hexdigest() == source_hash
    report = dict(
        rows=len(plotted),
        numeric_data_unchanged=True,
        colors_unchanged=True,
        markers_unchanged=True,
        source_sha256=source_hash,
        size_inches=list(fig.get_size_inches()),
        png_dpi=600,
    )
    (source.parent / '62_review_1008_figure2_style_validation.json').write_text(
        json.dumps(report, indent=2) + '\n'
    )
    plt.close(fig)
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()

"""Gene-collapsed external log-OR weighting; returned scores stay controlled."""
import numpy as np
import pandas as pd

def read_weights(path, genes):
    table = pd.read_csv(path)
    assert {'gene', 'OR'} <= set(table), 'Weight table must contain gene and OR'
    assert not table.gene.duplicated().any(), 'Duplicate gene weights'
    assert set(genes) <= set(table.gene), 'Missing panel gene weight'
    values = table.set_index('gene').loc[genes, 'OR'].astype(float)
    assert np.isfinite(values).all() and values.gt(0).all(), 'OR must be finite and positive'
    return np.log(values)

def score_panel(frame, genes, logweights, version):
    gene_cols = [f'{g}_{version}' for g in genes]
    ptv_cols = [f'{g}_PTV_{version}' for g in genes]
    miss_cols = [f'{g}_missense_{version}' for g in genes]
    assert frame[gene_cols + ptv_cols + miss_cols].isin([0, 1]).all().all()
    for g in genes:
        assert frame[f'{g}_{version}'].equals(frame[[f'{g}_PTV_{version}', f'{g}_missense_{version}']].max(axis=1)), 'Gene union must count once'
    raw = frame[gene_cols].to_numpy(dtype=float) @ logweights.loc[genes].to_numpy()
    mean = float(raw.mean())
    sd = float(raw.std(ddof=0))
    assert np.isfinite(raw).all() and sd > 0, 'Cannot standardize zero-variance burden'
    ptv = frame[ptv_cols].max(axis=1).astype(bool)
    miss = frame[miss_cols].max(axis=1).astype(bool)
    tag = np.select([ptv & miss, ptv & ~miss, miss & ~ptv], ['both', 'PTV-only', 'missense-only'], default='non-carrier')
    return (pd.DataFrame({'B_raw': raw, 'B_standardized': (raw - mean) / sd, 'carried_gene_count': frame[gene_cols].sum(axis=1), 'PTV_gene_count': frame[ptv_cols].sum(axis=1), 'missense_gene_count': frame[miss_cols].sum(axis=1), 'carrier_class': tag}), {'mean': mean, 'SD': sd, 'ddof': 0})

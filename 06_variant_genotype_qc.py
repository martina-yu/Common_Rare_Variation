"""Genotype selection and one-per-person-per-gene carrier collapse in memory."""
import pandas as pd
from functions.genotype import passing_call, normalized_native_sites, local_to_global_gt


def filter_calls(calls, candidates):
    """Caller must normalize against GRCh38 and join complete alleles, not position.

    candidates: gene, normalized_key, variant_class (PTV/missense), eligible.
    calls: sample_id, gene, normalized_key, gt, ad, gq, target_index, n_alts,
    ft, site_filters. Original GT/AD indexing must be retained through normalization.
    A single allele may belong to more than one gene; join includes gene.
    """
    required = ['sample_id', 'gene', 'normalized_key']
    if calls.duplicated(required).any() or candidates.duplicated(['gene', 'normalized_key']).any():
        raise ValueError('Duplicate allele/gene calls require reconciliation')
    selected = calls.merge(candidates.loc[candidates.eligible], on=['gene', 'normalized_key'],
                           how='inner', validate='many_to_one')
    flags = [passing_call(r.gt, r.ad, r.gq, int(r.target_index), int(r.n_alts), r.gene,
                          r.ft, r.site_filters) for r in selected.itertuples(index=False)]
    for key in ('main', 'strict_vaf'):
        selected[key] = [flag[key] for flag in flags]
    return selected


def collapse_gene_carriers(calls, sample_ids, genes, rule='main'):
    """Multiple qualifying alleles in one gene contribute one carrier indicator."""
    if len(set(sample_ids)) != len(sample_ids):
        raise ValueError('Duplicate cohort identifiers')
    if not set(calls.sample_id) <= set(sample_ids):
        raise ValueError('Calls include samples outside the frozen cohort')
    result = pd.DataFrame(index=pd.Index(sample_ids, name='sample_id'))
    for gene in genes:
        for kind in ('PTV', 'missense'):
            mask = calls[rule].eq(True) & calls.gene.eq(gene) & calls.variant_class.eq(kind)
            result[f'{gene}_{kind}_{rule}'] = result.index.isin(calls.loc[mask, 'sample_id']).astype(int)
        result[f'{gene}_{rule}'] = result[[f'{gene}_PTV_{rule}', f'{gene}_missense_{rule}']].max(axis=1)
    return result

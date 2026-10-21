"""Gene panels, same-gene transcript collapse and candidate annotation rules."""
import re
import numpy as np
import pandas as pd

PANELS = {'main_13': ['BRCA1', 'BRCA2', 'PALB2', 'CHEK2', 'ATM', 'BARD1', 'RAD51C', 'RAD51D', 'TP53', 'PTEN', 'NF1', 'STK11', 'CDH1'], 'without_stk11_12': ['BRCA1', 'BRCA2', 'PALB2', 'CHEK2', 'ATM', 'BARD1', 'RAD51C', 'RAD51D', 'TP53', 'PTEN', 'NF1', 'CDH1'], 'baseline_9': ['ATM', 'BARD1', 'BRCA1', 'BRCA2', 'CHEK2', 'PALB2', 'RAD51C', 'RAD51D', 'TP53']}
EXTERNAL_OR = {'BRCA1': 8.73, 'BRCA2': 5.68, 'PALB2': 4.3, 'CHEK2': 2.4, 'ATM': 2.16, 'BARD1': 2.34, 'RAD51C': 1.53, 'RAD51D': 1.76, 'TP53': 3.62, 'PTEN': 2.63, 'NF1': 2.01, 'STK11': 1.6, 'CDH1': 2.01}
WEIGHT_DOI = '10.1016/j.annonc.2024.07.244'
PTV = {'stop_gained', 'frameshift_variant', 'splice_donor_variant', 'splice_acceptor_variant'}


def consequence_terms(value):
    return set(re.split('[,&]', str(value)))


def collapse_transcripts(transcripts):
    """Input includes normalized_key from reference-based allele normalization.

    PTV support requires HC and PTV on the SAME retained same-gene transcript.
    Representative fields use canonical rows where available. Missense needs
    unambiguous canonical consequence. Uncertain/missing AF is not filled with 0.
    """
    rows = []
    for (gene, allele), group in transcripts.groupby(['gene_symbol', 'normalized_key'], sort=False):
        canonical = group.loc[group.is_canonical_transcript.astype(str).str.lower().eq('true')]
        ptv_hc = group.apply(lambda r: bool(consequence_terms(r.consequence) & PTV) and r.LoF == 'HC', axis=1)
        representative = (canonical if len(canonical) else group).iloc[0]
        af = pd.to_numeric(group.gvs_eur_af, errors='coerce').dropna().unique()
        missense = (not canonical.empty and canonical.consequence.nunique() == 1
                    and 'missense_variant' in consequence_terms(canonical.iloc[0].consequence))
        rows.append({**representative.to_dict(), 'gene_symbol': gene, 'normalized_key': allele,
                     'ptv_hc_any_transcript': bool(ptv_hc.any()), 'canonical_missense': missense,
                     'af_known_consistent': len(af) == 1,
                     'af_value': float(af[0]) if len(af) == 1 else np.nan,
                     'transcript_review_required': canonical.empty or canonical.transcript.nunique() > 1
                         or group.consequence.nunique() > 1 or group.LoF.nunique() > 1})
    return pd.DataFrame(rows)


def gene_coordinates(gtf_lines, genes):
    """Ensembl GRCh38 release 110 GTF gene coordinates with 2kb query padding.

    Frozen prior gene-coordinate tables can be passed by the caller instead.
    Padding is for region lookup and does not replace transcript eligibility.
    """
    result = []
    for line in gtf_lines:
        if line.startswith('#'):
            continue
        fields = line.rstrip().split('\t')
        if len(fields) < 9 or fields[2] != 'gene':
            continue
        match = re.search(r'gene_name "([^";]+)"', fields[8])
        if match and match.group(1) in genes:
            result.append({'gene': match.group(1), 'contig': 'chr' + fields[0].removeprefix('chr'),
                           'start': int(fields[3]), 'end': int(fields[4]),
                           'query_start': max(1, int(fields[3]) - 2000),
                           'query_end': int(fields[4]) + 2000, 'strand': fields[6]})
    table = pd.DataFrame(result)
    if table.empty or table.gene.duplicated().any() or set(table.gene) != set(genes):
        raise ValueError('Incomplete or ambiguous gene coordinates')
    return table

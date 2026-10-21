"""PTV/LOFTEE and ClinVar germline P/LP missense eligibility.

Variant definitions remain provisional pending final manual adjudication. No
additional final-exon or DP threshold is silently introduced.
"""
import numpy as np
from functions.clinvar import clean, exclusions, PLP, STARS, review_variant

CLINVAR_RELEASE = '2026-09-13'


def rare_consistent_af(values, upper=0.01):
    numeric = []
    for value in values:
        try:
            number = float(value)
        except (ValueError, TypeError):
            continue
        if np.isfinite(number):
            numeric.append(number)
    unique = set(numeric)
    return len(unique) == 1 and 0 <= next(iter(unique)) < upper


def qualifying_ptv(ptv_hc_any_transcript, af_values):
    return bool(ptv_hc_any_transcript and rare_consistent_af(af_values))


def qualifying_missense(canonical_missense, af_values, reference_rows, xml_record,
                        gene, min_stars=2):
    """reference_rows are normalized ClinVar allele records, not cohort counts.

    Every matching record contributes to blanket conflict/risk exclusion before
    allowed exact germline P/LP, stars, historical XML, disease and inheritance.
    """
    if not canonical_missense or not rare_consistent_af(af_values):
        return False, 'consequence_or_AF_ineligible'
    all_text = [v for row in reference_rows for v in (row.get('sig', ''), row.get('review', ''))]
    if any(row.get('conflict', '') not in {'', '.'} for row in reference_rows):
        all_text.append('conflicting classifications')
    block = exclusions(all_text)
    if block:
        return False, block
    allowed = [row for row in reference_rows if clean(row.get('sig', '')) in PLP
               and STARS.get(clean(row.get('review', '')), 0) >= min_stars]
    if not allowed or xml_record is None:
        return False, 'public_reference_or_review_unavailable'
    if len({(clean(row['sig']), clean(row['review'])) for row in allowed}) != 1:
        return False, 'reference_aggregate_disagreement_requires_review'
    ok, reason, _, _ = review_variant(xml_record, gene, allowed[0]['sig'], allowed[0]['review'], min_stars)
    return ok, reason


def ptv_position_flags(position, exon_number, exons, consequence, lof_flags='', lof_info=''):
    """Review flags only: terminal exons and <=55nt final-junction proximity.

    exons are (start, end, strand), one transcript. These are position proxies,
    not altered-transcript/NMD predictions, and do not exclude an eligible PTV.
    """
    flags = []
    parts = str(exon_number).split('/')
    if len(parts) == 2 and all(x.isdigit() for x in parts) and parts[0] == parts[1]:
        flags.append('terminal_exon_annotation')
    if not exons:
        flags.append('transcript_unresolved')
    else:
        ordered = sorted(set(exons), reverse=exons[0][2] == '-')
        if len(ordered) == 1:
            flags.append('single_exon_possible_NMD_escape')
        else:
            final, previous = ordered[-1], ordered[-2]
            if final[0] <= position <= final[1]:
                flags.append('terminal_exon_possible_NMD_escape')
            if previous[0] <= position <= previous[1]:
                distance = previous[1] - position if previous[2] == '+' else position - previous[0]
                if distance <= 55:
                    flags.append('within_55nt_final_junction_position_proxy')
    if set(str(consequence).split(',')) & {'frameshift_variant', 'splice_donor_variant', 'splice_acceptor_variant'}:
        flags.append('altered_transcript_review_required')
    if lof_flags or lof_info:
        flags.append('LOFTEE_flags_review_required')
    return flags

"""Disclosure support, not an automatic authorization to release a table.

Protect counts 1 through 20 inclusive. All results remain controlled by default.
No generic small-cell rule certifies a collection of overlapping releases.
"""
import numbers

MIN_PUBLIC_N = 21
SUPPRESSED = 'suppressed'


def format_count(value):
    if not isinstance(value, numbers.Integral) or value < 0:
        raise ValueError('Counts must be nonnegative integers')
    return str(value) if value == 0 or value >= MIN_PUBLIC_N else SUPPRESSED


def mask_related_counts(cells, equations=(), already_public=()):
    """Conservatively withhold connected arithmetic components containing small cells.

    cells maps opaque keys to counts. Each equation lists ALL nonzero-coefficient
    cells in a known exact relationship, including published totals. Propagate
    masking through the full component, rather than masking only one partner.
    A mixed component with already-public values fails closed. Equations must
    include all known cross-table and margin relationships. Unlisted relationships,
    percentages, confidence intervals and external releases require human review.
    """
    for value in cells.values():
        format_count(value)
    known = set(already_public)
    if not known <= set(cells):
        raise ValueError('Unknown published cell')
    groups = [set(equation) for equation in equations]
    if any(not group <= set(cells) for group in groups):
        raise ValueError('Unknown related cell')
    masked = {key for key, value in cells.items() if 0 < value < MIN_PUBLIC_N}
    changed = True
    while changed:
        old = set(masked)
        for group in groups:
            if group & masked:
                masked.update(group)
        changed = masked != old
    if masked & known:
        raise ValueError('Prior release requires a separate disclosure review')
    return {key: SUPPRESSED if key in masked else str(value)
            for key, value in cells.items()}

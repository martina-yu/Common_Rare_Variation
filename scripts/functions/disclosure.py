"""Joint small-cell disclosure across related aggregate tables (CLAUDE.md 6, 13).

Per-table binning is not enough once two tables share an exact arithmetic
relation: each table can be individually compliant while the pair is not. This
module decides suppression over the whole published grid at once.

The grid is {level label -> {column key -> count}}. A "level" is one category of
one variable (e.g. "Menopause: premenopausal"); a "column" is one published
population slice (e.g. ("keep", "early")). Partitions are lists of grid cells
that sum to a quantity a reader already knows, so that all-but-one determines
the last one. Three come for free from the caller's declarations:

  * columns    -- cells summing along a declared column partition, e.g.
                  pre = keep + excl, or all = early + late;
  * levels     -- all levels of one variable in one column sum to that column's n;
  * aggregates -- a declared parent level equals the sum of its children
                  (a coarsened category, or a two-way cross-tab collapsing to
                  its margin).

Rules, iterated to a fixed point:
  1. bin any published cell whose count is 1-19;
  2. if exactly one cell of a partition is binned, bin a second one -- the
     smallest cell still open -- or the first returns by subtraction.
A count of exactly 0 is not a disclosure and never forces a partner. Cells
listed as already public are never spent as a partner, because suppressing
something that is already published buys nothing.

render() writes a decided cell. A genuinely small cell gets its bin; a cell that
is not small but had to be withheld to protect a partner is written
"suppressed (secondary)" -- printing a bin label there would state something
false about its size.
"""

MIN_N, LOW_N = 20, 5
SECONDARY = "suppressed (secondary)"


def decide(cells, var_of, published, column_partitions=(), aggregates=(), public=frozenset()):
    """Return {level -> set of masked column keys}.

    cells       {level -> {col -> int}}
    var_of      {level -> variable name}
    published   {level -> set of cols actually printed somewhere}
    column_partitions  iterable of column-key tuples that sum exactly
                       (e.g. ("pre", "keep", "excl") where pre = keep + excl)
    aggregates  iterable of (parent_level, [child_levels]); the relation is
                checked in every column both sides publish
    public      set of (level, col) already published elsewhere
    """
    grid = {(lv, c): cells[lv][c] for lv in cells for c in published[lv]}
    masked = {g for g, v in grid.items() if 0 < v < MIN_N}

    parts = []
    for lv in cells:
        for cp in column_partitions:
            p = [(lv, c) for c in cp if (lv, c) in grid]
            if len(p) > 1:
                parts.append(p)
    by_var = {}
    for lv in cells:
        by_var.setdefault(var_of[lv], []).append(lv)
    for lvs in by_var.values():
        cols = set.intersection(*[set(published[lv]) for lv in lvs])
        for c in cols:
            parts.append([(lv, c) for lv in lvs])
    for parent, children in aggregates:
        cols = set(published[parent]).intersection(*[set(published[ch]) for ch in children])
        for c in cols:
            parts.append([(parent, c)] + [(ch, c) for ch in children])

    changed = True
    while changed:
        changed = False
        for part in parts:
            if all(g in public for g in part):
                continue  # nothing here can still be protected
            if len([g for g in part if g in masked]) != 1:
                continue
            openg = [g for g in part if g not in masked and grid[g] > 0 and g not in public]
            if openg:
                masked.add(min(openg, key=lambda g: grid[g]))
                changed = True

    out = {lv: set() for lv in cells}
    for lv, c in masked:
        out[lv].add(c)
    return out


def render(count, is_masked, denom):
    """One published cell: exact count with percentage, a bin, or a withholding."""
    if not is_masked:
        return f"{count} {100 * count / denom:.1f}%" if denom else str(count)
    if count == 0:
        return "0"
    if count < LOW_N:
        return "<5"
    if count < MIN_N:
        return "5-19"
    return SECONDARY


def format_count(count, style='cohort'):
    """Preserve the existing cohort/QC count displays, including the n=20 boundary."""
    count = int(count)
    if count == 0 or count >= 20:
        return f'{count:,d}' if style == 'sample_qc' else str(count)
    return ('<5' if count < 5 else '5-19') if style == 'sample_qc' else '<20'

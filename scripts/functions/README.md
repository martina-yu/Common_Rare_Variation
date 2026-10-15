# Shared analysis functions

Scientific selection, model specification and output assembly live in numbered
scripts. These modules retain existing, genuinely reusable responsibilities:

- `paths`, `io`, `provenance`, `bigquery`, `cli`: paths, hashes, source/input locks, queries and explicit CLI selection.
- `disclosure`: existing joint suppression and explicitly distinct count-display policies.
- `stats`, `burden`, `frailty`: reused statistics, scoring kernels and frailty mathematics.
- `prs`, `variants`, `clinvar`, `site_qc`, `genotype_qc`, `sample_variant_audit`: low-level allele/genotype and assertion mechanics.
- `cohort`: release-checked ancestry-file extraction.
- `age_group_contrasts`, `ehr_age_audit`, `nine_gene_missense_age`: existing helpers shared with historical/exploratory analyses.
- `plotting`: shared rendering support.

Domain-specific existing modules are retained where merging them into one
`variants.py` or `stats.py` would obscure their purpose. Whole-stage implementation
modules and lazy task/forwarding modules from the interrupted pass are archived.

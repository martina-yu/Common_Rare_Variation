# Revision log: 1006 manuscript to thirteen-gene main score

User instruction: use the research skill to revise manuscript 1006 against the thirteen-gene revision outline. Original DOCX and analysis data were preserved; no models were rerun.

| Outline commitment | Completed revision |
|---|---|
| Sole thirteen-gene burden in main | Replaced construction, Results, Abstract and Discussion; removed obsolete main-panel methods and values. |
| Accurate construction | Added all thirteen genes, PTV retained-transcript/LOFTEE rule, contextual eligible missense rules, actual genotype QC and once-per-gene weighted formula. |
| Correct PRS | Distinguished 313 intended, 304 matched and 302 retained loci in PRS311; retained dosage-imputation description. |
| Updated results | Main ORs 1.554/1.269; interaction OR 0.972; mean AUCs 0.7452/0.7682/0.7588/0.7815; mean joint-minus-PRS increment 0.01328. |
| Diagnosis age | PRS −0.43 years/SD, burden −1.82 years/SD, R² 0.02628; conclusions restricted to recorded diagnosis age among cases. |
| Historical comparisons in supplement | All five historical stages and four sequential increments included with descriptive names; historical increment 0.001905 distinguished from the main contrast. |
| Figure consistency | Generated a branching flow, current-only ROC/AUC and observed joint proportions; moved full sensitivity/legacy forest figures to supplement. |
| Negative interaction | Retained exploratory gene interaction table, box plots, robust audit and Barnes 2020 citation; no biological antagonism or positive synergy claim. |
| Literature comparison | Added verified Guare 2025 comparison, with continuous-model efficiency stated only as a possible explanation. |
| Status and limits | Reported chronology, pending adjudication, external-weight mismatch, EHR ascertainment, case conditioning and ancestry scope. |

Independent evidence and draft reviews identified missing historical tables, multipanel figure/caption mismatches and an overstated missingness measurement; these were corrected. The chromosome-21 approximately 1,500-locus proxy is now stated explicitly.

The clean DOCX body was rebuilt using the original document package and styles. `main_revision.diff` is a unified text diff against extracted source text; it is not Word tracked changes or an ARS block-patch artifact. `main_13gene_revised.md` and `supplement_13gene_revised.md` provide readable review copies. Numerical source hashes and package/XML checks are recorded in `revision_validation.json`. The deterministic acronym checker was run; gene symbols, score names, software/filter labels and mathematical subscripts create advisory false positives/coverage limits, recorded separately.

Remaining author/publication work: original correspondence, acknowledgements, ethics, funding, competing interests and contribution placeholders need completion; variant adjudication and family-history phenotype coding remain unresolved. No prospective calibration or clinical-use validation is claimed. Full Word pagination was not rendered because a Word/LibreOffice renderer was unavailable.

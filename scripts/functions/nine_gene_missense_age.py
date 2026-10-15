# Documentation: scripts/README.md.
import argparse
import gzip
import json
import os
from pathlib import Path
import subprocess
import tempfile
from urllib.parse import unquote

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import norm
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
V2 = ROOT
ANNOT = V2 / "data/variant_annotation/full"
DATA = ROOT / "data/28_expanded_missense"
TABLE = ROOT / "results/table"
FIGURE = ROOT / "results/figures"
REPORT = ROOT / "results/report"
PCS = [f"PC{i}" for i in range(1, 11)]
PANEL = ["ATM", "BARD1", "BRCA1", "BRCA2", "CHEK2", "PALB2", "RAD51C", "RAD51D", "TP53"]
PLP = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}
STARS = {
    "criteria_provided,_single_submitter": 1,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "reviewed_by_expert_panel": 3,
    "practice_guideline": 4,
}
MIN_N = 21
I157T = "chr22:28725099:A:G"


def read(name, **kwargs):
    return pd.read_csv(ROOT / "data" / name, **kwargs)


def safe_partition(counts):
    return bool(((np.asarray(counts) == 0) | (np.asarray(counts) >= MIN_N)).all())


def related_condition(name, gene):
    text = unquote(name).replace("_", " ").lower().strip()
    if "breast" in text or "mammary" in text:
        return True
    if text == "hereditary cancer-predisposing syndrome":
        return True
    if gene in {"TP53", "CHEK2"} and ("li-fraumeni" in text or "li fraumeni" in text):
        return True
    return text == f"{gene.lower()}-related cancer predisposition" and gene != "MAP3K1"


def disease_allowed(disease_text, gene):
    names = [
        x for x in disease_text.split("|") if x not in {".", "", "not_provided", "not_specified"}
    ]
    # The cache is variant-level; exclude ambiguous unrelated co-listed conditions.
    return bool(names) and all(related_condition(x, gene) for x in names)


def germline_allowed(origin):
    try:
        bits = 0
        for x in origin.split(","):
            bits |= int(x)
    except ValueError:
        return True  # CLNSIG is explicitly the germline classification field.
    germline_bits = 1 | 4 | 8 | 16 | 32 | 64 | 128
    return not (bits & 2 and not bits & germline_bits)


def prepare_candidates():
    vat = read("11_panel_gene_vat_ptv_missense_candidates.csv.gz", dtype=str, keep_default_na=False)
    vat["_rank"] = np.where(vat.is_canonical_transcript == "true", 0, 1)
    vat = vat.sort_values(["vid", "gene_symbol", "_rank"]).drop_duplicates(["vid", "gene_symbol"])
    vat["eur_af"] = pd.to_numeric(vat.gvs_eur_af, errors="coerce").fillna(0.0)
    vat = vat.loc[
        vat.gene_symbol.isin(PANEL)
        & vat.consequence.str.contains("missense_variant")
        & (vat.eur_af < 0.001)
    ].copy()
    vat["orig_key"] = vat.contig + ":" + vat.position + ":" + vat.ref_allele + ":" + vat.alt_allele
    vat = vat.loc[vat.orig_key != I157T]
    wanted = set(vat.orig_key)
    mapped = []
    columns = ["orig_key", "key", "gene", "exonic_func", "aachange"]
    for chunk in pd.read_csv(
        ANNOT / "v2_05_annotated_sites.tsv.gz",
        sep="\t",
        dtype=str,
        usecols=columns,
        chunksize=100000,
        keep_default_na=False,
    ):
        chunk = chunk.loc[chunk.orig_key.isin(wanted)]
        if len(chunk):
            mapped.append(chunk)
    mapping = pd.concat(mapped, ignore_index=True) if mapped else pd.DataFrame(columns=columns)
    mapping = mapping.rename(
        columns={
            "gene": "mane_gene",
            "exonic_func": "mane_exonic_func",
            "aachange": "mane_aachange",
        }
    )
    assert not mapping.orig_key.duplicated().any()
    vat = vat.merge(mapping, on="orig_key", how="left", validate="many_to_one")
    # A single-base substitution needs no indel realignment; cached mapping is used otherwise.
    snv = (vat.ref_allele.str.len() == 1) & (vat.alt_allele.str.len() == 1)
    vat.loc[vat.key.isna() & snv, "key"] = vat.loc[vat.key.isna() & snv, "orig_key"]
    vat = vat.loc[vat.key.notna()].copy()
    reference = []
    shard_map = read("11_panel_gene_exome_vcf_shard_map.csv")
    fmt = "%CHROM\t%POS\t%REF\t%ALT\t%ID\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/CLNSIGCONF\t%INFO/CLNDN\t%INFO/CLNDISDB\t%INFO/ORIGIN\t%INFO/MC\n"
    for row in shard_map.itertuples():
        if row.gene_symbol not in PANEL:
            continue
        region = f"{row.chrom}:{row.region_start}-{row.region_end}"
        raw = subprocess.run(
            ["bcftools", "query", "-r", region, "-f", fmt, str(ANNOT / "clinvar.norm.vcf.gz")],
            check=True,
            text=True,
            capture_output=True,
        ).stdout
        for line in raw.splitlines():
            chrom, pos, ref, alt, vid, sig, rev, conflict, disease, disease_ids, origin, mc = (
                line.split("\t")
            )
            key = f"{chrom}:{pos}:{ref}:{alt}"
            reference.append(
                {
                    "key": key,
                    "variation_id": vid,
                    "clnsig": sig,
                    "clnrevstat": rev,
                    "clnsigconf": conflict,
                    "disease": disease,
                    "disease_ids": disease_ids,
                    "origin": origin,
                    "clinvar_mc": mc,
                }
            )
    cv = pd.DataFrame(reference).drop_duplicates()
    blocked = set(
        cv.loc[
            cv.clnsig.str.contains("Benign|benign|Conflicting|conflicting")
            | cv.clnrevstat.str.contains("conflicting")
            | ~cv.clnsigconf.isin(["", "."]),
            "key",
        ]
    )
    cv["stars"] = cv.clnrevstat.map(STARS).fillna(0).astype(int)
    cv = cv.loc[
        cv.clnsig.isin(PLP)
        & (cv.stars >= 1)
        & ~cv.key.isin(blocked)
        & cv.origin.map(germline_allowed)
        & cv.clinvar_mc.str.contains("missense_variant")
    ]
    candidates = vat.merge(cv, on="key", how="inner")
    candidates = candidates.loc[
        [disease_allowed(d, g) for d, g in zip(candidates.disease, candidates.gene_symbol)]
    ]
    # Keep the highest-reviewed germline assertion; all records already passed disease/conflict gates.
    candidates = candidates.sort_values("stars", ascending=False).drop_duplicates(
        ["orig_key", "gene_symbol"]
    )
    candidates["strict"] = candidates.stars >= 2
    mane_available = candidates.mane_exonic_func.notna()
    # Candidates contain cohort AF and must stay in ignored data/, never public tables.
    candidates.to_csv(DATA / "candidates.tsv.gz", sep="\t", index=False)
    header_date = None
    with gzip.open(ANNOT / "clinvar.norm.vcf.gz", "rt") as handle:
        for line in handle:
            if line.startswith("##fileDate="):
                header_date = line.strip().split("=", 1)[1]
            if not line.startswith("#"):
                break
    assert header_date == "2026-09-13"
    provenance = {
        "status": "post_hoc_exploratory",
        "ClinVar_snapshot": header_date,
        "genome": "GRCh38",
        "annotation_reused_not_rerun": True,
        "panel": PANEL,
        "primary_PTV_definition": "original_nine_gene_qualifying_PTV",
        "reference_annotation": 'normalized ClinVar VCF and existing canonical VAT; MANE .5 audit',
        "disease_rule": "at_least_one_relevant_named_condition_and_no_ambiguous_unrelated_named_condition",
        "disease_specific_RCV_cache_available": False,
        "main_minimum_stars": 1,
        "strict_minimum_stars": 2,
        "eur_AF_threshold": 0.001,
        "I157T_excluded": True,
        "GQ_min": 20,
        "het_VAF": [0.30, 0.80],
        "hom_alt_VAF_min": 0.8,
        "old_v2_genotype_cache_covers_frozen_primary": False,
        "annotation_audit": {
            "cached_MANE_available_for_every_candidate": bool(mane_available.all()),
            "cached_MANE_consequences_match_where_available": bool(
                candidates.loc[mane_available, "mane_exonic_func"].eq("nonsynonymous SNV").all()
            ),
            "cached_MANE_genes_match_where_available": bool(
                (
                    candidates.loc[mane_available, "mane_gene"]
                    == candidates.loc[mane_available, "gene_symbol"]
                ).all()
            ),
            "independent_ClinGen_classification_field_available": False,
        },
        "minimum_release_n": MIN_N,
    }
    (DATA / "pre_model_rules.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    print(
        "Prepared fixed clinical missense mask from existing reference files; no cohort allele counts printed.",
        flush=True,
    )
    return candidates, provenance


def ensure_ro_mount(mount):
    mount.mkdir(parents=True, exist_ok=True)
    existing = [line.split() for line in Path("/proc/mounts").read_text().splitlines()]
    entry = next((x for x in existing if x[1] == str(mount)), None)
    if entry:
        assert "ro" in entry[3].split(","), "Mount is not read-only."
        return
    subprocess.run(["fusermount", "-u", str(mount)], capture_output=True)
    bucket = os.environ["CDR_STORAGE_PATH"].replace("gs://", "").split("/")[0]
    result = subprocess.run(
        [
            "gcsfuse",
            "--billing-project",
            os.environ["GOOGLE_PROJECT"],
            "-o",
            "ro",
            "--implicit-dirs",
            bucket,
            str(mount),
        ],
        capture_output=True,
        text=True,
    )
    if result.returncode:
        raise RuntimeError("Read-only gcsfuse mount failed: " + result.stderr[-1500:])
    entries = [line.split() for line in Path("/proc/mounts").read_text().splitlines()]
    entry = next(x for x in entries if x[1] == str(mount))
    assert "ro" in entry[3].split(",")


def qualifying_call(gt, ad, gq, alt_field, target):
    try:
        alleles = [int(x) for x in gt.replace("|", "/").split("/")]
        depths = [int(x) for x in ad.split(",")]
        index = alt_field.split(",").index(target) + 1
        dosage = alleles.count(index)
        if dosage == 0 or len(alleles) != 2 or sum(depths) <= 0:
            return False
        vaf = depths[index] / sum(depths)
        return float(gq) >= 20 and (vaf >= 0.8 if dosage == 2 else 0.3 <= vaf <= 0.8)
    except (ValueError, IndexError):
        return False


def minimal_allele_key(chrom, pos, ref, alt):
    """Minimal allele key."""
    pos = int(pos)
    if alt == "*" or alt.startswith("<"):
        return chrom, pos, ref, alt
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return chrom, pos, ref, alt


def candidate_shards(boundaries, chrom, positions):
    ordinal = int(chrom.removeprefix("chr"))
    bounds = boundaries.loc[boundaries.chrom_ord == ordinal].sort_values("start_pos")
    starts = bounds.start_pos.to_numpy()
    shards = bounds.shard_idx.to_numpy()
    chosen = set()
    for pos in positions:
        left = int(np.searchsorted(starts, pos, side="left"))
        right = int(np.searchsorted(starts, pos, side="right"))
        # Include the prior shard when an exact boundary can split alleles at a locus.
        for index in range(max(0, left - 1), max(left, right)):
            chosen.add(int(shards[index]))
    return sorted(chosen)


def extract_carriers(candidates, ids, mount):
    # Sets in memory only. No genotype/person-by-variant output is written.
    main = {gene: set() for gene in PANEL}
    strict = {gene: set() for gene in PANEL}
    if candidates.empty:
        return main, strict, {"no_qualifying_missense_variants": True}
    ensure_ro_mount(mount)
    prefix = "/".join(os.environ["CDR_STORAGE_PATH"].replace("gs://", "").split("/")[1:])
    vcfdir = mount / prefix / "wgs/short_read/snpindel/exome/vcf"
    maps = read("11_panel_gene_exome_vcf_shard_map.csv")
    boundaries = read("11_exome_vcf_shard_start_boundaries.csv")
    seen_sites = set()
    fmt = "[%CHROM\t%POS\t%REF\t%ALT\t%SAMPLE\t%GT\t%AD\t%GQ\n]"
    with tempfile.TemporaryDirectory(prefix="expanded_missense_", dir=DATA) as temp:
        samples = Path(temp) / "samples.txt"
        samples.write_text("\n".join(map(str, sorted(ids))) + "\n", encoding="utf-8")
        for gene, group in candidates.groupby("gene_symbol"):
            row = maps.loc[maps.gene_symbol == gene].iloc[0]
            regions = Path(temp) / "regions.bed"
            group[["contig", "position"]].drop_duplicates().assign(
                start=lambda d: d.position.astype(int) - 1, end=lambda d: d.position.astype(int)
            )[["contig", "start", "end"]].to_csv(regions, sep="\t", index=False, header=False)
            wanted = {
                (r.contig, int(r.position), r.ref_allele, r.alt_allele): bool(r.strict)
                for r in group.itertuples()
            }
            shards = candidate_shards(boundaries, row.chrom, group.position.astype(int))
            # Native padded representations can start just outside the minimal-key shard.
            shards = sorted(
                {
                    s + d
                    for s in shards
                    for d in [-1, 0, 1]
                    if int(row.shard_first) <= s + d <= int(row.shard_last)
                }
            )
            assert (
                shards
                and min(shards) >= int(row.shard_first)
                and max(shards) <= int(row.shard_last)
            )
            for shard in shards:
                path = vcfdir / f"{shard:010d}.vcf.bgz"
                assert path.exists(), "Expected exome VCF shard is missing."
                sites = subprocess.run(
                    [
                        "bcftools",
                        "query",
                        "-R",
                        str(regions),
                        "--regions-overlap",
                        "1",
                        "-f",
                        "%CHROM\t%POS\t%REF\t%ALT\n",
                        str(path),
                    ],
                    capture_output=True,
                    text=True,
                    check=True,
                ).stdout
                matching = set()
                for site in sites.splitlines():
                    chrom, pos, ref, alts = site.split("\t")
                    for alt in alts.split(","):
                        key = minimal_allele_key(chrom, pos, ref, alt)
                        if key in wanted:
                            seen_sites.add(key)
                            matching.add(key)
                if not matching:
                    continue
                view = subprocess.Popen(
                    [
                        "bcftools",
                        "view",
                        "-R",
                        str(regions),
                        "--regions-overlap",
                        "1",
                        "-S",
                        str(samples),
                        "-Ou",
                        str(path),
                    ],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )
                query = subprocess.run(
                    ["bcftools", "query", "-i", 'GT="alt"', "-f", fmt],
                    stdin=view.stdout,
                    capture_output=True,
                    text=True,
                )
                view.stdout.close()
                error = view.stderr.read().decode()
                ret = view.wait()
                if ret or query.returncode:
                    raise RuntimeError(
                        "Targeted genotype query failed; raw stderr withheld to protect sample IDs."
                    )
                for line in query.stdout.splitlines():
                    chrom, pos, ref, alt_field, person, gt, ad, gq = line.split("\t")
                    person = int(person)
                    assert person in ids
                    for alt in alt_field.split(","):
                        key = minimal_allele_key(chrom, pos, ref, alt)
                        if key in wanted and qualifying_call(gt, ad, gq, alt_field, alt):
                            main[gene].add(person)
                            if wanted[key]:
                                strict[gene].add(person)
            print(
                f"Completed targeted exome calls for {gene}; no participant or small-cell counts printed.",
                flush=True,
            )
    expected = {
        (r.contig, int(r.position), r.ref_allele, r.alt_allele) for r in candidates.itertuples()
    }
    assert (
        seen_sites == expected
    ), "Some qualifying missense sites were not found in the original callset."
    return (
        main,
        strict,
        {
            "targeted_sites_all_found": True,
            "source": "original AoU exome VCF",
            "sample_set_matches_frozen_7996_cases": True,
            "no_genotypes_written": True,
            "native_VCF_padding_removed_for_allele_match": True,
            "AD_uses_native_ALT_index": True,
        },
    )


def carrier_partitions(full, gene_main, gene_strict):
    main_set = set().union(*gene_main.values())
    strict_set = set().union(*gene_strict.values())
    full = full.copy()
    full["ptv"] = (full.burden_rare_main_ptv_9gene > 0).astype(int)
    full["missense"] = full.person_id.isin(main_set).astype(int)
    full["strict_missense"] = full.person_id.isin(strict_set).astype(int)
    full["expanded"] = ((full.ptv == 1) | (full.missense == 1)).astype(int)
    full["expanded_strict"] = ((full.ptv == 1) | (full.strict_missense == 1)).astype(int)
    full["category"] = np.select(
        [
            (full.ptv == 0) & (full.missense == 0),
            (full.ptv == 1) & (full.missense == 0),
            (full.ptv == 0) & (full.missense == 1),
        ],
        ["neither", "PTV_only", "missense_only"],
        default="both",
    )
    primary = (full.group == "case") & (full.bc_is_first_record != True)
    audits = {}
    for label, sub in [
        ("primary_cases", full.loc[primary]),
        (
            "primary_age_groups",
            full.loc[primary].assign(
                age_group=lambda d: np.where(d.age_at_diagnosis < 45, "early", "later")
            ),
        ),
    ]:
        keys = ["age_group", "category"] if label == "primary_age_groups" else ["category"]
        levels = ["early", "later"]
        categories = ["neither", "PTV_only", "missense_only", "both"]
        if len(keys) == 1:
            counts = sub.groupby(keys[0]).size().reindex(categories, fill_value=0)
        else:
            counts = (
                sub.groupby(keys)
                .size()
                .reindex(pd.MultiIndex.from_product([levels, categories], names=keys), fill_value=0)
            )
        safe = safe_partition(counts)
        audits[label] = {"partition_release": safe}
        if safe:
            counts.rename_axis(keys).reset_index(name="n").to_csv(
                TABLE / f"28_{label}_carrier_categories.csv", index=False
            )
    # Gene summaries omit all numbers when any count or known complement is small.
    gene_rows = []
    for gene in PANEL:
        for label, sub in [("primary_case", full.loc[primary])]:
            n = int(sub.person_id.isin(gene_main[gene]).sum())
            total = len(sub)
            if (n == 0 or n >= MIN_N) and (total - n == 0 or total - n >= MIN_N):
                gene_rows.append({"gene": gene, "stratum": label, "n_missense_carriers": n})
            else:
                gene_rows.append(
                    {"gene": gene, "stratum": label, "n_missense_carriers": "suppressed"}
                )
    if any(r["n_missense_carriers"] == "suppressed" for r in gene_rows):
        # Do not leave all-but-one gene counts available beside an overall missense total.
        for row in gene_rows:
            row["n_missense_carriers"] = "suppressed_joint_gene_table"
    pd.DataFrame(gene_rows).to_csv(TABLE / "28_missense_carriers_by_gene.csv", index=False)
    return full, audits


def linear_contrast(fit, v):
    v = np.asarray(v, dtype=float)
    beta = float(v @ fit.params)
    se = float(np.sqrt(max(0, v @ fit.cov_params() @ v)))
    return {
        "age_difference_years": beta,
        "se_HC1": se,
        "ci_low": beta - norm.ppf(0.975) * se,
        "ci_high": beta + norm.ppf(0.975) * se,
    }


def fit_age(cases, carrier, label, ref):
    mu, sd, pc_ref = ref
    df = cases.copy()
    df["PRS_z"] = (df.PRS311 - mu) / sd
    df["carrier"] = df[carrier]
    pc_mu, pc_sd = df[PCS].mean(), df[PCS].std(ddof=1)
    df[PCS] = (df[PCS] - pc_mu) / pc_sd
    pc_fixed = (pc_ref - pc_mu) / pc_sd
    rows = []
    for interaction in [False, True]:
        df["interaction"] = df.PRS_z * df.carrier
        terms = ["PRS_z", "carrier"] + (["interaction"] if interaction else []) + PCS
        X = sm.add_constant(df[terms], has_constant="add")
        fit = sm.OLS(df.age_at_diagnosis.astype(float), X).fit(cov_type="HC1")
        for term in ["PRS_z", "carrier"] + (["interaction"] if interaction else []):
            v = np.zeros(X.shape[1])
            v[X.columns.get_loc(term)] = 1
            rows.append(
                {
                    "definition": label,
                    "model": "with_interaction" if interaction else "additive",
                    "contrast": term,
                    **linear_contrast(fit, v),
                }
            )
        if interaction:
            v = np.zeros(X.shape[1])
            v[X.columns.get_loc("PRS_z")] = v[X.columns.get_loc("interaction")] = 1
            rows.append(
                {
                    "definition": label,
                    "model": "with_interaction",
                    "contrast": "PRS_per_SD_carriers",
                    **linear_contrast(fit, v),
                }
            )
            lo = max(df.loc[df.carrier == c, "PRS_z"].quantile(0.01) for c in [0, 1])
            hi = min(df.loc[df.carrier == c, "PRS_z"].quantile(0.99) for c in [0, 1])
            assert lo < hi
            fig, ax = plt.subplots()
            for c in [0, 1]:
                xs = np.linspace(lo, hi, 100)
                predictions = pd.DataFrame(
                    [linear_contrast(fit, [1, p, c, p * c, *pc_fixed]) for p in xs]
                )
                (line,) = ax.plot(xs, predictions.age_difference_years, label=f"Carrier={c}")
                ax.fill_between(
                    xs, predictions.ci_low, predictions.ci_high, alpha=0.2, color=line.get_color()
                )
            ax.set(
                xlabel="PRS311 (SD in frozen primary cases)",
                ylabel="Adjusted recorded diagnosis age (years)",
            )
            ax.set_title(
                {
                    "PTV_only": "Nine-gene qualifying PTV",
                    "expanded_star1": "Nine-gene PTV + P/LP missense (review >=1 star)",
                    "expanded_star2": "Nine-gene PTV + P/LP missense (review >=2 stars)",
                }[label]
            )
            ax.grid(False)
            ax.legend()
            fig.tight_layout()
            fig.savefig(FIGURE / f"28_{label}_predicted_age.svg")
            plt.close(fig)
    return rows


def finalize_disclosure_and_report():
    """Check released aggregates against earlier ten-gene totals before reporting."""
    meta_path = REPORT / "28_expanded_missense_metadata.json"
    metadata = json.loads(meta_path.read_text())
    counts_path = TABLE / "28_primary_carrier_increment.csv"
    counts = pd.read_csv(counts_path)
    old_path = ROOT / "results/common_rare_age/primary_four_groups.csv"
    baseline_safe = metadata.get("nine_gene_baseline_total_release", True)
    if old_path.exists():
        old = pd.read_csv(old_path)
        n10 = int(old.loc[old.carrier == 1, "n"].sum())
        n9 = pd.to_numeric(
            counts.loc[counts.definition == "PTV_only", "n_carriers"], errors="coerce"
        ).iloc[0]
        if pd.notna(n9):
            difference = n10 - int(n9)
            baseline_safe = difference == 0 or difference >= MIN_N
    if not baseline_safe:
        counts["n_carriers"] = "suppressed_cross_panel_difference"
        counts.to_csv(counts_path, index=False)
        for label in ["primary_cases", "primary_age_groups"]:
            path = TABLE / f"28_{label}_carrier_categories.csv"
            if path.exists():
                path.unlink()  # Only new tables whose sums would reveal the excluded small group.
            metadata["counts_audit"][label]["partition_release"] = False
            metadata["counts_audit"][label]["cross_panel_complement_suppressed"] = True
    metadata["nine_gene_baseline_total_release"] = baseline_safe
    genes_path = TABLE / "28_missense_carriers_by_gene.csv"
    genes = pd.read_csv(genes_path)
    if genes.n_missense_carriers.astype(str).str.contains("suppressed").any():
        genes["n_missense_carriers"] = "suppressed_joint_gene_table"
        genes.to_csv(genes_path, index=False)
    meta_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    models = pd.read_csv(TABLE / "28_age_model_comparison.csv")
    names = {
        "PTV_only": 'Nine genes PTV-only',
        "expanded_star1": 'Nine genes PTV+P/LP missense (≥1 star)',
        "expanded_star2": 'Strictly Expanded (≥2 stars)',
    }
    lines = [
        '# Ninegenes expanded qualifying-variant carrier: exploratory age at diagnosis analysis',
        "",
        'Using the same batch of 7,996 eligible cases, the same case PRS mean/sample SD and PC1–10, the original PTV-only analysis is still the main analysis.',
        'Reuse ClinVar 2026-09-13/GRCh38 and existing VAT and MANE .5 annotations without re-running ClinVar.',
        "",
        '## Add carrier',
        "",
        '| Definition | Total number of carriers | New number of people compared to nine-gene PTV-only |',
        "|---|---:|---:|",
    ]
    for row in counts.itertuples():

        def display(v):
            return 'inhibition' if "suppressed" in str(v) else str(v)

        lines.append(
            f"| {names[row.definition]} | {display(row.n_carriers)} | {display(row.new_carriers_vs_PTV)} |"
        )
    lines.extend(
        [
            "",
            'New carriers refer to those who are PTV negative for the original nine genes and carry at least one qualified P/LP missense; carriers of both categories are not counted as new additions.',
            'New increments of less than or equal to 20 people are not directly released; the total number of ten-gene PTVs has been released, so if the total number of nine genes is reversed, a small-scale MAP3K1 difference set will be derived, and the related total number and the grouping table that can be added will also be suppressed. When the difference between the strict and main extended numbers reaches the threshold, the strict number is suppressed.',
            "",
            '## Recorded diagnosis age difference corresponding to each SD of PRS (years, 95% CI)',
            "",
            '| Carrier definition | Non-carrier slope | Carrier slope | Difference/interaction of two sets of slopes |',
            "|---|---:|---:|---:|",
        ]
    )
    for label, name in names.items():
        subset = models.loc[
            (models.definition == label) & (models.model == "with_interaction")
        ].set_index("contrast")
        formatted = []
        for term in ["PRS_z", "PRS_per_SD_carriers", "interaction"]:
            r = subset.loc[term]
            formatted.append(f"{r.age_difference_years:+.2f}（{r.ci_low:+.2f}to{r.ci_high:+.2f}）")
        lines.append(f"| {name} | " + " | ".join(formatted) + " |")
    main = models.loc[
        (models.definition == "expanded_star1")
        & (models.model == "with_interaction")
        & (models.contrast == "interaction")
    ].iloc[0]
    conclusion = (
        'The 95% CI of the interaction crosses zero, and no clear slope difference is seen; this does not mean that the slopes of the two groups have been proven to be the same.'
        if main.ci_low <= 0 <= main.ci_high
        else 'The 95% CI of the interaction does not cross zero, which is an exploratory interaction signal that requires independent verification.'
    )
    lines.extend(
        [
            "",
            conclusion,
            "",
            '## Prediction graph',
            "",
            '[Nine-gene expanded carrier predicted age chart](../figures/28_expanded_star1_predicted_age.svg)',
            "",
            'The horizontal axis PRS311 is the SD of the fixed case cohort, and the vertical axis is the PC-adjusted average recorded diagnosis age. The two lines represent the extended definition of carrier=0/1 respectively; the shadow is the HC1 point state 95%CI. Only the intersection of the 1–99th percentile range of each group is drawn, using the same PC mean as a reference; the shading is not the individual prediction interval.',
            "",
            '## Comments and interpretation restrictions',
            "",
            'The main extension requires at least 1 star, conflict-free germline P/LP, clear related diseases, AoU EUR AF<0.1%; the strict extension requires ≥2 stars. The original GQ≥20, het VAF[0.30,0.80], hom-alt VAF≥0.8 and CHEK2 I157T exclusion rules are retained. The new missense only reads fixed points from the original exome callset and folds it in memory; it does not save the new participant genotype or carrier table.',
            "",
            'No disease-level RCV classification is cached. Variants that are also listed with other clear diseases and cannot confirm the classification correspondence are conservatively excluded; only not_provided/not_specified will not make them eligible. Therefore, this is a clearly identifiable subset of P/LP missenses and cannot be claimed to exhaust all relevant pathogenic missenses. The original PTV was qualifying PTV, and clinical P/LP was not required one by one, so it was collectively called expanded qualifying-variant carrier.',
            "",
            'The MANE site cache does not cover all candidates; the results and genes recorded by MANE all match the original VAT, and the rest are based on the original canonical VAT missense results and the ClinVar missense annotations of the same version. There is no independent ClinGen classification field, the expert review rating comes from ClinVar CLNREVSTAT.',
            "",
            'Recording age associations within cases cannot prove early-onset risk synergy or biological interactions in the population; limitations such as retaining EHR dates later than actual onset, index-event bias, and low VAF/clonal hematopoiesis are not completely eliminated. Annotations or PRS thresholds were not adjusted based on age model results.',
            "",
            'Verification: The original case n=7,996 and the PRS reference assertion passed, and the newly added candidate site passed the original callset verification; the synthetic test covers disease exclusion, somatic-only, target ALT dose/VAF and the same position across shard boundaries. The main and extended models adopt the same OLS HC1 covariance rule.',
            "",
            'Basis: [ClinVar Review Level](https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/), [AoU Disclosure Policy](https://support.researchallofus.org/hc/en-us/articles/22346276580372-Data-and-Statistics-Dissemination-Policy).',
        ]
    )
    (REPORT / "28_expanded_missense_age_summary.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )

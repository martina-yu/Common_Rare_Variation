# Documentation: scripts/README.md.
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data/29_rowlands_extension"
CONFIG = ROOT / "config/rowlands_13gene_extension.yaml"
PTV = {"stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant"}
FIELDS = [
    "vid",
    "contig",
    "position",
    "ref_allele",
    "alt_allele",
    "gene_symbol",
    "transcript",
    "transcript_source",
    "consequence",
    "is_canonical_transcript",
    "gvs_eur_af",
    "gvs_eur_ac",
    "gvs_eur_an",
    "gnomad_nfe_af",
    "exon_number",
    "intron_number",
    "aa_change",
    "dna_change_in_transcript",
    "mane_select_name",
    "mane_plus_clinical_name",
    "LoF",
    "LoF_filter",
    "LoF_flags",
    "LoF_info",
]


def terms(value):
    return set(value.split(","))


def key(row):
    return f"{row.contig}:{row.position}:{row.ref_allele}:{row.alt_allele}"


def coordinates(panel):
    old = pd.read_csv(ROOT / "data/11_panel_gene_coordinates_hg38.csv")
    rows = old.loc[old.gene_symbol.isin(panel)].to_dict("records")
    for row in rows:
        row["coordinate_source"] = "prespecified frozen gene coordinates"
    missing = set(panel) - set(old.gene_symbol)
    gtf = ROOT / "data/29_ptv_audit/Homo_sapiens.GRCh38.110.gtf.gz"
    with gzip.open(gtf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            v = line.rstrip().split("\t")
            if v[2] != "gene":
                continue
            m = re.search(r'gene_name "([^";]+)"', v[8])
            if not m or m.group(1) not in missing:
                continue
            gene = m.group(1)
            rows.append(
                dict(
                    gene_symbol=gene,
                    chrom="chr" + v[0],
                    start_hg38=int(v[3]),
                    end_hg38=int(v[4]),
                    strand=1 if v[6] == "+" else -1,
                    query_region_start=int(v[3]) - 2000,
                    query_region_end=int(v[4]) + 2000,
                    coordinate_source="cached Ensembl GRCh38 release 110 GTF; region lookup only",
                )
            )
    assert {r["gene_symbol"] for r in rows} == set(panel) and len(rows) == len(panel)
    return pd.DataFrame(rows)


def normalize(frame, fasta):
    unique = frame[["orig_key", "contig", "position", "ref_allele", "alt_allele"]].drop_duplicates(
        "orig_key"
    )
    assert not unique.alt_allele.str.contains(",", regex=False).any()
    assert unique.ref_allele.str.fullmatch("[ACGTN]+", case=False).all()
    assert unique.alt_allele.str.fullmatch(
        "[ACGTN]+", case=False
    ).all(), "Unsupported ALT; inspect rather than silently dropping"
    keys = {}
    path = OUT / "normalization_input.vcf"
    with path.open("w") as f:
        f.write("##fileformat=VCFv4.2\n##reference=GRCh38\n")
        for line in Path(str(fasta) + ".fai").read_text().splitlines():
            chrom, length, *_ = line.split("\t")
            f.write(f"##contig=<ID={chrom},length={length}>\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for i, r in enumerate(unique.itertuples()):
            token = f"audit{i}"
            keys[token] = r.orig_key
            f.write(f"{r.contig}\t{r.position}\t{token}\t{r.ref_allele}\t{r.alt_allele}\t.\t.\t.\n")
    result = subprocess.run(
        ["bcftools", "norm", "-f", str(fasta), "-c", "e", "-m", "-any", "-Ov", str(path)],
        capture_output=True,
        text=True,
        check=True,
    )
    mapping = {}
    for line in result.stdout.splitlines():
        if line.startswith("#"):
            continue
        chrom, pos, token, ref, alt, *_ = line.split("\t")
        assert token in keys and keys[token] not in mapping
        mapping[keys[token]] = f"{chrom}:{pos}:{ref}:{alt}"
    assert len(mapping) == len(unique)
    path.unlink()
    return mapping

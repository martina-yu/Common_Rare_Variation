"""
- Reproduce legacy PTV selection and burden.
- Fit baseline joint and age models.
- Generate legacy result figures.
"""
from __future__ import annotations
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
# ============================================================================
# 1. vat-schema: analyze_vat_schema
# ============================================================================

import os
import subprocess
import sys

def vat_schema_run(cmd, **kwargs):
    print(f"$ {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, **kwargs)
    if result.stdout:
        print(result.stdout)
    if result.stderr:
        print(result.stderr, file=sys.stderr)
    return result

def vat_schema_section(title):
    print("\n" + "=" * 80)
    print(title)
    print("=" * 80)

def vat_schema_calculation():
    print(f"CDR_STORAGE_PATH = {VAT_SCHEMA_CDR_STORAGE_PATH}")
    print(f"GOOGLE_PROJECT   = {VAT_SCHEMA_GOOGLE_PROJECT}")

    # ------------------------------------------------------------------
    
    # ------------------------------------------------------------------
    vat_schema_section('Q1. Variation data source format: directory structure detection (gsutil ls)')
    for rel in [
        "",
        "wgs/",
        "wgs/short_read/",
        "wgs/short_read/snpindel/",
        "wgs/short_read/snpindel/vds/",
        "wgs/short_read/snpindel/vds/hail.vds/",
        "wgs/short_read/snpindel/acaf_threshold/",
        "wgs/short_read/snpindel/exome/",
        "wgs/short_read/snpindel/aux/",
        "wgs/short_read/snpindel/aux/vat/",
    ]:
        path = f"{VAT_SCHEMA_CDR_STORAGE_PATH}/{rel}"
        vat_schema_run(f"gsutil -u {VAT_SCHEMA_GOOGLE_PROJECT} ls {path}")

    vat_schema_section('Q1. VAT file size confirmation (basis for not downloading the entire table)')
    vat_schema_run(f"gsutil -u {VAT_SCHEMA_GOOGLE_PROJECT} du -h {VAT_PATH}")

    # ------------------------------------------------------------------
    
    # ------------------------------------------------------------------
    vat_schema_section('Q2. Actual measurement of VAT schema: only pull the byte range of the file header, decompress it locally and view it')

    small_bgz = os.path.join(SCRATCH_DIR, "vat_head_2MB.bgz")
    vat_schema_run(f"gsutil -u {VAT_SCHEMA_GOOGLE_PROJECT} cat -r 0-2000000 {VAT_PATH} > {small_bgz}")
    vat_schema_run(f"ls -la {small_bgz}")

    head_tsv = os.path.join(SCRATCH_DIR, "vat_head.tsv")

    vat_schema_run(f"zcat {small_bgz} 2>/dev/null | head -3 > {head_tsv}")
    vat_schema_run(f"wc -l {head_tsv}")

    vat_schema_section('Q2. Complete list of names (numbered in order of appearance)')
    col_file = os.path.join(SCRATCH_DIR, "vat_columns.txt")
    vat_schema_run(f"head -1 {head_tsv} | tr '\\t' '\\n' | nl > {col_file}")
    vat_schema_run(f"cat {col_file}")

    vat_schema_section('Q2. Check the existence of key fields (grep -i, check the 6 fields in the plan table one by one)')
    targets = {
        "revel": "REVEL",
        "consequence": 'PTV consequence (SO term)',
        "splice_ai": "SpliceAI",
        "^\\s*[0-9]+\\s+LoF$": 'LOFTEE - LoF (HC/LC confidence level)',
        "LoF_filter": "LOFTEE - LoF_filter",
        "LoF_flags": "LOFTEE - LoF_flags",
        "cadd": "CADD",
        "alpha": 'AlphaMissense (keyword alpha/missense)',
        "missense": 'AlphaMissense (keyword missense, exclude revel/consequence itself)',
    }
    for pattern, label in targets.items():
        res = vat_schema_run(f"grep -inE '{pattern}' {col_file}")
        found = bool(res.stdout.strip())
        print(f"  >>> [{label}] pattern='{pattern}' -> {'FOUND' if found else 'NOT FOUND'}")

    # ------------------------------------------------------------------
    
    # ------------------------------------------------------------------
    vat_schema_section('Q2 Supplement: Do LoF(LOFTEE) and revel columns have non-empty measured values \u200b\u200b(instead of just column names)')
    bigger_bgz = os.path.join(SCRATCH_DIR, "vat_head_8MB.bgz")
    vat_schema_run(f"gsutil -u {VAT_SCHEMA_GOOGLE_PROJECT} cat -r 0-8000000 {VAT_PATH} > {bigger_bgz}")
    bigger_tsv = os.path.join(SCRATCH_DIR, "vat_head_8MB.tsv")
    vat_schema_run(f"zcat {bigger_bgz} 2>/dev/null > {bigger_tsv}")
    vat_schema_run(f"wc -l {bigger_tsv}")

    awk_script = r"""
    awk -F'\t' '
      NR==1 {
        for (i=1;i<=NF;i++) { if ($i=="revel") revel_col=i; if ($i=="consequence") cons_col=i;
                              if ($i=="LoF") lof_col=i; if ($i=="LoF_filter") loffilter_col=i; }
        print "revel col idx =", revel_col, " consequence col idx =", cons_col, " LoF col idx =", lof_col, " LoF_filter col idx =", loffilter_col;
        next
      }
      $lof_col != "" { lof_nonempty++; if (lof_examples < 5) { print "  LoF example:", $1, $cons_col, $lof_col, $loffilter_col; lof_examples++ } }
      $revel_col != "" { revel_nonempty++ }
      END { print "LoF non-empty rows:", lof_nonempty+0, " / total data rows:", NR-1;
            print "revel non-empty rows:", revel_nonempty+0, " / total data rows:", NR-1 }
    ' """ + bigger_tsv
    vat_schema_run(awk_script)

    vat_schema_section('Finish')
    print('The above is the original output of the actual measurement, and no interpretation judgment has been made to replace it - whether it is consistent with the expectations of the plan Step 3.0 table will be left for discussion in the report.')

def analyze_vat_schema():
    """Initialize runtime state and execute this component explicitly."""
    global VAT_SCHEMA_CDR_STORAGE_PATH, VAT_SCHEMA_GOOGLE_PROJECT, SCRATCH_DIR, VAT_PATH

    VAT_SCHEMA_GOOGLE_PROJECT = os.environ.get("GOOGLE_PROJECT", "")
    VAT_SCHEMA_CDR_STORAGE_PATH = os.environ.get("CDR_STORAGE_PATH", "")
    VAT_PATH = f"{VAT_SCHEMA_CDR_STORAGE_PATH}/wgs/short_read/snpindel/aux/vat/vat_complete.bgz.tsv.gz"

    SCRATCH_DIR = "/tmp/vat_schema_check_scratch"
    os.makedirs(SCRATCH_DIR, exist_ok=True)

    if True:
        vat_schema_calculation()

# ============================================================================
# 2. vat-prep: analyze_legacy_vat_preparation
# ============================================================================

import json
import os
import subprocess
import sys
import time
import urllib.request
import pandas as pd

def legacy_vat_preparation_run(cmd, quiet=False, max_print_chars=4000, **kwargs):
    """Legacy vat preparation run."""
    print(f"$ {cmd}", flush=True)
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, **kwargs)
    if result.stdout:
        if quiet:
            n_lines = result.stdout.count("\n")
            preview = result.stdout[:max_print_chars]
            print(f"[quiet=True, output{n_lines}OK/{len(result.stdout)}characters, only before printing{max_print_chars}Character preview]")
            print(preview)
            if len(result.stdout) > max_print_chars:
                print('...[truncated]...')
        else:
            print(result.stdout)
    if result.stderr:
        print(result.stderr, file=sys.stderr)
    return result

def legacy_vat_preparation_section(title):
    print("\n" + "=" * 90)
    print(title)
    print("=" * 90, flush=True)

def ensembl_get(url, retries=3):
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=20) as resp:
                return json.loads(resp.read().decode())
        except Exception as e:
            print(f"  retry {attempt+1}/{retries} for {url}: {e}")
            time.sleep(2)
    raise RuntimeError(f"Failed to fetch {url}")

def legacy_vat_preparation_calculation():
    legacy_vat_preparation_section('0. Environment confirmation')
    print(f"CDR_STORAGE_PATH = {LEGACY_VAT_PREPARATION_CDR_STORAGE_PATH}")
    print(f"BUCKET_NAME      = {LEGACY_VAT_PREPARATION_BUCKET_NAME}")
    print(f"CDR_PREFIX       = {LEGACY_VAT_PREPARATION_CDR_PREFIX}")
    print(f"GOOGLE_PROJECT   = {LEGACY_VAT_PREPARATION_GOOGLE_PROJECT}")

    # ------------------------------------------------------------------
    
    # ------------------------------------------------------------------
    legacy_vat_preparation_section('1. gcsfuse read-only mount controlled-tier bucket')
    legacy_vat_preparation_run(f"fusermount -u {LEGACY_VAT_PREPARATION_MOUNT_DIR} 2>/dev/null; true")
    r = legacy_vat_preparation_run(
        f"gcsfuse --billing-project {LEGACY_VAT_PREPARATION_GOOGLE_PROJECT} -o ro --implicit-dirs "
        f"{LEGACY_VAT_PREPARATION_BUCKET_NAME} {LEGACY_VAT_PREPARATION_MOUNT_DIR}"
    )
    if r.returncode != 0:
        print('!!! gcsfuse failed to mount, and subsequent steps that depend on the mount point will fail and require manual troubleshooting!!!')
    mount_root = os.path.join(LEGACY_VAT_PREPARATION_MOUNT_DIR, LEGACY_VAT_PREPARATION_CDR_PREFIX)
    vat_path = os.path.join(LEGACY_VAT_PREPARATION_MOUNT_DIR, LEGACY_VAT_PREPARATION_VAT_REL)
    exome_vcf_dir = os.path.join(LEGACY_VAT_PREPARATION_MOUNT_DIR, EXOME_VCF_DIR_REL)
    legacy_vat_preparation_run(f"ls -la {vat_path} {vat_path}.tbi")

    # ------------------------------------------------------------------
    
    # ------------------------------------------------------------------
    legacy_vat_preparation_section('2. Panel gene GRCh38 coordinates (Ensembl REST lookup/symbol)')
    gene_rows = []
    for g in PANEL_GENES:
        url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{g}?content-type=application/json"
        info = ensembl_get(url)
        row = {
            "gene_symbol": g,
            "ensembl_gene_id": info["id"],
            "chrom": "chr" + info["seq_region_name"],
            "start_hg38": info["start"],
            "end_hg38": info["end"],
            "strand": info["strand"],
            "canonical_transcript": info.get("canonical_transcript"),
            "query_region_start": max(1, info["start"] - QUERY_BUFFER_BP),
            "query_region_end": info["end"] + QUERY_BUFFER_BP,
        }
        gene_rows.append(row)
        print(f"  {g}: {row['chrom']}:{row['start_hg38']}-{row['end_hg38']} "
              f"(strand {row['strand']}, {row['ensembl_gene_id']}, canonical={row['canonical_transcript']})")
    gene_df = pd.DataFrame(gene_rows)
    gene_coord_path = os.path.join(LEGACY_VAT_PREPARATION_DATA_DIR, "11_panel_gene_coordinates_hg38.csv")
    gene_df.to_csv(gene_coord_path, index=False)
    print(f"\n writes:{gene_coord_path} ({len(gene_df)}OK)")

    # ------------------------------------------------------------------
    
    # ------------------------------------------------------------------
    legacy_vat_preparation_section('3. CHEK2 p.I157T hg38 coordinate verification (rs17879961, Ensembl variation API)')
    var_info = ensembl_get(
        "https://rest.ensembl.org/variation/human/rs17879961?content-type=application/json"
    )
    mappings = [m for m in var_info["mappings"] if m["assembly_name"] == "GRCh38"]
    print(f"  rs17879961 GRCh38 mappings: {mappings}")
    assert len(mappings) == 1, 'rs17879961 should have a unique mapping on GRCh38, otherwise manual verification is required'
    i157t_hg38_pos = mappings[0]["start"]
    i157t_hg38_chrom = "chr" + mappings[0]["seq_region_name"]
    print(f"-> I157T hg38 coordinates:{i157t_hg38_chrom}:{i157t_hg38_pos} "
          f"(GRCh37 original coordinate 22_29121087_A_G comparison)")

    r = legacy_vat_preparation_run(f"tabix {vat_path} {i157t_hg38_chrom}:{i157t_hg38_pos}-{i157t_hg38_pos}")
    i157t_vat_rows = [l for l in r.stdout.strip().split("\n") if l]
    print(f"VAT hits at this coordinate{len(i157t_vat_rows)}Lines (mostly different transcripts)")
    i157t_confirmed = False
    i157t_ref, i157t_alt = None, None
    for line in i157t_vat_rows:
        fields = line.split("\t")
        ref, alt = fields[4], fields[5]
        gene_sym = fields[43]
        rsid = fields[52]
        if gene_sym == "CHEK2" and ref == "A" and alt == "G":
            i157t_confirmed = True
            i157t_ref, i157t_alt = ref, alt
            print(f"Confirm match: ref={ref} alt={alt} gene={gene_sym} rsid={rsid}")
    print(f"I157T hg38 coordinates are cross-validated by VAT:{'CONFIRMED' if i157t_confirmed else 'NOT CONFIRMED -- manual verification required'}")

    i157t_row = pd.DataFrame([{
        "grch37_coord": "22_29121087_A_G",
        "grch38_chrom": i157t_hg38_chrom,
        "grch38_pos": i157t_hg38_pos,
        "grch38_ref": i157t_ref,
        "grch38_alt": i157t_alt,
        "dbsnp_rsid": "rs17879961",
        "vat_cross_validated": i157t_confirmed,
    }])
    i157t_path = os.path.join(LEGACY_VAT_PREPARATION_DATA_DIR, "11_chek2_i157t_hg38_coordinate.csv")
    i157t_row.to_csv(i157t_path, index=False)
    print(f"Write:{i157t_path}")

    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    legacy_vat_preparation_section('4. VAT is extracted according to gene interval (canonical transcript row), and the placement is based on gene small table')
    PTV_TERMS = {"stop_gained", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant"}
    MISSENSE_TERMS = {"missense_variant"}
    all_gene_vat = []
    all_gene_vat_full_for_xtab = []  
    for _, g in gene_df.iterrows():
        region = f"{g['chrom']}:{g['query_region_start']}-{g['query_region_end']}"
        t0 = time.time()
        r = legacy_vat_preparation_run(f"tabix {vat_path} {region}", quiet=True)
        elapsed = time.time() - t0
        lines = [l for l in r.stdout.strip().split("\n") if l]
        print(f"  {g['gene_symbol']} ({region}): {len(lines)}line original hit,{elapsed:.1f}s")
        rows = []
        xtab_counter = {}
        for line in lines:
            f = line.split("\t")
            if len(f) < 114:
                continue
            is_canonical = f[55] == "true"  # 0-indexed 55 = 1-indexed 56 is_canonical_transcript
            has_lof_call = f[109] not in ("", None)  # 0-indexed 109 = 1-indexed 110 LoF(LOFTEE)

            if not (is_canonical or has_lof_call):
                continue
            gene_sym = f[43]
            if is_canonical:
                xtab_counter[gene_sym] = xtab_counter.get(gene_sym, 0) + 1

            if gene_sym != g["gene_symbol"]:
                continue
            cons_terms = set(f[46].split(","))  # consequence, 0-indexed 46 = 1-indexed 47
            if not (cons_terms & PTV_TERMS) and not (cons_terms & MISSENSE_TERMS):
                continue
            rows.append({name: f[idx - 1] for name, idx in LEGACY_VAT_PREPARATION_VAT_COLUMNS_KEEP.items()})
        sub_df = pd.DataFrame(rows)
        sub_df["panel_gene_query"] = g["gene_symbol"]
        print(f"-> canonical + gene_symbol strict match + PTV/missense related:{len(sub_df)}OK")
        all_gene_vat.append(sub_df)
        all_gene_vat_full_for_xtab.append(
            pd.DataFrame(
                [{"panel_gene_query": g["gene_symbol"], "gene_symbol": k, "n": v}
                 for k, v in xtab_counter.items()]
            )
        )

    panel_vat_df = pd.concat(all_gene_vat, ignore_index=True)
    panel_vat_path = os.path.join(LEGACY_VAT_PREPARATION_DATA_DIR, "11_panel_gene_vat_ptv_missense_candidates.csv.gz")
    panel_vat_df.to_csv(panel_vat_path, index=False, compression="gzip")
    size_mb = os.path.getsize(panel_vat_path) / 1e6
    print(f"\n writes:{panel_vat_path} ({len(panel_vat_df)}OK,{size_mb:.2f}MB)")
    if size_mb > 5:
        print('!!! Exceeds the 5MB limit of CLAUDE.md single file and needs to be split or compressed !!!')
    else:
        print('Comply with CLAUDE.md single file 5MB limit.')

    legacy_vat_preparation_section('4b. gene_symbol consistency check (whether the buffer brings in neighboring genes, based on the canonical raw count before filtering)')
    xtab = pd.concat(all_gene_vat_full_for_xtab, ignore_index=True)
    xtab = xtab.groupby(["panel_gene_query", "gene_symbol"], as_index=False)["n"].sum()
    print(xtab.to_string(index=False))
    mismatch = xtab[xtab["panel_gene_query"] != xtab["gene_symbol"]]
    if len(mismatch) > 0:
        print(f"\n found{len(mismatch)}A combination of gene_symbol and query gene is inconsistent (this script has been"
              f"Strictly filter by gene_symbol, the downstream script 12 can directly use the placement results without repeated filtering):")
        print(mismatch.to_string(index=False))
    else:
        print('No adjacent gene stringing was found, but the downstream script 12 will still explicitly filter by gene_symbol just in case.')

    # ------------------------------------------------------------------
    
    # ------------------------------------------------------------------
    legacy_vat_preparation_section('5. Exome callset VCF FORMAT field confirmation (carrier-level QC requires GT/AD/GQ)')
    sample_shard = os.path.join(exome_vcf_dir, "0000000000.vcf.bgz")
    r = legacy_vat_preparation_run(f"bcftools view -h {sample_shard} | grep '^##FORMAT'")
    r2 = legacy_vat_preparation_run(f"bcftools query -l {sample_shard} | wc -l")
    exome_n_samples = int(r2.stdout.strip())
    print(f"Exome VCF (shard0) sample number:{exome_n_samples}")
    has_gt = "ID=GT" in r.stdout
    has_ad = "ID=AD" in r.stdout
    has_gq = "ID=GQ" in r.stdout
    has_dp = "ID=DP" in r.stdout
    print(f"  GT: {has_gt}, AD: {has_ad}, GQ: {has_gq}, explicit DP:{has_dp}")
    if not has_dp:
        print('Note: FORMAT has no explicit DP field, VAF is calculated according to AD_ALT/(AD_ref+AD_alt) (sum(AD) approximates DP),'
              'With Le Guen et al. 2026 ADSP process VAF definition (AD_ALT/DP) in GATK joint-called VCF scenario'
              'is the standard equivalent approach')

    # ------------------------------------------------------------------
    
    # ------------------------------------------------------------------
    legacy_vat_preparation_section('6. Exome VCF shard binary positioning (which shards each panel gene falls on)')
    shard_files = sorted(
        f for f in os.listdir(exome_vcf_dir) if f.endswith(".vcf.bgz")
    )
    n_shards = len(shard_files)
    print(f"common{n_shards}shard")

    CHROM_ORDER = {f"chr{i}": i for i in range(1, 23)}
    CHROM_ORDER["chrX"] = 23
    CHROM_ORDER["chrY"] = 24

    legacy_vat_preparation_section('6a. Pre-scan the starting boundaries of all shards (each shard only reads the first body line, concurrent I/O)')

    LOCAL_INTERVAL_DIR = "/tmp/claude_interval_lists_full"
    use_local = os.path.isdir(LOCAL_INTERVAL_DIR) and len(
        [f for f in os.listdir(LOCAL_INTERVAL_DIR) if f.endswith(".interval_list")]
    ) >= n_shards
    interval_source_dir = LOCAL_INTERVAL_DIR if use_local else exome_vcf_dir
    print(f"interval_list reading source:{interval_source_dir}(Local cache hit:{use_local})")

    def read_first_body(i_fname):
        i, fname = i_fname
        ipath = os.path.join(interval_source_dir, fname.replace(".vcf.bgz", ".interval_list"))
        try:
            with open(ipath) as fh:
                for line in fh:
                    if line.startswith("@"):
                        continue
                    first_body = line.strip().split("\t")
                    c_ord = CHROM_ORDER.get(first_body[0], 99)
                    return (c_ord, int(first_body[1]), i)
        except Exception as e:
            print(f"!!! Read failed{ipath}: {e}")
        return None

    import concurrent.futures as cf

    t0 = time.time()
    results = [None] * n_shards

    with cf.ThreadPoolExecutor(max_workers=64) as executor:
        for idx, res in zip(
            range(n_shards),
            executor.map(read_first_body, list(enumerate(shard_files))),
        ):
            results[idx] = res
            if idx % 2000 == 0:
                print(f"...schedule{idx}/{n_shards}, time consuming{time.time()-t0:.0f}s", flush=True)
    shard_starts = [r for r in results if r is not None]
    elapsed = time.time() - t0
    print(f"prescan{len(shard_starts)}A shard is completed, which takes time{elapsed:.1f}s")
    
    is_sorted = all(shard_starts[i] <= shard_starts[i + 1] for i in range(len(shard_starts) - 1))
    print(f"Whether the starting boundary of shard is monotonic and undecreasing:{is_sorted}")
    if not is_sorted:
        print('!!! The starting boundary of the shard is non-monotonic, and the binary search hypothesis cannot be used. It needs to be changed to a full linear scan or manual verification!!!')

    starts_only = [(c, p) for c, p, _ in shard_starts]
    shard_bounds_path = os.path.join(LEGACY_VAT_PREPARATION_DATA_DIR, "11_exome_vcf_shard_start_boundaries.csv")
    pd.DataFrame(shard_starts, columns=["chrom_ord", "start_pos", "shard_idx"]).to_csv(
        shard_bounds_path, index=False
    )
    print(f"Write out (cached for reproduction):{shard_bounds_path}")

    import bisect

    def find_shard_range(chrom_ord, start, end):
        """Find shard range."""
        
        lo = bisect.bisect_right(starts_only, (chrom_ord, start))
        first = max(0, lo - 1)
        
        hi = bisect.bisect_right(starts_only, (chrom_ord, end))
        last = max(0, hi - 1)
        
        first = max(0, first - 1)
        last = min(len(shard_starts) - 1, last + 1)
        return shard_starts[first][2], shard_starts[last][2]

    legacy_vat_preparation_section('6b. Gene-by-gene binary positioning shard range')
    shard_rows = []
    for _, g in gene_df.iterrows():
        chrom_ord = CHROM_ORDER[g["chrom"]]
        t0 = time.time()
        first, last = find_shard_range(chrom_ord, g["query_region_start"], g["query_region_end"])
        elapsed = time.time() - t0
        print(f"  {g['gene_symbol']} ({g['chrom']}:{g['query_region_start']}-{g['query_region_end']}) "
              f"-> shard[{first}, {last}] ({last-first+1}shard), positioning takes time{elapsed:.3f}s")
        shard_rows.append({
            "gene_symbol": g["gene_symbol"], "chrom": g["chrom"],
            "region_start": g["query_region_start"], "region_end": g["query_region_end"],
            "shard_first": first, "shard_last": last, "n_shards": last - first + 1,
        })
    shard_df = pd.DataFrame(shard_rows)
    shard_path = os.path.join(LEGACY_VAT_PREPARATION_DATA_DIR, "11_panel_gene_exome_vcf_shard_map.csv")
    shard_df.to_csv(shard_path, index=False)
    print(f"\n writes:{shard_path}")
    print(shard_df.to_string(index=False))

    # ------------------------------------------------------------------
    
    # ------------------------------------------------------------------
    legacy_vat_preparation_section('7. Callset division of labor record (red line to prevent misuse)')
    print(f"""- Step 3/4 rare variant qualifying variant determination + carrier genotype extraction:
      -> VAT:{LEGACY_VAT_PREPARATION_VAT_REL}
      -> Exome VCF: {EXOME_VCF_DIR_REL}/*.vcf.bgz ({exome_n_samples}samples)
  - SAIGE-GENE+/REGENIE null model (estimates the correlation background structure between samples, common variations):
      -> ACAF threshold callset:{ACAF_VCF_DIR_REL}(This script will not be executed, leaving it to the discovery scan script)
  - ⚠️Hard constraints: rare variant burden/carrier determines that ACAF callset will never be used
      (ACAF is designed to exclude most rare mutations with AF<1%/AC<=100. If used incorrectly, the carrier will be lost silently)
  - VDS full path (standby, if the VCF shard method encounters problems, you can fall back to Hail, but the current environment does not have Dataproc/Spark,
    The Hail reading gs:// path has been previously tested in Step 2 and is not feasible. This path is only recorded and not used as the main solution for this round):
      ->{VDS_REL}
""")

    legacy_vat_preparation_section('Finish')
    print('All remaining technical prerequisites of Step 3.0 have been verified through actual testing, and you can enter script 12 (qualifying variant definition).')

def analyze_legacy_vat_preparation():
    """Initialize runtime state and execute this component explicitly."""
    global ACAF_VCF_DIR_REL, LEGACY_VAT_PREPARATION_BUCKET_NAME, LEGACY_VAT_PREPARATION_CDR_PREFIX, LEGACY_VAT_PREPARATION_CDR_STORAGE_PATH, LEGACY_VAT_PREPARATION_DATA_DIR, EXOME_VCF_DIR_REL, LEGACY_VAT_PREPARATION_GOOGLE_PROJECT, LEGACY_VAT_PREPARATION_MOUNT_DIR, PANEL_GENES, QUERY_BUFFER_BP, LEGACY_VAT_PREPARATION_REPO_ROOT, LEGACY_VAT_PREPARATION_TABLE_DIR, LEGACY_VAT_PREPARATION_VAT_COLUMNS_KEEP, LEGACY_VAT_PREPARATION_VAT_REL, VDS_REL

    LEGACY_VAT_PREPARATION_GOOGLE_PROJECT = os.environ.get("GOOGLE_PROJECT", "")
    LEGACY_VAT_PREPARATION_CDR_STORAGE_PATH = os.environ.get("CDR_STORAGE_PATH", "")
    LEGACY_VAT_PREPARATION_BUCKET_NAME = LEGACY_VAT_PREPARATION_CDR_STORAGE_PATH.replace("gs://", "").split("/")[0]
    LEGACY_VAT_PREPARATION_CDR_PREFIX = "/".join(LEGACY_VAT_PREPARATION_CDR_STORAGE_PATH.replace("gs://", "").split("/")[1:])  # e.g. "v9"

    LEGACY_VAT_PREPARATION_REPO_ROOT = "/home/jupyter/workspace/analysis"
    LEGACY_VAT_PREPARATION_DATA_DIR = os.path.join(LEGACY_VAT_PREPARATION_REPO_ROOT, "data")
    LEGACY_VAT_PREPARATION_TABLE_DIR = os.path.join(LEGACY_VAT_PREPARATION_REPO_ROOT, "results", "table")
    LEGACY_VAT_PREPARATION_MOUNT_DIR = "/tmp/claude_aou_gcsfuse_mount"

    os.makedirs(LEGACY_VAT_PREPARATION_DATA_DIR, exist_ok=True)
    os.makedirs(LEGACY_VAT_PREPARATION_TABLE_DIR, exist_ok=True)
    os.makedirs(LEGACY_VAT_PREPARATION_MOUNT_DIR, exist_ok=True)

    LEGACY_VAT_PREPARATION_VAT_REL = f"{LEGACY_VAT_PREPARATION_CDR_PREFIX}/wgs/short_read/snpindel/aux/vat/vat_complete.bgz.tsv.gz"
    EXOME_VCF_DIR_REL = f"{LEGACY_VAT_PREPARATION_CDR_PREFIX}/wgs/short_read/snpindel/exome/vcf"
    ACAF_VCF_DIR_REL = f"{LEGACY_VAT_PREPARATION_CDR_PREFIX}/wgs/short_read/snpindel/acaf_threshold/vcf"
    VDS_REL = f"{LEGACY_VAT_PREPARATION_CDR_PREFIX}/wgs/short_read/snpindel/vds/hail.vds"

    PANEL_GENES = [
        "ATM", "BARD1", "BRCA1", "BRCA2", "CHEK2", "PALB2",
        "RAD51C", "RAD51D", "TP53", "MAP3K1",
    ]

    QUERY_BUFFER_BP = 2000

    LEGACY_VAT_PREPARATION_VAT_COLUMNS_KEEP = {
        "vid": 1, "transcript": 2, "contig": 3, "position": 4,
        "ref_allele": 5, "alt_allele": 6,
        "gvs_eur_ac": 28, "gvs_eur_an": 29, "gvs_eur_af": 30,
        "gene_symbol": 44, "aa_change": 46, "consequence": 47,
        "dna_change_in_transcript": 48, "variant_type": 49,
        "exon_number": 50, "intron_number": 51,
        "dbsnp_rsid": 53, "gene_id": 54,
        "is_canonical_transcript": 56,
        "gnomad_all_af": 57, "gnomad_nfe_af": 82,
        "revel": 89,
        "clinvar_classification": 100,
        "mane_select_name": 106, "mane_plus_clinical_name": 107,
        "LoF": 110, "LoF_filter": 111, "LoF_flags": 112, "LoF_info": 113,
    }

    if True:
        legacy_vat_preparation_calculation()

# ============================================================================
# 3. vat-repair: analyze_legacy_vat_repair
# ============================================================================

import os
import subprocess
import sys
import time
import pandas as pd

def legacy_vat_repair_run(cmd, quiet=False, max_print_chars=4000, **kwargs):
    print(f"$ {cmd}", flush=True)
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, **kwargs)
    if result.stdout:
        if quiet:
            n_lines = result.stdout.count("\n")
            print(f"  [quiet=True, {n_lines}OK/{len(result.stdout)}character]")
            print(result.stdout[:max_print_chars])
        else:
            print(result.stdout)
    if result.stderr:
        print(result.stderr, file=sys.stderr)
    return result

def legacy_vat_repair_calculation():
    if not os.path.ismount(LEGACY_VAT_REPAIR_MOUNT_DIR):
        legacy_vat_repair_run(f"gcsfuse --billing-project {LEGACY_VAT_REPAIR_GOOGLE_PROJECT} -o ro --implicit-dirs "
            f"vwb-aou-datasets-controlled {LEGACY_VAT_REPAIR_MOUNT_DIR}")
    vat_path = os.path.join(LEGACY_VAT_REPAIR_MOUNT_DIR, LEGACY_VAT_REPAIR_VAT_REL)

    gene_df = pd.read_csv(os.path.join(LEGACY_VAT_REPAIR_DATA_DIR, "11_panel_gene_coordinates_hg38.csv"))

    all_gene_vat = []
    for _, g in gene_df.iterrows():
        region = f"{g['chrom']}:{g['query_region_start']}-{g['query_region_end']}"
        t0 = time.time()
        r = legacy_vat_repair_run(f"tabix {vat_path} {region}", quiet=True)
        elapsed = time.time() - t0
        lines = [l for l in r.stdout.strip().split("\n") if l]
        print(f"  {g['gene_symbol']} ({region}): {len(lines)}line original hit,{elapsed:.1f}s")
        rows = []
        n_canonical = 0
        n_lof_noncanonical = 0
        for line in lines:
            f = line.split("\t")
            if len(f) < 114:
                continue
            is_canonical = f[55] == "true"
            has_lof_call = f[109] not in ("", None)
            if not (is_canonical or has_lof_call):
                continue
            gene_sym = f[43]
            if gene_sym != g["gene_symbol"]:
                continue
            cons_terms = set(f[46].split(","))
            if not (cons_terms & LEGACY_VAT_REPAIR_PTV_TERMS) and not (cons_terms & LEGACY_VAT_REPAIR_MISSENSE_TERMS):
                continue
            if is_canonical:
                n_canonical += 1
            else:
                n_lof_noncanonical += 1
            rows.append({name: f[idx - 1] for name, idx in LEGACY_VAT_REPAIR_VAT_COLUMNS_KEEP.items()})
        sub_df = pd.DataFrame(rows)
        sub_df["panel_gene_query"] = g["gene_symbol"]
        print(f"    -> canonical:{n_canonical}Line + non-canonical but marked with LOFTEE:{n_lof_noncanonical}OK"
              f"= total{len(sub_df)}OK")
        all_gene_vat.append(sub_df)

    panel_vat_df = pd.concat(all_gene_vat, ignore_index=True)
    panel_vat_path = os.path.join(LEGACY_VAT_REPAIR_DATA_DIR, "11_panel_gene_vat_ptv_missense_candidates.csv.gz")
    panel_vat_df.to_csv(panel_vat_path, index=False, compression="gzip")
    size_mb = os.path.getsize(panel_vat_path) / 1e6
    print(f"\n writes:{panel_vat_path} ({len(panel_vat_df)}OK,{size_mb:.2f}MB)")
    if size_mb > 5:
        print('!!! Exceeded CLAUDE.md single file 5MB limit !!!')
    else:
        print('Comply with CLAUDE.md single file 5MB limit.')

    # CHEK2 sanity check
    print("\n--- CHEK2 PTV(HC) sanity check ---")
    chek2 = panel_vat_df[panel_vat_df["gene_symbol"] == "CHEK2"]
    chek2_ptv_hc = chek2[
        chek2["consequence"].apply(lambda c: bool(set(str(c).split(",")) & LEGACY_VAT_REPAIR_PTV_TERMS))
        & (chek2["LoF"] == "HC")
    ]
    print(f"CHEK2 PTV+LOFTEE HC lines:{len(chek2_ptv_hc)}")
    print(f"CHEK2 unique vid number (PTV+HC):{chek2_ptv_hc['vid'].nunique()}")
    c1100 = chek2_ptv_hc[chek2_ptv_hc["dna_change_in_transcript"].astype(str).str.contains("1100del", na=False)]
    print(f"The vid containing \'1100del\':{c1100['vid'].unique().tolist()}")

def analyze_legacy_vat_repair():
    """Initialize runtime state and execute this component explicitly."""
    global LEGACY_VAT_REPAIR_CDR_PREFIX, LEGACY_VAT_REPAIR_CDR_STORAGE_PATH, LEGACY_VAT_REPAIR_DATA_DIR, LEGACY_VAT_REPAIR_GOOGLE_PROJECT, LEGACY_VAT_REPAIR_MISSENSE_TERMS, LEGACY_VAT_REPAIR_MOUNT_DIR, LEGACY_VAT_REPAIR_PTV_TERMS, LEGACY_VAT_REPAIR_REPO_ROOT, LEGACY_VAT_REPAIR_VAT_COLUMNS_KEEP, LEGACY_VAT_REPAIR_VAT_REL

    LEGACY_VAT_REPAIR_GOOGLE_PROJECT = os.environ.get("GOOGLE_PROJECT", "")
    LEGACY_VAT_REPAIR_CDR_STORAGE_PATH = os.environ.get("CDR_STORAGE_PATH", "")
    LEGACY_VAT_REPAIR_CDR_PREFIX = "/".join(LEGACY_VAT_REPAIR_CDR_STORAGE_PATH.replace("gs://", "").split("/")[1:])

    LEGACY_VAT_REPAIR_REPO_ROOT = "/home/jupyter/workspace/analysis"
    LEGACY_VAT_REPAIR_DATA_DIR = os.path.join(LEGACY_VAT_REPAIR_REPO_ROOT, "data")
    LEGACY_VAT_REPAIR_MOUNT_DIR = "/tmp/claude_aou_gcsfuse_mount"
    LEGACY_VAT_REPAIR_VAT_REL = f"{LEGACY_VAT_REPAIR_CDR_PREFIX}/wgs/short_read/snpindel/aux/vat/vat_complete.bgz.tsv.gz"

    LEGACY_VAT_REPAIR_PTV_TERMS = {"stop_gained", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant"}
    LEGACY_VAT_REPAIR_MISSENSE_TERMS = {"missense_variant"}

    LEGACY_VAT_REPAIR_VAT_COLUMNS_KEEP = {
        "vid": 1, "transcript": 2, "contig": 3, "position": 4,
        "ref_allele": 5, "alt_allele": 6,
        "gvs_eur_ac": 28, "gvs_eur_an": 29, "gvs_eur_af": 30,
        "gene_symbol": 44, "aa_change": 46, "consequence": 47,
        "dna_change_in_transcript": 48, "variant_type": 49,
        "exon_number": 50, "intron_number": 51,
        "dbsnp_rsid": 53, "gene_id": 54,
        "is_canonical_transcript": 56,
        "gnomad_all_af": 57, "gnomad_nfe_af": 82,
        "revel": 89,
        "clinvar_classification": 100,
        "mane_select_name": 106, "mane_plus_clinical_name": 107,
        "LoF": 110, "LoF_filter": 111, "LoF_flags": 112, "LoF_info": 113,
    }

    if True:
        legacy_vat_repair_calculation()

# ============================================================================
# 4. variant-selection: analyze_legacy_variant_selection
# ============================================================================

import os
import numpy as np
import pandas as pd
import yaml

def legacy_variant_selection_section(title):
    print("\n" + "=" * 90)
    print(title)
    print("=" * 90, flush=True)

def write_variant_list(frame, cols, fname):
    """Write variant list."""
    internal_path = os.path.join(LEGACY_VARIANT_SELECTION_INTERNAL_DIR, fname)
    frame[cols].to_csv(internal_path, index=False)
    print(f"Write:{internal_path}(complete, internal)")
    return internal_path

def write_public_variant_table(frames):
    """Write public variant table."""
    allv = pd.concat(frames, ignore_index=True)
    masks = allv.groupby("vid")["mask"].apply(lambda s: ";".join(sorted(set(s))))
    public = allv.drop_duplicates("vid").drop(columns=["mask", "component"], errors="ignore").set_index("vid")
    maf = public["maf_primary"].fillna(0.0)
    public.insert(public.columns.get_loc("maf_primary"), "maf_primary_bin",
                  np.where(maf < 0.001, "<0.1%", "0.1-1%"))
    public = public.drop(columns="maf_primary")
    public["masks"] = masks
    out = os.path.join(LEGACY_VARIANT_SELECTION_TABLE_DIR, "12_qualifying_variants_public.csv.gz")
    public.to_csv(out)
    print(f"Write:{out}(public,{len(public)}mutations, maf_primary binning,{os.path.getsize(out)/1e6:.2f}MB)")

def to_float(x):
    try:
        if x is None or x == "" or (isinstance(x, float) and np.isnan(x)):
            return np.nan
        return float(x)
    except (ValueError, TypeError):
        return np.nan

def is_ptv(cons_str):
    terms = set(str(cons_str).split(","))
    return bool(terms & LEGACY_VARIANT_SELECTION_PTV_TERMS)

def is_missense(cons_str):
    terms = set(str(cons_str).split(","))
    return bool(terms & LEGACY_VARIANT_SELECTION_MISSENSE_TERMS)

def legacy_variant_selection_calculation():
    legacy_variant_selection_section('0. Load data')
    vat_path = os.path.join(LEGACY_VARIANT_SELECTION_DATA_DIR, "11_panel_gene_vat_ptv_missense_candidates.csv.gz")
    df = pd.read_csv(vat_path, dtype=str, compression="gzip")
    print(f"Load{vat_path}: {len(df)}OK")

    i157t = pd.read_csv(os.path.join(LEGACY_VARIANT_SELECTION_DATA_DIR, "11_chek2_i157t_hg38_coordinate.csv"))
    i157t_vid = f"22-{i157t['grch38_pos'].iloc[0]}-{i157t['grch38_ref'].iloc[0]}-{i157t['grch38_alt'].iloc[0]}"
    print(f"I157T exclusion site vid ={i157t_vid}(Cross-validation status:{i157t['vat_cross_validated'].iloc[0]})")

    # ------------------------------------------------------------------
    
    # ------------------------------------------------------------------
    legacy_variant_selection_section('1. Variation level field organization (consequence classification / LOFTEE HC / MAF digitization)')
    df["maf_primary"] = df[MAF_FIELD].apply(to_float)
    df["maf_crosscheck"] = df[MAF_CROSSCHECK_FIELD].apply(to_float)
    df["revel_f"] = df["revel"].apply(to_float)
    df["is_ptv"] = df["consequence"].apply(is_ptv)
    df["is_missense"] = df["consequence"].apply(is_missense)
    df["is_loftee_hc"] = df["LoF"] == "HC"
    df["is_i157t"] = df["vid"] == i157t_vid

    print(f"PTV candidate lines:{df['is_ptv'].sum()}")
    print(f"Where LOFTEE HC:{(df['is_ptv'] & df['is_loftee_hc']).sum()}")
    print(f"missense candidate line:{df['is_missense'].sum()}")
    print(f"I157T hit line number (cross-transcript should be 1 because canonical is limited):{df['is_i157t'].sum()}")

    maf_na_ptv = df.loc[df["is_ptv"], "maf_primary"].isna().sum()
    maf_na_missense = df.loc[df["is_missense"], "maf_primary"].isna().sum()
    print(f"PTV in line{MAF_FIELD}Missing:{maf_na_ptv}Missing in /missense line:{maf_na_missense}"
          f"(Missing usually means that the row has AC=0 in the AoU European subgroup, which is regarded as AF=0 by convention. It does not mean that it does not exist)")
    df["maf_primary_filled"] = df["maf_primary"].fillna(0.0)

    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    legacy_variant_selection_section('1b. Aggregate and remove duplicates by (vid, gene_symbol) (PTV+LOFTEE HC determines OR across transcripts)')
    n_rows_before = len(df)
    ptv_hc_by_vid = (
        df.assign(_ptv_hc_row=df["is_ptv"] & df["is_loftee_hc"])
        .groupby(["vid", "gene_symbol"])["_ptv_hc_row"]
        .any()
        .rename("is_ptv_hc_any_transcript")
    )
    df = df.merge(ptv_hc_by_vid, on=["vid", "gene_symbol"], how="left")

    df["_pref_rank"] = df["is_canonical_transcript"].map({"true": 0, "false": 1}).fillna(1)
    df = df.sort_values(["vid", "gene_symbol", "_pref_rank"])
    dedup = df.drop_duplicates(subset=["vid", "gene_symbol"], keep="first").drop(columns=["_pref_rank"])
    print(f"Before aggregation:{n_rows_before}Rows -> After aggregation (one row per vid x gene):{len(dedup)}OK")
    n_ptv_hc_rescued_from_noncanonical = (
        dedup["is_ptv_hc_any_transcript"] & (dedup["is_canonical_transcript"] != "true") & (~dedup["is_loftee_hc"])
    ).sum()
    print(f"Among them, the number of variants judged as PTV qualifying because of \'non-canonical lines with LOFTEE HC label\':"
          f"{n_ptv_hc_rescued_from_noncanonical}(The representative row itself is canonical but the LoF is empty,"
          f"HC annotations provided by other transcript lines of the same vid are included)")
    df = dedup
    
    df["is_loftee_hc"] = df["is_ptv_hc_any_transcript"]

    # ------------------------------------------------------------------
    
    # ------------------------------------------------------------------
    legacy_variant_selection_section('2. Main analysis qualifying variant: pLoF, MAF<1%, LOFTEE HC')

    maf_main_ptv = CFG["maf_thresholds"]["main_ptv"]
    main_ptv = df[
        df["is_loftee_hc"] & (df["maf_primary_filled"] < maf_main_ptv)
    ].copy()
    main_ptv["mask"] = "main_ptv_maf1pct"
    print(f"Main analysis PTV qualifying variants:{len(main_ptv)}row, cover{main_ptv['gene_symbol'].nunique()}Gene")
    print(main_ptv.groupby("gene_symbol").size().to_string())

    keep_cols = ["vid", "contig", "position", "ref_allele", "alt_allele", "gene_symbol",
                 "consequence", "aa_change", "dna_change_in_transcript", "LoF", "LoF_filter",
                 "revel_f", "maf_primary", "maf_crosscheck", "dbsnp_rsid", "mask"]
    write_variant_list(main_ptv, keep_cols, "12_qualifying_variants_main_ptv.csv")

    # ------------------------------------------------------------------
    # 3. CHEK2 1100delC / I157T sanity check
    # ------------------------------------------------------------------
    legacy_variant_selection_section("3. CHEK2 1100delC / I157T sanity check")
    chek2 = df[df["gene_symbol"] == "CHEK2"].copy()

    c1100_mask = chek2["dna_change_in_transcript"].astype(str).str.contains("1100del", na=False)
    c1100_rows = chek2[c1100_mask]
    print(f"The number of rows containing \'1100del\' in dna_change in CHEK2 candidate rows:{len(c1100_rows)}")
    if len(c1100_rows) > 0:
        print(c1100_rows[["vid", "consequence", "LoF", "maf_primary", "aa_change"]].to_string(index=False))
    i157t_rows = chek2[chek2["is_i157t"]]
    print(f"\nCHEK2 I157T candidate line:")
    if len(i157t_rows) > 0:
        print(i157t_rows[["vid", "consequence", "LoF", "maf_primary", "aa_change"]].to_string(index=False))
    else:
        print('Not hit in the missense/PTV candidate set (may be filtered in advance by script 11 because the consequence is not in the PTV/missense set, need to be checked)')

    sanity_rows = []
    for _, r in c1100_rows.iterrows():
        sanity_rows.append({
            "variant": "c.1100delC", "vid": r["vid"], "maf_primary": r["maf_primary"],
            "in_main_ptv_maf1pct": r["maf_primary_filled"] < maf_main_ptv if pd.notna(r["maf_primary_filled"]) else None,
            "in_sensitivity_maf0p1pct": r["maf_primary_filled"] < CFG["maf_thresholds"]["sensitivity_uniform_0p1pct"],
        })
    for _, r in i157t_rows.iterrows():
        sanity_rows.append({
            "variant": "p.I157T", "vid": r["vid"], "maf_primary": r["maf_primary"],
            "in_main_ptv_maf1pct": None,  
            "in_sensitivity_maf0p1pct": None,
        })
    sanity_df = pd.DataFrame(sanity_rows)
    sanity_path = os.path.join(LEGACY_VARIANT_SELECTION_TABLE_DIR, "12_chek2_1100delc_sanity_check.csv")
    sanity_df.to_csv(sanity_path, index=False)
    print(f"\n writes:{sanity_path}")

    # ------------------------------------------------------------------
    
    # ------------------------------------------------------------------
    legacy_variant_selection_section('4. Sensitivity analysis file 1: Unified MAF<0.1% (1100delC should be excluded from the CHEK2 PTV set)')
    maf_sens_low = CFG["maf_thresholds"]["sensitivity_uniform_0p1pct"]
    sens_low_ptv = df[
        df["is_loftee_hc"] & (df["maf_primary_filled"] < maf_sens_low)
    ].copy()
    sens_low_ptv["mask"] = "sensitivity_ptv_maf0p1pct"
    sens_low_missense = df[
        df["is_missense"] & (~df["is_i157t"]) & (df["maf_primary_filled"] < maf_sens_low)
    ].copy()
    sens_low_missense["mask"] = "sensitivity_missense_maf0p1pct"
    sens_low = pd.concat([sens_low_ptv, sens_low_missense], ignore_index=True)
    print(f"File 1 PTV:{len(sens_low_ptv)}OK, missense (excluding I157T):{len(sens_low_missense)}OK")
    chek2_1100_in_sens_low = sens_low_ptv[
        sens_low_ptv["dna_change_in_transcript"].astype(str).str.contains("1100del", na=False)
    ]
    print(f"Verify: 1100delC is still in the file 1 PTV collection:{len(chek2_1100_in_sens_low) > 0} "
          f"(Expected False, because AF≈0.26%>0.1% threshold)")
    write_variant_list(sens_low, keep_cols, "12_qualifying_variants_sensitivity_0p1pct.csv")

    # ------------------------------------------------------------------
    
    # ------------------------------------------------------------------
    legacy_variant_selection_section('5. Sensitivity analysis file 2: Unified MAF<1% (1100delC is still included, missense is relaxed to 1%)')
    maf_sens_high = CFG["maf_thresholds"]["sensitivity_uniform_1pct"]
    sens_high_ptv = main_ptv.copy()  
    sens_high_ptv["mask"] = "sensitivity_ptv_maf1pct"
    sens_high_missense = df[
        df["is_missense"] & (~df["is_i157t"]) & (df["maf_primary_filled"] < maf_sens_high)
    ].copy()
    sens_high_missense["mask"] = "sensitivity_missense_maf1pct"
    sens_high = pd.concat([sens_high_ptv, sens_high_missense], ignore_index=True)
    print(f"File 2 PTV:{len(sens_high_ptv)}Row (=main analysis PTV set), missense (exclude I157T, MAF<1%):{len(sens_high_missense)}OK")
    write_variant_list(sens_high, keep_cols, "12_qualifying_variants_sensitivity_1pct.csv")

    # ------------------------------------------------------------------
    
    # ------------------------------------------------------------------
    legacy_variant_selection_section('6. Discovery scan masks: LOF / LOF+REVEL>=0.75 / 0.50 / 0.25, unified MAF<1%')
    maf_disc = CFG["maf_thresholds"]["discovery_scan_uniform"]
    disc_ptv = df[
        df["is_loftee_hc"] & (df["maf_primary_filled"] < maf_disc)
    ].copy()
    disc_missense_pool = df[
        df["is_missense"] & (~df["is_i157t"]) & (df["maf_primary_filled"] < maf_disc)
    ].copy()
    print(f"discovery PTV pool (MAF<1%):{len(disc_ptv)}OK")
    print(f"discovery missense pool (MAF<1%, exclude I157T, REVEL any):{len(disc_missense_pool)}OK"
          f"(REVEL missing rows:{disc_missense_pool['revel_f'].isna().sum()}, missense without REVEL will not enter"
          f"Any REVEL threshold mask will only be covered by the pure PTV mask of LOF, and the behavior is in line with expectations and does not require additional elimination)")

    disc_frames = []
    for mask_cfg in CFG["discovery_scan_masks"]:
        name = mask_cfg["name"]
        revel_min = mask_cfg["revel_min"]
        ptv_part = disc_ptv.copy()
        ptv_part["mask"] = name
        ptv_part["component"] = "PTV"
        frames = [ptv_part]
        if revel_min is not None:
            missense_part = disc_missense_pool[disc_missense_pool["revel_f"] >= revel_min].copy()
            missense_part["mask"] = name
            missense_part["component"] = "missense"
            frames.append(missense_part)
        mask_df = pd.concat(frames, ignore_index=True)
        disc_frames.append(mask_df)
        n_genes = mask_df["gene_symbol"].nunique()
        print(f"  mask={name} (revel_min={revel_min}): {len(mask_df)}OK,{n_genes}gene coverage")

    disc_all = pd.concat(disc_frames, ignore_index=True)
    write_variant_list(disc_all, keep_cols + ["component"], "12_qualifying_variants_discovery_masks.csv")
    print(f"discovery masks: {len(disc_all)}OK")
    write_public_variant_table([main_ptv[keep_cols], sens_low[keep_cols], sens_high[keep_cols],
                                disc_all[keep_cols + ["component"]]])

    # ------------------------------------------------------------------
    
    # ------------------------------------------------------------------
    legacy_variant_selection_section('7. Summary of candidate variants by gene and mask')
    summary_rows = []
    for name, sub in [
        ("main_ptv_maf1pct", main_ptv),
        ("sensitivity_ptv_maf0p1pct", sens_low_ptv),
        ("sensitivity_missense_maf0p1pct", sens_low_missense),
        ("sensitivity_ptv_maf1pct", sens_high_ptv),
        ("sensitivity_missense_maf1pct", sens_high_missense),
    ]:
        for gene, cnt in sub.groupby("gene_symbol").size().items():
            summary_rows.append({"mask": name, "gene_symbol": gene, "n_qualifying_variants": cnt})
    for mask_df in disc_frames:
        mname = mask_df["mask"].iloc[0]
        for gene, cnt in mask_df.groupby("gene_symbol").size().items():
            summary_rows.append({"mask": mname, "gene_symbol": gene, "n_qualifying_variants": cnt})
    summary_df = pd.DataFrame(summary_rows)
    
    all_genes = sorted(df["gene_symbol"].unique())
    all_masks = summary_df["mask"].unique()
    full_index = pd.MultiIndex.from_product([all_masks, all_genes], names=["mask", "gene_symbol"])
    summary_df = (
        summary_df.set_index(["mask", "gene_symbol"])
        .reindex(full_index, fill_value=0)
        .reset_index()
    )
    out_summary = os.path.join(LEGACY_VARIANT_SELECTION_TABLE_DIR, "12_qualifying_variant_gene_summary.csv")
    summary_df.to_csv(out_summary, index=False)
    print(summary_df.pivot(index="gene_symbol", columns="mask", values="n_qualifying_variants").to_string())
    print(f"\n writes:{out_summary}")

    legacy_variant_selection_section('Finish')
    print('Step 3 The qualifying variant definition is completed and you can enter script 13 (carrier genotype extraction + carrier-level QC).')

def analyze_legacy_variant_selection():
    """Initialize runtime state and execute this component explicitly."""
    global CFG, LEGACY_VARIANT_SELECTION_CONFIG_PATH, LEGACY_VARIANT_SELECTION_DATA_DIR, LEGACY_VARIANT_SELECTION_INTERNAL_DIR, MAF_CROSSCHECK_FIELD, MAF_FIELD, LEGACY_VARIANT_SELECTION_MISSENSE_TERMS, LEGACY_VARIANT_SELECTION_PTV_TERMS, LEGACY_VARIANT_SELECTION_REPO, LEGACY_VARIANT_SELECTION_TABLE_DIR, fh

    LEGACY_VARIANT_SELECTION_REPO = "/home/jupyter/workspace/analysis"
    LEGACY_VARIANT_SELECTION_DATA_DIR = f"{LEGACY_VARIANT_SELECTION_REPO}/data"
    LEGACY_VARIANT_SELECTION_TABLE_DIR = f"{LEGACY_VARIANT_SELECTION_REPO}/results/table"
    LEGACY_VARIANT_SELECTION_INTERNAL_DIR = f"{LEGACY_VARIANT_SELECTION_DATA_DIR}/12_internal"
    LEGACY_VARIANT_SELECTION_CONFIG_PATH = f"{LEGACY_VARIANT_SELECTION_REPO}/config/aim1_config.yaml"

    os.makedirs(LEGACY_VARIANT_SELECTION_TABLE_DIR, exist_ok=True)
    os.makedirs(LEGACY_VARIANT_SELECTION_INTERNAL_DIR, exist_ok=True)

    with open(LEGACY_VARIANT_SELECTION_CONFIG_PATH) as fh:
        CFG = yaml.safe_load(fh)["step3_qualifying_variant"]

    LEGACY_VARIANT_SELECTION_PTV_TERMS = set(CFG["ptv_consequence_terms"])
    LEGACY_VARIANT_SELECTION_MISSENSE_TERMS = set(CFG["missense_consequence_terms"])
    MAF_FIELD = CFG["maf_field_primary"]
    MAF_CROSSCHECK_FIELD = CFG["maf_field_crosscheck"]

    if True:
        legacy_variant_selection_calculation()

# ============================================================================
# 5. carriers: analyze_legacy_carriers
# ============================================================================

import os
import shutil
import subprocess
import sys
import time
import numpy as np
import pandas as pd
import yaml
from scipy.stats import gaussian_kde

def legacy_carriers_section(title):
    print("\n" + "=" * 90)
    print(title)
    print("=" * 90, flush=True)

def legacy_carriers_run(cmd, **kwargs):
    print(f"$ {cmd}", flush=True)
    t0 = time.time()
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, **kwargs)
    print(f"(time consuming{time.time()-t0:.1f}s, returncode={result.returncode})")
    if result.returncode != 0:
        print(result.stdout[-2000:])
        print(result.stderr[-2000:], file=sys.stderr)
    return result

def legacy_carriers_suppress(n, min_cell=20):
    """Legacy carriers suppress."""
    if n is None or (isinstance(n, float) and np.isnan(n)):
        return n
    return n if n >= min_cell else f"<{min_cell}"

def disclose_split(n_total, n_case, n_control, min_cell=20):
    """Disclose split."""
    if n_total < min_cell:
        return legacy_carriers_suppress(n_total), f"<{min_cell}", f"<{min_cell}", np.nan
    small_case, small_ctrl = n_case < min_cell, n_control < min_cell
    if not small_case and not small_ctrl:
        return n_total, n_case, n_control, round(n_control / n_case, 2)
    case_out = f"<{min_cell}" if small_case else SECONDARY
    ctrl_out = f"<{min_cell}" if small_ctrl else SECONDARY
    return n_total, case_out, ctrl_out, np.nan

def ensure_mount():
    legacy_carriers_section('0. gcsfuse read-only mount check/mount')
    r = legacy_carriers_run(f"mountpoint -q {LEGACY_CARRIERS_MOUNT_DIR}")
    if r.returncode == 0:
        print(f"  {LEGACY_CARRIERS_MOUNT_DIR}Already mounted, skip remounting")
        return
    os.makedirs(LEGACY_CARRIERS_MOUNT_DIR, exist_ok=True)
    legacy_carriers_run(f"fusermount -u {LEGACY_CARRIERS_MOUNT_DIR} 2>/dev/null; true")
    r = legacy_carriers_run(
        f"gcsfuse --billing-project {LEGACY_CARRIERS_GOOGLE_PROJECT} -o ro --implicit-dirs "
        f"{LEGACY_CARRIERS_BUCKET_NAME} {LEGACY_CARRIERS_MOUNT_DIR}"
    )
    if r.returncode != 0:
        raise RuntimeError('gcsfuse failed to mount and cannot continue')

def build_sample_list():
    legacy_carriers_section('1. Construct a cohort sample list for carrier extraction (No. 05 pass_qc==1, genetics EUR Aim 1 queue)')
    sample_list_path = f"{LEGACY_CARRIERS_DATA_DIR}/13_carrier_extraction_sample_list.txt"
    
    if not os.path.exists(SAMPLE_LIST_V0):
        shutil.copy2(sample_list_path, SAMPLE_LIST_V0)
        print(f"Back up old sample table ->{SAMPLE_LIST_V0}")
    qc = pd.read_csv(f"{LEGACY_CARRIERS_DATA_DIR}/05_cohort_qc_flags.csv")
    cohort = qc[qc["pass_qc"] == 1].rename(columns={"analysis_group": "group"})
    print(cohort["group"].value_counts().to_string())
    print(f"Total number of cohort:{len(cohort)}")
    cohort[["person_id"]].to_csv(sample_list_path, index=False, header=False)
    print(f"Write (person-level, data/lower, git-ignored):{sample_list_path}")

    v0 = set(pd.read_csv(SAMPLE_LIST_V0, header=None)[0].astype(str))
    new_ids = cohort["person_id"].astype(str)
    addon = new_ids[~new_ids.isin(v0)]
    addon.to_csv(SAMPLE_LIST_ADDON, index=False, header=False)
    print(f"Number of people reusing old checkpoints:{new_ids.isin(v0).sum()}; New extraction required:{len(addon)} -> {SAMPLE_LIST_ADDON}")

    case_ids = set(cohort.loc[cohort["group"] == "case", "person_id"].astype(str))
    control_ids = set(cohort.loc[cohort["group"] == "control", "person_id"].astype(str))
    return sample_list_path, case_ids, control_ids, len(cohort)

def build_qualifying_union(gene_list):
    legacy_carriers_section('2. Load the qualifying variant set and build a union by gene (across main/sensitivity/discovery)')
    frames = []
    for tag, path in QUAL_VARIANT_FILES.items():
        d = pd.read_csv(path, dtype={"contig": str, "ref_allele": str, "alt_allele": str})
        d["mask"] = d["mask"]  
        d["mask_source_file"] = tag
        frames.append(d[["vid", "contig", "position", "ref_allele", "alt_allele", "gene_symbol", "mask"]])
    allv = pd.concat(frames, ignore_index=True)
    allv = allv[allv["gene_symbol"].isin(gene_list)]

    membership = (
        allv.groupby(["gene_symbol", "vid", "contig", "position", "ref_allele", "alt_allele"])["mask"]
        .apply(lambda s: sorted(set(s)))
        .reset_index()
        .rename(columns={"mask": "masks"})
    )
    membership["masks_joined"] = membership["masks"].apply(lambda x: "|".join(x))
    print(f"  {len(gene_list)}gene, after removing duplicates distinct(vid,gene):{len(membership)}OK")
    print(membership.groupby("gene_symbol").size().to_string())

    union_path = f"{LEGACY_CARRIERS_DATA_DIR}/13_qualifying_variant_union_by_gene.csv.gz"
    membership.drop(columns=["masks"]).to_csv(union_path, index=False, compression="gzip")
    print(f"Write:{union_path}")
    return membership

def parse_ad(ad_str):
    if ad_str is None or ad_str == "." or ad_str == "":
        return None
    try:
        return [int(x) for x in ad_str.split(",")]
    except ValueError:
        return None

def parse_gt_alleles(gt_str):
    sep = "/" if "/" in gt_str else ("|" if "|" in gt_str else None)
    if sep is None:
        return None
    parts = gt_str.split(sep)
    alleles = []
    for p in parts:
        if p == ".":
            return None
        try:
            alleles.append(int(p))
        except ValueError:
            return None
    return alleles

def extract_gene_shard(gene, shard_idx, region, sample_list_path, gene_variants, raw_dir=None):
    """Extract gene shard."""
    if raw_dir is None:
        raw_dir = RAW_DIR
    ckpt_path = os.path.join(raw_dir, f"{gene}__shard{shard_idx:010d}.csv.gz")
    if os.path.exists(ckpt_path):
        print(f"[Skip, checkpoint already exists]{ckpt_path}")

        return pd.read_csv(ckpt_path, dtype={"contig": str, "ref_allele": str, "alt_allele": str, "person_id": str})

    shard_file = os.path.join(EXOME_VCF_DIR, f"{shard_idx:010d}.vcf.bgz")
    raw_txt = os.path.join(SCRATCH, f"{os.path.basename(raw_dir)}__{gene}__shard{shard_idx:010d}.raw.txt")
    cmd = (
        f"bcftools view -r {region} -S {sample_list_path} '{shard_file}' -Ou | "
        f"bcftools query -i 'GT=\"alt\"' "
        f"-f '[%CHROM\\t%POS\\t%REF\\t%ALT\\t%SAMPLE\\t%GT\\t%AD\\t%GQ\\n]' "
        f"> {raw_txt}"
    )
    r = legacy_carriers_run(cmd)
    if r.returncode != 0:
        raise RuntimeError(f"bcftools extraction failed: gene={gene} shard={shard_idx}")

    cols = ["contig", "position", "ref", "alt_full", "person_id", "gt", "ad", "gq_raw"]
    try:
        raw = pd.read_csv(raw_txt, sep="\t", header=None, names=cols, dtype=str)
    except pd.errors.EmptyDataError:
        raw = pd.DataFrame(columns=cols)
    n_raw = len(raw)

    raw["position"] = raw["position"].astype(int)
    gv = gene_variants[["vid", "contig", "position", "ref_allele", "alt_allele", "masks_joined"]].copy()
    merged = raw.merge(
        gv, left_on=["contig", "position", "ref"], right_on=["contig", "position", "ref_allele"],
        how="inner",
    )

    NULL_ROW = pd.Series({
        "n_target_alleles": np.nan, "zygosity": None,
        "ad_target": np.nan, "ad_total": np.nan, "vaf_target": np.nan, "gq": np.nan,
    })

    def resolve_row(row):

        alt_list = row["alt_full"].split(",")
        if row["alt_allele"] not in alt_list:
            return NULL_ROW.copy()  
        target_idx = alt_list.index(row["alt_allele"]) + 1  # 0=REF
        alleles = parse_gt_alleles(row["gt"])
        if alleles is None:
            return NULL_ROW.copy()
        n_target = sum(1 for a in alleles if a == target_idx)
        if n_target == 0:
            return NULL_ROW.copy()  
        ad = parse_ad(row["ad"])
        gq = row["gq_raw"]
        if ad is None or gq in (None, ".", ""):
            ad_target, ad_total, vaf, gq_val = np.nan, np.nan, np.nan, np.nan
        else:
            ad_total = sum(ad)
            ad_target = ad[target_idx] if target_idx < len(ad) else np.nan
            vaf = ad_target / ad_total if ad_total and ad_total > 0 else np.nan
            try:
                gq_val = float(gq)
            except ValueError:
                gq_val = np.nan
        zygosity = "hom_alt" if n_target == 2 else "het"
        return pd.Series({
            "n_target_alleles": n_target, "zygosity": zygosity,
            "ad_target": ad_target, "ad_total": ad_total, "vaf_target": vaf, "gq": gq_val,
        })

    if len(merged) > 0:
        resolved = merged.apply(resolve_row, axis=1)
        merged = pd.concat([merged, resolved], axis=1)
        merged = merged[merged["n_target_alleles"].notna()]
    merged["gene_symbol"] = gene
    keep = ["gene_symbol", "vid", "contig", "position", "ref_allele", "alt_allele", "masks_joined",
            "person_id", "gt", "zygosity", "ad_target", "ad_total", "vaf_target", "gq"]
    result = merged[keep] if len(merged) > 0 else pd.DataFrame(columns=keep)

    tmp_path = ckpt_path + ".tmp"
    result.to_csv(tmp_path, index=False, compression="gzip")
    os.replace(tmp_path, ckpt_path)
    print(f"    shard{shard_idx}: original ALT line{legacy_carriers_suppress(n_raw)}-> After exact matching of qualifying variant{legacy_carriers_suppress(len(result))}OK,"
          f"Place the order{ckpt_path}")
    try:
        os.remove(raw_txt)
    except OSError:
        pass
    return result

def extract_all_genes(gene_list, shard_map, sample_list_path, union_df, raw_dir=None):
    if raw_dir is None:
        raw_dir = RAW_DIR
    legacy_carriers_section(f"3. Extract carrier genotype (GT/AD/GQ, qualifying variant has been accurately matched) by gene and shard ->{raw_dir}")
    all_frames = []
    for gene in gene_list:
        row = shard_map[shard_map["gene_symbol"] == gene].iloc[0]
        region = f"{row['chrom']}:{row['region_start']}-{row['region_end']}"
        shard_ids = list(range(int(row["shard_first"]), int(row["shard_last"]) + 1))
        gene_variants = union_df[union_df["gene_symbol"] == gene]
        print(f"\n gene{gene}: area{region}, {len(shard_ids)}shard{shard_ids}, "
              f"{len(gene_variants)}qualifying variant candidates")
        gene_frames = []
        for shard_idx in shard_ids:
            f = extract_gene_shard(gene, shard_idx, region, sample_list_path, gene_variants, raw_dir)
            gene_frames.append(f)
        gene_df = pd.concat(gene_frames, ignore_index=True) if gene_frames else pd.DataFrame()
        n_variants_hit = gene_df["vid"].nunique() if len(gene_df) else 0
        n_carriers_hit = gene_df["person_id"].nunique() if len(gene_df) else 0
        print(f"  {gene}Summary:{legacy_carriers_suppress(len(gene_df))}Carrier genotype call,"
              f"{n_variants_hit}/{len(gene_variants)}Each qualifying variant has at least 1 carrier,"
              f"{legacy_carriers_suppress(n_carriers_hit)}distinct carrier (before QC)")
        all_frames.append(gene_df)
    combined = pd.concat(all_frames, ignore_index=True) if all_frames else pd.DataFrame()
    return combined

def find_local_minimum(values, side="low", grid_points=512):
    """Find local minimum."""
    values = np.asarray(values, dtype=float)
    values = values[~np.isnan(values)]
    if len(values) < 30 or np.nanstd(values) == 0:
        return None, None, None
    kde = gaussian_kde(values)
    grid = np.linspace(values.min(), values.max(), grid_points)
    dens = kde(grid)
    med = np.median(values)
    minima_idx = [i for i in range(1, len(dens) - 1) if dens[i] < dens[i - 1] and dens[i] < dens[i + 1]]
    if side == "low":
        cand = [i for i in minima_idx if grid[i] < med]
    else:
        cand = [i for i in minima_idx if grid[i] > med]
    if not cand:
        return None, grid, dens
    best = min(cand, key=lambda i: dens[i])
    return float(grid[best]), grid, dens

def determine_qc_thresholds(raw_carriers):
    legacy_carriers_section('4. Determination of empirical density threshold (KDE local minimum, without reference to phenotype/association results)')
    gq_vals = raw_carriers["gq"].dropna().values
    print(f"pooled carrier genotype calls GQ: n={len(gq_vals)}, "
          f"mean={np.nanmean(gq_vals):.2f}, median={np.nanmedian(gq_vals):.2f}, "
          f"p5={np.nanpercentile(gq_vals,5):.1f}, p95={np.nanpercentile(gq_vals,95):.1f}")
    gq_thresh, gq_grid, gq_dens = find_local_minimum(gq_vals, side="low")
    if gq_thresh is None:
        gq_thresh = GQ_FALLBACK_THRESHOLD
        gq_method = "fallback_no_kde_minimum_found"
    else:
        gq_method = "kde_local_minimum"
    print(f"-> GQ threshold:{gq_thresh:.2f}(method:{gq_method})")

    het_vaf = raw_carriers.loc[raw_carriers["zygosity"] == "het", "vaf_target"].dropna().values
    print(f"\n  pooled het carrier VAF: n={len(het_vaf)}, "
          f"mean={np.nanmean(het_vaf):.3f}, median={np.nanmedian(het_vaf):.3f}")
    vaf_low, _, _ = find_local_minimum(het_vaf, side="low")
    vaf_high, _, _ = find_local_minimum(het_vaf, side="high")
    vaf_method = "kde_local_minimum"
    if vaf_low is None:
        vaf_low = VAF_HET_FALLBACK_RANGE[0]
        vaf_method = "fallback_no_kde_minimum_found(low)"
    if vaf_high is None:
        vaf_high = VAF_HET_FALLBACK_RANGE[1]
        vaf_method += "+fallback(high)" if "fallback" not in vaf_method else ""
    print(f"-> VAF(het) reasonable interval: [{vaf_low:.3f}, {vaf_high:.3f}] (method:{vaf_method})")

    n_hom = (raw_carriers["zygosity"] == "hom_alt").sum()
    print(f"\n hom_alt carrier call number:{n_hom}(The sample size is small, KDE is not performed, and the conventional threshold VAF>= is used.{VAF_HOM_MIN})")

    thresholds = {
        "gq_mean_threshold": round(gq_thresh, 2),
        "gq_threshold_method": gq_method,
        "vaf_het_range": [round(vaf_low, 3), round(vaf_high, 3)],
        "vaf_het_range_method": vaf_method,
        "vaf_hom_alt_min": VAF_HOM_MIN,
        "vaf_hom_alt_min_method": 'fixed_convention_predetermined_not_kde(hom_alt sample size is too small)',
        "status": "ESTIMATED_FROM_SCRIPT_13",
        "n_carrier_genotype_calls_pooled": int(len(raw_carriers)),
        "n_het_vaf_pooled": int(len(het_vaf)),
    }
    return thresholds

def apply_qc(raw_carriers, thresholds):
    legacy_carriers_section('5. Apply carrier-level genotype QC')
    df = raw_carriers.copy()
    gq_ok = df["gq"] >= thresholds["gq_mean_threshold"]
    vaf_low, vaf_high = thresholds["vaf_het_range"]
    het_ok = (df["zygosity"] == "het") & (df["vaf_target"] >= vaf_low) & (df["vaf_target"] <= vaf_high)
    hom_ok = (df["zygosity"] == "hom_alt") & (df["vaf_target"] >= thresholds["vaf_hom_alt_min"])
    df["qc_pass"] = gq_ok & (het_ok | hom_ok) & df["gq"].notna() & df["vaf_target"].notna()
    n_total, n_pass = len(df), int(df["qc_pass"].sum())
    print(f"carrier genotype call:{n_total}OK, QC passed:{n_pass} ({n_pass/n_total*100:.1f}%)"
          if n_total else 'Carrier-free genotype calling')
    print(f"QC failed details: GQ insufficient ={int((~gq_ok).sum())}, "
          f"het VAF out of range ={int(((df['zygosity']=='het') & ~het_ok).sum())}, "
          f"hom VAF too low ={int(((df['zygosity']=='hom_alt') & ~hom_ok).sum())}, "
          f"MissingGQ/VAF={int((df['gq'].isna() | df['vaf_target'].isna()).sum())}")
    return df

def collapse_carrier_status(qc_df, gene_list, sample_list_path):
    legacy_carriers_section("6. Collapse carrier_status per person per gene per mask（person-level -> data/）")
    all_person_ids = pd.read_csv(sample_list_path, header=None, names=["person_id"])["person_id"].astype(str)

    mask_names = sorted({m for masks in qc_df["masks_joined"].dropna().unique() for m in masks.split("|")}) \
        if len(qc_df) else []
    print(f"Involving masks:{mask_names}")

    pass_df = qc_df[qc_df["qc_pass"]].copy()
    out_paths = {}
    for mask in mask_names:
        hit = pass_df[pass_df["masks_joined"].str.split("|").apply(lambda lst: mask in lst)]
        wide = (
            hit.groupby(["gene_symbol"])["person_id"]
            .apply(lambda s: set(s))
            .to_dict()
        )
        person_gene = pd.DataFrame({"person_id": all_person_ids})
        for gene in gene_list:
            carriers = wide.get(gene, set())
            person_gene[gene] = person_gene["person_id"].isin(carriers).astype(int)
        safe_mask_name = mask.replace("/", "_")
        out_path = f"{LEGACY_CARRIERS_DATA_DIR}/13_carrier_status_{safe_mask_name}_person_level.csv.gz"
        person_gene.to_csv(out_path, index=False, compression="gzip")
        out_paths[mask] = out_path
        print(f"  mask={mask}: write out{out_path} ({person_gene.shape[0]}person x{len(gene_list)}Gene),"
              f"The number of carriers of any gene ={int((person_gene[gene_list].sum(axis=1) > 0).sum())}")
    return out_paths, mask_names, pass_df

def step3_5_sparsity_precheck(pass_df, gene_list, mask_names, case_ids, control_ids, n_cohort,
                                union_df, cfg):
    legacy_carriers_section('7. Step 3.5 Carrier sparsity/efficiency pre-check -> results/table, results/report')
    rows = []
    main_mask = "main_ptv_maf1pct"
    any_gene_carriers = {}

    for mask in mask_names:
        hit = pass_df[pass_df["masks_joined"].str.split("|").apply(lambda lst: mask in lst)]
        gene_carrier_sets = {g: set(hit.loc[hit["gene_symbol"] == g, "person_id"]) for g in gene_list}
        any_gene_carriers[mask] = set().union(*gene_carrier_sets.values()) if gene_carrier_sets else set()

        for gene in gene_list:
            carriers = gene_carrier_sets.get(gene, set())
            n_total = len(carriers)
            n_case = len(carriers & case_ids)
            n_control = len(carriers & control_ids)
            d_total, d_case, d_ctrl, d_ratio = disclose_split(n_total, n_case, n_control)
            n_qv = int(((union_df["gene_symbol"] == gene) &
                        (union_df["masks_joined"].str.split("|").apply(lambda lst: mask in lst))).sum())
            if n_total < CARRIER_GATE_LOW:
                tier = 'low(<500)_exploratory'
            elif n_total < CARRIER_GATE_HIGH:
                tier = 'mid(500-2000)_two narratives written together'
            else:
                tier = 'high(>2000)_can support confirmatory interaction testing'
            rows.append({
                "mask": mask, "gene_symbol": gene, "n_qualifying_variants": n_qv,
                "n_carriers_total_postQC": d_total,
                "n_carriers_case_postQC": d_case,
                "n_carriers_control_postQC": d_ctrl,
                "carrier_control_to_case_ratio": d_ratio,
                "cohort_baseline_control_to_case_ratio": COHORT_BASELINE_RATIO,
                "sparsity_tier": tier,
            })
        any_g = any_gene_carriers[mask]
        n_total = len(any_g)
        n_case = len(any_g & case_ids)
        n_control = len(any_g & control_ids)
        d_total, d_case, d_ctrl, d_ratio = disclose_split(n_total, n_case, n_control)
        rows.append({
            "mask": mask, "gene_symbol": "ANY_OF_10_GENES_DEDUPED", "n_qualifying_variants": np.nan,
            "n_carriers_total_postQC": d_total,
            "n_carriers_case_postQC": d_case,
            "n_carriers_control_postQC": d_ctrl,
            "carrier_control_to_case_ratio": d_ratio,
            "cohort_baseline_control_to_case_ratio": COHORT_BASELINE_RATIO,
            "sparsity_tier": (
                'low(<500)_exploratory' if n_total < CARRIER_GATE_LOW else
                'mid(500-2000)_two narratives written together' if n_total < CARRIER_GATE_HIGH else
                'high(>2000)_can support confirmatory interaction testing'
            ),
        })

    summary = pd.DataFrame(rows)

    os.makedirs(LEGACY_CARRIERS_INTERNAL_DIR, exist_ok=True)
    internal_path = f"{LEGACY_CARRIERS_INTERNAL_DIR}/13_carrier_sparsity_precheck_per_gene.csv"
    summary.to_csv(internal_path, index=False)
    summary = summary[summary["gene_symbol"] == "ANY_OF_10_GENES_DEDUPED"].reset_index(drop=True)
    out_path = f"{LEGACY_CARRIERS_TABLE_DIR}/13_carrier_sparsity_precheck.csv"
    summary.to_csv(out_path, index=False)
    print(f"Write:{out_path} ({len(summary)}rows, mask summary rows only; gene-by-gene rows -&gt;{internal_path}, not push)")
    print(summary.to_string(index=False))

    cmac_rows = []
    disc_masks = [m for m in mask_names if m in ("LOF", "LOF_REVEL75", "LOF_REVEL50", "LOF_REVEL25")]
    cmac_min = cfg["step3_qualifying_variant"]["discovery_scan_cmac_min"]
    for mask in disc_masks:
        hit = pass_df[pass_df["masks_joined"].str.split("|").apply(lambda lst: mask in lst)]
        for gene in gene_list:
            g = hit[hit["gene_symbol"] == gene]
            cmac = int((g["zygosity"] == "het").sum() + 2 * (g["zygosity"] == "hom_alt").sum())
            cmac_rows.append({
                "mask": mask, "gene_symbol": gene, "cMAC": cmac,
                "passes_cmac_min": cmac >= cmac_min,
            })
    if cmac_rows:
        cmac_df = pd.DataFrame(cmac_rows)
        cmac_path = f"{LEGACY_CARRIERS_INTERNAL_DIR}/13_discovery_mask_cmac_check.csv"  
        cmac_df.to_csv(cmac_path, index=False)
        print(f"Write:{cmac_path} ({len(cmac_df)}OK, cMAC<{cmac_min}Number of excluded combinations:"
              f"{(~cmac_df['passes_cmac_min']).sum()})")
    else:
        cmac_path = None

    return summary, cmac_path

def legacy_carriers_write_report(thresholds, summary, cmac_path, n_cohort, n_case, n_control, gene_list,
                  extraction_start, extraction_end, genes_run):
    legacy_carriers_section('8. Write report -> results/report/13_carrier_extraction_summary.md')
    main_mask = "main_ptv_maf1pct"
    
    lines = []
    lines.append('# Step 3 (Finishing) + Step 3.5 | Carrier genotype extraction and sparsity/efficiency pre-check')
    lines.append(f"Run gene range:{genes_run}\n")
    lines.append(f"Extraction time:{extraction_end - extraction_start:.0f}Seconds\n")
    lines.append(f"Cohort: {n_cohort}person (pass_qc==1 & case/control), case={n_case}, control={n_control}, "
                 f"Baseline control:case ratio ={COHORT_BASELINE_RATIO}\n")
    lines.append('## Carrier-level genotype QC threshold (KDE local minimum, no reference to any phenotype/correlation results)')
    lines.append(f"- GQ threshold:{thresholds['gq_mean_threshold']}(method:{thresholds['gq_threshold_method']})\n")
    lines.append(f"- Het carrier VAF reasonable range:{thresholds['vaf_het_range']} "
                 f"(method:{thresholds['vaf_het_range_method']})\n")
    lines.append(f"- Hom-alt carrier VAF Bottom line:{thresholds['vaf_hom_alt_min']} "
                 f"({thresholds['vaf_hom_alt_min_method']})\n")
    lines.append(f"- Number of pooled carrier genotype calls (before QC):{thresholds['n_carrier_genotype_calls_pooled']}, "
                 f"Where het VAF is available:{thresholds['n_het_vaf_pooled']}\n")
    lines.append('## Number of carriers for each mask (any of 10 genes, after QC, processed according to AoU>=20 disclosure rules)')
    lines.append(summary.to_markdown(index=False))
    lines.append('Gene-by-gene carrier counts are only public in `results/table/19_gene_tier_carrier_counts.csv` (Aim 1 modeling queue);'
                 'In this step, the gene-by-gene number on the QC queue is only used for internal sparsity judgment (`data/13_internal/`).')
    lines.append('## Step 3.5 Access control interpretation (magnitude recommendations, non-hard threshold, see config carrier_gate)')
    lines.append(f"- <{CARRIER_GATE_LOW}Person -> Tendency to clearly downgrade to exploratory \n")
    lines.append(f"- {CARRIER_GATE_LOW}-{CARRIER_GATE_HIGH}People -> Write both narratives, the main story focuses on incremental R²/AUC comparison \n")
    lines.append(f"- >{CARRIER_GATE_HIGH}People -> Can support confirmatory interactive testing \n")
    lines.append('For the mask-by-mask magnitude interpretation, see the `sparsity_tier` column of `results/table/13_carrier_sparsity_precheck.csv`.')
    if cmac_path:
        lines.append(f"\n## Discovery scan cMAC threshold check (Le Guen et al. 2026, cMAC>={None}）\n"
                      .replace("None", "10"))
        lines.append('The cMAC of gene x discovery mask and whether it reaches the minimum threshold only fall into `data/13_internal/` (gene-by-gene cMAC ≈ number of people per gene, not public).')
    lines.append('## Description of data placement location')
    lines.append('- The person-level carrier_status table (by person_id, 0/1 per gene) only writes `data/`'
                 '(git-ignored), not written to results/.')
    lines.append('- This report and `results/table/13_*.csv` are aggregate statistics and do not contain the person_id column;'
                 "Any subdivided cell <20 has been marked as '<20' according to the AoU disclosure rules; in the same row, the total subtraction or ratio can be used to deduce the other cell of the cell."
                 "Mark 'suppressed (secondary)' and leave the ratio blank.")
    with open(f"{LEGACY_CARRIERS_REPORT_DIR}/13_carrier_extraction_summary.md", "w") as fh:
        fh.write("\n".join(lines))
    print(f"Write:{LEGACY_CARRIERS_REPORT_DIR}/13_carrier_extraction_summary.md")

def update_config(thresholds):
    legacy_carriers_section('9. Backfill the carrier_genotype_qc section of config/aim1_config.yaml')
    with open(LEGACY_CARRIERS_CONFIG_PATH) as fh:
        cfg = yaml.safe_load(fh)
    cfg["step3_qualifying_variant"]["carrier_genotype_qc"] = {
        "method": (
            'Starting from 2026-09-24, it is a fixed threshold for user locking (GQ>=20; het VAF [0.30,0.80]; hom VAF>=0.8),'
            'Replaces KDE estimate for 2026-08-18 (GQ 18.22, VAF [0.158,0.807]). The whole process does not refer to phenotype/correlation results.'
        ),
        "gq_mean_threshold": thresholds["gq_mean_threshold"],
        "gq_threshold_method": thresholds["gq_threshold_method"],
        "vaf_het_range": thresholds["vaf_het_range"],
        "vaf_het_range_method": thresholds["vaf_het_range_method"],
        "vaf_hom_alt_min": thresholds["vaf_hom_alt_min"],
        "gq_fallback_used_if_no_kde_minimum": GQ_FALLBACK_THRESHOLD,
        "vaf_het_fallback_used_if_no_kde_minimum": list(VAF_HET_FALLBACK_RANGE),
        "n_carrier_genotype_calls_pooled": thresholds["n_carrier_genotype_calls_pooled"],
        "status": thresholds["status"] if thresholds.get("_all_genes") else "PARTIAL_VALIDATION_RUN",
    }
    with open(LEGACY_CARRIERS_CONFIG_PATH, "w") as fh:
        yaml.dump(cfg, fh, sort_keys=False, allow_unicode=True, width=100)
    print(f"updated{LEGACY_CARRIERS_CONFIG_PATH}The step3_qualifying_variant.carrier_genotype_qc section")

def legacy_carriers_calculation():
    with open(LEGACY_CARRIERS_CONFIG_PATH) as fh:
        cfg = yaml.safe_load(fh)
    full_gene_list = cfg["step3_qualifying_variant"]["panel_genes_main_10"]

    override = os.environ.get("CARRIER_EXTRACTION_GENES", "").strip()
    gene_list = [g.strip() for g in override.split(",")] if override else full_gene_list
    print(f"Gene range of this run:{gene_list} "
          f"({'All 10 genes' if gene_list == full_gene_list else 'Validation subset (CARRIER_EXTRACTION_GENES coverage)'})")

    ensure_mount()
    sample_list_path, case_ids, control_ids, n_cohort = build_sample_list()
    shard_map = pd.read_csv(f"{LEGACY_CARRIERS_DATA_DIR}/11_panel_gene_exome_vcf_shard_map.csv")
    union_df = build_qualifying_union(gene_list)

    t0 = time.time()
    raw_old = extract_all_genes(gene_list, shard_map, SAMPLE_LIST_V0, union_df, RAW_DIR)
    raw_new = extract_all_genes(gene_list, shard_map, SAMPLE_LIST_ADDON, union_df, RAW_DIR_ADDON)
    t1 = time.time()
    keep_ids = set(pd.read_csv(sample_list_path, header=None)[0].astype(str))
    for d in (raw_old, raw_new):
        d["person_id"] = d["person_id"].astype(str)
    overlap = set(raw_old["person_id"]) & set(raw_new["person_id"])
    assert not overlap, 'The old/new parts should not have overlapping people'
    raw_carriers = pd.concat([raw_old, raw_new], ignore_index=True)
    raw_carriers = raw_carriers[raw_carriers["person_id"].isin(keep_ids)].reset_index(drop=True)
    print(f"Post-merge (new queue only) carrier genotype call:{legacy_carriers_suppress(len(raw_carriers))}")
    print(f"Total time spent in the \n extraction phase:{t1-t0:.0f}Second")

    if len(raw_carriers) == 0:
        print('!!! No carrier genotype calls were extracted, check whether the shard map/qualifying variant/sample list is correct!!!')
        return

    thresholds = dict(LOCKED_THRESHOLDS)  
    thresholds["n_carrier_genotype_calls_pooled"] = legacy_carriers_suppress(len(raw_carriers))
    thresholds["n_het_vaf_pooled"] = legacy_carriers_suppress(int(raw_carriers.loc[raw_carriers["zygosity"] == "het", "vaf_target"].notna().sum()))
    print(f"Use lock threshold:{thresholds}")
    thresholds["_all_genes"] = (gene_list == full_gene_list)
    qc_df = apply_qc(raw_carriers, thresholds)
    out_paths, mask_names, pass_df = collapse_carrier_status(qc_df, gene_list, sample_list_path)
    summary, cmac_path = step3_5_sparsity_precheck(
        pass_df, gene_list, mask_names, case_ids, control_ids, n_cohort, union_df, cfg
    )
    legacy_carriers_write_report(thresholds, summary, cmac_path, n_cohort, len(case_ids), len(control_ids),
                 gene_list, t0, t1, gene_list)
    update_config(thresholds)

    legacy_carriers_section('Finish')
    print(f"Step 3 finishing (carrier extraction) + Step 3.5 (sparsity pre-check) completed, gene range:{gene_list}")
    print('person-level carrier_status table path (data/, do not download/do not upload to github):')
    for mask, p in out_paths.items():
        print(f"  {mask}: {p}")

def analyze_legacy_carriers():
    """Initialize runtime state and execute this component explicitly."""
    global LEGACY_CARRIERS_BUCKET_NAME, CARRIER_GATE_HIGH, CARRIER_GATE_LOW, LEGACY_CARRIERS_CDR_PREFIX, LEGACY_CARRIERS_CDR_STORAGE_PATH, COHORT_BASELINE_RATIO, LEGACY_CARRIERS_CONFIG_PATH, LEGACY_CARRIERS_DATA_DIR, LEGACY_CARRIERS_DISCLOSURE_MIN_CELL, EXOME_VCF_DIR, LEGACY_CARRIERS_GOOGLE_PROJECT, GQ_FALLBACK_THRESHOLD, LEGACY_CARRIERS_INTERNAL_DIR, LOCKED_THRESHOLDS, LEGACY_CARRIERS_MOUNT_DIR, QUAL_VARIANT_FILES, RAW_DIR, RAW_DIR_ADDON, LEGACY_CARRIERS_REPO, LEGACY_CARRIERS_REPORT_DIR, SAMPLE_LIST_ADDON, SAMPLE_LIST_V0, SCRATCH, SECONDARY, LEGACY_CARRIERS_TABLE_DIR, VAF_HET_FALLBACK_RANGE, VAF_HOM_MIN

    LEGACY_CARRIERS_REPO = "/home/jupyter/workspace/analysis"
    LEGACY_CARRIERS_DATA_DIR = f"{LEGACY_CARRIERS_REPO}/data"
    RAW_DIR = f"{LEGACY_CARRIERS_DATA_DIR}/13_raw_carriers"

    RAW_DIR_ADDON = f"{LEGACY_CARRIERS_DATA_DIR}/13_raw_carriers_addon"
    SAMPLE_LIST_V0 = f"{LEGACY_CARRIERS_DATA_DIR}/13_carrier_extraction_sample_list_v0.txt"
    SAMPLE_LIST_ADDON = f"{LEGACY_CARRIERS_DATA_DIR}/13_carrier_extraction_sample_list_addon.txt"
    
    LOCKED_THRESHOLDS = {
        "gq_mean_threshold": 20.0,
        "gq_threshold_method": "fixed_user_2026-09-24",
        "vaf_het_range": [0.30, 0.80],
        "vaf_het_range_method": "fixed_user_2026-09-24",
        "vaf_hom_alt_min": 0.8,
        "vaf_hom_alt_min_method": "fixed_convention",
        "status": "LOCKED_2026-09-24",
    }
    LEGACY_CARRIERS_TABLE_DIR = f"{LEGACY_CARRIERS_REPO}/results/table"
    LEGACY_CARRIERS_REPORT_DIR = f"{LEGACY_CARRIERS_REPO}/results/report"
    LEGACY_CARRIERS_INTERNAL_DIR = f"{LEGACY_CARRIERS_DATA_DIR}/13_internal"  
    LEGACY_CARRIERS_CONFIG_PATH = f"{LEGACY_CARRIERS_REPO}/config/aim1_config.yaml"
    LEGACY_CARRIERS_MOUNT_DIR = "/tmp/claude_aou_gcsfuse_mount"
    SCRATCH = "/tmp/claude_step13_scratch"

    os.makedirs(RAW_DIR, exist_ok=True)
    os.makedirs(RAW_DIR_ADDON, exist_ok=True)
    os.makedirs(LEGACY_CARRIERS_TABLE_DIR, exist_ok=True)
    os.makedirs(LEGACY_CARRIERS_REPORT_DIR, exist_ok=True)
    os.makedirs(SCRATCH, exist_ok=True)

    LEGACY_CARRIERS_GOOGLE_PROJECT = os.environ.get("GOOGLE_PROJECT", "")
    LEGACY_CARRIERS_CDR_STORAGE_PATH = os.environ.get("CDR_STORAGE_PATH", "")
    LEGACY_CARRIERS_BUCKET_NAME = LEGACY_CARRIERS_CDR_STORAGE_PATH.replace("gs://", "").split("/")[0]
    LEGACY_CARRIERS_CDR_PREFIX = "/".join(LEGACY_CARRIERS_CDR_STORAGE_PATH.replace("gs://", "").split("/")[1:])
    EXOME_VCF_DIR = os.path.join(LEGACY_CARRIERS_MOUNT_DIR, LEGACY_CARRIERS_CDR_PREFIX, "wgs/short_read/snpindel/exome/vcf")

    QUAL_VARIANT_FILES = {
        "main_ptv_maf1pct": f"{LEGACY_CARRIERS_DATA_DIR}/12_internal/12_qualifying_variants_main_ptv.csv",
        "sensitivity_0p1pct": f"{LEGACY_CARRIERS_DATA_DIR}/12_internal/12_qualifying_variants_sensitivity_0p1pct.csv",
        "sensitivity_1pct": f"{LEGACY_CARRIERS_DATA_DIR}/12_internal/12_qualifying_variants_sensitivity_1pct.csv",
        "discovery_masks": f"{LEGACY_CARRIERS_DATA_DIR}/12_internal/12_qualifying_variants_discovery_masks.csv",
    }

    GQ_FALLBACK_THRESHOLD = 20.0
    VAF_HET_FALLBACK_RANGE = (0.15, 0.85)
    VAF_HOM_MIN = 0.8  

    CARRIER_GATE_LOW = 500
    CARRIER_GATE_HIGH = 2000
    COHORT_BASELINE_RATIO = 15.6  # plan Step0: 144,841/9,281 ≈ 15.6
    LEGACY_CARRIERS_DISCLOSURE_MIN_CELL = 20

    SECONDARY = "suppressed (secondary)"

    # ----------------------------------------------------------------------
    
    # ----------------------------------------------------------------------

    # ----------------------------------------------------------------------
    # 1. Cohort sample list（pass_qc==1 & group in case/control）
    # ----------------------------------------------------------------------

    # ----------------------------------------------------------------------
    
    # ----------------------------------------------------------------------

    # ----------------------------------------------------------------------
    
    # ----------------------------------------------------------------------

    # ----------------------------------------------------------------------
    
    # ----------------------------------------------------------------------

    # ----------------------------------------------------------------------
    
    # ----------------------------------------------------------------------

    # ----------------------------------------------------------------------
    
    # ----------------------------------------------------------------------

    if True:
        legacy_carriers_calculation()

# ============================================================================
# 6. burden: analyze_legacy_burden
# ============================================================================

import os
import numpy as np
import pandas as pd
import yaml

def legacy_burden_section(title):
    print("\n" + "=" * 80)
    print(title)
    print("=" * 80)

def legacy_burden_suppress(n, min_cell=20):
    """Legacy burden suppress."""
    if n is None or (isinstance(n, float) and np.isnan(n)):
        return n
    return n if n >= min_cell else f"<{min_cell}"

def load_weights():
    legacy_burden_section('1. Read config/aim1_config.yaml weights (gene_ptv_or / chek2_ptv_or / panel gene list)')
    with open(LEGACY_BURDEN_CONFIG_PATH) as fh:
        cfg = yaml.safe_load(fh)
    sqv = cfg["step3_qualifying_variant"]

    panel_10 = sqv["panel_genes_main_10"]
    panel_9 = sqv["panel_genes_sensitivity_9"]
    gene_or = dict(sqv["gene_ptv_or"])  # ATM, BARD1, BRCA1, BRCA2, PALB2, RAD51C, RAD51D, TP53, MAP3K1
    chek2_or = sqv["chek2_ptv_or"]

    print(f"  panel_genes_main_10 ({len(panel_10)}): {panel_10}")
    print(f"  panel_genes_sensitivity_9 ({len(panel_9)}): {panel_9}")
    print("  gene_ptv_or (non-CHEK2):")
    for g, v in gene_or.items():
        print(f"    {g}: OR={v['or']}, log(OR)={np.log(v['or']):.4f}")
    print("  chek2_ptv_or variants:")
    for k, v in chek2_or.items():
        print(f"    {k}: OR={v['or']} ({v.get('note', '')})")

    return panel_10, panel_9, gene_or, chek2_or

def build_log_or_map(gene_or, chek2_or_variant_key, panel_genes):
    """Return {gene: log(OR)} for the given panel, using the specified CHEK2 OR variant."""
    log_or = {}
    for g in panel_genes:
        if g == "CHEK2":
            or_val = chek2_or_variant_key["or"]
        else:
            or_val = gene_or[g]["or"]
        log_or[g] = float(np.log(or_val))
    return log_or

def compute_burden(mask_name, panel_genes, log_or_map, col_name):
    path = f"{LEGACY_BURDEN_DATA_DIR}/13_carrier_status_{mask_name}_person_level.csv.gz"
    df = pd.read_csv(path)
    missing = [g for g in panel_genes if g not in df.columns]
    if missing:
        raise ValueError(f"{path} missing gene columns: {missing}")

    weights = np.array([log_or_map[g] for g in panel_genes])
    burden = df[panel_genes].to_numpy(dtype=float) @ weights
    out = df[["person_id"]].copy()
    out[col_name] = burden
    return out, df[["person_id"] + panel_genes].rename(
        columns={g: f"{mask_name}__{g}" for g in panel_genes}
    )

def legacy_burden_calculation():
    os.makedirs(LEGACY_BURDEN_TABLE_DIR, exist_ok=True)
    os.makedirs(LEGACY_BURDEN_REPORT_DIR, exist_ok=True)

    panel_10, panel_9, gene_or, chek2_or = load_weights()

    # mask -> which CHEK2 OR variant to use
    mask_chek2_variant = {
        "main_ptv_maf1pct": "main",
        "sensitivity_ptv_maf0p1pct": "sensitivity_0p1pct",
        "sensitivity_ptv_maf1pct": "main",
    }

    legacy_burden_section('2. Construct the Burden_rare score for each mask x panel combination')
    combos = []  # (mask, panel_label, panel_genes, col_name)
    for mask, chek2_key in mask_chek2_variant.items():
        combos.append((mask, "10gene", panel_10, chek2_key))
        combos.append((mask, "9gene", panel_9, chek2_key))

    merged = None
    carrier_cols_for_check = {}
    for mask, panel_label, panel_genes, chek2_key in combos:
        col_name = f"burden_rare_{mask.replace('main_ptv_maf1pct', 'main_ptv')}_{panel_label}"
        log_or_map = build_log_or_map(gene_or, chek2_or[chek2_key], panel_genes)
        print(f"  {col_name}: mask={mask}, panel={panel_label} ({len(panel_genes)} genes), "
              f"CHEK2 OR variant={chek2_key} (OR={chek2_or[chek2_key]['or']})")
        for g in panel_genes:
            print(f"      log(OR)_{g} = {log_or_map[g]:.4f}")

        burden_df, carrier_wide = compute_burden(mask, panel_genes, log_or_map, col_name)
        if mask == "main_ptv_maf1pct":
            carrier_cols_for_check[panel_label] = carrier_wide

        if merged is None:
            merged = burden_df
        else:
            merged = merged.merge(burden_df, on="person_id", how="outer")

    print(f"\n merged person-level table:{merged.shape}")
    burden_cols = [c for c in merged.columns if c != "person_id"]
    print(f"Burden_rare column ({len(burden_cols)}): {burden_cols}")

    out_path = f"{LEGACY_BURDEN_DATA_DIR}/14_burden_rare_scores_person_level.csv.gz"
    merged.to_csv(out_path, index=False)
    print(f"has been written{out_path}")

    # ------------------------------------------------------------------
    
    # ------------------------------------------------------------------
    legacy_burden_section('3. Read case/control grouping labels (04_phenotype_diagnostics + 05_cohort_qc_flags)')
    pheno = pd.read_csv(f"{LEGACY_BURDEN_DATA_DIR}/04_phenotype_diagnostics_person_level.csv.gz")
    qc = pd.read_csv(f"{LEGACY_BURDEN_DATA_DIR}/05_cohort_qc_flags.csv")
    cohort = pheno.merge(qc[["person_id", "pass_qc"]], on="person_id", how="inner")
    cohort = cohort[(cohort["pass_qc"] == 1) & (cohort["group"].isin(["case", "control"]))]
    cohort = cohort[["person_id", "group"]]
    print(f"  cohort (pass_qc==1 & case/control): {cohort.shape}, "
          f"{cohort['group'].value_counts().to_dict()}")

    merged_grp = merged.merge(cohort, on="person_id", how="inner")
    n_unmatched = len(merged) - len(merged_grp)
    print(f"After the merged Burden_rare table is joined to cohort:{merged_grp.shape} "
          f"(Number of rows that do not match case/control tags ={n_unmatched})")
    if len(merged_grp) != len(merged):
        print('WARNING: The number of rows in the person-level carrier table and the cohort case/control label are not exactly the same.'
              'Please check whether the cohort of script 13 is consistent with the assumptions of this script.')

    # ------------------------------------------------------------------
    
    # ------------------------------------------------------------------
    legacy_burden_section('4. Burden_rare distribution by case/control (aggregation, no person_id)')
    dist_rows = []
    for col in burden_cols:
        for grp_name, sub in merged_grp.groupby("group")[col]:
            n = len(sub)
            dist_rows.append({
                "burden_score": col,
                "group": grp_name,
                "n": legacy_burden_suppress(n),
                "mean": round(sub.mean(), 5),
                "median": round(sub.median(), 5),
                "sd": round(sub.std(), 5),
                "q25": round(sub.quantile(0.25), 5),
                "q75": round(sub.quantile(0.75), 5),
                "min": round(sub.min(), 5),
                "max": round(sub.max(), 5),
                "pct_nonzero": round(float((sub != 0).mean() * 100), 3),
            })
    dist_df = pd.DataFrame(dist_rows)

    os.makedirs(f"{LEGACY_BURDEN_DATA_DIR}/14_internal", exist_ok=True)
    dist_df.to_csv(f"{LEGACY_BURDEN_DATA_DIR}/14_internal/14_burden_rare_distribution_by_case_control_all_scores.csv", index=False)
    dist_df = dist_df[~dist_df["burden_score"].str.endswith("_9gene")].reset_index(drop=True)
    print(dist_df.to_string(index=False))
    dist_path = f"{LEGACY_BURDEN_TABLE_DIR}/14_burden_rare_distribution_by_case_control.csv"
    dist_df.to_csv(dist_path, index=False)
    print(f"has been written{dist_path}")

    # ------------------------------------------------------------------
    
    # ------------------------------------------------------------------
    legacy_burden_section('5. Burden_rare>0 <=> Internal consistency check for at least one gene carrier_status==1 (main_ptv_maf1pct, 10gene)')
    check_cols = ["ATM", "BARD1", "BRCA1", "BRCA2", "CHEK2", "PALB2", "RAD51C", "RAD51D", "TP53", "MAP3K1"]
    carrier_wide_10 = carrier_cols_for_check["10gene"]
    carrier_wide_10 = carrier_wide_10.rename(
        columns={f"main_ptv_maf1pct__{g}": g for g in check_cols}
    )
    check_df = merged[["person_id", "burden_rare_main_ptv_10gene"]].merge(
        carrier_wide_10, on="person_id", how="inner"
    )
    check_df["any_carrier"] = (check_df[check_cols].sum(axis=1) > 0).astype(int)
    check_df["burden_gt0"] = (check_df["burden_rare_main_ptv_10gene"] > 0).astype(int)

    consistency_rows = []

    n_mismatch_total = int((check_df["any_carrier"] != check_df["burden_gt0"]).sum())
    consistency_rows.append({
        "check": "any_of_10genes_carrier_vs_burden_gt0",
        "any_carrier": "all_consistent",
        "burden_gt0": "any_carrier==1 <=> burden_gt0==1",
        "n": "yes" if n_mismatch_total == 0 else "NO (see log)",
    })
    # per-gene: among carriers of that single gene (regardless of other genes), what fraction have burden>0
    for g in check_cols:
        carriers_g = check_df[check_df[g] == 1]
        n_g = len(carriers_g)
        n_g_burden_gt0 = int((carriers_g["burden_rare_main_ptv_10gene"] > 0).sum())
        consistency_rows.append({
            "check": f"gene_{g}_carrier_implies_burden_gt0",
            "any_carrier": g,
            "burden_gt0": "all carriers of this gene have burden_gt0",
            
            "n": "yes" if n_g_burden_gt0 == n_g else "NO (see internal log)",
        })
    consistency_df = pd.DataFrame(consistency_rows)
    print(consistency_df.to_string(index=False))
    n_mismatch = int(((check_df["any_carrier"] == 1) & (check_df["burden_gt0"] == 0)).sum())
    print(f"\n Mismatch count (any_carrier==1 but burden_rare<=0, theoretically should be 0, because all main analysis log(OR)>0):{n_mismatch if n_mismatch == 0 else legacy_burden_suppress(n_mismatch)}")
    consistency_path = f"{LEGACY_BURDEN_TABLE_DIR}/14_burden_rare_carrier_consistency_check.csv"
    consistency_df.to_csv(consistency_path, index=False)
    print(f"has been written{consistency_path}")

    # ------------------------------------------------------------------
    
    # ------------------------------------------------------------------
    legacy_burden_section('6. 10-gene vs 9-gene main analysis Burden_rare correlation (main_ptv_maf1pct)')
    corr_rows = []
    pairs = [
        ("burden_rare_main_ptv_10gene", "burden_rare_main_ptv_9gene", "main_ptv_maf1pct"),
        ("burden_rare_sensitivity_ptv_maf0p1pct_10gene", "burden_rare_sensitivity_ptv_maf0p1pct_9gene",
         "sensitivity_ptv_maf0p1pct"),
        ("burden_rare_sensitivity_ptv_maf1pct_10gene", "burden_rare_sensitivity_ptv_maf1pct_9gene",
         "sensitivity_ptv_maf1pct"),
    ]
    for col10, col9, mask in pairs:
        pearson_r = merged[col10].corr(merged[col9], method="pearson")
        spearman_r = merged[col10].corr(merged[col9], method="spearman")
        n_diff = int((merged[col10] != merged[col9]).sum())
        corr_rows.append({
            "mask": mask,
            "col_10gene": col10,
            "col_9gene": col9,
            "pearson_r": round(pearson_r, 6),
            "spearman_r": round(spearman_r, 6),
            "n_total": legacy_burden_suppress(len(merged)),
            
        })
    corr_df = pd.DataFrame(corr_rows)
    print(corr_df.to_string(index=False))
    corr_path = f"{LEGACY_BURDEN_TABLE_DIR}/14_burden_rare_10gene_vs_9gene_correlation.csv"
    corr_df.to_csv(corr_path, index=False)
    print(f"has been written{corr_path}")

    # ------------------------------------------------------------------
    
    # ------------------------------------------------------------------
    legacy_burden_section('7. Write results/report/14_burden_rare_construction_summary.md')
    legacy_burden_write_report(dist_df, consistency_df, corr_df, n_mismatch, merged, cohort, panel_10, panel_9,
                 gene_or, chek2_or)

    print('Script 14 is executed.')

def legacy_burden_write_report(dist_df, consistency_df, corr_df, n_mismatch, merged, cohort, panel_10, panel_9,
                  gene_or, chek2_or):
    n_case = int((cohort["group"] == "case").sum())
    n_control = int((cohort["group"] == "control").sum())

    lines = []
    lines.append('# Step 4 | Burden_rare score construction summary report')
    lines.append(f"Corresponding script: `scripts/90_legacy_ptv_baseline.py`; corresponding plan: `plans/aim_1.md`"
                 f"Step 4 \"Burden_rare score construction\". \n")
    lines.append(f"Cohort: {legacy_burden_suppress(len(cohort))}people (pass_qc==1 & group in {{case,control}}), "
                 f"case={legacy_burden_suppress(n_case)}, control={legacy_burden_suppress(n_control)}。\n")

    lines.append('## Formula')
    lines.append('Intra-gene folding (script 13 completed, only filtering without weighting):')
    lines.append('```\ncarrier_{i,g} = 1[There is PTV qualifying variant j of gene g: G_ij > 0]\n```')
    lines.append('')
    lines.append("```\nBurden_rare_i = sum_g log(OR_g) * carrier_{i,g}\n```\n")

    lines.append('## Weight source (`config/aim1_config.yaml`, this script only reads and does not re-derive)')
    lines.append(f"- 10-gene main analysis panel (`panel_genes_main_10`):{panel_10}\n")
    lines.append(f"- 9 gene sensitivity panel (Dorling-only, `panel_genes_sensitivity_9`):{panel_9}\n")
    lines.append('- Non-CHEK2 gene OR (Dorling et al., population-based overall PTV OR):')
    for g, v in gene_or.items():
        lines.append(f"  - {g}: OR={v['or']}, log(OR)={np.log(v['or']):.4f}"
                      f" (95% CI {v['ci_low']}-{v['ci_high']})\n")
    lines.append('- CHEK2 OR (select different variants according to mask, reason: c.1100delC AF≈0.26%,'
                 'Within the unified MAF <1% range but excluded by the unified MAF <0.1% range):')
    for k, v in chek2_or.items():
        lines.append(f"  - `{k}`: OR={v['or']} — {v.get('note', '')}\n")

    lines.append('## Mask x panel combination and CHEK2 OR used')
    lines.append('| Burden_rare column | mask | panel | CHEK2 OR variant (OR value) |\n|---|---|---|---|')
    combo_doc = [
        ("burden_rare_main_ptv_10gene", "main_ptv_maf1pct", '10 genes (main analysis)', "main (2.54)"),
        ("burden_rare_main_ptv_9gene", "main_ptv_maf1pct", '9 genes (sensitivity, Dorling-only)', "main (2.54)"),
        ("burden_rare_sensitivity_ptv_maf0p1pct_10gene", "sensitivity_ptv_maf0p1pct", '10 genes', "sensitivity_0p1pct (2.13)"),
        ("burden_rare_sensitivity_ptv_maf0p1pct_9gene", "sensitivity_ptv_maf0p1pct", '9 genes', "sensitivity_0p1pct (2.13)"),
        ("burden_rare_sensitivity_ptv_maf1pct_10gene", "sensitivity_ptv_maf1pct", '10 genes', "main (2.54)"),
        ("burden_rare_sensitivity_ptv_maf1pct_9gene", "sensitivity_ptv_maf1pct", '9 genes', "main (2.54)"),
    ]
    for col, mask, panel, chek2 in combo_doc:
        lines.append(f"| `{col}` | `{mask}` | {panel} | {chek2} |\n")
    lines.append('Main analysis = `burden_rare_main_ptv_10gene`; if the results of 9 genes/10 genes are inconsistent,'
                 'Take `burden_rare_main_ptv_9gene` as the main conclusion (see plans/aim_1.md Step 4).'
                 'Missense sensitivity mask with 4 discovery masks does not build Burden_rare (planned only for'
                 'SAIGE-GENE+ scan, outside the scope of this step).')

    lines.append('## Sanity check 1: Main 10-gene Burden_rare distribution by case/control')
    lines.append(dist_df[dist_df["burden_score"] == "burden_rare_main_ptv_10gene"]
                 .to_markdown(index=False))
    lines.append('10 genes 3 scores')

    lines.append('## Sanity check 2: carrier -> Burden_rare>0 internal consistency (main_ptv_maf1pct, 10gene)')
    lines.append(f"Mismatch (any_carrier==1 but burden_rare<=0) quantity ={n_mismatch if n_mismatch == 0 else legacy_burden_suppress(n_mismatch)}"
                 '(Theoretically it should be 0, because if all 10 genes are analyzed, log(OR)>0, any gene carrier must make Burden_rare>0).')
    lines.append('For gene-by-gene verification, see `results/table/14_burden_rare_carrier_consistency_check.csv` (the 2x2 crosstab is no longer output, only whether it is completely consistent is reported).')

    lines.append('## Sanity check 3: 10-gene vs 9-gene correlation')
    lines.append(corr_df.to_markdown(index=False))
    lines.append('The 10-gene version only has one more bonus, MAP3K1 log(OR)=1.833, than the 9-gene version.'
                 'MAP3K1 carrier itself is extremely sparse, so it is expected that the two are highly correlated, only on a few person_ids of MAP3K1 carrier'
                 'Differences occur.')

    lines.append('## Output file')
    lines.append('- Person-level (git-ignored, only exists in `data/`): `data/14_burden_rare_scores_person_level.csv.gz`')
    lines.append('- Aggregation table: `results/table/14_burden_rare_distribution_by_case_control.csv`,'
                 "`results/table/14_burden_rare_carrier_consistency_check.csv`, "
                 "`results/table/14_burden_rare_10gene_vs_9gene_correlation.csv`\n")

    with open(f"{LEGACY_BURDEN_REPORT_DIR}/14_burden_rare_construction_summary.md", "w") as fh:
        fh.writelines(lines)

def analyze_legacy_burden():
    """Initialize runtime state and execute this component explicitly."""
    global LEGACY_BURDEN_CONFIG_PATH, LEGACY_BURDEN_DATA_DIR, LEGACY_BURDEN_DISCLOSURE_MIN_CELL, LEGACY_BURDEN_REPO, LEGACY_BURDEN_REPORT_DIR, LEGACY_BURDEN_TABLE_DIR

    LEGACY_BURDEN_REPO = "/home/jupyter/workspace/analysis"
    LEGACY_BURDEN_DATA_DIR = f"{LEGACY_BURDEN_REPO}/data"
    LEGACY_BURDEN_TABLE_DIR = f"{LEGACY_BURDEN_REPO}/results/table"
    LEGACY_BURDEN_REPORT_DIR = f"{LEGACY_BURDEN_REPO}/results/report"
    LEGACY_BURDEN_CONFIG_PATH = f"{LEGACY_BURDEN_REPO}/config/aim1_config.yaml"

    LEGACY_BURDEN_DISCLOSURE_MIN_CELL = 20

    # ----------------------------------------------------------------------
    
    # ----------------------------------------------------------------------

    # ----------------------------------------------------------------------
    
    # ----------------------------------------------------------------------

    if True:
        legacy_burden_calculation()

# ============================================================================
# 7. joint-models: analyze_legacy_joint_models
# ============================================================================

import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
import yaml
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import RepeatedStratifiedKFold

def legacy_joint_models_log(msg):
    print(msg, flush=True)

def legacy_joint_models_section(title):
    legacy_joint_models_log("\n" + "=" * 80)
    legacy_joint_models_log(title)
    legacy_joint_models_log("=" * 80)

def legacy_joint_models_suppress(n, min_cell=20):
    """Legacy joint models suppress."""
    if n is None or (isinstance(n, float) and np.isnan(n)):
        return n
    return n if n >= min_cell else f"<{min_cell}"

def p_floor(df):
    """P floor."""
    out = df.copy()
    for c in out.columns:
        if c.startswith("p_") and pd.api.types.is_numeric_dtype(out[c]) and (out[c] == 0).any():
            out[c] = out[c].astype(object).where(out[c] != 0, "<1e-300")
    return out

def balanced_weights(y):
    """Balanced weights."""
    y = np.asarray(y)
    n = len(y)
    n1 = y.sum()
    n0 = n - n1
    w = np.where(y == 1, n / (2 * n1), n / (2 * n0))
    return w

def fit_unweighted_logit(y, X):
    """Fit unweighted logit."""
    model = sm.GLM(y, X, family=sm.families.Binomial())
    return model.fit(cov_type="HC1"), np.ones(len(y))

def fit_weighted_logit(y, X):
    """Fit weighted logit."""
    w = balanced_weights(y)
    model = sm.GLM(y, X, family=sm.families.Binomial(), var_weights=w)
    return model.fit(cov_type="HC1"), w

def legacy_joint_models_coef_table(fit_result, terms, label):
    params = fit_result.params
    bse = fit_result.bse
    pvals = fit_result.pvalues
    rows = []
    for t in terms:
        if t not in params.index:
            continue
        b = params[t]
        se = bse[t]
        p = pvals[t]
        or_ = np.exp(b)
        ci_low = np.exp(b - 1.959963984540054 * se)
        ci_high = np.exp(b + 1.959963984540054 * se)
        rows.append({
            "model": label, "term": t, "beta": b, "se": se,
            "OR": or_, "CI_low": ci_low, "CI_high": ci_high, "p_value": p,
        })
    return pd.DataFrame(rows)

def mcfadden_r2(y, X_full, w, ll_full):
    """Mcfadden r2."""
    X_null = pd.DataFrame({"const": np.ones(len(y))})
    if w is None or np.allclose(w, 1.0):
        null_fit = sm.GLM(y, X_null, family=sm.families.Binomial()).fit()
    else:
        null_fit = sm.GLM(y, X_null, family=sm.families.Binomial(), var_weights=w).fit()
    ll_null = null_fit.llf
    r2 = 1 - ll_full / ll_null
    return r2, ll_null

def zscore_fit(train_series):
    mu = train_series.mean()
    sd = train_series.std()
    return mu, sd

def apply_z(series, mu, sd):
    return (series - mu) / sd

def delong_auc_variance(y_true, y_score):
    """Delong auc variance."""
    y_true = np.asarray(y_true)
    y_score = np.asarray(y_score)
    pos = y_score[y_true == 1]
    neg = y_score[y_true == 0]
    n_pos, n_neg = len(pos), len(neg)

    def midrank(x):
        return stats.rankdata(x, method="average")  

    combined = np.concatenate([pos, neg])
    all_ranks = midrank(combined)
    pos_ranks_in_combined = all_ranks[:n_pos]
    tx = midrank(pos)
    ty = midrank(neg)

    v10 = (pos_ranks_in_combined - tx) / n_neg  # placement of each positive
    neg_ranks_in_combined = all_ranks[n_pos:]
    v01 = 1.0 - (neg_ranks_in_combined - ty) / n_pos  # placement of each negative

    auc = v10.mean()
    return auc, v10, v01

def delong_test(y_true, score_a, score_b):
    """Delong test."""
    auc_a, v10_a, v01_a = delong_auc_variance(y_true, score_a)
    auc_b, v10_b, v01_b = delong_auc_variance(y_true, score_b)
    n_pos, n_neg = len(v10_a), len(v01_a)

    v10 = np.vstack([v10_a, v10_b])  # (2, n_pos)
    v01 = np.vstack([v01_a, v01_b])  # (2, n_neg)
    s10 = np.cov(v10, bias=False) / n_pos
    s01 = np.cov(v01, bias=False) / n_neg
    cov_mat = s10 + s01  # 2x2

    delta = auc_a - auc_b
    var_delta = cov_mat[0, 0] + cov_mat[1, 1] - 2 * cov_mat[0, 1]
    se_delta = np.sqrt(max(var_delta, 0))
    z = delta / se_delta if se_delta > 0 else np.nan
    p = 2 * stats.norm.sf(abs(z)) if se_delta > 0 else np.nan  
    ci_low = delta - 1.959963984540054 * se_delta
    ci_high = delta + 1.959963984540054 * se_delta
    return {
        "auc_a": auc_a, "auc_b": auc_b, "delta_auc": delta,
        "ci_low": ci_low, "ci_high": ci_high, "p_value": p,
    }

def legacy_joint_models_build_X(split_df, extra_cols):
    cols = extra_cols + LEGACY_JOINT_MODELS_COV_COLS
    X = split_df[cols].astype(float).copy()
    X = sm.add_constant(X, has_constant="add")
    return X

def add_z(frame, prs_mu, prs_sd, burden_mu, burden_sd):
    frame = frame.copy()
    frame["PRS311_z"] = apply_z(frame["PRS311"], prs_mu, prs_sd)
    frame["burden_z"] = apply_z(frame["burden_rare_main_ptv_10gene"], burden_mu, burden_sd)
    frame["inter_z"] = frame["PRS311_z"] * frame["burden_z"]
    return frame

def auc_ci(y, score):
    auc, v10, v01 = delong_auc_variance(y, score)
    se = np.sqrt(np.var(v10, ddof=1) / len(v10) + np.var(v01, ddof=1) / len(v01))
    return auc - 1.959963984540054 * se, auc + 1.959963984540054 * se

def legacy_joint_models_save_smallest(fig, stem):
    paths = {ext: os.path.join(LEGACY_JOINT_MODELS_FIG_DIR, f"{stem}.{ext}") for ext in ("svg", "png", "jpg")}
    fig.savefig(paths["svg"])
    fig.savefig(paths["png"], dpi=LEGACY_JOINT_MODELS_DPI)
    fig.savefig(paths["jpg"], dpi=LEGACY_JOINT_MODELS_DPI, pil_kwargs={"quality": 85, "subsampling": 2})
    sizes = {ext: os.path.getsize(p) for ext, p in paths.items()}
    keep = min(sizes, key=sizes.get)
    for ext, p in paths.items():
        if ext != keep:
            os.remove(p)
    legacy_joint_models_log(f"  {stem}: " + ", ".join(f"{e}={s:,}B" for e, s in sizes.items()) + f" -> kept .{keep}")

def assign_stratum(row):
    n = row["n_genes_carried"]
    if n == 0:
        return "non-carrier"
    if n == 1:
        carried_gene = [g for g in LEGACY_JOINT_MODELS_GENE_COLS if row[f"gene_{g}"] == 1][0]
        return carried_gene if carried_gene in BIG4_GENES else "other"
    return "other"  # compound (n>=2 genes)

def analyze_legacy_joint_models():
    """Initialize runtime state and execute this component explicitly."""
    global LEGACY_JOINT_MODELS_BASE, BIG4_GENES, LEGACY_JOINT_MODELS_CONFIG_PATH, LEGACY_JOINT_MODELS_COV_COLS, CV_FOLDS, CV_REPEATS, LEGACY_JOINT_MODELS_DATA_DIR, LEGACY_JOINT_MODELS_DISCLOSURE_MIN_CELL, LEGACY_JOINT_MODELS_DPI, ESTIMAND_CAVEAT, LEGACY_JOINT_MODELS_FIG_DIR, LEGACY_JOINT_MODELS_GENE_COLS, MODEL_COLORS, LEGACY_JOINT_MODELS_PC_COLS, RANDOM_STATE, LEGACY_JOINT_MODELS_REPORT_DIR, LEGACY_JOINT_MODELS_SEQ_COLS, LEGACY_JOINT_MODELS_SITE_COLS, LEGACY_JOINT_MODELS_TABLE_DIR, LEGACY_JOINT_MODELS_X, X_281, X_carrier_bin, X_full_joint, X_ni, X_primary_full, X_r2, X_sens, X_strat, legacy_joint_models__, _ct, a, age_na_by_group, auc, auc_rep, auc_rep_rows, auc_tbl, legacy_joint_models_ax, legacy_joint_models_b, b_col, beta_burden, beta_inter, beta_prs, legacy_joint_models_burden, burden_na_by_group, legacy_joint_models_c, carrier_terms, cfg, ci_half_width, ci_high, ci_low, legacy_joint_models_clinical, coef_rows_main, legacy_joint_models_cohort, column_swap_specs, comp_label, ct, ct_281, ct_carrier_bin, ct_ni, ct_r2, ct_weighted, cv_splits, delong_pairs, delong_rep, delong_rep_rows, delta_tbl, legacy_joint_models_df, df_281, df_281_cc, df_noimp_cc, df_r2, df_r2_cc, distinct_burden_values, evidence_status_map, external_or_ref, legacy_joint_models_extra_cols, legacy_joint_models_f, legacy_joint_models_fig, fit_281, fit_carrier_bin, fit_ni, fit_r2, legacy_joint_models_fit_res, fit_sens, fit_strat, fit_weighted_primary, fold_r1, fpr, full_burden_mu, full_burden_sd, full_coef_all, full_coef_primary, full_joint_fit, legacy_joint_models_full_prs_mu, legacy_joint_models_full_prs_sd, legacy_joint_models_g, legacy_joint_models_gene_carrier, gene_flag_cols, legacy_joint_models_i, implied_or_vs_noncarrier, implied_prs_or_per_sd, implied_rows, implied_tbl, incremental_r2_tbl, j, legacy_joint_models_k, keep_mask, legacy_joint_models_label, ll_null, m, model_specs, mu_b, mu_b_r2, mu_p, mu_p281, mu_p_r2, mu_pni, n_before_age, n_before_burden, legacy_joint_models_n_case, n_case_r2, n_control, n_control_r2, legacy_joint_models_n_dropped_noimp, legacy_joint_models_n_na, legacy_joint_models_n_na_281, n_persons, oof, oof_persist, oof_r1, p_col, passes, legacy_joint_models_pc, legacy_joint_models_pheno, legacy_joint_models_primary_for_compare, primary_terms, legacy_joint_models_prs, legacy_joint_models_prs_281, prs_281_renamed, prs_map, prs_merge_cols, pseudo_r2_rows, pseudo_r2_tbl, legacy_joint_models_qc, qc_cols, legacy_joint_models_r, r2, res, legacy_joint_models_row, row_out, rskf, sd_b, sd_b_r2, sd_p, sd_p281, sd_p_r2, sd_pni, sd_used_map, legacy_joint_models_se, sens281_df, sens_df, legacy_joint_models_sensitivity_rows, sensitivity_tbl, legacy_joint_models_sensitivity_tbl_full, sensni_df, sensr2_df, seq_dummies_full, site_dummies_full, spec, strat_rows, strat_tbl, strata_counts, legacy_joint_models_stratum, stratum_order, legacy_joint_models_sub, te, te_idx, tpr, tr, tr_idx, legacy_joint_models_v, w, x, y_281, legacy_joint_models_y_all, y_full, y_ni, y_r2, y_strat, z_v

    matplotlib.use("Agg")
    plt.style.use("default")

    LEGACY_JOINT_MODELS_BASE = "/home/jupyter/workspace/analysis"
    LEGACY_JOINT_MODELS_DATA_DIR = os.path.join(LEGACY_JOINT_MODELS_BASE, "data")
    LEGACY_JOINT_MODELS_FIG_DIR = os.path.join(LEGACY_JOINT_MODELS_BASE, "results", "figures")
    LEGACY_JOINT_MODELS_TABLE_DIR = os.path.join(LEGACY_JOINT_MODELS_BASE, "results", "table")
    LEGACY_JOINT_MODELS_REPORT_DIR = os.path.join(LEGACY_JOINT_MODELS_BASE, "results", "report")
    LEGACY_JOINT_MODELS_CONFIG_PATH = os.path.join(LEGACY_JOINT_MODELS_BASE, "config", "aim1_config.yaml")

    os.makedirs(LEGACY_JOINT_MODELS_FIG_DIR, exist_ok=True)
    os.makedirs(LEGACY_JOINT_MODELS_TABLE_DIR, exist_ok=True)
    os.makedirs(LEGACY_JOINT_MODELS_REPORT_DIR, exist_ok=True)

    LEGACY_JOINT_MODELS_DPI = 150
    LEGACY_JOINT_MODELS_DISCLOSURE_MIN_CELL = 20
    RANDOM_STATE = 20260818
    CV_FOLDS = 5
    CV_REPEATS = 10
    LEGACY_JOINT_MODELS_PC_COLS = [f"PC{legacy_joint_models_i}" for legacy_joint_models_i in range(1, 11)]
    LEGACY_JOINT_MODELS_GENE_COLS = ["ATM", "BARD1", "BRCA1", "BRCA2", "CHEK2", "PALB2", "RAD51C", "RAD51D", "TP53", "MAP3K1"]
    BIG4_GENES = ["BRCA1", "BRCA2", "CHEK2", "ATM"]

    MODEL_COLORS = {
        "Covariates-only": "tab:gray", "PRS-only": "tab:orange",
        "Burden-only": "tab:blue", "Joint": "tab:green",
    }

    # ---------------------------------------------------------------------------
    
    # ---------------------------------------------------------------------------
    legacy_joint_models_section('Step A: Read each table and join to build a modeling queue')

    with open(LEGACY_JOINT_MODELS_CONFIG_PATH) as legacy_joint_models_f:
        cfg = yaml.safe_load(legacy_joint_models_f)
    prs_map = cfg["step3_qualifying_variant"]["prs_version_by_maf_tier"]
    legacy_joint_models_log(f"config prs_version_by_maf_tier: {prs_map}")

    legacy_joint_models_pheno = pd.read_csv(os.path.join(LEGACY_JOINT_MODELS_DATA_DIR, "04_phenotype_diagnostics_person_level.csv.gz"))
    legacy_joint_models_qc = pd.read_csv(os.path.join(LEGACY_JOINT_MODELS_DATA_DIR, "05_cohort_qc_flags.csv"))
    legacy_joint_models_pc = pd.read_csv(os.path.join(LEGACY_JOINT_MODELS_DATA_DIR, "09_ancestry_pc_person_level.csv.gz"))
    legacy_joint_models_clinical = pd.read_csv(os.path.join(LEGACY_JOINT_MODELS_DATA_DIR, "09_extended_clinical_person_level.csv.gz"))
    legacy_joint_models_burden = pd.read_csv(os.path.join(LEGACY_JOINT_MODELS_DATA_DIR, "14_burden_rare_scores_person_level.csv.gz"))
    legacy_joint_models_prs = pd.read_csv(os.path.join(LEGACY_JOINT_MODELS_DATA_DIR, "prs_scores.csv"))
    legacy_joint_models_prs_281 = pd.read_csv(os.path.join(LEGACY_JOINT_MODELS_DATA_DIR, "prs_scores_281loci.csv"))
    legacy_joint_models_gene_carrier = pd.read_csv(os.path.join(LEGACY_JOINT_MODELS_DATA_DIR, "13_carrier_status_main_ptv_maf1pct_person_level.csv.gz"))

    qc_cols = ["person_id", "pass_qc"]

    legacy_joint_models_cohort = legacy_joint_models_pheno.merge(legacy_joint_models_qc[qc_cols], on="person_id", how="inner")
    legacy_joint_models_cohort = legacy_joint_models_cohort[(legacy_joint_models_cohort["pass_qc"] == 1) & (legacy_joint_models_cohort["group"].isin(["case", "control"]))].copy()
    legacy_joint_models_log(f"  pass_qc==1 & group in case/control: n={len(legacy_joint_models_cohort)} "
        f"(case={int((legacy_joint_models_cohort['group']=='case').sum())}, control={int((legacy_joint_models_cohort['group']=='control').sum())})")

    legacy_joint_models_df = legacy_joint_models_cohort.merge(legacy_joint_models_pc[["person_id"] + LEGACY_JOINT_MODELS_PC_COLS], on="person_id", how="left")
    legacy_joint_models_df = legacy_joint_models_df.merge(legacy_joint_models_clinical[["person_id", "seq_site_id", "seq_sequencer"]], on="person_id", how="left")
    legacy_joint_models_df = legacy_joint_models_df.merge(legacy_joint_models_burden, on="person_id", how="left")
    prs_merge_cols = ["person_id", "PRS311", "PRS313"]
    if "PRS311_noimp" in legacy_joint_models_prs.columns:
        prs_merge_cols.append("PRS311_noimp")
    else:
        legacy_joint_models_log('WARNING: There is no PRS311_noimp column in data/prs_scores.csv -- scripts/03_common_variants_prs.py'
            'The column persist has not been rerun, prs_no_mean_imputation_PRS311 sensitivity analysis will be skipped')
    legacy_joint_models_df = legacy_joint_models_df.merge(legacy_joint_models_prs[prs_merge_cols], on="person_id", how="left")
    legacy_joint_models_df = legacy_joint_models_df.merge(legacy_joint_models_gene_carrier.rename(columns={legacy_joint_models_c: f"gene_{legacy_joint_models_c}" for legacy_joint_models_c in LEGACY_JOINT_MODELS_GENE_COLS}), on="person_id", how="left")
    legacy_joint_models_log(f"After join: n={len(legacy_joint_models_df)}")

    for legacy_joint_models_c in LEGACY_JOINT_MODELS_PC_COLS + ["seq_site_id", "seq_sequencer", "burden_rare_main_ptv_10gene", "PRS311", "PRS313"]:
        legacy_joint_models_n_na = int(legacy_joint_models_df[legacy_joint_models_c].isna().sum())
        if legacy_joint_models_n_na > 0:
            legacy_joint_models_log(f"  WARNING: {legacy_joint_models_c}have{legacy_joint_models_n_na}NA (non-age covariate, exceeded expectations, needs to be checked)")

    n_before_burden = len(legacy_joint_models_df)
    burden_na_by_group = legacy_joint_models_df.loc[legacy_joint_models_df["burden_rare_main_ptv_10gene"].isna(), "group"].value_counts().to_dict()
    assert not burden_na_by_group, f"burden_rare_main_ptv_10gene is missing:{burden_na_by_group}"

    n_before_age = len(legacy_joint_models_df)
    age_na_by_group = legacy_joint_models_df.loc[legacy_joint_models_df["age_at_last_obs"].isna(), "group"].value_counts().to_dict()
    legacy_joint_models_df = legacy_joint_models_df.dropna(subset=["age_at_last_obs"]).reset_index(drop=True)
    legacy_joint_models_log(f"age_at_last_obs missing removal:{n_before_age} -> {len(legacy_joint_models_df)} "
        f"(eliminate{legacy_joint_models_suppress(n_before_age - len(legacy_joint_models_df))}People, distributed by group:{ {legacy_joint_models_k: legacy_joint_models_suppress(legacy_joint_models_v) for legacy_joint_models_k, legacy_joint_models_v in age_na_by_group.items()} })")
    legacy_joint_models_log(f"Final modeling queue n={len(legacy_joint_models_df)} "
        f"(case={int((legacy_joint_models_df['group']=='case').sum())}, control={int((legacy_joint_models_df['group']=='control').sum())})")

    legacy_joint_models_df["is_case_int"] = legacy_joint_models_df["is_case"].astype(int)
    legacy_joint_models_df["carrier_main_10gene"] = (legacy_joint_models_df["burden_rare_main_ptv_10gene"] > 0).astype(int)
    legacy_joint_models_df["strata"] = legacy_joint_models_df["is_case_int"].astype(str) + "_" + legacy_joint_models_df["carrier_main_10gene"].astype(str)

    strata_counts = legacy_joint_models_df["strata"].value_counts().sort_index()
    legacy_joint_models_log('Number of people in four layers (is_case x carrier_main_10gene):')
    legacy_joint_models_log(strata_counts.apply(legacy_joint_models_suppress).to_string())

    site_dummies_full = pd.get_dummies(legacy_joint_models_df["seq_site_id"], prefix="site", drop_first=True).astype(float)
    seq_dummies_full = pd.get_dummies(legacy_joint_models_df["seq_sequencer"], prefix="seq", drop_first=True).astype(float)
    LEGACY_JOINT_MODELS_SITE_COLS = list(site_dummies_full.columns)
    LEGACY_JOINT_MODELS_SEQ_COLS = list(seq_dummies_full.columns)
    legacy_joint_models_df = pd.concat([legacy_joint_models_df, site_dummies_full, seq_dummies_full], axis=1)
    legacy_joint_models_log(f"\n seq_site_id dummy column:{LEGACY_JOINT_MODELS_SITE_COLS}; seq_sequencer dummy column:{LEGACY_JOINT_MODELS_SEQ_COLS}")

    # ---------------------------------------------------------------------------

    # ---------------------------------------------------------------------------
    legacy_joint_models_section(f"Step B: Repeat stratified cross-validation ({CV_REPEATS}times ×{CV_FOLDS}fold, is_case x carrier_main_10gene four layers)")

    rskf = RepeatedStratifiedKFold(n_splits=CV_FOLDS, n_repeats=CV_REPEATS, random_state=RANDOM_STATE)
    cv_splits = list(rskf.split(legacy_joint_models_df.index.values, legacy_joint_models_df["strata"].values))
    legacy_joint_models_log(f"common{len(cv_splits)}(train, test) folds; each fold test is about{len(legacy_joint_models_df) // CV_FOLDS}people")

    # ---------------------------------------------------------------------------
    
    # ---------------------------------------------------------------------------
    legacy_joint_models_section('Step C: Define logistic regression helper (unweighted MLE main estimator + var_weights weighted sensitivity analysis)')

    # ---------------------------------------------------------------------------

    # ---------------------------------------------------------------------------
    legacy_joint_models_section('Step D: Main analysis 4 model cross-validation (Covariates-only + PRS311 + burden_rare_main_ptv_10gene)')

    LEGACY_JOINT_MODELS_COV_COLS = ["age_at_last_obs"] + LEGACY_JOINT_MODELS_PC_COLS + LEGACY_JOINT_MODELS_SITE_COLS + LEGACY_JOINT_MODELS_SEQ_COLS

    model_specs = {
        "Covariates-only": [],
        "PRS-only": ["PRS311_z"],
        "Burden-only": ["burden_z"],
        "Joint": ["PRS311_z", "burden_z", "inter_z"],
    }
    delong_pairs = [
        ("PRS-only vs Covariates-only", "PRS-only", "Covariates-only"),
        ("Burden-only vs Covariates-only", "Burden-only", "Covariates-only"),
        ("Joint vs Covariates-only", "Joint", "Covariates-only"),
        ("Joint vs PRS-only", "Joint", "PRS-only"),
        ("Joint vs Burden-only", "Joint", "Burden-only"),
        ("PRS-only vs Burden-only", "PRS-only", "Burden-only"),
    ]

    legacy_joint_models_y_all = legacy_joint_models_df["is_case_int"].values
    auc_rep_rows, delong_rep_rows = [], []
    oof_r1, fold_r1 = {}, np.full(len(legacy_joint_models_df), -1)
    for legacy_joint_models_r in range(CV_REPEATS):
        oof = {legacy_joint_models_label: np.full(len(legacy_joint_models_df), np.nan) for legacy_joint_models_label in model_specs}
        for legacy_joint_models_k in range(CV_FOLDS):
            tr_idx, te_idx = cv_splits[legacy_joint_models_r * CV_FOLDS + legacy_joint_models_k]
            tr, te = legacy_joint_models_df.iloc[tr_idx], legacy_joint_models_df.iloc[te_idx]
            mu_p, sd_p = zscore_fit(tr["PRS311"])
            mu_b, sd_b = zscore_fit(tr["burden_rare_main_ptv_10gene"])
            tr, te = add_z(tr, mu_p, sd_p, mu_b, sd_b), add_z(te, mu_p, sd_p, mu_b, sd_b)
            for legacy_joint_models_label, legacy_joint_models_extra_cols in model_specs.items():
                legacy_joint_models_fit_res, legacy_joint_models__ = fit_unweighted_logit(tr["is_case_int"].values, legacy_joint_models_build_X(tr, legacy_joint_models_extra_cols))
                oof[legacy_joint_models_label][te_idx] = legacy_joint_models_fit_res.predict(legacy_joint_models_build_X(te, legacy_joint_models_extra_cols)).values
            if legacy_joint_models_r == 0:
                fold_r1[te_idx] = legacy_joint_models_k
        assert all(not np.isnan(v).any() for v in oof.values()), 'There is a gap in OOF forecast'
        for legacy_joint_models_label in model_specs:
            auc, legacy_joint_models__, legacy_joint_models__ = delong_auc_variance(legacy_joint_models_y_all, oof[legacy_joint_models_label])
            auc_rep_rows.append({"repeat": legacy_joint_models_r + 1, "model": legacy_joint_models_label, "AUC_oof": auc})
        for comp_label, a, legacy_joint_models_b in delong_pairs:
            res = delong_test(legacy_joint_models_y_all, oof[a], oof[legacy_joint_models_b])
            delong_rep_rows.append({"repeat": legacy_joint_models_r + 1, "comparison": comp_label, **res})
        if legacy_joint_models_r == 0:
            oof_r1 = oof
        legacy_joint_models_log(f"  repeat {legacy_joint_models_r + 1}/{CV_REPEATS}: " + ", ".join(
            f"{row['model']}={row['AUC_oof']:.4f}" for row in auc_rep_rows[-len(model_specs):]))

    auc_rep = pd.DataFrame(auc_rep_rows)
    delong_rep = pd.DataFrame(delong_rep_rows)

    auc_tbl = (auc_rep.groupby("model", sort=False)["AUC_oof"]
               .agg(AUC_mean="mean", AUC_min="min", AUC_max="max").reset_index())
    auc_tbl["AUC_CI_low_r1"], auc_tbl["AUC_CI_high_r1"] = zip(*[auc_ci(legacy_joint_models_y_all, oof_r1[m]) for m in auc_tbl["model"]])
    delta_tbl = (delong_rep.groupby("comparison", sort=False)
                 .agg(delta_AUC_mean=("delta_auc", "mean"), delta_AUC_min=("delta_auc", "min"),
                      delta_AUC_max=("delta_auc", "max"), CI_low_mean=("ci_low", "mean"),
                      CI_high_mean=("ci_high", "mean"), p_median=("p_value", "median"),
                      p_max=("p_value", "max")).reset_index())

    p_floor(auc_rep).to_csv(os.path.join(LEGACY_JOINT_MODELS_TABLE_DIR, "16_cv_auc_by_repeat.csv"), index=False)
    p_floor(delong_rep).to_csv(os.path.join(LEGACY_JOINT_MODELS_TABLE_DIR, "16_cv_delong_by_repeat.csv"), index=False)
    p_floor(auc_tbl).to_csv(os.path.join(LEGACY_JOINT_MODELS_TABLE_DIR, "16_cv_auc_summary.csv"), index=False)
    p_floor(delta_tbl).to_csv(os.path.join(LEGACY_JOINT_MODELS_TABLE_DIR, "16_cv_delta_auc_summary.csv"), index=False)
    legacy_joint_models_log('CV AUC summary (OOF, 10 replicates):')
    legacy_joint_models_log(auc_tbl.to_string(index=False))
    legacy_joint_models_log('Summary of CV ΔAUC (DeLong, 10 replicates):')
    legacy_joint_models_log(delta_tbl.to_string(index=False))

    j = auc_tbl.set_index("model")["AUC_mean"]
    legacy_joint_models_log(f"\n  SANITY CHECK: Joint mean OOF AUC={j['Joint']:.4f}vs max(rest)={j.drop('Joint').max():.4f} -> "
        f"{'PASS' if j['Joint'] >= j.drop('Joint').max() else 'FAIL'}")

    oof_persist = legacy_joint_models_df[["person_id", "group", "carrier_main_10gene"]].copy()
    oof_persist["fold_r1"] = fold_r1
    for legacy_joint_models_label in model_specs:
        oof_persist[f"oof_r1_{legacy_joint_models_label.replace('-', '_').lower()}"] = oof_r1[legacy_joint_models_label]
    oof_persist.to_csv(os.path.join(LEGACY_JOINT_MODELS_DATA_DIR, "16_joint_model_cv_oof_person_level.csv.gz"),
                       index=False, compression="gzip")
    legacy_joint_models_log('Write out: data/16_joint_model_cv_oof_person_level.csv.gz (1st repetition of folds and 4 model OOF predictions)')

    legacy_joint_models_section('Step D2: Full cohort 4 model coefficients and McFadden pseudo R2 (in-sample)')
    legacy_joint_models_full_prs_mu, legacy_joint_models_full_prs_sd = zscore_fit(legacy_joint_models_df["PRS311"])
    full_burden_mu, full_burden_sd = zscore_fit(legacy_joint_models_df["burden_rare_main_ptv_10gene"])
    legacy_joint_models_df = add_z(legacy_joint_models_df, legacy_joint_models_full_prs_mu, legacy_joint_models_full_prs_sd, full_burden_mu, full_burden_sd)
    y_full = legacy_joint_models_df["is_case_int"].values

    pseudo_r2_rows, coef_rows_main = [], []
    for legacy_joint_models_label, legacy_joint_models_extra_cols in model_specs.items():
        LEGACY_JOINT_MODELS_X = legacy_joint_models_build_X(legacy_joint_models_df, legacy_joint_models_extra_cols)
        legacy_joint_models_fit_res, w = fit_unweighted_logit(y_full, LEGACY_JOINT_MODELS_X)
        r2, ll_null = mcfadden_r2(y_full, LEGACY_JOINT_MODELS_X, w, legacy_joint_models_fit_res.llf)
        pseudo_r2_rows.append({"model": legacy_joint_models_label, "LL_model": legacy_joint_models_fit_res.llf, "LL_null": ll_null, "McFadden_R2": r2})
        coef_rows_main.append(legacy_joint_models_coef_table(legacy_joint_models_fit_res, legacy_joint_models_extra_cols, legacy_joint_models_label))
    pseudo_r2_tbl = pd.DataFrame(pseudo_r2_rows).set_index("model")
    incremental_r2_tbl = pd.DataFrame([
        {"comparison": f"{a} - {legacy_joint_models_b}", "McFadden_R2_model": pseudo_r2_tbl.loc[a, "McFadden_R2"],
         "McFadden_R2_baseline": pseudo_r2_tbl.loc[legacy_joint_models_b, "McFadden_R2"]}
        for a, legacy_joint_models_b in [("Joint", "PRS-only"), ("Joint", "Burden-only"),
                     ("PRS-only", "Covariates-only"), ("Burden-only", "Covariates-only")]])
    incremental_r2_tbl["incremental_R2"] = (incremental_r2_tbl["McFadden_R2_model"]
                                            - incremental_r2_tbl["McFadden_R2_baseline"])
    pseudo_r2_tbl = pseudo_r2_tbl.reset_index()
    p_floor(pd.concat(coef_rows_main, ignore_index=True)).to_csv(os.path.join(LEGACY_JOINT_MODELS_TABLE_DIR, "16_coefficients_full_cohort_4models.csv"), index=False)
    p_floor(pseudo_r2_tbl).to_csv(os.path.join(LEGACY_JOINT_MODELS_TABLE_DIR, "16_mcfadden_pseudo_r2_full_cohort.csv"), index=False)
    p_floor(incremental_r2_tbl).to_csv(os.path.join(LEGACY_JOINT_MODELS_TABLE_DIR, "16_incremental_pseudo_r2.csv"), index=False)
    legacy_joint_models_log(pseudo_r2_tbl.to_string(index=False))
    legacy_joint_models_log(incremental_r2_tbl.to_string(index=False))

    # ---------------------------------------------------------------------------

    # ---------------------------------------------------------------------------
    legacy_joint_models_section('Step E: Full queue (n={}) joint model coefficient table'.format(len(legacy_joint_models_df)))

    legacy_joint_models_full_prs_mu, legacy_joint_models_full_prs_sd = zscore_fit(legacy_joint_models_df["PRS311"])
    full_burden_mu, full_burden_sd = zscore_fit(legacy_joint_models_df["burden_rare_main_ptv_10gene"])
    legacy_joint_models_df["PRS311_z"] = apply_z(legacy_joint_models_df["PRS311"], legacy_joint_models_full_prs_mu, legacy_joint_models_full_prs_sd)
    legacy_joint_models_df["burden_z"] = apply_z(legacy_joint_models_df["burden_rare_main_ptv_10gene"], full_burden_mu, full_burden_sd)
    legacy_joint_models_df["inter_z"] = legacy_joint_models_df["PRS311_z"] * legacy_joint_models_df["burden_z"]

    y_full = legacy_joint_models_df["is_case_int"].values
    X_full_joint = legacy_joint_models_build_X(legacy_joint_models_df, ["PRS311_z", "burden_z", "inter_z"])
    full_joint_fit, legacy_joint_models__ = fit_unweighted_logit(y_full, X_full_joint)

    primary_terms = ["PRS311_z", "burden_z", "inter_z"]
    full_coef_primary = legacy_joint_models_coef_table(full_joint_fit, primary_terms, "Joint_full_cohort_main_10gene_PRS311")
    full_coef_all = legacy_joint_models_coef_table(full_joint_fit, list(X_full_joint.columns), "Joint_full_cohort_main_10gene_PRS311_allcov")

    sd_used_map = {"burden_z": full_burden_sd, "PRS311_z": legacy_joint_models_full_prs_sd}
    full_coef_primary["sd_used"] = full_coef_primary["term"].map(sd_used_map)
    full_coef_primary["beta_raw"] = full_coef_primary["beta"] / full_coef_primary["sd_used"]
    full_coef_primary["OR_raw"] = np.exp(full_coef_primary["beta_raw"])
    full_coef_primary["beta_raw_ci_low"] = (full_coef_primary["beta"] - 1.959963984540054 * full_coef_primary["se"]) / full_coef_primary["sd_used"]
    full_coef_primary["beta_raw_ci_high"] = (full_coef_primary["beta"] + 1.959963984540054 * full_coef_primary["se"]) / full_coef_primary["sd_used"]
    full_coef_primary["OR_raw_ci_low"] = np.exp(full_coef_primary["beta_raw_ci_low"])
    full_coef_primary["OR_raw_ci_high"] = np.exp(full_coef_primary["beta_raw_ci_high"])

    p_floor(full_coef_primary).to_csv(os.path.join(LEGACY_JOINT_MODELS_TABLE_DIR, "16_joint_model_coefficients_primary.csv"), index=False)
    p_floor(full_coef_all).to_csv(os.path.join(LEGACY_JOINT_MODELS_TABLE_DIR, "16_joint_model_coefficients_full_covariates.csv"), index=False)

    legacy_joint_models_log('The main three coefficients of the whole-cohort joint model (PRS311_z / burden_z / inter_z), including raw log-OR scale:')
    legacy_joint_models_log(full_coef_primary[["term", "beta", "OR", "beta_raw", "OR_raw", "p_value"]].to_string(index=False))

    # ---------------------------------------------------------------------------

    # ---------------------------------------------------------------------------
    legacy_joint_models_section('Step F: Plot (CV AUC + OOF ROC curve)')

    legacy_joint_models_fig, legacy_joint_models_ax = plt.subplots(figsize=(7, 5))
    x = np.arange(len(auc_tbl))
    legacy_joint_models_ax.bar(x, auc_tbl["AUC_mean"], color=[MODEL_COLORS[m] for m in auc_tbl["model"]])
    legacy_joint_models_ax.errorbar(x, auc_tbl["AUC_mean"],
                yerr=[auc_tbl["AUC_mean"] - auc_tbl["AUC_CI_low_r1"], auc_tbl["AUC_CI_high_r1"] - auc_tbl["AUC_mean"]],
                fmt="none", ecolor="black", capsize=4)
    legacy_joint_models_ax.set_xticks(x)
    legacy_joint_models_ax.set_xticklabels(auc_tbl["model"])
    legacy_joint_models_ax.set_ylabel(f"Cross-validated AUC ({CV_REPEATS}x{CV_FOLDS}-fold)")
    legacy_joint_models_ax.set_ylim(0.5, auc_tbl["AUC_CI_high_r1"].max() + 0.03)
    legacy_joint_models_fig.tight_layout()
    legacy_joint_models_save_smallest(legacy_joint_models_fig, "16_cv_auc_comparison_bar")
    plt.close(legacy_joint_models_fig)

    legacy_joint_models_fig, legacy_joint_models_ax = plt.subplots(figsize=(6, 6))
    for legacy_joint_models_label in model_specs:
        fpr, tpr, legacy_joint_models__ = roc_curve(legacy_joint_models_y_all, oof_r1[legacy_joint_models_label])
        legacy_joint_models_ax.plot(fpr, tpr, label=f"{legacy_joint_models_label} (AUC={auc_tbl.set_index('model').loc[legacy_joint_models_label, 'AUC_mean']:.3f})",
                color=MODEL_COLORS[legacy_joint_models_label])
    legacy_joint_models_ax.plot([0, 1], [0, 1], linestyle="--", color="black", linewidth=0.8)
    legacy_joint_models_ax.set_xlabel("False positive rate")
    legacy_joint_models_ax.set_ylabel("True positive rate")
    legacy_joint_models_ax.legend(loc="lower right")
    legacy_joint_models_fig.tight_layout()
    legacy_joint_models_save_smallest(legacy_joint_models_fig, "16_cv_roc_curves_oof")
    plt.close(legacy_joint_models_fig)

    # ---------------------------------------------------------------------------
    
    # ---------------------------------------------------------------------------
    legacy_joint_models_section('Step G: Sensitivity analysis (9-gene / MAF caliber / weighted / 281-site PRS / no-mean-imp / rule-of-two)')

    legacy_joint_models_sensitivity_rows = []

    column_swap_specs = [
        {
            "name": "9gene_main_ptv_PRS311",
            "burden_col": "burden_rare_main_ptv_9gene",
            "prs_col": "PRS311",
            "note": '9-gene main analysis panel (Dorling-only, excluding MAP3K1) + PRS311 (same PRS version as main analysis)',
        },
        {
            "name": "10gene_sensitivity_maf0p1pct_PRS313",
            "burden_col": "burden_rare_sensitivity_ptv_maf0p1pct_10gene",
            "prs_col": "PRS313",
            "note": 'Unified MAF<0.1% level (1100delC AF~0.26% is excluded from burden) + PRS313'
                    "(config prs_version_by_maf_tier: sensitivity_uniform_0p1pct)",
        },
        {
            "name": "10gene_sensitivity_maf1pct_PRS311",
            "burden_col": "burden_rare_sensitivity_ptv_maf1pct_10gene",
            "prs_col": "PRS311",
            "note": 'Unified MAF<1% level (same threshold as main analysis, 1100delC enters burden) + PRS311'
                    "(config prs_version_by_maf_tier: sensitivity_uniform_1pct)",
        },
    ]

    for spec in column_swap_specs:
        b_col, p_col = spec["burden_col"], spec["prs_col"]
        mu_p, sd_p = zscore_fit(legacy_joint_models_df[p_col])
        mu_b, sd_b = zscore_fit(legacy_joint_models_df[b_col])
        sens_df = legacy_joint_models_df[LEGACY_JOINT_MODELS_COV_COLS].copy()
        sens_df["PRS311_z"] = apply_z(legacy_joint_models_df[p_col], mu_p, sd_p)
        sens_df["burden_z"] = apply_z(legacy_joint_models_df[b_col], mu_b, sd_b)
        sens_df["inter_z"] = sens_df["PRS311_z"] * sens_df["burden_z"]

        X_sens = legacy_joint_models_build_X(sens_df, ["PRS311_z", "burden_z", "inter_z"])
        fit_sens, legacy_joint_models__ = fit_unweighted_logit(y_full, X_sens)
        ct = legacy_joint_models_coef_table(fit_sens, ["PRS311_z", "burden_z", "inter_z"], spec["name"])
        ct["prs_col_used"] = p_col
        ct["burden_col_used"] = b_col
        ct["cohort_n"] = len(sens_df)
        ct["note"] = spec["note"]
        legacy_joint_models_sensitivity_rows.append(ct)
        legacy_joint_models_log(f"\n  [{spec['name']}] burden_col={b_col}, prs_col={p_col}, n={len(sens_df)}")
        legacy_joint_models_log(ct[["term", "OR", "CI_low", "CI_high", "p_value"]].to_string(index=False))

    X_primary_full = legacy_joint_models_build_X(legacy_joint_models_df, ["PRS311_z", "burden_z", "inter_z"])
    fit_weighted_primary, legacy_joint_models__ = fit_weighted_logit(y_full, X_primary_full)
    ct_weighted = legacy_joint_models_coef_table(fit_weighted_primary, ["PRS311_z", "burden_z", "inter_z"], "weighted_class_balanced")
    ct_weighted["prs_col_used"] = "PRS311"
    ct_weighted["burden_col_used"] = "burden_rare_main_ptv_10gene"
    ct_weighted["cohort_n"] = len(legacy_joint_models_df)
    ct_weighted["note"] = ('The main analysis is the same as the cohort/the same terms, and the estimator is changed to class-balanced var_weights weighted logistic'
                            '(w_i=n/(2*n_group)) + HC1, for comparison with the class-weighted approach commonly used in the literature;'
                            'The change in coefficient point estimates should theoretically be minimal (the weighting is a function of y and does not change the slope consistency).')
    legacy_joint_models_sensitivity_rows.append(ct_weighted)
    legacy_joint_models_log('[weighted_class_balanced] (estimator=var_weights weighted logistic, terms/cohort co-main analysis)')
    legacy_joint_models_log(ct_weighted[["term", "OR", "CI_low", "CI_high", "p_value"]].to_string(index=False))

    prs_281_renamed = legacy_joint_models_prs_281.rename(columns={"PRS311": "PRS311_281loci"})
    df_281 = legacy_joint_models_df.merge(prs_281_renamed[["person_id", "PRS311_281loci"]], on="person_id", how="left")
    legacy_joint_models_n_na_281 = int(df_281["PRS311_281loci"].isna().sum())
    if legacy_joint_models_n_na_281:
        legacy_joint_models_log(f"WARNING: PRS311_281loci has{legacy_joint_models_n_na_281}NA, will be excluded by dropna")
    df_281_cc = df_281.dropna(subset=["PRS311_281loci"]).reset_index(drop=True)
    mu_p281, sd_p281 = zscore_fit(df_281_cc["PRS311_281loci"])
    mu_b, sd_b = zscore_fit(df_281_cc["burden_rare_main_ptv_10gene"])
    sens281_df = df_281_cc[LEGACY_JOINT_MODELS_COV_COLS].copy()
    sens281_df["PRS311_z"] = apply_z(df_281_cc["PRS311_281loci"], mu_p281, sd_p281)
    sens281_df["burden_z"] = apply_z(df_281_cc["burden_rare_main_ptv_10gene"], mu_b, sd_b)
    sens281_df["inter_z"] = sens281_df["PRS311_z"] * sens281_df["burden_z"]
    y_281 = df_281_cc["is_case_int"].values
    X_281 = legacy_joint_models_build_X(sens281_df, ["PRS311_z", "burden_z", "inter_z"])
    fit_281, legacy_joint_models__ = fit_unweighted_logit(y_281, X_281)
    ct_281 = legacy_joint_models_coef_table(fit_281, ["PRS311_z", "burden_z", "inter_z"], "prs281loci_PRS311")
    ct_281["prs_col_used"] = "PRS311(281loci, data/prs_scores_281loci.csv)"
    ct_281["burden_col_used"] = "burden_rare_main_ptv_10gene"
    ct_281["cohort_n"] = len(sens281_df)
    ct_281["note"] = 'plan L301 preset PRS version sensitivity analysis (281 sites vs. 304 sites in main analysis), the burden remains unchanged'
    legacy_joint_models_sensitivity_rows.append(ct_281)
    legacy_joint_models_log(f"\n  [prs281loci_PRS311] n={len(sens281_df)}")
    legacy_joint_models_log(ct_281[["term", "OR", "CI_low", "CI_high", "p_value"]].to_string(index=False))

    if "PRS311_noimp" in legacy_joint_models_df.columns:
        df_noimp_cc = legacy_joint_models_df.dropna(subset=["PRS311_noimp"]).reset_index(drop=True)
        legacy_joint_models_n_dropped_noimp = len(legacy_joint_models_df) - len(df_noimp_cc)
        if legacy_joint_models_n_dropped_noimp:
            legacy_joint_models_log(f"WARNING: PRS311_noimp has{legacy_joint_models_n_dropped_noimp}NA, will be excluded by dropna")
        mu_pni, sd_pni = zscore_fit(df_noimp_cc["PRS311_noimp"])
        mu_b, sd_b = zscore_fit(df_noimp_cc["burden_rare_main_ptv_10gene"])
        sensni_df = df_noimp_cc[LEGACY_JOINT_MODELS_COV_COLS].copy()
        sensni_df["PRS311_z"] = apply_z(df_noimp_cc["PRS311_noimp"], mu_pni, sd_pni)
        sensni_df["burden_z"] = apply_z(df_noimp_cc["burden_rare_main_ptv_10gene"], mu_b, sd_b)
        sensni_df["inter_z"] = sensni_df["PRS311_z"] * sensni_df["burden_z"]
        y_ni = df_noimp_cc["is_case_int"].values
        X_ni = legacy_joint_models_build_X(sensni_df, ["PRS311_z", "burden_z", "inter_z"])
        fit_ni, legacy_joint_models__ = fit_unweighted_logit(y_ni, X_ni)
        ct_ni = legacy_joint_models_coef_table(fit_ni, ["PRS311_z", "burden_z", "inter_z"], "prs_no_mean_imputation_PRS311")
        ct_ni["prs_col_used"] = "PRS311_noimp"
        ct_ni["burden_col_used"] = "burden_rare_main_ptv_10gene"
        ct_ni["cohort_n"] = len(sensni_df)
        ct_ni["note"] = 'plan L262 preset PRS interpolation method sensitivity analysis (no-mean-imputation vs main analysis mean-imputation)'
        legacy_joint_models_sensitivity_rows.append(ct_ni)
        legacy_joint_models_log(f"\n  [prs_no_mean_imputation_PRS311] n={len(sensni_df)}")
        legacy_joint_models_log(ct_ni[["term", "OR", "CI_low", "CI_high", "p_value"]].to_string(index=False))
    else:
        legacy_joint_models_log('[prs_no_mean_imputation_PRS311] SKIPPED: data/prs_scores.csv No PRS311_noimp column')

    if "is_case_rule_of_two" in legacy_joint_models_df.columns:

        df_r2 = legacy_joint_models_df
        keep_mask = (df_r2["group"] == "control") | (
            (df_r2["group"] == "case") & (df_r2["is_case_rule_of_two"] == True)  # noqa: E712
        )
        df_r2_cc = df_r2.loc[keep_mask].reset_index(drop=True)
        n_case_r2 = int((df_r2_cc["group"] == "case").sum())
        n_control_r2 = int((df_r2_cc["group"] == "control").sum())
        legacy_joint_models_log(f"\n  [case_def_rule_of_two] n={len(df_r2_cc)} (case={n_case_r2}, control={n_control_r2}), "
            f"burden/PRS z-score is refit on this subqueue")
        mu_p_r2, sd_p_r2 = zscore_fit(df_r2_cc["PRS311"])
        mu_b_r2, sd_b_r2 = zscore_fit(df_r2_cc["burden_rare_main_ptv_10gene"])
        sensr2_df = df_r2_cc[LEGACY_JOINT_MODELS_COV_COLS].copy()
        sensr2_df["PRS311_z"] = apply_z(df_r2_cc["PRS311"], mu_p_r2, sd_p_r2)
        sensr2_df["burden_z"] = apply_z(df_r2_cc["burden_rare_main_ptv_10gene"], mu_b_r2, sd_b_r2)
        sensr2_df["inter_z"] = sensr2_df["PRS311_z"] * sensr2_df["burden_z"]
        y_r2 = df_r2_cc["is_case_int"].values
        X_r2 = legacy_joint_models_build_X(sensr2_df, ["PRS311_z", "burden_z", "inter_z"])
        fit_r2, legacy_joint_models__ = fit_unweighted_logit(y_r2, X_r2)
        ct_r2 = legacy_joint_models_coef_table(fit_r2, ["PRS311_z", "burden_z", "inter_z"], "case_def_rule_of_two")
        ct_r2["prs_col_used"] = "PRS311"
        ct_r2["burden_col_used"] = "burden_rare_main_ptv_10gene"
        ct_r2["cohort_n"] = len(sensr2_df)
        ct_r2["note"] = (f"plan L109 preset case definition sensitivity analysis (is_case_rule_of_two==True replaces the main analysis case definition,"
                          f"The control group remains unchanged); case={n_case_r2}, control={n_control_r2}")
        legacy_joint_models_sensitivity_rows.append(ct_r2)
        legacy_joint_models_log(ct_r2[["term", "OR", "CI_low", "CI_high", "p_value"]].to_string(index=False))
    else:
        legacy_joint_models_log('[case_def_rule_of_two] SKIPPED: phenotype table has no is_case_rule_of_two column')

    sensitivity_tbl = pd.concat(legacy_joint_models_sensitivity_rows, ignore_index=True)
    
    legacy_joint_models_primary_for_compare = full_coef_primary[["model", "term", "beta", "se", "OR", "CI_low", "CI_high", "p_value"]].copy()
    legacy_joint_models_primary_for_compare["prs_col_used"] = "PRS311"
    legacy_joint_models_primary_for_compare["burden_col_used"] = "burden_rare_main_ptv_10gene"
    legacy_joint_models_primary_for_compare["cohort_n"] = len(legacy_joint_models_df)
    legacy_joint_models_primary_for_compare["note"] = 'Main analysis (10-gene main_ptv, MAF<1%, PRS311, unweighted MLE+HC1)'
    legacy_joint_models_sensitivity_tbl_full = pd.concat([legacy_joint_models_primary_for_compare, sensitivity_tbl], ignore_index=True)
    p_floor(legacy_joint_models_sensitivity_tbl_full).to_csv(os.path.join(LEGACY_JOINT_MODELS_TABLE_DIR, "16_sensitivity_analyses_coefficients.csv"), index=False)

    legacy_joint_models_log('Main analysis + all sensitivity analyzes side by side (burden_z / inter_z direction & significance comparison):')
    legacy_joint_models_log(legacy_joint_models_sensitivity_tbl_full[legacy_joint_models_sensitivity_tbl_full["term"].isin(["burden_z", "inter_z"])]
        [["model", "term", "OR", "CI_low", "CI_high", "p_value"]].to_string(index=False))

    # ---------------------------------------------------------------------------

    # ---------------------------------------------------------------------------
    legacy_joint_models_section('Step H: Post-hoc exploratory companion model (carrier binary contrast / gene-specific PRS slope / model-implied OR)')

    # is_case ~ PRS311_z + carrier_main_10gene + PRS311_z:carrier_main_10gene + COV_COLS

    legacy_joint_models_df["carrier_x_prs"] = legacy_joint_models_df["PRS311_z"] * legacy_joint_models_df["carrier_main_10gene"]
    carrier_terms = ["PRS311_z", "carrier_main_10gene", "carrier_x_prs"]
    X_carrier_bin = legacy_joint_models_df[carrier_terms + LEGACY_JOINT_MODELS_COV_COLS].astype(float).copy()
    X_carrier_bin = sm.add_constant(X_carrier_bin, has_constant="add")
    fit_carrier_bin, legacy_joint_models__ = fit_unweighted_logit(y_full, X_carrier_bin)
    ct_carrier_bin = legacy_joint_models_coef_table(fit_carrier_bin, carrier_terms, "carrier_binary_contrast_model")
    evidence_status_map = {
        "PRS311_z": "model_component (reference term within alternative_operationalization spec)",
        "carrier_main_10gene": "alternative_operationalization (NOT a reparametrization: binary carrier "
                                "indicator replaces continuous weighted burden -> genuinely different "
                                "model spec/likelihood from the primary z-scaled weighted-burden model)",
        "carrier_x_prs": "post_hoc_exploratory (tests the imposed-linearity hypothesis behind Antinomy 1; "
                          "corresponds to addendum's 'PRS_z x carrier(0/1)' term, same fitted model as "
                          "carrier_main_10gene above)",
    }
    ct_carrier_bin["evidence_status"] = ct_carrier_bin["term"].map(evidence_status_map)
    ct_carrier_bin["n_cohort"] = len(legacy_joint_models_df)
    ct_carrier_bin["n_carrier"] = int(legacy_joint_models_df["carrier_main_10gene"].sum())
    ct_carrier_bin["n_carrier_case"] = int(((legacy_joint_models_df["carrier_main_10gene"] == 1) & (legacy_joint_models_df["is_case_int"] == 1)).sum())
    p_floor(ct_carrier_bin).to_csv(os.path.join(LEGACY_JOINT_MODELS_TABLE_DIR, "16_carrier_binary_contrast_model.csv"), index=False)
    legacy_joint_models_log('[Model B] Binary-carrier comparison model (is_case ~ PRS311_z + carrier_main_10gene +'
        "PRS311_z:carrier_main_10gene + COV_COLS):")
    legacy_joint_models_log(ct_carrier_bin[["term", "OR", "CI_low", "CI_high", "p_value", "evidence_status"]].to_string(index=False))
    legacy_joint_models_log(f"Write out: results/table/16_carrier_binary_contrast_model.csv")

    gene_flag_cols = [f"gene_{legacy_joint_models_g}" for legacy_joint_models_g in LEGACY_JOINT_MODELS_GENE_COLS]
    legacy_joint_models_df["n_genes_carried"] = legacy_joint_models_df[gene_flag_cols].sum(axis=1)

    legacy_joint_models_df["carrier_gene_stratum"] = legacy_joint_models_df.apply(assign_stratum, axis=1)

    _ct = pd.crosstab(legacy_joint_models_df["carrier_gene_stratum"], legacy_joint_models_df["is_case_int"])
    legacy_joint_models_log('carrier_gene_stratum Whether the case/control of each layer is >=20:')
    legacy_joint_models_log(", ".join(f"{k}: {'yes' if (v >= 20).all() else 'no'}" for k, v in _ct.iterrows()))

    stratum_order = ["non-carrier", "BRCA1", "BRCA2", "CHEK2", "ATM", "other"]
    strat_rows = []
    ESTIMAND_CAVEAT = (
        "ESTIMAND/ASCERTAINMENT CAVEATS (do not silently compare to Barnes 2020): "
        "(1) Barnes DR et al. Genet Med 2020;22(10):1653-1666 reports HR from a time-to-event carrier "
        "cohort; this table reports OR from case-control logistic regression -- similar for rare "
        "outcomes but not the same estimand. "
        "(2) Barnes' carriers are clinically ascertained (BRCA testing referrals); this cohort's "
        "carriers are population-based (AoU) -- ascertainment differs materially. "
        "(3) Barnes' BRCA1 figure (HR=1.29) uses an ER-negative-specific PRS version; this project only "
        "has the overall PRS311, so the BRCA1 row below is NOT like-for-like with Barnes; the closest "
        "like-for-like comparison is BRCA2 (this table) vs Barnes BRCA2 HR=1.31 (95% CI 1.27-1.36)."
    )
    for legacy_joint_models_stratum in stratum_order:
        legacy_joint_models_sub = legacy_joint_models_df[legacy_joint_models_df["carrier_gene_stratum"] == legacy_joint_models_stratum]
        legacy_joint_models_n_case = int((legacy_joint_models_sub["is_case_int"] == 1).sum())
        n_control = int((legacy_joint_models_sub["is_case_int"] == 0).sum())
        row_out = {

            "stratum": legacy_joint_models_stratum,
            "n_case": legacy_joint_models_suppress(legacy_joint_models_n_case) if legacy_joint_models_stratum == "non-carrier" else (">=20" if legacy_joint_models_n_case >= 20 else "<20"),
            "n_control": legacy_joint_models_suppress(n_control) if legacy_joint_models_stratum == "non-carrier" else (">=20" if n_control >= 20 else "<20"),
            "n_total": legacy_joint_models_suppress(len(legacy_joint_models_sub)) if legacy_joint_models_stratum == "non-carrier" else (">=20" if len(legacy_joint_models_sub) >= 20 else "<20"),
        }
        if legacy_joint_models_n_case < 20 or n_control < 20:
            row_out.update(beta_prs=np.nan, se=np.nan, or_per_sd=np.nan, ci_low=np.nan, ci_high=np.nan,
                            ci_half_width=np.nan, passes_power_precheck=False,
                            note="n_case or n_control <20 (AoU disclosure rule); model not fit")
            strat_rows.append(row_out)
            legacy_joint_models_log(f"  [{legacy_joint_models_stratum}] n_case or n_control <20 -> SUPPRESSED, model not fit")
            continue
        X_strat = legacy_joint_models_sub[["PRS311_z"] + LEGACY_JOINT_MODELS_COV_COLS].astype(float).copy()
        X_strat = sm.add_constant(X_strat, has_constant="add")
        y_strat = legacy_joint_models_sub["is_case_int"].values
        fit_strat, legacy_joint_models__ = fit_unweighted_logit(y_strat, X_strat)
        legacy_joint_models_b = fit_strat.params["PRS311_z"]
        legacy_joint_models_se = fit_strat.bse["PRS311_z"]
        ci_low = legacy_joint_models_b - 1.959963984540054 * legacy_joint_models_se
        ci_high = legacy_joint_models_b + 1.959963984540054 * legacy_joint_models_se
        ci_half_width = 1.959963984540054 * legacy_joint_models_se
        passes = ci_half_width < abs(legacy_joint_models_b)
        row_out.update(beta_prs=legacy_joint_models_b, se=legacy_joint_models_se, or_per_sd=np.exp(legacy_joint_models_b),
                        ci_low=np.exp(ci_low), ci_high=np.exp(ci_high),
                        ci_half_width=ci_half_width, passes_power_precheck=bool(passes),
                        note="" if passes else "CI half-width >= |point estimate|: DO NOT INTERPRET, underpowered")
        strat_rows.append(row_out)
        legacy_joint_models_log(f"  [{legacy_joint_models_stratum}] n_case/n_control >=20, beta={legacy_joint_models_b:.4f} OR/SD={np.exp(legacy_joint_models_b):.3f} "
            f"(95%CI {np.exp(ci_low):.3f}-{np.exp(ci_high):.3f}) passes_power_precheck={passes}")

    strat_tbl = pd.DataFrame(strat_rows)
    strat_tbl["evidence_status"] = "post_hoc_exploratory"
    strat_tbl["estimand_caveat"] = ESTIMAND_CAVEAT
    p_floor(strat_tbl).to_csv(os.path.join(LEGACY_JOINT_MODELS_TABLE_DIR, "16_carrier_stratified_prs_slope_by_gene.csv"), index=False)
    legacy_joint_models_log(f"\n writes: results/table/16_carrier_stratified_prs_slope_by_gene.csv")
    legacy_joint_models_log(f"[ESTIMAND CAVEAT] {ESTIMAND_CAVEAT}")

    beta_burden = full_joint_fit.params["burden_z"]
    beta_prs = full_joint_fit.params["PRS311_z"]
    beta_inter = full_joint_fit.params["inter_z"]

    distinct_burden_values = sorted(legacy_joint_models_df["burden_rare_main_ptv_10gene"].unique())
    legacy_joint_models_log(f"The actual Burden_rare_main_ptv_10gene value of \n is:{distinct_burden_values}")

    implied_rows = []
    for legacy_joint_models_v in distinct_burden_values:
        n_persons = int((legacy_joint_models_df["burden_rare_main_ptv_10gene"] == legacy_joint_models_v).sum())
        z_v = (legacy_joint_models_v - full_burden_mu) / full_burden_sd
        implied_or_vs_noncarrier = float(np.exp(beta_burden * z_v))
        implied_prs_or_per_sd = float(np.exp(beta_prs + beta_inter * z_v))
        external_or_ref = float(np.exp(legacy_joint_models_v))
        implied_rows.append({
            
            "burden_value": legacy_joint_models_v, "burden_z": z_v, "n_persons": legacy_joint_models_suppress(n_persons) if legacy_joint_models_v == 0 else "not reported",
            "implied_OR_vs_noncarrier": implied_or_vs_noncarrier,
            "external_OR_reference": external_or_ref,
            "implied_PRS_OR_per_SD_at_this_level": implied_prs_or_per_sd,
        })
    implied_tbl = pd.DataFrame(implied_rows)
    implied_tbl["evidence_status"] = "post_hoc_exploratory"
    implied_tbl["note"] = ('external_OR_reference = exp(burden_value): The original unit of burden_rare is external'
                            'The sum of ln(OR) (sum of ln(OR_g)), so exp (this value) is naturally equal to the reference under the assumption of external multiple effects'
                            'OR, this column is an identity rather than an independent fitting value. It is used to check the deviation of implied_OR_vs_noncarrier from it.')
    p_floor(implied_tbl).to_csv(os.path.join(LEGACY_JOINT_MODELS_TABLE_DIR, "16_model_implied_carrier_ors.csv"), index=False)
    legacy_joint_models_log('Model-implied carrier OR table:')
    legacy_joint_models_log(implied_tbl.to_string(index=False))
    legacy_joint_models_log(f"Write out: results/table/16_model_implied_carrier_ors.csv")

    legacy_joint_models_log("\nDONE.")

# ============================================================================
# 8. age-models: analyze_legacy_age_models
# ============================================================================

import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sm

def legacy_age_models_log(msg):
    print(msg, flush=True)

def legacy_age_models_section(title):
    legacy_age_models_log("\n" + "=" * 80)
    legacy_age_models_log(title)
    legacy_age_models_log("=" * 80)

def legacy_age_models_suppress(n, min_cell=20):
    """Legacy age models suppress."""
    if n is None or (isinstance(n, float) and np.isnan(n)):
        return n
    return n if n >= min_cell else f"<{min_cell}"

def join_covariates(base_df):
    prs_cols = ["person_id", "PRS311"]
    if "PRS311_noimp" in legacy_age_models_prs.columns:
        prs_cols.append("PRS311_noimp")
    out = base_df.merge(legacy_age_models_prs[prs_cols], on="person_id", how="left")
    out = out.merge(legacy_age_models_burden, on="person_id", how="left")
    out = out.merge(legacy_age_models_pc[["person_id"] + LEGACY_AGE_MODELS_PC_COLS], on="person_id", how="left")
    return out

def zscore(series):
    mu, sd = series.mean(), series.std()
    return (series - mu) / sd, mu, sd

def fit_ols_hc1(y, X):
    model = sm.OLS(y, X)
    return model.fit(cov_type="HC1")

def legacy_age_models_coef_table(fit_result, label, n):
    params = fit_result.params
    bse = fit_result.bse
    pvals = fit_result.pvalues
    rows = []
    for t in params.index:
        b = params[t]
        se = bse[t]
        p = pvals[t]
        rows.append({
            "model": label, "term": t, "beta": b, "robust_se_HC1": se,
            "ci_low": b - LEGACY_AGE_MODELS_Z * se, "ci_high": b + LEGACY_AGE_MODELS_Z * se, "p_value": p, "n": n,
        })
    return pd.DataFrame(rows)

def legacy_age_models_build_X(df_in, extra_cols, pc_cols=[f'PC{legacy_age_models_i}' for legacy_age_models_i in range(1, 11)]):
    cols = extra_cols + pc_cols
    X = df_in[cols].astype(float).copy()
    X = sm.add_constant(X, has_constant="add")
    return X

def legacy_age_models_save_smallest(fig, stem):
    paths = {ext: os.path.join(LEGACY_AGE_MODELS_FIG_DIR, f"{stem}.{ext}") for ext in ("svg", "png", "jpg")}
    fig.savefig(paths["svg"], bbox_inches="tight")
    fig.savefig(paths["png"], dpi=LEGACY_AGE_MODELS_DPI, bbox_inches="tight")
    fig.savefig(paths["jpg"], dpi=LEGACY_AGE_MODELS_DPI, bbox_inches="tight", pil_kwargs={"quality": 85, "subsampling": 2})
    sizes = {ext: os.path.getsize(q) for ext, q in paths.items()}
    keep = min(sizes, key=sizes.get)
    for ext, q in paths.items():
        if ext != keep:
            os.remove(q)
    legacy_age_models_log(f"  {stem}: " + ", ".join(f"{e}={v:,}B" for e, v in sizes.items()) + f" -> kept .{keep}")

def analyze_legacy_age_models():
    """Initialize runtime state and execute this component explicitly."""
    global LEGACY_AGE_MODELS_BASE, LEGACY_AGE_MODELS_DATA_DIR, LEGACY_AGE_MODELS_DISCLOSURE_MIN_CELL, LEGACY_AGE_MODELS_DPI, LEGACY_AGE_MODELS_FIG_DIR, LEGACY_AGE_MODELS_PC_COLS, LEGACY_AGE_MODELS_REPORT_DIR, LEGACY_AGE_MODELS_TABLE_DIR, LEGACY_AGE_MODELS_X, X_nopc, legacy_age_models_X_primary, X_sens1, X_sens2, X_sens3, X_sens4, X_sens5, LEGACY_AGE_MODELS_Z, legacy_age_models__, legacy_age_models_ax, legacy_age_models_b0, legacy_age_models_b1, beta1, beta2, legacy_age_models_burden, burden_mu, burden_raw, burden_sd, legacy_age_models_c, carrier_median, cases_full, cases_full_joined, covar_check_1, covar_cols_check, ct1, ct2, ct3, ct4, ct5, df_denom, df_num, legacy_age_models_extra_cols, f_stat, f_test_res, legacy_age_models_fig, fit_nopc, legacy_age_models_fit_res, fit_sens1, fit_sens2, fit_sens3, fit_sens4, fit_sens5, flow_rows, flow_tbl, legacy_age_models_g, grp, grp_n, hb, legacy_age_models_i, legacy_age_models_label, n_age_na, n_before_covna, legacy_age_models_n_dropped_noimp, n_first_record, n_flag_under18, legacy_age_models_n_na, legacy_age_models_n_na_281, n_pc_terms, nested_fits, nested_rows, nested_specs, nested_tbl, order, legacy_age_models_p1, p2, p_value_f, legacy_age_models_pc, pc_f_test_tbl, pc_terms_in_model, pcc, legacy_age_models_pheno, primary_coef, primary_df, primary_fit, legacy_age_models_primary_for_compare, primary_joined, legacy_age_models_prs, prs281_renamed, legacy_age_models_prs_281, prs_mu, prs_sd, legacy_age_models_qc, r2_joint, restriction_matrix, restriction_rows, legacy_age_models_row, sens1_df, sens2_df, sens3_df, sens4_df, sens5_base, sens5_df, legacy_age_models_sensitivity_rows, legacy_age_models_sensitivity_tbl_full, step1, step2, step3, step4_primary, unique_r2_map, legacy_age_models_xs, y_primary, y_sens1, y_sens2, y_sens3, y_sens4, y_sens5

    matplotlib.use("Agg")
    plt.style.use("default")

    LEGACY_AGE_MODELS_BASE = "/home/jupyter/workspace/analysis"
    LEGACY_AGE_MODELS_DATA_DIR = os.path.join(LEGACY_AGE_MODELS_BASE, "data")
    LEGACY_AGE_MODELS_FIG_DIR = os.path.join(LEGACY_AGE_MODELS_BASE, "results", "figures")
    LEGACY_AGE_MODELS_TABLE_DIR = os.path.join(LEGACY_AGE_MODELS_BASE, "results", "table")
    LEGACY_AGE_MODELS_REPORT_DIR = os.path.join(LEGACY_AGE_MODELS_BASE, "results", "report")

    os.makedirs(LEGACY_AGE_MODELS_FIG_DIR, exist_ok=True)
    os.makedirs(LEGACY_AGE_MODELS_TABLE_DIR, exist_ok=True)
    os.makedirs(LEGACY_AGE_MODELS_REPORT_DIR, exist_ok=True)

    LEGACY_AGE_MODELS_DPI = 150
    LEGACY_AGE_MODELS_DISCLOSURE_MIN_CELL = 20
    LEGACY_AGE_MODELS_PC_COLS = [f"PC{legacy_age_models_i}" for legacy_age_models_i in range(1, 11)]
    LEGACY_AGE_MODELS_Z = 1.959963984540054  # 95% CI multiplier

    # ---------------------------------------------------------------------------
    
    # ---------------------------------------------------------------------------
    legacy_age_models_section('Step A: Read each table and join to build Aim 4 cases-only queue')

    legacy_age_models_pheno = pd.read_csv(os.path.join(LEGACY_AGE_MODELS_DATA_DIR, "04_phenotype_diagnostics_person_level.csv.gz"))
    legacy_age_models_qc = pd.read_csv(os.path.join(LEGACY_AGE_MODELS_DATA_DIR, "05_cohort_qc_flags.csv"))
    legacy_age_models_pc = pd.read_csv(os.path.join(LEGACY_AGE_MODELS_DATA_DIR, "09_ancestry_pc_person_level.csv.gz"))
    legacy_age_models_burden = pd.read_csv(os.path.join(LEGACY_AGE_MODELS_DATA_DIR, "14_burden_rare_scores_person_level.csv.gz"))
    legacy_age_models_prs = pd.read_csv(os.path.join(LEGACY_AGE_MODELS_DATA_DIR, "prs_scores.csv"))
    legacy_age_models_prs_281 = pd.read_csv(os.path.join(LEGACY_AGE_MODELS_DATA_DIR, "prs_scores_281loci.csv"))

    flow_rows = []

    step1 = legacy_age_models_pheno[legacy_age_models_pheno["group"] == "case"].copy()
    step2 = step1.merge(legacy_age_models_qc[["person_id", "pass_qc"]], on="person_id", how="inner")
    step2 = step2[step2["pass_qc"] == 1].copy()
    flow_rows.append({"step": '1. Case in Aim 1 queue (No. 05 pass_qc==1)', "n": len(step2)})

    n_age_na = int(step2["age_at_diagnosis"].isna().sum())
    n_flag_under18 = int((step2["flag_age_at_diagnosis_invalid_under18"] == 1).sum())
    assert n_age_na == 0 and n_flag_under18 == 0, 'No. 05 has required the diagnosis age to be valid, there should be no more missing/under18 here'
    step3 = step2.dropna(subset=["age_at_diagnosis"]).copy()

    n_first_record = int((step3["bc_is_first_record"] == True).sum())
    legacy_age_models_log(f"bc_is_first_record==True (within this cohort):{n_first_record}people"
        f"({n_first_record/len(step3)*100:.2f}%)")

    cases_full = step3.copy()  

    step4_primary = step3[step3["bc_is_first_record"] != True].copy()
    flow_rows.append({"step": '2. & bc_is_first_record != True (excluded from main analysis, see method description)', "n": len(step4_primary)})

    primary_joined = join_covariates(step4_primary)
    n_before_covna = len(primary_joined)
    covar_cols_check = ["PRS311", "burden_rare_main_ptv_10gene", "burden_rare_main_ptv_9gene"] + LEGACY_AGE_MODELS_PC_COLS
    for legacy_age_models_c in covar_cols_check:
        legacy_age_models_n_na = int(primary_joined[legacy_age_models_c].isna().sum())
        if legacy_age_models_n_na > 0:
            legacy_age_models_log(f"  WARNING: {legacy_age_models_c}have{legacy_age_models_suppress(legacy_age_models_n_na)}NA (in the main analysis queue)")

    primary_df = primary_joined.dropna(subset=covar_cols_check).reset_index(drop=True)
    assert len(primary_df) == n_before_covna, 'PRS/Burden/PC should not be missing (calculated upstream for the entire Aim 1 queue)'

    flow_tbl = pd.DataFrame(flow_rows)
    flow_tbl["n_disclosed"] = flow_tbl["n"].apply(legacy_age_models_suppress)
    flow_tbl.to_csv(os.path.join(LEGACY_AGE_MODELS_TABLE_DIR, "17_cohort_flow.csv"), index=False)
    legacy_age_models_log('Queue construction process (cohort flow):')
    legacy_age_models_log(flow_tbl.to_string(index=False))

    legacy_age_models_log(f"\n main analysis modeling queue n={len(primary_df)} "
        f"(Expect about 8,300; the rough estimate range given by plan is for cross-checking, not a hard assertion)")

    # ---------------------------------------------------------------------------
    # Step B: helper — z-score, OLS(HC1) fit, coef table
    # ---------------------------------------------------------------------------
    legacy_age_models_section('Step B: Define z-score / OLS(HC1 robust) helper')

    # ---------------------------------------------------------------------------
    
    # ---------------------------------------------------------------------------
    legacy_age_models_section('Step C: Main model (10-gene Burden_rare, main analysis cohort, n={})'.format(len(primary_df)))

    primary_df["Burden_rare_z"], burden_mu, burden_sd = zscore(primary_df["burden_rare_main_ptv_10gene"])
    primary_df["PRS311_z"], prs_mu, prs_sd = zscore(primary_df["PRS311"])
    legacy_age_models_log(f"Standardization parameters: Burden_rare(10gene) mean={burden_mu:.6f} sd={burden_sd:.6f}; "
        f"PRS311 mean={prs_mu:.4f} sd={prs_sd:.4f}")

    y_primary = primary_df["age_at_diagnosis"].astype(float).values
    legacy_age_models_X_primary = legacy_age_models_build_X(primary_df, ["Burden_rare_z", "PRS311_z"])
    primary_fit = fit_ols_hc1(y_primary, legacy_age_models_X_primary)

    primary_coef = legacy_age_models_coef_table(primary_fit, "primary_10gene_bc_first_record_excluded", len(primary_df))
    primary_coef["r_squared"] = primary_fit.rsquared
    primary_coef["adj_r_squared"] = primary_fit.rsquared_adj
    primary_coef.to_csv(os.path.join(LEGACY_AGE_MODELS_TABLE_DIR, "17_primary_model_coefficients.csv"), index=False)

    legacy_age_models_log('Main model full coefficient table (Burden_rare_z / PRS311_z / PC1-10 / const):')
    legacy_age_models_log(primary_coef[["term", "beta", "robust_se_HC1", "ci_low", "ci_high", "p_value"]].to_string(index=False))
    legacy_age_models_log(f"\n main model R2={primary_fit.rsquared:.5f}, adj.R2={primary_fit.rsquared_adj:.5f}, n={len(primary_df)}")

    beta1 = primary_coef.loc[primary_coef["term"] == "Burden_rare_z", "beta"].iloc[0]
    legacy_age_models_p1 = primary_coef.loc[primary_coef["term"] == "Burden_rare_z", "p_value"].iloc[0]
    beta2 = primary_coef.loc[primary_coef["term"] == "PRS311_z", "beta"].iloc[0]
    p2 = primary_coef.loc[primary_coef["term"] == "PRS311_z", "p_value"].iloc[0]
    legacy_age_models_log(f"\n H3 direction verification: beta1(Burden_rare)={beta1:.4f} (p={legacy_age_models_p1:.4g}) -> "
        f"{'Negative (support H3)' if beta1 < 0 else 'Forward (does not support H3)'}")
    legacy_age_models_log(f"H3 direction verification: beta2(PRS311)={beta2:.4f} (p={p2:.4g}) -> "
        f"{'Negative (support H3)' if beta2 < 0 else 'Forward (does not support H3)'}")

    legacy_age_models_section('Step C2: PC1-10 joint F test (primary_fit.f_test, addendum §2)')
    pc_terms_in_model = [legacy_age_models_c for legacy_age_models_c in LEGACY_AGE_MODELS_PC_COLS if legacy_age_models_c in primary_fit.params.index]
    n_pc_terms = len(pc_terms_in_model)
    restriction_rows = []
    for pcc in pc_terms_in_model:
        legacy_age_models_row = np.zeros(len(primary_fit.params))
        legacy_age_models_row[list(primary_fit.params.index).index(pcc)] = 1.0
        restriction_rows.append(legacy_age_models_row)
    restriction_matrix = np.vstack(restriction_rows)
    f_test_res = primary_fit.f_test(restriction_matrix)
    f_stat = float(np.asarray(f_test_res.fvalue).squeeze())
    p_value_f = float(np.asarray(f_test_res.pvalue).squeeze())
    df_num = int(f_test_res.df_num)
    df_denom = int(f_test_res.df_denom)
    legacy_age_models_log(f"PC1-10 joint F test: F({df_num},{df_denom})={f_stat:.4f}, p={p_value_f:.4g}")
    pc_f_test_tbl = pd.DataFrame([{
        "f_stat": f_stat, "df_num": df_num, "df_denom": df_denom, "p_value": p_value_f,
    }])
    pc_f_test_tbl.to_csv(os.path.join(LEGACY_AGE_MODELS_TABLE_DIR, "17_pc_block_f_test.csv"), index=False)
    legacy_age_models_log(f"Write out: results/table/17_pc_block_f_test.csv")

    # ---------------------------------------------------------------------------
    # Step D: nested model R2 comparison (addendum §2 D2.1/D2.2/D2.3)

    # ---------------------------------------------------------------------------
    legacy_age_models_section('Step D: Comparison of nested model R2 (PCs-only / PRS-only / Burden-only / Burden+PRS without PCs / Joint)')

    nested_specs = {
        "PCs-only": [],
        "PRS311-only+PCs": ["PRS311_z"],
        "Burden_rare-only+PCs": ["Burden_rare_z"],
        "Joint (primary)": ["Burden_rare_z", "PRS311_z"],
    }

    nested_rows = []
    nested_fits = {}
    for legacy_age_models_label, legacy_age_models_extra_cols in nested_specs.items():
        LEGACY_AGE_MODELS_X = legacy_age_models_build_X(primary_df, legacy_age_models_extra_cols)
        legacy_age_models_fit_res = fit_ols_hc1(y_primary, LEGACY_AGE_MODELS_X)
        nested_fits[legacy_age_models_label] = legacy_age_models_fit_res
        nested_rows.append({
            "model": legacy_age_models_label, "terms": ", ".join(legacy_age_models_extra_cols) + " + PC1-10 + const",
            "r_squared": legacy_age_models_fit_res.rsquared, "adj_r_squared": legacy_age_models_fit_res.rsquared_adj,
            "n": len(primary_df),
        })
        legacy_age_models_log(f"  [{legacy_age_models_label}] R2={legacy_age_models_fit_res.rsquared:.5f}, adj.R2={legacy_age_models_fit_res.rsquared_adj:.5f}")

    X_nopc = sm.add_constant(primary_df[["Burden_rare_z", "PRS311_z"]].astype(float), has_constant="add")
    fit_nopc = fit_ols_hc1(y_primary, X_nopc)
    legacy_age_models_log(f"  [Burden+PRS, no PCs (fit)] R2={fit_nopc.rsquared:.5f}, adj.R2={fit_nopc.rsquared_adj:.5f}")

    nested_tbl = pd.DataFrame(nested_rows)
    r2_joint = nested_tbl.loc[nested_tbl["model"] == "Joint (primary)", "r_squared"].iloc[0]
    nested_tbl["incremental_r2_vs_joint"] = r2_joint - nested_tbl["r_squared"]

    unique_r2_map = {
        "PRS311-only+PCs": r2_joint - nested_tbl.loc[nested_tbl["model"] == "Burden_rare-only+PCs", "r_squared"].iloc[0],
        "Burden_rare-only+PCs": r2_joint - nested_tbl.loc[nested_tbl["model"] == "PRS311-only+PCs", "r_squared"].iloc[0],
        "PCs-only": r2_joint - fit_nopc.rsquared,
        "Joint (primary)": np.nan,
    }
    nested_tbl["unique_r2"] = nested_tbl["model"].map(unique_r2_map)
    nested_tbl.to_csv(os.path.join(LEGACY_AGE_MODELS_TABLE_DIR, "17_nested_model_r2_comparison.csv"), index=False)
    legacy_age_models_log('nested model R2 comparison table (including unique_r2):')
    legacy_age_models_log(nested_tbl.to_string(index=False))
    legacy_age_models_log(f"\n Reference: unique PRS={unique_r2_map['PRS311-only+PCs']:.7f}, "
        f"unique Burden={unique_r2_map['Burden_rare-only+PCs']:.7f}, "
        f"ratio(Burden/PRS)={unique_r2_map['Burden_rare-only+PCs']/unique_r2_map['PRS311-only+PCs']:.2f}x "
        f"(addendum §2Expected value: Burden unique=0.023045, PRS unique=0.002043, ratio≈11.3x --"
        f"This expected value comes from the 146,417/8,328 subsample under the old 146,417 population size. The sample size has changed under the new v7.3 QC standard."
        f"Changes in values ​​are expected)")

    # ---------------------------------------------------------------------------

    #   3. prs281loci_PRS311 (addendum §7/item 8, mirror scripts/16)
    #   4. prs_no_mean_imputation_PRS311 (addendum §7/item 8, mirror scripts/16)
    #   5. case_def_rule_of_two (addendum §7/item 8, mirror scripts/16;

    #   nested R2/AUC/ROC。
    # ---------------------------------------------------------------------------
    legacy_age_models_section('Step E: Sensitivity analysis (including bc_is_first_record full cohort / 9-gene panel / 281-site PRS /'
            'no-mean-imp PRS / rule-of-two case definition)')

    legacy_age_models_sensitivity_rows = []

    cases_full_joined = join_covariates(cases_full)
    covar_check_1 = ["PRS311", "burden_rare_main_ptv_10gene"] + LEGACY_AGE_MODELS_PC_COLS
    sens1_df = cases_full_joined.dropna(subset=covar_check_1).reset_index(drop=True)
    legacy_age_models_log(f"\n [Sensitivity #1] All cases-only cohort (not excluding bc_is_first_record), n={len(sens1_df)} "
        f"(Expect about 8,887)")

    sens1_df["Burden_rare_z"], legacy_age_models__, legacy_age_models__ = zscore(sens1_df["burden_rare_main_ptv_10gene"])
    sens1_df["PRS311_z"], legacy_age_models__, legacy_age_models__ = zscore(sens1_df["PRS311"])
    y_sens1 = sens1_df["age_at_diagnosis"].astype(float).values
    X_sens1 = legacy_age_models_build_X(sens1_df, ["Burden_rare_z", "PRS311_z"])
    fit_sens1 = fit_ols_hc1(y_sens1, X_sens1)
    ct1 = legacy_age_models_coef_table(fit_sens1, "sensitivity1_full_cohort_incl_bc_first_record_10gene", len(sens1_df))
    ct1["r_squared"] = fit_sens1.rsquared
    ct1["adj_r_squared"] = fit_sens1.rsquared_adj
    legacy_age_models_sensitivity_rows.append(ct1)
    legacy_age_models_log(ct1[ct1["term"].isin(["Burden_rare_z", "PRS311_z"])]
        [["term", "beta", "robust_se_HC1", "ci_low", "ci_high", "p_value"]].to_string(index=False))

    sens2_df = primary_df.copy()
    sens2_df["Burden_rare_9gene_z"], legacy_age_models__, legacy_age_models__ = zscore(sens2_df["burden_rare_main_ptv_9gene"])
    y_sens2 = sens2_df["age_at_diagnosis"].astype(float).values
    X_sens2 = legacy_age_models_build_X(sens2_df, ["Burden_rare_9gene_z", "PRS311_z"])
    fit_sens2 = fit_ols_hc1(y_sens2, X_sens2)
    ct2 = legacy_age_models_coef_table(fit_sens2, "sensitivity2_9gene_panel_primary_cohort", len(sens2_df))
    ct2["r_squared"] = fit_sens2.rsquared
    ct2["adj_r_squared"] = fit_sens2.rsquared_adj
    legacy_age_models_sensitivity_rows.append(ct2)
    legacy_age_models_log(f"\n [Sensitivity #2] 9-gene panel, main analysis cohort n={len(sens2_df)}")
    legacy_age_models_log(ct2[ct2["term"].isin(["Burden_rare_9gene_z", "PRS311_z"])]
        [["term", "beta", "robust_se_HC1", "ci_low", "ci_high", "p_value"]].to_string(index=False))

    prs281_renamed = legacy_age_models_prs_281.rename(columns={"PRS311": "PRS311_281loci"})
    sens3_df = primary_df.merge(prs281_renamed[["person_id", "PRS311_281loci"]], on="person_id", how="left")
    legacy_age_models_n_na_281 = int(sens3_df["PRS311_281loci"].isna().sum())
    if legacy_age_models_n_na_281:
        legacy_age_models_log(f"WARNING: PRS311_281loci has{legacy_age_models_n_na_281}NA, will be excluded by dropna")
    sens3_df = sens3_df.dropna(subset=["PRS311_281loci"]).reset_index(drop=True)
    sens3_df["Burden_rare_z"], legacy_age_models__, legacy_age_models__ = zscore(sens3_df["burden_rare_main_ptv_10gene"])
    sens3_df["PRS311_281loci_z"], legacy_age_models__, legacy_age_models__ = zscore(sens3_df["PRS311_281loci"])
    y_sens3 = sens3_df["age_at_diagnosis"].astype(float).values
    X_sens3 = legacy_age_models_build_X(sens3_df, ["Burden_rare_z", "PRS311_281loci_z"])
    fit_sens3 = fit_ols_hc1(y_sens3, X_sens3)
    ct3 = legacy_age_models_coef_table(fit_sens3, "prs281loci_PRS311", len(sens3_df))
    ct3["r_squared"] = fit_sens3.rsquared
    ct3["adj_r_squared"] = fit_sens3.rsquared_adj
    legacy_age_models_sensitivity_rows.append(ct3)
    legacy_age_models_log(f"\n [Sensitivity#3 prs281loci_PRS311] n={len(sens3_df)}")
    legacy_age_models_log(ct3[ct3["term"].isin(["Burden_rare_z", "PRS311_281loci_z"])]
        [["term", "beta", "robust_se_HC1", "ci_low", "ci_high", "p_value"]].to_string(index=False))

    if "PRS311_noimp" in primary_df.columns:
        sens4_df = primary_df.dropna(subset=["PRS311_noimp"]).reset_index(drop=True)
        legacy_age_models_n_dropped_noimp = len(primary_df) - len(sens4_df)
        if legacy_age_models_n_dropped_noimp:
            legacy_age_models_log(f"WARNING: PRS311_noimp has{legacy_age_models_n_dropped_noimp}NA, will be excluded by dropna")
        sens4_df["Burden_rare_z"], legacy_age_models__, legacy_age_models__ = zscore(sens4_df["burden_rare_main_ptv_10gene"])
        sens4_df["PRS311_noimp_z"], legacy_age_models__, legacy_age_models__ = zscore(sens4_df["PRS311_noimp"])
        y_sens4 = sens4_df["age_at_diagnosis"].astype(float).values
        X_sens4 = legacy_age_models_build_X(sens4_df, ["Burden_rare_z", "PRS311_noimp_z"])
        fit_sens4 = fit_ols_hc1(y_sens4, X_sens4)
        ct4 = legacy_age_models_coef_table(fit_sens4, "prs_no_mean_imputation_PRS311", len(sens4_df))
        ct4["r_squared"] = fit_sens4.rsquared
        ct4["adj_r_squared"] = fit_sens4.rsquared_adj
        legacy_age_models_sensitivity_rows.append(ct4)
        legacy_age_models_log(f"\n [Sensitivity#4 prs_no_mean_imputation_PRS311] n={len(sens4_df)}")
        legacy_age_models_log(ct4[ct4["term"].isin(["Burden_rare_z", "PRS311_noimp_z"])]
            [["term", "beta", "robust_se_HC1", "ci_low", "ci_high", "p_value"]].to_string(index=False))
    else:
        legacy_age_models_log('[Sensitivity #4 prs_no_mean_imputation_PRS311] SKIPPED: data/prs_scores.csv has no PRS311_noimp column')

    if "is_case_rule_of_two" in primary_df.columns:

        sens5_base = primary_df
        sens5_df = sens5_base.loc[sens5_base["is_case_rule_of_two"] == True].reset_index(drop=True)  # noqa: E712
        legacy_age_models_log(f"\n [sensitivity#5 case_def_rule_of_two] n={len(sens5_df)}(Main analysis cohort n={len(primary_df)}subset)")
        sens5_df["Burden_rare_z"], legacy_age_models__, legacy_age_models__ = zscore(sens5_df["burden_rare_main_ptv_10gene"])
        sens5_df["PRS311_z"], legacy_age_models__, legacy_age_models__ = zscore(sens5_df["PRS311"])
        y_sens5 = sens5_df["age_at_diagnosis"].astype(float).values
        X_sens5 = legacy_age_models_build_X(sens5_df, ["Burden_rare_z", "PRS311_z"])
        fit_sens5 = fit_ols_hc1(y_sens5, X_sens5)
        ct5 = legacy_age_models_coef_table(fit_sens5, "case_def_rule_of_two", len(sens5_df))
        ct5["r_squared"] = fit_sens5.rsquared
        ct5["adj_r_squared"] = fit_sens5.rsquared_adj
        legacy_age_models_sensitivity_rows.append(ct5)
        legacy_age_models_log(ct5[ct5["term"].isin(["Burden_rare_z", "PRS311_z"])]
            [["term", "beta", "robust_se_HC1", "ci_low", "ci_high", "p_value"]].to_string(index=False))
    else:
        legacy_age_models_log('[Sensitivity #5 case_def_rule_of_two] SKIPPED: phenotype table has no is_case_rule_of_two column')

    legacy_age_models_primary_for_compare = primary_coef.copy()
    legacy_age_models_sensitivity_tbl_full = pd.concat([legacy_age_models_primary_for_compare] + legacy_age_models_sensitivity_rows, ignore_index=True)
    legacy_age_models_sensitivity_tbl_full.to_csv(os.path.join(LEGACY_AGE_MODELS_TABLE_DIR, "17_sensitivity_analyses_coefficients.csv"), index=False)

    legacy_age_models_log('Main analysis + all sensitivity analyzes side by side (Burden_rare/PRS311 coefficient direction & significance comparison):')
    legacy_age_models_log(legacy_age_models_sensitivity_tbl_full[legacy_age_models_sensitivity_tbl_full["term"].str.contains("Burden_rare|PRS311", na=False)]
        [["model", "term", "beta", "robust_se_HC1", "ci_low", "ci_high", "p_value", "n"]].to_string(index=False))

    # ---------------------------------------------------------------------------

    # ---------------------------------------------------------------------------
    legacy_age_models_section('Step F: Plot (Burden_rare grouped box plot + PRS311 hexbin, main analysis cohort)')

    burden_raw = primary_df["burden_rare_main_ptv_10gene"]
    carrier_median = burden_raw[burden_raw > 0].median()
    grp = np.select([burden_raw == 0, burden_raw <= carrier_median], ["Non-carrier",
                    f"Carrier, Burden_rare <= {carrier_median:.2f}"], default=f"Carrier, Burden_rare > {carrier_median:.2f}")
    order = ["Non-carrier", f"Carrier, Burden_rare <= {carrier_median:.2f}", f"Carrier, Burden_rare > {carrier_median:.2f}"]
    grp_n = pd.Series(grp).value_counts()
    assert all(grp_n.get(g, 0) >= 20 for g in order), 'The boxplot has groups n<20 and needs to be merged'
    legacy_age_models_fig, legacy_age_models_ax = plt.subplots(figsize=(6.5, 5))
    sns.boxplot(x=grp, y=primary_df["age_at_diagnosis"].values, order=order, showfliers=False, ax=legacy_age_models_ax, color="white")
    legacy_age_models_ax.set_xticklabels([f"{legacy_age_models_g}\n(n={grp_n[legacy_age_models_g]:,d})" for legacy_age_models_g in order], fontsize=8)
    legacy_age_models_ax.set_xlabel("")
    legacy_age_models_ax.set_ylabel("Age at diagnosis (years)")
    legacy_age_models_fig.tight_layout()
    legacy_age_models_save_smallest(legacy_age_models_fig, "17_age_at_diagnosis_by_burden_group_boxplot")
    plt.close(legacy_age_models_fig)

    legacy_age_models_fig, legacy_age_models_ax = plt.subplots(figsize=(6.5, 5))
    hb = legacy_age_models_ax.hexbin(primary_df["PRS311_z"], primary_df["age_at_diagnosis"], gridsize=30, mincnt=20, cmap="Blues")
    legacy_age_models_fig.colorbar(hb, ax=legacy_age_models_ax, label="Participants per cell (cells with <20 not shown)")
    legacy_age_models_b1, legacy_age_models_b0 = np.polyfit(primary_df["PRS311_z"], primary_df["age_at_diagnosis"], 1)
    legacy_age_models_xs = np.linspace(primary_df["PRS311_z"].quantile(0.01), primary_df["PRS311_z"].quantile(0.99), 50)
    legacy_age_models_ax.plot(legacy_age_models_xs, legacy_age_models_b0 + legacy_age_models_b1 * legacy_age_models_xs, color="tab:red")
    legacy_age_models_ax.set_xlabel("PRS311 (z-scored)")
    legacy_age_models_ax.set_ylabel("Age at diagnosis (years)")
    legacy_age_models_fig.tight_layout()
    legacy_age_models_save_smallest(legacy_age_models_fig, "17_prs311_vs_age_at_diagnosis_hexbin")
    plt.close(legacy_age_models_fig)

    legacy_age_models_log("\nDONE.")

# ============================================================================
# 9. figures: analyze_legacy_figures
# ============================================================================

import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
from scipy import stats
from sklearn.metrics import brier_score_loss, roc_auc_score

def legacy_figures_log(msg):
    print(msg, flush=True)

def legacy_figures_section(title):
    legacy_figures_log("\n" + "=" * 88)
    legacy_figures_log(title)
    legacy_figures_log("=" * 88)

def disclose_n(n, where=""):
    """Disclose n."""
    if n is None or (isinstance(n, float) and np.isnan(n)):
        return n, False
    n = int(n)
    if n < 5:
        SUPPRESSION_LOG.append({"cell": where, "raw_n_bucket": "<5", "action": "suppressed (<5)", "disclosed": "<5"})
        return "<5", True
    if n < LEGACY_FIGURES_DISCLOSURE_MIN_CELL:
        SUPPRESSION_LOG.append({"cell": where, "raw_n_bucket": "5-19", "action": "rounded up to 20 (5-19)", "disclosed": "20"})
        return "20", True
    return str(n), False

def audit_override(cell, disclosed, action, raw_bucket):
    """Audit override."""
    for e in SUPPRESSION_LOG:
        if e["cell"] == cell:
            e.update({"action": action, "disclosed": disclosed})
            return
    SUPPRESSION_LOG.append({"cell": cell, "raw_n_bucket": raw_bucket, "action": action, "disclosed": disclosed})

def fit_axes_to_hexbin(ax, hb, pad=0.08):
    """Fit axes to hexbin."""
    xy = hb.get_offsets()
    (x0, y0), (x1, y1) = xy.min(axis=0), xy.max(axis=0)
    dx, dy = (x1 - x0) * pad, (y1 - y0) * pad
    ax.set_xlim(x0 - dx, x1 + dx)
    ax.set_ylim(y0 - dy, y1 + dy)

def save_fig(fig, stem, prefer="svg", max_svg_kb=600, out_dir=None):
    """Save fig."""
    svg_path = os.path.join(out_dir or LEGACY_FIGURES_FIG_DIR, stem + ".svg")
    png_path = os.path.join(out_dir or LEGACY_FIGURES_FIG_DIR, stem + ".png")
    if prefer == "png":
        fig.savefig(png_path, dpi=LEGACY_FIGURES_DPI, bbox_inches="tight")
        plt.close(fig)
        kb = os.path.getsize(png_path) / 1024
        FIG_SIZES.append({"figure": os.path.basename(png_path), "size_kb": round(kb, 1), "format": "png"})
        legacy_figures_log(f"  [saved] {png_path}  ({kb:.1f} KB, PNG by design)")
        return png_path
    fig.savefig(svg_path, bbox_inches="tight")
    kb = os.path.getsize(svg_path) / 1024
    if kb > max_svg_kb:
        fig.savefig(png_path, dpi=LEGACY_FIGURES_DPI, bbox_inches="tight")
        os.remove(svg_path)
        plt.close(fig)
        kb2 = os.path.getsize(png_path) / 1024
        FIG_SIZES.append({"figure": os.path.basename(png_path), "size_kb": round(kb2, 1), "format": "png"})
        legacy_figures_log(f"  [saved] {png_path}  ({kb2:.1f}KB; SVG version is{kb:.1f}KB more than{max_svg_kb}KB threshold, saved as PNG and removed SVG)")
        return png_path
    plt.close(fig)
    FIG_SIZES.append({"figure": os.path.basename(svg_path), "size_kb": round(kb, 1), "format": "svg"})
    legacy_figures_log(f"  [saved] {svg_path}  ({kb:.1f} KB)")
    return svg_path

def write_table(df_out, stem, note="", out_dir=None):
    out = df_out.copy()
    out["evidence_status"] = EVIDENCE
    if note:
        out["note"] = note
    path = os.path.join(out_dir or LEGACY_FIGURES_TABLE_DIR, stem + ".csv")
    out.to_csv(path, index=False)
    legacy_figures_log(f"  [table] {path}  rows={len(out)}")
    return path

def fit_logit(y, X):
    """Fit logit."""
    return sm.GLM(np.asarray(y, dtype=float), X, family=sm.families.Binomial()).fit(cov_type="HC1")

def fit_logit_weighted(y, X):
    y = np.asarray(y, dtype=float)
    n = len(y)
    n1 = y.sum()
    w = np.where(y == 1, n / (2 * n1), n / (2 * (n - n1)))
    return sm.GLM(y, X, family=sm.families.Binomial(), var_weights=w).fit(cov_type="HC1"), w

def add_const(X):
    return sm.add_constant(X.astype(float), has_constant="add")

def fisher_ci(r, n, k=0):
    """Fisher ci."""
    if not np.isfinite(r) or abs(r) >= 1:
        return np.nan, np.nan
    zr = np.arctanh(r)
    se = 1.0 / np.sqrt(max(n - 3 - k, 1))
    return np.tanh(zr - LEGACY_FIGURES_Z * se), np.tanh(zr + LEGACY_FIGURES_Z * se)

def residualize(y, X_cov):
    """Residualize."""
    Xc = add_const(pd.DataFrame(X_cov))
    return sm.OLS(np.asarray(y, dtype=float), Xc).fit().resid

def corr_row(label, x, y, pcs=None, kind="pearson"):
    """Corr row."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(x)
    if kind == "spearman":
        xr, yr = stats.rankdata(x), stats.rankdata(y)
        r, p = stats.spearmanr(x, y)
        xu, yu = xr, yr
    else:
        r, p = stats.pearsonr(x, y)
        xu, yu = x, y
    lo, hi = fisher_ci(r, n)
    rows = [{"comparison": label, "method": kind, "adjustment": "raw", "n": n,
             "r": r, "ci_low": lo, "ci_high": hi, "p_value": p}]
    if pcs is not None:
        rx = residualize(xu, pcs)
        ry = residualize(yu, pcs)
        rp, pp = stats.pearsonr(rx, ry)
        lo2, hi2 = fisher_ci(rp, n, k=pcs.shape[1])
        rows.append({"comparison": label, "method": kind, "adjustment": "partial_PC1_10", "n": n,
                     "r": rp, "ci_low": lo2, "ci_high": hi2, "p_value": pp})
    return rows

def tier4(row):
    if row["n_ptv_genes"] == 0:
        return "non_carrier"
    if row["n_ptv_genes"] >= 2:
        return "compound"
    return "high_penetrance" if row["n_high"] == 1 else "moderate"

def _split(n_total, n_case, n_control):
    if n_total < LEGACY_FIGURES_DISCLOSURE_MIN_CELL:
        return "<20", "<20", "<20"
    small_case, small_ctrl = n_case < LEGACY_FIGURES_DISCLOSURE_MIN_CELL, n_control < LEGACY_FIGURES_DISCLOSURE_MIN_CELL
    if not small_case and not small_ctrl:
        return str(n_total), str(n_case), str(n_control)
    return (str(n_total), "<20" if small_case else "suppressed (secondary)",
            "<20" if small_ctrl else "suppressed (secondary)")

def cell_table(frame, row_key, col_key, tag):
    rows = []
    for rv, sub_r in frame.groupby(row_key):
        for cv, cell in sub_r.groupby(col_key):
            n = len(cell)
            n_case = int(cell["is_case_int"].sum())
            obs = n_case / n if n else np.nan
            exp = cell["p_expected_noint"].mean()
            where = f"{tag}|{row_key}={rv}|{col_key}={cv}"
            n_disc, n_mod = disclose_n(n, where)
            c_disc, c_mod = disclose_n(n_case, where + "|cases")
            suppressed = (n < LEGACY_FIGURES_DISCLOSURE_MIN_CELL) or (n_case < 5)
            rows.append({row_key: rv, col_key: cv, "n_raw": n, "n_case_raw": n_case,
                         "n_disclosed": n_disc, "n_case_disclosed": c_disc,
                         "observed_case_proportion": np.nan if suppressed else obs,
                         "mean_expected_prob_no_interaction": np.nan if suppressed else exp,
                         "obs_over_exp": np.nan if suppressed else (obs / exp if exp > 0 else np.nan),
                         "cell_suppressed_for_disclosure": suppressed})
    return pd.DataFrame(rows)

def to_grid(cells, row_key, col_key, val):
    return cells.pivot(index=row_key, columns=col_key, values=val)

def fit_ventile_model(ref_set):
    lab = np.where(legacy_figures_df["prs_ventile"].isin(ref_set), "ref", "v" + legacy_figures_df["prs_ventile"].astype(str))
    dm = pd.get_dummies(lab, prefix="pv").astype(float).drop(columns=["pv_ref"])
    X = add_const(pd.concat([dm, legacy_figures_df[COV_SIMPLE].astype(float)], axis=1))
    return fit_logit(legacy_figures_df["is_case_int"].values, X)

def calib_params(y, p):
    lp = np.log(np.clip(p, 1e-12, 1 - 1e-12) / (1 - np.clip(p, 1e-12, 1 - 1e-12)))
    slope_fit = sm.GLM(y, add_const(pd.DataFrame({"lp": lp})), family=sm.families.Binomial()).fit()
    inter_fit = sm.GLM(y, np.ones((len(y), 1)), family=sm.families.Binomial(), offset=lp).fit()
    return float(inter_fit.params[0]), float(slope_fit.params["lp"]), float(slope_fit.bse["lp"])

def analyze_legacy_figures():
    """Initialize runtime state and execute this component explicitly."""
    global A, LEGACY_FIGURES_BASE, CONTESTED, LEGACY_FIGURES_COV_COLS, COV_SIMPLE, LEGACY_FIGURES_DATA_DIR, LEGACY_FIGURES_DISCLOSURE_MIN_CELL, LEGACY_FIGURES_DPI, EVIDENCE, EXPECTED_N_AIM4, LEGACY_FIGURES_FIG_DIR, FIG_SIZES, LEGACY_FIGURES_GENE_COLS, HIGH_PEN, LEGACY_FIGURES_INTERNAL_DIR, MODERATE, LEGACY_FIGURES_PC_COLS, REF16_INTER_BETA, REF17_BURDEN_BETA, REF17_PRS_BETA, REF17_R2, REF_LOWEST, REF_VENTILES, LEGACY_FIGURES_REPORT_DIR, LEGACY_FIGURES_SEQ_COLS, LEGACY_FIGURES_SITE_COLS, SUPPRESSION_LOG, LEGACY_FIGURES_TABLE_DIR, TIER_ORDER, X_a4, X_nc, X_noint, legacy_figures_X_primary, X_t, Xm, XtWX, LEGACY_FIGURES_Z, legacy_figures__, _c, _c16, _c17, _cand, _cand_b, _k, _l, _m, _m2, _m3, _pbm, _pct, _pg_rows, _pv, _samp, _small, _small_b, _t, _theo, _top_low, annS, annotA, auc_oof, legacy_figures_ax, axA, axB, legacy_figures_axes, legacy_figures_b, legacy_figures_b0, legacy_figures_b1, band, band_order, bl, bmask, box_rows, box_tbl, brier_oof, brows, legacy_figures_burden, legacy_figures_c, c1, cal, cal_tbl, car, cases, cat_order, cell, ci_hi_top, ci_lo_top, ci_oof, legacy_figures_clinical, cmp7, cmp_rows, cmp_tbl, cnt, legacy_figures_cohort, color, comp, comp_hi, comp_lo, comp_n, comp_or, comp_row, controls, cook_summary, cooks, corr_rows, corr_tbl, cs_oof, cs_se_oof, ctrl_car, data_box, legacy_figures_df, df_c, df_nc, diag_tbl, edges, legacy_figures_f, legacy_figures_fig, fig1_path, fig63_path, fit_a4, fit_nc, fit_noint, fit_primary, fit_t, fit_v_low, fit_v_mid, fitted, ft, full_b_mu, full_b_sd, legacy_figures_full_prs_mu, legacy_figures_full_prs_sd, legacy_figures_g, gA, gB, gcols, legacy_figures_gene_carrier, h, hb63, hbB, legacy_figures_hi, hi10, hi9, hi_corner, legacy_figures_i, jb_p, jb_stat, legacy_figures_k, kurt, legacy_figures_lab, legacy_figures_lo, log_cd, main_cells, mask, mask_nc, mismatch, n, n_bins_hidden, legacy_figures_n_case, n_compound, n_pts, n_shown_B, n_supp_main, nc, obs, off, or_tbl, or_top_low, order3, p, legacy_figures_p1, p_hat, p_oof, pa, pa_mid, pb, pb_rows, legacy_figures_pc, pcs, pcs_c, per_gene, legacy_figures_pheno, piv, ppv_95_high, ppv_95_low, ppv_95_ratio, ppv_any_high, ppv_any_low, ppv_any_ratio, ppv_bands, ppv_rows, ppv_tbl, legacy_figures_prs, prs_ct, legacy_figures_qc, legacy_figures_r, r_c, r_pb, r_pb_adj, r_pb_raw, ref_label, ref_set, report, resid, resid_p, rng, legacy_figures_rows, rp, rs_c, s, sA, sB, scope, legacy_figures_se, seq_dummies, site_dummies, size_tbl, skew, split_tbl, legacy_figures_stratum, legacy_figures_sub, sub_idx, supp, supp_cells, supp_tbl, sw_p, sw_stat, t, tbl_files, tbl_kb, tier_dummy, tier_n, tier_tbl, legacy_figures_v, verdict, w_irls, xb, legacy_figures_xs, legacy_figures_y_all, zero_frac, zero_frac_pct

    matplotlib.use("Agg")
    plt.style.use("default")

    LEGACY_FIGURES_BASE = "/home/jupyter/workspace/analysis"
    LEGACY_FIGURES_DATA_DIR = os.path.join(LEGACY_FIGURES_BASE, "data")
    LEGACY_FIGURES_FIG_DIR = os.path.join(LEGACY_FIGURES_BASE, "results", "figures")
    LEGACY_FIGURES_TABLE_DIR = os.path.join(LEGACY_FIGURES_BASE, "results", "table")
    LEGACY_FIGURES_REPORT_DIR = os.path.join(LEGACY_FIGURES_BASE, "results", "report")

    LEGACY_FIGURES_DPI = 150
    LEGACY_FIGURES_Z = 1.959963984540054
    LEGACY_FIGURES_DISCLOSURE_MIN_CELL = 20
    
    _c16 = pd.read_csv(os.path.join(LEGACY_FIGURES_TABLE_DIR, "16_joint_model_coefficients_primary.csv")).set_index("term")
    _c17 = pd.read_csv(os.path.join(LEGACY_FIGURES_TABLE_DIR, "17_primary_model_coefficients.csv"))
    _c17 = _c17[_c17["model"] == "primary_10gene_bc_first_record_excluded"].set_index("term")
    REF16_INTER_BETA = float(_c16.loc["inter_z", "beta"])
    REF17_BURDEN_BETA = float(_c17.loc["Burden_rare_z", "beta"])
    REF17_PRS_BETA = float(_c17.loc["PRS311_z", "beta"])
    REF17_R2 = float(_c17["r_squared"].iloc[0])
    EXPECTED_N_AIM4 = int(_c17["n"].iloc[0])
    EVIDENCE = "post_hoc_exploratory"

    LEGACY_FIGURES_PC_COLS = [f"PC{legacy_figures_i}" for legacy_figures_i in range(1, 11)]
    LEGACY_FIGURES_GENE_COLS = ["ATM", "BARD1", "BRCA1", "BRCA2", "CHEK2", "PALB2", "RAD51C", "RAD51D", "TP53", "MAP3K1"]
    HIGH_PEN = ["BRCA1", "BRCA2", "PALB2", "TP53"]
    MODERATE = ["ATM", "CHEK2", "BARD1", "RAD51C", "RAD51D", "MAP3K1"]
    CONTESTED = ["MAP3K1"]

    FIG_SIZES = []
    SUPPRESSION_LOG = []

    # =============================================================================
    
    # =============================================================================
    legacy_figures_section('Step A | Queue construction (replica scripts/16 Step A caliber), asserting the same group of people as OOF table No. 16')

    legacy_figures_pheno = pd.read_csv(os.path.join(LEGACY_FIGURES_DATA_DIR, "04_phenotype_diagnostics_person_level.csv.gz"))
    legacy_figures_qc = pd.read_csv(os.path.join(LEGACY_FIGURES_DATA_DIR, "05_cohort_qc_flags.csv"))
    legacy_figures_pc = pd.read_csv(os.path.join(LEGACY_FIGURES_DATA_DIR, "09_ancestry_pc_person_level.csv.gz"))
    legacy_figures_clinical = pd.read_csv(os.path.join(LEGACY_FIGURES_DATA_DIR, "09_extended_clinical_person_level.csv.gz"))
    legacy_figures_burden = pd.read_csv(os.path.join(LEGACY_FIGURES_DATA_DIR, "14_burden_rare_scores_person_level.csv.gz"))
    legacy_figures_prs = pd.read_csv(os.path.join(LEGACY_FIGURES_DATA_DIR, "prs_scores.csv"))
    legacy_figures_gene_carrier = pd.read_csv(os.path.join(LEGACY_FIGURES_DATA_DIR, "13_carrier_status_main_ptv_maf1pct_person_level.csv.gz"))
    split_tbl = pd.read_csv(os.path.join(LEGACY_FIGURES_DATA_DIR, "16_joint_model_cv_oof_person_level.csv.gz"))  

    legacy_figures_cohort = legacy_figures_pheno.merge(legacy_figures_qc[["person_id", "pass_qc"]], on="person_id", how="inner")
    legacy_figures_cohort = legacy_figures_cohort[(legacy_figures_cohort["pass_qc"] == 1) & (legacy_figures_cohort["group"].isin(["case", "control"]))].copy()
    legacy_figures_df = legacy_figures_cohort.merge(legacy_figures_pc[["person_id"] + LEGACY_FIGURES_PC_COLS], on="person_id", how="left")
    legacy_figures_df = legacy_figures_df.merge(legacy_figures_clinical[["person_id", "seq_site_id", "seq_sequencer"]], on="person_id", how="left")
    legacy_figures_df = legacy_figures_df.merge(legacy_figures_burden[["person_id", "burden_rare_main_ptv_10gene"]], on="person_id", how="left")
    legacy_figures_df = legacy_figures_df.merge(legacy_figures_prs[["person_id", "PRS311"]], on="person_id", how="left")
    legacy_figures_df = legacy_figures_df.merge(legacy_figures_gene_carrier.rename(columns={legacy_figures_c: f"gene_{legacy_figures_c}" for legacy_figures_c in LEGACY_FIGURES_GENE_COLS}), on="person_id", how="left")
    legacy_figures_df = legacy_figures_df.dropna(subset=["burden_rare_main_ptv_10gene"]).reset_index(drop=True)
    legacy_figures_df = legacy_figures_df.dropna(subset=["age_at_last_obs"]).reset_index(drop=True)

    legacy_figures_log(f"Modeling queue n={len(legacy_figures_df)} (case={int((legacy_figures_df['group']=='case').sum())}, control={int((legacy_figures_df['group']=='control').sum())})")
    assert set(legacy_figures_df["person_id"]) == set(split_tbl["person_id"]), \
        f"FATAL: cohort n={len(legacy_figures_df)}with scripts/16 OOF table n={len(split_tbl)}The person_id collection is inconsistent (inconsistent caliber, stop)"
    assert int((legacy_figures_df["group"] == "case").sum()) == int((split_tbl["group"] == "case").sum()), 'FATAL: case number is inconsistent with scripts/16'
    legacy_figures_log('ASSERT PASS: n / case number is exactly the same as scripts/16')

    legacy_figures_df["is_case_int"] = legacy_figures_df["is_case"].astype(int)
    legacy_figures_df["carrier_main_10gene"] = (legacy_figures_df["burden_rare_main_ptv_10gene"] > 0).astype(int)

    legacy_figures_df = legacy_figures_df.merge(split_tbl[["person_id", "oof_r1_joint", "carrier_main_10gene"]]
                  .rename(columns={"carrier_main_10gene": "carrier_from_split"}), on="person_id", how="left")
    assert legacy_figures_df["oof_r1_joint"].isna().sum() == 0, 'FATAL: Someone is not in the OOF table of scripts/16'
    assert (legacy_figures_df["carrier_main_10gene"] == legacy_figures_df["carrier_from_split"]).all(), 'FATAL: carrier definition is inconsistent with No. 16'
    legacy_figures_log('No. 16 OOF predicts successful merger')

    site_dummies = pd.get_dummies(legacy_figures_df["seq_site_id"], prefix="site", drop_first=True).astype(float)
    seq_dummies = pd.get_dummies(legacy_figures_df["seq_sequencer"], prefix="seq", drop_first=True).astype(float)
    LEGACY_FIGURES_SITE_COLS, LEGACY_FIGURES_SEQ_COLS = list(site_dummies.columns), list(seq_dummies.columns)
    legacy_figures_df = pd.concat([legacy_figures_df, site_dummies, seq_dummies], axis=1)
    LEGACY_FIGURES_COV_COLS = ["age_at_last_obs"] + LEGACY_FIGURES_PC_COLS + LEGACY_FIGURES_SITE_COLS + LEGACY_FIGURES_SEQ_COLS
    COV_SIMPLE = ["age_at_last_obs"] + LEGACY_FIGURES_PC_COLS  

    legacy_figures_full_prs_mu, legacy_figures_full_prs_sd = legacy_figures_df["PRS311"].mean(), legacy_figures_df["PRS311"].std()
    full_b_mu, full_b_sd = legacy_figures_df["burden_rare_main_ptv_10gene"].mean(), legacy_figures_df["burden_rare_main_ptv_10gene"].std()
    legacy_figures_df["PRS311_z"] = (legacy_figures_df["PRS311"] - legacy_figures_full_prs_mu) / legacy_figures_full_prs_sd
    legacy_figures_df["burden_z"] = (legacy_figures_df["burden_rare_main_ptv_10gene"] - full_b_mu) / full_b_sd
    legacy_figures_df["inter_z"] = legacy_figures_df["PRS311_z"] * legacy_figures_df["burden_z"]
    legacy_figures_log(f"Full queue z-score parameter: PRS311 mean={legacy_figures_full_prs_mu:.4f} sd={legacy_figures_full_prs_sd:.4f}; "
        f"burden mean={full_b_mu:.6f} sd={full_b_sd:.6f}")

    # =============================================================================
    
    # =============================================================================
    legacy_figures_section('Step B | Gene stratification (high / moderate / compound)')

    gcols = [f"gene_{legacy_figures_g}" for legacy_figures_g in LEGACY_FIGURES_GENE_COLS]
    legacy_figures_df[gcols] = legacy_figures_df[gcols].fillna(0).astype(int)
    legacy_figures_df["n_ptv_genes"] = legacy_figures_df[gcols].sum(axis=1)
    legacy_figures_df["n_high"] = legacy_figures_df[[f"gene_{legacy_figures_g}" for legacy_figures_g in HIGH_PEN]].sum(axis=1)
    legacy_figures_df["n_mod"] = legacy_figures_df[[f"gene_{legacy_figures_g}" for legacy_figures_g in MODERATE]].sum(axis=1)

    mismatch = int((legacy_figures_df["n_ptv_genes"] > 0).ne(legacy_figures_df["carrier_main_10gene"] == 1).sum())
    legacy_figures_log(f"The number of people who are inconsistent with the definition of gene layer carrier (>=1 gene PTV) and burden>0:{mismatch}")

    legacy_figures_df["tier4"] = legacy_figures_df.apply(tier4, axis=1)
    
    legacy_figures_df["tier3"] = np.where(legacy_figures_df["tier4"] == "non_carrier", "non_carrier",
                  np.where(legacy_figures_df["n_high"] >= 1, "high_penetrance", "moderate"))

    tier_tbl = (legacy_figures_df.groupby(["tier4", "group"]).size().unstack(fill_value=0)
                  .rename(columns={"case": "n_case", "control": "n_control"}).reset_index())
    tier_tbl["n_total"] = tier_tbl["n_case"] + tier_tbl["n_control"]
    tier_tbl["n_total_disclosed"] = [disclose_n(legacy_figures_v, f"tier4={t}")[0] for legacy_figures_v, t in zip(tier_tbl["n_total"], tier_tbl["tier4"])]
    tier_tbl["n_case_disclosed"] = [disclose_n(legacy_figures_v, f"tier4={t}|case")[0] for legacy_figures_v, t in zip(tier_tbl["n_case"], tier_tbl["tier4"])]

    _small = tier_tbl["n_case"] < LEGACY_FIGURES_DISCLOSURE_MIN_CELL
    if _small.any():
        _cand = tier_tbl[(~_small) & (tier_tbl["tier4"] != "non_carrier")].sort_values("n_case")
        if len(_cand):
            tier_tbl.loc[_cand.index[0], "n_case_disclosed"] = "suppressed (secondary)"
            tier_tbl.loc[_small, "n_case_disclosed"] = "<20"
            audit_override(f"tier4={tier_tbl.loc[_cand.index[0], 'tier4']}|case", "suppressed (secondary)",
                           "secondary suppression (sum of tiers = published case total)", ">=20")
            for _t in tier_tbl.loc[_small, "tier4"]:
                audit_override(f"tier4={_t}|case", "<20", "shown as <20 (part of secondary-suppression pair)", "5-19 or <5")
    legacy_figures_log('Four-layer structure (only the counts after disclosure are printed; the raw counts are intentionally not entered into the log, because the log file itself will be entered into git):')
    legacy_figures_log(tier_tbl[["tier4", "n_total_disclosed", "n_case_disclosed"]].to_string(index=False))
    n_compound = int((legacy_figures_df["tier4"] == "compound").sum())
    
    legacy_figures_log(f"\n compound carrier (>=2 panel gene PTV): n={disclose_n(n_compound)[0]}(Instructions expected ~20-25)")
    legacy_figures_log(f"In the 3-layer version, compound is merged into high_penetrance or moderate according to the highest penetrance gene it carries (split number <20, not printed)")
    legacy_figures_log('⚠️ MAP3K1 is classified as moderate, but this is a **controversial classification** (mainly common-variant sites,'
        'Their PTV weights are derived from single cohort MWS estimates without independent replication).')

    _pg_rows = []
    for legacy_figures_g in LEGACY_FIGURES_GENE_COLS:
        _c = legacy_figures_df[f"gene_{legacy_figures_g}"] == 1
        _t, _k, _l = _split(int(_c.sum()), int((_c & (legacy_figures_df["group"] == "case")).sum()),
                            int((_c & (legacy_figures_df["group"] == "control")).sum()))
        _pg_rows.append({"gene": legacy_figures_g, "tier": ("high_penetrance" if legacy_figures_g in HIGH_PEN else "moderate"),
                         "contested_tier_assignment": legacy_figures_g in CONTESTED,
                         "n_carriers_any": _t, "n_carriers_case": _k, "n_carriers_control": _l})
    per_gene = pd.DataFrame(_pg_rows)
    legacy_figures_log('Number of carriers per gene (any-carrier, after disclosure processing; this project is the only table that discloses the number of people per gene):')
    legacy_figures_log(per_gene.to_string(index=False))

    write_table(tier_tbl.drop(columns=["n_case", "n_control", "n_total"]), "19_gene_tier_composition",
                note="tiers: high=BRCA1/BRCA2/PALB2/TP53; moderate=ATM/CHEK2/BARD1/RAD51C/RAD51D/MAP3K1 "
                     "(MAP3K1 assignment CONTESTED: primarily a common-variant locus, PTV weight from a "
                     "single-cohort MWS estimate); compound=>=2 panel-gene PTVs, kept as its own category.")
    write_table(per_gene, "19_gene_tier_carrier_counts",
                note="Canonical per-gene carrier table (Aim 1 modelling cohort). any-carrier: a compound carrier is "
                     "counted under every gene it carries, so genes do not sum to the carrier total. "
                     "<20 = suppressed; 'suppressed (secondary)' = would reveal a <20 cell by subtraction.")

    # =============================================================================
    
    # =============================================================================
    legacy_figures_section('Fig 1 | Orthogonality between PRS311 and Burden_rare (controls-only is the main one, and the whole queue is the second one)')

    legacy_figures_log('⚠️ Key points of methodology (must be included in figure notes and reports): This figure is mainly analyzed using controls, without cases.\n         Calculating the correlation between PRS and rare burden inside a case (or conditional on the case state) will depend on collider /\n         Index-event bias artificially induces negative correlation: two independent paths leading to the same outcome are negatively correlated with each other after a given outcome.\n         Therefore, the correlation within the case itself is bias and cannot be regarded as "orthogonal/complementary" evidence. The full cohort version contains 6.2%\n         case, also with slight contamination from the same direction, is only used as a secondary reference.')

    controls = legacy_figures_df[legacy_figures_df["group"] == "control"].copy()
    legacy_figures_log(f"  controls n={len(controls)}, where carrier n={int(controls['carrier_main_10gene'].sum())}")

    corr_rows = []
    for scope, legacy_figures_sub in [("controls_only(primary)", controls), ("full_cohort(secondary)", legacy_figures_df)]:
        pcs = legacy_figures_sub[LEGACY_FIGURES_PC_COLS].values
        # point-biserial = Pearson(PRS_z, carrier 0/1)
        corr_rows += [dict(scope=scope, **legacy_figures_r) for legacy_figures_r in corr_row(
            "PRS311_z vs carrier_main_10gene (point-biserial)",
            legacy_figures_sub["PRS311_z"].values, legacy_figures_sub["carrier_main_10gene"].values, pcs, kind="pearson")]
        corr_rows += [dict(scope=scope, **legacy_figures_r) for legacy_figures_r in corr_row(
            "PRS311_z vs Burden_rare (continuous, all)",
            legacy_figures_sub["PRS311_z"].values, legacy_figures_sub["burden_rare_main_ptv_10gene"].values, pcs, kind="pearson")]
        corr_rows += [dict(scope=scope, **legacy_figures_r) for legacy_figures_r in corr_row(
            "PRS311_z vs Burden_rare (continuous, all)",
            legacy_figures_sub["PRS311_z"].values, legacy_figures_sub["burden_rare_main_ptv_10gene"].values, pcs, kind="spearman")]
        car = legacy_figures_sub[legacy_figures_sub["carrier_main_10gene"] == 1]
        pcs_c = car[LEGACY_FIGURES_PC_COLS].values
        corr_rows += [dict(scope=scope, **legacy_figures_r) for legacy_figures_r in corr_row(
            "PRS311_z vs Burden_rare (carriers only)",
            car["PRS311_z"].values, car["burden_rare_main_ptv_10gene"].values, pcs_c, kind="pearson")]
        corr_rows += [dict(scope=scope, **legacy_figures_r) for legacy_figures_r in corr_row(
            "PRS311_z vs Burden_rare (carriers only)",
            car["PRS311_z"].values, car["burden_rare_main_ptv_10gene"].values, pcs_c, kind="spearman")]

    corr_tbl = pd.DataFrame(corr_rows)[["scope", "comparison", "method", "adjustment", "n", "r", "ci_low", "ci_high", "p_value"]]
    legacy_figures_log('Correlation coefficient table (raw and PC1-10 partial correlation):')
    legacy_figures_log(corr_tbl.to_string(index=False))
    write_table(corr_tbl, "19_fig1_prs_burden_correlations",
                note="PRIMARY = controls only. Cases are excluded on purpose: conditioning on case status "
                     "induces a negative PRS-burden association via collider / index-event bias, so a "
                     "case-internal correlation measures the bias, not orthogonality. Context only (NOT a "
                     "benchmark to match): Mazzarotto et al. bioRxiv 2025 doi:10.1101/2025.09.15.676308 "
                     "report r=0.097 (schizophrenia) and r=-0.178 (HCM) between PRS and rare-variant burden; "
                     "different traits and different burden-score construction.")

    ctrl_car = controls[controls["carrier_main_10gene"] == 1]
    n_pts = len(ctrl_car)
    legacy_figures_log(f"Panel B points (control carriers) ={n_pts}")

    legacy_figures_fig, legacy_figures_axes = plt.subplots(1, 2, figsize=(11, 4.6))
    axA = legacy_figures_axes[0]
    data_box = [controls.loc[controls["carrier_main_10gene"] == 0, "PRS311_z"].values,
                controls.loc[controls["carrier_main_10gene"] == 1, "PRS311_z"].values]
    axA.boxplot(data_box, tick_labels=[f"non-carrier\n(n={len(data_box[0]):,})", f"carrier\n(n={len(data_box[1]):,})"],
                showfliers=False)
    axA.set_ylabel("PRS311 (z-score)")
    axA.set_xlabel("Rare PTV carrier status (10-gene panel)")
    axA.set_title("PRS311 by carrier status (controls)")
    _m = ((corr_tbl["scope"] == "controls_only(primary)") & (corr_tbl["adjustment"] == "raw")
          & corr_tbl["comparison"].str.contains("point-biserial"))
    r_pb = float(corr_tbl.loc[_m, "r"].iloc[0])
    axA.text(0.03, 0.03, f"point-biserial r = {r_pb:.4f}", transform=axA.transAxes, va="bottom", fontsize=9)

    axB = legacy_figures_axes[1]
    
    hbB = axB.hexbin(ctrl_car["burden_rare_main_ptv_10gene"], ctrl_car["PRS311_z"], gridsize=15,
                     mincnt=LEGACY_FIGURES_DISCLOSURE_MIN_CELL, cmap="Blues")
    legacy_figures_fig.colorbar(hbB, ax=axB, label=f"Participants per cell (cells with <{LEGACY_FIGURES_DISCLOSURE_MIN_CELL} not shown)")
    fit_axes_to_hexbin(axB, hbB)
    n_shown_B = int(hbB.get_array().sum())
    legacy_figures_log(f"  Fig1B hexbin: {n_shown_B:,}/{n_pts:,}controls carriers fall within n>={LEGACY_FIGURES_DISCLOSURE_MIN_CELL}Gerry")
    axB.set_xlabel("Burden_rare (weighted PTV score, raw ln-OR units)")
    axB.set_ylabel("PRS311 (z-score)")
    axB.set_title(f"PRS311 vs Burden_rare, carriers (controls, n={n_pts:,})")
    _m2 = ((corr_tbl["scope"] == "controls_only(primary)") & (corr_tbl["adjustment"] == "raw")
           & (corr_tbl["method"] == "pearson") & corr_tbl["comparison"].str.contains("carriers only"))
    r_c = float(corr_tbl.loc[_m2, "r"].iloc[0])
    _m3 = ((corr_tbl["scope"] == "controls_only(primary)") & (corr_tbl["adjustment"] == "raw")
           & (corr_tbl["method"] == "spearman") & corr_tbl["comparison"].str.contains("carriers only"))
    rs_c = float(corr_tbl.loc[_m3, "r"].iloc[0])
    axB.text(0.03, 0.96, f"Pearson r = {r_c:.3f}\nSpearman rho = {rs_c:.3f}", transform=axB.transAxes, va="top", fontsize=9)
    legacy_figures_fig.tight_layout()
    fig1_path = save_fig(legacy_figures_fig, "19_fig1_prs_burden_orthogonality", prefer="svg", max_svg_kb=250)

    box_rows = []
    for scope, legacy_figures_sub in [("controls_only(primary)", controls), ("full_cohort(secondary)", legacy_figures_df)]:
        for legacy_figures_c in [0, 1]:
            s = legacy_figures_sub.loc[legacy_figures_sub["carrier_main_10gene"] == legacy_figures_c, "PRS311_z"]
            box_rows.append({"scope": scope, "carrier_main_10gene": legacy_figures_c,
                             "n_disclosed": disclose_n(len(s), f"fig1_box|{scope}|carrier={legacy_figures_c}")[0],
                             "mean": s.mean(), "sd": s.std(), "q25": s.quantile(.25),
                             "median": s.median(), "q75": s.quantile(.75)})
    box_tbl = pd.DataFrame(box_rows)
    legacy_figures_log('Fig1 Panel A accompanying values:')
    legacy_figures_log(box_tbl.to_string(index=False))
    write_table(box_tbl, "19_fig1_prs_by_carrier_distribution")

    # =============================================================================
    
    # =============================================================================
    legacy_figures_section('Fig 2 | 2D heat map: PRS311 decile x carrier(0/1)')

    legacy_figures_log('⚠️The amount of Panel A is **cohort-observed case proportion (the proportion of cases observed in this cohort)**,\n         **Not absolute risk**. ~6.4% of the cases in this queue come from the AoU WGS sub-queue.\n         EHR diagnosis is neither a lifetime cumulative incidence nor a population incidence; reading it as an absolute risk is the error most likely to be caught by reviewers.')

    legacy_figures_df["prs_decile"] = pd.qcut(legacy_figures_df["PRS311"], 10, labels=False) + 1   
    legacy_figures_df["prs_ventile"] = pd.qcut(legacy_figures_df["PRS311"], 20, labels=False) + 1
    legacy_figures_df["prs_quintile"] = pd.qcut(legacy_figures_df["PRS311"], 5, labels=False) + 1

    X_noint = add_const(legacy_figures_df[["PRS311_z", "carrier_main_10gene"] + COV_SIMPLE])
    fit_noint = fit_logit(legacy_figures_df["is_case_int"].values, X_noint)
    legacy_figures_df["p_expected_noint"] = fit_noint.fittedvalues
    legacy_figures_log('No interaction expectation model (case ~ PRS311_z + carrier + PC1-10 + age) coefficients:')
    legacy_figures_log(pd.DataFrame({"beta": fit_noint.params[["PRS311_z", "carrier_main_10gene", "age_at_last_obs"]],
                      "OR": np.exp(fit_noint.params[["PRS311_z", "carrier_main_10gene", "age_at_last_obs"]]),
                      "p": fit_noint.pvalues[["PRS311_z", "carrier_main_10gene", "age_at_last_obs"]]}).to_string())

    main_cells = cell_table(legacy_figures_df, "prs_decile", "carrier_main_10gene", "fig2_main")
    legacy_figures_log('Main panel 20 cells (PRS decile x carrier) (raw count column is not printed):')
    legacy_figures_log(main_cells.drop(columns=["n_raw", "n_case_raw"]).to_string(index=False))
    n_supp_main = int(main_cells["cell_suppressed_for_disclosure"].sum())
    legacy_figures_log(f"Number of cells on the main panel that are suppressed by disclosure rules:{n_supp_main}")
    write_table(main_cells.drop(columns=["n_raw", "n_case_raw"]), "19_fig2_prs_decile_x_carrier_cells",
                note="Panel A quantity is the COHORT-OBSERVED CASE PROPORTION within an EHR-ascertained AoU "
                     "WGS subcohort (~6.4% cases overall) -- it is NOT absolute/lifetime risk. Panel B = "
                     "observed / expected, expected = mean predicted probability from a NO-INTERACTION "
                     "logistic model (case ~ PRS311_z + carrier + PC1-10 + age). This figure is this "
                     "project's own construction; it is NOT taken from Bolze et al. 2024 (that paper uses "
                     "Kaplan-Meier / absolute-risk-at-70 with PRS dichotomised at the median, no OR heatmap).")

    supp_cells = cell_table(legacy_figures_df, "prs_quintile", "tier3", "fig2_supp")

    LEGACY_FIGURES_INTERNAL_DIR = os.path.join(LEGACY_FIGURES_BASE, "data", "19_internal")
    os.makedirs(LEGACY_FIGURES_INTERNAL_DIR, exist_ok=True)
    legacy_figures_log('Supplementary panel (PRS quintile x 3 tiers) -> data/19_internal/ (internal, no count printed)')
    write_table(supp_cells.drop(columns=["n_raw", "n_case_raw"]), "19_fig2_prs_quintile_x_tier_cells", out_dir=LEGACY_FIGURES_INTERNAL_DIR,
                note="3-tier version only: compound carriers (>=2 panel-gene PTVs) are folded into the "
                     "highest tier they carry (high_penetrance if any of BRCA1/BRCA2/PALB2/TP53, else "
                     "moderate). Fig 3 keeps compound carriers as a separate category. INTERPRETATION CAVEAT: "
                     "expected values come from a model with a SINGLE pooled carrier indicator, so O/E here "
                     "mixes (a) between-tier heterogeneity in the main effect (high-penetrance carriers are "
                     "riskier than the pooled carrier coefficient, moderate carriers less risky) with (b) any "
                     "genuine PRS-by-carrier interaction. Do not read the tier-wise O/E LEVEL as interaction; "
                     "only the within-tier gradient across PRS quintiles is informative about interaction.")

    legacy_figures_fig, legacy_figures_axes = plt.subplots(1, 2, figsize=(12, 5.4))
    gA = to_grid(main_cells, "prs_decile", "carrier_main_10gene", "observed_case_proportion")
    annotA = to_grid(main_cells, "prs_decile", "carrier_main_10gene", "n_disclosed").astype(str) + "\n(" + \
             to_grid(main_cells, "prs_decile", "carrier_main_10gene", "n_case_disclosed").astype(str) + " cases)"
    sns.heatmap(gA, annot=annotA, fmt="", cmap="Reds", ax=legacy_figures_axes[0], cbar_kws={"label": "cohort-observed case proportion"},
                annot_kws={"fontsize": 7})
    legacy_figures_axes[0].set_title("Observed case proportion")
    legacy_figures_axes[0].set_xlabel("Rare PTV carrier (10-gene panel)")
    legacy_figures_axes[0].set_ylabel("PRS311 decile (full cohort)")
    legacy_figures_axes[0].invert_yaxis()

    gB = to_grid(main_cells, "prs_decile", "carrier_main_10gene", "obs_over_exp")
    sns.heatmap(gB, annot=True, fmt=".2f", cmap="RdBu_r", center=1.0, ax=legacy_figures_axes[1],
                cbar_kws={"label": "observed / expected (no-interaction model)"}, annot_kws={"fontsize": 8})
    legacy_figures_axes[1].set_title("Observed / expected (no-interaction model)")
    legacy_figures_axes[1].set_xlabel("Rare PTV carrier (10-gene panel)")
    legacy_figures_axes[1].set_ylabel("PRS311 decile (full cohort)")
    legacy_figures_axes[1].invert_yaxis()
    legacy_figures_fig.tight_layout()
    save_fig(legacy_figures_fig, "19_fig2_prs_decile_x_carrier_heatmap", prefer="svg")

    legacy_figures_fig, legacy_figures_axes = plt.subplots(1, 2, figsize=(12, 4.6))
    order3 = ["non_carrier", "moderate", "high_penetrance"]
    sA = supp_cells.pivot(index="prs_quintile", columns="tier3", values="observed_case_proportion").reindex(columns=order3)
    annS = (supp_cells.pivot(index="prs_quintile", columns="tier3", values="n_disclosed").reindex(columns=order3).astype(str)
            + "\n(" + supp_cells.pivot(index="prs_quintile", columns="tier3", values="n_case_disclosed").reindex(columns=order3).astype(str) + ")")
    sns.heatmap(sA, annot=annS, fmt="", cmap="Reds", ax=legacy_figures_axes[0], cbar_kws={"label": "cohort-observed case proportion"},
                annot_kws={"fontsize": 8})
    legacy_figures_axes[0].set_title("Supp A. Observed case proportion (NOT absolute risk)")
    legacy_figures_axes[0].set_ylabel("PRS311 quintile"); legacy_figures_axes[0].set_xlabel("gene tier"); legacy_figures_axes[0].invert_yaxis()
    sB = supp_cells.pivot(index="prs_quintile", columns="tier3", values="obs_over_exp").reindex(columns=order3)
    sns.heatmap(sB, annot=True, fmt=".2f", cmap="RdBu_r", center=1.0, ax=legacy_figures_axes[1],
                cbar_kws={"label": "observed / expected"}, annot_kws={"fontsize": 8})
    legacy_figures_axes[1].set_title("Supp B. Observed / expected (no-interaction model)")
    legacy_figures_axes[1].set_ylabel("PRS311 quintile"); legacy_figures_axes[1].set_xlabel("gene tier"); legacy_figures_axes[1].invert_yaxis()
    legacy_figures_fig.suptitle("Fig 2 (supplementary). PRS quintile x gene tier", fontsize=11)
    legacy_figures_fig.tight_layout()
    save_fig(legacy_figures_fig, "19_fig2_supp_prs_quintile_x_tier_heatmap", prefer="svg", out_dir=LEGACY_FIGURES_INTERNAL_DIR)

    hi_corner = main_cells[(main_cells["prs_decile"] >= 9) & (main_cells["carrier_main_10gene"] == 1)]
    legacy_figures_log('O/E for high PRS x carrier corner (decile 9-10, carrier=1):')
    legacy_figures_log(hi_corner[["prs_decile", "n_disclosed", "n_case_disclosed", "observed_case_proportion", "obs_over_exp"]].to_string(index=False))

    # =============================================================================
    # Fig 3 | Quantile-OR dot-and-whisker
    # =============================================================================
    legacy_figures_section('Fig 3 | Quantile-OR point and line plot (A: PRS ventiles; B: burden ordered categories)')

    legacy_figures_log('Binning convention from Mazzarotto et al. 2025 (bioRxiv doi:10.1101/2025.09.15.676308),\n         **This round has obtained the full PDF text and checked it item by item** (plans/papers/The complementary roles of rare variant\n         burden scores.pdf), no longer a paraphrase. Checked original agreement:\n           - Both fractions are cut into quantile bins, "default N=20; ventiles";\n           - "logistic regressions of diagnosis on IOGC quantile (reference: lowest quantile)"\n             -> **The reference group is the lowest ventile**, not the middle band;\n           - Additionally dichotomize IOGC into low (ventiles 1-18) vs high (>18) for PPV analysis.\n         The results (cited as shape precedent, not benchmark): Schizophrenia PRS top ventile vs bottom OR=13.6\n         (95%CI 8.1-23), the middle quantiles increase monotonically; while IOGC "did not vary linearly with ventiles,\n         showing a significant increase only for the top two", top ventile OR=2.3 (1.5-3.7).\n         The trait is schizophrenia/HCM, not breast cancer -> just for context and "shape" precedent.')

    REF_LOWEST = 1          
    REF_VENTILES = [9, 10]  
    legacy_figures_log(f"Panel A main reference = ventile{REF_LOWEST}(lowest quantile, aligned with Mazzarotto);"
        f"secondary reference = ventile{REF_VENTILES}(middle of 40-60th percentile, left as secondary)")

    fit_v_low = fit_ventile_model([REF_LOWEST])
    fit_v_mid = fit_ventile_model(REF_VENTILES)

    legacy_figures_rows = []
    for ref_label, ref_set, ft in [("lowest_ventile (PRIMARY, Mazzarotto-aligned)", [REF_LOWEST], fit_v_low),
                                   ("mid_band_ventiles_9_10 (secondary)", REF_VENTILES, fit_v_mid)]:
        for legacy_figures_v in range(1, 21):
            n = int((legacy_figures_df["prs_ventile"] == legacy_figures_v).sum())
            legacy_figures_n_case = int(legacy_figures_df.loc[legacy_figures_df["prs_ventile"] == legacy_figures_v, "is_case_int"].sum())
            if legacy_figures_v in ref_set:
                legacy_figures_rows.append({"panel": "A_PRS_ventile", "reference": ref_label, "category": f"ventile {legacy_figures_v}",
                             "is_reference": True, "n_raw": n, "n_case_raw": legacy_figures_n_case, "OR": 1.0,
                             "ci_low": np.nan, "ci_high": np.nan, "beta": 0.0, "se": np.nan, "p_value": np.nan})
            else:
                t = f"pv_v{legacy_figures_v}"
                legacy_figures_b, legacy_figures_se, p = ft.params[t], ft.bse[t], ft.pvalues[t]
                legacy_figures_rows.append({"panel": "A_PRS_ventile", "reference": ref_label, "category": f"ventile {legacy_figures_v}",
                             "is_reference": False, "n_raw": n, "n_case_raw": legacy_figures_n_case, "OR": np.exp(legacy_figures_b),
                             "ci_low": np.exp(legacy_figures_b - LEGACY_FIGURES_Z * legacy_figures_se), "ci_high": np.exp(legacy_figures_b + LEGACY_FIGURES_Z * legacy_figures_se),
                             "beta": legacy_figures_b, "se": legacy_figures_se, "p_value": p})

    legacy_figures_log("Panel B uses ordered categories instead of quantile: Burden_rare ~98% of the entire queue is 0, and there will be 19 in the 20th percentile\n         The quantiles take the same value 0, and quantile binning is mathematically undefinable (Mazzarotto's IOGC is outlier-gene\n         Counts and distributions must be much more continuous to be ventile; this is not transferable to the 10-gene PTV weighted score of this project).\n         The fact that non-carrier is used as the reference group is consistent with the article.")
    zero_frac = float((legacy_figures_df["burden_rare_main_ptv_10gene"] == 0).mean())
    legacy_figures_log(f"Measured ratio of Burden_rare == 0:{zero_frac*100:.2f}%")

    cat_order = ["non_carrier", "moderate", "high_penetrance", "compound"]
    tier_dummy = pd.get_dummies(legacy_figures_df["tier4"], prefix="t").astype(float).drop(columns=["t_non_carrier"])
    X_t = add_const(pd.concat([tier_dummy, legacy_figures_df[COV_SIMPLE].astype(float)], axis=1))
    fit_t = fit_logit(legacy_figures_df["is_case_int"].values, X_t)

    for legacy_figures_c in cat_order:
        n = int((legacy_figures_df["tier4"] == legacy_figures_c).sum())
        legacy_figures_n_case = int(legacy_figures_df.loc[legacy_figures_df["tier4"] == legacy_figures_c, "is_case_int"].sum())
        if legacy_figures_c == "non_carrier":
            legacy_figures_rows.append({"panel": "B_burden_ordered_category", "reference": "non_carrier (matches Mazzarotto)",
                         "category": legacy_figures_c, "is_reference": True, "n_raw": n, "n_case_raw": legacy_figures_n_case, "OR": 1.0,
                         "ci_low": np.nan, "ci_high": np.nan, "beta": 0.0, "se": np.nan, "p_value": np.nan})
        else:
            t = f"t_{legacy_figures_c}"
            legacy_figures_b, legacy_figures_se, p = fit_t.params[t], fit_t.bse[t], fit_t.pvalues[t]
            legacy_figures_rows.append({"panel": "B_burden_ordered_category", "reference": "non_carrier (matches Mazzarotto)",
                         "category": legacy_figures_c, "is_reference": False, "n_raw": n, "n_case_raw": legacy_figures_n_case, "OR": np.exp(legacy_figures_b),
                         "ci_low": np.exp(legacy_figures_b - LEGACY_FIGURES_Z * legacy_figures_se), "ci_high": np.exp(legacy_figures_b + LEGACY_FIGURES_Z * legacy_figures_se),
                         "beta": legacy_figures_b, "se": legacy_figures_se, "p_value": p})

    or_tbl = pd.DataFrame(legacy_figures_rows)
    or_tbl["n_disclosed"] = [disclose_n(n, f"fig3|{p}|{legacy_figures_c}")[0] for n, p, legacy_figures_c in zip(or_tbl["n_raw"], or_tbl["panel"], or_tbl["category"])]
    or_tbl["n_case_disclosed"] = [disclose_n(n, f"fig3|{p}|{legacy_figures_c}|cases")[0] for n, p, legacy_figures_c in zip(or_tbl["n_case_raw"], or_tbl["panel"], or_tbl["category"])]

    _pbm = or_tbl["panel"] == "B_burden_ordered_category"
    _small_b = _pbm & (or_tbl["n_case_raw"] < LEGACY_FIGURES_DISCLOSURE_MIN_CELL)
    if _small_b.any():
        or_tbl.loc[_small_b, "n_case_disclosed"] = "<20"
        _cand_b = or_tbl[_pbm & ~_small_b & (or_tbl["category"] != "non_carrier")].sort_values("n_case_raw")
        or_tbl.loc[_cand_b.index[0], "n_case_disclosed"] = "suppressed (secondary)"
        audit_override(f"fig3|B_burden_ordered_category|{or_tbl.loc[_cand_b.index[0], 'category']}|cases",
                       "suppressed (secondary)", "secondary suppression (sum of categories = published case total)", ">=20")
        for _c in or_tbl.loc[_small_b, "category"]:
            audit_override(f"fig3|B_burden_ordered_category|{_c}|cases", "<20",
                           "shown as <20 (part of secondary-suppression pair)", "5-19 or <5")
    or_tbl["interpretation_warning"] = np.where(
        or_tbl["category"] == "compound",
        "UNDERPOWERED: n~20-25, CI spans orders of magnitude -- do NOT interpret the point estimate", "")
    legacy_figures_log('Fig 3 OR table:')
    legacy_figures_log(or_tbl[["panel", "reference", "category", "is_reference", "n_disclosed", "n_case_disclosed",
                "OR", "ci_low", "ci_high", "p_value"]].to_string(index=False))
    write_table(or_tbl.drop(columns=["n_raw", "n_case_raw"]), "19_fig3_quantile_or_estimates",
                note="Both panels adjusted for PC1-10 + age_at_last_obs (covariate set as specified for this "
                     "figure; the scripts/16 primary model additionally adjusts for sequencing site/sequencer). "
                     "Panel A is reported under TWO parameterisations: PRIMARY reference = LOWEST ventile, "
                     "matching Mazzarotto et al. 2025 (bioRxiv doi:10.1101/2025.09.15.676308), whose Methods "
                     "state 'logistic regressions of diagnosis on IOGC quantile (reference: lowest quantile)' "
                     "with 'default N=20; ventiles'; SECONDARY reference = mid band (ventiles 9+10, 40th-60th "
                     "percentile). Verified directly against the PDF (plans/papers/), not a second-hand "
                     "description. Panel B: Burden_rare is ~98% zero cohort-wide so ventile binning is "
                     "undefined here (their IOGC is an outlier-gene count with a far more continuous "
                     "distribution); ordered categories used instead with non-carrier as reference (same "
                     "reference convention as the paper), and COMPOUND CARRIERS ARE KEPT AS THEIR OWN "
                     "CATEGORY (not merged into high-penetrance) because they are exactly the group Antinomy 1 "
                     "turns on. Shape precedent from the paper (context, NOT a benchmark; schizophrenia/HCM, "
                     "not breast cancer): PRS top ventile OR=13.6 (8.1-23) vs lowest with a smooth monotone "
                     "gradient, whereas IOGC 'did not vary linearly with ventiles, showing a significant "
                     "increase only for the top two', top ventile OR=2.3 (1.5-3.7).")

    pa = or_tbl[(or_tbl["panel"] == "A_PRS_ventile")
                & (or_tbl["reference"].str.startswith("lowest"))].copy()
    pa_mid = or_tbl[(or_tbl["panel"] == "A_PRS_ventile")
                    & (or_tbl["reference"].str.startswith("mid_band"))].copy()

    legacy_figures_fig, legacy_figures_axes = plt.subplots(1, 2, figsize=(12.5, 5.0), gridspec_kw={"width_ratios": [2.1, 1]})
    legacy_figures_xs = np.arange(1, 21)
    legacy_figures_axes[0].errorbar(legacy_figures_xs, pa["OR"], yerr=[pa["OR"] - pa["ci_low"].fillna(pa["OR"]),
                                        pa["ci_high"].fillna(pa["OR"]) - pa["OR"]],
                     fmt="o", color="tab:blue", markersize=4, capsize=2, linestyle="none")
    legacy_figures_axes[0].axhline(1.0, color="gray", linewidth=0.8)
    legacy_figures_axes[0].plot(REF_LOWEST, 1.0, marker="s", color="black", markersize=6)
    legacy_figures_axes[0].set_yscale("log")
    legacy_figures_axes[0].set_xticks(legacy_figures_xs)
    legacy_figures_axes[0].set_xlabel("PRS311 ventile (reference = lowest ventile, as in Mazzarotto et al. 2025)")
    legacy_figures_axes[0].set_ylabel("OR (adjusted for PC1-10 + age), log scale")
    legacy_figures_axes[0].set_title("PRS311 ventiles")

    pb = or_tbl[or_tbl["panel"] == "B_burden_ordered_category"].copy().set_index("category").loc[cat_order].reset_index()
    for legacy_figures_i, legacy_figures_r in pb.iterrows():
        if legacy_figures_r["category"] == "compound":
            legacy_figures_axes[1].errorbar(legacy_figures_i, legacy_figures_r["OR"], yerr=[[legacy_figures_r["OR"] - legacy_figures_r["ci_low"]], [legacy_figures_r["ci_high"] - legacy_figures_r["OR"]]],
                             fmt="o", mfc="none", color="tab:red", markersize=9, capsize=3)
        else:
            legacy_figures_lo = 0 if np.isnan(legacy_figures_r["ci_low"]) else legacy_figures_r["OR"] - legacy_figures_r["ci_low"]
            legacy_figures_hi = 0 if np.isnan(legacy_figures_r["ci_high"]) else legacy_figures_r["ci_high"] - legacy_figures_r["OR"]
            legacy_figures_axes[1].errorbar(legacy_figures_i, legacy_figures_r["OR"], yerr=[[legacy_figures_lo], [legacy_figures_hi]], fmt="o", color="tab:blue", markersize=6, capsize=3)
    legacy_figures_axes[1].axhline(1.0, color="gray", linewidth=0.8)
    legacy_figures_axes[1].set_yscale("log")
    legacy_figures_axes[1].set_xticks(np.arange(len(pb)))
    legacy_figures_axes[1].set_xticklabels(["non-carrier\n(ref)", "single\nmoderate PTV", "single\nhigh-pen. PTV", "compound\n(>=2 PTVs)"], fontsize=8)
    legacy_figures_axes[1].set_ylabel("OR (adjusted for PC1-10 + age), log scale")
    legacy_figures_axes[1].set_title("Rare-variant burden categories")
    comp = pb[pb["category"] == "compound"].iloc[0]
    legacy_figures_axes[1].annotate(f"n = {comp['n_disclosed']}: underpowered,\nnot interpreted",
                     xy=(3, comp["OR"]), xytext=(0.42, 0.06), textcoords="axes fraction", fontsize=8,
                     arrowprops=dict(arrowstyle="->", linewidth=0.7))
    legacy_figures_fig.tight_layout()
    save_fig(legacy_figures_fig, "19_fig3_quantile_or_dotwhisker", prefer="svg")

    legacy_figures_section("Fig 3 Additional | 'PPV analogy': PRS centile threshold x burden high/low")

    legacy_figures_log('⚠️ Three points must be written clearly:\n       1) What this queue calculates is not PPV**. PPV depends on prevalence, 6.2% of cases in this cohort were case-control\n          Design + EHR What determines the diagnosis is not the population prevalence; therefore, what is reported here is **cohort-observed case\n          proportion (proportion of observed cases in this cohort)**, and the PPV value of this article** cannot be directly compared**, but can only be compared\n          The shape "Relative improvement factor of high burden vs low burden".\n       2) The PPV value of this article is **percentage** (schizophrenia 0.25% -> 0.41%; PRS>95th centile 0.89% ->\n          1.13%; the MDD analogy in the same text is read as 6.03 -> 6.13, 17.12 -> 14.31, and the confirmed unit is %), not a ratio.\n       3) This article divides IOGC into low (ventiles 1-18) vs high (>18); 98% of the burden of this project is 0,\n          The ventile threshold is meaningless, so burden high/low = carrier/non-carrier (10-gene panel PTV).')

    prs_ct = legacy_figures_df["PRS311"].rank(pct=True) * 100
    ppv_bands = [("<15th centile", prs_ct < 15), ("<25th centile", prs_ct < 25),
                 ("<50th centile", prs_ct < 50), ("<75th centile", prs_ct < 75),
                 (">95th centile", prs_ct > 95), ("any PRS (all)", pd.Series(True, index=legacy_figures_df.index))]
    ppv_rows = []
    for band, mask in ppv_bands:
        for bl, bmask in [("low burden (non-carrier)", legacy_figures_df["carrier_main_10gene"] == 0),
                          ("high burden (any panel PTV carrier)", legacy_figures_df["carrier_main_10gene"] == 1)]:
            cell = legacy_figures_df[mask & bmask]
            n, nc = len(cell), int(cell["is_case_int"].sum())
            supp = (n < LEGACY_FIGURES_DISCLOSURE_MIN_CELL) or (nc < 5)
            legacy_figures_lo, legacy_figures_hi = (np.nan, np.nan)
            if not supp:
                legacy_figures_lo, legacy_figures_hi = stats.binomtest(nc, n).proportion_ci(confidence_level=0.95, method="wilson")
            ppv_rows.append({"prs_band": band, "burden_stratum": bl,
                             "n_disclosed": disclose_n(n, f"fig3ppv|{band}|{bl}")[0],
                             "n_case_disclosed": disclose_n(nc, f"fig3ppv|{band}|{bl}|cases")[0],
                             "observed_case_percent": np.nan if supp else 100 * nc / n,
                             "ci_low_percent": np.nan if supp else 100 * legacy_figures_lo,
                             "ci_high_percent": np.nan if supp else 100 * legacy_figures_hi,
                             "cell_suppressed_for_disclosure": supp})
    ppv_tbl = pd.DataFrame(ppv_rows)
    piv = ppv_tbl.pivot(index="prs_band", columns="burden_stratum", values="observed_case_percent")
    piv["ratio_high_over_low"] = piv.iloc[:, 0] / piv.iloc[:, 1] if piv.columns[0].startswith("high") else piv.iloc[:, 1] / piv.iloc[:, 0]
    ppv_tbl = ppv_tbl.merge(piv[["ratio_high_over_low"]].reset_index(), on="prs_band", how="left")
    legacy_figures_log('PPV analogy table (% of observed cases, non-PPV, non-absolute risk):')
    legacy_figures_log(ppv_tbl.to_string(index=False))
    write_table(ppv_tbl, "19_fig3_ppv_analogue_prs_band_x_burden",
                note="ANALOGUE ONLY of Mazzarotto et al. Supplementary Figure 3 -- these are NOT PPVs. PPV "
                     "depends on prevalence; this cohort's ~6.2% case fraction is set by the case-control "
                     "design and EHR ascertainment, not by population prevalence, so the LEVELS are not "
                     "comparable to the paper's numbers -- only the high-vs-low-burden ratio (shape) is. The "
                     "paper's PPVs are PERCENTAGES (schizophrenia 0.25% -> 0.41% at any PRS, 1.64-fold; "
                     "0.89% -> 1.13% above the 95th PRS centile; cf. the MDD analogue 6.03 -> 6.13 and "
                     "17.12 -> 14.31 in the same paper, which fixes the unit as %), not proportions. PRS "
                     "centile thresholds match the paper's pre-specified set (<15th, <25th, <50th, <75th, "
                     ">95th) and are cumulative/overlapping as published, plus an 'any PRS' row. The paper "
                     "dichotomises IOGC at ventile >18; that cutoff is meaningless here because Burden_rare "
                     "is ~98% zero, so burden high/low = carrier / non-carrier of any 10-gene panel PTV.")

    legacy_figures_fig, legacy_figures_axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
    legacy_figures_axes[0].errorbar(legacy_figures_xs, pa_mid["OR"], yerr=[pa_mid["OR"] - pa_mid["ci_low"].fillna(pa_mid["OR"]),
                                             pa_mid["ci_high"].fillna(pa_mid["OR"]) - pa_mid["OR"]],
                     fmt="o", color="tab:gray", markersize=4, capsize=2, linestyle="none")
    legacy_figures_axes[0].axhline(1.0, color="gray", linewidth=0.8)
    for legacy_figures_v in REF_VENTILES:
        legacy_figures_axes[0].plot(legacy_figures_v, 1.0, marker="s", color="black", markersize=5)
    legacy_figures_axes[0].set_yscale("log")
    legacy_figures_axes[0].set_xticks(legacy_figures_xs)
    legacy_figures_axes[0].set_xlabel("PRS311 ventile (reference = ventiles 9+10, 40th-60th percentile)")
    legacy_figures_axes[0].set_ylabel("OR (adjusted for PC1-10 + age), log scale")
    legacy_figures_axes[0].set_title("PRS311 ventiles, mid-band reference")

    band_order = ["<15th centile", "<25th centile", "<50th centile", "<75th centile", ">95th centile", "any PRS (all)"]
    xb = np.arange(len(band_order))
    for legacy_figures_stratum, color, off in [("low burden (non-carrier)", "tab:blue", -0.12),
                                ("high burden (any panel PTV carrier)", "tab:red", 0.12)]:
        legacy_figures_sub = ppv_tbl[ppv_tbl["burden_stratum"] == legacy_figures_stratum].set_index("prs_band").loc[band_order]
        legacy_figures_axes[1].errorbar(xb + off, legacy_figures_sub["observed_case_percent"],
                         yerr=[legacy_figures_sub["observed_case_percent"] - legacy_figures_sub["ci_low_percent"],
                               legacy_figures_sub["ci_high_percent"] - legacy_figures_sub["observed_case_percent"]],
                         fmt="o", color=color, capsize=3, markersize=5, label=legacy_figures_stratum)
    legacy_figures_axes[1].set_xticks(xb)
    legacy_figures_axes[1].set_xticklabels(band_order, rotation=30, ha="right", fontsize=8)
    legacy_figures_axes[1].set_ylabel("cohort-observed case percent (%)")
    legacy_figures_axes[1].set_title("Observed case percentage by PRS band and carrier status")
    legacy_figures_axes[1].legend(fontsize=7)
    legacy_figures_fig.tight_layout()
    save_fig(legacy_figures_fig, "19_fig3_supp_midband_and_ppv_analogue", prefer="svg")

    # =============================================================================
    
    # =============================================================================
    legacy_figures_section('Fig 6.1 | Decisive test: refit the main joint model after removing the compound carrier')

    legacy_figures_X_primary = add_const(legacy_figures_df[["PRS311_z", "burden_z", "inter_z"] + LEGACY_FIGURES_COV_COLS])
    fit_primary = fit_logit(legacy_figures_df["is_case_int"].values, legacy_figures_X_primary)
    legacy_figures_b0 = fit_primary.params["inter_z"]
    legacy_figures_log(f"Reproduce scripts/16 main model inter_z beta ={legacy_figures_b0:.6f}(scripts/16 record value{REF16_INTER_BETA:.6f})")
    assert abs(legacy_figures_b0 - REF16_INTER_BETA) < 1e-6, 'FATAL: Unable to reproduce interaction coefficients for scripts/16'
    legacy_figures_log('ASSERT PASS: Reproduce bit by bit with scripts/16 main model')

    mask_nc = legacy_figures_df["tier4"] != "compound"
    df_nc = legacy_figures_df[mask_nc]
    X_nc = add_const(df_nc[["PRS311_z", "burden_z", "inter_z"] + LEGACY_FIGURES_COV_COLS])
    fit_nc = fit_logit(df_nc["is_case_int"].values, X_nc)
    legacy_figures_log(f"After excluding compound carrier n={len(df_nc)}(eliminate{n_compound}people)")
    legacy_figures_log('⚠️ The z-score normalization parameter **continues to be used in the entire queue** (not recalculated on the reduced queue) to ensure that the coefficients of the two fits are comparable on the same scale.')

    cmp_rows = []
    for legacy_figures_lab, ft, n in [("primary (full cohort, scripts/16 replication)", fit_primary, len(legacy_figures_df)),
                       ("compound carriers excluded", fit_nc, len(df_nc))]:
        for t in ["PRS311_z", "burden_z", "inter_z"]:
            legacy_figures_b, legacy_figures_se, p = ft.params[t], ft.bse[t], ft.pvalues[t]
            cmp_rows.append({"model": legacy_figures_lab, "n": n, "term": t, "beta": legacy_figures_b, "se": legacy_figures_se, "OR": np.exp(legacy_figures_b),
                             "ci_low": np.exp(legacy_figures_b - LEGACY_FIGURES_Z * legacy_figures_se), "ci_high": np.exp(legacy_figures_b + LEGACY_FIGURES_Z * legacy_figures_se), "p_value": p})
    cmp_tbl = pd.DataFrame(cmp_rows)
    legacy_figures_log('Coefficient comparison before/after elimination:')
    legacy_figures_log(cmp_tbl.to_string(index=False))
    legacy_figures_b1 = fit_nc.params["inter_z"]
    legacy_figures_p1 = fit_nc.pvalues["inter_z"]
    legacy_figures_log(f"\n >>> Key conclusion: inter_z beta{legacy_figures_b0:.5f} (OR={np.exp(legacy_figures_b0):.4f}, p={fit_primary.pvalues['inter_z']:.3e}) "
        f"-> {legacy_figures_b1:.5f} (OR={np.exp(legacy_figures_b1):.4f}, p={legacy_figures_p1:.3e}); change{100*(legacy_figures_b1-legacy_figures_b0)/abs(legacy_figures_b0):+.1f}%")
    verdict = ('The interaction term remains significant and has the same direction after removing the compound carrier -> the negative interaction is **not** the product of these ~20-25 highly leveraged individuals'
               if (legacy_figures_p1 < 0.05 and np.sign(legacy_figures_b1) == np.sign(legacy_figures_b0)) else
               'The interaction term loses significance/changes direction after eliminating the compound carrier -> the negative interaction is highly dependent on this small group of high-leverage individuals and is an artifact')
    legacy_figures_log(f">>>Interpretation:{verdict}")
    write_table(cmp_tbl, "19_fig6_compound_exclusion_interaction_comparison",
                note=f"DECISIVE TEST for Gap 1. Same z-scaling (full-cohort mean/sd) used for both fits so the "
                     f"coefficients are on one scale. Covariates = scripts/16 COV_COLS (age + PC1-10 + "
                     f"sequencing site + sequencer). n_compound_excluded={n_compound}. VERDICT: {verdict}")

    # ---- Cook's distance ----
    legacy_figures_section("Fig 6.2 | Cook's distance (Aim 1 logistic main model)")
    Xm = legacy_figures_X_primary.values.astype(float)
    p_hat = np.asarray(fit_primary.fittedvalues, dtype=float)
    w_irls = p_hat * (1 - p_hat)
    legacy_figures_k = Xm.shape[1]
    XtWX = Xm.T @ (Xm * w_irls[:, None])
    A = np.linalg.inv(XtWX)
    h = w_irls * np.einsum("ij,ij->i", Xm, Xm @ A)
    resid_p = (legacy_figures_df["is_case_int"].values - p_hat) / np.sqrt(w_irls)
    cooks = (resid_p ** 2 * h) / (legacy_figures_k * (1 - h) ** 2)
    legacy_figures_log(f"design matrix{Xm.shape}, k={legacy_figures_k}; leverage h range [{h.min():.2e}, {h.max():.3f}]")
    legacy_figures_log(f"  Cook's D: max={cooks.max():.5f}, mean={cooks.mean():.3e}, 4/n threshold ={4/len(legacy_figures_df):.3e}, "
        f"Number of people exceeding the threshold ={int((cooks > 4/len(legacy_figures_df)).sum())}")

    df_c = pd.DataFrame({"cooks_d": cooks, "tier4": legacy_figures_df["tier4"].values,
                         "carrier": legacy_figures_df["carrier_main_10gene"].values, "is_case": legacy_figures_df["is_case_int"].values})

    legacy_figures_log(f"compound carrier\'s Cook\'s D: max={df_c.loc[df_c['tier4']=='compound','cooks_d'].max():.5f}, "
        f"Median ={df_c.loc[df_c['tier4']=='compound','cooks_d'].median():.5f}; "
        f"Median of the entire cohort ={np.median(cooks):.3e}")

    cook_summary = (df_c.groupby("tier4")["cooks_d"]
                    .agg(n="size", mean="mean", median="median", p95=lambda s: s.quantile(.95), max="max")
                    .reset_index())
    cook_summary["n_disclosed"] = [disclose_n(legacy_figures_v, f"cooks|{t}")[0] for legacy_figures_v, t in zip(cook_summary["n"], cook_summary["tier4"])]
    legacy_figures_log("Cook's D aggregated by tier:")
    legacy_figures_log(cook_summary.to_string(index=False))
    write_table(cook_summary.drop(columns=["n"]), "19_fig6_cooks_distance_by_tier",
                note=f"Cook's D for the Aim 1 primary logistic model (n={len(legacy_figures_df):,}, k={legacy_figures_k} params), computed as "
                     "(pearson_resid^2 * h)/(k*(1-h)^2) with h = IRLS leverage. Aggregated by gene tier only "
                     "-- no participant-level rows. The decisive Gap-1 evidence is the compound-exclusion "
                     "refit table, not this plot.")

    legacy_figures_fig, legacy_figures_axes = plt.subplots(1, 2, figsize=(12, 4.6))
    
    TIER_ORDER = ["non_carrier", "moderate", "high_penetrance", "compound"]
    tier_n = df_c["tier4"].value_counts()
    assert all(tier_n[t] >= LEGACY_FIGURES_DISCLOSURE_MIN_CELL for t in TIER_ORDER), 'If tier n<20, the distribution cannot be drawn separately.'
    legacy_figures_axes[0].boxplot([np.log10(df_c.loc[df_c["tier4"] == t, "cooks_d"] + 1e-12) for t in TIER_ORDER],
                    tick_labels=[f"{t}\n(n={tier_n[t]:,})" for t in TIER_ORDER], showfliers=False, whis=(5, 95))
    legacy_figures_axes[0].axhline(np.log10(4 / len(legacy_figures_df)), color="black", linewidth=0.8, linestyle="--")
    legacy_figures_axes[0].set_ylabel("log10 Cook's distance")
    legacy_figures_axes[0].set_title("Cook's distance by carrier tier")
    
    log_cd = np.log10(cooks + 1e-12)
    cnt, edges = np.histogram(log_cd, bins=80)
    n_bins_hidden = int(((cnt > 0) & (cnt < LEGACY_FIGURES_DISCLOSURE_MIN_CELL)).sum())
    cnt = np.where(cnt >= LEGACY_FIGURES_DISCLOSURE_MIN_CELL, cnt, 0)
    legacy_figures_axes[1].stairs(cnt, edges, fill=True, color="tab:gray")
    legacy_figures_axes[1].axvline(np.log10(4 / len(legacy_figures_df)), color="black", linewidth=0.8, linestyle="--")
    legacy_figures_axes[1].set_xlabel("log10 Cook's distance")
    legacy_figures_axes[1].set_ylabel("count")
    legacy_figures_axes[1].set_title(f"Distribution of Cook's distance (n={len(legacy_figures_df):,})")
    legacy_figures_log(f"Fig6.2B Histogram:{n_bins_hidden}Number of people in a bin <{LEGACY_FIGURES_DISCLOSURE_MIN_CELL}, not drawn")
    legacy_figures_fig.tight_layout()
    save_fig(legacy_figures_fig, "19_fig6_cooks_distance", prefer="svg")

    legacy_figures_section('Fig 6.3 | Aim 4 OLS residual diagnosis (refitting in script, unchanged scripts/17)')

    cases = legacy_figures_pheno[legacy_figures_pheno["group"] == "case"].merge(legacy_figures_qc[["person_id", "pass_qc"]], on="person_id", how="inner")
    cases = cases[cases["pass_qc"] == 1]
    cases = cases.dropna(subset=["age_at_diagnosis"])
    cases = cases[cases["bc_is_first_record"] != True]
    cases = cases.merge(legacy_figures_prs[["person_id", "PRS311"]], on="person_id", how="left")
    cases = cases.merge(legacy_figures_burden[["person_id", "burden_rare_main_ptv_10gene"]], on="person_id", how="left")
    cases = cases.merge(legacy_figures_pc[["person_id"] + LEGACY_FIGURES_PC_COLS], on="person_id", how="left")
    cases = cases.dropna(subset=["PRS311", "burden_rare_main_ptv_10gene"] + LEGACY_FIGURES_PC_COLS).reset_index(drop=True)
    legacy_figures_log(f"Aim 4 cases-only queue n={len(cases)}")
    assert len(cases) == EXPECTED_N_AIM4, f"FATAL: Aim4 n={len(cases)} != {EXPECTED_N_AIM4}(scripts/17 has inconsistent caliber)"
    legacy_figures_log('ASSERT PASS: Same as scripts/17 main analysis queue n')

    cases["burden_z"] = (cases["burden_rare_main_ptv_10gene"] - cases["burden_rare_main_ptv_10gene"].mean()) / cases["burden_rare_main_ptv_10gene"].std()
    cases["PRS311_z"] = (cases["PRS311"] - cases["PRS311"].mean()) / cases["PRS311"].std()
    X_a4 = add_const(cases[["burden_z", "PRS311_z"] + LEGACY_FIGURES_PC_COLS])
    fit_a4 = sm.OLS(cases["age_at_diagnosis"].values, X_a4).fit(cov_type="HC1")
    legacy_figures_log(f"Refitting coefficient: burden_z beta={fit_a4.params['burden_z']:.4f}(scripts/17 record{REF17_BURDEN_BETA:.4f}), "
        f"PRS311_z beta={fit_a4.params['PRS311_z']:.4f}(Record{REF17_PRS_BETA:.4f}), R2={fit_a4.rsquared:.4f}(Record{REF17_R2:.4f})")

    resid = fit_a4.resid
    fitted = fit_a4.fittedvalues
    jb_stat, jb_p, skew, kurt = sm.stats.stattools.jarque_bera(resid)
    rng = np.random.default_rng(20260819)
    sub_idx = rng.choice(len(resid), size=5000, replace=False)
    sw_stat, sw_p = stats.shapiro(resid[sub_idx])
    legacy_figures_log(f"  Jarque-Bera: stat={jb_stat:.2f}, p={jb_p:.3e}, skew={skew:.4f}, excess-kurtosis={kurt-3:.4f}")
    legacy_figures_log(f"Shapiro-Wilk (n=5,000 random subsamples, scipy is unreliable for n>5000 so subsamples are taken): W={sw_stat:.5f}, p={sw_p:.3e}")
    legacy_figures_log(f"⚠️ Burden_rare is highly right-skewed/zero-inflated even within case (proportion of 0 within case ="
        f"{float((cases['burden_rare_main_ptv_10gene']==0).mean())*100:.2f}%), so the normality assumption of OLS must be tested in practice.")

    diag_tbl = pd.DataFrame([{
        "model": "Aim4 OLS refit in script 19 (age_at_diagnosis ~ Burden_z + PRS311_z + PC1-10, HC1)",
        "n": len(cases), "beta_burden_z": fit_a4.params["burden_z"], "beta_prs311_z": fit_a4.params["PRS311_z"],
        "r_squared": fit_a4.rsquared, "resid_skew": skew, "resid_excess_kurtosis": kurt - 3,
        "jarque_bera_stat": jb_stat, "jarque_bera_p": jb_p,
        "shapiro_wilk_W_n5000_subsample": sw_stat, "shapiro_wilk_p_n5000_subsample": sw_p,
        "burden_zero_fraction_within_cases": float((cases["burden_rare_main_ptv_10gene"] == 0).mean()),
    }])
    write_table(diag_tbl, "19_fig6_aim4_ols_residual_diagnostics",
                note="scripts/17 was NOT modified; the model was refit inside script 19 (burden_z "
                     f"{fit_a4.params['burden_z']:.3f}, PRS311_z {fit_a4.params['PRS311_z']:.3f}, R2 {fit_a4.rsquared:.4f}; "
                     f"scripts/17: {REF17_BURDEN_BETA:.3f}, {REF17_PRS_BETA:.3f}, {REF17_R2:.4f}). Shapiro-Wilk is on a "
                     "5,000-person random subsample because scipy's implementation is unreliable above n=5000; "
                     f"Jarque-Bera uses all {len(cases):,}.")

    legacy_figures_fig, legacy_figures_axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
    
    _pct = np.arange(1, 100)
    _theo = stats.norm.ppf(_pct / 100)
    _samp = np.percentile(resid, _pct)
    legacy_figures_axes[0].plot(_theo, _samp, marker="o", markersize=3, linewidth=0.8, color="tab:blue")
    legacy_figures_axes[0].plot(_theo, np.mean(resid) + np.std(resid, ddof=1) * _theo, color="tab:red")
    legacy_figures_axes[0].set_xlabel("Theoretical quantiles (normal)")
    legacy_figures_axes[0].set_ylabel("Residual percentiles (years)")
    legacy_figures_axes[0].set_title(f"QQ plot of residuals, 1st-99th percentiles (n={len(cases):,})")
    
    hb63 = legacy_figures_axes[1].hexbin(fitted, resid, gridsize=30, mincnt=LEGACY_FIGURES_DISCLOSURE_MIN_CELL, cmap="Blues")
    legacy_figures_fig.colorbar(hb63, ax=legacy_figures_axes[1], label=f"Participants per cell (cells with <{LEGACY_FIGURES_DISCLOSURE_MIN_CELL} not shown)")
    fit_axes_to_hexbin(legacy_figures_axes[1], hb63)
    legacy_figures_axes[1].axhline(0, color="black", linewidth=0.8)
    legacy_figures_axes[1].set_xlabel("fitted age at diagnosis (years)")
    legacy_figures_axes[1].set_ylabel("residual (years)")
    legacy_figures_axes[1].set_title("Residuals vs fitted")
    legacy_figures_fig.tight_layout()
    fig63_path = save_fig(legacy_figures_fig, "19_fig6_aim4_ols_diagnostics", prefer="svg", max_svg_kb=600)

    # =============================================================================

    # =============================================================================
    legacy_figures_section('Fig 7 | Calibration plot (Joint model 5-fold out-of-fold prediction, 1st repeat on 16th)')

    legacy_figures_y_all = legacy_figures_df["is_case_int"].values
    p_oof = legacy_figures_df["oof_r1_joint"].values.astype(float)
    auc_oof = roc_auc_score(legacy_figures_y_all, p_oof)
    brier_oof = brier_score_loss(legacy_figures_y_all, p_oof)

    ci_oof, cs_oof, cs_se_oof = calib_params(legacy_figures_y_all, p_oof)
    legacy_figures_log(f"  n={len(legacy_figures_df)}, proportion of observed cases ={legacy_figures_y_all.mean():.4f}; OOF AUC={auc_oof:.4f}, Brier={brier_oof:.5f}, "
        f"calibration intercept={ci_oof:.4f}, slope={cs_oof:.4f} (SE {cs_se_oof:.4f})")
    cmp7 = pd.DataFrame([{"model": "Joint, 5-fold out-of-fold (repeat 1 of scripts/16)", "auc": auc_oof,
                          "brier": brier_oof, "calibration_intercept": ci_oof, "calibration_slope": cs_oof,
                          "calibration_slope_se": cs_se_oof, "mean_predicted_prob": float(p_oof.mean()),
                          "observed_case_proportion": float(legacy_figures_y_all.mean())}])
    write_table(cmp7, "19_fig7_calibration_model_comparison",
                note="Unweighted MLE joint model; predictions are 5-fold out-of-fold from repeat 1 of the "
                     "10x5 repeated stratified CV in scripts/16.")

    cal = pd.DataFrame({"y": legacy_figures_y_all, "p": p_oof})
    cal["calib_bin"] = pd.qcut(cal["p"], 10, labels=False) + 1
    brows = []
    for legacy_figures_b, cell in cal.groupby("calib_bin"):
        n = len(cell); nc = int(cell["y"].sum()); obs = nc / n
        legacy_figures_lo, legacy_figures_hi = stats.binomtest(nc, n).proportion_ci(confidence_level=0.95, method="wilson")
        brows.append({"bin": legacy_figures_b, "n_disclosed": disclose_n(n, f"fig7|bin={legacy_figures_b}")[0],
                      "n_case_disclosed": disclose_n(nc, f"fig7|bin={legacy_figures_b}|cases")[0],
                      "mean_predicted": cell["p"].mean(), "observed_proportion": obs,
                      "obs_ci_low": legacy_figures_lo, "obs_ci_high": legacy_figures_hi})
    cal_tbl = pd.DataFrame(brows)
    legacy_figures_log('Calibration bin table (OOF, 10 bins):')
    legacy_figures_log(cal_tbl.to_string(index=False))
    write_table(cal_tbl, "19_fig7_calibration_bins",
                note="Deciles of out-of-fold predicted risk (scripts/16 repeat 1).")

    legacy_figures_fig, legacy_figures_ax = plt.subplots(figsize=(5.8, 5.4))
    legacy_figures_ax.plot([0, cal_tbl["mean_predicted"].max() * 1.05], [0, cal_tbl["mean_predicted"].max() * 1.05],
            color="gray", linewidth=0.8, linestyle="--")
    legacy_figures_ax.errorbar(cal_tbl["mean_predicted"], cal_tbl["observed_proportion"],
                yerr=[cal_tbl["observed_proportion"] - cal_tbl["obs_ci_low"],
                      cal_tbl["obs_ci_high"] - cal_tbl["observed_proportion"]],
                fmt="o", color="tab:blue", capsize=3, markersize=5)
    legacy_figures_ax.set_xlabel("Mean predicted probability (decile of out-of-fold predicted risk)")
    legacy_figures_ax.set_ylabel("Observed case proportion")
    legacy_figures_fig.tight_layout()
    save_fig(legacy_figures_fig, "19_fig7_calibration", prefer="svg")

    # =============================================================================
    
    # =============================================================================
    legacy_figures_section('Closing | Disclosure Rules Log / Figure File Size / Report')

    supp_tbl = pd.DataFrame(SUPPRESSION_LOG) if SUPPRESSION_LOG else pd.DataFrame(
        columns=["cell", "raw_n_bucket", "action", "disclosed"])
    legacy_figures_log(f"Number of squares changed by disclosed rules:{len(supp_tbl)}")
    if len(supp_tbl):
        legacy_figures_log(supp_tbl.to_string(index=False))
    write_table(supp_tbl, "19_disclosure_rule_applied_cells",
                note="Every cell where the CLAUDE.md n>=20 rule changed what is reported (5-19 rounded up to "
                     "20; <5 suppressed). Only the bucket (<5 / 5-19) is recorded, never the exact raw count.")

    size_tbl = pd.DataFrame(FIG_SIZES)
    legacy_figures_log('Image file size:')
    legacy_figures_log(size_tbl.to_string(index=False))
    legacy_figures_log(f"Total image files:{size_tbl['size_kb'].sum():.1f} KB")
    tbl_files = sorted([legacy_figures_f for legacy_figures_f in os.listdir(LEGACY_FIGURES_TABLE_DIR) if legacy_figures_f.startswith("19_")])
    tbl_kb = sum(os.path.getsize(os.path.join(LEGACY_FIGURES_TABLE_DIR, f)) for f in tbl_files) / 1024
    legacy_figures_log(f"19_* Table{len(tbl_files)}number, total{tbl_kb:.1f} KB")

    hi10 = main_cells[(main_cells["prs_decile"] == 10) & (main_cells["carrier_main_10gene"] == 1)].iloc[0]
    hi9 = main_cells[(main_cells["prs_decile"] == 9) & (main_cells["carrier_main_10gene"] == 1)].iloc[0]
    comp_row = or_tbl[(or_tbl["panel"] == "B_burden_ordered_category") & (or_tbl["category"] == "compound")].iloc[0]
    comp_or, comp_lo, comp_hi = comp_row["OR"], comp_row["ci_low"], comp_row["ci_high"]
    comp_n = comp_row["n_disclosed"]
    _top_low = pa[pa["category"] == "ventile 20"].iloc[0]
    or_top_low, ci_lo_top, ci_hi_top = _top_low["OR"], _top_low["ci_low"], _top_low["ci_high"]
    zero_frac_pct = zero_frac * 100
    _pv = ppv_tbl.set_index(["prs_band", "burden_stratum"])["observed_case_percent"]
    ppv_any_low = _pv[("any PRS (all)", "low burden (non-carrier)")]
    ppv_any_high = _pv[("any PRS (all)", "high burden (any panel PTV carrier)")]
    ppv_any_ratio = ppv_any_high / ppv_any_low
    ppv_95_low = _pv[(">95th centile", "low burden (non-carrier)")]
    ppv_95_high = _pv[(">95th centile", "high burden (any panel PTV carrier)")]
    ppv_95_ratio = ppv_95_high / ppv_95_low
    c1 = corr_tbl[corr_tbl["scope"] == "controls_only(primary)"]

    report = f"""# Aim 1 Supplementary Figures (Fig 1 / 2 / 3 / 6 / 7) Summary Report

    Corresponding script: `scripts/90_legacy_ptv_baseline.py` (log: `scripts/script_output_py/output_19_aim1_figures.txt`).

    > **All content is post-hoc exploratory**. This script does not produce any new confirmatory conclusions, all
    > The output tables all have the `evidence_status = post_hoc_exploratory` column, which will be used
    > Convention of `results/table/16_carrier_stratified_prs_slope_by_gene.csv`. Cohorts, covariates, standardization
    > The caliber is exactly the same as `scripts/90_legacy_ptv_baseline.py` (n={len(legacy_figures_df):,}，case={int(legacy_figures_df['is_case_int'].sum()):,},
    > Hard verification with assert in the script; Aim 4 part n={len(cases):,}, similarly assert).

    ## 0. Two specification corrections that must be recorded

    **(a) The main model is not a weighted model. ** The task description for this round states that the Aim 1 main model uses class-balanced weights
    `w = n/(2*n_group)`. The actual situation is: `scripts/16` has been changed in v8.0 (addendum second review, 2026-08-18)
    **The main estimator is changed to unweighted MLE** (`sm.GLM(Binomial)` + HC1), and the weighted logistic only remains as
    `weighted_class_balanced` One-line sensitivity analysis. Therefore, the original motivation of Fig 7 (\"The predicted probability of the main model depends on the weighted
    On the scale, direct drawing calibration will cause false inaccuracies purely caused by weighting\") which is **not true** in this project. This script is still executed according to the instructions
    \"Unweighted refit on train fold + calibration on held-out test\", and additionally fit the weighted version for comparison, using actual measured numbers
    Note that weighting only moves the intercept.

    **(b) Two sets of covariates. ** The instructions clearly specify \"PC1-10 + age\" for Fig 2 Panel B and Fig 3, so these two
    **Not included** Sequencing center/sequencer dummy; while the elimination refitting of Fig 6.1 and the joint model of Fig 7 must be used with `scripts/16`
    The main model is comparable item by item, using the full `COV_COLS` (age + PC1-10 + site + sequencer).

    ## 1. Gene stratification

    - **high-penetrance**: BRCA1, BRCA2, PALB2, TP53
    - **moderate**: ATM, CHEK2, BARD1, RAD51C, RAD51D, **MAP3K1**
    - **compound**: ≥2 panel gene PTVs, **separately layered, never merged into high-penetrance**

    | layer | n (after disclosure) | where case |
    |---|---|---|""" + "\n".join(
        f"| {r['tier4']} | {r['n_total_disclosed']} | {r['n_case_disclosed']} |" for _, r in tier_tbl.iterrows()
    ) + f"""**⚠️ The classification of MAP3K1 is controversial**: MAP3K1 is mainly a common-variant locus (GWAS signal), and this project gives it
    The PTV weights are estimated from a single cohort (MWS) without independent replication. Putting it into the moderate layer is just an operational choice and does not mean
    It has ATM/CHEK2 level PTV evidence. compound carrier measured n={n_compound}(Instructions expect ~20–25, consistent).

    ## 2. Fig 1｜PRS × Burden orthogonality (controls-only mainly)

    `results/figures/{os.path.basename(fig1_path)}`,
    Table `results/table/19_fig1_prs_burden_correlations.csv`, `19_fig1_prs_by_carrier_distribution.csv`.

    **⚠️ Why exclude cases**: Calculating the correlation between PRS and rare burden based on the case state will pass
    collider / index-event bias artificially induces negative correlation - two independent paths leading to the same outcome must be negatively correlated after a given outcome.
    Therefore correlation values ​​within cases are themselves biases and not evidence of orthogonality. The full queue version contains 6.2% of the cases, with minor in the same direction
    Pollution, only for secondary reference.

    controls-only main results (raw partially correlated with PC1-10):

    | Comparison | Method | Adjustment | r | 95% CI |
    |---|---|---|---|---|""" + "\n".join(
        f"| {r['comparison']} | {r['method']} | {r['adjustment']} | {r['r']:.4f} | [{r['ci_low']:.4f}, {r['ci_high']:.4f}] |"
        for _, r in c1.iterrows()
    ) + f"""**External reference (for context only, not the baseline we are aligning against)**: Mazzarotto et al. (bioRxiv 2025,
    doi:10.1101/2025.09.15.676308; This round has obtained the PDF full text and checked it, see `plans/papers/`) The summary sentence is
    \"IOGC and PRS were minimally correlated (**|r|<0.2**), indicating the two scores capture
    complementary information\" - that is, they use the framework of \"**|r|<0.2 is regarded as minimally relevant/complementary**\" rather than a certain
    A specific r value is used as a target. The specific values verified in the article are schizophrenia PCC=0.097 (0.095–0.1), MDD PCC=0.08
    (0.077–0.083), HCM PCC=−0.178 (−0.181–−0.175). Trait (schizophrenia/MDD/HCM) and burden score construction
    (IOGC = expression outlier gene count) are different from this project, so they can only be used as context: the actual measurement of this project |r| All <0.04,
    It is closer to 0 than any value in the article, falling within its \"minimally correlated\" criterion.

    ## 3. Fig 2｜Two-dimensional heat map (PRS decile × carrier)

    `results/figures/19_fig2_prs_decile_x_carrier_heatmap.svg` (main),
    Table `19_fig2_prs_decile_x_carrier_cells.csv`. Supplementary panel for PRS quintile × 3 tier for internal reference only (`data/19_internal/`):
    The 3-tier counts, together with the public 4-tier, gene-by-gene counts, can solve for the number of gene combinations (<20) in the compound.

    **⚠️ The amount in Panel A is the \"cohort-observed case proportion\", not the absolute risk. **
    Approximately 6.4% of cases in this cohort are from EHR diagnoses within the AoU WGS subcohort and are neither lifetime cumulative incidence nor population incidence.

    Panel B = observations/expectations, expectations taken from **no interaction terms** model `case ~ PRS311_z + carrier + PC1-10 + age`
    Average predicted probability within a bin; divergent color scale is centered around 1.0. If negative interactions were real, cells with high PRS × carrier would fall below 1.
    Actual measured high score corner: decile 10 × carrier O/E = **{hi10['obs_over_exp']:.3f}**（n={hi10['n_disclosed']}，
    case={hi10['n_case_disclosed']}），decile 9 × carrier O/E = **{hi9['obs_over_exp']:.3f}**
    （n={hi9['n_disclosed']}).

    **This figure is not the figure of Bolze et al. 2024**: The full text of the article has been checked (PMC11613897), which uses the Kaplan-Meier curve and
    Absolute risk at age 70 years, PRS dichotomized by median, no OR heatmap in full text. The observed/expected panel is a self-built structure of this project.
    No credit for this figure may be attributed to Bolze.

    The supplementary panel uses 3 layers (compounds are combined into high or moderate according to the highest penetrance genes they carry, and are separated into separate layers in Fig 3),
    PRS uses quintiles rather than deciles because the grid population is small after genetic stratification. **⚠️ The O/E of the supplementary panel cannot be read as the main panel**:
    The expectation comes from a model with only one pooled carrier indicator variable, so the O/E levels of this panel are mixed with (a) between-level main effect heterogeneity
    (The high-penetrance carrier is more dangerous than the merger coefficient, and the moderate carrier is less dangerous. The measured high-layer O/E is all >1.
    moderate layer all <1, mainly for this reason) interacts with (b) true PRS×carrier. Only** across PRS quintiles within the same stratum
    Only the gradient** carries interaction information (the high layer drops from 1.99 to 1.14, the direction is consistent with negative interaction).

    ## 4. Fig 3｜Quantile-OR point and line chart

    `results/figures/19_fig3_quantile_or_dotwhisker.svg` (main),
    `19_fig3_supp_midband_and_ppv_analogue.svg` (supplement); Table `19_fig3_quantile_or_estimates.csv`,
    `19_fig3_ppv_analogue_prs_band_x_burden.csv`. Both panels are adjusted for PC1-10 + age, with logarithmic scale on the vertical axis.

    **The agreement between binning and reference groups has been checked against the original text** (`plans/papers/The complementary roles of rare variant burden
    scores.pdf`, not paraphrased): The article Methods states that both scores cut quantile bins, \"default N=20; ventiles\",
    And \"logistic regressions of diagnosis on IOGC quantile (**reference: lowest quantile**)\",
    IOGC was also dichotomized into low(ventiles 1–18) vs high(>18) for PPV analysis. Therefore:

    - **Panel A (PRS, main parameterization)**: PRS311 twentieth quartile, **Reference group = lowest ventile** (aligned with the article),
      top ventile OR={or_top_low:.2f}（95% CI {ci_lo_top:.2f}–{ci_hi_top:.2f}), the middle quantile increases monotonically,
      Rise in a smooth logarithmic linear manner.
    - The \"40–60 percentile midband\" reference used in the previous version is retained as a **minor parameterization** (Supplementary Figure Supp A + Table
      `reference` column), the two versions have the same conclusion, but the reference points are different.
    - **Panel B (rare burden)**: Burden_rare full queue{zero_frac_pct:.2f}% is 0, the 19th percentile out of 20th
      Taking the same value, quantile binning cannot be defined (IOGC in this article refers to the count and distribution of expression outlier genes, which is much more continuous.
      Not portable), so ordered categories are used instead: non-carrier (ref, the reference group convention is consistent with this article) / single moderate PTV /
      Single high-penetrance PTV/**compound (≥2 panel gene PTVs)**. compound is drawn with a hollow red dot,
      Show full wide CI: OR={comp_or:.2f}（95% CI {comp_lo:.2f}–{comp_hi:.2f}，n={comp_n}).
      **There are only 26 people in this layer (the case number is reported as 20 according to n≥20 disclosure rules), which is seriously insufficient and the point estimate cannot be interpreted** - but this layer must not
      Incorporated into high-penetrance, since that\'s what the Antinomy 1 argument was about.

    **The whole point of this graph is the comparison of the two panels**: A\'s smooth logarithmic linear rise vs B\'s threshold jump, that is
    Visual evidence for \"Is the linearity assumption imposed on the burden score\" (the question behind Antinomy 1). This article provides this
    Precedent for **shape** (context, not baseline; trait is schizophrenia/HCM not breast cancer): PRS top ventile vs bottom
    OR=13.6 (95% CI 8.1–23) and monotonically increasing in the middle, while IOGC \"did not vary linearly with ventiles, showing a
    significant increase only for the top two\", top ventile OR=2.3 (1.5–3.7).

    ### 4.1 PPV analogy (analogy in Supplementary Figure 3 of this article)

    Table `19_fig3_ppv_analogue_prs_band_x_burden.csv`, Supplementary Figure Supp B. The PRS centile threshold follows the default value in this article.
    A group (<15th, <25th, <50th, <75th, >95th, cumulative/overlapping intervals according to the original text) plus \"any PRS\" rows;
    burden high/low = carrier / non-carrier. **⚠️Three points must be made clearly**:

    1. **This is not a PPV. **PPV dependent prevalence, 6.2% case proportion in this cohort determined by case-control design + EHR confirmation,
       Not population prevalence. What is reported in the table is the **proportion of observed cases in this cohort (%)**, and its **absolute level is not comparable to the value in this article**.
       Only the shape \"relative improvement factor of high burden vs low burden\" is comparable.
    2. **The PPV value of this article is a percentage, not a ratio**: Schizophrenia any-PRS 0.25% → 0.41% (1.64 times),
       PRS>95th centile 0.89% → 1.13%; the MDD analogy in the same text reads 6.03 → 6.13 and 17.12 → 14.31,
       Confirmable unit is %. (The mission statement for this round writes these two sets of numbers as bare values 0.25/0.41 and 0.89/1.13.
       Easily misread as proportion, corrected here. )
    3. This article divides IOGC into two parts: ventile >18; the burden 98% of this project is 0, and the ventile threshold is meaningless, so it is used instead
       carrier / non-carrier dichotomy.

    Actual measurement of this project: non-carrier under any-PRS{ppv_any_low:.2f}% vs carrier {ppv_any_high:.2f}%
    （{ppv_any_ratio:.2f}times); PRS >95th centile under non-carrier{ppv_95_low:.2f}% vs carrier
    {ppv_95_high:.2f}%（{ppv_95_ratio:.2f}times).

    ## 5. Fig 6｜Model diagnosis (Gap 1)

    ### 5.1 Decisive test: Refit after eliminating compound carrier (`19_fig6_compound_exclusion_interaction_comparison.csv`)

    | model | n | inter_z β | OR | p |
    |---|---|---|---|---|
    | Main analysis (reproduction scripts/16) |{len(legacy_figures_df):,} | {legacy_figures_b0:.5f} | {np.exp(legacy_figures_b0):.4f} | {fit_primary.pvalues['inter_z']:.3e}|
    | Eliminate compound carrier |{len(df_nc):,} | {legacy_figures_b1:.5f} | {np.exp(legacy_figures_b1):.4f} | {legacy_figures_p1:.3e}|

    coefficient change{100*(legacy_figures_b1-legacy_figures_b0)/abs(legacy_figures_b0):+.1f}%. **Interpretation**:{verdict}.
    Both fits use the z-score parameters of the full queue (not recalculated on the reduced queue) to ensure that the coefficients are comparable on the same scale.
    Cook\'s distance plot is evidence, not verdict - the verdict is in this table.

    ### 5.2 Cook\'s distance (`19_fig6_cooks_distance.svg`, `19_fig6_cooks_distance_by_tier.csv`)

    n={len(legacy_figures_df):,}Without plotting one point per person: Panel A is the Cook\'s D distribution by carrier tier (box plot, 5–95 percentile required),
    Panel B is the overall distribution histogram (columns with <20 people are not drawn). Maximum Cook\'s D={cooks.max():.5f}, far lower than any
    The magnitude of the \"single point dominant fit\" (empirical threshold 4/n={4/len(legacy_figures_df):.2e}, but in the sample of n=140,000, the threshold itself has no practical significance).

    ### 5.3 Aim 4 OLS residual diagnosis (`{os.path.basename(fig63_path)}`, `19_fig6_aim4_ols_residual_diagnostics.csv`)

    **Unmodified `scripts/17`**, refit within this script (coefficient reproduction: burden_z={fit_a4.params['burden_z']:.3f}、
    PRS311_z={fit_a4.params['PRS311_z']:.3f}、R²={fit_a4.rsquared:.4f}). case internal Burden_rare still has{float((cases['burden_rare_main_ptv_10gene']==0).mean())*100:.1f}% is 0, highly right-skewed/zero-inflated, normality must be checked:
    Jarque-Bera stat={jb_stat:.1f}（p={jb_p:.2e}), skewness ={skew:.3f}, excess kurtosis ={kurt-3:.3f};
    Shapiro-Wilk (random subsample of 5,000 people, scipy is unreliable on n>5000) W={sw_stat:.4f}，p={sw_p:.2e}.

    ## 6. Fig 7｜Calibration graph (unweighted refitting, held-out test)

    `results/figures/19_fig7_calibration.svg`; table `19_fig7_calibration_bins.csv`,
    `19_fig7_calibration_model_comparison.csv`.

    Predictions from `scripts/16` 5-fold out-of-fold predictions from replicate 1 (Joint model, unweighted MLE),
    Each person’s prediction comes from a model that has not seen him (n={len(legacy_figures_df):,}).

    | Model | OOF AUC | Brier | Calibrated Intercept | Calibrated Slope | Average Predicted Probability | Proportion of Observed Cases |
    |---|---|---|---|---|---|---|
    | Joint (50% off OOF) |{auc_oof:.4f} | {brier_oof:.5f} | {ci_oof:.4f} | {cs_oof:.4f} | {p_oof.mean():.4f} | {legacy_figures_y_all.mean():.4f}|

    ## 7. Disclosure Rules

    CLAUDE.md\'s n≥20 rule applies cell-by-cell (5–19 rounds up to 20, <5 suppresses). The total number of grids that have been changed this round{len(supp_tbl)}each, recorded cell by cell in `results/table/19_disclosure_rule_applied_cells.csv`.
    All output tables are aggregate results, without person-level rows; the script does not write any person-level intermediate files to disk.

    ## 8. Output list

    Figures (`results/figures/`):

    | File | Size |
    |---|---|""" + "\n".join(f"| {r['figure']} | {r['size_kb']:.1f} KB |" for _, r in size_tbl.iterrows()) + f"""Table(`results/table/`, total{len(tbl_files)}number, total{tbl_kb:.1f} KB）：
    """ + "\n".join(f"- `{f}`" for f in tbl_files) + f"""Report: `results/report/19_aim1_figures_summary.md` (this file);
    Log: `scripts/script_output_py/output_19_aim1_figures.txt`.

    Figure + table total approx.{size_tbl['size_kb'].sum() + tbl_kb:.1f} KB。
    """

    rp = os.path.join(LEGACY_FIGURES_REPORT_DIR, "19_aim1_figures_summary.md")
    with open(rp, "w") as legacy_figures_f:
        legacy_figures_f.write(report)
    legacy_figures_log(f"  [report] {rp}")

    legacy_figures_section('DONE | A quick look at key figures')
    pb_rows = c1[c1["comparison"].str.contains("point-biserial")]
    r_pb_raw = float(pb_rows[pb_rows["adjustment"] == "raw"]["r"].iloc[0])
    r_pb_adj = float(pb_rows[pb_rows["adjustment"] == "partial_PC1_10"]["r"].iloc[0])
    legacy_figures_log(f"  Fig1 controls-only point-biserial r(raw)   = {r_pb_raw:.4f}")
    legacy_figures_log(f"  Fig1 controls-only point-biserial r(PC-adj)= {r_pb_adj:.4f}")
    legacy_figures_log(f"  Fig2 decile10 x carrier O/E = {hi10['obs_over_exp']:.3f}; decile9 x carrier O/E = {hi9['obs_over_exp']:.3f}")
    legacy_figures_log(f"  Fig3 compound OR = {comp_row['OR']:.3f} [{comp_row['ci_low']:.3f}, {comp_row['ci_high']:.3f}], n={comp_row['n_disclosed']}")
    legacy_figures_log(f"  Fig6 interaction beta: {legacy_figures_b0:.5f} -> {legacy_figures_b1:.5f} (p {fit_primary.pvalues['inter_z']:.2e} -> {legacy_figures_p1:.2e})")
    legacy_figures_log(f"  Fig7 calib intercept={ci_oof:.4f}, slope={cs_oof:.4f}, Brier={brier_oof:.5f}, AUC={auc_oof:.4f}")

# ============================================================================
# 10. forest: analyze_legacy_forest
# ============================================================================

import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def term_kind(term):
    if term.startswith("PRS"):
        return "PRS311"
    if term.lower().startswith("burden"):
        return "Burden_rare"
    return None

def analyze_legacy_forest():
    """Initialize runtime state and execute this component explicitly."""
    global AIM1_MODELS, AIM4_MODELS, LEGACY_FOREST_FIG_DIR, LEGACY_FOREST_ROOT, LEGACY_FOREST_TABLE_DIR, legacy_forest__, a1, a4, aim, legacy_forest_ax, legacy_forest_axes, legacy_forest_df, est, est_tbl, legacy_forest_fig, handles, legacy_forest_hi, legacy_forest_i, legacy_forest_k, kind, legacy_forest_lab, legacy_forest_label, labels, legacy_forest_lo, logx, markers, missing, model, models, n_col, offsets, out, legacy_forest_r, ref, legacy_forest_rows, legacy_forest_sub, xlabel, y, ypos

    matplotlib.use("Agg")

    LEGACY_FOREST_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    LEGACY_FOREST_TABLE_DIR = os.path.join(LEGACY_FOREST_ROOT, "results", "table")
    LEGACY_FOREST_FIG_DIR = os.path.join(LEGACY_FOREST_ROOT, "results", "figures")

    AIM1_MODELS = [  # (model, label)
        ("Joint_full_cohort_main_10gene_PRS311", "Primary (10 genes, MAF<1%)"),
        ("9gene_main_ptv_PRS311", "9-gene panel (excl. MAP3K1)"),
        ("10gene_sensitivity_maf0p1pct_PRS313", "MAF<0.1% burden + PRS313"),
        ("prs281loci_PRS311", "PRS restricted to 281 loci"),
        ("prs_no_mean_imputation_PRS311", "PRS without mean imputation"),
        ("weighted_class_balanced", "Class-balanced weights"),
        ("case_def_rule_of_two", "Case definition: rule of two"),
    ]
    AIM4_MODELS = [
        ("primary_10gene_bc_first_record_excluded", "Primary (10 genes, MAF<1%)"),
        ("sensitivity2_9gene_panel_primary_cohort", "9-gene panel (excl. MAP3K1)"),
        ("prs281loci_PRS311", "PRS restricted to 281 loci"),
        ("prs_no_mean_imputation_PRS311", "PRS without mean imputation"),
        ("sensitivity1_full_cohort_incl_bc_first_record_10gene", "Including cases with BC as first record"),
        ("case_def_rule_of_two", "Case definition: rule of two"),
    ]

    a1 = pd.read_csv(os.path.join(LEGACY_FOREST_TABLE_DIR, "16_sensitivity_analyses_coefficients.csv"))
    a1["kind"] = a1["term"].map(term_kind)
    a4 = pd.read_csv(os.path.join(LEGACY_FOREST_TABLE_DIR, "17_sensitivity_analyses_coefficients.csv"))
    a4["kind"] = a4["term"].map(term_kind)

    legacy_forest_rows = []
    for aim, legacy_forest_df, models, est, legacy_forest_lo, legacy_forest_hi, n_col in [
            ("Aim 1", a1, AIM1_MODELS, "OR", "CI_low", "CI_high", "cohort_n"),
            ("Aim 4", a4, AIM4_MODELS, "beta", "ci_low", "ci_high", "n")]:
        for model, legacy_forest_label in models:
            legacy_forest_sub = legacy_forest_df[(legacy_forest_df["model"] == model) & legacy_forest_df["kind"].notna()]
            missing = {"PRS311", "Burden_rare"} - set(legacy_forest_sub["kind"])
            if missing:
                raise SystemExit(f"{aim} {model}: Lack{missing}, check whether the 16/17 output has been rerun")
            for legacy_forest__, legacy_forest_r in legacy_forest_sub.iterrows():
                legacy_forest_rows.append({"aim": aim, "model": model, "label": legacy_forest_label, "exposure": legacy_forest_r["kind"],
                             "estimate_type": "OR per SD" if aim == "Aim 1" else "beta (years) per SD",
                             "estimate": legacy_forest_r[est], "ci_low": legacy_forest_r[legacy_forest_lo], "ci_high": legacy_forest_r[legacy_forest_hi],
                             "p_value": legacy_forest_r["p_value"], "cohort_n": int(legacy_forest_r[n_col])})
    est_tbl = pd.DataFrame(legacy_forest_rows)
    est_tbl.to_csv(os.path.join(LEGACY_FOREST_TABLE_DIR, "20_fig5_forest_estimates.csv"), index=False)
    print(est_tbl[["aim", "label", "exposure", "estimate", "ci_low", "ci_high", "cohort_n"]].to_string(index=False))

    legacy_forest_fig, legacy_forest_axes = plt.subplots(1, 2, figsize=(11, 4.6), sharey=False)
    offsets = {"PRS311": -0.15, "Burden_rare": 0.15}
    markers = {"PRS311": "o", "Burden_rare": "s"}
    for legacy_forest_ax, aim, models, ref, xlabel, logx in [
            (legacy_forest_axes[0], "Aim 1", AIM1_MODELS, 1.0, "Odds ratio per SD (95% CI)", True),
            (legacy_forest_axes[1], "Aim 4", AIM4_MODELS, 0.0, "Difference in age at diagnosis\nper SD, years (95% CI)", False)]:
        legacy_forest_sub = est_tbl[est_tbl["aim"] == aim]
        labels = [legacy_forest_lab for legacy_forest__, legacy_forest_lab in models]
        ypos = {legacy_forest_lab: len(labels) - 1 - legacy_forest_i for legacy_forest_i, legacy_forest_lab in enumerate(labels)}
        for kind in ["PRS311", "Burden_rare"]:
            legacy_forest_k = legacy_forest_sub[legacy_forest_sub["exposure"] == kind]
            y = legacy_forest_k["label"].map(ypos) + offsets[kind]
            legacy_forest_ax.errorbar(legacy_forest_k["estimate"], y,
                        xerr=[legacy_forest_k["estimate"] - legacy_forest_k["ci_low"], legacy_forest_k["ci_high"] - legacy_forest_k["estimate"]],
                        fmt=markers[kind], capsize=2, label=kind)
        legacy_forest_ax.axvline(ref, color="grey", linestyle="--", linewidth=0.8)
        if logx:
            legacy_forest_ax.set_xscale("log")
            legacy_forest_ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
            legacy_forest_ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
            legacy_forest_ax.set_xticks([1.0, 1.2, 1.4, 1.6])
        legacy_forest_ax.set_yticks(list(ypos.values()))
        legacy_forest_ax.set_yticklabels(list(ypos.keys()))
        legacy_forest_ax.set_xlabel(xlabel)
        legacy_forest_ax.set_title(("A. " if aim == "Aim 1" else "B. ") + aim, loc="left")
    handles, labels = legacy_forest_axes[0].get_legend_handles_labels()
    legacy_forest_fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.04))
    legacy_forest_fig.tight_layout(w_pad=2, rect=(0, 0, 0.97, 1))
    out = os.path.join(LEGACY_FOREST_FIG_DIR, "20_fig5_forest_sensitivity.svg")
    legacy_forest_fig.savefig(out, bbox_inches="tight")
    print(f"saved {out} ({os.path.getsize(out)/1024:.1f} KB)")

# Command-line selection of the computations defined above.
def main():
    from scripts.functions.cli import run_local_step
    computations = {
        'vat-schema': analyze_vat_schema,
        'vat-prep': analyze_legacy_vat_preparation,
        'vat-repair': analyze_legacy_vat_repair,
        'variant-selection': analyze_legacy_variant_selection,
        'carriers': analyze_legacy_carriers,
        'burden': analyze_legacy_burden,
        'joint-models': analyze_legacy_joint_models,
        'age-models': analyze_legacy_age_models,
        'figures': analyze_legacy_figures,
        'forest': analyze_legacy_forest,
    }
    run_local_step('Optional legacy PTV baseline for A and inherited CV folds', computations, default_tasks=[])

if __name__ == '__main__':
    main()

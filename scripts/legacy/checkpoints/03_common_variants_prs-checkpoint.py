"""
- Match PRS variants and orient alleles.
- Calculate weighted PRS.
- Audit variant matches and scoring.
"""
from __future__ import annotations
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
# ============================================================================
# 1. prs: calculate_common_variant_prs
# ============================================================================

import concurrent.futures
import os
import subprocess
import sys
import time
import numpy as np
import pandas as pd
import yaml

def run(cmd, **kw):
    print(f"$ {cmd}")
    res = subprocess.run(cmd, shell=True, capture_output=True, text=True, **kw)
    if res.stdout:
        print(res.stdout)
    if res.stderr:
        print(res.stderr, file=sys.stderr)
    if res.returncode != 0:
        print(f"!!! exit code {res.returncode}", file=sys.stderr)
    return res

def section(title):
    print("\n" + "=" * 90)
    print(title)
    print("=" * 90)

def download_scoring_files():
    section('Step 1 | Download PGS000004 harmonized scoring file (PGS Catalog FTP, no manual liftover)')
    for url, fname in [
        (HM38_URL, f"{PGS_ID}_hmPOS_GRCh38.txt.gz"),
        (HM37_URL, f"{PGS_ID}_hmPOS_GRCh37.txt.gz"),
    ]:
        dest = f"{PRS_RAW_DIR}/{fname}"
        if os.path.exists(dest) and os.path.getsize(dest) > 0:
            print(f"already present: {dest} ({os.path.getsize(dest)} bytes)")
            continue
        run(f"wget -q '{url}' -O {dest}")
        print(f"downloaded {dest}: {os.path.getsize(dest)} bytes")

def load_cohort():
    section('Step 2 | Load the genetics EUR queue (BigQuery eur_female_breast_cancer_cohort, output 01)')
    client = bigquery.Client()
    q = f"""
        SELECT person_id, is_case
        FROM `{GOOGLE_PROJECT}.my_analysis_results.eur_female_breast_cancer_cohort`
    """
    df = client.query(q).to_dataframe()
    df["person_id"] = df["person_id"].astype("int64")
    df = df.sort_values("person_id").reset_index(drop=True)
    print(f"cohort n={len(df)}  case={int((df.is_case==1).sum())}  control={int((df.is_case==0).sum())}")
    assert len(df) == 170890, f"The number of people did not match the expected 170,890:{len(df)}"
    return df

def build_target_table():
    section('Step 3 | Construct PRS313 target site table')
    pgs = pd.read_csv(f"{PRS_RAW_DIR}/{PGS_ID}_hmPOS_GRCh38.txt.gz", sep="\t", comment="#")
    assert len(pgs) == 313, f"Expected 313 sites, actually read{len(pgs)}indivual"
    pgs["hm_chr"] = pgs["hm_chr"].astype(str)
    pgs["hm_pos"] = pgs["hm_pos"].astype(int)
    pgs["chr_name"] = pgs["chr_name"].astype(str)
    pgs["chr_position"] = pgs["chr_position"].astype(int)
    pgs["is_indel"] = (pgs["effect_allele"].str.len() > 1) | (pgs["other_allele"].str.len() > 1)

    def _key(r):
        try:
            return (int(r["chr_name"]), int(r["chr_position"]))
        except ValueError:
            return (None, None)

    keys = pgs.apply(_key, axis=1)
    pgs["is_chek2_tagsnp"] = keys.isin(CHEK2_TAGSNP_HG19)
    pgs["is_chek2_i157t"] = keys.map(lambda t: t == CHEK2_I157T_HG19)
    pgs["variant_uid"] = np.arange(len(pgs))

    n_indel = int(pgs["is_indel"].sum())
    n_snp = len(pgs) - n_indel
    print(f"Position 313: SNP={n_snp}, indel(effect or other allele length>1)={n_indel}")
    print(f"CHEK2 tag SNP (PRS311 exclusion, hg19 22:29203724/22:29551872):"
          f"{pgs.loc[pgs.is_chek2_tagsnp, ['chr_name','chr_position','hm_chr','hm_pos','effect_allele','other_allele']].to_string(index=False)}")
    print(f"CHEK2 p.I157T (mark only, PRS313/311 are reserved):"
          f"{pgs.loc[pgs.is_chek2_i157t, ['chr_name','chr_position','hm_chr','hm_pos','effect_allele','other_allele']].to_string(index=False)}")
    assert pgs["is_chek2_tagsnp"].sum() == 2, 'Should match exactly 2 CHEK2 tag SNPs'
    return pgs

def match_variants_to_aou(pgs):
    section('Step 4 | Chromosome-by-chromosome streaming matching AoU ACAF plink_bed .bim (match by chr:pos:ref:alt exact ID)')
    match_records = []
    for chrom in range(1, 23):
        sub = pgs[pgs.hm_chr == str(chrom)]
        if sub.empty:
            continue
        id_to_uids = {}
        for _, r in sub.iterrows():
            id_ref_is_other = f"chr{chrom}:{r.hm_pos}:{r.other_allele}:{r.effect_allele}"  # ref=other, alt=effect
            id_ref_is_effect = f"chr{chrom}:{r.hm_pos}:{r.effect_allele}:{r.other_allele}"  # ref=effect, alt=other
            id_to_uids.setdefault(id_ref_is_other, []).append(int(r.variant_uid))
            id_to_uids.setdefault(id_ref_is_effect, []).append(int(r.variant_uid))

        ids_file = f"{SCRATCH}/ids_chr{chrom}.txt"
        with open(ids_file, "w") as f:
            for k in id_to_uids:
                f.write(k + "\n")

        bim_remote = f"{BED_DIR}/acaf_threshold.chr{chrom}.bim"
        awk_cmd = (
            f"gsutil -u {GOOGLE_PROJECT} cat {bim_remote} | "
            f"awk -F'\\t' -v idsfile={ids_file} '"
            f"BEGIN {{ while ((getline line < idsfile) > 0) ids[line]=1 }} "
            f'{{ if ($2 in ids) print (NR-1)"\\t"$0 }}\''
        )
        t0 = time.time()
        res = subprocess.run(awk_cmd, shell=True, capture_output=True, text=True)
        elapsed = time.time() - t0
        n_lines = len([l for l in res.stdout.strip().split("\n") if l])
        print(f"  chr{chrom}: {len(sub)} target(s), {elapsed:.1f}s, {n_lines} candidate row(s) matched"
              + (f"  [stderr: {res.stderr[:200]}]" if res.stderr else ""))
        for line in res.stdout.strip().split("\n"):
            if not line:
                continue
            row_idx_s, bim_chrom, bim_id, cm, pos, a1, a2 = line.split("\t")
            for uid in id_to_uids.get(bim_id, []):
                match_records.append(dict(
                    variant_uid=uid, chrom=chrom, row_idx=int(row_idx_s),
                    bim_id=bim_id, pos=int(pos), A1=a1, A2=a2,
                ))
    match_df = pd.DataFrame(match_records)
    return match_df

def _orient_effect_allele(row):
    """ orient effect allele."""
    if row["effect_allele"] == row["A1"]:
        return "A1"
    if row["effect_allele"] == row["A2"]:
        return "A2"
    return "MISMATCH"

def resolve_matches(pgs, match_df):
    """Resolve matches."""
    section('Step 4b | Matching result disambiguation (clean / ambiguous / unmatched)')
    if len(match_df) == 0:
        clean = match_df.copy()
        ambiguous_uids, unmatched_uids = set(), set(pgs.variant_uid)
    else:
        n_distinct_rows = match_df.groupby("variant_uid")["row_idx"].nunique()
        clean_uids = set(n_distinct_rows[n_distinct_rows == 1].index)
        ambiguous_uids = set(n_distinct_rows[n_distinct_rows > 1].index)
        matched_uids = set(match_df["variant_uid"].unique())
        unmatched_uids = set(pgs["variant_uid"]) - matched_uids
        clean = match_df[match_df.variant_uid.isin(clean_uids)].drop_duplicates("variant_uid").copy()

    clean = clean.merge(
        pgs[["variant_uid", "effect_allele", "other_allele", "effect_weight", "is_indel", "is_chek2_tagsnp"]],
        on="variant_uid", how="left",
    )

    clean["effect_is"] = clean.apply(_orient_effect_allele, axis=1)
    n_mismatch = int((clean["effect_is"] == "MISMATCH").sum())
    if n_mismatch:
        print(f"WARNING: {n_mismatch}A clean matching allele is neither equal to A1 nor A2, and is eliminated.")
        clean = clean[clean["effect_is"] != "MISMATCH"].copy()

    print(f"clean matches:{len(clean)} / 313")
    print(f"Ambiguous (Multiple ref/alt representations hit at the same time, cannot be uniquely determined, AF comparison disambiguation will be tried):{len(ambiguous_uids)} / 313"
          + (f"  variant_uid={sorted(ambiguous_uids)}" if ambiguous_uids else ""))
    print(f"unmatched (not found at all in AoU ACAF callset):{len(unmatched_uids)} / 313"
          + (f"  variant_uid={sorted(unmatched_uids)}" if unmatched_uids else ""))
    return clean, ambiguous_uids, unmatched_uids

def load_fam_and_cohort_index(cohort_df):
    section('Step 5 | Get the sample order (.fam) and create a queue subset index')
    fam_local = f"{SCRATCH}/acaf_threshold.chr1.fam"
    if not os.path.exists(fam_local):
        run(f"gsutil -u {GOOGLE_PROJECT} cat {BED_DIR}/acaf_threshold.chr1.fam > {fam_local}")
    fam = pd.read_csv(fam_local, sep=r"\s+", header=None, names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"])
    N = len(fam)
    fam["IID"] = fam["IID"].astype("int64")
    print(f"The total number of ACAF callset samples N={N}(Including all AoU WGS participants, not just the 174,381 cohort of this project)")

    row_of_iid = pd.Series(np.arange(N), index=fam["IID"].values)
    in_bed = cohort_df["person_id"].isin(row_of_iid.index)
    n_missing = int((~in_bed).sum())
    if n_missing:
        print(f"WARNING: {n_missing}Queue person_id is not in the ACAF .fam sample list and will be removed from the PRS output")
    cohort_df = cohort_df[in_bed].reset_index(drop=True)
    cohort_row_idx = row_of_iid.reindex(cohort_df["person_id"].values).values.astype(int)
    print(f"Final queue size used for PRS calculation:{len(cohort_df)}")
    return N, cohort_row_idx, cohort_df

def _decode_bed_block(raw_bytes, N):
    byte_arr = np.frombuffer(raw_bytes, dtype=np.uint8)
    codes = np.zeros((len(byte_arr), 4), dtype=np.uint8)
    for k in range(4):
        codes[:, k] = (byte_arr >> (2 * k)) & 0b11
    codes = codes.reshape(-1)[:N]
    
    dosage_A1 = np.select([codes == 0, codes == 2, codes == 3], [2.0, 1.0, 0.0], default=np.nan)
    return dosage_A1

def fetch_and_decode_all(rows_df, N, key_col="variant_uid", label='site'):
    """Fetch and decode all."""
    section(f"Step 6 | Parallel byte-range precise crawling{len(rows_df)}indivual{label}Genotype byte chunks and decoding")
    record_size = -(-N // 4)  # ceil(N/4)
    print(f"record_size = ceil({N}/4) = {record_size} bytes/variant; "
          f"{len(rows_df)}indivual{label}Estimated total crawl volume ≈{len(rows_df)*record_size/1e6:.1f} MB")

    def worker(row):
        byte_start = 3 + row.row_idx * record_size
        byte_end = byte_start + record_size - 1
        bed_remote = f"{BED_DIR}/acaf_threshold.chr{row.chrom}.bed"
        cmd = f"gsutil -u {GOOGLE_PROJECT} cat -r {byte_start}-{byte_end} {bed_remote}"
        res = subprocess.run(cmd, shell=True, capture_output=True)
        key = getattr(row, key_col)
        if len(res.stdout) != record_size:
            print(f"WARNING {key_col}={key}: Number of bytes to grab{len(res.stdout)}≠Expectation{record_size}")
        dosage_A1 = _decode_bed_block(res.stdout, N)
        if row.effect_is == "A1":
            dosage_eff = dosage_A1
        else:  # "A2"
            dosage_eff = 2.0 - dosage_A1
        return key, dosage_eff

    results = {}
    t0 = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_FETCH_WORKERS) as ex:
        futs = [ex.submit(worker, row) for row in rows_df.itertuples()]
        for i, fut in enumerate(concurrent.futures.as_completed(futs)):
            key, dosage = fut.result()
            results[key] = dosage
            if (i + 1) % 50 == 0:
                print(f"  {i+1}/{len(rows_df)}Fetched, elapsed={time.time()-t0:.0f}s")
    print(f"All crawling completed:{len(results)}indivual{label}, time used{time.time()-t0:.0f}s")
    return results

def resolve_ambiguous_by_af(pgs, match_df, ambiguous_uids, N, cohort_row_idx):
    """Resolve ambiguous by af."""
    section('Step 6b | Constrained AF alignment disambiguation of 24 ambiguous indels (new in , R1-R4)')
    cfg = AF_DISAMBIG_CFG
    resolved_schema = ["variant_uid", "chrom", "row_idx", "bim_id", "pos", "A1", "A2",
                        "effect_allele", "other_allele", "effect_weight", "is_indel",
                        "is_chek2_tagsnp", "effect_is"]
    if not ambiguous_uids:
        return pd.DataFrame(columns=resolved_schema), pd.DataFrame()

    cand = match_df[match_df.variant_uid.isin(ambiguous_uids)].copy()
    cand = cand.merge(
        pgs[["variant_uid", "effect_allele", "other_allele", "effect_weight",
             "is_indel", "is_chek2_tagsnp", "allelefrequency_effect", "hm_rsID"]],
        on="variant_uid", how="left",
    )
    cand["effect_is"] = cand.apply(_orient_effect_allele, axis=1)
    n_mismatch = int((cand["effect_is"] == "MISMATCH").sum())
    if n_mismatch:
        print(f"WARNING: {n_mismatch}The allele of an ambiguous candidate row is neither equal to A1 nor A2, and the candidate row is eliminated.")
        cand = cand[cand["effect_is"] != "MISMATCH"].copy()
    cand = cand.reset_index(drop=True)
    cand["cand_key"] = cand.index  

    print(f"Number of ambiguous bits:{len(ambiguous_uids)}, total number of candidate rows (after removing MISMATCH):{len(cand)}")
    cand_group_sizes = cand.groupby("variant_uid").size()
    if (cand_group_sizes != 2).any():
        odd = cand_group_sizes[cand_group_sizes != 2]
        print(f"WARNING: {len(odd)}The number of candidate rows for an ambiguous site ≠ 2 (R3-3 mirror/complementary decision is only defined for exactly 2 candidates,"
              f"These sites will skip the R3-3 check):{odd.to_dict()}")

    geno_dict = fetch_and_decode_all(cand, N, key_col="cand_key", label='ambiguous candidate lines')

    af_rows = []
    for row in cand.itertuples():
        full_dosage = geno_dict[row.cand_key]
        cohort_dosage = full_dosage[cohort_row_idx]
        n_nonmissing = int(np.sum(~np.isnan(cohort_dosage)))
        af = float(np.nansum(cohort_dosage)) / (2 * n_nonmissing) if n_nonmissing > 0 else np.nan
        missing_rate = float(np.mean(np.isnan(cohort_dosage)))
        af_rows.append(dict(cand_key=row.cand_key, af_cohort=af,
                             missing_rate_cohort=missing_rate, n_nonmissing=n_nonmissing))
    cand = cand.merge(pd.DataFrame(af_rows), on="cand_key", how="left")

    diag_rows = []
    resolved_rows_list = []
    for uid, g in cand.groupby("variant_uid"):
        g = g.reset_index(drop=True)
        af_ref = g["allelefrequency_effect"].iloc[0]
        row_out = dict(variant_uid=int(uid), n_candidate_rows=len(g), af_ref=af_ref)

        d = (g["af_cohort"] - af_ref).abs() if pd.notna(af_ref) else pd.Series([np.nan] * len(g))
        order = d.sort_values(na_position="last").index
        g_sorted = g.loc[order].reset_index(drop=True)
        d_sorted = d.loc[order].reset_index(drop=True)

        for i in range(min(len(g_sorted), 2)):
            row_out[f"cand{i+1}_id"] = g_sorted.loc[i, "bim_id"]
            row_out[f"cand{i+1}_af"] = g_sorted.loc[i, "af_cohort"]
            row_out[f"cand{i+1}_missing_rate"] = g_sorted.loc[i, "missing_rate_cohort"]

        if pd.isna(af_ref):
            
            row_out.update(d_min=np.nan, d_2nd=np.nan, resolved_by="none",
                            resolved_row_id=None, high_missing_flag=False,
                            reject_reason="R2_af_ref_missing")
            diag_rows.append(row_out)
            continue

        d_min = float(d_sorted.iloc[0])
        row_out["d_min"] = d_min
        if len(d_sorted) < 2:
            row_out.update(d_2nd=np.nan, resolved_by="none", resolved_row_id=None,
                            high_missing_flag=False,
                            reject_reason="only_one_valid_candidate_no_discriminability_check")
            diag_rows.append(row_out)
            continue
        d_2nd = float(d_sorted.iloc[1])
        row_out["d_2nd"] = d_2nd

        degenerate = False
        if len(g_sorted) == 2:
            af_sum_dev = abs(g_sorted.loc[0, "af_cohort"] + g_sorted.loc[1, "af_cohort"] - 1.0)
            is_complementary = af_sum_dev <= cfg["complementary_af_tol"]
            is_af_ref_mid = cfg["degenerate_af_ref_low"] <= af_ref <= cfg["degenerate_af_ref_high"]
            degenerate = is_complementary and is_af_ref_mid
            row_out["af_sum_deviation_from_1"] = af_sum_dev
            row_out["is_complementary_pair"] = is_complementary

        pass_r1 = d_min <= cfg["d_min_max"]
        
        pass_ratio = True if d_min < cfg["d_min_ratio_floor"] else (d_2nd >= cfg["d_ratio_min"] * d_min)
        pass_gap = (d_2nd - d_min) >= cfg["d_gap_min"]
        pass_r2 = pass_gap and pass_ratio

        if degenerate:
            row_out.update(resolved_by="none", resolved_row_id=None, high_missing_flag=False,
                            reject_reason="R3-3_degenerate_mirror_zone")
        elif not pass_r1:
            row_out.update(resolved_by="none", resolved_row_id=None, high_missing_flag=False,
                            reject_reason="R3-1_d_min_exceeds_0.10")
        elif not pass_r2:
            row_out.update(resolved_by="none", resolved_row_id=None, high_missing_flag=False,
                            reject_reason="R3-2_insufficient_discriminative_power")
        else:
            accepted = g_sorted.iloc[0]
            high_missing = bool(accepted["missing_rate_cohort"] > cfg["max_missing_rate_flag"])
            row_out.update(resolved_by="af", resolved_row_id=accepted["bim_id"],
                            high_missing_flag=high_missing, reject_reason=None)
            resolved_rows_list.append(accepted)
        diag_rows.append(row_out)

    ambig_diag = pd.DataFrame(diag_rows)
    if resolved_rows_list:
        resolved_rows = pd.DataFrame(resolved_rows_list)[resolved_schema].reset_index(drop=True)
    else:
        resolved_rows = pd.DataFrame(columns=resolved_schema)

    n_resolved = len(resolved_rows)
    n_unresolved = len(ambiguous_uids) - n_resolved
    print(f"Disambiguation results:{n_resolved}/{len(ambiguous_uids)}ambiguous sites were successfully disambiguated (ambiguous_resolved_af),"
          f"{n_unresolved}remains skipped (ambiguous_unresolved)")
    if not ambig_diag.empty:
        reasons = ambig_diag.loc[ambig_diag["resolved_by"] != "af", "reject_reason"].value_counts()
        if len(reasons):
            print('Skip reason distribution:')
            print(reasons.to_string())
    n_high_missing = int(ambig_diag.get("high_missing_flag", pd.Series(dtype=bool)).sum())
    if n_high_missing:
        print(f"⚠️ {n_high_missing}Accepted locus deletion rate >{cfg['max_missing_rate_flag']:.0%}"
              f"(Still retained but marked with high_missing_flag, incorporated into existing mean-imputation sensitivity analysis)")

    return resolved_rows, ambig_diag

def _encode_bed_bytes(dosage_mat):
    """ encode bed bytes."""
    n_var, n_samp = dosage_mat.shape
    codes = np.ones(dosage_mat.shape, dtype=np.uint8)  
    codes[dosage_mat == 2] = 0
    codes[dosage_mat == 1] = 2
    codes[dosage_mat == 0] = 3
    record_size = -(-n_samp // 4)
    pad = (-n_samp) % 4
    if pad:
        codes = np.concatenate([codes, np.zeros((n_var, pad), dtype=np.uint8)], axis=1)
    codes = codes.reshape(n_var, -1, 4)
    byte_vals = (codes[:, :, 0] | (codes[:, :, 1] << 2) | (codes[:, :, 2] << 4) | (codes[:, :, 3] << 6)).astype(np.uint8)
    out = bytearray(3 + n_var * record_size)
    out[0:3] = bytes([0x6C, 0x1B, 0x01])
    out[3:] = byte_vals.tobytes()
    return bytes(out), record_size

def plink2_cross_check(pgs, clean, geno_dict, cohort_df, cohort_row_idx):
    section('Step 7 | Rebuild small .bed/.bim/.fam, cross-validate with real plink2 --score')
    uids_sorted = sorted(geno_dict.keys())
    meta = clean.set_index("variant_uid").loc[uids_sorted]
    mat_full = np.vstack([geno_dict[u] for u in uids_sorted])
    mat_cohort = mat_full[:, cohort_row_idx]  # (n_matched, n_cohort), effect-allele dosage

    mini_prefix = f"{SCRATCH}/mini_prs"
    bed_bytes, record_size = _encode_bed_bytes(mat_cohort)
    with open(f"{mini_prefix}.bed", "wb") as f:
        f.write(bed_bytes)
    print(f"mini.bed writes:{len(bed_bytes)/1e6:.1f} MB ({mat_cohort.shape[0]} variants x {mat_cohort.shape[1]} samples)")

    var_ids = [f"prs_{u}" for u in uids_sorted]
    bim = pd.DataFrame({
        "chrom": meta["chrom"].values,
        "id": var_ids,
        "cm": 0,
        "pos": meta["pos"].values,
        "A1": meta["effect_allele"].values,   
        "A2": meta["other_allele"].values,
    })
    bim.to_csv(f"{mini_prefix}.bim", sep="\t", header=False, index=False)

    fam = pd.DataFrame({
        "FID": 0, "IID": cohort_df["person_id"].values,
        "PAT": 0, "MAT": 0, "SEX": 0, "PHENO": -9,
    })
    fam.to_csv(f"{mini_prefix}.fam", sep="\t", header=False, index=False)

    is_tagsnp = meta["is_chek2_tagsnp"].values
    score313 = pd.DataFrame({"id": var_ids, "A1": meta["effect_allele"].values, "weight": meta["effect_weight"].values})
    score313.to_csv(f"{SCRATCH}/score313.txt", sep="\t", header=True, index=False)
    score311 = score313[~is_tagsnp]
    score311.to_csv(f"{SCRATCH}/score311.txt", sep="\t", header=True, index=False)

    plink_results = {}
    for tag, scorefile in [("313", "score313.txt"), ("311", "score311.txt")]:
        for imp_tag, extra_flag in [("meanimp", ""), ("noimp", "no-mean-imputation")]:
            out_prefix = f"{SCRATCH}/plink_out_{tag}_{imp_tag}"
            cmd = (
                f"{PLINK2} --bfile {mini_prefix} --score {SCRATCH}/{scorefile} 1 2 3 header cols=+scoresums "
                f"{extra_flag} --out {out_prefix}"
            )
            res = run(cmd)
            sscore_path = f"{out_prefix}.sscore"
            if not os.path.exists(sscore_path):
                print(f"WARNING: plink2 not generated{sscore_path}, stderr is as above, this cross-validation is skipped")
                continue
            sscore = pd.read_csv(sscore_path, sep="\t")
            iid_col = "#IID" if "#IID" in sscore.columns else "IID"
            sscore = sscore.rename(columns={iid_col: "person_id"})

            score_col = [c for c in sscore.columns if c.startswith("SCORE") and c.endswith("_SUM")][0]
            plink_results[(tag, imp_tag)] = sscore[["person_id", score_col]].rename(columns={score_col: f"PRS{tag}_{imp_tag}_plink2"})

    return plink_results, mat_cohort, uids_sorted, meta

def compute_prs_numpy(mat_cohort, uids_sorted, meta, cohort_df):
    section('Step 8 | numpy self-calculates 313/311 x mean-imputation/no-mean-imputation (cross-validation plink2 output)')
    missing_mask = np.isnan(mat_cohort)
    missing_rate = missing_mask.mean(axis=1)
    col_mean = np.nanmean(mat_cohort, axis=1)
    mean_imputed = np.where(missing_mask, col_mean[:, None], mat_cohort)
    zero_imputed = np.nan_to_num(mat_cohort, nan=0.0)

    weights = meta["effect_weight"].values.astype(float)
    is_tagsnp = meta["is_chek2_tagsnp"].values
    w311 = weights.copy()
    w311[is_tagsnp] = 0.0

    out = pd.DataFrame({"person_id": cohort_df["person_id"].values})
    out["PRS313_meanimp_np"] = weights @ mean_imputed
    out["PRS313_noimp_np"] = weights @ zero_imputed
    out["PRS311_meanimp_np"] = w311 @ mean_imputed
    out["PRS311_noimp_np"] = w311 @ zero_imputed

    diag = pd.DataFrame({
        "variant_uid": uids_sorted,
        "missing_rate_in_cohort": missing_rate,
    })
    return out, diag

def common_variants_calculation():
    download_scoring_files()
    cohort_df = load_cohort()
    pgs = build_target_table()
    match_df = match_variants_to_aou(pgs)
    clean, ambiguous_uids, unmatched_uids = resolve_matches(pgs, match_df)

    N, cohort_row_idx, cohort_df = load_fam_and_cohort_index(cohort_df)

    resolved_rows, ambig_diag = resolve_ambiguous_by_af(pgs, match_df, ambiguous_uids, N, cohort_row_idx)
    clean_extended = pd.concat([clean, resolved_rows], ignore_index=True) if len(resolved_rows) else clean

    geno_dict = fetch_and_decode_all(clean_extended, N)

    plink_results, mat_cohort, uids_sorted, meta = plink2_cross_check(pgs, clean_extended, geno_dict, cohort_df, cohort_row_idx)
    np_scores, variant_diag = compute_prs_numpy(mat_cohort, uids_sorted, meta, cohort_df)

    clean_mask = np.isin(np.asarray(uids_sorted), clean["variant_uid"].values)
    mm = mat_cohort[clean_mask]
    mm = np.where(np.isnan(mm), np.nanmean(mm, axis=1)[:, None], mm)
    w = meta["effect_weight"].values.astype(float)[clean_mask]
    w311 = np.where(meta["is_chek2_tagsnp"].values[clean_mask], 0.0, w)
    clean_scores = pd.DataFrame({"person_id": cohort_df["person_id"].values,
                                 "PRS313": w @ mm, "PRS311": w311 @ mm})
    clean_path = f"{DATA_DIR}/prs_scores_281loci.csv"
    clean_scores.to_csv(clean_path, index=False)
    print(f"Clean only PRS of matching sites:{int(clean_mask.sum())}single point ->{clean_path} ({len(clean_scores)}OK)")

    section('Step 9 | plink2 vs numpy cross check')
    merged = np_scores.copy()
    for (tag, imp_tag), df in plink_results.items():
        merged = merged.merge(df, on="person_id", how="left")

    for tag in ["313", "311"]:
        for imp_tag in ["meanimp", "noimp"]:
            np_col = f"PRS{tag}_{imp_tag}_np"
            pk_col = f"PRS{tag}_{imp_tag}_plink2"
            if pk_col in merged.columns:
                diff = (merged[np_col] - merged[pk_col]).abs()
                print(f"PRS{tag} {imp_tag}: numpy vs plink2 max|diff|={diff.max():.6g}, "
                      f"corr={merged[np_col].corr(merged[pk_col]):.6f}")
            else:
                print(f"PRS{tag} {imp_tag}: plink2 cross-validation results are missing, only numpy version reported")

    section('Step 10 | Write the main output data/prs_scores.csv')
    final = merged[["person_id"]].copy()
    
    for tag in ["313", "311"]:
        pk_col = f"PRS{tag}_meanimp_plink2"
        np_col = f"PRS{tag}_meanimp_np"
        final[f"PRS{tag}"] = merged[pk_col] if pk_col in merged.columns else merged[np_col]
        pk_col_noimp = f"PRS{tag}_noimp_plink2"
        np_col_noimp = f"PRS{tag}_noimp_np"
        final[f"PRS{tag}_noimp"] = merged[pk_col_noimp] if pk_col_noimp in merged.columns else merged[np_col_noimp]

    final = cohort_df[["person_id"]].merge(final, on="person_id", how="left")
    out_path = f"{DATA_DIR}/prs_scores.csv"
    final.to_csv(out_path, index=False)
    size_mb = os.path.getsize(out_path) / 1e6
    print(f"Write:{out_path}  ({len(final)}OK,{size_mb:.2f}MB, column={list(final.columns)})")

    if size_mb >= 5:
        print(f"  NOTE: {out_path}volume{size_mb:.2f}MB is more than 5MB, but the entire data/ directory is already in .gitignore,"
              f"Never enter git push, do not violate the push volume rules of CLAUDE.md, only record without interruption.")

    section('Step 11 | Write a site hit rate diagnosis table (: add ambiguous AF disambiguation column)')
    resolved_uids = set(resolved_rows["variant_uid"]) if len(resolved_rows) else set()
    unresolved_ambiguous_uids = ambiguous_uids - resolved_uids

    pgs_diag = pgs.copy()
    pgs_diag["match_status"] = "unmatched"
    pgs_diag.loc[pgs_diag.variant_uid.isin(unresolved_ambiguous_uids), "match_status"] = "ambiguous_unresolved"
    pgs_diag.loc[pgs_diag.variant_uid.isin(resolved_uids), "match_status"] = "ambiguous_resolved_af"
    pgs_diag.loc[pgs_diag.variant_uid.isin(clean.variant_uid), "match_status"] = "matched"
    pgs_diag = pgs_diag.merge(variant_diag, on="variant_uid", how="left")

    ambig_diag_cols = [
        "variant_uid", "n_candidate_rows", "cand1_id", "cand1_af", "cand1_missing_rate",
        "cand2_id", "cand2_af", "cand2_missing_rate", "af_ref", "d_min", "d_2nd",
        "is_complementary_pair", "af_sum_deviation_from_1", "resolved_by", "resolved_row_id",
        "high_missing_flag", "reject_reason",
    ]
    if not ambig_diag.empty:
        pgs_diag = pgs_diag.merge(
            ambig_diag[[c for c in ambig_diag_cols if c in ambig_diag.columns]],
            on="variant_uid", how="left",
        )
    for c in ambig_diag_cols:
        if c != "variant_uid" and c not in pgs_diag.columns:
            pgs_diag[c] = np.nan

    diag_cols = [
        "variant_uid", "chr_name", "chr_position", "hm_chr", "hm_pos",
        "effect_allele", "other_allele", "effect_weight", "is_indel",
        "is_chek2_tagsnp", "is_chek2_i157t", "match_status", "missing_rate_in_cohort",
    ] + ambig_diag_cols[1:]
    diag_out_path = f"{COMMON_VARIANTS_TABLE_DIR}/07_prs313_variant_match_diagnostics.csv"

    n_c = len(cohort_df)
    pub = pgs_diag[diag_cols].copy()
    n_banded = 0
    for col, miss_col, min_count, both_tails in [
            ("missing_rate_in_cohort", None, 20, False),
            ("cand1_missing_rate", None, 20, False), ("cand2_missing_rate", None, 20, False),
            ("cand1_af", "cand1_missing_rate", 40, True), ("cand2_af", "cand2_missing_rate", 40, True)]:
        v = pd.to_numeric(pub[col], errors="coerce")
        if miss_col is None:
            denom = pd.Series(float(n_c), index=v.index)
        else:
            denom = 2 * n_c * (1 - pd.to_numeric(pgs_diag[miss_col], errors="coerce").fillna(0))
        lo = (v * denom).round()
        small = (lo > 0) & (lo < min_count)
        if both_tails:
            hi = ((1 - v) * denom).round()
            small_hi = (hi > 0) & (hi < min_count)
        else:
            small_hi = pd.Series(False, index=v.index)
        pub[col] = pub[col].astype(object)
        pub.loc[small, col] = [f"<{min_count / d:.1e}" for d in denom[small]]
        pub.loc[small_hi, col] = [f">{1 - min_count / d:.6f}" for d in denom[small_hi]]
        n_banded += int(small.sum() + small_hi.sum())
    print(f"Disclosure Processing:{n_banded}The proportion grid is rewritten as a gear (the number of people in reverse calculation falls in 1-19 / the number of alleles falls in 1-39)")
    pub.to_csv(diag_out_path, index=False)
    print(f"Write:{diag_out_path}  ({len(pgs_diag)}OK)")

    section('Step 12 | Summary Statistics')
    n_total = len(pgs_diag)
    n_matched = int((pgs_diag.match_status == "matched").sum())
    n_ambig_resolved = int((pgs_diag.match_status == "ambiguous_resolved_af").sum())
    n_ambig_unresolved = int((pgs_diag.match_status == "ambiguous_unresolved").sum())
    n_ambiguous = n_ambig_resolved + n_ambig_unresolved
    n_unmatched = int((pgs_diag.match_status == "unmatched").sum())
    n_used = n_matched + n_ambig_resolved  
    n_indel_total = int(pgs_diag.is_indel.sum())
    n_snp_total = n_total - n_indel_total
    n_used_mask = pgs_diag.match_status.isin(["matched", "ambiguous_resolved_af"])
    n_indel_used = int((n_used_mask & pgs_diag.is_indel).sum())
    n_snp_used = int((n_used_mask & (~pgs_diag.is_indel)).sum())

    print(f"Total number of points:{n_total}")
    print(f"matched (native clean):{n_matched}ambiguous_resolved_af (new disambiguation recovery):{n_ambig_resolved}  "
          f"ambiguous_unresolved (still skipped):{n_ambig_unresolved}  unmatched: {n_unmatched}")
    print(f"The total number of sites actually entered into the PRS calculation:{n_used}/{n_total} ({n_used/n_total:.1%})")
    print(f"SNP: {n_snp_used}/{n_snp_total} used ({n_snp_used/n_snp_total:.1%})")
    print(f"indel: {n_indel_used}/{n_indel_total} used ({n_indel_used/n_indel_total:.1%})")

    for tag in ["313", "311"]:
        col = f"PRS{tag}"
        print(f"{col}: mean={final[col].mean():.4f}, sd={final[col].std():.4f}, "
              f"min={final[col].min():.4f}, max={final[col].max():.4f}")

    for tag in ["313", "311"]:
        meanc, noimpc = f"PRS{tag}_meanimp_np", f"PRS{tag}_noimp_np"
        corr = merged[meanc].corr(merged[noimpc])
        diff = (merged[meanc] - merged[noimpc])
        print(f"PRS{tag} mean-imputation vs no-mean-imputation: corr={corr:.6f}, "
              f"mean|diff|={diff.abs().mean():.6g}, max|diff|={diff.abs().max():.6g}")

    summary = dict(
        n_total=n_total, n_matched=n_matched,
        n_ambig_resolved=n_ambig_resolved, n_ambig_unresolved=n_ambig_unresolved,
        n_ambiguous=n_ambiguous, n_unmatched=n_unmatched, n_used=n_used,
        n_snp_total=n_snp_total, n_snp_used=n_snp_used,
        n_indel_total=n_indel_total, n_indel_used=n_indel_used,
        cohort_n=len(final),
    )
    return summary, merged, final, pgs_diag

def calculate_common_variant_prs():
    """Initialize runtime state and execute this component explicitly."""
    global AF_DISAMBIG_CFG, AIM1_CONFIG, BED_DIR, CDR_STORAGE_PATH, CHEK2_I157T_HG19, CHEK2_TAGSNP_HG19, CONFIG_PATH, DATA_DIR, GOOGLE_PROJECT, HM37_URL, HM38_URL, MAX_FETCH_WORKERS, PGS_ID, PLINK2, PRS_RAW_DIR, REPO, REPORT_DIR, SCRATCH, COMMON_VARIANTS_TABLE_DIR, _f, bigquery, d

    from google.cloud import bigquery

    GOOGLE_PROJECT = os.environ["GOOGLE_PROJECT"]
    CDR_STORAGE_PATH = os.environ["CDR_STORAGE_PATH"]

    REPO = "/home/jupyter/workspace/analysis"
    DATA_DIR = f"{REPO}/data"
    PRS_RAW_DIR = f"{DATA_DIR}/prs_raw"
    COMMON_VARIANTS_TABLE_DIR = f"{REPO}/results/table"
    REPORT_DIR = f"{REPO}/results/report"
    CONFIG_PATH = f"{REPO}/config/aim1_config.yaml"
    SCRATCH = "/tmp/claude-1001/-home-jupyter-workspace-analysis/4bfcbf4d-1cb3-4b4d-b35e-84140637b300/scratchpad/prs_scratch"

    for d in (PRS_RAW_DIR, COMMON_VARIANTS_TABLE_DIR, REPORT_DIR, SCRATCH):
        os.makedirs(d, exist_ok=True)

    with open(CONFIG_PATH) as _f:
        AIM1_CONFIG = yaml.safe_load(_f)
    AF_DISAMBIG_CFG = AIM1_CONFIG["prs_indel_af_disambiguation"]

    BED_DIR = f"{CDR_STORAGE_PATH}/wgs/short_read/snpindel/acaf_threshold/plink_bed"
    PLINK2 = "/opt/workbench-tools/binaries/bin/plink2"

    PGS_ID = "PGS000004"
    HM38_URL = f"https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{PGS_ID}/ScoringFiles/Harmonized/{PGS_ID}_hmPOS_GRCh38.txt.gz"
    HM37_URL = f"https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{PGS_ID}/ScoringFiles/Harmonized/{PGS_ID}_hmPOS_GRCh37.txt.gz"

    CHEK2_TAGSNP_HG19 = {(22, 29203724), (22, 29551872)}
    
    CHEK2_I157T_HG19 = (22, 29121087)

    MAX_FETCH_WORKERS = 16

    if True:
        common_variants_calculation()

# ============================================================================
# 2. match-audit: audit_prs_matches
# ============================================================================

import os
import pandas as pd

def audit_prs_matches():
    """Initialize runtime state and execute this component explicitly."""
    global BASE, PRS_MATCH_AUDIT_TABLE_DIR, consistency_ok, dbsnp_column_finding, diag, n_ambiguous, n_failed_only_one_candidate, n_failed_r2_missing_af, n_failed_r31, n_failed_r32_discrim, n_failed_r33_degenerate, n_flagged_high_missing, n_matched_clean, n_passed_r32_via_waiver, n_resolved_by_af, n_resolved_by_dbsnp, n_total, n_unmatched, n_unresolved, out, out_path, reject_reason, resolved_af_mask, resolved_by_values

    BASE = "/home/jupyter/workspace/analysis"
    PRS_MATCH_AUDIT_TABLE_DIR = os.path.join(BASE, "results", "table")

    diag = pd.read_csv(os.path.join(PRS_MATCH_AUDIT_TABLE_DIR, "07_prs313_variant_match_diagnostics.csv"))
    print(f"Read diagnostic table:{diag.shape}, column={diag.columns.tolist()}")

    n_total = len(diag)
    n_ambiguous = int((diag["match_status"].isin(["ambiguous_resolved_af", "ambiguous_unresolved"])).sum())
    n_matched_clean = int((diag["match_status"] == "matched").sum())
    n_unmatched = int((diag["match_status"] == "unmatched").sum())

    reject_reason = diag["reject_reason"]
    n_failed_r2_missing_af = int((reject_reason == "R2_af_ref_missing").sum())
    n_failed_r31 = int((reject_reason == "R3-1_d_min_exceeds_0.10").sum())
    n_failed_r32_discrim = int((reject_reason == "R3-2_insufficient_discriminative_power").sum())
    n_failed_r33_degenerate = int((reject_reason == "R3-3_degenerate_mirror_zone").sum())
    n_failed_only_one_candidate = int((reject_reason == "only_one_valid_candidate_no_discriminability_check").sum())

    resolved_af_mask = diag["resolved_by"] == "af"
    n_resolved_by_af = int(resolved_af_mask.sum())

    if "d_min" in diag.columns:
        n_passed_r32_via_waiver = int((resolved_af_mask & (diag["d_min"] < 0.01)).sum())
    else:
        n_passed_r32_via_waiver = None

    n_flagged_high_missing = int((diag.get("high_missing_flag") == True).sum()) if "high_missing_flag" in diag.columns else None  # noqa: E712

    if "resolved_by" in diag.columns:
        resolved_by_values = set(diag["resolved_by"].dropna().unique())
        n_resolved_by_dbsnp = int((diag["resolved_by"] == "dbsnp").sum())
        dbsnp_column_finding = (
            "The resolved_by column exists and has been read, and the measured value set = {} -- does not include 'dbsnp'. This is consistent with the script header comment:"
            'R6 (dbSNP/Ensembl orthogonality) has never been implemented in (not missing, but the disambiguation path itself is not implemented/run),'
            'This 0 value is reported truthfully, not back-fill.'.format(sorted(resolved_by_values))
        )
    else:
        n_resolved_by_dbsnp = None
        dbsnp_column_finding = "The resolved_by column itself does not exist in the diagnostic table (a more serious deletion than 'no dbsnp in the value')."

    n_unresolved = int((diag["match_status"] == "ambiguous_unresolved").sum())

    print("\n--- Gate-by-gate attrition ---")
    print(f"n_ambiguous (total number of ambiguous sites to be disambiguated):{n_ambiguous}")
    print(f"n_matched_clean (non-ambiguous, native clean matching):{n_matched_clean}")
    print(f"n_unmatched (not found at all in AoU ACAF callset):{n_unmatched}")
    print(f"n_failed_R2_missing_AF (if_ref is missing, discard it):{n_failed_r2_missing_af}")
    print(f"n_failed_only_one_valid_candidate (only 1 valid candidate, unable to do discriminant test and skipped):{n_failed_only_one_candidate}")
    print(f"n_failed_R3-1 (d_min>0.10): {n_failed_r31}")
    print(f"n_failed_R3-2 (discrimination): {n_failed_r32_discrim}")
    print(f"n_passed_R3-2_via_d_min<0.01_auto_waiver (among accepted sites, passed d_min<0.01 exemption discrimination test):"
          f"{n_passed_r32_via_waiver}")
    print(f"n_failed_R3-3 (degenerate mirror zone): {n_failed_r33_degenerate}")
    print(f"n_flagged_high_missing (accepted but missing rate > threshold, retained after marking):{n_flagged_high_missing}")
    print(f"n_resolved_by=af: {n_resolved_by_af}")
    print(f"n_resolved_by=dbsnp (R6): {n_resolved_by_dbsnp}")
    print(dbsnp_column_finding)
    print(f"n_unresolved (still skipped, does not enter PRS calculation):{n_unresolved}")

    consistency_ok = n_ambiguous == (n_resolved_by_af + n_unresolved)
    print(f"\n self-consistency check: n_ambiguous({n_ambiguous}) == n_resolved_by_af({n_resolved_by_af}) + "
          f"n_unresolved({n_unresolved}) -> {'PASS' if consistency_ok else 'FAIL'}")

    out = pd.DataFrame([{
        "n_total_loci": n_total,
        "n_matched_clean": n_matched_clean,
        "n_unmatched": n_unmatched,
        "n_ambiguous": n_ambiguous,
        "n_failed_R2_missing_AF": n_failed_r2_missing_af,
        "n_failed_only_one_valid_candidate_no_discriminability_check": n_failed_only_one_candidate,
        "n_failed_R3-1_d_min_gt_0.10": n_failed_r31,
        "n_failed_R3-2_discrimination": n_failed_r32_discrim,
        "n_passed_R3-2_via_d_min_lt_0.01_auto_waiver": n_passed_r32_via_waiver,
        "n_failed_R3-3_degenerate_mirror_zone": n_failed_r33_degenerate,
        "n_flagged_high_missing": n_flagged_high_missing,
        "n_resolved_by=af": n_resolved_by_af,
        "n_resolved_by=dbsnp": n_resolved_by_dbsnp,
        "n_resolved_by=dbsnp_finding": dbsnp_column_finding,
        "n_unresolved": n_unresolved,
        "consistency_check_n_ambiguous_eq_resolved_plus_unresolved": consistency_ok,
    }])
    out_path = os.path.join(PRS_MATCH_AUDIT_TABLE_DIR, "07_r3_gate_attrition.csv")
    out.to_csv(out_path, index=False)
    print(f"\n writes:{out_path}")

# Command-line selection of the computations defined above.
def main():
    from scripts.functions.cli import run_local_step
    computations = {
        'prs': calculate_common_variant_prs,
        'match-audit': audit_prs_matches,
    }
    run_local_step('Common SNP harmonization, selection and PRS calculation', computations, default_tasks=['prs', 'match-audit'])

if __name__ == '__main__':
    main()

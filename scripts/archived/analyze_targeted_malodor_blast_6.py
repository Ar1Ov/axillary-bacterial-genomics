import re
from pathlib import Path

import numpy as np
import pandas as pd

# Six-reference targeted malodor homology analysis.
# This preserves the original best-hit logic and descriptive thresholds.
# Only the BLAST input/output names and QUERY_MAP/ref list are expanded.

BLAST_FILE = "malodor_targeted_blast_6_out.txt"
PANEL_FILE = "malodor_targeted_genomes.tsv"

OUT_ALL = "malodor_blast_hits_annotated_6.tsv"
OUT_BEST = "malodor_best_hit_per_genome_6.tsv"
OUT_RELEVANT = "malodor_best_hits_relevant_genus_6.tsv"
OUT_SUMMARY = "malodor_homology_group_summary_6.tsv"
OUT_SUMMARY_ALL = "malodor_homology_group_summary_all_genomes_6.tsv"
OUT_TOP = "malodor_top_hits_by_reference_6.tsv"
OUT_TOP_ALL = "malodor_top_hits_by_reference_all_genomes_6.tsv"

for f in [BLAST_FILE, PANEL_FILE]:
    if not Path(f).exists():
        raise FileNotFoundError(f)

cols = [
    "qseqid","sseqid","pident","length","mismatch","gapopen",
    "qstart","qend","sstart","send","evalue","bitscore","qlen","slen"
]
blast = pd.read_csv(
    BLAST_FILE, sep="\t", header=None, names=cols,
    dtype=str, keep_default_na=False, low_memory=False
)
panel = pd.read_csv(
    PANEL_FILE, sep="\t", dtype=str, keep_default_na=False, low_memory=False
)

required_panel = {
    "genome_id","species","bioproject","site","source_group","analysis_roles"
}
missing = required_panel - set(panel.columns)
if missing:
    raise ValueError(f"Panel file missing columns: {sorted(missing)}")
if panel["genome_id"].duplicated().any():
    raise ValueError("Duplicate genome IDs in malodor panel.")

QUERY_MAP = {
    "WP_119633472.1": {
        "reference":"PatB","reference_accession":"WP_119633472.1",
        "relevant_genus":"Staphylococcus",
        "biological_role":"3M3SH-releasing C-S lyase",
    },
    "sp|Q8GGD4|AGAA_CORST": {
        "reference":"AgaA","reference_accession":"Q8GGD4",
        "relevant_genus":"Corynebacterium",
        "biological_role":"N(alpha)-acyl-glutamine aminoacylase",
    },
    "sp|B2KZE7|TPDA_CORST": {
        "reference":"TpdA","reference_accession":"B2KZE7",
        "relevant_genus":"Corynebacterium",
        "biological_role":"Cysteinylglycine-S-conjugate dipeptidase",
    },
    "tr|Q4JWQ6|Q4JWQ6_CORJK": {
        "reference":"AecD","reference_accession":"Q4JWQ6",
        "relevant_genus":"Corynebacterium",
        "biological_role":"Cysteine-S-conjugate beta-lyase",
    },
    "WP_002449073.1": {
        "reference":"DtpT","reference_accession":"WP_002449073.1",
        "relevant_genus":"Staphylococcus",
        "biological_role":"S-Cys-Gly-3M3SH uptake peptide MFS transporter",
    },
    "A0A657M1L9": {
        "reference":"PepV","reference_accession":"A0A657M1L9",
        "relevant_genus":"Staphylococcus",
        "biological_role":"Cys-Gly-3M3SH-processing dipeptidase",
    },
}

unexpected = sorted(set(blast["qseqid"]) - set(QUERY_MAP))
if unexpected:
    raise ValueError(f"Unexpected BLAST query IDs: {unexpected}")

for c in ["pident","length","qstart","qend","sstart","send","evalue","bitscore","qlen","slen"]:
    blast[c] = pd.to_numeric(blast[c], errors="coerce")

if blast[["pident","length","bitscore","qlen","slen"]].isna().any().any():
    raise ValueError("Could not parse one or more required BLAST numeric fields.")

pat = re.compile(r"^fig\|(.+)\.peg\.\d+$")
blast["genome_id"] = blast["sseqid"].map(
    lambda s: pat.match(str(s)).group(1) if pat.match(str(s)) else ""
)
if (blast["genome_id"] == "").any():
    bad = blast.loc[blast["genome_id"] == "", "sseqid"].head(10).tolist()
    raise ValueError(f"Could not parse genome ID from subject IDs: {bad}")

for field in ["reference","reference_accession","relevant_genus","biological_role"]:
    blast[field] = blast["qseqid"].map(lambda q: QUERY_MAP[q][field])

blast["query_coverage_pct"] = (100.0 * blast["length"] / blast["qlen"]).clip(upper=100)
blast["subject_coverage_pct"] = (100.0 * blast["length"] / blast["slen"]).clip(upper=100)

blast = blast.merge(panel, on="genome_id", how="left", validate="many_to_one")
if blast["species"].eq("").any():
    bad = blast.loc[blast["species"].eq(""), ["sseqid","genome_id"]].head(10)
    raise ValueError("BLAST hit(s) outside fixed panel:\n" + bad.to_string(index=False))

blast["genus"] = blast["species"].str.split().str[0]
blast["relevant_genus_match"] = blast["genus"].eq(blast["relevant_genus"])

def similarity_band(row):
    pid = row["pident"]
    qcov = row["query_coverage_pct"]
    if pid >= 70 and qcov >= 80:
        return "HIGH_SIMILARITY"
    if pid >= 50 and qcov >= 80:
        return "MODERATE_SIMILARITY"
    if pid >= 30 and qcov >= 70:
        return "DISTANT_FULL_LENGTH_HOMOLOGY"
    return "WEAK_OR_PARTIAL"

blast["similarity_band"] = blast.apply(similarity_band, axis=1)
blast.to_csv(OUT_ALL, sep="\t", index=False, encoding="utf-8-sig")

best_hit = (
    blast.sort_values(
        ["reference","genome_id","bitscore","query_coverage_pct","pident","evalue"],
        ascending=[True,True,False,False,False,True]
    )
    .drop_duplicates(["reference","genome_id"], keep="first")
    .copy()
)

refs = ["PatB","AgaA","TpdA","AecD","DtpT","PepV"]
matrix = pd.MultiIndex.from_product(
    [panel["genome_id"].tolist(), refs],
    names=["genome_id","reference"]
).to_frame(index=False)

ref_meta = pd.DataFrame([
    {
        "reference":v["reference"],
        "reference_accession":v["reference_accession"],
        "relevant_genus":v["relevant_genus"],
        "biological_role":v["biological_role"],
    }
    for v in QUERY_MAP.values()
]).drop_duplicates("reference")

matrix = matrix.merge(panel, on="genome_id", how="left", validate="many_to_one")
matrix = matrix.merge(ref_meta, on="reference", how="left", validate="many_to_one")

keep_best = [
    "genome_id","reference","sseqid","pident","length","evalue","bitscore",
    "qlen","slen","query_coverage_pct","subject_coverage_pct","similarity_band"
]
matrix = matrix.merge(
    best_hit[keep_best],
    on=["genome_id","reference"],
    how="left",
    validate="one_to_one"
)

matrix["detected_hit"] = matrix["sseqid"].notna()
matrix["genus"] = matrix["species"].str.split().str[0]
matrix["relevant_genus_match"] = matrix["genus"].eq(matrix["relevant_genus"])
matrix["similarity_band"] = matrix["similarity_band"].fillna("NO_DETECTABLE_HIT")

# This all-genome table is the main source for cross-genus figures and phylogeny selection.
matrix.to_csv(OUT_BEST, sep="\t", index=False, encoding="utf-8-sig")

# Retain the original relevant-genus subset for direct comparability with the 4-protein analysis.
relevant = matrix[matrix["relevant_genus_match"]].copy()
relevant.to_csv(OUT_RELEVANT, sep="\t", index=False, encoding="utf-8-sig")

def make_summary(df):
    summary_rows = []
    for (ref,species,group), g in df.groupby(
        ["reference","species","source_group"], dropna=False
    ):
        hit = g[g["detected_hit"]].copy()
        row = {
            "reference":ref,
            "species":species,
            "source_group":group,
            "genomes":len(g),
            "genomes_with_detectable_hit":int(g["detected_hit"].sum()),
            "detectable_hit_fraction":float(g["detected_hit"].mean()),
            "high_similarity_n":int((g["similarity_band"]=="HIGH_SIMILARITY").sum()),
            "moderate_similarity_n":int((g["similarity_band"]=="MODERATE_SIMILARITY").sum()),
            "distant_full_length_n":int((g["similarity_band"]=="DISTANT_FULL_LENGTH_HOMOLOGY").sum()),
            "weak_partial_n":int((g["similarity_band"]=="WEAK_OR_PARTIAL").sum()),
            "no_detectable_hit_n":int((g["similarity_band"]=="NO_DETECTABLE_HIT").sum()),
        }
        if len(hit):
            row.update({
                "median_best_identity_pct":float(hit["pident"].median()),
                "min_best_identity_pct":float(hit["pident"].min()),
                "max_best_identity_pct":float(hit["pident"].max()),
                "median_query_coverage_pct":float(hit["query_coverage_pct"].median()),
                "median_bitscore":float(hit["bitscore"].median()),
            })
        else:
            row.update({
                "median_best_identity_pct":np.nan,
                "min_best_identity_pct":np.nan,
                "max_best_identity_pct":np.nan,
                "median_query_coverage_pct":np.nan,
                "median_bitscore":np.nan,
            })
        summary_rows.append(row)

    return pd.DataFrame(summary_rows).sort_values(
        ["reference","species","source_group"]
    )

summary = make_summary(relevant)
summary.to_csv(OUT_SUMMARY, sep="\t", index=False, encoding="utf-8-sig")

summary_all = make_summary(matrix)
summary_all.to_csv(OUT_SUMMARY_ALL, sep="\t", index=False, encoding="utf-8-sig")

top = (
    relevant[relevant["detected_hit"]]
    .sort_values(
        ["reference","pident","query_coverage_pct","bitscore"],
        ascending=[True,False,False,False]
    )
    .groupby("reference", group_keys=False)
    .head(25)
    .copy()
)
top.to_csv(OUT_TOP, sep="\t", index=False, encoding="utf-8-sig")

top_all = (
    matrix[matrix["detected_hit"]]
    .sort_values(
        ["reference","pident","query_coverage_pct","bitscore"],
        ascending=[True,False,False,False]
    )
    .groupby("reference", group_keys=False)
    .head(25)
    .copy()
)
top_all.to_csv(OUT_TOP_ALL, sep="\t", index=False, encoding="utf-8-sig")

print("================================================")
print("SIX-REFERENCE TARGETED MALODOR HOMOLOGY ANALYSIS")
print("================================================")
print("Raw BLAST hits:", len(blast))
print("Fixed panel genomes:", len(panel))
print()

for ref in refs:
    x = matrix[matrix["reference"] == ref]
    print(f"{ref}: all panel genomes = {len(x)}")
    print(x["similarity_band"].value_counts(dropna=False).to_string())
    print()

print("Relevant-genus species / body-site summary:")
show_cols = [
    "reference","species","source_group","genomes","genomes_with_detectable_hit",
    "high_similarity_n","moderate_similarity_n","distant_full_length_n",
    "weak_partial_n","no_detectable_hit_n","median_best_identity_pct",
    "min_best_identity_pct","max_best_identity_pct","median_query_coverage_pct"
]
print(summary[show_cols].to_string(index=False))

print("\nTop 10 ALL-GENOME hits per reference:")
print(
    top_all.groupby("reference", group_keys=False).head(10)[
        ["reference","species","source_group","genome_id","sseqid","pident",
         "query_coverage_pct","subject_coverage_pct","evalue","bitscore",
         "similarity_band"]
    ].to_string(index=False)
)

print("\nSaved:")
for f in [
    OUT_ALL, OUT_BEST, OUT_RELEVANT, OUT_SUMMARY, OUT_SUMMARY_ALL,
    OUT_TOP, OUT_TOP_ALL
]:
    print(" ", f)

print(
    "\nIMPORTANT:"
    "\n- Similarity bands are descriptive sequence-similarity categories, not validated functional cutoffs."
    "\n- Interpret identity and coverage continuously."
    "\n- The all-genome best-hit table should be used for the cross-genus malodor figures."
    "\n- S. hominis and Cory species lacking adequate controls remain descriptive mechanistic follow-up only."
)

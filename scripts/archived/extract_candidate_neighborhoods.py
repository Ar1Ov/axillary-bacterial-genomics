import re
from pathlib import Path

import numpy as np
import pandas as pd

# ============================================================
# EXTRACT +/- 5-CDS GENOMIC NEIGHBORHOODS AROUND CANDIDATE PGFAMS
# ============================================================
# Inputs:
#   candidate_neighborhood_features_with_coords.tsv
#   candidate_neighborhood_carriers.tsv
#
# Outputs:
#   candidate_neighborhoods_long.tsv
#   candidate_neighborhood_summary.tsv
#   candidate_target_colocalization.tsv
#   candidate_mobile_context_summary.tsv
#
# Each target occurrence gets the 5 CDS immediately upstream and
# 5 CDS immediately downstream ON THE SAME CONTIG.
# ============================================================

FEATURE_FILE = "candidate_neighborhood_features_with_coords.tsv"
CARRIER_FILE = "candidate_neighborhood_carriers.tsv"

OUT_LONG = "candidate_neighborhoods_long.tsv"
OUT_SUMMARY = "candidate_neighborhood_summary.tsv"
OUT_COLOCAL = "candidate_target_colocalization.tsv"
OUT_MOBILE = "candidate_mobile_context_summary.tsv"

WINDOW = 5

for f in [FEATURE_FILE, CARRIER_FILE]:
    if not Path(f).exists():
        raise FileNotFoundError(f)

features = pd.read_csv(
    FEATURE_FILE,
    sep="\t",
    dtype=str,
    keep_default_na=False,
    low_memory=False,
)

carriers = pd.read_csv(
    CARRIER_FILE,
    sep="\t",
    dtype=str,
    keep_default_na=False,
    low_memory=False,
)

required_f = {
    "genome_id",
    "feature.patric_id",
    "feature.pgfam_id",
    "feature.product",
    "feature.sequence_id",
    "feature.start",
    "feature.end",
    "feature.strand",
}
missing = required_f - set(features.columns)
if missing:
    raise ValueError(f"Feature file missing columns: {sorted(missing)}")

required_c = {
    "analysis_arm",
    "candidate_system",
    "feature.pgfam_id",
    "genome_id",
    "group",
}
missing = required_c - set(carriers.columns)
if missing:
    raise ValueError(f"Carrier file missing columns: {sorted(missing)}")

if features["feature.patric_id"].duplicated().any():
    raise ValueError("Duplicate feature IDs found.")

# Numeric coordinates
features["start_num"] = pd.to_numeric(features["feature.start"], errors="coerce")
features["end_num"] = pd.to_numeric(features["feature.end"], errors="coerce")

if features["start_num"].isna().any() or features["end_num"].isna().any():
    raise ValueError("Non-numeric start/end coordinates found.")

# Sort genes in genomic order within each contig.
features = features.sort_values(
    ["genome_id", "feature.sequence_id", "start_num", "end_num", "feature.patric_id"]
).reset_index(drop=True)

# Give each CDS an index within its contig.
features["contig_gene_index"] = (
    features.groupby(["genome_id", "feature.sequence_id"]).cumcount()
)

# Lookup structure by genome + contig
contig_groups = {
    key: sub.reset_index(drop=True)
    for key, sub in features.groupby(["genome_id", "feature.sequence_id"], sort=False)
}

# Target-system lookup for recognizing another shortlisted PGFam in the same neighborhood.
target_map = (
    carriers[["feature.pgfam_id", "candidate_system"]]
    .drop_duplicates()
    .groupby("feature.pgfam_id")["candidate_system"]
    .agg(lambda s: " || ".join(sorted(set(s))))
    .to_dict()
)

MOBILE_RE = re.compile(
    r"phage|prophage|capsid|portal|terminase|integrase|transposase|"
    r"insertion sequence|mobile element|mobilization|pathogenicity island|"
    r"\bSaPI\b|plasmid partition|recombinase",
    re.I,
)

TRANSPORT_RE = re.compile(
    r"transporter|transport system|permease|symporter|antiporter|efflux|"
    r"channel|uptake|ABC|ECF|MFS|oligopeptide|dipeptide",
    re.I,
)

METABOLIC_RE = re.compile(
    r"peptidase|dipeptidase|amidohydrolase|hydrolase|lyase|lipase|esterase|"
    r"thioesterase|oxidase|reductase|dehydrogenase|transferase|synthetase|"
    r"synthase|ligase|deaminase|decarboxylase",
    re.I,
)

rows = []
target_occurrence_counter = 0

# Remove duplicate carrier rows that refer to the same arm/system/pgfam/genome.
carriers = carriers.drop_duplicates(
    ["analysis_arm", "candidate_system", "feature.pgfam_id", "genome_id"]
)

for _, car in carriers.iterrows():
    genome_id = car["genome_id"]
    target_pgfam = car["feature.pgfam_id"]
    candidate_system = car["candidate_system"]

    genome_features = features[
        (features["genome_id"] == genome_id)
        & (features["feature.pgfam_id"] == target_pgfam)
    ].copy()

    if genome_features.empty:
        # This should never happen after previous QC.
        rows.append({
            "analysis_arm": car["analysis_arm"],
            "candidate_system": candidate_system,
            "target_pgfam_id": target_pgfam,
            "genome_id": genome_id,
            "group": car["group"],
            "target_occurrence_id": "MISSING",
            "neighbor_offset": np.nan,
            "neighbor_patric_id": "",
            "neighbor_pgfam_id": "",
            "neighbor_gene": "",
            "neighbor_product": "TARGET NOT FOUND",
            "neighbor_sequence_id": "",
            "neighbor_start": "",
            "neighbor_end": "",
            "neighbor_strand": "",
            "is_target": False,
            "is_other_shortlisted_target": False,
            "other_shortlisted_system": "",
            "mobile_phage_flag": False,
            "transport_flag": False,
            "metabolic_enzyme_flag": False,
        })
        continue

    # One PGFam can have multiple copies in one genome.
    for copy_num, (_, target) in enumerate(genome_features.iterrows(), start=1):
        target_occurrence_counter += 1
        occurrence_id = (
            f"{genome_id}|{target_pgfam}|copy{copy_num}|"
            f"{target['feature.patric_id']}"
        )

        contig_key = (genome_id, target["feature.sequence_id"])
        contig = contig_groups[contig_key]

        matches = contig.index[
            contig["feature.patric_id"] == target["feature.patric_id"]
        ].tolist()

        if len(matches) != 1:
            raise RuntimeError(
                f"Could not uniquely locate {target['feature.patric_id']} in contig."
            )

        center = matches[0]
        left = max(0, center - WINDOW)
        right = min(len(contig) - 1, center + WINDOW)

        neighborhood = contig.iloc[left:right + 1].copy()

        for local_i, (_, gene_row) in enumerate(neighborhood.iterrows()):
            absolute_i = left + local_i
            offset = absolute_i - center

            product = gene_row["feature.product"]
            pgfam = gene_row["feature.pgfam_id"]

            is_target = gene_row["feature.patric_id"] == target["feature.patric_id"]
            is_other_target = (
                (not is_target)
                and bool(pgfam)
                and pgfam in target_map
            )

            base = {
                "analysis_arm": car["analysis_arm"],
                "candidate_system": candidate_system,
                "target_pgfam_id": target_pgfam,
                "genome_id": genome_id,
                "group": car["group"],
                "target_occurrence_id": occurrence_id,
                "target_copy_number": copy_num,
                "target_sequence_id": target["feature.sequence_id"],
                "target_patric_id": target["feature.patric_id"],
                "target_start": int(target["start_num"]),
                "target_end": int(target["end_num"]),
                "target_strand": target["feature.strand"],
                "neighbor_offset": offset,
                "neighbor_patric_id": gene_row["feature.patric_id"],
                "neighbor_pgfam_id": pgfam,
                "neighbor_gene": gene_row["feature.gene"] if "feature.gene" in gene_row else "",
                "neighbor_product": product,
                "neighbor_sequence_id": gene_row["feature.sequence_id"],
                "neighbor_start": int(gene_row["start_num"]),
                "neighbor_end": int(gene_row["end_num"]),
                "neighbor_strand": gene_row["feature.strand"],
                "is_target": is_target,
                "is_other_shortlisted_target": is_other_target,
                "other_shortlisted_system": target_map.get(pgfam, "") if is_other_target else "",
                "mobile_phage_flag": bool(MOBILE_RE.search(product)),
                "transport_flag": bool(TRANSPORT_RE.search(product)),
                "metabolic_enzyme_flag": bool(METABOLIC_RE.search(product)),
            }

            # Preserve useful cohort metadata when available.
            for col in ["species", "bioproject", "site"]:
                if col in carriers.columns:
                    base[col] = car[col]

            rows.append(base)

long = pd.DataFrame(rows)

if (long["target_occurrence_id"] == "MISSING").any():
    missing_n = int((long["target_occurrence_id"] == "MISSING").sum())
    raise RuntimeError(f"{missing_n} target carrier rows were not found in coordinate features.")

long = long.sort_values(
    ["candidate_system", "analysis_arm", "genome_id", "target_occurrence_id", "neighbor_offset"]
)

long.to_csv(
    OUT_LONG,
    sep="\t",
    index=False,
    encoding="utf-8-sig",
)

# ============================================================
# OCCURRENCE-LEVEL SUMMARY
# ============================================================

summary_rows = []

for occ, g in long.groupby("target_occurrence_id", sort=False):
    target = g[g["is_target"]].iloc[0]

    non_target = g[~g["is_target"]]

    other_targets = sorted(
        set(
            non_target.loc[
                non_target["is_other_shortlisted_target"],
                "other_shortlisted_system",
            ].replace("", np.nan).dropna()
        )
    )

    mobile_neighbors = int(non_target["mobile_phage_flag"].sum())
    transport_neighbors = int(non_target["transport_flag"].sum())
    metabolic_neighbors = int(non_target["metabolic_enzyme_flag"].sum())

    summary_rows.append({
        "target_occurrence_id": occ,
        "analysis_arm": target["analysis_arm"],
        "candidate_system": target["candidate_system"],
        "target_pgfam_id": target["target_pgfam_id"],
        "genome_id": target["genome_id"],
        "group": target["group"],
        "species": target.get("species", ""),
        "bioproject": target.get("bioproject", ""),
        "site": target.get("site", ""),
        "target_patric_id": target["target_patric_id"],
        "target_sequence_id": target["target_sequence_id"],
        "target_start": target["target_start"],
        "target_end": target["target_end"],
        "target_strand": target["target_strand"],
        "genes_in_window": len(g),
        "other_shortlisted_targets_in_window": len(other_targets),
        "other_shortlisted_systems": " || ".join(other_targets),
        "mobile_phage_neighbors": mobile_neighbors,
        "transport_neighbors": transport_neighbors,
        "metabolic_enzyme_neighbors": metabolic_neighbors,
        "mobile_context": mobile_neighbors > 0,
    })

summary = pd.DataFrame(summary_rows)

summary.to_csv(
    OUT_SUMMARY,
    sep="\t",
    index=False,
    encoding="utf-8-sig",
)

# ============================================================
# TARGET CO-LOCALIZATION SUMMARY
# ============================================================

coloc = (
    summary[
        summary["other_shortlisted_targets_in_window"] > 0
    ][
        [
            "candidate_system",
            "target_pgfam_id",
            "analysis_arm",
            "group",
            "genome_id",
            "target_patric_id",
            "other_shortlisted_systems",
            "mobile_context",
        ]
    ]
    .copy()
)

coloc.to_csv(
    OUT_COLOCAL,
    sep="\t",
    index=False,
    encoding="utf-8-sig",
)

# ============================================================
# MOBILE CONTEXT BY CANDIDATE
# ============================================================

mobile_summary = (
    summary.groupby(
        ["candidate_system", "target_pgfam_id", "analysis_arm", "group"],
        dropna=False,
    )
    .agg(
        target_occurrences=("target_occurrence_id", "nunique"),
        genomes=("genome_id", "nunique"),
        occurrences_with_mobile_neighbor=("mobile_context", "sum"),
        mean_mobile_neighbors=("mobile_phage_neighbors", "mean"),
        mean_transport_neighbors=("transport_neighbors", "mean"),
        mean_metabolic_neighbors=("metabolic_enzyme_neighbors", "mean"),
        occurrences_with_other_shortlisted_target=(
            "other_shortlisted_targets_in_window",
            lambda x: int((x > 0).sum()),
        ),
    )
    .reset_index()
)

mobile_summary["fraction_mobile_context"] = (
    mobile_summary["occurrences_with_mobile_neighbor"]
    / mobile_summary["target_occurrences"]
)

mobile_summary["fraction_colocalized_with_shortlisted_target"] = (
    mobile_summary["occurrences_with_other_shortlisted_target"]
    / mobile_summary["target_occurrences"]
)

mobile_summary.to_csv(
    OUT_MOBILE,
    sep="\t",
    index=False,
    encoding="utf-8-sig",
)

# ============================================================
# CONSOLE SUMMARY
# ============================================================

print("================================================")
print("CANDIDATE GENOMIC NEIGHBORHOOD EXTRACTION")
print("================================================")
print("Neighborhood window: +/-", WINDOW, "CDS")
print("Target carrier rows:", len(carriers))
print("Target feature occurrences:", summary["target_occurrence_id"].nunique())
print("Neighborhood rows:", len(long))
print("Unique genomes:", summary["genome_id"].nunique())
print()

print("Candidates with the most co-localization to another shortlisted target:")
tmp = (
    mobile_summary.groupby(["candidate_system", "target_pgfam_id"])
    .agg(
        occurrences=("target_occurrences", "sum"),
        coloc_occurrences=("occurrences_with_other_shortlisted_target", "sum"),
        mobile_occurrences=("occurrences_with_mobile_neighbor", "sum"),
    )
    .reset_index()
)
tmp["coloc_fraction"] = tmp["coloc_occurrences"] / tmp["occurrences"]
tmp["mobile_fraction"] = tmp["mobile_occurrences"] / tmp["occurrences"]
tmp = tmp.sort_values(
    ["coloc_fraction", "occurrences"],
    ascending=[False, False],
)

print(
    tmp[
        [
            "candidate_system",
            "target_pgfam_id",
            "occurrences",
            "coloc_occurrences",
            "coloc_fraction",
            "mobile_occurrences",
            "mobile_fraction",
        ]
    ].head(25).to_string(index=False)
)

print()
print("Saved:")
print(" ", OUT_LONG)
print(" ", OUT_SUMMARY)
print(" ", OUT_COLOCAL)
print(" ", OUT_MOBILE)

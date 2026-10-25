import pandas as pd
import numpy as np
from pathlib import Path

# ============================================================
# FINAL CANDIDATE NEIGHBORHOOD AUDIT
# ============================================================
# This script does NOT create a new significance test.
# It summarizes genomic context for the shortlisted systems so we can
# decide which loci deserve final follow-up.
# ============================================================

LONG_FILE = "candidate_neighborhoods_long.tsv"
SUMMARY_FILE = "candidate_neighborhood_summary.tsv"
CARRIER_SUMMARY_FILE = "candidate_neighborhood_carrier_summary.tsv"

OUT_AUDIT = "final_candidate_neighborhood_audit.tsv"
OUT_NEIGHBORS = "final_candidate_common_neighbors.tsv"
OUT_COHERENCE = "final_candidate_system_coherence.tsv"

for f in [LONG_FILE, SUMMARY_FILE, CARRIER_SUMMARY_FILE]:
    if not Path(f).exists():
        raise FileNotFoundError(f)

long = pd.read_csv(LONG_FILE, sep="\t", dtype=str, keep_default_na=False, low_memory=False)
summary = pd.read_csv(SUMMARY_FILE, sep="\t", dtype=str, keep_default_na=False, low_memory=False)
carrier = pd.read_csv(CARRIER_SUMMARY_FILE, sep="\t", dtype=str, keep_default_na=False, low_memory=False)

# Candidate systems to audit.
FOCUS = {
    "ABC antimicrobial-peptide transport ATPase",
    "Arginine/ornithine antiporter ArcD",
    "Carnitine-utilization associated thioesterase",
    "ECF transporter - ATPase component",
    "ECF transporter - transmembrane component",
    "Esterase/lipase",
    "Glycine betaine/L-proline transporter ProW",
    "N-acyl-L-amino-acid amidohydrolase",
    "Opp oligopeptide transport - OppA",
    "Opp oligopeptide transport - OppB",
    "Opp oligopeptide transport - OppC",
    "Opp oligopeptide transport - OppD",
    "Opp oligopeptide transport - OppF",
    "Uncharacterized MFS-type transporter",
}

# Known controls retained only as contextual baselines.
CONTROLS = {
    "Known control - DtpT",
    "Known control - PepV",
    "Known control - Staphylococcus MetC-like lyase",
    "Known control - Corynebacterium AecD/MetC-like lyase",
}

USE = FOCUS | CONTROLS

long = long[long["candidate_system"].isin(USE)].copy()
summary = summary[summary["candidate_system"].isin(USE)].copy()
carrier = carrier[carrier["candidate_system"].isin(USE)].copy()

# Numeric conversions
for c in [
    "mobile_phage_neighbors",
    "transport_neighbors",
    "metabolic_enzyme_neighbors",
    "other_shortlisted_targets_in_window",
]:
    if c in summary.columns:
        summary[c] = pd.to_numeric(summary[c], errors="coerce").fillna(0)

summary["mobile_context_bool"] = (
    summary["mobile_context"].astype(str).str.lower().eq("true")
)

# ============================================================
# 1. Candidate audit by analysis arm + group
# ============================================================

audit_rows = []

for keys, g in summary.groupby(
    ["candidate_system", "target_pgfam_id", "analysis_arm", "group"],
    dropna=False
):
    candidate_system, pgfam, arm, group = keys

    carrier_match = carrier[
        (carrier["candidate_system"] == candidate_system)
        & (carrier["feature.pgfam_id"] == pgfam)
        & (carrier["analysis_arm"] == arm)
        & (carrier["group"] == group)
    ]

    carrier_genomes = (
        pd.to_numeric(carrier_match["carrier_genomes"], errors="coerce").sum()
        if len(carrier_match)
        else g["genome_id"].nunique()
    )

    audit_rows.append({
        "candidate_system": candidate_system,
        "pgfam_id": pgfam,
        "role": "CONTROL" if candidate_system in CONTROLS else "CANDIDATE",
        "analysis_arm": arm,
        "group": group,
        "carrier_genomes": int(carrier_genomes),
        "target_occurrences": g["target_occurrence_id"].nunique(),
        "mobile_context_occurrences": int(g["mobile_context_bool"].sum()),
        "mobile_context_fraction": float(g["mobile_context_bool"].mean()),
        "colocalized_shortlist_occurrences": int(
            (g["other_shortlisted_targets_in_window"] > 0).sum()
        ),
        "colocalized_shortlist_fraction": float(
            (g["other_shortlisted_targets_in_window"] > 0).mean()
        ),
        "mean_transport_neighbors": float(g["transport_neighbors"].mean()),
        "mean_metabolic_neighbors": float(g["metabolic_enzyme_neighbors"].mean()),
    })

audit = pd.DataFrame(audit_rows).sort_values(
    ["role", "candidate_system", "analysis_arm", "group"]
)

audit.to_csv(OUT_AUDIT, sep="\t", index=False, encoding="utf-8-sig")

# ============================================================
# 2. Most recurrent neighboring genes around each candidate
# ============================================================

neighbor_rows = []

# Exclude target row itself.
neighbors = long[
    (long["is_target"].astype(str).str.lower() != "true")
].copy()

# Deduplicate same PGFam within one target occurrence.
neighbors = neighbors[
    neighbors["neighbor_pgfam_id"].ne("")
].drop_duplicates(
    [
        "candidate_system",
        "target_pgfam_id",
        "analysis_arm",
        "group",
        "target_occurrence_id",
        "neighbor_pgfam_id",
    ]
)

for keys, g in neighbors.groupby(
    ["candidate_system", "target_pgfam_id", "analysis_arm", "group"],
    dropna=False
):
    candidate_system, target_pgfam, arm, group = keys

    denom = summary[
        (summary["candidate_system"] == candidate_system)
        & (summary["target_pgfam_id"] == target_pgfam)
        & (summary["analysis_arm"] == arm)
        & (summary["group"] == group)
    ]["target_occurrence_id"].nunique()

    if denom == 0:
        continue

    counts = (
        g.groupby(["neighbor_pgfam_id", "neighbor_product"])
        ["target_occurrence_id"]
        .nunique()
        .rename("occurrences_with_neighbor")
        .reset_index()
    )
    counts["neighbor_occurrence_fraction"] = (
        counts["occurrences_with_neighbor"] / denom
    )

    counts = counts.sort_values(
        ["neighbor_occurrence_fraction", "occurrences_with_neighbor"],
        ascending=[False, False]
    ).head(12)

    for _, r in counts.iterrows():
        neighbor_rows.append({
            "candidate_system": candidate_system,
            "target_pgfam_id": target_pgfam,
            "analysis_arm": arm,
            "group": group,
            "target_occurrences": denom,
            "neighbor_pgfam_id": r["neighbor_pgfam_id"],
            "neighbor_product": r["neighbor_product"],
            "occurrences_with_neighbor": int(r["occurrences_with_neighbor"]),
            "neighbor_occurrence_fraction": float(r["neighbor_occurrence_fraction"]),
        })

common = pd.DataFrame(neighbor_rows)
common.to_csv(OUT_NEIGHBORS, sep="\t", index=False, encoding="utf-8-sig")

# ============================================================
# 3. Coherence of multi-component systems
# ============================================================

SYSTEMS = {
    "OPP_LOCUS": {
        "PGF_02147779",  # OppA
        "PGF_00977978",  # OppB
        "PGF_10569020",  # OppC
        "PGF_09990804",  # OppD
        "PGF_07131324",  # OppF
    },
    "ECF_LOCUS": {
        "PGF_03285992",  # ATPase
        "PGF_00267745",  # transmembrane
    },
}

coherence_rows = []

# long already contains each target occurrence +/-5 CDS.
for system_name, components in SYSTEMS.items():
    sub = long[long["target_pgfam_id"].isin(components)].copy()

    for occ, g in sub.groupby("target_occurrence_id", sort=False):
        target = g[g["is_target"].astype(str).str.lower().eq("true")].iloc[0]

        present_components = set(
            g.loc[
                g["neighbor_pgfam_id"].isin(components),
                "neighbor_pgfam_id"
            ]
        )

        coherence_rows.append({
            "system": system_name,
            "analysis_arm": target["analysis_arm"],
            "group": target["group"],
            "genome_id": target["genome_id"],
            "target_occurrence_id": occ,
            "target_component": target["target_pgfam_id"],
            "components_detected_within_window": len(present_components),
            "component_set_within_window": " || ".join(sorted(present_components)),
            "all_components_within_window": components.issubset(present_components),
        })

coherence = pd.DataFrame(coherence_rows)

coherence_summary = (
    coherence.groupby(
        ["system", "analysis_arm", "group", "target_component"],
        dropna=False
    )
    .agg(
        target_occurrences=("target_occurrence_id", "nunique"),
        mean_components_within_window=("components_detected_within_window", "mean"),
        all_components_within_window_occurrences=("all_components_within_window", "sum"),
    )
    .reset_index()
)

coherence_summary["all_components_within_window_fraction"] = (
    coherence_summary["all_components_within_window_occurrences"]
    / coherence_summary["target_occurrences"]
)

coherence_summary.to_csv(
    OUT_COHERENCE, sep="\t", index=False, encoding="utf-8-sig"
)

# ============================================================
# Console summary
# ============================================================

print("================================================")
print("FINAL CANDIDATE NEIGHBORHOOD AUDIT")
print("================================================")

print("\nCandidate context summary:")
candidate_audit = audit[audit["role"] == "CANDIDATE"].copy()
print(
    candidate_audit[
        [
            "candidate_system",
            "analysis_arm",
            "group",
            "carrier_genomes",
            "mobile_context_fraction",
            "colocalized_shortlist_fraction",
            "mean_transport_neighbors",
            "mean_metabolic_neighbors",
        ]
    ].to_string(index=False)
)

print("\nMulti-component locus coherence:")
print(coherence_summary.to_string(index=False))

print("\nMost recurrent neighbors for ECF + ABC antimicrobial-peptide transport:")
focus_names = {
    "ECF transporter - ATPase component",
    "ECF transporter - transmembrane component",
    "ABC antimicrobial-peptide transport ATPase",
}
focus_neighbors = common[common["candidate_system"].isin(focus_names)].copy()
print(
    focus_neighbors[
        [
            "candidate_system",
            "analysis_arm",
            "group",
            "neighbor_pgfam_id",
            "neighbor_product",
            "neighbor_occurrence_fraction",
        ]
    ].head(40).to_string(index=False)
)

print("\nSaved:")
print(" ", OUT_AUDIT)
print(" ", OUT_NEIGHBORS)
print(" ", OUT_COHERENCE)

import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# ============================================================
# PRIMARY S. EPIDERMIDIS STRICT AXILLA ANALYSIS
# ============================================================
# Cohort:
#   5 STRICT_AXILLA genomes
#   117 explicit OTHER_SKIN genomes
#
# Statistical design:
#   - PGFam presence/absence per genome
#   - unbiased prefilter: family present in >=3 genomes AND absent in >=3 genomes
#   - two-sided Fisher exact test
#   - Benjamini-Hochberg FDR across all tested PGFams
#   - axilla enrichment is assigned only AFTER testing, based on prevalence
#   - independent axilla BioProject recurrence is reported as a replication metric
#
# BioProject is NOT "adjusted away" because strict AXILLA and OTHER_SKIN genomes
# do not share BioProjects. Instead, project recurrence is used as an explicit
# robustness criterion, and a same-study PRJNA803478 sensitivity analysis should
# be performed separately.
# ============================================================

COHORT_FILE = "staph_epi_strict_primary_cohort.tsv"
PRESENCE_FILE = "staph_epi_strict_primary_pgfam_presence_long.tsv"
FEATURE_FILE = "staph_epi_strict_primary_features.tsv"

OUT_ALL = "staph_epi_strict_primary_pgfam_results.tsv"
OUT_FDR = "staph_epi_strict_primary_fdr_axilla_candidates.tsv"
OUT_REPLICATED = "staph_epi_strict_primary_replicated_candidates.tsv"

MIN_PRESENT = 3
MIN_ABSENT = 3
FDR_CUTOFF = 0.05
MIN_AXILLA_PROJECTS = 3

# ============================================================
# LOAD + VALIDATE COHORT
# ============================================================

cohort = pd.read_csv(
    COHORT_FILE, sep="\t", dtype=str, keep_default_na=False
)

required_cohort = {"genome_id", "bioproject", "site", "group"}
missing = required_cohort - set(cohort.columns)
if missing:
    raise ValueError(f"Missing cohort columns: {sorted(missing)}")

if cohort["genome_id"].duplicated().any():
    raise ValueError("Duplicate genome_id values found in cohort.")

expected_groups = {"AXILLA", "OTHER_SKIN"}
if set(cohort["group"]) != expected_groups:
    raise ValueError(
        f"Unexpected groups: {sorted(set(cohort['group']))}"
    )

group_counts = cohort["group"].value_counts().to_dict()

if group_counts.get("AXILLA", 0) != 5:
    raise ValueError(
        f"Expected 5 AXILLA genomes, found {group_counts.get('AXILLA', 0)}"
    )

if group_counts.get("OTHER_SKIN", 0) != 117:
    raise ValueError(
        f"Expected 117 OTHER_SKIN genomes, found "
        f"{group_counts.get('OTHER_SKIN', 0)}"
    )

ax_projects = sorted(
    cohort.loc[cohort["group"] == "AXILLA", "bioproject"]
    .replace("", np.nan)
    .dropna()
    .unique()
)

control_projects = sorted(
    cohort.loc[cohort["group"] == "OTHER_SKIN", "bioproject"]
    .replace("", np.nan)
    .dropna()
    .unique()
)

if len(ax_projects) != 4:
    raise ValueError(
        f"Expected 4 AXILLA BioProjects, found {len(ax_projects)}"
    )

if len(control_projects) != 15:
    raise ValueError(
        f"Expected 15 OTHER_SKIN BioProjects, found {len(control_projects)}"
    )

shared_projects = sorted(set(ax_projects) & set(control_projects))
if shared_projects:
    raise ValueError(
        "Strict cohort unexpectedly contains BioProjects shared between "
        f"AXILLA and OTHER_SKIN: {shared_projects}"
    )

# ============================================================
# LOAD + VALIDATE PGFAM PRESENCE
# ============================================================

presence = pd.read_csv(
    PRESENCE_FILE, sep="\t", dtype=str, keep_default_na=False
)

required_presence = {
    "genome_id", "feature.pgfam_id", "bioproject", "site", "group"
}
missing = required_presence - set(presence.columns)
if missing:
    raise ValueError(f"Missing presence columns: {sorted(missing)}")

presence = presence[presence["feature.pgfam_id"].ne("")].copy()

# Presence table must be one row per genome-PGFam.
if presence.duplicated(["genome_id", "feature.pgfam_id"]).any():
    raise ValueError(
        "Duplicate genome-PGFam pairs found. Presence table is not deduplicated."
    )

cohort_ids = set(cohort["genome_id"])
presence_ids = set(presence["genome_id"])

extra_ids = sorted(presence_ids - cohort_ids)
if extra_ids:
    raise ValueError(
        f"Presence file contains genomes outside cohort. Examples: {extra_ids[:10]}"
    )

missing_ids = sorted(cohort_ids - presence_ids)
if missing_ids:
    raise ValueError(
        f"Cohort genomes missing entirely from presence file: {missing_ids[:10]}"
    )

# Verify labels in the presence file match the locked cohort.
label_check = presence[
    ["genome_id", "bioproject", "site", "group"]
].drop_duplicates()

merged_check = label_check.merge(
    cohort,
    on="genome_id",
    how="left",
    suffixes=("_presence", "_cohort"),
    validate="one_to_one",
)

for col in ["bioproject", "site", "group"]:
    bad = merged_check[
        merged_check[f"{col}_presence"] != merged_check[f"{col}_cohort"]
    ]
    if not bad.empty:
        raise ValueError(
            f"Presence/cohort mismatch in {col}. "
            f"Example genome: {bad.iloc[0]['genome_id']}"
        )

# ============================================================
# PRODUCT ANNOTATION SUMMARY
# ============================================================

features = pd.read_csv(
    FEATURE_FILE,
    sep="\t",
    dtype=str,
    keep_default_na=False,
    low_memory=False,
)

required_features = {"feature.pgfam_id", "feature.product"}
missing = required_features - set(features.columns)
if missing:
    raise ValueError(f"Missing feature columns: {sorted(missing)}")

annot = features[
    features["feature.pgfam_id"].ne("")
][["feature.pgfam_id", "feature.product"]].copy()

def summarize_products(s):
    vals = [x.strip() for x in s if str(x).strip()]
    if not vals:
        return pd.Series(
            {
                "representative_product": "",
                "n_distinct_products": 0,
            }
        )

    vc = pd.Series(vals).value_counts()
    return pd.Series(
        {
            "representative_product": vc.index[0],
            "n_distinct_products": int(vc.size),
        }
    )

product_summary = (
    annot.groupby("feature.pgfam_id")["feature.product"]
    .apply(summarize_products)
    .unstack()
    .reset_index()
)

product_summary["n_distinct_products"] = (
    pd.to_numeric(product_summary["n_distinct_products"], errors="coerce")
    .fillna(0)
    .astype(int)
)

# ============================================================
# BUILD COUNTS
# ============================================================

ax_ids = set(
    cohort.loc[cohort["group"] == "AXILLA", "genome_id"]
)
control_ids = set(
    cohort.loc[cohort["group"] == "OTHER_SKIN", "genome_id"]
)

n_ax = len(ax_ids)
n_control = len(control_ids)
n_total = n_ax + n_control

# Sets make project support calculations easy and reproducible.
genome_to_project = cohort.set_index("genome_id")["bioproject"].to_dict()

rows = []

for pgfam, g in presence.groupby("feature.pgfam_id", sort=False):
    carriers = set(g["genome_id"])

    ax_present = len(carriers & ax_ids)
    control_present = len(carriers & control_ids)

    ax_absent = n_ax - ax_present
    control_absent = n_control - control_present

    total_present = ax_present + control_present
    total_absent = n_total - total_present

    ax_prev = ax_present / n_ax
    control_prev = control_present / n_control
    prevalence_diff = ax_prev - control_prev

    ax_carrier_projects = {
        genome_to_project[x]
        for x in (carriers & ax_ids)
        if genome_to_project.get(x, "")
    }

    control_carrier_projects = {
        genome_to_project[x]
        for x in (carriers & control_ids)
        if genome_to_project.get(x, "")
    }

    tested = (
        total_present >= MIN_PRESENT
        and total_absent >= MIN_ABSENT
    )

    if tested:
        table = [
            [ax_present, ax_absent],
            [control_present, control_absent],
        ]
        crude_or, p = fisher_exact(table, alternative="two-sided")
    else:
        crude_or = np.nan
        p = np.nan

    # Haldane-Anscombe corrected OR + Wald CI, useful when a cell is zero.
    aa = ax_present + 0.5
    bb = ax_absent + 0.5
    cc = control_present + 0.5
    dd = control_absent + 0.5

    corrected_or = (aa * dd) / (bb * cc)
    se_log_or = math.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    log_or = math.log(corrected_or)
    ci_low = math.exp(log_or - 1.96 * se_log_or)
    ci_high = math.exp(log_or + 1.96 * se_log_or)

    rows.append(
        {
            "pgfam_id": pgfam,
            "AXILLA_present": ax_present,
            "AXILLA_absent": ax_absent,
            "OTHER_SKIN_present": control_present,
            "OTHER_SKIN_absent": control_absent,
            "AXILLA_prevalence": ax_prev,
            "OTHER_SKIN_prevalence": control_prev,
            "prevalence_difference": prevalence_diff,
            "total_present": total_present,
            "total_absent": total_absent,
            "tested": tested,
            "fisher_odds_ratio": crude_or,
            "haldane_anscombe_or": corrected_or,
            "or_95ci_low": ci_low,
            "or_95ci_high": ci_high,
            "p_value": p,
            "AXILLA_projects_with_gene": len(ax_carrier_projects),
            "AXILLA_projects_total": len(ax_projects),
            "AXILLA_project_fraction": (
                len(ax_carrier_projects) / len(ax_projects)
            ),
            "OTHER_SKIN_projects_with_gene": len(control_carrier_projects),
            "OTHER_SKIN_projects_total": len(control_projects),
            "OTHER_SKIN_project_fraction": (
                len(control_carrier_projects) / len(control_projects)
            ),
        }
    )

results = pd.DataFrame(rows)

# ============================================================
# MULTIPLE-TEST CORRECTION
# ============================================================

results["fdr_bh"] = np.nan

mask = results["tested"] & results["p_value"].notna()

if mask.sum() == 0:
    raise RuntimeError("No PGFams passed the independent prefilter.")

results.loc[mask, "fdr_bh"] = multipletests(
    results.loc[mask, "p_value"].astype(float),
    alpha=FDR_CUTOFF,
    method="fdr_bh",
)[1]

results["direction"] = np.select(
    [
        results["prevalence_difference"] > 0,
        results["prevalence_difference"] < 0,
    ],
    [
        "AXILLA_ENRICHED",
        "OTHER_SKIN_ENRICHED",
    ],
    default="EQUAL",
)

results["fdr_significant"] = (
    results["tested"]
    & results["fdr_bh"].notna()
    & (results["fdr_bh"] < FDR_CUTOFF)
)

results["axilla_project_replicated"] = (
    results["AXILLA_projects_with_gene"] >= MIN_AXILLA_PROJECTS
)

# Merge representative BV-BRC product annotations.
results = results.merge(
    product_summary,
    left_on="pgfam_id",
    right_on="feature.pgfam_id",
    how="left",
)

results = results.drop(columns=["feature.pgfam_id"], errors="ignore")

results["representative_product"] = (
    results["representative_product"].fillna("")
)
results["n_distinct_products"] = (
    results["n_distinct_products"].fillna(0).astype(int)
)

# Sort with tested families first, then significance/effect.
results = results.sort_values(
    [
        "tested",
        "fdr_bh",
        "prevalence_difference",
        "AXILLA_projects_with_gene",
    ],
    ascending=[False, True, False, False],
    na_position="last",
)

results.to_csv(
    OUT_ALL, sep="\t", index=False, encoding="utf-8-sig"
)

fdr_hits = results[
    results["fdr_significant"]
    & results["direction"].eq("AXILLA_ENRICHED")
].copy()

fdr_hits = fdr_hits.sort_values(
    [
        "fdr_bh",
        "AXILLA_projects_with_gene",
        "prevalence_difference",
    ],
    ascending=[True, False, False],
)

fdr_hits.to_csv(
    OUT_FDR, sep="\t", index=False, encoding="utf-8-sig"
)

replicated = fdr_hits[
    fdr_hits["axilla_project_replicated"]
].copy()

replicated.to_csv(
    OUT_REPLICATED, sep="\t", index=False, encoding="utf-8-sig"
)

# ============================================================
# SUMMARY
# ============================================================

print("================================================")
print("S. EPIDERMIDIS STRICT PRIMARY PGFAM ANALYSIS")
print("================================================")
print("Cohort:", n_ax, "AXILLA /", n_control, "OTHER_SKIN")
print("AXILLA BioProjects:", len(ax_projects))
print("OTHER_SKIN BioProjects:", len(control_projects))
print("Shared BioProjects:", len(shared_projects))
print("PGFams observed:", len(results))
print(
    f"PGFams tested (present >= {MIN_PRESENT}, absent >= {MIN_ABSENT}):",
    int(results["tested"].sum()),
)
print(
    "FDR < 0.05, AXILLA-enriched:",
    len(fdr_hits),
)
print(
    f"FDR < 0.05, AXILLA-enriched, present in >= "
    f"{MIN_AXILLA_PROJECTS}/4 AXILLA projects:",
    len(replicated),
)

if len(replicated):
    show = replicated[
        [
            "pgfam_id",
            "representative_product",
            "AXILLA_present",
            "OTHER_SKIN_present",
            "AXILLA_prevalence",
            "OTHER_SKIN_prevalence",
            "prevalence_difference",
            "haldane_anscombe_or",
            "or_95ci_low",
            "or_95ci_high",
            "p_value",
            "fdr_bh",
            "AXILLA_projects_with_gene",
            "OTHER_SKIN_projects_with_gene",
        ]
    ].head(30)

    print("\nTop replicated AXILLA-enriched candidates:")
    print(show.to_string(index=False))
else:
    print(
        "\nNo candidate met the conservative FDR + "
        "3-project replication rule."
    )

print("\nSaved:")
print(" ", OUT_ALL)
print(" ", OUT_FDR)
print(" ", OUT_REPLICATED)

print(
    "\nINTERPRETATION NOTE:"
    "\nThe strict AXILLA and OTHER_SKIN genomes come from different BioProjects, "
    "so these associations cannot by themselves separate body site from study/"
    "lineage effects. The >=3/4 AXILLA-project criterion tests recurrence across "
    "independent axilla projects, and the PRJNA803478 same-study sensitivity "
    "analysis should be used next to assess robustness."
)

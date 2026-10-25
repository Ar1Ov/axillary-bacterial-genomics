import math
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# ============================================================
# EXPLORATORY E. COLI AXILLA VS OTHER-SKIN ANALYSIS
# ============================================================
# Final cohort:
#   3 AXILLA (axillary swab; PRJNA1238219)
#   5 OTHER_SKIN (umbilicus/lower-leg/hand; 3 BioProjects)
#
# IMPORTANT:
#   Case/control status is completely separated by BioProject.
#   n is very small (3 vs 5), so this analysis is descriptive/
#   exploratory and is intended only for cross-species support.
#
# Statistical approach:
#   - PGFam presence/absence per genome
#   - unbiased prefilter: present in >=2 genomes and absent in >=2 genomes
#   - two-sided Fisher exact test
#   - Benjamini-Hochberg FDR across all tested PGFams
#   - effect sizes + prevalence differences retained for synthesis
# ============================================================

COHORT_FILE = "ecoli_final_cohort.tsv"
PRESENCE_FILE = "ecoli_final_pgfam_presence_long.tsv"
FEATURE_FILE = "ecoli_final_features.tsv"

OUT_ALL = "ecoli_exploratory_pgfam_results.tsv"
OUT_FDR = "ecoli_exploratory_fdr_axilla_candidates.tsv"
OUT_NOMINAL = "ecoli_exploratory_nominal_axilla.tsv"
OUT_TOP = "ecoli_exploratory_top_directional.tsv"

MIN_PRESENT = 2
MIN_ABSENT = 2
FDR_CUTOFF = 0.05
NOMINAL_P = 0.05

# ----------------------------
# Load + validate cohort
# ----------------------------

cohort = pd.read_csv(COHORT_FILE, sep="\t", dtype=str, keep_default_na=False)

required = {"genome_id", "bioproject", "site", "group"}
missing = required - set(cohort.columns)
if missing:
    raise ValueError(f"Missing cohort columns: {sorted(missing)}")

if cohort["genome_id"].duplicated().any():
    raise ValueError("Duplicate genome IDs found in cohort.")

counts = cohort["group"].value_counts().to_dict()

if counts.get("AXILLA", 0) != 3:
    raise ValueError(f"Expected 3 AXILLA genomes, found {counts.get('AXILLA', 0)}")

if counts.get("OTHER_SKIN", 0) != 5:
    raise ValueError(f"Expected 5 OTHER_SKIN genomes, found {counts.get('OTHER_SKIN', 0)}")

ax_projects = sorted(
    cohort.loc[cohort["group"] == "AXILLA", "bioproject"].unique()
)
skin_projects = sorted(
    cohort.loc[cohort["group"] == "OTHER_SKIN", "bioproject"].unique()
)

shared = sorted(set(ax_projects) & set(skin_projects))

if shared:
    raise ValueError(f"Unexpected shared BioProjects between groups: {shared}")

# ----------------------------
# Load + validate presence
# ----------------------------

presence = pd.read_csv(
    PRESENCE_FILE, sep="\t", dtype=str, keep_default_na=False
)

required_p = {"genome_id", "feature.pgfam_id", "bioproject", "site", "group"}
missing = required_p - set(presence.columns)
if missing:
    raise ValueError(f"Missing presence columns: {sorted(missing)}")

presence = presence[presence["feature.pgfam_id"].ne("")].copy()

if presence.duplicated(["genome_id", "feature.pgfam_id"]).any():
    raise ValueError("Duplicate genome-PGFam pairs found.")

cohort_ids = set(cohort["genome_id"])
presence_ids = set(presence["genome_id"])

if presence_ids != cohort_ids:
    raise ValueError(
        f"Presence/cohort mismatch. Missing: {sorted(cohort_ids-presence_ids)[:10]}; "
        f"extra: {sorted(presence_ids-cohort_ids)[:10]}"
    )

# ----------------------------
# Product annotation summary
# ----------------------------

features = pd.read_csv(
    FEATURE_FILE,
    sep="\t",
    dtype=str,
    keep_default_na=False,
    low_memory=False,
)

if "feature.pgfam_id" not in features or "feature.product" not in features:
    raise ValueError("Feature file missing PGFam or product columns.")

annot = features.loc[
    features["feature.pgfam_id"].ne(""),
    ["feature.pgfam_id", "feature.product"],
].copy()

def top_product(s):
    vals = pd.Series([x.strip() for x in s if str(x).strip()])
    if vals.empty:
        return ""
    return vals.value_counts().index[0]

product_summary = (
    annot.groupby("feature.pgfam_id")["feature.product"]
    .agg(representative_product=top_product, n_distinct_products="nunique")
    .reset_index()
)

# ----------------------------
# Fisher tests
# ----------------------------

ax_ids = set(cohort.loc[cohort["group"] == "AXILLA", "genome_id"])
skin_ids = set(cohort.loc[cohort["group"] == "OTHER_SKIN", "genome_id"])

n_ax = len(ax_ids)
n_skin = len(skin_ids)
n_total = n_ax + n_skin

rows = []

for pgfam, g in presence.groupby("feature.pgfam_id", sort=False):
    carriers = set(g["genome_id"])

    ax_present = len(carriers & ax_ids)
    skin_present = len(carriers & skin_ids)

    ax_absent = n_ax - ax_present
    skin_absent = n_skin - skin_present

    total_present = ax_present + skin_present
    total_absent = n_total - total_present

    ax_prev = ax_present / n_ax
    skin_prev = skin_present / n_skin
    diff = ax_prev - skin_prev

    tested = (
        total_present >= MIN_PRESENT
        and total_absent >= MIN_ABSENT
    )

    if tested:
        crude_or, p = fisher_exact(
            [[ax_present, ax_absent],
             [skin_present, skin_absent]],
            alternative="two-sided",
        )
    else:
        crude_or, p = np.nan, np.nan

    # Haldane-Anscombe corrected OR + approximate 95% CI
    aa = ax_present + 0.5
    bb = ax_absent + 0.5
    cc = skin_present + 0.5
    dd = skin_absent + 0.5

    corrected_or = (aa * dd) / (bb * cc)
    se = math.sqrt(1/aa + 1/bb + 1/cc + 1/dd)
    log_or = math.log(corrected_or)

    ci_low = math.exp(log_or - 1.96 * se)
    ci_high = math.exp(log_or + 1.96 * se)

    rows.append({
        "pgfam_id": pgfam,
        "AXILLA_present": ax_present,
        "AXILLA_absent": ax_absent,
        "OTHER_SKIN_present": skin_present,
        "OTHER_SKIN_absent": skin_absent,
        "AXILLA_prevalence": ax_prev,
        "OTHER_SKIN_prevalence": skin_prev,
        "prevalence_difference": diff,
        "total_present": total_present,
        "total_absent": total_absent,
        "tested": tested,
        "fisher_odds_ratio": crude_or,
        "haldane_anscombe_or": corrected_or,
        "or_95ci_low": ci_low,
        "or_95ci_high": ci_high,
        "p_value": p,
    })

results = pd.DataFrame(rows)

# ----------------------------
# Multiple-testing correction
# ----------------------------

results["fdr_bh"] = np.nan

mask = results["tested"] & results["p_value"].notna()

if mask.sum() == 0:
    raise RuntimeError("No PGFams passed the independent prefilter.")

results.loc[mask, "fdr_bh"] = multipletests(
    results.loc[mask, "p_value"].astype(float),
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

# Add BV-BRC product annotations.
results = results.merge(
    product_summary,
    left_on="pgfam_id",
    right_on="feature.pgfam_id",
    how="left",
).drop(columns=["feature.pgfam_id"], errors="ignore")

results["representative_product"] = (
    results["representative_product"].fillna("")
)

results["n_distinct_products"] = pd.to_numeric(
    results["n_distinct_products"], errors="coerce"
).fillna(0).astype(int)

# ----------------------------
# Save result tables
# ----------------------------

results = results.sort_values(
    ["tested", "fdr_bh", "prevalence_difference"],
    ascending=[False, True, False],
    na_position="last",
)

results.to_csv(
    OUT_ALL, sep="\t", index=False, encoding="utf-8-sig"
)

fdr_hits = results[
    results["fdr_significant"]
    & results["direction"].eq("AXILLA_ENRICHED")
].copy()

fdr_hits.to_csv(
    OUT_FDR, sep="\t", index=False, encoding="utf-8-sig"
)

nominal = results[
    results["tested"]
    & results["direction"].eq("AXILLA_ENRICHED")
    & (results["p_value"] < NOMINAL_P)
].copy()

nominal.to_csv(
    OUT_NOMINAL, sep="\t", index=False, encoding="utf-8-sig"
)

# Descriptive ranking only; not a significance threshold.
top = results[
    results["tested"]
    & results["direction"].eq("AXILLA_ENRICHED")
].copy()

top = top.sort_values(
    ["p_value", "prevalence_difference", "haldane_anscombe_or"],
    ascending=[True, False, False],
).head(50)

top.to_csv(
    OUT_TOP, sep="\t", index=False, encoding="utf-8-sig"
)

# ----------------------------
# Console summary
# ----------------------------

print("================================================")
print("E. COLI EXPLORATORY AXILLA ANALYSIS")
print("================================================")
print("Cohort:", n_ax, "AXILLA /", n_skin, "OTHER_SKIN")
print("AXILLA BioProjects:", ", ".join(ax_projects))
print("OTHER_SKIN BioProjects:", ", ".join(skin_projects))
print("Shared BioProjects:", len(shared))
print("PGFams observed:", len(results))
print(
    f"PGFams tested (present >= {MIN_PRESENT}, absent >= {MIN_ABSENT}):",
    int(results["tested"].sum()),
)
print("FDR < 0.05, AXILLA-enriched:", len(fdr_hits))
print("Nominal p < 0.05, AXILLA-enriched:", len(nominal))

show_cols = [
    "pgfam_id",
    "representative_product",
    "AXILLA_present",
    "OTHER_SKIN_present",
    "AXILLA_prevalence",
    "OTHER_SKIN_prevalence",
    "prevalence_difference",
    "haldane_anscombe_or",
    "p_value",
    "fdr_bh",
]

if len(fdr_hits):
    print("\nFDR-significant AXILLA-enriched candidates:")
    print(fdr_hits[show_cols].head(30).to_string(index=False))

if len(nominal):
    print("\nNominal AXILLA-enriched candidates:")
    print(nominal[show_cols].head(30).to_string(index=False))
else:
    print("\nNo nominal p<0.05 AXILLA-enriched candidates.")
    print("Top directional signals:")
    print(top[show_cols].head(20).to_string(index=False))

print("\nSaved:")
for f in [OUT_ALL, OUT_FDR, OUT_NOMINAL, OUT_TOP]:
    print(" ", f)

print(
    "\nINTERPRETATION NOTE:"
    "\nThis is a very small exploratory analysis (3 vs 5). "
    "AXILLA genomes all come from PRJNA1238219, whereas controls come from "
    "three different BioProjects. Therefore any apparent association can reflect "
    "lineage/study differences as well as body site. Use these results only as "
    "supporting evidence in the later cross-species functional synthesis."
)

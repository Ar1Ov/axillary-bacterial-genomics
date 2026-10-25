import math
from functools import lru_cache

import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

# ============================================================
# FINAL EXACT-SITE CORYNEBACTERIUM EXPLORATORY ANALYSIS
# ============================================================
# Final cohort:
#   C. amycolatum          2 AXILLA / 3 OTHER_SKIN
#   C. kefirresidentii     5 AXILLA / 4 OTHER_SKIN
#   C. tuberculostearicum  9 AXILLA / 1 OTHER_SKIN
#
# WHY THIS METHOD:
# The three species have very different case/control ratios, so a naive pooled
# Fisher test would be misleading. We therefore condition on species and test
# enrichment by an exact stratified test.
#
# For each PGFam, within each species:
#   X_s ~ Hypergeometric(N_s, K_s, n_axilla_s)
# under the null that body-site labels are exchangeable WITHIN species.
#
# The exact null distribution of total AXILLA carriers is obtained by
# convolving the three species-specific hypergeometric distributions.
# This automatically controls species composition and remains valid with the
# sparse 9-vs-1 C. tuberculostearicum stratum.
#
# FDR is applied to the exact TWO-SIDED p-values across all tested PGFams.
# One-sided enrichment p-values are also reported for ranking only.
#
# This remains EXPLORATORY because AXILLA and OTHER_SKIN genomes come from
# different BioProjects and several strata are small.
# ============================================================

COHORT_FILE = "cory_exact_site_final_cohort.tsv"
PRESENCE_FILE = "cory_exact_site_final_pgfam_presence_long.tsv"
FEATURE_FILE = "cory_expanded_features.tsv"

OUT_ALL = "cory_exact_site_stratified_results.tsv"
OUT_FDR = "cory_exact_site_fdr_axilla_candidates.tsv"
OUT_NOMINAL = "cory_exact_site_nominal_axilla.tsv"
OUT_CROSS = "cory_exact_site_cross_species_candidates.tsv"
OUT_TOP = "cory_exact_site_top_directional.tsv"

MIN_PRESENT = 2
MIN_ABSENT = 2
FDR_CUTOFF = 0.05
NOMINAL_P = 0.05

SPECIES_ORDER = [
    "Corynebacterium amycolatum",
    "Corynebacterium kefirresidentii",
    "Corynebacterium tuberculostearicum",
]

SHORT = {
    "Corynebacterium amycolatum": "amycolatum",
    "Corynebacterium kefirresidentii": "kefirresidentii",
    "Corynebacterium tuberculostearicum": "tuberculostearicum",
}

EXPECTED = {
    "Corynebacterium amycolatum": {"AXILLA": 2, "OTHER_SKIN": 3},
    "Corynebacterium kefirresidentii": {"AXILLA": 5, "OTHER_SKIN": 4},
    "Corynebacterium tuberculostearicum": {"AXILLA": 9, "OTHER_SKIN": 1},
}

# ============================================================
# LOAD + VALIDATE COHORT
# ============================================================

cohort = pd.read_csv(
    COHORT_FILE, sep="\t", dtype=str, keep_default_na=False
)

required = {"genome_id", "species", "bioproject", "site", "group"}
missing = required - set(cohort.columns)
if missing:
    raise ValueError(f"Missing cohort columns: {sorted(missing)}")

if cohort["genome_id"].duplicated().any():
    raise ValueError("Duplicate genome IDs found in cohort.")

if set(cohort["species"]) != set(SPECIES_ORDER):
    raise ValueError(
        f"Unexpected species: {sorted(set(cohort['species']))}"
    )

if set(cohort["group"]) != {"AXILLA", "OTHER_SKIN"}:
    raise ValueError(
        f"Unexpected groups: {sorted(set(cohort['group']))}"
    )

for species in SPECIES_ORDER:
    got = (
        cohort.loc[cohort["species"] == species, "group"]
        .value_counts()
        .to_dict()
    )
    for group, expected_n in EXPECTED[species].items():
        if got.get(group, 0) != expected_n:
            raise ValueError(
                f"{species}: expected {expected_n} {group}, "
                f"found {got.get(group, 0)}"
            )

if len(cohort) != 24:
    raise ValueError(f"Expected 24 genomes, found {len(cohort)}")

# BioProject structure
project_summary = (
    cohort.groupby(["species", "group"])["bioproject"]
    .agg(lambda x: sorted(set(x)))
)

# ============================================================
# LOAD + VALIDATE PRESENCE TABLE
# ============================================================

presence = pd.read_csv(
    PRESENCE_FILE, sep="\t", dtype=str, keep_default_na=False
)

required_p = {
    "genome_id", "feature.pgfam_id", "species", "bioproject", "site", "group"
}
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

if "feature.pgfam_id" not in features or "feature.product" not in features:
    raise ValueError("Feature file missing PGFam or product columns.")

cohort_feature_rows = features[
    features["genome_id"].isin(cohort_ids)
    & features["feature.pgfam_id"].ne("")
][["feature.pgfam_id", "feature.product"]].copy()

def top_product(s):
    vals = pd.Series([x.strip() for x in s if str(x).strip()])
    if vals.empty:
        return ""
    return vals.value_counts().index[0]

product_summary = (
    cohort_feature_rows.groupby("feature.pgfam_id")["feature.product"]
    .agg(representative_product=top_product, n_distinct_products="nunique")
    .reset_index()
)

# ============================================================
# EXACT STRATIFIED NULL DISTRIBUTION
# ============================================================

species_meta = {}

for species in SPECIES_ORDER:
    sub = cohort[cohort["species"] == species]
    species_meta[species] = {
        "N": len(sub),
        "n_ax": int((sub["group"] == "AXILLA").sum()),
        "n_skin": int((sub["group"] == "OTHER_SKIN").sum()),
        "ax_ids": set(sub.loc[sub["group"] == "AXILLA", "genome_id"]),
        "skin_ids": set(sub.loc[sub["group"] == "OTHER_SKIN", "genome_id"]),
    }

@lru_cache(maxsize=None)
def exact_null_distribution(k_tuple):
    """
    k_tuple = total carriers in each species, in SPECIES_ORDER.
    Returns a dict {total_axilla_carriers: probability} under the null.
    """
    pmf_total = np.array([1.0])

    for species, K in zip(SPECIES_ORDER, k_tuple):
        N = species_meta[species]["N"]
        n_ax = species_meta[species]["n_ax"]

        xmin = max(0, n_ax - (N - K))
        xmax = min(n_ax, K)
        xs = np.arange(xmin, xmax + 1)

        probs = hypergeom.pmf(xs, N, K, n_ax)

        species_vec = np.zeros(xmax + 1, dtype=float)
        species_vec[xs] = probs

        pmf_total = np.convolve(pmf_total, species_vec)

    # numerical normalization safeguard
    total = pmf_total.sum()
    if not np.isclose(total, 1.0, atol=1e-10):
        pmf_total = pmf_total / total

    return pmf_total

def exact_p_values(k_tuple, observed_ax_total):
    pmf = exact_null_distribution(k_tuple)

    if observed_ax_total >= len(pmf):
        raise RuntimeError("Observed AXILLA count outside exact null support.")

    p_obs = pmf[observed_ax_total]

    # Fisher-style exact two-sided p:
    # sum probabilities no larger than the probability of the observed table.
    p_two = float(pmf[pmf <= p_obs + 1e-12].sum())

    # Directional enrichment tail, for descriptive ranking only.
    p_enrich = float(pmf[observed_ax_total:].sum())

    # Depletion tail, retained for completeness.
    p_deplete = float(pmf[:observed_ax_total + 1].sum())

    return min(p_two, 1.0), min(p_enrich, 1.0), min(p_deplete, 1.0)

# ============================================================
# TEST EVERY PGFAM
# ============================================================

n_ax_total = int((cohort["group"] == "AXILLA").sum())
n_skin_total = int((cohort["group"] == "OTHER_SKIN").sum())
n_total = len(cohort)

rows = []

for pgfam, g in presence.groupby("feature.pgfam_id", sort=False):
    carriers = set(g["genome_id"])

    row = {"pgfam_id": pgfam}

    total_present = len(carriers)
    total_absent = n_total - total_present

    overall_ax_present = 0
    overall_skin_present = 0

    k_tuple = []
    positive_species = 0
    negative_species = 0
    equal_species = 0

    mh_num = 0.0
    mh_den = 0.0
    mh_num_corr = 0.0
    mh_den_corr = 0.0

    observed_ax_total = 0

    for species in SPECIES_ORDER:
        meta = species_meta[species]
        short = SHORT[species]

        ax_present = len(carriers & meta["ax_ids"])
        skin_present = len(carriers & meta["skin_ids"])

        ax_absent = meta["n_ax"] - ax_present
        skin_absent = meta["n_skin"] - skin_present

        ax_prev = ax_present / meta["n_ax"]
        skin_prev = skin_present / meta["n_skin"]
        diff = ax_prev - skin_prev

        row[f"{short}_AXILLA_present"] = ax_present
        row[f"{short}_OTHER_SKIN_present"] = skin_present
        row[f"{short}_AXILLA_prevalence"] = ax_prev
        row[f"{short}_OTHER_SKIN_prevalence"] = skin_prev
        row[f"{short}_prevalence_difference"] = diff

        if diff > 0:
            positive_species += 1
        elif diff < 0:
            negative_species += 1
        else:
            equal_species += 1

        overall_ax_present += ax_present
        overall_skin_present += skin_present
        observed_ax_total += ax_present

        K = ax_present + skin_present
        k_tuple.append(K)

        N = meta["N"]

        # Mantel-Haenszel common OR, uncorrected.
        mh_num += (ax_present * skin_absent) / N
        mh_den += (ax_absent * skin_present) / N

        # 0.5 continuity-corrected descriptive common OR.
        aa = ax_present + 0.5
        bb = ax_absent + 0.5
        cc = skin_present + 0.5
        dd = skin_absent + 0.5
        N_corr = N + 2.0

        mh_num_corr += (aa * dd) / N_corr
        mh_den_corr += (bb * cc) / N_corr

    tested = total_present >= MIN_PRESENT and total_absent >= MIN_ABSENT

    if tested:
        p_two, p_enrich, p_deplete = exact_p_values(
            tuple(k_tuple), observed_ax_total
        )
    else:
        p_two = np.nan
        p_enrich = np.nan
        p_deplete = np.nan

    overall_ax_prev = overall_ax_present / n_ax_total
    overall_skin_prev = overall_skin_present / n_skin_total
    crude_diff = overall_ax_prev - overall_skin_prev

    if mh_den == 0:
        mh_or = np.inf if mh_num > 0 else np.nan
    else:
        mh_or = mh_num / mh_den

    mh_or_corrected = (
        mh_num_corr / mh_den_corr if mh_den_corr > 0 else np.nan
    )

    row.update({
        "AXILLA_present_total": overall_ax_present,
        "OTHER_SKIN_present_total": overall_skin_present,
        "AXILLA_prevalence_total": overall_ax_prev,
        "OTHER_SKIN_prevalence_total": overall_skin_prev,
        "crude_prevalence_difference": crude_diff,
        "total_present": total_present,
        "total_absent": total_absent,
        "positive_species_count": positive_species,
        "negative_species_count": negative_species,
        "equal_species_count": equal_species,
        "tested": tested,
        "exact_two_sided_p": p_two,
        "exact_axilla_enrichment_p": p_enrich,
        "exact_axilla_depletion_p": p_deplete,
        "mh_common_or": mh_or,
        "mh_common_or_corrected": mh_or_corrected,
    })

    rows.append(row)

results = pd.DataFrame(rows)

# ============================================================
# MULTIPLE-TEST CORRECTION
# ============================================================

results["fdr_bh"] = np.nan

mask = results["tested"] & results["exact_two_sided_p"].notna()

if mask.sum() == 0:
    raise RuntimeError("No PGFams passed the independent prefilter.")

results.loc[mask, "fdr_bh"] = multipletests(
    results.loc[mask, "exact_two_sided_p"].astype(float),
    method="fdr_bh",
)[1]

# Overall direction is based on the species-stratified evidence pattern,
# not merely the crude pooled prevalence.
results["direction"] = np.select(
    [
        (results["positive_species_count"] > 0)
        & (results["negative_species_count"] == 0),
        (results["negative_species_count"] > 0)
        & (results["positive_species_count"] == 0),
    ],
    [
        "AXILLA_ENRICHED_CONSISTENT",
        "OTHER_SKIN_ENRICHED_CONSISTENT",
    ],
    default="MIXED_OR_EQUAL",
)

results["fdr_significant"] = (
    results["tested"]
    & results["fdr_bh"].notna()
    & (results["fdr_bh"] < FDR_CUTOFF)
)

results["cross_species_consistent"] = (
    (results["positive_species_count"] >= 2)
    & (results["negative_species_count"] == 0)
)

# Add BV-BRC product annotation.
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

results = results.sort_values(
    [
        "tested",
        "fdr_bh",
        "exact_axilla_enrichment_p",
        "positive_species_count",
        "mh_common_or_corrected",
    ],
    ascending=[False, True, True, False, False],
    na_position="last",
)

results.to_csv(
    OUT_ALL, sep="\t", index=False, encoding="utf-8-sig"
)

fdr_hits = results[
    results["fdr_significant"]
    & (results["positive_species_count"] > 0)
    & (results["negative_species_count"] == 0)
].copy()

fdr_hits.to_csv(
    OUT_FDR, sep="\t", index=False, encoding="utf-8-sig"
)

nominal = results[
    results["tested"]
    & (results["exact_two_sided_p"] < NOMINAL_P)
    & (results["positive_species_count"] > 0)
    & (results["negative_species_count"] == 0)
].copy()

nominal.to_csv(
    OUT_NOMINAL, sep="\t", index=False, encoding="utf-8-sig"
)

cross = results[
    results["tested"]
    & results["cross_species_consistent"]
].copy()

cross = cross.sort_values(
    [
        "fdr_bh",
        "exact_axilla_enrichment_p",
        "positive_species_count",
        "mh_common_or_corrected",
    ],
    ascending=[True, True, False, False],
    na_position="last",
)

cross.to_csv(
    OUT_CROSS, sep="\t", index=False, encoding="utf-8-sig"
)

top = results[
    results["tested"]
    & (results["positive_species_count"] > 0)
    & (results["negative_species_count"] == 0)
].head(75)

top.to_csv(
    OUT_TOP, sep="\t", index=False, encoding="utf-8-sig"
)

# ============================================================
# SUMMARY
# ============================================================

print("================================================")
print("CORYNEBACTERIUM EXACT-SITE STRATIFIED ANALYSIS")
print("================================================")
print("Cohort: 16 AXILLA / 8 OTHER_SKIN across 3 species")
print()
for species in SPECIES_ORDER:
    a = EXPECTED[species]["AXILLA"]
    s = EXPECTED[species]["OTHER_SKIN"]
    print(f"{species}: {a} AXILLA / {s} OTHER_SKIN")

print("\nBioProject structure:")
for species in SPECIES_ORDER:
    print(
        species,
        "| AXILLA:", project_summary.loc[(species, "AXILLA")],
        "| OTHER_SKIN:", project_summary.loc[(species, "OTHER_SKIN")]
    )

print("\nPGFams observed:", len(results))
print(
    f"PGFams tested (present >= {MIN_PRESENT}, absent >= {MIN_ABSENT}):",
    int(results["tested"].sum()),
)
print(
    "FDR < 0.05, consistently AXILLA-enriched:",
    len(fdr_hits),
)
print(
    "Nominal exact two-sided p < 0.05, consistently AXILLA-enriched:",
    len(nominal),
)
print(
    "Cross-species-consistent directional candidates "
    "(positive in >=2 species, negative in 0):",
    int(results["cross_species_consistent"].sum()),
)

show_cols = [
    "pgfam_id",
    "representative_product",
    "AXILLA_present_total",
    "OTHER_SKIN_present_total",
    "positive_species_count",
    "negative_species_count",
    "mh_common_or_corrected",
    "exact_two_sided_p",
    "exact_axilla_enrichment_p",
    "fdr_bh",
    "amycolatum_prevalence_difference",
    "kefirresidentii_prevalence_difference",
    "tuberculostearicum_prevalence_difference",
]

if len(fdr_hits):
    print("\nFDR-significant consistently AXILLA-enriched candidates:")
    print(fdr_hits[show_cols].head(30).to_string(index=False))

if len(nominal):
    print("\nTop nominal consistently AXILLA-enriched candidates:")
    print(nominal[show_cols].head(30).to_string(index=False))
else:
    print("\nNo nominal exact two-sided p<0.05 consistently AXILLA-enriched candidates.")

print("\nTop cross-species-consistent directional candidates:")
print(cross[show_cols].head(30).to_string(index=False))

print("\nSaved:")
for f in [OUT_ALL, OUT_FDR, OUT_NOMINAL, OUT_CROSS, OUT_TOP]:
    print(" ", f)

print(
    "\nINTERPRETATION NOTE:"
    "\nThis exact test controls SPECIES by permuting body-site labels only within "
    "species, but it does not remove BioProject/lineage confounding because "
    "AXILLA and OTHER_SKIN genomes are project-separated. The Cory analysis "
    "therefore remains a secondary exploratory arm."
)

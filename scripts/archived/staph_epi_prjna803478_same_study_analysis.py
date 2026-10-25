import math
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

COHORT_FILE = "staph_epi_prjna803478_final_cohort.tsv"
PRESENCE_FILE = "staph_epi_prjna803478_final_pgfam_presence_long.tsv"
FEATURE_FILE = "staph_epi_prjna803478_features.tsv"
STRICT_RESULTS_FILE = "staph_epi_strict_primary_pgfam_results.tsv"

OUT_ALL = "staph_epi_prjna803478_same_study_results.tsv"
OUT_FDR = "staph_epi_prjna803478_same_study_fdr_axilla_candidates.tsv"
OUT_NOMINAL = "staph_epi_prjna803478_same_study_nominal_axilla.tsv"
OUT_CONCORDANT = "staph_epi_cross_cohort_concordant_candidates.tsv"

MIN_PRESENT = 2
MIN_ABSENT = 2
FDR_CUTOFF = 0.05
NOMINAL_P = 0.05

cohort = pd.read_csv(COHORT_FILE, sep="\t", dtype=str, keep_default_na=False)

required = {"genome_id", "bioproject", "site", "group"}
missing = required - set(cohort.columns)
if missing:
    raise ValueError(f"Missing cohort columns: {sorted(missing)}")

if cohort["genome_id"].duplicated().any():
    raise ValueError("Duplicate genome IDs found in cohort.")

counts = cohort["group"].value_counts().to_dict()
if counts.get("AXILLARY_GROIN", 0) != 4:
    raise ValueError(f"Expected 4 AXILLARY_GROIN genomes, found {counts.get('AXILLARY_GROIN', 0)}")
if counts.get("OTHER_SKIN", 0) != 22:
    raise ValueError(f"Expected 22 OTHER_SKIN genomes, found {counts.get('OTHER_SKIN', 0)}")

if set(cohort["bioproject"]) != {"PRJNA803478"}:
    raise ValueError(f"Expected only PRJNA803478, found {sorted(set(cohort['bioproject']))}")

presence = pd.read_csv(PRESENCE_FILE, sep="\t", dtype=str, keep_default_na=False)

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
        f"Presence/cohort genome mismatch. Missing: {sorted(cohort_ids-presence_ids)[:10]}; "
        f"extra: {sorted(presence_ids-cohort_ids)[:10]}"
    )

features = pd.read_csv(
    FEATURE_FILE, sep="\t", dtype=str, keep_default_na=False, low_memory=False
)

if "feature.pgfam_id" not in features or "feature.product" not in features:
    raise ValueError("Feature file is missing PGFam or product columns.")

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

ax_ids = set(cohort.loc[cohort["group"] == "AXILLARY_GROIN", "genome_id"])
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

    tested = total_present >= MIN_PRESENT and total_absent >= MIN_ABSENT

    if tested:
        crude_or, p = fisher_exact(
            [[ax_present, ax_absent], [skin_present, skin_absent]],
            alternative="two-sided",
        )
    else:
        crude_or, p = np.nan, np.nan

    aa, bb, cc, dd = ax_present + 0.5, ax_absent + 0.5, skin_present + 0.5, skin_absent + 0.5
    corrected_or = (aa * dd) / (bb * cc)
    se = math.sqrt(1/aa + 1/bb + 1/cc + 1/dd)
    log_or = math.log(corrected_or)
    ci_low = math.exp(log_or - 1.96 * se)
    ci_high = math.exp(log_or + 1.96 * se)

    rows.append({
        "pgfam_id": pgfam,
        "AXILLARY_GROIN_present": ax_present,
        "AXILLARY_GROIN_absent": ax_absent,
        "OTHER_SKIN_present": skin_present,
        "OTHER_SKIN_absent": skin_absent,
        "AXILLARY_GROIN_prevalence": ax_prev,
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
results["fdr_bh"] = np.nan

mask = results["tested"] & results["p_value"].notna()
if mask.sum() == 0:
    raise RuntimeError("No PGFams passed the independent prefilter.")

results.loc[mask, "fdr_bh"] = multipletests(
    results.loc[mask, "p_value"].astype(float),
    method="fdr_bh",
)[1]

results["direction"] = np.select(
    [results["prevalence_difference"] > 0, results["prevalence_difference"] < 0],
    ["AXILLARY_GROIN_ENRICHED", "OTHER_SKIN_ENRICHED"],
    default="EQUAL",
)

results["fdr_significant"] = (
    results["tested"]
    & results["fdr_bh"].notna()
    & (results["fdr_bh"] < FDR_CUTOFF)
)

results = results.merge(
    product_summary,
    left_on="pgfam_id",
    right_on="feature.pgfam_id",
    how="left",
).drop(columns=["feature.pgfam_id"], errors="ignore")

results["representative_product"] = results["representative_product"].fillna("")
results["n_distinct_products"] = pd.to_numeric(
    results["n_distinct_products"], errors="coerce"
).fillna(0).astype(int)

strict = pd.read_csv(STRICT_RESULTS_FILE, sep="\t", dtype=str, keep_default_na=False)

strict_cols = [
    "pgfam_id",
    "AXILLA_present",
    "OTHER_SKIN_present",
    "AXILLA_prevalence",
    "OTHER_SKIN_prevalence",
    "prevalence_difference",
    "haldane_anscombe_or",
    "p_value",
    "fdr_bh",
    "AXILLA_projects_with_gene",
    "OTHER_SKIN_projects_with_gene",
    "direction",
]

missing = set(strict_cols) - set(strict.columns)
if missing:
    raise ValueError(f"Strict-results file missing columns: {sorted(missing)}")

strict = strict[strict_cols].copy()
strict = strict.rename(columns={c: f"strict_{c}" for c in strict_cols if c != "pgfam_id"})

for c in strict.columns:
    if c != "pgfam_id" and c != "strict_direction":
        strict[c] = pd.to_numeric(strict[c], errors="coerce")

results = results.merge(strict, on="pgfam_id", how="left")

results["same_direction_strict"] = (
    (results["direction"] == "AXILLARY_GROIN_ENRICHED")
    & (results["strict_prevalence_difference"] > 0)
)

results["strict_recurred_3plus_axilla_projects"] = (
    results["strict_AXILLA_projects_with_gene"] >= 3
)

results = results.sort_values(
    ["tested", "fdr_bh", "prevalence_difference"],
    ascending=[False, True, False],
    na_position="last",
)

results.to_csv(OUT_ALL, sep="\t", index=False, encoding="utf-8-sig")

fdr_hits = results[
    results["fdr_significant"]
    & results["direction"].eq("AXILLARY_GROIN_ENRICHED")
].copy()
fdr_hits.to_csv(OUT_FDR, sep="\t", index=False, encoding="utf-8-sig")

nominal = results[
    results["tested"]
    & results["direction"].eq("AXILLARY_GROIN_ENRICHED")
    & (results["p_value"] < NOMINAL_P)
].copy()
nominal.to_csv(OUT_NOMINAL, sep="\t", index=False, encoding="utf-8-sig")

concordant = nominal[nominal["same_direction_strict"]].copy()
concordant = concordant.sort_values(
    ["p_value", "strict_p_value", "prevalence_difference"],
    ascending=[True, True, False],
)
concordant.to_csv(OUT_CONCORDANT, sep="\t", index=False, encoding="utf-8-sig")

print("================================================")
print("S. EPIDERMIDIS PRJNA803478 SAME-STUDY ANALYSIS")
print("================================================")
print("Cohort:", n_ax, "AXILLARY_GROIN /", n_skin, "OTHER_SKIN")
print("BioProject: PRJNA803478")
print("PGFams observed:", len(results))
print(f"PGFams tested (present >= {MIN_PRESENT}, absent >= {MIN_ABSENT}):", int(results["tested"].sum()))
print("FDR < 0.05, AXILLARY_GROIN-enriched:", len(fdr_hits))
print("Nominal p < 0.05, AXILLARY_GROIN-enriched:", len(nominal))
print("Nominal same-study + same-direction strict cohort:", len(concordant))

show_cols = [
    "pgfam_id",
    "representative_product",
    "AXILLARY_GROIN_present",
    "OTHER_SKIN_present",
    "AXILLARY_GROIN_prevalence",
    "OTHER_SKIN_prevalence",
    "haldane_anscombe_or",
    "p_value",
    "fdr_bh",
    "strict_AXILLA_present",
    "strict_OTHER_SKIN_present",
    "strict_prevalence_difference",
    "strict_p_value",
    "strict_fdr_bh",
    "strict_AXILLA_projects_with_gene",
]

if len(fdr_hits):
    print("\nFDR-significant AXILLARY_GROIN-enriched candidates:")
    print(fdr_hits[show_cols].head(30).to_string(index=False))

if len(concordant):
    print("\nExploratory cross-cohort concordant candidates:")
    print(concordant[show_cols].head(30).to_string(index=False))
elif len(nominal):
    print("\nTop nominal same-study AXILLARY_GROIN-enriched candidates:")
    print(nominal[show_cols].head(30).to_string(index=False))
else:
    top = results[
        results["tested"]
        & results["direction"].eq("AXILLARY_GROIN_ENRICHED")
    ].head(20)
    print("\nNo nominal p<0.05 AXILLARY_GROIN hits. Top directional signals:")
    print(top[show_cols].to_string(index=False))

print("\nSaved:")
for f in [OUT_ALL, OUT_FDR, OUT_NOMINAL, OUT_CONCORDANT]:
    print(" ", f)

print(
    "\nINTERPRETATION NOTE:"
    "\nFDR significance remains the genome-wide significance criterion. "
    "The cross-cohort concordant file is exploratory only: nominal p<0.05 "
    "is not treated as genome-wide significance."
)

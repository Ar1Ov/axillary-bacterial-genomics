import re
import unicodedata
from pathlib import Path

import pandas as pd

# ============================================================
# PURPOSE
# ============================================================
# Objective, all-human-bacteria species-selection screen for:
#
# "Within human skin bacterial species, which genes are enriched
#  in axillary isolates compared with non-axillary skin isolates?"
#
# This script DOES NOT select Corynebacterium or Staphylococcus in
# advance. It begins with the all-human-bacterial-genome BV-BRC
# master file and determines which named species have enough
# STRICT_AXILLA and explicitly NON-AXILLARY OTHER_SKIN genomes.
#
# IMPORTANT DESIGN CHOICE:
# Generic metadata such as "skin" or "cutaneous" is NOT counted as
# OTHER_SKIN, because it does not prove the sample was non-axillary.
# It is placed in SKIN_UNSPECIFIED.
#
# Comments are not allowed to silently override structured site
# metadata. Conflicts are flagged for manual review.
# ============================================================

INPUT = "human_bacteria_master.tsv"

OUT_ALL = "human_bacteria_master_classified_v4.tsv"
OUT_SPECIES = "human_bacteria_species_eligibility_v4.tsv"
OUT_PROJECTS = "human_bacteria_species_project_counts_v4.tsv"
OUT_AXILLA_REVIEW = "human_bacteria_axilla_candidates_review_v4.tsv"
OUT_AXILLA_LABELS = "human_bacteria_axilla_unique_labels_v4.tsv"
OUT_ELIGIBLE_AUDIT = "human_bacteria_eligible_species_audit_v4.tsv"
OUT_DUPLICATES = "human_bacteria_exact_duplicate_rows_v4.tsv"

QC_COMPLETENESS = 95.0
QC_CONTAMINATION = 5.0

# Pragmatic pre-specified screening thresholds.
# These are eligibility thresholds, NOT guarantees of statistical power.
PRIMARY_MIN_STRICT_AXILLA = 5
PRIMARY_MIN_OTHER_SKIN = 10

EXPLORATORY_MIN_STRICT_AXILLA = 3
EXPLORATORY_MIN_OTHER_SKIN = 5

STRUCTURED_SITE_COLS = [
    "genome.body_sample_site",
    "genome.body_sample_subsite",
    "genome.isolation_site",
    "genome.isolation_source",
]
COMMENTS_COL = "genome.isolation_comments"

REQUIRED = {
    "genome.genome_id",
    "genome.bioproject_accession",
    "genome.species",
    *STRUCTURED_SITE_COLS,
    COMMENTS_COL,
    "genome.checkm_completeness",
    "genome.checkm_contamination",
}


# ============================================================
# NORMALIZATION
# ============================================================

def norm_text(x):
    if x is None or pd.isna(x):
        return ""
    s = str(x)
    s = unicodedata.normalize("NFKC", s)
    s = (
        s.replace("\u2013", "-")
         .replace("\u2014", "-")
         .replace("\u2212", "-")
    )
    s = s.lower().strip()
    s = re.sub(r"\s+", " ", s)
    return s


# ============================================================
# BODY-SITE DEFINITIONS — V4 CONSERVATIVE RULES
# ============================================================

AXILLA_RE = re.compile(
    r"\b(?:axilla|axillae|axillary|armpit|underarm)\b",
    re.I,
)

MAXILLARY_RE = re.compile(
    r"\b(?:maxilla|maxillary|maxillofacial)\b",
    re.I,
)

MIXED_AXILLA_RE = re.compile(
    r"\b(?:"
    r"neck|groin|inguinal|umbilicus|umbilical|"
    r"perineum|perineal|multisite|multi-site"
    r")\b",
    re.I,
)

CLINICAL_RE = re.compile(
    r"\b(?:"
    r"abscess|pus|purulent|lesion|wound|ulcer|ulceration|"
    r"dermatitis|eczema|cellulitis|infection|infected|"
    r"lymph\s*node|lymphnode|soft\s*tissue|"
    r"surgical|postoperative|post-operative|incision|"
    r"biopsy|tumou?r|cyst|"
    r"catheter|transplant|bone\s*marrow"
    r")\b",
    re.I,
)

# Specimen types that are not evidence of a normal skin-surface niche.
NON_SURFACE_SPECIMEN_RE = re.compile(
    r"\b(?:"
    r"intra-abdominal|intraabdominal|"
    r"fluid|secretion|exudate|"
    r"gauze|drain|jp\s*drain|"
    r"aspirate|aspiration|"
    r"tissue|deep\s*tissue"
    r")\b",
    re.I,
)

NON_SKIN_RE = re.compile(
    r"\b(?:"
    r"blood|serum|plasma|urine|urinary|urethra|urethral|"
    r"vagina|vaginal|cervix|cervical|rectum|rectal|"
    r"feces|faeces|fecal|faecal|stool|"
    r"mouth|oral|saliva|tongue|gingiva|gingival|tooth|dental|"
    r"nasal|nares|nostril|nasopharyn\w*|oropharyn\w*|"
    r"pharyn\w*|tonsil\w*|sputum|trache\w*|lung|bronch\w*|"
    r"sinus|"
    r"eye|conjunctiva|tear|ear\s*canal|"
    r"csf|cerebrospinal|brain|heart|cardiac|"
    r"gut|intestinal|intestine|colon|ileum|jejunum|duodenum|"
    r"semen|prostate|placenta|amniotic|"
    r"bone|joint|synovial|liver|kidney|renal"
    r")\b",
    re.I,
)

OTHER_APOCRINE_RE = re.compile(
    r"\b(?:groin|inguinal|perineum|perineal)\b",
    re.I,
)

# Highly recognizable skin-surface niches. These can count as OTHER_SKIN
# even when the word "skin" is omitted.
SPECIFIC_SKIN_NICHE_RE = re.compile(
    r"\b(?:"
    r"forehead|"
    r"antecubital(?:\s+fossa)?|"
    r"occiput|scalp|"
    r"toe\s*web(?:\s*space)?|interdigital|"
    r"alar\s+crease|retroauricular|"
    r"umbilicus|umbilical|navel"
    r")\b",
    re.I,
)

# Body regions that can represent skin, but need surface-sampling evidence
# to count automatically as explicit OTHER_SKIN.
SURFACE_REGION_RE = re.compile(
    r"\b(?:"
    r"forearm|cheek|face|facial|temple|"
    r"hand|palm|finger|wrist|"
    r"arm|elbow|shoulder|"
    r"leg|thigh|knee|shin|calf|"
    r"foot|heel|plantar|sole|"
    r"buttock|gluteal"
    r")\b",
    re.I,
)

# These regions are especially ambiguous because they may refer to internal,
# surgical, lymphatic, or other non-surface specimens.
AMBIGUOUS_REGION_RE = re.compile(
    r"\b(?:"
    r"abdomen|abdominal|chest|back|neck|torso|trunk"
    r")\b",
    re.I,
)

SURFACE_EVIDENCE_RE = re.compile(
    r"\b(?:"
    r"skin|cutaneous|swab|surface|scraping|wipe|tape\s*strip"
    r")\b",
    re.I,
)

GENERIC_SKIN_RE = re.compile(
    r"\b(?:skin|cutaneous)\b",
    re.I,
)


# ============================================================
# SPECIES VALIDITY
# ============================================================

INVALID_SPECIES_RE = re.compile(
    r"(?:"
    r"\bsp\.(?:\s|$)|"
    r"\buncultured\b|"
    r"\bunclassified\b|"
    r"\bmetagenome\b|"
    r"\bbacterium\b|"
    r"\bgroup\b"
    r")",
    re.I,
)


def valid_named_species(species):
    s = norm_text(species)
    if not s:
        return False
    if INVALID_SPECIES_RE.search(s):
        return False
    # Named bacterial species should contain at least genus + species epithet.
    return len(s.split()) >= 2


# ============================================================
# CLASSIFIER — V4
# ============================================================

def classify_text(text):
    """
    Conservative body-site classification for the study's phenotype.
    """
    text = norm_text(text)

    if not text:
        return "UNKNOWN"

    has_axilla = bool(AXILLA_RE.search(text))
    has_maxillary = bool(MAXILLARY_RE.search(text))
    has_clinical = bool(CLINICAL_RE.search(text))
    has_non_surface = bool(NON_SURFACE_SPECIMEN_RE.search(text))
    has_non_skin = bool(NON_SKIN_RE.search(text))
    has_other_apocrine = bool(OTHER_APOCRINE_RE.search(text))
    has_specific_skin = bool(SPECIFIC_SKIN_NICHE_RE.search(text))
    has_surface_region = bool(SURFACE_REGION_RE.search(text))
    has_ambiguous_region = bool(AMBIGUOUS_REGION_RE.search(text))
    has_surface_evidence = bool(SURFACE_EVIDENCE_RE.search(text))
    has_generic_skin = bool(GENERIC_SKIN_RE.search(text))

    if has_maxillary and not has_axilla:
        return "OTHER_BODY_SITE"

    if has_axilla:
        if has_clinical:
            return "EXCLUDE_CLINICAL"
        if has_non_surface:
            return "EXCLUDE_NON_SURFACE"
        if has_non_skin:
            return "AMBIGUOUS"
        if MIXED_AXILLA_RE.search(text):
            return "BROAD_AXILLA"
        return "STRICT_AXILLA"

    if has_clinical:
        return "EXCLUDE_CLINICAL"

    if has_non_surface:
        return "EXCLUDE_NON_SURFACE"

    if has_other_apocrine:
        return "OTHER_APOCRINE_SITE"

    if has_specific_skin:
        if has_non_skin:
            return "AMBIGUOUS"
        return "OTHER_SKIN"

    if has_surface_region:
        if has_non_skin:
            return "AMBIGUOUS"
        if has_surface_evidence:
            return "OTHER_SKIN"
        return "SKIN_SITE_UNCONFIRMED"

    if has_ambiguous_region:
        if has_non_skin:
            return "OTHER_BODY_SITE"
        if has_surface_evidence:
            return "OTHER_SKIN_POSSIBLE"
        return "SKIN_SITE_UNCONFIRMED"

    if has_generic_skin:
        if has_non_skin:
            return "AMBIGUOUS"
        return "SKIN_UNSPECIFIED"

    return "OTHER_BODY_SITE"


# ============================================================
# BUILT-IN SELF TESTS
# ============================================================

def run_self_tests():
    expected = {
        "axilla": "STRICT_AXILLA",
        "left axilla swab": "STRICT_AXILLA",
        "bilateral axilla": "STRICT_AXILLA",
        "armpit swab": "STRICT_AXILLA",
        "Axillary-groin Swab": "BROAD_AXILLA",
        "Neck/Axilla": "BROAD_AXILLA",
        "axilla/inguinal swab": "BROAD_AXILLA",
        "maxillary sinus": "OTHER_BODY_SITE",
        "tropical superficial axillary dermatitis": "EXCLUDE_CLINICAL",
        "Axillary lymph node": "EXCLUDE_CLINICAL",
        "Pus left axilla": "EXCLUDE_CLINICAL",
        "axillary blood": "AMBIGUOUS",
        "forehead swab": "OTHER_SKIN",
        "toe web space swab": "OTHER_SKIN",
        "alar crease swab": "OTHER_SKIN",
        "umbilicus": "OTHER_SKIN",
        "hand swab": "OTHER_SKIN",
        "skin of the hand": "OTHER_SKIN",
        "skin lower leg": "OTHER_SKIN",
        "farm worker hand swab": "OTHER_SKIN",
        "right arm": "SKIN_SITE_UNCONFIRMED",
        "neck": "SKIN_SITE_UNCONFIRMED",
        "neck swab": "OTHER_SKIN_POSSIBLE",
        "swab, abdomen": "OTHER_SKIN_POSSIBLE",
        "abdominal": "SKIN_SITE_UNCONFIRMED",
        "abdominal gauze": "EXCLUDE_NON_SURFACE",
        "abdominal fluid": "EXCLUDE_NON_SURFACE",
        "intra-abdominal": "EXCLUDE_NON_SURFACE",
        "abdominal fluid jp drain": "EXCLUDE_NON_SURFACE",
        "swab, abdominal incision": "EXCLUDE_CLINICAL",
        "foot secretion": "EXCLUDE_NON_SURFACE",
        "leg secretion": "EXCLUDE_NON_SURFACE",
        "skin": "SKIN_UNSPECIFIED",
        "cutaneous": "SKIN_UNSPECIFIED",
        "skin swab": "SKIN_UNSPECIFIED",
        "groin swab": "OTHER_APOCRINE_SITE",
        "inguinal swab": "OTHER_APOCRINE_SITE",
        "nares swab": "OTHER_BODY_SITE",
        "blood": "OTHER_BODY_SITE",
        "skin wound": "EXCLUDE_CLINICAL",
        "forearm blood": "AMBIGUOUS",
        "": "UNKNOWN",
    }

    failures = []
    for label, want in expected.items():
        got = classify_text(label)
        if got != want:
            failures.append((label, want, got))

    if failures:
        lines = ["Built-in classifier self-test FAILED:"]
        for label, want, got in failures:
            lines.append(f"  {label!r}: expected {want}, got {got}")
        raise RuntimeError("\n".join(lines))

    print(f"Classifier self-tests passed: {len(expected)} / {len(expected)}")


# ============================================================
# LOAD DATA
# ============================================================

run_self_tests()

df = pd.read_csv(
    INPUT,
    sep="\t",
    dtype=str,
    keep_default_na=False,
    low_memory=False,
    encoding="utf-8",
    encoding_errors="replace",
)

missing = REQUIRED - set(df.columns)
if missing:
    raise ValueError(f"Missing required columns: {sorted(missing)}")

print("Rows loaded:", len(df))

# Genome IDs should represent unique genomes. BV-BRC exports can contain
# repeated identical rows, so handle them conservatively:
#   1) identify every repeated genome ID;
#   2) verify that all metadata are identical within each repeated ID;
#   3) save those rows for audit;
#   4) collapse only exact duplicates;
#   5) STOP if any repeated genome ID has conflicting metadata.
dups = df[df["genome.genome_id"].duplicated(keep=False)].copy()

if not dups.empty:
    duplicate_id_count = dups["genome.genome_id"].nunique()

    # Count how many distinct values each column has within each duplicated genome ID.
    per_id_nunique = (
        dups.groupby("genome.genome_id", dropna=False)
        .nunique(dropna=False)
    )

    conflicting_ids = per_id_nunique.index[
        per_id_nunique.max(axis=1) > 1
    ].tolist()

    if conflicting_ids:
        raise ValueError(
            f"Found {len(conflicting_ids)} duplicated genome IDs with conflicting "
            "metadata. Stop and inspect them before species counting. "
            f"Example IDs: {conflicting_ids[:10]}"
        )

    dups.to_csv(
        OUT_DUPLICATES,
        sep="\t",
        index=False,
        encoding="utf-8-sig",
    )

    before = len(df)
    df = df.drop_duplicates().copy()
    removed = before - len(df)

    # Final assertion: after collapsing exact duplicate rows, genome IDs must be unique.
    if df["genome.genome_id"].duplicated().any():
        raise RuntimeError(
            "Exact-row deduplication did not produce unique genome IDs. "
            "Stop before continuing."
        )

    print(
        f"Exact duplicate audit: {duplicate_id_count} duplicated genome IDs; "
        f"{removed} repeated rows safely collapsed."
    )
    print("Unique genome IDs after deduplication:", len(df))
else:
    print("No duplicate genome IDs found.")



# ============================================================
# QC
# ============================================================

df["checkm_completeness_num"] = pd.to_numeric(
    df["genome.checkm_completeness"], errors="coerce"
)
df["checkm_contamination_num"] = pd.to_numeric(
    df["genome.checkm_contamination"], errors="coerce"
)

df["qc_pass"] = (
    (df["checkm_completeness_num"] >= QC_COMPLETENESS)
    & (df["checkm_contamination_num"] <= QC_CONTAMINATION)
)

df["valid_named_species"] = df["genome.species"].map(valid_named_species)

print("QC-pass genomes:", int(df["qc_pass"].sum()))
print("QC-fail/missing-quality genomes:", int((~df["qc_pass"]).sum()))
print("QC-pass genomes with named species:",
      int((df["qc_pass"] & df["valid_named_species"]).sum()))


# ============================================================
# NORMALIZE STRUCTURED FIELDS + COMMENTS
# ============================================================

normalized = {}
for col in STRUCTURED_SITE_COLS + [COMMENTS_COL]:
    normalized[col] = df[col].map(norm_text)

# Keep the structured metadata separate from free-text comments.
df["structured_site_text"] = (
    normalized[STRUCTURED_SITE_COLS[0]]
    + " || " + normalized[STRUCTURED_SITE_COLS[1]]
    + " || " + normalized[STRUCTURED_SITE_COLS[2]]
    + " || " + normalized[STRUCTURED_SITE_COLS[3]]
).str.strip(" |")

df["comments_text"] = normalized[COMMENTS_COL]

structured_nonempty = df["structured_site_text"].str.replace(
    "|", "", regex=False
).str.strip().ne("")

# Structured metadata is authoritative when available.
df["classification_source"] = "STRUCTURED"
df.loc[~structured_nonempty & df["comments_text"].ne(""),
       "classification_source"] = "COMMENTS_ONLY"
df.loc[~structured_nonempty & df["comments_text"].eq(""),
       "classification_source"] = "NO_SITE_METADATA"

df["site_class"] = [
    classify_text(structured if has_struct else comments)
    for structured, comments, has_struct in zip(
        df["structured_site_text"],
        df["comments_text"],
        structured_nonempty,
    )
]

# Independently classify comments to detect contradictions that should be reviewed.
df["comments_class"] = df["comments_text"].map(classify_text)

df["comments_have_axilla"] = df["comments_text"].str.contains(
    AXILLA_RE, regex=True, na=False
)
df["structured_have_axilla"] = df["structured_site_text"].str.contains(
    AXILLA_RE, regex=True, na=False
)

# Manual-review triggers:
# 1) classification depended only on comments;
# 2) structured metadata and comments disagree in a biologically meaningful way;
# 3) structured "good" skin/axilla record has clinical wording in comments.
meaningful_classes = {
    "STRICT_AXILLA",
    "BROAD_AXILLA",
    "OTHER_SKIN",
    "OTHER_SKIN_POSSIBLE",
    "SKIN_SITE_UNCONFIRMED",
    "SKIN_UNSPECIFIED",
    "OTHER_APOCRINE_SITE",
    "EXCLUDE_CLINICAL",
    "EXCLUDE_NON_SURFACE",
    "AMBIGUOUS",
}

df["needs_manual_review"] = (
    df["classification_source"].eq("COMMENTS_ONLY")
    | (
        df["comments_text"].ne("")
        & df["comments_class"].isin(meaningful_classes)
        & df["site_class"].isin(meaningful_classes)
        & df["comments_class"].ne(df["site_class"])
    )
    | (
        df["site_class"].isin(["STRICT_AXILLA", "BROAD_AXILLA", "OTHER_SKIN"])
        & df["comments_text"].str.contains(CLINICAL_RE, regex=True, na=False)
    )
)

# Any axilla wording anywhere should be available for audit.
df["any_axilla_term"] = (
    df["structured_have_axilla"] | df["comments_have_axilla"]
)


# ============================================================
# SAVE FULL CLASSIFIED MASTER
# ============================================================

df.to_csv(
    OUT_ALL,
    sep="\t",
    index=False,
    encoding="utf-8-sig",
)


# ============================================================
# AXILLA AUDIT
# ============================================================

ax_review = df[
    df["qc_pass"]
    & df["valid_named_species"]
    & df["any_axilla_term"]
].copy()

ax_review.to_csv(
    OUT_AXILLA_REVIEW,
    sep="\t",
    index=False,
    encoding="utf-8-sig",
)

ax_labels = (
    ax_review.groupby(
        [
            "site_class",
            "classification_source",
            "structured_site_text",
            "comments_text",
        ],
        dropna=False,
    )
    .size()
    .reset_index(name="count")
    .sort_values(
        ["site_class", "count"],
        ascending=[True, False],
    )
)

ax_labels.to_csv(
    OUT_AXILLA_LABELS,
    sep="\t",
    index=False,
    encoding="utf-8-sig",
)


# ============================================================
# AUTO-USABLE RECORDS FOR OBJECTIVE SPECIES SCREEN
# ============================================================
# We use only:
#   - QC-passed genomes
#   - confidently named species
#   - records not requiring manual review
#
# Manual-review records remain available in the audit files and may be
# restored only after explicit inspection.

usable = df[
    df["qc_pass"]
    & df["valid_named_species"]
    & ~df["needs_manual_review"]
].copy()

classes = [
    "STRICT_AXILLA",
    "BROAD_AXILLA",
    "OTHER_SKIN",
    "OTHER_SKIN_POSSIBLE",
    "SKIN_SITE_UNCONFIRMED",
    "SKIN_UNSPECIFIED",
    "OTHER_APOCRINE_SITE",
    "EXCLUDE_CLINICAL",
    "EXCLUDE_NON_SURFACE",
    "AMBIGUOUS",
    "OTHER_BODY_SITE",
    "UNKNOWN",
]

species_counts = (
    usable.groupby(["genome.species", "site_class"])
    .size()
    .unstack(fill_value=0)
)

for c in classes:
    if c not in species_counts.columns:
        species_counts[c] = 0

species_counts = species_counts[classes].copy()


# ============================================================
# PROJECT-LEVEL STRUCTURE
# ============================================================

project_usable = usable[
    usable["genome.bioproject_accession"].str.strip().ne("")
].copy()

project_counts = (
    project_usable.groupby(
        ["genome.species", "genome.bioproject_accession", "site_class"]
    )
    .size()
    .unstack(fill_value=0)
)

for c in ["STRICT_AXILLA", "BROAD_AXILLA", "OTHER_SKIN"]:
    if c not in project_counts.columns:
        project_counts[c] = 0

project_counts = (
    project_counts[
        ["STRICT_AXILLA", "BROAD_AXILLA", "OTHER_SKIN"]
    ]
    .reset_index()
    .rename(
        columns={
            "genome.species": "species",
            "genome.bioproject_accession": "bioproject",
        }
    )
)

project_counts["STRICT_SHARED"] = (
    (project_counts["STRICT_AXILLA"] > 0)
    & (project_counts["OTHER_SKIN"] > 0)
)

project_counts["STRICT_SHARED_STRONG"] = (
    (project_counts["STRICT_AXILLA"] >= 3)
    & (project_counts["OTHER_SKIN"] >= 5)
)

project_counts["BROAD_SHARED"] = (
    (project_counts["BROAD_AXILLA"] > 0)
    & (project_counts["OTHER_SKIN"] > 0)
)

project_counts["BROAD_SHARED_STRONG"] = (
    (project_counts["BROAD_AXILLA"] >= 3)
    & (project_counts["OTHER_SKIN"] >= 5)
)

project_counts.to_csv(
    OUT_PROJECTS,
    sep="\t",
    index=False,
    encoding="utf-8-sig",
)


def project_summary_for_species(g):
    return pd.Series(
        {
            "STRICT_AXILLA_PROJECTS": int((g["STRICT_AXILLA"] > 0).sum()),
            "OTHER_SKIN_PROJECTS": int((g["OTHER_SKIN"] > 0).sum()),
            "STRICT_SHARED_PROJECTS": int(g["STRICT_SHARED"].sum()),
            "STRICT_STRONG_SHARED_PROJECTS": int(
                g["STRICT_SHARED_STRONG"].sum()
            ),
            "BROAD_AXILLA_PROJECTS": int((g["BROAD_AXILLA"] > 0).sum()),
            "BROAD_SHARED_PROJECTS": int(g["BROAD_SHARED"].sum()),
            "BROAD_STRONG_SHARED_PROJECTS": int(
                g["BROAD_SHARED_STRONG"].sum()
            ),
        }
    )


if project_counts.empty:
    project_summary = pd.DataFrame()
else:
    project_summary = (
        project_counts.groupby("species")
        .apply(project_summary_for_species, include_groups=False)
    )


# ============================================================
# SPECIES ELIGIBILITY
# ============================================================

species_counts["PRIMARY_COUNT_ELIGIBLE"] = (
    (species_counts["STRICT_AXILLA"] >= PRIMARY_MIN_STRICT_AXILLA)
    & (species_counts["OTHER_SKIN"] >= PRIMARY_MIN_OTHER_SKIN)
)

species_counts["EXPLORATORY_COUNT_ELIGIBLE"] = (
    (species_counts["STRICT_AXILLA"] >= EXPLORATORY_MIN_STRICT_AXILLA)
    & (species_counts["OTHER_SKIN"] >= EXPLORATORY_MIN_OTHER_SKIN)
)

species_counts = (
    species_counts.reset_index()
    .rename(columns={"genome.species": "species"})
)

if not project_summary.empty:
    species_counts = species_counts.merge(
        project_summary,
        left_on="species",
        right_index=True,
        how="left",
    )

project_cols = [
    "STRICT_AXILLA_PROJECTS",
    "OTHER_SKIN_PROJECTS",
    "STRICT_SHARED_PROJECTS",
    "STRICT_STRONG_SHARED_PROJECTS",
    "BROAD_AXILLA_PROJECTS",
    "BROAD_SHARED_PROJECTS",
    "BROAD_STRONG_SHARED_PROJECTS",
]

for c in project_cols:
    if c not in species_counts:
        species_counts[c] = 0
    species_counts[c] = species_counts[c].fillna(0).astype(int)


def design_flag(row):
    if not row["EXPLORATORY_COUNT_ELIGIBLE"]:
        return "NOT_COUNT_ELIGIBLE"

    if row["STRICT_STRONG_SHARED_PROJECTS"] >= 1:
        return "STRONG_SAME_PROJECT_CONTRAST"

    if row["STRICT_SHARED_PROJECTS"] >= 1:
        return "SAME_PROJECT_CONTRAST_AVAILABLE"

    if (
        row["STRICT_AXILLA_PROJECTS"] >= 2
        and row["OTHER_SKIN_PROJECTS"] >= 2
    ):
        return "MULTIPROJECT_BUT_NOT_MATCHED"

    return "PROJECT_SEPARATED_OR_SINGLE_PROJECT"


species_counts["PROJECT_DESIGN_FLAG"] = species_counts.apply(
    design_flag, axis=1
)

species_counts = species_counts.sort_values(
    [
        "PRIMARY_COUNT_ELIGIBLE",
        "EXPLORATORY_COUNT_ELIGIBLE",
        "STRICT_AXILLA",
        "OTHER_SKIN",
    ],
    ascending=[False, False, False, False],
)

species_counts.to_csv(
    OUT_SPECIES,
    sep="\t",
    index=False,
    encoding="utf-8-sig",
)


# ============================================================
# AUDIT RECORDS FOR SPECIES THAT COULD ENTER THE STUDY
# ============================================================

eligible_species = set(
    species_counts.loc[
        species_counts["EXPLORATORY_COUNT_ELIGIBLE"],
        "species",
    ]
)

eligible_audit = df[
    df["genome.species"].isin(eligible_species)
    & df["qc_pass"]
    & df["site_class"].isin(
        [
            "STRICT_AXILLA",
            "BROAD_AXILLA",
            "OTHER_SKIN",
            "OTHER_SKIN_POSSIBLE",
            "SKIN_SITE_UNCONFIRMED",
            "SKIN_UNSPECIFIED",
            "OTHER_APOCRINE_SITE",
            "EXCLUDE_CLINICAL",
            "EXCLUDE_NON_SURFACE",
            "AMBIGUOUS",
        ]
    )
].copy()

eligible_audit.to_csv(
    OUT_ELIGIBLE_AUDIT,
    sep="\t",
    index=False,
    encoding="utf-8-sig",
)


# ============================================================
# CONSOLE SUMMARY
# ============================================================

print("\n================================================")
print("ALL-HUMAN-BACTERIA SPECIES-SELECTION SCREEN V4")
print("================================================")

print("\nAUTO-USABLE QC-passed named-species genomes by class:")
print(usable["site_class"].value_counts().to_string())

print(
    "\nSpecies with >=1 auto-usable STRICT_AXILLA genome:",
    int((species_counts["STRICT_AXILLA"] >= 1).sum()),
)

print(
    f"PRIMARY count-eligible "
    f"(>= {PRIMARY_MIN_STRICT_AXILLA} STRICT_AXILLA and "
    f">= {PRIMARY_MIN_OTHER_SKIN} explicit OTHER_SKIN):",
    int(species_counts["PRIMARY_COUNT_ELIGIBLE"].sum()),
)

print(
    f"EXPLORATORY count-eligible "
    f"(>= {EXPLORATORY_MIN_STRICT_AXILLA} STRICT_AXILLA and "
    f">= {EXPLORATORY_MIN_OTHER_SKIN} explicit OTHER_SKIN):",
    int(species_counts["EXPLORATORY_COUNT_ELIGIBLE"].sum()),
)

show_cols = [
    "species",
    "STRICT_AXILLA",
    "STRICT_AXILLA_PROJECTS",
    "OTHER_SKIN",
    "OTHER_SKIN_PROJECTS",
    "STRICT_SHARED_PROJECTS",
    "STRICT_STRONG_SHARED_PROJECTS",
    "BROAD_AXILLA",
    "BROAD_SHARED_PROJECTS",
    "PRIMARY_COUNT_ELIGIBLE",
    "EXPLORATORY_COUNT_ELIGIBLE",
    "PROJECT_DESIGN_FLAG",
]

print("\nTop species by STRICT_AXILLA count:")
print(
    species_counts.loc[
        species_counts["STRICT_AXILLA"] > 0,
        show_cols,
    ]
    .head(50)
    .to_string(index=False)
)

print("\nSaved:")
for f in [
    OUT_ALL,
    OUT_SPECIES,
    OUT_PROJECTS,
    OUT_AXILLA_REVIEW,
    OUT_AXILLA_LABELS,
    OUT_ELIGIBLE_AUDIT,
]:
    print(" ", f)

print(
    "\nIMPORTANT:"
    "\n1) Only high-confidence explicit non-axillary skin records count as OTHER_SKIN controls."
    "\n2) Groin/inguinal records are NOT counted as ordinary skin controls."
    "\n3) Species selection is based on all human bacterial genomes, not "
    "preselected genera."
    "\n4) Count eligibility and BioProject design quality are reported separately."
    "\n5) Before final inclusion, manually audit the eligible-species records and "
    "axilla labels. Free-text metadata cannot be made 100% error-proof by regex."
)

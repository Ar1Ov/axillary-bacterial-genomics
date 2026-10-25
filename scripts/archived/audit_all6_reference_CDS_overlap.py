import itertools
import pandas as pd

REFS = ["DtpT", "PepV", "PatB", "AgaA", "TpdA", "AecD"]

INPUT = "Supplementary_Malodor_best_hits_annotated_186.tsv"
OUT_LONG = "Supplementary_Malodor_reference_pair_CDS_overlap.tsv"
OUT_MATRIX = "Malodor_reference_CDS_overlap_matrix.tsv"

d = pd.read_csv(INPUT, sep="\t", dtype=str).fillna("")

required = {"genome_id", "reference", "best_hit_bvbrc_feature_id"}
missing = required - set(d.columns)
if missing:
    raise ValueError(f"Missing required columns: {sorted(missing)}")

if d.duplicated(["genome_id", "reference"]).any():
    raise ValueError("Duplicate genome/reference rows detected.")

w = (
    d.pivot(
        index="genome_id",
        columns="reference",
        values="best_hit_bvbrc_feature_id"
    )
    .reindex(columns=REFS)
    .fillna("")
)

rows = []

same_matrix = pd.DataFrame(index=REFS, columns=REFS, dtype=object)
denom_matrix = pd.DataFrame(index=REFS, columns=REFS, dtype=object)
display_matrix = pd.DataFrame(index=REFS, columns=REFS, dtype=object)

for ref in REFS:
    detected = int((w[ref] != "").sum())
    same_matrix.loc[ref, ref] = detected
    denom_matrix.loc[ref, ref] = detected
    display_matrix.loc[ref, ref] = f"{detected}/{detected} (100.0%)"

for a, b in itertools.combinations(REFS, 2):
    both = (w[a] != "") & (w[b] != "")
    same = both & (w[a] == w[b])

    both_n = int(both.sum())
    same_n = int(same.sum())
    pct = 100.0 * same_n / both_n if both_n else 0.0

    rows.append({
        "reference_1": a,
        "reference_2": b,
        "both_detected": both_n,
        "same_exact_CDS": same_n,
        "pct_same_CDS_among_both": round(pct, 1),
    })

    same_matrix.loc[a, b] = same_n
    same_matrix.loc[b, a] = same_n

    denom_matrix.loc[a, b] = both_n
    denom_matrix.loc[b, a] = both_n

    text = f"{same_n}/{both_n} ({pct:.1f}%)"
    display_matrix.loc[a, b] = text
    display_matrix.loc[b, a] = text

out = pd.DataFrame(rows)

out.to_csv(
    OUT_LONG,
    sep="\t",
    index=False
)

display_matrix.to_csv(
    OUT_MATRIX,
    sep="\t",
    index=True
)

print("PAIRWISE EXACT-CDS OVERLAP")
print(out.to_string(index=False))

print("\nMATRIX")
print(display_matrix.to_string())

print("\nSaved:")
print(" ", OUT_LONG)
print(" ", OUT_MATRIX)

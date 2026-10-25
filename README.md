# Cross-species comparative genomics of the human axillary niche and malodor-related proteins

Supporting data, analysis scripts, and figures for the study by **Archana Ovitigala**.

This study examines two related questions: whether bacterial protein families are associated with axillary isolation within species, and how similarity to six malodor-related reference proteins varies across bacterial lineages. Genome metadata, annotations, BLAST output, and organismal trees originate from the Bacterial and Viral Bioinformatics Resource Center (BV-BRC).

The repository supports regeneration of Figures 1–7, postprocessing of the archived six-reference BLAST output, and repetition of the four-control sensitivity analysis. It also preserves additional original analysis scripts and their inputs. It does **not** reconstruct every original database download, manual cohort decision, or remote BV-BRC job.

## Quick start

Download and extract the repository, or clone it. Open a terminal in the repository root: the folder containing `README.md`, `requirements.txt`, `data`, and `scripts`. Extract the files before running Python.

Python 3.10 or later is recommended. The package was previously checked with Python 3.12. Create an isolated environment:

**Windows PowerShell**

```powershell
py -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe scripts/generate_figures.py --font-dir C:/Windows/Fonts --output regenerated_figures
```

**macOS or Linux**

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/generate_figures.py --output regenerated_figures
```

Expected output: `Figure1.pdf` through `Figure7.pdf`, together with seven corresponding JPEGs, in `regenerated_figures/`.

The examples below use `python`. On Windows, substitute `.\.venv\Scripts\python.exe` to use the environment without activating it. Run commands from the repository root. The plotting script resolves its input directory relative to the script itself.

Dependencies include NumPy, pandas, SciPy, statsmodels, Matplotlib, Biopython, PyMuPDF, and Pillow. Installation requires internet access; the tested downstream commands use local bundled inputs and do not require a BV-BRC account. Requirements specify minimum versions rather than a locked environment, so identical pixels across platforms or future dependency versions are not guaranteed.

## Repository contents

| Location | Contents |
| --- | --- |
| `data/cohort_manifest.tsv` | Genome IDs and available cohort/project/site information. The same genome can appear in several analysis arms; this is not a list of independent samples. |
| `data/figures/` | The exact tables and Newick trees read by the figure-generation script. |
| `data/archive/` | Selected archived metadata, feature and presence tables, cohorts, BLAST alignments, reference sequences, neighborhood results, and analysis outputs. Larger text files are compressed as `.gz`. |
| `scripts/generate_figures.py` | Generates all seven figures directly from the bundled tables and trees. |
| `scripts/run_analysis.py` | Prepares an isolated workspace and invokes an analysis. |
| `scripts/run_control_sensitivity.py` | Checks the primary result and repeats the analysis with four additional controls. |
| `scripts/archived/` | Original analysis scripts retained for provenance and further inspection. |
| `figures/` | Reference JPEG and PDF figure outputs. |
| `results/` | Checked sensitivity outputs and the strict-cohort control-count audit. |
| `FILE_MANIFEST.tsv` | Relative paths, file sizes, and SHA-256 checksums for the packaged files, excluding the manifest itself. |

## Regenerate the figures

```bash
python scripts/generate_figures.py --output regenerated_figures
```

The script reads numerical inputs; it does not trace or load existing figure images. Some inputs are archived statistical summaries rather than raw features, so rendering a figure does not automatically rerun its underlying association test.

| Figure | Main inputs in `data/figures/` | What is reproduced |
| --- | --- | --- |
| 1 | `workflow_counts.tsv`, `malodor_targeted_genomes.tsv` | Recorded workflow and cohort counts. This does not rerun the initial metadata screen. |
| 2 | `staph_epi_prjna803478_same_study_fdr_axilla_candidates.tsv` and the corresponding same-study/strict cohort tables | Prevalence of the 14 same-study families in the two overlapping comparisons. |
| 3 | `cory_exact_site_stratified_results.tsv`, `cory_exact_site_final_cohort.tsv`, `candidate_mobile_context_summary.tsv` | Candidate prevalence, saved exact-test results, and neighborhood annotation summaries. |
| 4 | `malodor_best_hit_per_genome_6.tsv` | Individual best-hit similarities and species/site-group medians. |
| 5 | `cory_malodor_codon_tree_tree.nwk`, `SupplementaryCory.tsv` | Saved Corynebacterium tree and aligned protein-similarity tracks. |
| 6 | `staph_malodor_codon_tree_tree.nwk`, `SupplementaryStaph.tsv` | Saved Staphylococcus tree and aligned protein-similarity tracks. |
| 7 | `Figure7_pairwise_exact_CDS_overlap_v3.tsv`, `Figure7B_PatB_AecD_shared_species_summary_v4.tsv` | Exact-CDS overlap and shared-feature reference similarities from saved audit summaries. |

Figure 4 displays all detected hits. Medians use hits with at least 70% query coverage, separately by species and plotted site group; AXILLA and AXILLARY_GROIN are pooled for this display. No-hit entries are excluded from medians and are not treated as zero identity. Median lines appear above points; no IQR or error bars are shown. Species summaries in the manuscript use all reported detected hits and therefore need not equal these filtered, site-specific display medians.

Figures 5–6 mark low-coverage hits with hatching and retain the saved tree-tip order. Figure 7 compares exact coding-sequence identifiers: two reference matches can represent the same genomic feature.

### Fonts and export

The supplied figures use Arial. To reproduce the typography, install Arial or supply a directory containing its TTF files:

```bash
python scripts/generate_figures.py --font-dir C:/Windows/Fonts --output regenerated_figures
```

Arial is not distributed with this repository. Without it, the script uses an available fallback font. Exported figures are 6.5 inches wide, with regular 10-point titles below the panels; JPEGs are 400 dpi and 2,600 pixels wide. Figure legends are supplied separately with the manuscript and are not included in this repository.

## Reproduce BLAST postprocessing

```bash
python scripts/run_analysis.py blast --workdir reruns/blast
```

This command decompresses/copies archived inputs into the specified workspace and runs `analyze_targeted_malodor_blast_6.py`. It uses the archived alignment output and the fixed 186-genome panel. It does not submit a new BLAST search.

Seven tables are produced in `reruns/blast/`:

- `malodor_blast_hits_annotated_6.tsv`
- `malodor_best_hit_per_genome_6.tsv`
- `malodor_best_hits_relevant_genus_6.tsv`
- `malodor_homology_group_summary_6.tsv`
- `malodor_homology_group_summary_all_genomes_6.tsv`
- `malodor_top_hits_by_reference_6.tsv`
- `malodor_top_hits_by_reference_all_genomes_6.tsv`

Expected totals are 1,981 archived alignment rows, 186 panel genomes, and 1,116 genome–reference combinations. The best-hit table contains 1,067 detected combinations. The panel comprises 145 Staphylococcus and 41 Corynebacterium genomes.

The six query sequences are supplied in `data/archive/malodor_reference_proteins_6.faa`: DtpT, PepV, PatB, AgaA, TpdA, and AecD. Reference labels describe the query proteins and do not establish the function of every matching sequence.

The recovered BV-BRC job was 23736298, using a 447,943-protein database, an E-value cutoff of 1e-5, a maximum of 1,000 target sequences per query, and no minimum query-coverage requirement. The postprocessor selects one alignment per genome/reference by descending bit score, coverage, and identity, then ascending E-value. The full search database and original BLAST executable version are not supplied; this package reproduces downstream processing of saved search output.

## Reproduce control sensitivity

```bash
python scripts/run_analysis.py sensitivity --workdir reruns/sensitivity
```

The script first checks the primary 5-axilla/117-other-skin result against the saved results. It then restores controls 1282.4413–1282.4416, yielding 5 axillary and 121 other-skin genomes.

| Expected result | Primary 5/117 | Expanded 5/121 |
| --- | ---: | ---: |
| Observed families | 4,410 | 4,465 |
| Tested families | 1,485 | 1,528 |
| FDR-significant axilla-enriched families | 0 | 0 |
| FDR-significant axilla-depleted families | 5 | 2 |

Outputs are written to `reruns/sensitivity/sensitivity_results/`. Inspect `run_status.json`, `baseline_recomputed.tsv`, and `expanded_5_vs_121.tsv`; reference copies are in `results/`. A successful check reports `baseline_matches_saved: true` and `expanded_run_completed: true`.

The notation “2–6/117” describes control carriage of the 14 same-study families in the strict cohort: two families occurred in two controls each, eleven in three controls each, and one in six controls. All fourteen occurred in one of five strict axillary genomes. See `results/Figure2StrictControlCounts.tsv`.

## Additional archived analyses

These entry points are supplied for inspection and rerunning. They were not newly executed during the previous package-validation pass:

```bash
python scripts/run_analysis.py strict --workdir reruns/strict
python scripts/run_analysis.py same-study --workdir reruns/same_study
python scripts/run_analysis.py cory --workdir reruns/cory
python scripts/run_analysis.py haemolyticus --workdir reruns/haemolyticus
python scripts/run_analysis.py ecoli --workdir reruns/ecoli
```

Each analysis needs a new or empty work directory. The runner intentionally refuses to overwrite a populated directory. It copies all archived inputs into each workspace, so allow additional disk space. Files already copied into a workspace are not evidence that they were newly generated; use the script's completion status and reported outputs to assess the run.

Changing an analysis input does not automatically update `data/figures/`. If performing a new analysis, inspect the new results and explicitly update the corresponding figure input tables before rendering; the packaged figure inputs preserve the manuscript's recorded analysis.

## Validation and interpretation

In the previous package check, all seven figures were generated from a clean copy containing only the plotting script and figure inputs. All seven BLAST postprocessing tables matched the prior verified outputs numerically, and the sensitivity check reproduced the counts above. This README revision does not represent a new rerun of every archived analysis.

The strict and same-study S. epidermidis cohorts overlap: all 22 same-study controls are included in the strict comparison. These are sensitivity comparisons, not independent replication. Body site and project/population structure remain confounded. Corynebacterium candidates did not survive FDR correction.

Protein similarity is descriptive evidence of sequence variation, not proof of odor production. PatB and AecD selected the same feature in 184 of 185 co-detected genomes; reference matches should not be counted automatically as separate genes. Organismal trees and reference similarities do not establish independent evolutionary origins, convergence, or horizontal transfer.

Initial metadata collection/filtering, the remote BLAST search, and codon-tree jobs are not reconstructed end to end. Workflow counts are archived summaries, and the original reason for omitting the four controls was not recorded. Historical filenames or script comments using “replicated” should be interpreted in light of the overlapping cohorts described above.

## Troubleshooting

- **Missing Python module:** install `requirements.txt` using the same Python executable that runs the script.
- **File not found:** extract the complete repository and preserve its folder structure. Do not move a script out of `scripts/`.
- **Work directory not empty:** choose a new `--workdir` rather than overwriting an existing rerun.
- **Different typography:** supply Arial with `--font-dir`; otherwise font substitution is expected.
- **Baseline check fails:** retain the error and inputs, confirm that files are unchanged, and report the Python/package versions. Do not disable the assertion to obtain an apparently successful run.

## Citation and reuse

When using this work, credit Archana Ovitigala and cite the accompanying manuscript and the exact repository commit or release used. The README does not imply that the manuscript has been accepted or published. Retain BV-BRC provenance and cite the original data and methods as appropriate.

No code license is included in this package. Consult the author about reuse permissions until a license is supplied; third-party source data retain their applicable terms. Report reproducibility problems through the repository's issue tracker, including the command, error output, operating system, and Python/package versions.

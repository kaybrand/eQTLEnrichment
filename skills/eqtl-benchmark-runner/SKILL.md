---
name: eqtl-benchmark-runner
description: Run, extend, and interpret the eQTL benchmark pipeline (enrichment/recall of regulatory-element predictions against eQTL Catalogue signal) in the IGVF eQTL-benchmark project. Use when the user wants to sync IGVF E2G predictions, build/validate predictions and biosample-tissue-pair configs, run the default scE2G/distanceToTSS benchmark, add a new model or pair set to a benchmark run, submit the benchmark via sbatch, or interpret enrichment/recall outputs. Synthesized from 6 past sessions across the eQTLEnrichment, CATlas, and IGVF benchmark projects (2026-10-15) — several sections are flagged 🔲 TODO / 📌 VERIFY IN IGVF and should be checked against the current IGVF branch before being treated as ground truth.
---

# eQTL Benchmark Runner (IGVF)

> Synthesized from sessions: b73cdace, 10cdeda8, 7d1019e3, 6a5dc004, b4cc6b41, 1336b414
> Source projects: eQTLEnrichment, CATlas, IGVF — generated 2026-10-15

## Overview

The eQTL benchmark measures how well a regulatory-element-to-gene prediction model's scored candidate elements overlap real eQTL signal in a matched tissue, relative to background SNPs — it is a way to evaluate regulatory element prediction models (scE2G, distance-to-TSS, and others) against an orthogonal ground truth (eQTL Catalogue variants) rather than against each other's own assumptions. This SKILL covers the full path: syncing IGVF E2G predictions, building and validating the predictions/config tables, running the default four-method benchmark with no changes, extending it with new methods or pairs, submitting via `sbatch`, and reading the resulting enrichment/recall outputs. The default comparison is **scE2G Multiome, scE2G scATAC, scABC, and Distance to TSS** (the baseline). **Fairness constraint, stated plainly: every method in a single benchmark run must be compared on the exact same set of biosample-tissue pairs — mixing pair sets across methods in one run invalidates the comparison.** The only exception is an explicit user request to do otherwise, and even then this SKILL must spell out exactly what becomes non-comparable.

---

## What the Benchmark Measures

### Enrichment

Enrichment asks: *are real eQTL variants more likely to fall inside this model's predicted (scored) candidate elements than background SNPs are?* Concretely, for a given method, tissue, and distance bin:

- `x1` = number of eQTL variants (in the matched tissue) that overlap a predicted element passing the method's score threshold
- `n1` = total number of eQTL variants in that tissue (within the distance bin / gene universe)
- `x2` = number of background SNPs (1000G) overlapping a predicted element passing the same threshold
- `n2` = total number of background SNPs considered

```
enrichment = (x1 / n1) / (x2 / n2)
SE(log enrichment) = sqrt( (n1-x1)/(x1*n1) + (n2-x2)/(x2*n2) )
CI_low  = exp(log(enrichment) - z * SE_log_enr)
CI_high = exp(log(enrichment) + z * SE_log_enr)
```

⚠️ **CONFLICT / fixed-bug history** — three distinct bugs were found and fixed across sessions, and the fixes are consistent with each other:
1. The SE formula was originally computed as `sqrt(a) + b` instead of `sqrt(a + b)` (background term left outside the square root), which always *understated* the true SE. Fixed in `counts_to_enrichment_distance.R` and `counts_to_enrichment_across_thresholds.R` (sessions 10cdeda8, b73cdace, 1336b414 all agree).
2. The numerator's "0–30000Kb" (all-distances) bin skipped the score threshold filter entirely (count_matrix_distance.R). Fixed by applying the same `score>=score.thresh` filter to that branch (commit `133acc0`, session 7d1019e3).
3. The denominator (background overlap count) was **never** thresholded by score, in any distance bin, for any method — so four methods with four different thresholds all reported the identical background-overlap count. Fixed by adding `dplyr::filter(commonVarPredInt, score>=score.thresh)` before counting (commit `ab9d27c`, session 7d1019e3).
   - Per 10cdeda8's `timestamp_signal`, commits `ab9d27c`/`133acc0`/`2af8c86` (bugs 2 and 3, plus the heatmap thresholded/candidateUniverse split) were **present on the IGVF branch already** (inherited via cherry-pick from the branch that originated them) but **not yet on `catlas-benchmark-2026-09`** as of that session. 📌 VERIFY IN IGVF that your checkout includes these three commits before trusting absolute enrichment numbers.

A high enrichment score (CI excluding 1) means the model's scored elements are genuinely concentrating real eQTL signal relative to background; a score near or below 1 means the model's threshold is not isolating true regulatory signal better than chance.

**Pooled vs. per-pair:** an individual (biosample, tissue, distance bin) row is one enrichment estimate. A pooled **"AllMatches"** row sums `x1,n1,x2,n2` across every biosample matched to that tissue *before* computing a single enrichment ratio — it is not an average of per-pair enrichments. 🔲 TODO: NEEDS VALIDATION / known limitation (documented in README, commit `23bf8a9`, not fixed) — pooling treats repeated/related eVariant-biosample pairs as statistically independent, which they are not; a bootstrap-based correction was drafted and explicitly parked ("backburner") by the user, not implemented.

### Recall

Recall asks: *of the real eQTL-linked (eGene, eVariant) pairs in a tissue, how many does the model's predicted set actually capture?* It is computed per method, per distance bin or per threshold, via `compute_prediction_table_distance.R` (by-distance) and the across-thresholds equivalent — both consistently threshold the score **before** branching on distance bin (session 7d1019e3 confirms this script was "never buggy"). Recall tables exist in the same thresholded-vs-**candidate-universe** (unthresholded) split as enrichment tables.

Universe grouping matters here because the **candidate-universe recall table is computed once per universe group, not once per method** — methods sharing an identical candidate-element file (e.g. scABC/distanceToTSS mirroring scE2G_multiome's file) read the *representative* member's table rather than recomputing their own. See **Universe Grouping** below.

### Reading Enrichment and Recall Together

Enrichment and recall are reported together because they trade off against the score threshold: a stricter threshold raises enrichment (fewer, higher-confidence calls) but lowers recall (fewer true links captured), and vice versa. The benchmark's **enrichment-recall curve** (`plot_enrichment_recall_curve`, one PDF+TSV+`.plotdata.rds` per tissue, plus an interactive per-tissue HTML on the IGVF branch) sweeps the threshold and plots enrichment against recall so models can be compared across the whole operating range, not just at one cutoff. A second output, **enrichment-at-fixed-recall** (`plot_enrichment_with_ci`, `enrichmentAtRecall/`), holds recall constant (e.g. recall=70%) across methods and compares enrichment at that matched recall, with a Bonferroni-corrected pairwise significance table (`pairwiseComparisons.Recall<r>.GTExTissue<tissue>.tsv`) — this is generally the more interpretable "which model wins" comparison since it controls for the recall/enrichment tradeoff. Models are not expected to be jointly optimal at every point; read the curve, not a single number.

---

## Environment Setup

Two separate conda environments are used, split by role (agreed across all 6 sessions):

| Env | Purpose |
|---|---|
| `eQTLEnv` | R/Python analysis, prediction-table generation, validation, all per-rule Snakemake `conda:` directives (`workflow/envs/eQTLEnv.yml`) |
| `run_snakemake9` | **Only** this env has `snakemake-executor-plugin-slurm` installed — it is required to drive the Snakemake **controller** process itself with `executor: slurm`. `eQTLEnv`'s own bundled snakemake (9.3.3) has no executor plugin and cannot submit Slurm jobs. |

Activation: do **not** use `conda activate` in batch scripts — call the interpreter by absolute path instead. Per `validate_predictions.sbatch`'s own comment: *"Absolute interpreter path rather than conda activate: auto-activating conda in a batch context is a documented source of held jobs on this cluster."*

```bash
source /oak/stanford/groups/engreitz/Users/kaybrand/software/miniforge3/etc/profile.d/conda.sh && conda activate eQTLEnv   # interactive use only
# OR, in scripts, call directly:
/oak/stanford/groups/engreitz/Users/kaybrand/software/miniforge3/envs/eQTLEnv/bin/Rscript ...
/oak/stanford/groups/engreitz/Users/kaybrand/software/miniforge3/envs/eQTLEnv/bin/python ...
/oak/stanford/groups/engreitz/Users/kaybrand/software/miniforge3/envs/run_snakemake9/bin/snakemake ...
```
📌 VERIFY IN IGVF: the exact `/oak/.../miniforge3/envs/...` paths are this user's Sherlock setup; confirm your own env locations.

Key dependencies: `snakemake 9.16.3` + `snakemake-executor-plugin-slurm 2.3.1` (`run_snakemake9`); R packages `dplyr, ggplot2, data.table, tidyr, egg, gridExtra, colorspace, plotly, htmlwidgets, htmltools, rmarkdown, optparse` (`eQTLEnv`); `bedtools`; `csvtk`; Python `click`, `synapseclient` (lazy-imported — table generation needs no Synapse credentials unless you use `--sync`), `pandas`.

---

## Directory Structure

```
<project_root>/                          # e.g. the eQTL_benchmark_IGVF worktree 📌 VERIFY IN IGVF
├── workflow/
│   ├── Snakefile                        # DAG entry; calls check_config_inputs(); computes THRESHOLDED_METHODS, UNIVERSE_GROUP_*
│   ├── rules/*.smk                      # preprocessing, recall, enrichment, visualization, utils, config_checks
│   ├── scripts/{preprocessing,enrichment_recall,visualization}/*.R, *.py
│   ├── IGVF_prediction_validation/      # make_prediction_config.py, main_validation_script.py, validators/, rescue_scripts/
│   ├── scripts/run_benchmark.sbatch     # controller submission script
│   ├── scripts/validate_predictions.sbatch
│   └── envs/eQTLEnv.yml
├── config/
│   ├── config_<name>.yaml               # run config: outDir, methodsTable, predictionsTable, methods list, etc.
│   ├── config_methods*.tsv              # model/methods table (tracked)
│   ├── predictions.tsv                  # example/committed predictions table
│   ├── eQTL_tissue_matches.tsv          # biosample-tissue pairs (tracked)
│   └── local/                           # GITIGNORED — derived predictions_*.tsv, synapse_manifests/, tissue-pairing worksheets
├── results/
│   └── <outDir>/                        # GITIGNORED — see Outputs section
└── .baton/                              # GITIGNORED — session handoff notes, <YYYY-MM-DD>-<slug>.md
```

Raw prediction files live outside the repo, synced from Synapse to a scratch/download directory, e.g. `$SCRATCH/<download_dir>/<dataset>/<biosample>/<biosample>_<Model>[.e2g].tsv.gz` — 📌 VERIFY IN IGVF: there is no default path anymore (the historical hardcoded default was deliberately removed, commit `8eadc42`); `--predictions-dir` is a required flag you must supply. `/tmp` (`$TMPDIR`) on Sherlock compute nodes is wiped between sessions — durable test configs belong under `results/test_configs/` (gitignored, on Oak), not `/tmp`.

---

## Configuration

### Predictions Table

**Format:** TSV. **Path pattern:** `config/predictions.tsv` (committed example) or `config/local/predictions_<name>.tsv` (derived/gitignored, regenerated by tooling — never hand-edit the generated table).

Required / common columns (agreed across all 6 sessions):

| Column | Description |
|---|---|
| `biosample` | Unique biosample ID. No leading/trailing whitespace; duplicates forbidden. Must be **dataset-prefixed** (`{dataset}_{biosample}`) — raw biosample folder names are **not** globally unique across IGVF datasets (e.g. `b_activated_d1` appears under both `igvf9` and `igvf12`); a bare name silently overwrites across datasets. |
| `GTExTissue` | Comma-separated list of eQTL tissues this biosample is paired to, e.g. `dendritic,myeloid,blood`. |
| `<method>` (one column per method) | Must match the method's `method` value in the methods table **exactly** (case-sensitive) — a mismatch produces a silent all-NaN column and the method silently disappears from the benchmark with no error. Cell value = path to that method's prediction file for this biosample, or blank/NaN if none exists. |
| `SampleTermName`, `CellAnnotation`, `SampleTermID` (optional) | Used as display-name / metadata fallbacks in plotting; not required for the benchmark to run. |

Example row shape (illustrative, paths abbreviated):
```
biosample            GTExTissue          scE2G_multiome                 scE2G_ATAC   scABC(mirror)  distanceToTSS(mirror)
igvf1_hl60           dendritic,myeloid,blood   .../igvf1_hl60_scE2G_multiome_powerlaw_v3.e2g.tsv.gz   .../igvf1_hl60_scE2G_scATAC....gz   (same file, score_col=ABC.Score)   (same file, score_col=distance)
```

### Building a Biosample-Tissue Prediction Config from a Synced Directory

1. **Sync predictions from Synapse** (IGVF-specific helper, optional if you already have a local tree):
   ```
   python workflow/IGVF_prediction_validation/make_prediction_config.py --sync -d $SCRATCH/<download_dir>
   ```
   Requires `synapseclient` + `~/.synapseConfig` credentials (lazy-imported, so this step is skippable if you already have files).

2. **Expected structure after sync:** `<predictions_dir>/<dataset>/<biosample>/<biosample>_<Model>[.e2g].tsv.gz`, one subdirectory per dataset (e.g. `igvf0`..`igvf18`, `GM12878_10XMultiome`), one subdirectory per biosample inside.
   - 🔲 TODO: NEEDS VALIDATION — some historical sync trees had excluded/legacy dataset directories (e.g. `(Archived)Y2_versions`, an "outdated prediction set" dataset); confirm your synced IGVF tree doesn't need an equivalent exclude list.

3. **Identify biosamples/tissues represented:** the tooling walks the synced tree and classifies each file by model using `extract_model_name()`, which anchors the match to the **trailing underscore-delimited token of the filename stem** (longest suffix first) — **not** a naive substring search. (Gotcha, fixed bug: an earlier unanchored substring match classified every file in biosample `b_quiescent`'s folder as `SCENT` because `"scent"` is a substring of `"quiescent"`.)

4. **Generate the predictions table:**
   ```
   python workflow/IGVF_prediction_validation/make_prediction_config.py \
     -d <predictions-dir> \
     -t config/<tissue_pairings>.tsv \
     -o config/local/predictions_<name>.tsv \
     [--matched-only] \
     [--models scE2G_multiome,scE2G_ATAC,scABC,distanceToTSS,...]
   ```
   - `--predictions-dir`/`-d` is **required** (no default).
   - `--matched-only` drops any biosample with no eQTL tissue pairing at all (empty/NaN `GTExTissue`).
   - `--models` restricts which method columns get built.
   - Mirrored/derived columns: some "methods" are not separate physical files but another method's file read with a different `score_col` — handled via a `MIRRORED_COLUMNS` / `SCE2G_DERIVED_BASELINES` mapping in the script (e.g. `scABC`, `distanceToTSS`, `Kendall` all mirror the `scE2G_multiome` file). 📌 VERIFY IN IGVF that the default four methods use this mirroring on your branch — it directly affects the intersection/fairness logic below, since "distanceToTSS has no prediction" may really mean "scE2G_multiome has no prediction for this biosample" if they mirror the same file.

5. **Missing predictions:** a biosample lacking a file for a given method is simply left blank/NaN in that column — this is not an error at this stage. It becomes relevant at run time because the pipeline (`utils.smk`) does `key[['biosample', method, 'GTExTissue']].dropna()` **per method**, silently excluding that biosample from that method's run. Whether blank cells get algebraically **excluded from the whole comparison** (fairness) vs. just from one method's own output depends on which pair-set-building step you run next — see **Default Configuration** and **Adding a New Method** below.

### Biosample-to-Tissue Pairs Config

**Format:** TSV. Two layouts have been seen and are both supported by the loader (`read_gtex_tissue_map()`, detected automatically by column presence):

| Layout | Columns |
|---|---|
| Legacy "wide" | `cluster` (biosample/cluster name), `tissue` (already comma-joined list) |
| 2026 "long" | `Dataset`, `Cluster_2026`, `eQTL_tissue` (**one row per biosample-tissue pair** — same biosample repeats once per matched tissue), `CellAnnotation`, `SampleTermID` |

**Path pattern:** 📌 VERIFY IN IGVF — eQTLEnrichment session examples show `config/local/preprare_2026_eQTL_tisue_matches/2026_eQTL_Tissue_pairs.tsv` (note: `preprare`/`tisue` are literal typos in that project's actual path, not a formatting error here); the IGVF-branch session (10cdeda8) instead lists a tracked `config/eQTL_tissue_matches.tsv`. Confirm which path your IGVF config actually points at.

**Resolution logic:** for the long format, the loader builds `key = f'{Dataset}_{Cluster_2026}'` (falls back to bare `Cluster_2026` if `Dataset` is blank) — this key format must match the dataset-prefixed biosample IDs produced by the prediction-config builder. Each row's tissue is appended into a comma-joined, de-duplicated `GTExTissue` list for that key. Downstream, the pipeline splits `GTExTissue` on commas again and emits one `(biosample, tissue)` pair per entry — so **N tissues for one biosample round-trips to N separate benchmark pairs.**

**Fairness note:** the pairs file itself is global, not per-method — but a given run's *active* pair set is determined by which rows survive the per-method `.dropna()` filtering (missing predictions) **and** by whichever intersection step you run (see Default Configuration / Adding a New Method). The pairs file does not need to change to add/remove a method from a comparison; the derived predictions-table subset does.

**File-resolution-logic example** (illustrative pattern from the eQTLEnrichment/Open4Gene session, `b4cc6b41` — this exact function is project-specific but the *pattern* of "try newest format, fall back, else raise rather than silently skip" is a generally sound approach worth replicating for any method with multiple file-format generations):
```python
def open4gene_file(b):
    v2  = f"{PRED}/{ds}/{cl}/{b}_Open4Gene.tsv.gz"                      # newest format
    re1 = f"{FILTERED}/{ds}/{cl}/{b}_open4gene_standard_zeroz_1Mb.tsv.gz"  # reissue, filtered
    re5 = f"{PRED}/{ds}/{cl}/{b}_open4gene_standard_zeroz.tsv.gz"          # reissue, unfiltered
    if os.path.exists(v2):  return v2
    if os.path.exists(re1): return re1
    if os.path.exists(re5): raise SystemExit(f"{b}: reissue present but not filtered yet")
    return None  # no file for this biosample -> excluded from this method's column
```
📌 VERIFY IN IGVF — this specific function belongs to the eQTLEnrichment/Open4Gene project; the four IGVF default methods do not need this multi-generation fallback as far as any session documents, but the general pattern applies if you add a method whose files have changed format over time.

### Model Config

**Path pattern:** `config/config_methods*.tsv` (e.g. `config_methods_IGVF_Annual_Meeting_2026.tsv` in the eQTLEnrichment examples — 📌 VERIFY IN IGVF for the IGVF-branch filename). One row per method.

**Required columns** (`REQUIRED_METHOD_COLUMNS`, enforced by `config_checks.smk`, agreed across all sessions):

| Field | Controls |
|---|---|
| `method` | Must exactly match (case-sensitive) both the predictions-table column header and the `config['methods']` list entry in the run YAML. Also the key used to look up a validator in `VALIDATOR_MAP`. |
| `boolean` | `TRUE`/`FALSE`. If `TRUE`, `threshold` must equal `1`. |
| `inverse_predictor` | `TRUE`/`FALSE`. If `TRUE`, lower raw score = better prediction; the Snakefile auto-negates the configured `threshold`'s sign (and `process_predictions.R` negates the score) so comparisons are uniformly `score >= threshold` internally. ⚠️ Fixed bug: the negation check must be an explicit string comparison (`value.strip().lower()=='true'`) — a naive Python truthiness check wrongly negated blank cells, because blank cells arrive as the truthy *string* `'None'` after `fillna('None')`. |
| `threshold` | Numeric, in the method's **own raw score units**. May be **blank** — a method with no threshold is excluded from thresholded outputs (heatmaps, distance comparisons, report tables) but still included in **all** threshold-free outputs (enrichment-recall curves, candidate-universe heatmaps/tables). Negative thresholds are allowed for inverse predictors. |
| `score_col` | Must be a literal column name **present in the actual prediction file's header** for every biosample using this method. |
| `geneUniverse` | ⚠️ **CONFIRMED DEAD CODE** per session 6a5dc004 (and independently re-confirmed from 2025 debug notes) — no rule actually reads this column; the real, global gene universe is `config['TSS_reference']` in the top-level run YAML. The column must still be **non-blank and point at a file that exists on disk** to pass `config_checks.smk`, even though its value is never used downstream. Do not expect a per-method universe override to work via this column. |

**Optional columns:**

| Field | Controls |
|---|---|
| `nonneg_only` | `TRUE`/`FALSE`, default `FALSE` if blank/absent. If `TRUE`, rows with raw score `< 0` are **dropped** (treated as "not predicted") *before* any score inversion — confirmed with submitters (Signac, Cicero) that a hard drop, not `abs()`, is the correct semantics for a signed correlation/z-statistic. |
| `pred_name_long` | Display label only; no pipeline effect. |
| `original_score_col` | Used only by certain validators' `rescue()` logic as the source column to recompute from (e.g. pgBoost, SCARlink, GEEES, Abs-correlation validators). |
| `color` | Hex string matching `#[0-9A-Fa-f]{6}`; validated if present. |
| `universe_group` | 🔲 TODO: NEEDS VALIDATION — this column only appears in session `1336b414` (eQTLEnrichment, added specifically to force Open4Gene's 11 scoring arms into one shared candidate universe). It is a manual override: methods given the same non-blank `universe_group` value are force-joined into one universe group, and declaring conflicting `universe_group` values for methods that already share a universe raises a `ValueError`. **Confirm whether this column exists on the IGVF branch** before using it — it may not be needed for the four default methods, which already group via file-path identity / `ALWAYS_SHARED_UNIVERSE` (see **Universe Grouping**). |

**Connections:** the `method` value is the single join key across three places — the predictions-table column header, the `config['methods']` list in the run YAML, and the `VALIDATOR_MAP` key in `main_validation_script.py`. Get any one of these three out of sync with the other two and the method silently has zero predictions (all-NaN column) or fails fuzzy-matched validation with a "did you mean" suggestion.

**Example config entry** (placeholder method `ModelX`, TSV row shown with explicit field labels for readability):
```
method           = ModelX
boolean          = FALSE
inverse_predictor= FALSE
nonneg_only      = FALSE
pred_name_long   = Model X
threshold        = 0.5
score_col        = Score
original_score_col =            # blank unless a validator's rescue() needs a source column
color            = #1F77B4
geneUniverse     = resources/genome_annotation/CollapsedGeneBounds.hg38.bed   # must exist on disk; value itself unused
universe_group   =              # blank unless ModelX must be force-grouped with another method (🔲 verify column exists on IGVF branch)
```

### Validator Package

#### Selecting an Existing Validator

Validators live in `workflow/IGVF_prediction_validation/validators/` and are dispatched by a `VALIDATOR_MAP` dict in `main_validation_script.py`, keyed on the exact `method` name. A method with no `VALIDATOR_MAP` entry falls back to `DefaultValidator` (float-parseable check only).

All prediction files, regardless of validator, must first pass a **core format check** (`_validate_core_format()`):
- Columns `ElementChr, ElementStart, ElementEnd, GeneSymbol` present (flexible order).
- `ElementChr` matches `^chr([0-9]+|[XYM])$`.
- `ElementStart`/`ElementEnd` are positive integers (no floats, no negatives).
- `GeneSymbol` matches `^[A-Za-z0-9\-]+$`.
- `GTExTissue` is explicitly exempted from path/score validation (treated as a plain string).

Existing validator classes and what they check:

| Validator | Model type | Checks |
|---|---|---|
| `BaseValidator` (abstract) | — | Defines `is_score_valid()` / `rescue()` / `pre_validate_and_reformat()` extension points. |
| `DefaultValidator` | Any method with no specific validator | Score column is float-parseable. |
| `pgBoostValidator` | pgBoost | `Percentile` column exists, numeric, in `[0,1]`. `Score==1` saturated rows must be excluded — **unconditionally**, in `pre_validate_and_reformat()`, not in `is_score_valid()` (see gotcha below). |
| `AbsCorrelationValidator` | Absolute-correlation-style scores | Target `abs_Score` column in `[0,1]`. |
| `DivergentAbsValidator` (subclasses `AbsCorrelationValidator`) | Signed z-statistics (e.g. `Open4Gene_abs`) | Overrides the bound check to numeric & `>=0`, **no upper bound** (a `|z|` is not a bounded correlation). |
| `SCARlinkValidator` | SCARlink | Boolean dtype check. |
| `GEEESValidator` / `scE2GValidator` | GEEES / scE2G family | Numeric float check. |

**Gotcha (pgBoost, high confidence):** `is_score_valid()` by default only reads the first ~100,000-row chunk. A check that is *conditional* on that chunk (e.g. "exclude `Score==1` rows if any are seen") can pass validation while saturated rows later in a multi-million-row file are silently never rescued. **Lesson: any exclusion that must hold for every row has to run unconditionally in `pre_validate_and_reformat()`**, which the orchestrator always runs regardless of chunking, not as a conditional inside `is_score_valid()`.

#### Building a New Validator

Build a new validator when: the model is a new type not covered above, has a non-standard score column, needs a percentile or absolute-value column checked, or needs custom numeric bounds.

Structure:
1. Subclass `BaseValidator` directly, or an existing validator if your model's semantics are close (e.g. subclass `AbsCorrelationValidator` the way `DivergentAbsValidator` does, just to relax/change the bound).
2. Implement `is_score_valid()` — the numeric/range/dtype check specific to this model's score semantics. Examples seen in source conversations:
   - Numeric, bounded `[0,1]`: `pgBoostValidator`/`AbsCorrelationValidator` pattern.
   - Presence of a `Percentile` column: `pgBoostValidator`.
   - Presence of an absolute-value column, unbounded above: `DivergentAbsValidator` (`numeric and >= 0`).
   - Boolean dtype: `SCARlinkValidator`.
3. Implement/override `rescue()` only if the model's files need format correction (e.g. recomputing a percentile from a raw score, dropping sentinel rows). `rescue()` writes a corrected file; if your rescue needs a source column, wire it through `original_score_col` in the methods table.
4. Put any **must-always-apply** exclusion rule in `pre_validate_and_reformat()`, not in `is_score_valid()`, if there's any chance of chunked/partial reads missing it later in the file (see the pgBoost gotcha above).
5. **Register** the new class in `VALIDATOR_MAP` in `main_validation_script.py`, keyed on the method's exact name as used in the methods table / predictions table.

#### Running Validation

```
python workflow/IGVF_prediction_validation/main_validation_script.py \
  -p <predictions_table.tsv> -m <config_methods.tsv> [--check-all-rows] [-v|--verbose]
```
- This is a **manual pre-step** — it is **not** wired into the Snakemake DAG and must be run and its output reviewed by a human before launching the pipeline.
- Writes `<table>_validated.tsv` next to the input. Cells that failed validation (and could not be rescued) are set to `NaN`.
- **Passing** looks like a per-model summary reporting "complete" or a small, expected "missing" list.
- **Failure interpretation:** if a model shows **every** biosample missing, that is almost always a method-name / `VALIDATOR_MAP` / `score_col` mismatch, not a real data problem — check the three-way join described in **Model Config** above first.
- An empty/NaN cell in the validated table **silently drops that biosample+method pair from the benchmark with no runtime error** — this is by design (manual-step convention) but means you must actually read the validation summary, not just check the exit code.

### Config Validation

This is a **separate, pipeline-level** step from the Validator package above — the Validator package checks individual prediction *files*; config validation checks that the *config/methods table/predictions table together* are internally consistent, and runs automatically as part of every dry run or real run (no separate manual invocation needed, though you can trigger it alone via a dry run).

Implemented in `workflow/rules/config_checks.smk`, function `check_config_inputs(methods_config, config)`, called directly from `workflow/Snakefile` **before** the DAG is built. It raises `snakemake.exceptions.WorkflowError` with a **numbered list of every problem found**, not just the first.

Rules enforced (most complete version documented in session 10cdeda8 / IGVF branch):
- All of `REQUIRED_METHOD_COLUMNS` present in the methods table.
- Methods-table index (`method`) has no duplicates.
- Every entry in `config['methods']` has a matching row in the methods table (fuzzy "did you mean" suggestion via `difflib` on mismatch).
- `threshold` may be blank, else must be float-parseable; negative values are allowed for `inverse_predictor=TRUE`.
- Boolean columns (`boolean`, `inverse_predictor`, `nonneg_only`) must be one of `{TRUE, FALSE, True, False}` — `'yes'`, `'1'`, `'T'` etc. are flagged as errors here precisely because the R side would otherwise silently read them as `FALSE`.
- If `boolean==TRUE`, `threshold` must equal `1`.
- `score_col` non-blank.
- `geneUniverse` file must exist on disk (even though its contents are unused downstream — see Model Config).
- `color`, if present, matches `#[0-9A-Fa-f]{6}`.
- The predictions table (`config['predictionsTable']`) must exist and have `biosample` and `GTExTissue` columns.
- Biosample names: no leading/trailing whitespace, no duplicates.
- Every method in `config['methods']` must have a column in the predictions table matching its name **exactly** (fuzzy suggestion on mismatch).
- Each method needs at least one biosample with both a non-null prediction file **and** a non-blank `GTExTissue`, or nothing would ever be benchmarked for it.
- Each prediction file's header must contain `ElementChr, ElementStart, ElementEnd, GeneSymbol` plus that method's `score_col`.

A separate **runtime** check (not at DAG-build time, inside the `process_predictions` rule): after other filters, it checks that **at least one row** passes the configured threshold (`score_col >= threshold`, or `<= -threshold` for inverse predictors) and fails loudly with the observed min/max score if not. Skipped entirely for methods with no threshold.

**Motivating incident** (2025-10-19, referenced across sessions): a mis-specified score column or an invalid/empty prediction file previously produced a silently-empty intermediate, which only surfaced as a cryptic empty-dataframe crash deep in the DAG after **~1 hour of compute**. Both the startup checks and the runtime zero-links check exist specifically to catch this class of error before any real compute happens.

**How to run manually:** just do a dry run — `check_config_inputs()` fires automatically:
```
snakemake -n -s workflow/Snakefile --configfile <cfg> --quiet rules
```

### Universe Grouping

No value is stored in a table column by default (see `universe_group` 🔲 TODO note in Model Config for an optional override). Grouping is **computed at Snakefile-parse time** by `determine_candidate_universe_groups()` in `workflow/rules/utils.smk`, producing three Python dicts: `UNIVERSE_GROUP_LABEL` (method → label), `UNIVERSE_GROUP_MEMBERS` (label → [methods]), `UNIVERSE_GROUP_REPRESENTATIVE` (label → the one member whose tables are actually read for the whole group).

⚠️ **CONFLICT — two generations of this logic exist; use the newer one.**
- **Deprecated** (seen in session 7d1019e3, pre-fix eQTLEnrichment): a 4-tier fingerprinting strategy — same file path (free) → hardcoded `KNOWN_SHARED_UNIVERSE_METHODS` set (free) → raw row count (cheap) → md5 hash of sorted element tuples (expensive, only on a row-count collision), cached to disk per `(method, biosample)`. This was **explicitly replaced** and should not be reintroduced.
- **Current** (sessions 10cdeda8/IGVF and b73cdace/CATlas, confirmed as the branch everyone converged on by `timestamp_signal` ordering): plain **union-find over predictions-table file paths** — two methods share a universe if, for every biosample they both cover, the predictions table points them at the **identical file path**, *or* if the pair appears in a hardcoded override list:
  ```python
  ALWAYS_SHARED_UNIVERSE = [('scE2G_multiome', 'scE2G_ATAC', 'scABC', 'distanceToTSS', 'ABC_distanceToTSS')]
  # ABC_distanceToTSS is the deprecated old name for distanceToTSS, kept for backward compatibility.
  ```
  Group label = `'-'.join(sorted(member methods))`. The representative member (the one covering the most biosamples) is the one whose candidate-universe tables are actually read for the entire group — this matters because all members are assumed to share **identical candidate element-gene coordinates**.

**Why it matters:** the candidate-universe (unthresholded) enrichment/recall heatmaps and tables are computed **once per universe group**, not once per method. Reading the wrong (thresholded) table here is exactly the bug that was fixed when the thresholded/candidate-universe split landed (see Enrichment section).

**Gotchas:**
- A method's group membership must be verified **consistent across every biosample it appears on in the run**, not assumed from a single biosample — if inconsistent, it falls back to its own solo group with a printed warning rather than a wrong merge.
- `distanceToTSS_Open4Gene` is **not** in `ALWAYS_SHARED_UNIVERSE` — it is scored against a genuinely different candidate element set (Open4Gene's own elements), so it should not be force-grouped with the scE2G family.
- For the four IGVF default methods specifically: `scABC` and `distanceToTSS` are documented (across eQTLEnrichment and CATlas sessions) as reading the **literal same raw file** as `scE2G_multiome` for a given biosample via the `MIRRORED_COLUMNS`/`SCE2G_DERIVED_BASELINES` mechanism — so under plain file-path grouping they would already group together even without the `ALWAYS_SHARED_UNIVERSE` override; the override exists for cases (e.g. a GM12878 override pointing at a different preprint file) where file-path identity alone would have missed it. 📌 VERIFY IN IGVF that your four default methods' predictions-table columns do in fact resolve to shared files, since that assumption underlies the "one shared candidate universe" claim in Default Configuration below.

### Default Configuration

The four default methods — **scE2G Multiome, scE2G scATAC, scABC, Distance to TSS** — map onto the methods table as four rows, each following the field layout in **Model Config** above (`scABC` and `distanceToTSS` are mirrored/derived columns of `scE2G_multiome`'s file per the mechanism above, not separate physical files). Per `ALWAYS_SHARED_UNIVERSE`, all four already belong to **one shared candidate universe** — candidate-universe enrichment/recall is computed once for the group, not four times.

**The default pair set is the intersection of pairs for which all four default methods have predictions** — this is the fairness constraint made concrete for the default run. Mechanism (`make_run_configs.py`'s "all_models" variant, documented in sessions 6a5dc004 and b4cc6b41):
```python
# keep only biosamples where EVERY listed model has a non-null prediction path
keep = df[["scE2G_multiome", "scE2G_ATAC", "scABC", "distanceToTSS"]].notna().all(axis=1)
```
Command pattern:
```
python workflow/scripts/generate_input/make_run_configs.py \
  -i config/local/predictions_<name>_validated.tsv \
  -o config/local \
  --all-models scE2G_multiome scE2G_ATAC scABC distanceToTSS \
  --candidate-models scE2G_multiome scE2G_ATAC scABC distanceToTSS \
  --prefix predictions_default
```
🔲 TODO: NEEDS VALIDATION — `make_run_configs.py` is documented primarily from eQTLEnrichment sessions; the IGVF session (10cdeda8) references equivalent `--include-biosample`/`--exclude-biosample`/`--models` flags but notes they were "not fully re-verified" on that branch. 📌 VERIFY IN IGVF the exact script path and flag names before relying on this command verbatim.

Because `scABC` and `distanceToTSS` mirror `scE2G_multiome`'s own file, in practice "all four have predictions" may reduce mostly to "`scE2G_multiome` **and** `scE2G_ATAC` both have predictions" for a biosample (since the mirrored columns inherit `scE2G_multiome`'s coverage automatically) — but `scE2G_ATAC` is a separately-generated file with potentially different coverage, so **do not assume this reduction; run the actual intersection check.**

---

## Running the Benchmark

### Pre-Run Checklist

Always complete these steps, in order, before submitting:

1. Predictions directory synced.
2. Biosample-tissue predictions config built from the synced directory.
3. Model config (methods table) updated with every method to be run.
4. A validator assigned (existing or newly built) for each method.
5. Prediction files validated — `main_validation_script.py` run, summary reviewed, no unexpected "missing" entries.
6. Config validation passes — i.e. a dry run completes without `WorkflowError`.
7. Dry run completes without errors and reports the expected job counts.

### Dry Run

**Always** run this before submitting to `sbatch` — it is also how config validation (above) actually gets triggered.
```
snakemake -n -s workflow/Snakefile --configfile <cfg> --quiet rules
```
equivalently, wrapped through the submission script in dry-run mode:
```
sbatch workflow/scripts/run_benchmark.sbatch <configfile.yaml> dryrun
```
**What to look for:** the per-rule job-count table prints with no `WorkflowError`; zero unexpected `WARNING` lines from the universe-grouping consistency check (a warning there means a method was force-split into its own solo group when you expected it to share a universe); job counts that reconcile with your expected biosample × method count (e.g. roughly `n_pairs_in_intersection × n_methods` for the per-biosample rules).

### Adding a New Method

Step-by-step for "run the eQTL benchmark with the addition of model X":

1. **Ensure predictions exist:** files for model X must exist under `<predictions_dir>/<dataset>/<biosample>/` following the `*_<biosample>_<ModelToken>[.e2g].tsv.gz` naming pattern, with header columns `ElementChr, ElementStart, ElementEnd, GeneSymbol` plus X's score column.
2. **Build/update the predictions table:** add X to `--models` when regenerating via `make_prediction_config.py` (see **Building a Biosample-Tissue Prediction Config** above), or if X's token isn't recognized, add an entry to the filename-suffix mapping used by `extract_model_name()` first.
3. **Update the model config:** add a row for `ModelX` to the methods table following the **Model Config** template above; add `ModelX` to the `methods:` list in the run YAML.
4. **Assign a validator:** use an existing one if X's score semantics match (see **Selecting an Existing Validator**), or build a new one (see **Building a New Validator**) and register it in `VALIDATOR_MAP`.
5. **Determine the pair set:** the active pair set becomes the intersection of pairs where **all** methods in this run — including X — have predictions. Rebuild the restricted predictions table with X added to `--all-models`:
   ```
   python workflow/scripts/generate_input/make_run_configs.py -i <validated_table> -o config/local \
     --all-models scE2G_multiome scE2G_ATAC scABC distanceToTSS ModelX --candidate-models ... --prefix predictions_withX
   ```
6. **Run validation:** `python workflow/IGVF_prediction_validation/main_validation_script.py -p config/local/predictions_withX.tsv -m <methods_table>`
7. **Dry run:** `snakemake -n -s workflow/Snakefile --configfile <cfg_pointing_at_predictions_withX> --quiet rules`
8. **Submit:** `sbatch workflow/scripts/run_benchmark.sbatch <cfg_pointing_at_predictions_withX> run`

⚠️ **FAIRNESS NOTE:** Adding model X changes the active pair set to the intersection of pairs where **all** methods (including X) have predictions. This run is **not** directly comparable to a prior default run unless the pair sets are identical — to check comparability, diff the biosample lists (or pair-count) between the two runs' predictions tables.

### Changing the Biosample-Tissue Pairs

Step-by-step for specifying a custom pair set. Per the task brief: **most users should be selecting a *subset* of the existing biosample–fine-grained-eQTL-tissue pairs** (wherever full coverage exists across the models being compared), not authoring brand-new pairs — new pairs require re-running the biology-based pairing process described in the Biosample-to-Tissue Pairs Config section, which is a much larger undertaking than this SKILL covers end-to-end. 🔲 TODO: NEEDS VALIDATION if you need to add genuinely new pairs, not just subset existing ones.

1. **Identify available pairs:** cross-reference the biosample-tissue pairs file against which biosamples actually have predictions for the method(s) you're comparing (see step 1 of **Combined** workflow below for the mechanics of checking file presence).
2. **Modify the pairs config:** for a *subset*, simply restrict the predictions table to the desired biosamples — e.g. via `make_run_configs.py --include-biosample <name>` (repeatable), which restricts both the `all_models`/`all_biosamples` outputs to only the named biosample(s). This does not require editing the tissue-pairs TSV itself.
3. **Validation:** re-run `main_validation_script.py` against the restricted table (cheap — it's a subset, but re-run it anyway in case the restriction changed which rows exist).
4. **Dry run:** confirm the reduced job count matches the expected subset size.
5. **Submit:** `sbatch workflow/scripts/run_benchmark.sbatch <cfg> run`.

⚠️ **FAIRNESS NOTE:** A run with a custom pair set is only comparable to other runs that used the **identical** pair set — not just the same number of pairs, the same specific biosample-tissue pairs.

### Combined: Adding Model X on Pairs Where X Has Predictions

Concrete walkthrough of "run the eQTL benchmark with the addition of model X on all pairs where predictions are present for model X":

1. **Identify which pairs have predictions for model X:** after generating X's column in the predictions table (step 2 of Adding a New Method), the biosamples with a non-null X entry are directly readable from that column — no separate helper script is documented beyond inspecting the table (e.g. `df[df['ModelX'].notna()]['biosample']`); 🔲 TODO: NEEDS VALIDATION whether a more direct CLI helper exists on the IGVF branch for this specific query.
2. **Set the pairs config to that set:** use `make_run_configs.py`'s `"all_biosamples"` variant with `--candidate-models` including `ModelX` — for each candidate model it checks `.notna().all()` across the retained biosamples and **drops any model that doesn't cover them all**, printing which biosamples caused the drop. In practice, run the `"all_models"` intersection (step 5 of Adding a New Method) including `ModelX` — this *is* "pairs where X has predictions" intersected with "pairs where the other methods being compared also have predictions."
3. **Add model X to the model config:** per **Model Config** above.
4. **Assign validator for model X:** per **Validator Package** above.
5. **Run validation:** `python workflow/IGVF_prediction_validation/main_validation_script.py -p <restricted_table> -m <methods_table>`
6. **Dry run:** `snakemake -n -s workflow/Snakefile --configfile <cfg> --quiet rules`
7. **Submit:** `sbatch workflow/scripts/run_benchmark.sbatch <cfg> run`

⚠️ **FAIRNESS NOTE:** This run compares all methods on the pair set defined by model X's prediction coverage intersected with the other methods'. Results are internally consistent but are **not** directly comparable to the default benchmark unless the default methods also have predictions for every pair in this set — to verify comparability, confirm this pair set is a **subset** of the default four-method intersection (see Default Configuration).

### Submitting with sbatch

```
sbatch workflow/scripts/run_benchmark.sbatch <configfile.yaml> [dryrun|run]
```
Internally, `run_benchmark.sbatch` builds and runs (parameters shown with their documented defaults):
```
snakemake --configfile <CONFIGFILE> --profile slurm.smk9 \
  --conda-prefix "$PWD/.snakemake/conda" \
  --jobs ${JOBS:-50} \
  --default-resources slurm_account=engreitz slurm_partition=owners,engreitz,normal \
      runtime="'${RUNTIME:-6h}'" slurm_extra="''" \
  --retries 3 --rerun-incomplete \
  [--dry-run --quiet rules   # dryrun mode
   | --keep-going]           # run mode — one failing method must not block generate_html_report
```
`JOBS` and `RUNTIME` are overridable env vars at submit time, e.g. `JOBS=50 RUNTIME=6h sbatch workflow/scripts/run_benchmark.sbatch <cfg> run`.

**Slurm parameters observed for the controller job** (📌 VERIFY IN IGVF — partition names are lab-specific): `--partition=engreitz,owners,normal --cpus-per-task=2 --mem=16GB --time=14:00:00` (controller only waits on children, doesn't do the compute itself). `validate_predictions.sbatch`: `--time=04:00:00 --cpus-per-task=1 --mem=32GB` (sized for loading a multi-million-row table into pandas during rescue).

**Monitoring:**
```
squeue --me --name=<jobname> --states=RUNNING --noheader -o "%i"
sacct -j <jobid> --format=JobID,JobName%14,State,Elapsed,MaxRSS,ExitCode --noheader
scontrol show job <jobid> | tr ' ' '\n' | grep -E "^(TimeLimit|Nice|Partition|NumCPUs|MinMemoryNode)="   # verify --default-resources overrides actually took effect
```

**Gotchas:**
- A bare `runtime=360` is parsed by Snakemake as **seconds**, not minutes — always pass a quoted unit string, `runtime="'6h'"`.
- `sbatch` snapshots the script content **at submit time**; editing the `.sbatch` file after submission does not affect an already-queued job.
- Sherlock's Slurm rejects jobs that are primarily a `sleep` call ("sleeper job detected"); any required grace period must use two separately `--begin`-scheduled jobs (see `stop_benchmarks.sbatch`'s two-stage signal/cancel pattern), not `sleep` inside one job.
- Two controllers sharing one `--conda-prefix` can **race** during first-time env creation (one process's env copy collides with the other's removal); once the env already exists, concurrent sharing is safe.
- `scancel`-ing a controller leaves a **stale lock** in `.snakemake/locks/`; the next invocation of the same config needs `snakemake --configfile <cfg> --unlock` first.
- Snakemake's lock is on the **overlapping output file set**, not the whole working directory — two configs with disjoint `outDir` values can run concurrently sharing one repo checkout.
- `--deadline` only gates whether Slurm will *start* a job it predicts can't finish in time — it does **not** terminate an already-running job; a hard stop needs `--time` (a real `TIMEOUT` kill) plus an explicit watchdog/cancel mechanism.

---

## CI Formula

See **Enrichment** above for full context. Repeated here for reference, consistent across sessions post-fix:
```
enrichment = (x1 / n1) / (x2 / n2)
SE(log enrichment) = sqrt( (n1-x1)/(x1*n1) + (n2-x2)/(x2*n2) )     # FIXED: previously sqrt(a) + b, understating SE
CI_low  = exp(log(enrichment) - z * SE_log_enr)
CI_high = exp(log(enrichment) + z * SE_log_enr)
p = phyper(x1, n1, n2, x1+x2, lower.tail=FALSE);  p.adjust(p, method="bonferroni")
```
Where `x1`=eQTL variants overlapping a predicted element, `n1`=total eQTL variants in tissue, `x2`=background SNPs overlapping, `n2`=total background SNPs. **Pooled ("AllMatches")** rows sum `x1,n1,x2,n2` across every matched biosample for a tissue *before* computing one ratio — not an average of per-pair ratios — and this pooling is a documented, unfixed, intentionally-parked limitation (treats dependent pairs as independent). See **Enrichment** for the full fix history (3 separate bugs, all fixed, with the two more serious ones — numerator/denominator thresholding — confirmed present-and-fixed on the IGVF branch per session 10cdeda8's timestamp data; 📌 VERIFY IN IGVF directly before trusting absolute numbers on your checkout).

---

## Outputs

All outputs live under `{outDir}` (set in the run YAML; gitignored, typically `results/<name>/`).

**Per-method subdirectories**, `{outDir}/{method}/`:
- `enrichmentTables/`: `enrichmentTable.{distance_min}to{distance_max}Kb.tsv` (5 files: 0–10, 10–100, 100–250, 250–1000, 0–30000 Kb — fixed-threshold enrichment) plus `enrichmentTable.acrossThresholds.tsv.gz` (full threshold-sweep, no distance dimension) plus `*.candidateUniverse.tsv` variants (unthresholded, one per universe-group representative).
- `recallTables/`: `recallTable.byDistance.tsv`, `recallTable.acrossThresholds.tsv.gz`, and candidate-universe variants.
- `intermediate/`: `GTExTissueBiosampleMap.tsv` (per-method biosample↔tissue pairing actually used), `basesPerEnhancerSet.tsv`, `thresholdSpan.tsv`, `nVariantsPerGTExTissue.tsv` (confirmed method-independent — safe as a shared tissue-axis ordering key).

**Cross-method plots/tables**, `{outDir}/plots/`:
- `enrichmentHeatmaps_thresholdedModel/{method}....pdf` and `enrichmentHeatmaps_candidateUniverse/{group-label}....pdf` (the thresholded-vs-unthresholded split — see Enrichment).
- `recallHeatmaps/` and `recallHeatmaps_candidateUniverse/`.
- `enrichmentRecall/enrichmentRecall.GTExTissue{tissue}.pdf` + matching `.tsv` + `.plotdata.rds`, and (IGVF-branch feature) an interactive `enrichmentRecall.GTExTissue{tissue}.html` per tissue with >1 biosample.
- `enrichmentAtRecall/enrichments.Recall{r}.GTExTissue{tissue}.pdf` + `pairwiseComparisons.Recall{r}.GTExTissue{tissue}.tsv` (Bonferroni-corrected pairwise significance — this is where to find "is model A significantly better than model B at matched recall").
- `allMatchedEnrichments.tsv`, `allMatchedPredictionMetrics.tsv` (pooled "AllMatches" summary tables), `thresholdedPerformanceComparison.pdf`.
- `colorPalette.tsv` (method → hex; ⚠️ gotcha — this is built from the methods table as a Snakemake **params** value, not an input, so editing a method's `color`/`pred_name_long` in the table does **not** auto-invalidate this rule; delete `colorPalette.tsv` manually to force regeneration).

**Top-level:** `{outDir}/benchmarkingReport.html` (built from `benchmarking_report.Rmd`) — includes its own embedded interactive `heatmaply` heatmap reading `enrichmentTable.0to30000Kb.tsv` directly. 🔲 TODO: NEEDS VALIDATION — per session 7d1019e3, this embedded heatmap automatically inherited the *numeric* enrichment fixes (same source file) but was **never rewired** to the thresholded/candidate-universe split, and its caption text was observed stale; confirm whether this was addressed on the IGVF branch before treating this particular widget as authoritative — use the dedicated `enrichmentHeatmaps_*` PDFs instead if in doubt.

Plotting/visual-design specifics (styling, canvas sizing, legend placement) are out of scope for this SKILL — see Known Issues below only insofar as a layout bug affected correctness of what's displayed, not aesthetics.

---

## Known Issues and Fixes

Ordered high → medium confidence. Each tagged `[session: project]`.

**High confidence:**
1. SE(log enrichment) computed as `sqrt(a) + b` instead of `sqrt(a + b)` — understated CI width. Fixed in `counts_to_enrichment_distance.R` / `counts_to_enrichment_across_thresholds.R`. `[10cdeda8/IGVF, b73cdace/CATlas, 1336b414/eQTLEnrichment]`
2. Enrichment numerator: the "0–30000Kb" (all-distances) bin skipped the score threshold filter entirely. Fixed by adding the same filter to that branch (commit `133acc0`). `[7d1019e3/eQTLEnrichment]`
3. Enrichment denominator: background overlap count was **never** thresholded, in any bin, for any method. Fixed by filtering `commonVarPredInt` by `score>=score.thresh` before counting (commit `ab9d27c`). Rated the more serious of the two threshold bugs. `[7d1019e3/eQTLEnrichment]`
4. Candidate-universe heatmap/table silently read the now-thresholded 0–30000Kb table after fixes #2/#3 landed, mislabeling thresholded numbers as "candidate universe." Fixed by adding dedicated `*CandidateUniverse`/`*candidateUniverse` rules with `threshold="-Inf"`. `[10cdeda8/IGVF, b73cdace/CATlas, 1336b414/eQTLEnrichment, 7d1019e3/eQTLEnrichment]`
5. Universe-grouping consistency check compared a per-biosample fingerprint **value** (hash/row-count) across a method's own biosamples — always "inconsistent" by construction even for correctly-shared universes, since the hash differs with genomic coordinates. Fixed by comparing the **frozenset of partner methods** sharing the tag on each biosample instead. `[7d1019e3/eQTLEnrichment]`
6. `process_predictions` used `set +o pipefail`; a `csvtk` column error exited 0 anyway, producing a silently-accepted 20-byte empty file. Fixed with explicit non-empty-output assertions. `[6a5dc004/eQTLEnrichment]`
7. `extract_model_name()` did an unanchored substring search, misclassifying files in biosample `b_quiescent`'s folder as model `SCENT`. Fixed by anchoring to the trailing filename-stem token. `[6a5dc004/eQTLEnrichment, 7d1019e3/eQTLEnrichment]`
8. Biosample folder names collide across IGVF datasets (e.g. `b_activated_d1` under two different datasets) — previously silently overwritten by a plain dict key. Fixed by requiring dataset-prefixed IDs, raising `DuplicateBiosampleConflict` on any remaining collision. `[6a5dc004/eQTLEnrichment, 7d1019e3/eQTLEnrichment, b4cc6b41/eQTLEnrichment]`
9. `inverse_predictor`/threshold-sign negation used Python truthiness on a config value; a blank cell (string `'None'` after `fillna`) is truthy, so blank cells were wrongly negated. Fixed with an explicit `value.strip().lower()=='true'` string check. `[10cdeda8/IGVF, 7d1019e3/eQTLEnrichment]`
10. A blank/missing `threshold` crashed the Snakefile with a bare `ValueError`. Fixed with a `has_threshold()` helper storing blank as `NaN` rather than raising; `THRESHOLDED_METHODS` built from that helper. `[10cdeda8/IGVF]`
11. `df[is.na(df)] = 0` coerced a missing threshold to literal `0`, drawing a spurious marker on the enrichment-recall plot. Fixed by selecting marker rows from the pre-coercion data with an explicit `is.na()` skip. `[10cdeda8/IGVF, b73cdace/CATlas]`
12. Single-tissue / single-biosample runs crashed (`print(str(tissue_dist))` referencing an undefined variable; `sd()` of one value returning `NA`, making a `filter(sd>0)` drop every row). Both fixed with explicit guards for `n==1`. `[10cdeda8/IGVF]`
13. `metadataTable` config key absent → `NULL` in R → `file.exists(NULL)` → `"argument is of length zero"` crash. Fixed with an explicit `!is.null() && length==1 && nzchar()` guard before `file.exists()`. `[10cdeda8/IGVF]`
14. `AbsCorrelationValidator.rescue()` required filenames ending `.e2g.tsv.gz`; Open4Gene's `*_Open4Gene.tsv.gz` (no `.e2g` infix) failed rescue entirely. Fixed by also accepting bare `.tsv.gz`. `[6a5dc004/eQTLEnrichment]`
15. `pgBoostValidator`'s rescue always recomputed `Percentile` via `rankdata(method='ordinal')` even when a correct `Percentile` already existed, disagreeing with the submitter's own `method='min'` convention on ~46/136,167 links at the 0.95 cutoff. Fixed to preserve an existing numeric column, only recomputing (with `method='min'`) when absent, and refusing to guess if partially numeric. `[6a5dc004/eQTLEnrichment]`
16. pgBoost `Score==1` saturated rows (~18–19% of links) retained a sentinel `Percentile` of `1.0` and were not excluded before thresholding. Fixed with an unconditional `pre_validate_and_reformat()` exclusion (see Validator gotcha above). `[6a5dc004/eQTLEnrichment, 7d1019e3/eQTLEnrichment]`

**Medium confidence:**
17. Fixed PDF canvas sizing (`ifelse(n>50,16,8)`) didn't adapt to actual matrix dimensions. Replaced with a `heatmap_canvas()` helper computing size from row/column counts plus fixed overhead. `[10cdeda8/IGVF, 6a5dc004/eQTLEnrichment]`
18. `htmlwidgets::saveWidget(selfcontained=TRUE)` left a stray `<name>_files` folder; plotly's `meta` field was an array, not a scalar, breaking a hover-highlight JS comparison. Both fixed. `[10cdeda8/IGVF]`
19. `render_report.py`'s pandoc fallback resolves relative image paths against the invoking shell's cwd, not the report's own directory — must `cd` into the report directory before rendering (workaround, not a code fix). `[7d1019e3/eQTLEnrichment]`
20. `validatedPredictionsTable` auto-derivation in the Snakefile used a prefix naming convention (`validated_<basename>`) that never matched `main_validation_script.py`'s actual suffix convention (`<basename>_validated`) — not fully resolved; point `config['predictionsTable']` directly at the validated table instead of relying on auto-derivation. `[6a5dc004/eQTLEnrichment]`

---

## Unresolved / Needs Validation

This is the user's action list for manual follow-up before trusting this SKILL blind.

**🔲 Needs validation (open questions in the source conversations themselves):**
- Whether `z>=1.96`-style threshold choices and model-specific sign conventions generalize correctly once the enrichment-formula fixes are applied — relevant if you add any model with a signed/divergent score, less relevant to the four IGVF defaults. `[7d1019e3/eQTLEnrichment]`
- Full end-to-end **real** (non-dry-run) pipeline execution had never been performed on the IGVF branch as of session 10cdeda8 — all verification there was via `snakemake -n` dry runs and hand-constructed mocked-R-script tests, never an actual `sbatch` run through the real DAG. **Run a real end-to-end dry-run-then-run on a small test config before trusting this SKILL's commands at scale.** `[10cdeda8/IGVF]`
- Interactive HTML enrichment-recall plots were only checked by parsing embedded plotly JSON programmatically, never opened and hovered in a real browser by a human. `[10cdeda8/IGVF]`
- `benchmarkingReport.html`'s embedded heatmaply widget was never rewired to the thresholded/candidate-universe split (see Outputs) — confirm current IGVF-branch status. `[7d1019e3/eQTLEnrichment]`
- Whether `make_run_configs.py` (all_models/all_biosamples intersection tooling) exists with the same flags on the IGVF branch, or needs porting — referenced but "not fully re-verified" there. `[10cdeda8/IGVF]`

**📌 VERIFY IN IGVF (collected from throughout this SKILL):**
- Exact conda env paths (`/oak/.../miniforge3/envs/eQTLEnv`, `/run_snakemake9`) are this user's Sherlock setup.
- Exact project-root / predictions-dir / config paths shown are illustrative examples from eQTLEnrichment/CATlas; confirm the IGVF branch's actual paths (e.g. `config/eQTL_tissue_matches.tsv` vs. the typo'd eQTLEnrichment path).
- Whether the `MIRRORED_COLUMNS`/`SCE2G_DERIVED_BASELINES` mechanism (scABC/distanceToTSS mirroring scE2G_multiome's file) is active for the four IGVF default methods — this directly affects what "all four have predictions" means for the default pair-set intersection.
- Whether the `universe_group` manual-override methods-table column (from eQTLEnrichment's Open4Gene work) exists on the IGVF branch at all; likely unnecessary for the four defaults.
- Slurm partition names (`engreitz,owners,normal`) are lab-specific.
- Whether commits `ab9d27c`/`133acc0`/`2af8c86` (the two serious enrichment-threshold fixes + heatmap split) are present on your actual IGVF checkout — confirmed present as of session 10cdeda8's timestamp, but verify against current `git log`.

---

## Source Attribution

| Section | Primary session(s) | Project(s) | Support |
|---|---|---|---|
| Overview / fairness framing | all 6 | all 3 | Well-supported — consistent across every session |
| Enrichment formula + fix history | 10cdeda8, 7d1019e3, b73cdace, 1336b414 | IGVF, eQTLEnrichment, CATlas | Well-supported, multiple agreeing sessions |
| Recall | 7d1019e3, 10cdeda8 | eQTLEnrichment, IGVF | Moderate — fewer independent confirmations than enrichment |
| Environment setup | 6a5dc004, 7d1019e3, 10cdeda8, b73cdace | eQTLEnrichment, IGVF, CATlas | Well-supported |
| Directory structure | 10cdeda8, 6a5dc004, 7d1019e3 | IGVF, eQTLEnrichment | Well-supported for eQTLEnrichment; IGVF paths are 📌 VERIFY |
| Predictions table format | all 6 | all 3 | Well-supported |
| Building prediction config from sync | 10cdeda8, 6a5dc004, b4cc6b41 | IGVF, eQTLEnrichment | Moderate — `make_prediction_config.py` mechanics mainly from eQTLEnrichment |
| Biosample-tissue pairs config | b4cc6b41, 1336b414, b73cdace, 6a5dc004 | eQTLEnrichment, CATlas | **Rests mainly on eQTLEnrichment/CATlas sessions — no IGVF-branch session documents this file in detail.** Treat with caution; resume a dedicated IGVF session if one exists. |
| Model config | 10cdeda8, 7d1019e3, 6a5dc004, 1336b414 | IGVF, eQTLEnrichment | Well-supported, consistent |
| Validator package (select + build) | 6a5dc004, 7d1019e3, 10cdeda8 | eQTLEnrichment, IGVF | Well-supported |
| Config validation (pipeline-level) | 10cdeda8, b73cdace, 7d1019e3 | IGVF, CATlas, eQTLEnrichment | Well-supported, multiple agreeing sessions |
| Universe grouping | 10cdeda8, 7d1019e3, b73cdace, 1336b414 | IGVF, eQTLEnrichment, CATlas | Well-supported on the *current* (file-path+override) logic; deprecated 4-tier logic is single-session (7d1019e3) and explicitly superseded |
| Default configuration / intersection logic | 6a5dc004, b4cc6b41 | eQTLEnrichment | **Rests on a single project (eQTLEnrichment) — no IGVF session confirms `make_run_configs.py` flags directly.** 📌 VERIFY IN IGVF before relying on the exact command. |
| Running the benchmark / sbatch | 6a5dc004, 7d1019e3, 10cdeda8, b73cdace | eQTLEnrichment, IGVF, CATlas | Well-supported |
| Known issues / bug fixes | all 6 | all 3 | Well-supported, cross-referenced |

---

## Synthesis Report

**Well-supported (multiple agreeing sessions):** the enrichment formula and its three-bug fix history, the thresholded-vs-candidate-universe split, universe grouping's current (file-path + `ALWAYS_SHARED_UNIVERSE`) logic, the model-config/validator three-way join, pipeline-level config validation, environment setup, and the sbatch/Slurm gotchas. These sections cross-check cleanly across eQTLEnrichment, CATlas, and IGVF sessions and can be trusted directly.

**Rests on a single session — treat with more caution:**
- **Biosample-tissue pairs config** (file format, path, resolution logic) comes almost entirely from eQTLEnrichment (`b4cc6b41`) and CATlas (`b73cdace`) sessions; no IGVF session in this batch documents the IGVF branch's own tissue-pairs file in comparable detail. If an IGVF-specific pairing session exists, **resume it** before trusting the path/column names here.
- **Default pair-set intersection tooling** (`make_run_configs.py`) is documented only from eQTLEnrichment (`6a5dc004`, `b4cc6b41`); the one IGVF session that touches it (`10cdeda8`) explicitly says it was "not fully re-verified" there. **Resume `10cdeda8` or a newer IGVF session** to confirm this script/flags exist as described before depending on it for the default run.
- The deprecated 4-tier universe-fingerprinting approach is single-session (`7d1019e3`) — correctly flagged as superseded, included only for historical context.

**Conflicts and how they were resolved:**
- Universe-grouping logic: two generations exist (4-tier fingerprinting vs. file-path+override). Resolved by timestamp ordering — the file-path+override version is confirmed present on both IGVF and CATlas branches in later sessions, so it is documented as current; the older one is marked deprecated.
- No other direct factual conflicts were found — the SE-formula fix, the two enrichment-threshold bugs, and the config-validation rules all agree across every session that touches them.

**Thin or missing — fill by resuming a specific session:**
- IGVF-branch-specific biosample-tissue pairs file path and `make_run_configs.py` flags — **resume `10cdeda8`** first (it's the one IGVF-project session in this batch), and if it doesn't have enough detail, a fresh IGVF session focused specifically on pair-set construction would need to be run.
- Recall section is thinner than Enrichment — if recall computation details matter for your use case beyond what's here, **resume `7d1019e3` or `10cdeda8`**, both of which touch `compute_prediction_table_distance.R` but not in exhaustive detail.

**`filtering_pairs_by_model` capability — gaps:** Partially documented. The *mechanism* (per-method `.dropna()` in `utils.smk`, and `make_run_configs.py`'s all_models/all_biosamples intersection logic) is clear and consistent across eQTLEnrichment sessions, but **no session gives a single, direct CLI command for "show me which pairs have predictions for model X"** — this SKILL had to compose that from adjacent primitives (inspecting the predictions-table column, or running the intersection with X in `--all-models`). 🔲 Before relying on the **Combined** workflow at scale, verify this composition actually produces the pair set you expect, ideally by resuming `10cdeda8` or `6a5dc004` to ask directly.

**Validator package — is it complete enough to build a new validator without resuming a source conversation?** Yes, with minor caveats. The extraction gives: the abstract base class's extension points (`is_score_valid()`, `rescue()`, `pre_validate_and_reformat()`), five worked concrete examples spanning the main check types (bounded-numeric, percentile-presence, unbounded-absolute, boolean-dtype, float-only), the `VALIDATOR_MAP` registration mechanism, and — critically — the one non-obvious gotcha (pgBoost chunked-read issue) that would otherwise cause a new validator to silently fail on large files. The one gap: no session shows the **literal file path** of `VALIDATOR_MAP` itself or the exact class-registration syntax (just that it's "a dict in `main_validation_script.py`") — a quick `grep -n VALIDATOR_MAP workflow/IGVF_prediction_validation/main_validation_script.py` on the actual IGVF checkout will resolve this in seconds and is not worth resuming a session for.


# eQTL benchmarking pipeline
We developed this "eQTL benchmarking pipeline” to assess whether predictive models of enhancer-gene regulatory interactions could (i) identify enhancers enriched for fine-mapped eQTL variants, and (ii) link those enhancers to their corresponding eQTL target genes. We have applied this pipeline using eQTLs fine-mapped with the Sum of Single Effects (SuSiE) model from the [GTEx V8 resource](https://gtexportal.org/home/downloads/adult-gtex/overview)  and from the [eQTL Catalogue release 7](https://www.ebi.ac.uk/eqtl/Data_access/).

For this analysis, we focus on distal noncoding eQTLs by by filtering out (i) coding sequences, 5′ and 3′ untranslated regions of protein-coding genes, and splice sites (within 10 bp of a intron–exon junction of a protein-coding gene) of protein-coding genes, and (ii) promoters (±250 bp from the gene TSS) of protein-coding genes. We define **enrichment** as (fraction of eQTL variants with PIP >= 50% that overlap predicted enhancers) / (fraction of all 1000G SNPs that overlap predicted enhancers). We define **recall (linking)** as as the fraction of eVariant-eGene pairs overlapping a predicted enhancer linked to the correct eGenes, and **recall (total)** as the fraction of eVariant-eGene pairs with the eVariant overlapping predicted a enhancer.

## Outputs
Running this pipeline will produce the following plots and corresponding tables, many of which are also included in a summary HTML report (`results/benchmarkingReport.html`). 
1. To compare comprehensively compare performance across predictors, we plot the enrichment and recall for user-specified eQTL biosample, prediction biosample matches, where predicted enhancers are defined as enhancer-gene links with a score greater than the provided threshold value.
	- `results/plots/thresholdedPerformanceComparison.pdf`: Boxplots for for enrichment, recall (total), and fraction of variants overlapping a predicted enhancer linked to the variant's eGene given that the variant overlaps an enhancer for eQTLs stratified by eGene-eVariant distance. Each datapoint represents variants from one eQTL biosample and predicted enhancers from the specified prediction biosample. Metrics for for eQTL biosamples with fewer than 20 variants are removed, and outliers are plotted as points.
	- `results/plots/allMatchedThresholdedPerformance.scatter.pdf`: Scatter plot of enrichment and recall (linking) for matched eQTL and prediction biosamples for variants with any eVariant-eGene distance. Points for eQTL biosamples with fewer than 20 variants are removed as above is performed.
2. To evaluate cell-type specificity of predictors, we plot enrichment and recall (linking) heatmaps across all eQTL biosamples, prediction biosample pairs. The eQTL tissue axis is ordered by number of eVariants in that tissue, fixed across every heatmap in a run so plots are comparable side by side (override with `heatmapTissueOrder: cluster` in the config to instead cluster each plot's own tissue axis by Pearson correlation, which will then vary between plots). The biosample axis is clustered by Pearson correlation. The number of base pairs predicted to be in enhancers for each biosample and the number of variants (enrichment) or eVariant-eGene pairs (recall) per eQTL biosample are plotted as well. Enrichment is reported two ways, split into two directories because they answer different questions and must not be compared to each other:
	- `results/plots/enrichmentHeatmaps_thresholdedModel/{predictor}.thresholdedModel.enrichmentHeatmap.withMetrics.pdf`: predicted enhancers are enhancer-gene links with a score above (or, for an inverse predictor, below) the predictor's own configured threshold — the plot title states the exact rule. Color scale is dynamic based on the range of enrichments for that predictor.
	- `results/plots/enrichmentHeatmaps_candidateUniverse/{universe}.candidateUniverse.enrichmentHeatmap.withMetrics.pdf`: no score threshold is applied — every candidate element-gene pair the predictor scored counts, answering how much eVariant signal sits in the candidate universe before any model scoring. One plot per *distinct* candidate universe. Two predictors share a universe when the predictions table points them at the same file for every biosample they both have (e.g. scABC and distanceToTSS, which read the scE2G multiome file); different files mean different universes, except that scE2G scATAC, scE2G multiome, scABC and distanceToTSS always share one.
	- `results/plots/recallHeatmaps/{predictor}.recallHeatmap.withMetrics.pdf`:  Color scale is constant across all predictors from 0 to 0.25. eQTL biosamples with constant recall for all prediction biosamples are removed. 
	- `results/plots/recallHeatmaps_candidateUniverse/{universe}.candidateUniverse.recallHeatmap.withMetrics.pdf` and `{universe}.candidateUniverse.recallTable.tsv`: recall with no score threshold, i.e. the most recall (linking) the candidate universe can reach, one per distinct universe as above. The table lists every biosample-eQTL tissue combination (all distances), whether it is a configured pairing, the number of eVariant-eGene pairs, how many have their eVariant in a candidate element and in one linked to the eGene, and the resulting recall (total) and recall (linking).
3. To evaluate predictors agnostic to a predefined threshold and account for the tradeoff between high enrichment and low recall, we calculate enrichment and recall (linking) curves across the full range of threshold values for specified eQTL biosample, prediction biosample matches. 
	- `results/plots.enrichmentRecall/enrichmentRecall.GTExTissue{tissue}.pdf`: One plot is created for each eQTL biosample with at least one specified prediction biosample match. One curve is plotted per prediction biosample matching the eQTL biosample, with straight lines connecting datapoints evaluated at thresholds yielding even steps across the recall range for all specified biosample matches. If multiple biosamples for the same predictor match one eQTL biosample, the curves will vary slightly in shade and the key will specify the specific match. Binary predictors are plotted as points corresponding to the positive predicted enhancers. The enrichment axis is cut off at 50 to eliminate spiking at very high thresholds.
	- `results/plots/enrichmentRecall/enrichmentRecall.GTExTissue{tissue}.html`: interactive version of the plot above, made only for eQTL tissues paired with more than one biosample. Hover a faint line to see which biosample (name and ID) and model it is, plus its threshold, recall and enrichment; hovering brightens that line. Click a model in the legend to hide it, double-click to show only that model. As in the PDF, the all-models plot is followed by one plot per model, all on the same axis ranges.
	- `results/plots.enrichmentRecall/enrichmentRecall.GTExTissueAllMatches.pdf`: Enrichment-recall curve as above where values are aggregated across all specified biosample matches for each predictor. Enrichment and recall are both aggregated by summing variant overlap counts across each biosample match before computing the metric.
	- `results/plots/normalizedEnrichmentRecall/normalizedEnrichmentRecall.GTExTissue{tissue}.pdf` (and `.tsv`), one per eQTL tissue plus `AllMatches`: enrichment-recall curves with enrichment divided by the enrichment of the model's own candidate universe (every element-gene pair the model scored, no score cutoff; as in the candidate-universe heatmaps). 1 (dashed line) = no better than calling every candidate, so the y-axis is the improvement the model's ranking adds over its universe; this makes models with differently enriched universes, and pairings with different universe enrichment, comparable. Thick lines pool counts over the tissue's pairings (all pairings for `AllMatches`) and divide by the universe enrichment pooled over the same pairings; faint lines are single pairings divided by their own universe enrichment. The dot with a bar is the configured threshold; its CI is the enrichment CI divided by the universe enrichment. The table has the pooled and per-pairing values (`level`), including the raw and universe enrichment.
	- `results/plots/enrichmentAtRecall/enrichmentAtRecall{recall}.GTExTissue{tissue}.pdf`: For user-specified recall (linking) values, a plot is generated for each eQTL biosample with at least one predictor biosample designated match. Each predictor biosample match is only included in the plot if they acheive a recall (linking) within 0.02 of the desired value; the exact recall corresponding to the threshold at which enrichment is plotted is shown in the legend. We calculate Bonferroni-adjusted p-values between all pairwise comparisons of included predictors using the approach of [comparing relative risks](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1125071/). Pairs with significant differences are not indicated on the plot because algorithimically calculating the y-position of the required bars is confusing. If fewer than two predictor biosample matches meet the designated recall, empty placeholder files are created. 

## Running the pipeline
1. Clone this repository
```
git@github.com:EngreitzLab/eQTLEnrichment.git
```
2. Activate a conda environment with [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) and [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) installed.
3. **Download the predictions.** Sync the Synapse E2G predictions project to a local
   directory. The expected layout is `<dataset>/<biosample>/<prediction files>`.
```
python workflow/IGVF_prediction_validation/make_prediction_config.py --sync -d $SCRATCH/<download_dir>
```
4. **Generate the predictions config table** from that download. Biosample IDs are
   emitted as `{dataset}_{biosample}`, because biosample folder names are *not*
   unique across datasets.
```
python workflow/IGVF_prediction_validation/make_prediction_config.py \
  -d $SCRATCH/<download_dir> \
  -t config/<tissue_pairings>.tsv \
  -o config/<predictions>.tsv
```
5. **Validate the table** (see below), then **review the result by hand** before
   launching. This step is deliberately manual: a single wrong cell silently removes
   a method from the benchmark.
6. Edit the main config file to point at your predictions table and methods table (see below)
7. Run pipeline
```
snakemake -j1 --configfile config/config_example.yml --use-conda
```

## Claude Code skill (`eqtl-benchmark-runner`)
A [Claude Code](https://docs.claude.com/en/docs/claude-code) skill for running, extending and interpreting this pipeline lives on the `skills` branch of [kaybrand/eQTLEnrichment](https://github.com/kaybrand/eQTLEnrichment), at `skills/eqtl-benchmark-runner/SKILL.md`. That branch holds only skills and is never merged into a code branch, so one copy serves every branch and checkout. You do not need any existing worktrees to use it.

1. **Get the `skills` branch into its own directory.**
   - No existing clone: clone just that branch (a few files, under 100 KB):
   ```
   git clone --branch skills --single-branch git@github.com:kaybrand/eQTLEnrichment.git eQTLEnrichment-skills
   ```
   - Already have a clone of the fork: add it as a worktree instead, which reuses the clone's history:
   ```
   git fetch origin skills
   git worktree add ../eQTLEnrichment-skills skills
   ```
   Never check out a different branch in this directory; Claude reads the skill from it.
2. **Make Claude Code load it.** Pick one route; using two loads the skill twice.
   - As its own plugin: start Claude with `--plugin-dir`, which can be repeated next to other plugins. The skill appears as `eQTLEnrichment-skills:eqtl-benchmark-runner`.
   ```
   claude --plugin-dir /path/to/eQTLEnrichment-skills
   ```
   - Through a local checkout of the lab's EngreitzLabAgents plugin: add a symlink and keep it out of that repo's git with a local-only exclude. Do not commit it to EngreitzLabAgents. The skill appears as `EngreitzLabAgents:eqtl-benchmark-runner`.
   ```
   cd /path/to/EngreitzLabAgents
   ln -s /path/to/eQTLEnrichment-skills/skills/eqtl-benchmark-runner skills/eqtl-benchmark-runner
   echo "/skills/eqtl-benchmark-runner" >> .git/info/exclude
   ```
   - As a personal skill, available in every session: symlink it into your Claude skills folder (`~/.claude/skills/`, or `$CLAUDE_CONFIG_DIR/skills/` if you set one).
   ```
   ln -s /path/to/eQTLEnrichment-skills/skills/eqtl-benchmark-runner ~/.claude/skills/eqtl-benchmark-runner
   ```
3. **Check:** start a new Claude session and run `/skills`; `eqtl-benchmark-runner` should be listed.
4. **Updates:** `git pull` in `eQTLEnrichment-skills/` to get changes. Make edits there too, and commit and push them on the `skills` branch.

## Prediction file format validation
Use the IGVF_prediction_validation package to verify and correct the format of predictions in the predictions configuration table before feeding it to the eQTL pipeline.
For models not already represented by the pgBoostValidator, SCARlinkValidator, AbsCorrelationValidator, scE2GValidator, or DefaultValidator, define your own class inheriting from the BaseValidator abstract class.

Validation is a **manual pre-step**; it is not wired into the Snakemake DAG.
```
python workflow/IGVF_prediction_validation/main_validation_script.py \
  -p config/<predictions>.tsv -m config/<methods>.tsv --check-all-rows
```
The validated table is written next to the input as `<name>_validated.tsv`, with a
`NaN` in place of any file that failed. **Inspect those `NaN`s before running**: an
empty cell means that biosample is dropped for that method, and the pipeline will not
warn you. Check in particular that every cell points at a file belonging to the right
biosample *and* the right method — a cell pasted into the wrong column is what caused
the 2025-10-19 run to lose all of one method's enrichment tables.

## Config files

**Config checks.** When the pipeline starts (including `snakemake -n`), it checks the three config files and every prediction file's header, and stops with a numbered list of all problems before running any job. It checks: each listed method has a methods-table row and a predictions-table column; thresholds are numbers or blank; `boolean`, `inverse_predictor` and `nonneg_only` are `TRUE`/`FALSE`; a `boolean` predictor's threshold is 1; `geneUniverse` files and prediction files exist; each prediction file has `ElementChr`, `ElementStart`, `ElementEnd`, `GeneSymbol` and its `score_col`; biosample names are unique; and each method has at least one biosample paired with an eQTL tissue. Whether a threshold keeps any links can only be known by reading the whole file, so that check runs in `process_predictions`, one of the first jobs, which stops with the observed score range if no link passes (skipped for methods with no threshold).
The pipeline requires three config files. The required inputs are outlined below.
1. **Main config file** - an example file is included at: `config/config_example.yml`. The following fields are required:
	- **outDir:** output directory name
	- **methodsTable:** configuration file with information about predictors (see below for specifications)
	- **methods:** list of methods from `methodsTable` to be benchmarked
	- **eQTLVariants:** gzipped, tab-separated file of eQTL variants with the following columns: `chr`, `start`, `end`, `varID_hg38` (variant ID), `gene_hgnc` (eGene), `tissue`, `pip`. Here, we use fine-mapped variants from GTEx or the eQTL Catalogue. For the GTEx dataset, we implement additional filtering for variants in a credible set and expressed with a TPM>1 in their respective tissue 
		- The unprocessed file of SuSiE-fine-mapped GTEx variants is available on Synapse [here](https://www.synapse.org/#!Synapse:syn52264297), and the expression data from GTEx is available [here](https://www.synapse.org/#!Synapse:syn52264240). The code to process these files into the format for this pipeline is under `workflow/scripts/generate_inputs`.
		- The code to generate an input file from the eQTL Catalogue is available upon request.
	- Genomic annotations **chrSizes**, **partition**, **TSS**, **GTExGeneUniverse** (all provided here in the `resources` directory)
	- **bgVariants:** a .bed file of 1000G SNPs  with columns (no header) `chr`, `start`, `end`, `rsid`. The list of background variants we use can be downloaded on Synapse [here](https://www.synapse.org/#!Synapse:syn52264319) and was collated from https://alkesgroup.broadinstitute.org/LDSCORE/baseline_v1.1_hg38_annots/
	- The following parameters:
		- **distances:** a list of values, in Kb, for eGene-eVariant distance-stratified benchmarking
		- **recalls:** a list of recall values at which to compare enrichments. May be left blank.
		- **thresholdPIP:** PIP threshold to filter variants
		- **thresholdPval:"** p-value threshold used in statistical tests
		- **nThresholdSteps:** number of thresholds to include when generating enrichment-recall curves. We recommend using 20-50 steps.
  
2. **Methods config table** - an example file is included at: `config/methods_config_example.tsv`. The following columns (with a header) are required, where each row represents one predictive method:
   - **method:** predictive method name, to match listed methods in the main config file outlined above
   - **pred_name_long:** a string with the full method name to be used in plots
   - **threshold:** score threshold value to be used to generate boxplots stratified by distance and enrichment heat maps (outputs (1) and (4)). We used the threshold values corresponding to 70% recall from our CRISPR benchmarking pipeline for each model. *Optional:* leave it blank if the predictor has no threshold. That method then still gets enrichment-recall curves (with no threshold dot) and enrichment-at-recall plots, but is left out of everything that needs a threshold: the thresholded enrichment heatmaps, recall heatmaps, distance-stratified comparisons and the report's thresholded tables. The report lists which methods had no threshold.
   - **score_col:** name of the column in prediction files with the prediction score
   - **color:** hex code specifying what color to use in plots for this method. May be left blank and will be automatically filled.
   - **inverse_predictor:** `TRUE` if high scores correspond to lower prediction confidence (eg. distance to TSS, or a p-value), otherwise `FALSE`
   - **boolean:** `TRUE` if this is a binary 0 or 1 predictor, otherwise `FALSE`
   - **nonneg_only:** *(optional)* `TRUE` to treat a negative score as "not predicted" and drop those rows, otherwise `FALSE`. Used for signed-correlation predictors such as Signac and Cicero, where an anti-correlated element-gene pair is not a weak prediction but no prediction at all. Without it those rows form the bottom of the threshold span and the enrichment-recall curve sweeps through them before reaching score 0. Omitting the column entirely is equivalent to `FALSE` for every method.

3. **Predictions config table** - an example file is included at: `config/predicitons_config_example.tsv`. The table contains the following columns (with a header), where each row represents one prediction biosample:
   - **biosample:** name of prediction biosample
   - **GTExTissue:** comma-separated list eQTL biosamples matching this prediciton biosample, where the eQTL biosample name corresponds to the `tissue` column in the eQTL variant input file
   - columns with names corresponding to each predictive method, with names matching the methods listed in the main config. These columns contain file paths to the enhancer-gene predictions for the corresponding method. The prediction files should be unthresholded and may be gzipped or uncompressed. They must contain the following columns at minimum, using the IGVF standard-format names: `ElementChr`, `ElementStart`, `ElementEnd` (corresponding to the enhancer elements), `GeneSymbol` (gene symbol), and a score column corresponding to the `score_col` indicated in the methods configuration table. Lines beginning with `#` are treated as a metadata header and skipped.
   - **The column name must match the `method` value exactly.** The pipeline stops at startup if a listed method has no matching column.
   - Biosample names must be unique. Prefer the `{dataset}_{biosample}` IDs produced by `make_prediction_config.py`: several biosample folder names recur across datasets.

## Known limitations

### Pooled statistics treat matched pairs as independent
The "AllMatches" results and the pooled per-tissue curves in the enrichment-recall plots sum counts over every matched (eQTL tissue, prediction biosample) pair at each score threshold. For each pair, `x1` is the tissue's eVariants that overlap the biosample's enhancers, `n1` is the tissue's eVariants, `x2` is the background SNPs that overlap the same enhancers, and `n2` is all background SNPs. Then:
- enrichment = (Σx1 / Σn1) / (Σx2 / Σn2)
- recall = Σ(eVariant-eGene pairs linked to the correct gene) / Σ(eVariant-eGene pairs)
- the 95% CI uses SE(log enrichment) = sqrt((n1−x1)/(x1·n1) + (n2−x2)/(x2·n2)) on the summed counts, the p-value is a hypergeometric test on the same sums, and model-to-model comparisons at a fixed recall use a z-test built from those SEs.

These formulas treat the summed counts as independent draws, and they are not. A tissue matched to k biosamples contributes the same eVariants k times. The same background SNPs are reused in every pair. Related biosamples share many enhancers, and eQTL tissues share eVariants. The effective sample size is therefore smaller than the summed counts, so pooled CIs are too narrow and pooled p-values too small. The pooled point estimates themselves are unaffected. Models scored on the same variants also give correlated estimates, which the model-to-model z-test ignores. Read pooled error bars and p-values as optimistic until this is resolved.

## Works-in-progress
- Edit terminology in configuration, code, file names to not be specific to GTEx (replace "GTEx tissue" with "eQTL biosample")
- Implement benchmarking of groups of prediction biosamples against a single eQTL biosamples to better map to hetergenous tissues
- Integrate computations from enrichment and recall into a single function
- Implement baseline predictors computed internally (e.g. distance to TSS, random expressed gene with 500kb)
- Review/improve the "enrichment at recall" analysis

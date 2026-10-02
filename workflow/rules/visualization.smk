def display_name(method):
	# blank methods-table cells arrive as the string "None" (see Snakefile)
	name = methods_config.loc[method, "pred_name_long"]
	return method if name == "None" else name

def threshold_rule_text(method):
	# methods_config thresholds are negated for inverse predictors (see Snakefile)
	col = methods_config.loc[method, "score_col"]
	t = methods_config.loc[method, "threshold"]
	if is_config_true(methods_config.loc[method, "inverse_predictor"]):
		return f"{col} <= {-t:g}"
	return f"{col} >= {t:g}"

# color palette
rule generate_color_palette:
	params:
		user_inputs = [methods_config.loc[method, "color"] for method in config["methods"]],
		names = config["methods"],
		methods_config = config["methodsTable"]
	output:
		colorPalette = os.path.join(config["outDir"], "plots", "colorPalette.tsv"),		
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script: 
		os.path.join(config["codeDir"], "visualization", "color_palette.R")

# gather data for enrichment recall curve per GTEx tissue/biosample match per method
rule gather_enrichment_recall:
	input: 
		predTable = os.path.join(config["outDir"], "{method}", "recallTables", "recallTable.acrossThresholds.tsv.gz"),
		enrichmentTable = os.path.join(config["outDir"], "{method}", "enrichmentTables", "enrichmentTable.acrossThresholds.tsv.gz"),
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	output:
		ERCurveTable = os.path.join(config["outDir"], "{method}", "enrichmentRecallTables", "GTExTissue{GTExTissue}.Biosample{biosample}.tsv")
	resources:
		mem_mb = determine_mem_mb
	script:
		os.path.join(config["codeDir"], "enrichment_recall", "gather_enrichment_recall.R")


# make table for "aggregate" enrichment/recall for all tissue/biosample matches per method
rule gather_enrichment_recall_aggregate:
	input: 
		predTable = os.path.join(config["outDir"], "{method}", "recallTables", "recallTable.acrossThresholds.tsv.gz"),
		enrichmentTable = os.path.join(config["outDir"], "{method}", "enrichmentTables", "enrichmentTable.acrossThresholds.tsv.gz"),
	conda: 
		os.path.join(config["envDir"], "eQTLEnv.yml")
	output:
		ERCurveTable = os.path.join(config["outDir"], "{method}", "enrichmentRecallTables", "AllMatches.tsv")
	resources:
		mem_mb = determine_mem_mb
	script:
		os.path.join(config["codeDir"], "enrichment_recall", "gather_aggregate_enrichment_recall.R")	

# combined enrichment recall curves-- one per GTEx tissue-- across methods
def get_table_files(GTExTissue):
	if (GTExTissue=="AllMatches"):
		tissue_files  = [os.path.join(config["outDir"], method, "enrichmentRecallTables", "AllMatches.tsv") for method in config["methods"]]
	else:
		substr = "GTExTissue" + GTExTissue + ".Biosample"
		all_files = flatten([[os.path.join(config["outDir"], method, "enrichmentRecallTables", f"GTExTissue{tissue}.Biosample{biosample}.tsv") for tissue, biosample in zip(methods_config.loc[method, "GTExTissue_map"], methods_config.loc[method, "biosample_map"])] for method in config['methods']])
		tissue_files =  [x for x in all_files if substr in x]
	return tissue_files

rule plot_enrichment_recall_curve:
	input:
		colorPalette = os.path.join(config["outDir"], "plots", "colorPalette.tsv"),
		enrichmentRecall_files = lambda wildcards: get_table_files(wildcards.GTExTissue),
		# Per-tissue view only: raw per-threshold counts for pooling biosamples,
		# unique eVariants per tissue, the pipeline's biosample-tissue pairings,
		# and one candidate-universe enrichment table per distinct universe.
		enrichmentTable_acrossThresholds = expand(os.path.join(config["outDir"], "{method}", "enrichmentTables", "enrichmentTable.acrossThresholds.tsv.gz"), method=config["methods"]),
		variantsFile = os.path.join(config["outDir"], "variants", "GTExVariants.filteredForUniverse.tsv.gz"),
		biosampleTissueMaps = expand(os.path.join(config["outDir"], "{method}", "intermediate", "GTExTissueBiosampleMap.tsv"), method=config["methods"]),
		candidateEnrichmentTables = [os.path.join(config["outDir"], UNIVERSE_GROUP_REPRESENTATIVE[label], "enrichmentTables", "enrichmentTable.0to30000Kb.candidateUniverse.tsv") for label in sorted(UNIVERSE_GROUP_MEMBERS)],
	params:
		score_thresholds = [methods_config.loc[method_name, "threshold"] for method_name in config["methods"]],
		methods = config["methods"],
		predictionsTable = config["predictionsTable"],
		candidateUniverseNames = [" / ".join(display_name(m) for m in sorted(UNIVERSE_GROUP_MEMBERS[label])) for label in sorted(UNIVERSE_GROUP_MEMBERS)]
	output:
		er_combined = os.path.join(config["outDir"],  "plots", "enrichmentRecall", "enrichmentRecall.GTExTissue{GTExTissue}.pdf"),
		er_combined_table = os.path.join(config["outDir"],  "plots", "enrichmentRecall", "enrichmentRecall.GTExTissue{GTExTissue}.tsv")
	resources:
		mem_mb = determine_mem_mb
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "visualization", "plot_enrichment_recall.R")
		
# plot enrichment with confidence intervals across methods, once per GTEx tissue and user-defined recalls
rule plot_enrichment_with_ci:
	input:
		colorPalette = os.path.join(config["outDir"], "plots", "colorPalette.tsv"),
		enrichmentRecall_files = lambda wildcards: get_table_files(wildcards.GTExTissue)
	params:
		thresholdPval = config["thresholdPval"]
	output:
		enr_at_recall = os.path.join(config["outDir"],  "plots", "enrichmentAtRecall", "enrichments.Recall{recall}.GTExTissue{GTExTissue}.pdf"),
		enr_at_recall_table = os.path.join(config["outDir"],  "plots", "enrichmentAtRecall", "enrichments.Recall{recall}.GTExTissue{GTExTissue}.tsv"),
		sign_table = os.path.join(config["outDir"],  "plots", "enrichmentAtRecall", "pairwiseComparisons.Recall{recall}.GTExTissue{GTExTissue}.tsv")
	resources:
		mem_mb = determine_mem_mb
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "visualization", "plot_enrichment_with_ci_at_recall.R")

# generate final comparison plot (triple boxplot)
# thresholded comparisons: only methods with a threshold (see Snakefile)
enrichmentTables_distance=  [os.path.join(config["outDir"], method, "enrichmentTables", f"enrichmentTable.{distance_min}to{distance_max}Kb.tsv") for distance_min, distance_max in zip(config["distances_min"], config["distances_max"]) for method in THRESHOLDED_METHODS]
predTables_distance = [os.path.join(config["outDir"], method, "recallTables", "recallTable.byDistance.tsv") for method in THRESHOLDED_METHODS]
maps = [os.path.join(config["outDir"], method, "intermediate", "GTExTissueBiosampleMap.tsv") for method in THRESHOLDED_METHODS]

rule plot_thresholded_matched_performance_comparison:
	input:
		enrichmentTable_files = enrichmentTables_distance,
		predTable_files = predTables_distance,
		map_files = maps,
		colorPalette = os.path.join(config["outDir"], "plots", "colorPalette.tsv"),
	params:
		distances_min = config["distances_min"],
		distances_max = config["distances_max"]
	output:
		outFile = os.path.join(config["outDir"], "plots", "thresholdedPerformanceComparison.pdf"),
		enrAllTable = os.path.join(config["outDir"], "plots", "allMatchedEnrichments.tsv"),
		predictionMetrics = os.path.join(config["outDir"], "plots", "allMatchedPredictionMetrics.tsv"),
		outScatter = os.path.join(config["outDir"], "plots", "allMatchedThresholdedPerformance.scatter.pdf"),
		outViolins = os.path.join(config["outDir"], "plots", "allMatchedThresholdedPerformance.violins.pdf")
	resources:
		mem_mb = determine_mem_mb
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "visualization", "plot_thresholded_performance_comparison.R")

# plot all thresholded enrichments + recalls
rule plot_all_thresholded_metrics:
	input:
		enrichmentTable_files = enrichmentTables_distance,
		predTable_files = predTables_distance,
		colorPalette = os.path.join(config["outDir"], "plots", "colorPalette.tsv"),
	params:
		distances_min = config["distances_min"],
		distances_max = config["distances_max"]
	output:
		outFile = os.path.join(config["outDir"], "plots", "allThresholdedEnrichmentRecall.byDistance.pdf"),
	resources:
		mem_mb = determine_mem_mb
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "visualization", "plot_all_thresholded_enrichment_recall.R")

# plot n variants per tissue
rule plot_variants_per_tissue:
	input:
		variantsPerTissue =  [os.path.join(config["outDir"], method, "intermediate", "nVariantsPerGTExTissue.tsv") for method in config["methods"]]
	params:
		distances_min = config["distances_min"],
		distances_max = config["distances_max"],
		tissues_matched = GTExTissues_matched,
		methods = config["methods"]
	output:
		out_plot = os.path.join(config["outDir"], "plots", "variantsPerTissue.pdf"),
		out_plot_matched = os.path.join(config["outDir"], "plots", "variantsPerTissue_matched.pdf"),
		out_table = os.path.join(config["outDir"], "plots", "avgVariantsPerTissue.tsv")
	resources:
		mem_mb = determine_mem_mb
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "visualization", "plot_number_variants.R")

# heatmaps -- split into two clearly-labeled products (see plot_enrichment_heatmap.R
# header and determine_candidate_universe_groups() in utils.smk):
#   thresholded model:  this method's own thresholded predictions, one per method.
#   candidate universe: the unthresholded candidate element-gene set, one per
#                        distinct universe (methods sharing one collapse into a
#                        single plot named for all of them).
rule plot_enrichment_heatmaps_thresholded:
	input:
		enrichmentTable = os.path.join(config["outDir"], "{method}", "enrichmentTables", "enrichmentTable.0to30000Kb.tsv"),
		enhancerSizes =  os.path.join(config["outDir"], "{method}", "intermediate", "basesPerEnhancerSet.tsv"),
		biosample_tissue_map = os.path.join(config["outDir"], "{method}", "intermediate", "GTExTissueBiosampleMap.tsv")
	params:
		mode = "thresholded",
		method = lambda wildcards: wildcards.method,
		methodsTable = config["methodsTable"],
		p_threshold = config["thresholdPval"],
		metadata = config.get("metadataTable"),
		tissueOrder = config.get("heatmapTissueOrder", "nVariants")
	output:
		outFile_combined =  os.path.join(config["outDir"], "plots", "enrichmentHeatmaps_thresholdedModel", "{method}.thresholdedModel.enrichmentHeatmap.withMetrics.pdf"),
		outFile_alone = os.path.join(config["outDir"], "plots", "enrichmentHeatmaps_thresholdedModel", "{method}.thresholdedModel.enrichmentHeatmap.pdf")
	resources:
		mem_mb = determine_mem_mb
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "visualization", "plot_enrichment_heatmap.R")

rule plot_enrichment_heatmaps_candidateUniverse:
	input:
		# Any group member works: by construction they share identical candidate
		# element-gene coordinates (and therefore identical enhancer-set bp), so
		# the representative's own tables already speak for the whole group.
		enrichmentTable = lambda wildcards: os.path.join(config["outDir"], UNIVERSE_GROUP_REPRESENTATIVE[wildcards.universeGroup], "enrichmentTables", "enrichmentTable.0to30000Kb.candidateUniverse.tsv"),
		enhancerSizes = lambda wildcards: os.path.join(config["outDir"], UNIVERSE_GROUP_REPRESENTATIVE[wildcards.universeGroup], "intermediate", "basesPerEnhancerSet.candidateUniverse.tsv"),
		biosample_tissue_map = lambda wildcards: os.path.join(config["outDir"], UNIVERSE_GROUP_REPRESENTATIVE[wildcards.universeGroup], "intermediate", "GTExTissueBiosampleMap.tsv")
	params:
		mode = "candidateUniverse",
		memberMethods = lambda wildcards: sorted(UNIVERSE_GROUP_MEMBERS[wildcards.universeGroup]),
		p_threshold = config["thresholdPval"],
		metadata = config.get("metadataTable"),
		tissueOrder = config.get("heatmapTissueOrder", "nVariants")
	output:
		outFile_combined = os.path.join(config["outDir"], "plots", "enrichmentHeatmaps_candidateUniverse", "{universeGroup}.candidateUniverse.enrichmentHeatmap.withMetrics.pdf"),
		outFile_alone = os.path.join(config["outDir"], "plots", "enrichmentHeatmaps_candidateUniverse", "{universeGroup}.candidateUniverse.enrichmentHeatmap.pdf")
	resources:
		mem_mb = determine_mem_mb
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "visualization", "plot_enrichment_heatmap.R")

rule plot_recall_heatmaps:
	input:
		recallTable = os.path.join(config["outDir"], "{method}", "recallTables", "recallTable.byDistance.tsv"),
		enhancerSizes =  os.path.join(config["outDir"], "{method}", "intermediate", "basesPerEnhancerSet.tsv"),
		biosample_tissue_map = os.path.join(config["outDir"], "{method}", "intermediate", "GTExTissueBiosampleMap.tsv")
	params:
		mode = "thresholded",
		title = lambda wildcards: f"Model: {display_name(wildcards.method)}  |  THRESHOLDED predictions",
		subtitle = lambda wildcards: ("recall = fraction of eVariant-eGene pairs whose eVariant is in a predicted enhancer linked to its eGene\n"
		                              f"predicted enhancer: {threshold_rule_text(wildcards.method)}")
	output:
		outFile_combined =  os.path.join(config["outDir"], "plots", "recallHeatmaps", "{method}.recallHeatmap.withMetrics.pdf"),
		outFile_alone = os.path.join(config["outDir"], "plots", "recallHeatmaps", "{method}.recallHeatmap.pdf")
	resources:
		mem_mb = determine_mem_mb
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "visualization", "plot_recall_heatmap.R")

# Recall with no score cutoff: the most recall each distinct candidate
# element-gene universe can reach. One plot and table per universe group.
rule plot_recall_heatmaps_candidateUniverse:
	input:
		recallTable = lambda wildcards: os.path.join(config["outDir"], UNIVERSE_GROUP_REPRESENTATIVE[wildcards.universeGroup], "recallTables", "recallTable.byDistance.candidateUniverse.tsv"),
		enhancerSizes = lambda wildcards: os.path.join(config["outDir"], UNIVERSE_GROUP_REPRESENTATIVE[wildcards.universeGroup], "intermediate", "basesPerEnhancerSet.candidateUniverse.tsv"),
		biosample_tissue_map = lambda wildcards: os.path.join(config["outDir"], UNIVERSE_GROUP_REPRESENTATIVE[wildcards.universeGroup], "intermediate", "GTExTissueBiosampleMap.tsv")
	params:
		mode = "candidateUniverse",
		title = "CANDIDATE E2G UNIVERSE (no score threshold applied)",
		subtitle = lambda wildcards: ("max. recall = fraction of eVariant-eGene pairs whose eVariant is in any candidate element linked to its eGene\n"
		                              "shared by: " + ", ".join(display_name(m) for m in sorted(UNIVERSE_GROUP_MEMBERS[wildcards.universeGroup])))
	output:
		outFile_combined = os.path.join(config["outDir"], "plots", "recallHeatmaps_candidateUniverse", "{universeGroup}.candidateUniverse.recallHeatmap.withMetrics.pdf"),
		outFile_alone = os.path.join(config["outDir"], "plots", "recallHeatmaps_candidateUniverse", "{universeGroup}.candidateUniverse.recallHeatmap.pdf"),
		table = os.path.join(config["outDir"], "plots", "recallHeatmaps_candidateUniverse", "{universeGroup}.candidateUniverse.recallTable.tsv")
	resources:
		mem_mb = determine_mem_mb
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "visualization", "plot_recall_heatmap.R")

# html report
rule generate_html_report:
	input:
		colorPalette = os.path.join(config["outDir"], "plots", "colorPalette.tsv"),
		# thresholded inputs are empty when no method has a threshold
		enrAllTable = [os.path.join(config["outDir"], "plots", "allMatchedEnrichments.tsv")] if THRESHOLDED_METHODS else [],
		predictionMetrics = [os.path.join(config["outDir"], "plots", "allMatchedPredictionMetrics.tsv")] if THRESHOLDED_METHODS else [],
		er_combined = expand(os.path.join(config["outDir"],  "plots", "enrichmentRecall", "enrichmentRecall.GTExTissue{GTExTissue}.tsv"), GTExTissue=GTExTissues_plus_all),
		enrMatrices_CRISPRthresh = expand(os.path.join(config["outDir"], "{method}", "enrichmentTables", "enrichmentTable.0to30000Kb.tsv"), method=THRESHOLDED_METHODS),
		recallTables_CRISPRthresh = expand(os.path.join(config["outDir"], "{method}", "recallTables", "recallTable.byDistance.tsv"), method=THRESHOLDED_METHODS),
		nVariants = os.path.join(config["outDir"], "plots", "avgVariantsPerTissue.tsv")
	params:
		methods = THRESHOLDED_METHODS,
		distances_min = config["distances_min"],
		distances_max = config["distances_max"],
		distances = config["distances"],
		GTExTissues_matched = GTExTissues_plus_all
	resources:
		mem_mb = determine_mem_mb
	output:
		htmlReport = os.path.join(config["outDir"], "benchmarkingReport.html")
	conda:
		os.path.join(config["envDir"], "eQTLEnv.yml")
	script:
		os.path.join(config["codeDir"], "visualization", "benchmarking_report.Rmd")







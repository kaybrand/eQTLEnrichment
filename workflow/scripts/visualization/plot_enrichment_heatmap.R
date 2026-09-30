suppressPackageStartupMessages({library(ggplot2)
        library(dplyr)
        library(scales)
        library(tidyr)
    	library(tibble)
        library(egg)
		library(data.table)})

# Produces one of two clearly-labeled enrichment heatmaps, selected by
# snakemake@params$mode:
#
#   "thresholded"      -- eQTL variants in THIS METHOD's thresholded predictions
#                         vs background in the same thresholded predictions.
#                         Title states the model and the exact threshold rule.
#   "candidateUniverse" -- eQTL variants in the UNTHRESHOLDED candidate
#                         element-gene universe vs background in the same
#                         unthresholded universe. One plot per *distinct*
#                         universe (methods sharing a universe -- see
#                         determine_candidate_universe_groups() in utils.smk --
#                         collapse into one plot named for all of them).
#
# Before this split, the single enrichmentHeatmap.pdf always read the
# "0-30000Kb" (all-distances) enrichmentTable, which -- prior to the two
# score-threshold fixes in this branch -- was implicitly a candidate-universe
# comparison despite README.md describing it as thresholded. The two modes
# here are that same code path, explicit about which one it is.

# Size the canvas from the actual matrix instead of a fixed 8x8 / 10x16.
#
# Rows are prediction biosamples and columns are eQTL tissues; both vary widely
# between runs (1 to ~140 biosamples, up to ~56 tissues). The legend, the 60-degree
# rotated tissue labels and the margins consume a roughly constant number of inches,
# so that overhead is ADDED to the per-cell budget rather than divided into it --
# otherwise small heatmaps lose their panel entirely to the legend and large ones
# still crush their rows. Bounds keep the output openable in a PDF viewer.
heatmap_canvas <- function(n_rows, n_cols,
                           row_in = 0.16, col_in = 0.16,
                           w_overhead = 3.2, h_overhead = 2.4,
                           min_w = 6, max_w = 24,
                           min_h = 4, max_h = 24) {
	w <- w_overhead + col_in * max(n_cols, 1)
	h <- h_overhead + row_in * max(n_rows, 1)
	c(width  = min(max(w, min_w), max_w),
	  height = min(max(h, min_h), max_h))
}

# Exact "what counts as predicted" clause from the methods table itself, so the
# plot cannot drift from the config that actually generated the numbers.
is_predicted_clause <- function(methods_table_file, method_name) {
	mt <- fread(methods_table_file, sep = "\t", header = TRUE)
	row <- dplyr::filter(mt, method == method_name)
	stopifnot(nrow(row) == 1)
	score_col <- row$score_col[1]
	thresh    <- row$threshold[1]
	inverse   <- toupper(as.character(row$inverse_predictor[1])) %in% c("TRUE", "TRUE.")
	op <- if (inverse) "<=" else ">="
	paste0(score_col, " ", op, " ", thresh)
}

main <- function() {
	mode = snakemake@params$mode  # "thresholded" | "candidateUniverse"
	stopifnot(mode %in% c("thresholded", "candidateUniverse"))

	# input data
	enrTableFile = snakemake@input$enrichmentTable # 0-30000kb table
	enhSizeFile = snakemake@input$enhancerSizes # per method; biosample / base pairs
	p_threshold = snakemake@params$p_threshold %>% as.numeric()
	outFile_combined = snakemake@output$outFile_combined
	outFile_solo = snakemake@output$outFile_alone
	metadata_file = snakemake@params$metadata
	biosample_tissue_file = snakemake@input$biosample_tissue_map
	# "nVariants" (default, fixed across every heatmap in a run) or "cluster"
	# (per-plot hclust on this table's own enrichment values -- the original
	# behavior; reshuffles tissue order between methods/modes because it is a
	# property of the enrichment numbers, not of the tissue set).
	tissue_order_by = snakemake@params$tissueOrder
	if (is.null(tissue_order_by) || !nzchar(tissue_order_by)) tissue_order_by = "nVariants"
	stopifnot(tissue_order_by %in% c("nVariants", "cluster"))
	highlight_pairs = TRUE

	enr = fread(enrTableFile, sep="\t", header=TRUE)
	enr = dplyr::filter(enr, nVariantsGTExTissue>20, is.finite(enrichment))
	enr$Biosample[enr$Biosample=="Cells_EBV-transformed_lymphocytes"] = "Cells_EBV_transformed_lymphocytes"
	enhSizes = fread(enhSizeFile, sep="\t", header=TRUE)
	colnames(enhSizes) = c("Biosample", "enhBp")

    # add base pairs per biosample
	enr = left_join(enr, enhSizes, by="Biosample")
	enr$enhMb = enr$enhBp/1e6

	# Add experiemental data source to enrichment table.
	# metadataTable is optional: config.get() yields None when the key is absent,
	# which arrives here as NULL, and file.exists(NULL) is logical(0) -- so a bare
	# if() on it fails with "argument is of length zero" rather than skipping the
	# block. Guard the length before testing existence.
	if (!is.null(metadata_file) && length(metadata_file) == 1 &&
	    nzchar(metadata_file) && file.exists(metadata_file)) {
		metadata = fread(metadata_file, sep=',', header=TRUE, select = c('folder_name', 'Dataset', 'SampleSummaryShort'))
		enr <- left_join(enr,
			metadata,
			by = c("Biosample" = "folder_name"),
			relationship = "many-to-one"
    	)
	}

    # pivot for biosample clustering (and, in "cluster" mode, tissue clustering)
    M = dplyr::select(enr, Biosample, GTExTissue, enrichment) %>% distinct() %>%
		pivot_wider(names_from=GTExTissue, values_from = enrichment) %>% column_to_rownames("Biosample") %>% drop_na()
	M[is.na(M)] <- 0

	# Tissue (column) order.
	if (tissue_order_by == "nVariants") {
		# Fixed across every heatmap in a run: nVariantsGTExTissue comes from the
		# eQTL variant file filtered only by distance bin, not from any method's
		# predictions, so it is identical for every method (and every mode) and
		# safe to use as one shared order.
		tissue_order_names = dplyr::select(enr, GTExTissue, nVariantsGTExTissue) %>% distinct() %>%
			arrange(desc(nVariantsGTExTissue), GTExTissue) %>% pull(GTExTissue)
		order_tissues = match(tissue_order_names, colnames(M))
		order_tissues = order_tissues[!is.na(order_tissues)]
	} else {
		n_tissues = ncol(M)
		if (n_tissues >= 2) {
			tissue_dist <- dist(1-cor(M))
			tissue_dist[is.na(tissue_dist)] <- 0
			order_tissues =  hclust(tissue_dist, method = "ward.D2")$order
		} else if (n_tissues == 1) {
			order_tissues = 1
		} else {
			print("There must be at least 1 tissue")
			order_tissues = 0
		}
	}

	# Get base biosample distance (disssimilarity of clusters) matrix
	biosample_dist <- dist(1-cor(t(M)))
	biosample_dist[is.na(biosample_dist)] <- 0
	biosample_dist_matrix <- as.matrix(biosample_dist) # Convert to a matrix for easier manipulation

	if ("Dataset" %in% colnames(enr)) {
		# Keep clusters from the same experiment close to each other.  Loop through each unique experiment
		for (experiment in unique(enr$Dataset)) {
		# Get the biosamples belonging to the current experiment
		biosamples_in_experiment <- enr$Biosample[enr$Dataset == experiment]

		# Loop through all pairs of biosamples within the experiment
		for (i in 1:length(biosamples_in_experiment)) {
			for (j in (i+1):length(biosamples_in_experiment)) {
			# Get the biosample names
			biosample1 <- rownames(biosample_dist_matrix)[rownames(biosample_dist_matrix) %in% biosamples_in_experiment[i]]
			biosample2 <- rownames(biosample_dist_matrix)[rownames(biosample_dist_matrix) %in% biosamples_in_experiment[j]]

			# Reduce the distance between these biosamples.  The amount of reduction
			# can be adjusted.  A value of 0 would force them to be clustered
			# perfectly together.  0.5 halves the distance between them.
			distance_reduction_factor <- 0.5
			biosample_dist_matrix[biosample1, biosample2] <- biosample_dist_matrix[biosample1, biosample2] * distance_reduction_factor
			biosample_dist_matrix[biosample2, biosample1] <- biosample_dist_matrix[biosample2, biosample1] * distance_reduction_factor # Ensure symmetry
			}
		}
		}

		# Convert the modified distance matrix back to a dist object
		constrained_biosample_dist <- as.dist(biosample_dist_matrix)

		#  Perform hierarchical clustering on the modified distance matrix
		if (attr(constrained_biosample_dist, "Size") >= 2) {
			order_biosamples = hclust(constrained_biosample_dist, method = "ward.D2")$order
		} else {
			order_biosamples = 1
		}
	} else {
		if (attr(biosample_dist, "Size") >= 2) {
			order_biosamples = hclust(biosample_dist, method = "ward.D2")$order
		} else {
			order_biosamples = 1
		}
	}
	# Turn the GTExTissue and Biosample columns into ordered categorical data
	enr$Biosample = factor(enr$Biosample, levels=rownames(M)[order_biosamples], ordered=TRUE)
	enr$GTExTissue = factor(enr$GTExTissue, levels=colnames(M)[order_tissues], ordered=TRUE)

	# set plotting params
	colors = c("#132b44","#154159","#18586b","#226f7a","#368785","#509f8e","#70b694","#95cd9a","#bee2a1","#ebf6ac") # dark blue to white yellow
	na_color = "#ffffff"

    # Define biosample-tissue pairs to highlight
    biosample_tissue_map <- fread(biosample_tissue_file, sep='\t', header=TRUE)
    colnames(biosample_tissue_map) = c("Biosample", "GTExTissue")
	# Create a data frame for the rectangles
    if (highlight_pairs) {
        rect_data <- biosample_tissue_map %>%
        left_join(
            data.frame(GTExTissue = levels(enr$GTExTissue), x = 1:length(levels(enr$GTExTissue))), # Get x positions
            by = "GTExTissue"
        ) %>%
        left_join(
            data.frame(Biosample = levels(enr$Biosample), y = 1:length(levels(enr$Biosample))), # Get y positions
            by = "Biosample"
        ) %>%
        mutate(
            xmin = x - 0.5,  # Adjust to center the rectangle on the cell
            xmax = x + 0.5,
            ymin = y - 0.5,
            ymax = y + 0.5
        )
    }

	# find max enrichment
	enr_lim <- dplyr::filter(enr, p_adjust_enr < p_threshold, nVariantsOverlappingEnhancers / nVariantsGTExTissue > 0.01)
	max_value <- round(quantile(enr_lim$enrichment, 0.99), 1)
	max_value <- max(1, max_value)
	lims = c(0, max_value)

	# Count what is actually drawn, not what survived into M: drop_na() above can
	# remove a biosample from M while it remains a factor level in enr and is still
	# rendered as a lane.
	n_rows_plotted = max(nlevels(enr$Biosample), nrow(M), 1)
	n_cols_plotted = max(nlevels(enr$GTExTissue), ncol(M), 1)
	canvas = heatmap_canvas(n_rows_plotted, n_cols_plotted)
	wd = canvas[["width"]]
	ht = canvas[["height"]]

	# Title/subtitle: exactly what is being measured, spelled out, plus which
	# model(s)/threshold this specific plot uses. Two lines, not one -- a single
	# long line ran off the edge of the canvas and got silently clipped.
	if (mode == "thresholded") {
		method_name = snakemake@params$method
		methods_table_file = snakemake@params$methodsTable
		predicted_clause = is_predicted_clause(methods_table_file, method_name)
		title_main = paste0("Model: ", method_name, "  |  THRESHOLDED predictions")
		title_sub = paste(
			"enrichment = (fraction of eVariants in predicted enhancers) / (fraction of common SNPs in predicted enhancers)",
			paste0("predicted enhancer: ", predicted_clause, " (eQTL variants AND background both restricted to this rule)"),
			sep = "\n"
		)
		banner_color = "#1a7a1a"
	} else {
		member_methods = snakemake@params$memberMethods
		title_main = "CANDIDATE UNIVERSE (no score threshold applied)"
		title_sub = paste(
			"candidate enrichment = (fraction of eVariants in candidate elements) / (fraction of common SNPs in candidate elements)",
			paste0("shared by: ", paste(member_methods, collapse = ", ")),
			sep = "\n"
		)
		banner_color = "#B2182B"
	}

    # heat map alone
	just_enr = ggplot(enr, aes(x=GTExTissue, y=Biosample, fill=enrichment)) +
		geom_tile() +
		scale_fill_gradientn(colors=colors, oob=scales::squish, na.value=na_color, limits=lims, name="Enrichment") +
		labs(title = title_main, subtitle = title_sub) +
		theme_minimal() + theme(axis.text = element_text(size = 7), axis.title = element_blank(), axis.text.x = element_text(angle=60, hjust=1),
			legend.position='top',  legend.direction='horizontal', legend.text=element_text(size=7), legend.title=element_text(size=7),
			plot.title = element_text(size = 10, face = "bold", color = banner_color),
			plot.subtitle = element_text(size = 7.5, color = banner_color)) +
        # Add black boxes around biosample-tissue pairs
		geom_rect(data = rect_data,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              color = "black",  # Customize the border color
              fill = NA,       # Make the rectangle transparent
              linewidth = 1)      # Customize the border thickness


	# plots for grid
	enr_grid  = ggplot(enr, aes(x=GTExTissue, y=Biosample, fill=enrichment)) +
		geom_tile() +
		scale_fill_gradientn(colors=colors, oob=scales::squish, na.value=na_color, limits=lims, name="Enrichment") +
		labs(title = title_main, subtitle = title_sub) +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_blank(), axis.text.x = element_blank(),
			legend.position='top',  legend.direction='horizontal', legend.text=element_text(size=7), legend.title=element_text(size=7),
			plot.title = element_text(size = 10, face = "bold", color = banner_color),
			plot.subtitle = element_text(size = 7.5, color = banner_color))
	if (highlight_pairs) {
		# Add black boxes around biosample-tissue pairs
        enr_grid = enr_grid +
        geom_rect(data = rect_data,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              color = "black",  # Customize the border color
              fill = NA,       # Make the rectangle transparent
              linewidth = 1)      # Customize the border thickness
    }

	nVar = dplyr::select(enr, GTExTissue, nVariantsGTExTissue) %>% distinct()
	var_count = ggplot(nVar, aes(x=GTExTissue, y=nVariantsGTExTissue)) +
		geom_bar(stat="identity", width=0.5) +
		ylab("# variants in tissue") +  xlab("") +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.text.x = element_text(angle=60, hjust=1))

	enhSizes = dplyr::select(enr, Biosample, enhMb) %>% distinct()
	enh_size = ggplot(enhSizes, aes(x=Biosample, y=enhMb)) +
		geom_bar(stat="identity", width=0.5) +
		ylab("Enhancer set size\n(Mb)") + xlab("") +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), axis.text.y = element_blank()) +
		coord_flip()

	blank = ggplot() + theme_void()

	# The bottom panel (variant count per tissue) needs enough absolute height for
	# its rotated x-axis tick labels, y-axis tick labels, and axis title to render
	# without crowding, independent of how tall the heatmap itself is. A fixed
	# ratio (heights=c(2, 0.2)) squeezes that panel to ~9% of the figure regardless
	# of absolute size -- fine at a small n_rows, crushes the y-axis labels once
	# the heatmap grows tall. Give it a fixed inch budget instead, added on top of
	# the heatmap's own height, so it never shrinks below what its labels need.
	bottom_panel_in = 1.3
	total_height = ht + bottom_panel_in
	assembled = egg::ggarrange(enr_grid, enh_size, var_count, blank, nrow=2, ncol=2, heights=c(ht, bottom_panel_in), widths=c(2, 0.3))

	ggsave(outFile_solo, just_enr, width=wd, height=ht, limitsize=FALSE)
	ggsave(outFile_combined, assembled, width=wd*1.15, height=total_height*1.05, limitsize=FALSE)

}

main()

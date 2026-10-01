suppressPackageStartupMessages({library(ggplot2)
        library(dplyr)
        library(scales)
        library(tidyr)
    	library(tibble)
        library(egg)
		library(data.table)})
    
# Size the canvas from the actual matrix instead of a fixed 8x8 / 10x16.
#
# Rows are prediction biosamples and columns are eQTL tissues; both vary widely
# between runs. The legend, the 60-degree rotated tissue labels and the margins
# consume a roughly constant number of inches, so that overhead is ADDED to the
# per-cell budget rather than divided into it. Bounds keep the output openable.
# Kept in step with the identical helper in plot_enrichment_heatmap.R.
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

cluster_tissues_biosamples <- function(res){
	M = dplyr::select(res, Biosample, GTExTissue, recall.linking) %>%
		pivot_wider(names_from=GTExTissue, values_from=recall.linking) %>% column_to_rownames("Biosample")
	M[is.na(M)] = 0
	
	n_tissues = ncol(M)
	if (n_tissues >= 2) { # at least 2 tissues
		tissue_dist = dist(1-cor(M))
		tissue_dist[is.na(tissue_dist)] = 0 
		order_tissues = hclust(tissue_dist, method = "ward.D2")$order
	} else if (n_tissues == 1) {
		# Nothing to cluster. This previously printed str(tissue_dist), a variable that
		# only exists in the n_tissues >= 2 branch, so a single-tissue run errored here.
		order_tissues = 1
	} else {
		print("There must be at least 1 tissue")
		order_tissues = 0
	}
	tissues_ordered = colnames(M)[order_tissues]

	n_biosamples = nrow(M)
	if (n_biosamples >=2) {
		biosample_dist = dist(1-cor(t(M)))
		biosample_dist[is.na(biosample_dist)] = 0
		order_biosamples = hclust(biosample_dist, method="ward.D2")$order
		biosamples_ordered = rownames(M)[order_biosamples]
	} else {
		biosamples_ordered = rownames(M)[1]
	}
	
	res$GTExTissue = factor(res$GTExTissue, levels=tissues_ordered, ordered=TRUE)
	res$Biosample = factor(res$Biosample, levels=biosamples_ordered, ordered=TRUE)
	return(res)
}

main <- function() {
	# input data
	recallTableFile = snakemake@input$recallTable # by distance
	enhSizeFile = snakemake@input$enhancerSizes # per method; biosample / base pairs
	#p_threshold = snakemake@params$p_threshold %>% as.numeric()
	outFile_combined = snakemake@output$outFile_combined
	outFile_solo = snakemake@output$outFile_alone

# 	if (grepl("ArchR", enhSizeFile)) {
#     file_name <- paste0("heatmap_ArchR", ".RData")
#     save.image(file = file_name)
#     quit(save = 'yes')
# 	}

	recall = data.table::fread(recallTableFile, sep="\t", header=TRUE)
	recall = dplyr::filter(recall, distance_min == 0, distance_max==30000,
		is.finite(recall.linking), !is.na(recall.linking), total.variants>20)
	# Remove tissues where recall is constant across biosamples -- uninformative for
	# clustering. Only meaningful with at least two biosamples to vary across: with a
	# single biosample sd() is NA for every tissue, NA > 0 is NA, and the filter would
	# silently drop every row and leave an empty table for the rest of the script.
	if (dplyr::n_distinct(recall$Biosample) > 1) {
		recall_sum = group_by(recall, GTExTissue) %>%
			summarize(sd_recall = sd(recall.linking)) %>%
			dplyr::filter(!is.na(sd_recall), sd_recall>0)
		recall = dplyr::filter(recall, GTExTissue %in% recall_sum$GTExTissue)
	}

	recall$Biosample[recall$Biosample=="Cells_EBV-transformed_lymphocytes"] = "Cells_EBV_transformed_lymphocytes"
	enhSizes = fread(enhSizeFile, sep="\t", header=TRUE) %>% distinct()
	colnames(enhSizes) = c("Biosample", "enhBp") 

    # add base pairs per biosample
	recall = left_join(recall, enhSizes, by="Biosample")
	recall$enhMb = recall$enhBp/1e6
	
    # cluster to get orders
   recall = cluster_tissues_biosamples(recall)

	# set plotting params
	colors = c("#edf8fb","#b3cde3", "#8c96c6", "#8856a7", "#810f7c")
	na_color = "#ffffff"
	lims = c(0, 0.25) 
	n_biosamples = recall$Biosample %>% unique() %>% length()
	n_tissues_plotted = recall$GTExTissue %>% unique() %>% length()
	canvas = heatmap_canvas(n_biosamples, n_tissues_plotted)
	wd = canvas[["width"]]
	ht = canvas[["height"]]

    # heat map alone
	just_recall = ggplot(recall, aes(x=GTExTissue, y=Biosample, fill=recall.linking)) + 
		geom_tile() +
		scale_fill_gradientn(colors=colors, oob=scales::squish, na.value="#FFFFFF", limits=lims, name="Recall (linked to correct eGene)") +
		theme_minimal() + theme(axis.text = element_text(size = 7), axis.title = element_blank(), axis.text.x = element_text(angle=60, hjust=1),
			legend.position='top',  legend.direction='horizontal', legend.text=element_text(size=7), legend.title=element_text(size=7))
		
	# plots for grid
	recall_grid  = ggplot(recall, aes(x=GTExTissue, y=Biosample, fill=recall.linking)) + 
		geom_tile() +
		scale_fill_gradientn(colors=colors, oob=scales::squish, na.value="#FFFFFF", limits=lims, name="Recall (linked to correct eGene)") +
		theme_minimal() + theme(axis.text = element_text(size = 7), axis.title = element_blank(), axis.text.x = element_blank(),
			legend.position='top',  legend.direction='horizontal', legend.text=element_text(size=7), legend.title=element_text(size=7))

	nVar = dplyr::select(recall, GTExTissue, total.variants) %>% distinct()
	var_count = ggplot(nVar, aes(x=GTExTissue, y=total.variants)) +
		geom_bar(stat="identity", width=0.5) +
		ylab("# eGene/eVariant pairs") +  xlab("") +
		theme_minimal() + theme(axis.text = element_text(size = 7), axis.text.x = element_text(angle=60, hjust=1))

	enhSizes = dplyr::select(recall, Biosample, enhMb) %>% distinct()
	enh_size = ggplot(enhSizes, aes(x=Biosample, y=enhMb)) +
		geom_bar(stat="identity", width=0.5) +
		ylab("Enhancer set size\n(Mb)") + xlab("") +
		theme_minimal() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), axis.text.y = element_blank()) + 
		coord_flip()

	blank = ggplot() + theme_void()

	# Same fixed inch budget as plot_enrichment_heatmap.R: a 2:0.2 height ratio
	# leaves the bottom panel too short for its tick labels and axis title,
	# especially when the heatmap has only a few rows.
	bottom_panel_in = 1.3
	total_height = ht + bottom_panel_in
	assembled = egg::ggarrange(recall_grid, enh_size, var_count, blank, nrow=2, ncol=2, heights=c(ht, bottom_panel_in), widths=c(2, 0.3))

	ggsave(outFile_solo, just_recall, width=wd, height=ht, limitsize=FALSE)
	# widths 2:0.3 adds the enhancer-size panel on the right; grow the canvas to match.
	ggsave(outFile_combined, assembled, width=wd*1.15, height=total_height*1.05, limitsize=FALSE)

}

main()

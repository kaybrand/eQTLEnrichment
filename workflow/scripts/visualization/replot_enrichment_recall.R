## Regenerate an enrichment-recall curve from a table a completed run already wrote,
## optionally dropping methods.
##
## Purpose: produce a figure variant (e.g. "same plot, without the scATAC curve")
## without re-running the benchmark and without hand-editing the PDF. Editing the
## vector output directly would leave the legend, the colour mapping and the y-axis
## limit describing curves that are no longer drawn; regenerating recomputes all
## three from the data.
##
## The plotting block is kept deliberately identical to plot_enrichment_recall.R so
## the variant is visually consistent with the pipeline's own output. If that script
## changes, mirror the change here.
##
## Usage:
##   Rscript workflow/scripts/visualization/replot_enrichment_recall.R \
##     --table results/<run>/plots/enrichmentRecall/enrichmentRecall.GTExTissueAllMatches.tsv \
##     --exclude scE2G_ATAC \
##     --out /path/to/out.pdf

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(colorspace)
  library(data.table)
  library(optparse)
})

option_list <- list(
  make_option("--table", type = "character", default = NA,
              help = "enrichmentRecall.GTExTissue<tissue>.tsv written by a pipeline run"),
  make_option("--out", type = "character", default = NA,
              help = "output PDF path"),
  make_option("--exclude", type = "character", default = "",
              help = "comma-separated method names to drop from the plot"),
  make_option("--tissue", type = "character", default = "AllMatches",
              help = "tissue label; only affects the x-axis caption wording"),
  make_option("--methods-table", type = "character", dest = "methods_table", default = NA,
              help = "config_methods TSV; refreshes pred_name_long so a display-name change reaches the legend without re-running the pipeline")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.na(opt$table) || is.na(opt$out)) stop("--table and --out are required")

df <- fread(opt$table, sep = "\t")

# pred_name_long is baked into the saved table, so a rename in the methods config
# would otherwise not appear here. Re-join it by method when a table is supplied.
if (!is.na(opt$methods_table)) {
  mc <- fread(opt$methods_table, sep = "\t", select = c("method", "pred_name_long"))
  changed <- df %>% dplyr::select(method, old = pred_name_long) %>% distinct() %>%
    dplyr::inner_join(mc, by = "method") %>% dplyr::filter(old != pred_name_long)
  if (nrow(changed) > 0) {
    for (i in seq_len(nrow(changed))) {
      cat("Renamed:", changed$old[i], "->", changed$pred_name_long[i], "\n")
    }
  }
  df <- dplyr::select(df, -pred_name_long) %>% dplyr::left_join(mc, by = "method")
}

exclude <- if (nzchar(opt$exclude)) trimws(strsplit(opt$exclude, ",")[[1]]) else character(0)
if (length(exclude) > 0) {
  absent <- setdiff(exclude, unique(df$method))
  if (length(absent) > 0) {
    stop("--exclude names not present in the table: ", paste(absent, collapse = ", "))
  }
  before <- nrow(df)
  df <- dplyr::filter(df, !method %in% exclude)
  cat("Dropped", paste(exclude, collapse = ", "),
      sprintf("(%d of %d rows removed)\n", before - nrow(df), before))
}
cat("Methods plotted:", paste(sort(unique(df$method)), collapse = ", "), "\n")

# Same legend rule as plot_enrichment_recall.R: the " (biosample)" suffix only
# disambiguates multi-biosample matches, and is noise for the all_matches aggregate.
df$key <- ifelse(df$Biosample == "all_matches",
                 df$pred_name_long,
                 paste0(df$pred_name_long, " (", df$Biosample, ")"))

df <- df[order(df$threshold), ]
df$key <- factor(df$key)

cp <- dplyr::select(df, method, key, hex) %>% distinct()
pred_colors <- cp$hex
names(pred_colors) <- cp$key

df_plot   <- dplyr::filter(df, recall.linking > 0)   # drop the (0,0) points
df_binary <- dplyr::filter(df_plot, nPoints == 2) %>% dplyr::filter(threshold == 1)
df_thresh <- dplyr::filter(df, threshold == score_threshold)  # the configured operating point
df_plot   <- dplyr::filter(df_plot, nPoints > 2)

# Recomputed after exclusion, so the axis and legend describe only what is drawn.
ylim <- min(50, max(df_plot$enrichment))
n_keys <- length(unique(cp$key))
n_legend_cols <- ceiling(n_keys / 12)

n_var_label <- if (opt$tissue == "AllMatches") {
  paste0(mean(df_plot$total.variants), " variant-biosample pairs")
} else {
  paste0(mean(df_plot$total.variants), " variants")
}
x_label <- paste0("Recall (variants overlapping prediction linked to eGene)\n", n_var_label)

g <- ggplot(data = df_plot, aes(x = recall.linking, y = enrichment, color = key)) +
  geom_line(linewidth = 0.75) +
  geom_point(data = df_binary, aes(x = recall.linking, y = enrichment, color = key),
             size = 3, shape = 16) +
  geom_linerange(data = df_thresh, aes(ymin = CI_enr_low, ymax = CI_enr_high), linewidth = 0.75) +
  geom_point(data = df_thresh, aes(x = recall.linking, y = enrichment, color = key),
             size = 3, shape = 16) +
  scale_color_manual(values = pred_colors) +
  ylab("Enrichment (eQTLs versus common variants)") + xlab(x_label) +
  labs(col = "Predictor") +
  coord_cartesian(ylim = c(0, ylim)) +
  theme_classic() + theme(axis.text = element_text(size = 7, color = "#000000"),
    axis.title = element_text(size = 8), axis.ticks = element_line(color = "#000000"),
    legend.text = element_text(size = 7), legend.title = element_text(size = 8),
    legend.position = "right", legend.direction = "vertical",
    aspect.ratio = 1) +
  guides(col = guide_legend(nrow = 12))

width <- 3 + 3 * n_legend_cols
ggsave(opt$out, g, width = width, height = 4)
cat("Wrote", opt$out, sprintf("(%.1f x 4.0 in)\n", width))

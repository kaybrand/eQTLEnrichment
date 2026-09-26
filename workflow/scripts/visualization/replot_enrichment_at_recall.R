## Regenerate an "enrichment at a fixed recall" bar plot from a table a completed
## run already wrote, optionally dropping methods and refreshing display names.
##
## Companion to replot_enrichment_recall.R; see that file for why regenerating beats
## editing the PDF. The plotting block mirrors plot_enrichment_with_ci_at_recall.R --
## keep the two in sync if that script changes.
##
## Dropping a method also changes the pairwise significance table, since the
## Bonferroni correction is over all pairs. That table is recomputed and rewritten
## alongside the plot rather than left describing a set of curves that no longer
## matches the figure.
##
## Usage:
##   Rscript workflow/scripts/visualization/replot_enrichment_at_recall.R \
##     --table results/<run>/plots/enrichmentAtRecall/enrichments.Recall0.05.GTExTissueAllMatches.tsv \
##     --exclude scE2G_ATAC --recall 0.05 \
##     --methods-table config/config_methods_IGVF_Annual_Meeting_2026.tsv \
##     --out <out.pdf> --out-sign <pairwise.tsv>

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(optparse)
})

option_list <- list(
  make_option("--table", type = "character", default = NA,
              help = "enrichments.Recall<r>.GTExTissue<tissue>.tsv from a pipeline run"),
  make_option("--out", type = "character", default = NA, help = "output PDF path"),
  make_option("--out-sign", type = "character", dest = "out_sign", default = NA,
              help = "optional: rewrite the pairwise-comparison TSV here"),
  make_option("--exclude", type = "character", default = "",
              help = "comma-separated method names to drop"),
  make_option("--recall", type = "double", default = NA,
              help = "the recall this table was built at; used for the axis label"),
  make_option("--methods-table", type = "character", dest = "methods_table", default = NA,
              help = "config_methods TSV; refreshes pred_name_long for the legend"),
  make_option("--p-threshold", type = "double", dest = "p_threshold", default = 0.05,
              help = "significance threshold for the pairwise table [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.na(opt$table) || is.na(opt$out)) stop("--table and --out are required")

df <- fread(opt$table, sep = "\t")
if (!"method" %in% colnames(df)) {
  stop("No 'method' column -- this table is the placeholder written when fewer than ",
       "two predictors reach the requested recall; nothing to replot.")
}

if (!is.na(opt$methods_table)) {
  mc <- fread(opt$methods_table, sep = "\t", select = c("method", "pred_name_long"))
  changed <- df %>% dplyr::select(method, old = pred_name_long) %>% distinct() %>%
    dplyr::inner_join(mc, by = "method") %>% dplyr::filter(old != pred_name_long)
  for (i in seq_len(nrow(changed))) {
    cat("Renamed:", changed$old[i], "->", changed$pred_name_long[i], "\n")
  }
  df <- dplyr::select(df, -pred_name_long) %>% dplyr::left_join(mc, by = "method")
}

exclude <- if (nzchar(opt$exclude)) trimws(strsplit(opt$exclude, ",")[[1]]) else character(0)
if (length(exclude) > 0) {
  absent <- setdiff(exclude, unique(df$method))
  if (length(absent) > 0) {
    stop("--exclude names not present in the table: ", paste(absent, collapse = ", "))
  }
  df <- dplyr::filter(df, !method %in% exclude)
  cat("Dropped", paste(exclude, collapse = ", "), "\n")
}
cat("Methods plotted:", paste(df$method, collapse = ", "), "\n")
if (nrow(df) < 2) stop("Fewer than two methods remain; a bar plot would be meaningless.")

# Same legend rule as the pipeline scripts.
df$key <- ifelse(df$Biosample == "all_matches",
                 df$pred_name_long,
                 paste0(df$pred_name_long, " (", df$Biosample, ")"))

# Pairwise comparison of relative risks, recomputed over the retained methods only.
df_p_val <- data.frame(t(combn(df$key, 2)))
colnames(df_p_val) <- c("group1", "group2")
df_p_val$p <- 0
for (i in seq_len(nrow(df_p_val))) {
  m1 <- dplyr::filter(df, key == df_p_val$group1[i])
  m2 <- dplyr::filter(df, key == df_p_val$group2[i])
  d <- log(m1$enrichment[1]) - log(m2$enrichment[1])
  SE_d <- sqrt(m1$SE_log_enr[1]**2 + m2$SE_log_enr[1]**2)
  df_p_val$p[i] <- pnorm(-(abs(d / SE_d))) * 2
}
df_p_val$p_adjust <- p.adjust(df_p_val$p, method = "bonferroni")
df_p_val$significant <- df_p_val$p_adjust < opt$p_threshold

df$recall.linking.rounded <- round(df$recall.linking, digits = 3)
df$plotting_label <- paste0(df$key, " (", df$recall.linking.rounded, ")")
df <- df[order(df$enrichment, decreasing = TRUE), ]
df$plotting_label <- factor(df$plotting_label, levels = unique(df$plotting_label), ordered = TRUE)
df$key <- factor(df$key, levels = unique(df$key), ordered = TRUE)

pred_colors <- df$hex
names(pred_colors) <- df$key
n_legend_cols <- ceiling(length(unique(df$key)) / 12)

recall_label <- if (is.na(opt$recall)) "" else paste0(" at recall ", opt$recall)

g <- ggplot(data = df, aes(x = key, y = enrichment)) +
  geom_col(aes(fill = key)) +
  geom_errorbar(aes(ymin = CI_enr_low, ymax = CI_enr_high), width = 0.2) +
  scale_fill_manual(values = pred_colors, labels = df$plotting_label) +
  labs(fill = "Predictor (exact recall)") +
  ylab(paste0("Enrichment (eQTLs versus common variants)", recall_label)) +
  theme_classic() + theme(axis.text = element_text(size = 7, color = "#000000"),
    axis.title = element_text(size = 8), axis.ticks = element_line(color = "#000000"),
    legend.text = element_text(size = 7), legend.title = element_text(size = 8)) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  guides(col = guide_legend(nrow = 12))

width <- 3 + 3 * n_legend_cols
ggsave(opt$out, g, width = width, height = 4)
cat("Wrote", opt$out, sprintf("(%.1f x 4.0 in)\n", width))

if (!is.na(opt$out_sign)) {
  write.table(df_p_val, opt$out_sign, quote = FALSE, col.names = TRUE,
              row.names = FALSE, sep = "\t")
  cat("Wrote", opt$out_sign, "\n")
}

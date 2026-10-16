## Enrichment-recall curves with enrichment divided by the enrichment of the model's own candidate
## universe -- every candidate element-gene pair the model scored, with no score cutoff (the
## candidate-universe heatmaps). 1 = no better than calling every candidate; the y-axis is the
## improvement the model's ranking adds over its universe. Use it to compare models whose
## universes enrich differently, and pairs whose universe enrichment differs.
##
## One plot per eQTL tissue and one for AllMatches. Per model:
##   thick line: pooled over the tissue's biosample pairings (all pairings for AllMatches). Counts
##     are summed over the pairings, as for the all_matches aggregate in
##     counts_to_enrichment_across_thresholds.R and compute_prediction_table_threshold.R; the
##     universe enrichment is pooled the same way over the same pairings.
##   faint lines: one per biosample-tissue pairing, divided by that pairing's universe enrichment.
##   dot with bar: the configured threshold, CI = enrichment CI / universe enrichment (the universe
##     enrichment is treated as fixed; pooled pairings are treated as independent).
## Filters as in plot_enrichment_recall.R: points with 20 or fewer eQTL pairs or recall.total of
## 0.001 or less are dropped (per pairing for faint lines, on the pooled counts for thick lines);
## curves drop recall.linking == 0; binary predictors are a point at threshold 1.
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
})

split_arg = function(x) unlist(strsplit(as.character(unlist(x)), " "))
method_names = split_arg(snakemake@params$methods)
score_thresholds = as.numeric(split_arg(snakemake@params$score_thresholds))  # negated for inverse predictors
universe_names = unlist(snakemake@params$universeNames)                      # one per method
recall_files = split_arg(snakemake@input$recallTables)
enrichment_files = split_arg(snakemake@input$enrichmentTables)
map_files = split_arg(snakemake@input$biosampleTissueMaps)
universe_files = split_arg(snakemake@input$candidateEnrichmentTables)        # one per method
z = qnorm(as.numeric(snakemake@params$thresholdPval) / 2, lower.tail = FALSE)
this_tissue = snakemake@wildcards$GTExTissue
out_plot = snakemake@output$plot
out_table = snakemake@output$table
stopifnot(length(recall_files) == length(method_names), length(enrichment_files) == length(method_names),
          length(map_files) == length(method_names), length(universe_files) == length(method_names),
          length(universe_names) == length(method_names))

cp = fread(snakemake@input$colorPalette, sep = "\t")
method_colors = setNames(cp$hex, cp$method)
method_labels = setNames(cp$pred_name_long, cp$method)

ratio_enrichment = function(x1, n1, x2, n2) (x1 / n1) / (x2 / n2)

curves = list()
for (i in seq_along(method_names)) {
  m = method_names[i]
  mp = fread(map_files[i], sep = "\t", colClasses = "character")
  if (this_tissue != "AllMatches") mp = mp[tissue == this_tissue]
  if (nrow(mp) == 0) next
  keys = unique(paste0(mp$tissue, "\t", mp$biosample))

  rec = fread(recall_files[i], sep = "\t", colClasses = list(character = "Biosample"))
  rec = rec[paste0(GTExTissue, "\t", Biosample) %in% keys,
            .(GTExTissue, Biosample, threshold, ov = nVariantGenePairsOverlappingEnhancers,
              lk = nVariantsOverlappingEnhancersCorrectGene, tv = total.variants)]
  enr = fread(enrichment_files[i], sep = "\t", colClasses = list(character = "Biosample"))
  enr = enr[paste0(GTExTissue, "\t", Biosample) %in% keys,
            .(GTExTissue, Biosample, threshold, x1 = nVariantsOverlappingEnhancers, n1 = nVariantsGTExTissue,
              x2 = nCommonVariantsOverlappingEnhancers, n2 = nCommonVariants)]
  cnt = merge(rec, enr, by = c("GTExTissue", "Biosample", "threshold"))
  if (nrow(cnt) != nrow(rec) || nrow(cnt) != nrow(enr))
    stop(m, ": recall and enrichment count tables do not line up for ", this_tissue)
  if (nrow(cnt) == 0) next

  uni = fread(universe_files[i], sep = "\t", colClasses = list(character = "Biosample"))
  uni = uni[paste0(GTExTissue, "\t", Biosample) %in% keys,
            .(GTExTissue, Biosample, u1 = nVariantsOverlappingEnhancers, un1 = nVariantsGTExTissue,
              u2 = nCommonVariantsOverlappingEnhancers, un2 = nCommonVariants)]
  missing = setdiff(unique(paste(cnt$GTExTissue, cnt$Biosample)), paste(uni$GTExTissue, uni$Biosample))
  if (length(missing) > 0) stop(m, ": no candidate-universe enrichment for ", paste(missing, collapse = ", "))
  cnt = merge(cnt, uni, by = c("GTExTissue", "Biosample"))
  binary = length(unique(cnt$threshold)) == 2

  # per pairing
  pair = cnt[, .(GTExTissue, Biosample, threshold, tv, rt = ov / tv, rl = lk / tv,
                 enrichment = ratio_enrichment(x1, n1, x2, n2),
                 universeEnrichment = ratio_enrichment(u1, un1, u2, un2))]
  pair[, level := "pairing"]

  # pooled over pairings
  pool = cnt[, .(ov = sum(ov), lk = sum(lk), tv = sum(tv), x1 = sum(x1), n1 = sum(n1), x2 = sum(x2), n2 = sum(n2)),
             by = threshold]
  u = unique(cnt[, .(GTExTissue, Biosample, u1, un1, u2, un2)])
  pooled_universe = ratio_enrichment(sum(u$u1), sum(u$un1), sum(u$u2), sum(u$un2))
  pool[, se := sqrt(((n1 - x1) / x1) / n1 + ((n2 - x2) / x2) / n2)]
  pool = pool[, .(GTExTissue = if (this_tissue == "AllMatches") "all_matches" else this_tissue,
                  Biosample = if (this_tissue == "AllMatches") "all_matches" else "pooled",
                  threshold, tv, rt = ov / tv, rl = lk / tv, enrichment = ratio_enrichment(x1, n1, x2, n2),
                  CI_low = exp(log(ratio_enrichment(x1, n1, x2, n2)) - z * se),
                  CI_high = exp(log(ratio_enrichment(x1, n1, x2, n2)) + z * se),
                  universeEnrichment = pooled_universe)]
  pool[, level := "pooled"]

  d = rbindlist(list(pool, pair), fill = TRUE)
  # As in plot_enrichment_recall.R: missing values (0/0) count as 0, then the small-count filter.
  for (col in c("rt", "rl", "enrichment", "CI_low", "CI_high")) set(d, which(is.na(d[[col]])), col, 0)
  d = d[tv > 20 & rt > 0.001]
  d[, `:=`(method = m, universe = universe_names[i], binary = binary,
           normalizedEnrichment = enrichment / universeEnrichment,
           normalizedCI_low = CI_low / universeEnrichment, normalizedCI_high = CI_high / universeEnrichment,
           isConfiguredThreshold = FALSE)]
  st = score_thresholds[i]
  if (!is.na(st) && !binary && any(d$level == "pooled")) {
    pl = which(d$level == "pooled")
    d$isConfiguredThreshold[pl[which.min(abs(d$threshold[pl] - st))]] = TRUE
  }
  curves[[m]] = d
}
d = rbindlist(curves, fill = TRUE)

out_cols = c("method", "universe", "level", "GTExTissue", "Biosample", "threshold", "tv", "rl", "rt",
             "enrichment", "universeEnrichment", "normalizedEnrichment", "normalizedCI_low", "normalizedCI_high",
             "isConfiguredThreshold")
if (nrow(d) == 0) {
  msg = paste0("No ", this_tissue, " pairing has more than 20 eQTL pairs; nothing to plot.")
  ggsave(out_plot, ggplot() + annotate("text", x = 0, y = 0, label = msg) + theme_void(), width = 6, height = 2)
  fwrite(data.table(matrix(ncol = length(out_cols), nrow = 0, dimnames = list(NULL, out_cols))), out_table, sep = "\t")
  quit(save = "no", status = 0)
}
tab = d[, ..out_cols]
setnames(tab, c("tv", "rl", "rt"), c("total.variants", "recall.linking", "recall.total"))
fwrite(tab[order(method, level, GTExTissue, Biosample, threshold)], out_table, sep = "\t")

models = names(method_colors)[names(method_colors) %in% d$method]
plot_d = d[rl > 0]
lines_pooled = plot_d[level == "pooled" & !binary][order(method, rl, threshold)]
lines_pair = plot_d[level == "pairing" & !binary][order(method, GTExTissue, Biosample, rl, threshold)]
points_pooled = plot_d[level == "pooled" & binary & threshold == 1]
points_pair = plot_d[level == "pairing" & binary & threshold == 1]
markers = d[isConfiguredThreshold == TRUE]

# The y range follows the pooled curves; faint per-pairing lines that spike at strict thresholds are clipped.
ylim = min(20, max(c(1.5, lines_pooled$normalizedEnrichment, points_pooled$normalizedEnrichment), na.rm = TRUE) * 1.1)
xlim = min(1, max(c(0.01, plot_d$rl), na.rm = TRUE) * 1.05)
n_pairings = length(unique(paste(d$GTExTissue[d$level == "pairing"], d$Biosample[d$level == "pairing"])))
x_label = paste0("Recall (eQTLs overlapping prediction linked to eGene)\n",
                 if (this_tissue == "AllMatches") paste0("all ", n_pairings, " biosample-tissue pairings")
                 else paste0(this_tissue, ", ", n_pairings, " biosample pairing", ifelse(n_pairings == 1, "", "s")))
# universeNames already name every model sharing the universe
caption = paste0("Thick: pooled over pairings. Faint: one per pairing, divided by its own universe enrichment.\n",
                 "Universe = the model's candidate element-gene pairs with no score cutoff, one per group of models sharing it:\n",
                 paste0("  ", unique(d$universe), collapse = "\n"))

g = ggplot() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_line(data = lines_pair, aes(rl, normalizedEnrichment, color = method, group = interaction(method, GTExTissue, Biosample)),
            linewidth = 0.35, alpha = 0.25) +
  geom_line(data = lines_pooled, aes(rl, normalizedEnrichment, color = method), linewidth = 1.3)
if (nrow(points_pair) > 0)
  g = g + geom_point(data = points_pair, aes(rl, normalizedEnrichment, color = method), size = 1.5, alpha = 0.35, show.legend = FALSE)
if (nrow(points_pooled) > 0)
  g = g + geom_point(data = points_pooled, aes(rl, normalizedEnrichment, color = method), size = 3, shape = 16)
if (nrow(markers) > 0)
  g = g + geom_linerange(data = markers, aes(x = rl, ymin = normalizedCI_low, ymax = normalizedCI_high, color = method),
                         linewidth = 0.75, show.legend = FALSE) +
          geom_point(data = markers, aes(rl, normalizedEnrichment, color = method), size = 3, shape = 16, show.legend = FALSE)
g = g +
  scale_color_manual(values = method_colors[models], labels = method_labels[models], breaks = models) +
  labs(x = x_label, y = "Enrichment / candidate-universe enrichment\n(1 = calling every candidate)",
       col = "Model", title = this_tissue, caption = caption) +
  coord_cartesian(xlim = c(0, xlim), ylim = c(0, ylim)) +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 10, color = "#000000"), axis.title = element_text(size = 11),
        axis.ticks = element_line(color = "#000000"),
        panel.grid.major = element_line(color = "grey92", linewidth = 0.3), panel.grid.minor = element_blank(),
        legend.text = element_text(size = 9), legend.title = element_text(size = 10),
        plot.caption = element_text(size = 7, color = "#555555", hjust = 0), plot.caption.position = "plot",
        aspect.ratio = 1) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.5, alpha = 1)))
ggsave(out_plot, g, width = 8, height = 7)

## NEW RECALL
## load libraries
library(dplyr)
library(ggplot2)
library(colorspace)
library(data.table)
library(tidyr)
library(gridExtra)
library(grid)

## INPUTS

files = snakemake@input$enrichmentRecall_files %>% strsplit(" ") %>% unlist()
score_thresholds = snakemake@params$score_thresholds %>% as.character() %>% strsplit(" ") %>% unlist() %>% as.numeric()
method_names = snakemake@params$methods %>% strsplit(" ") %>% unlist()
this_tissue = snakemake@wildcards$GTExTissue
out_plot = snakemake@output$er_combined
out_table = snakemake@output$er_combined_table
cpFilePlotting = snakemake@input$colorPalette
cp = fread(cpFilePlotting, sep="\t") # method, pred_name_long, hex

# Pristine method-level color/label lookup, captured before `cp` gets
# overwritten below for the legacy per-(method,biosample)-key table logic.
# The redesigned plot legend below always keys off THIS -- one entry per
# method, never per biosample -- regardless of how many biosamples this
# tissue matches.
method_colors = setNames(cp$hex, cp$method)
method_labels = setNames(cp$pred_name_long, cp$method)

# gather tables
for (i in 1:length(files)){
	temp = fread(files[i], sep="\t")
	temp$nPoints = nrow(temp)

	if (i==1){df = temp} else {
		df = rbind(df, temp)
	}
}

# organize for plotting
scores = data.frame(method=method_names, score_threshold=score_thresholds)
cp = left_join(cp, scores, by="method")

df = dplyr::left_join(df, cp, by="method")
# The " (biosample)" suffix exists to disambiguate a predictor that matches one
# eQTL tissue through more than one prediction biosample. In the AllMatches
# aggregate every row carries the synthetic biosample "all_matches", so the
# suffix adds nothing and appears in every legend entry as noise.
df$key = ifelse(df$Biosample == "all_matches",
                df$pred_name_long,
                paste0(df$pred_name_long, " (", df$Biosample, ")"))
df[is.na(df)] = 0
df = df[order(df$threshold),]

df <- dplyr::filter(df, total.variants > 20, recall.total > 0.001)

# A tissue whose pairings all have too few eQTLs leaves nothing to draw. Write a
# placeholder that says so instead of failing on empty data further down.
if (nrow(df) == 0) {
  msg = paste0("No ", this_tissue, " pairing has more than 20 eQTL pairs; nothing to plot.")
  ggsave(out_plot, ggplot() + annotate("text", x = 0, y = 0, label = msg) + theme_void(), width = 6, height = 2)
  write.table(df, out_table, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
  quit(save = "no", status = 0)
}

# handle multiple matches (kept for out_table's own key/hex columns -- these
# per-biosample-key colors are no longer used by the plot itself, which now
# always colors by method via method_colors/method_labels above)
cp = dplyr::select(df, method, key, hex) %>% distinct()
methods = unique(cp$method)
if (length(methods) < length(unique(cp$key))) {
	# set new hex codes
	for (i in 1:length(methods)){
		cp.this =  dplyr::filter(cp, method==methods[i])
		n = nrow(cp.this)
		if (n>1){
			hex_base = cp.this$hex[1]
			cols = colorRampPalette(c(hex_base, '#000000'))(n+2) # range of colors from hex_base to black (where black is 2 above the base)
			for (k in 1:nrow(cp.this)){
				cp$hex[cp$key==cp.this$key[k]] = cols[k]  # set corresponding hex code
			}
		}
	}
}

df = dplyr::select(df, -hex) %>% left_join(cp, by=c("key", "method")) # for saving later!

# data for binary predictors to be plotted as points
df$key = factor(df$key)
df_plot = dplyr::filter(df, recall.linking>0) # remove (0,0) points

df_binary = dplyr::filter(df_plot, nPoints==2) %>% dplyr::filter(threshold==1)
# df[is.na(df)] = 0 above turned a missing threshold into 0, so pick marker rows
# only for methods that have a threshold.
thresholded_methods = scores$method[!is.na(scores$score_threshold)]
df_thresh = dplyr::filter(df, method %in% thresholded_methods, threshold==score_threshold) # suggested thresholds
df_plot = dplyr::filter(df_plot, nPoints>2)

# plotting params
ylim = min(50, max(df_plot$enrichment))

write.table(df, out_table, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

if (this_tissue == "AllMatches") {

  ## AllMatches is already one pre-pooled row-set per method (Biosample ==
  ## "all_matches") -- there is no per-biosample clutter to fix here, so this
  ## path is unchanged from the original design other than coloring by the
  ## pristine method_colors/method_labels captured above instead of `key`.
  n_var_label = paste0(mean(df_plot$total.variants), " variant-biosample pairs")
  x_label = paste0("Recall (variants overlapping prediction linked to eGene)\n", n_var_label)

  g = ggplot(data=df_plot, aes(x=recall.linking, y=enrichment, color=method)) +
    geom_line(linewidth=0.75) +
    geom_point(data=df_binary, aes(x=recall.linking, y=enrichment, color=method), size=3, shape = 16) +
    geom_linerange(data=df_thresh, aes(ymin=CI_enr_low, ymax=CI_enr_high, color=method), linewidth = 0.75) +
    geom_point(data=df_thresh, aes(x=recall.linking, y=enrichment, color=method), size=3, shape = 16) +
    scale_color_manual(values=method_colors, labels=method_labels, breaks=names(method_colors)) +
    ylab("Enrichment (eQTLs versus common variants)") + xlab(x_label) +
    labs(col="Predictor") +
    coord_cartesian(ylim=c(0,ylim)) +
    theme_classic() + theme(axis.text = element_text(size = 7, color = "#000000"),
    	axis.title = element_text(size = 8), axis.ticks = element_line(color = "#000000"),
    	legend.text = element_text(size= 7), legend.title=element_text(size=8), legend.position="right", legend.direction="vertical",
    	aspect.ratio = 1)

  ggsave(out_plot, g, width=6, height=4)

} else {

  ## Per-tissue view. The input files hold one row-set per (method, biosample)
  ## match, so a tissue paired with many biosamples used to get one legend entry
  ## per pair. Instead: each biosample-method curve is thin and faint, the
  ## biosamples are pooled into one thick curve per method, and the legend has
  ## one entry per method. Pooling uses the same sum-of-counts convention as the
  ## all_matches aggregate in counts_to_enrichment_across_thresholds.R.

  models_here = names(method_colors)[names(method_colors) %in% unique(df$method)]
  binary_models = intersect(models_here, unique(df$method[df$nPoints == 2]))
  curve_models = setdiff(models_here, binary_models)

  counts_files = snakemake@input$enrichmentTable_acrossThresholds %>% strsplit(" ") %>% unlist()
  agg_list = list()
  for (i in seq_along(method_names)) {
    m = method_names[i]
    if (!(m %in% models_here)) next
    biosamples_m = unique(df$Biosample[df$method == m])
    counts = fread(counts_files[i], sep = "\t")
    counts = dplyr::filter(counts, GTExTissue == this_tissue, Biosample %in% biosamples_m)

    # As with all_matches, this CI treats the pooled pairs as independent; they
    # are not (README: Known limitations).
    # n1 and n2 are summed along with x1 and x2, as in the all_matches aggregate.
    # Pairing one biosample's n1 with a summed x1 lets x1 exceed n1, and the SE
    # then goes NaN.
    pooled = counts %>%
      group_by(threshold) %>%
      summarize(
        x1 = sum(nVariantsOverlappingEnhancers),
        n1 = sum(nVariantsGTExTissue),
        x2 = sum(nCommonVariantsOverlappingEnhancers),
        n2 = sum(nCommonVariants),
        .groups = "drop"
      ) %>%
      mutate(
        enrichment = (x1 / n1) / (x2 / n2),
        SE_log_enr = sqrt((n1 - x1) / x1 / n1 + (n2 - x2) / x2 / n2),
        CI_enr_low  = exp(log(enrichment) - 1.96 * SE_log_enr),
        CI_enr_high = exp(log(enrichment) + 1.96 * SE_log_enr)
      )

    recall_by_threshold = dplyr::filter(df, method == m) %>%
      group_by(threshold) %>%
      summarize(recall.linking = mean(recall.linking), .groups = "drop")

    agg_list[[m]] = dplyr::inner_join(pooled, recall_by_threshold, by = "threshold") %>%
      dplyr::filter(is.finite(enrichment), recall.linking > 0) %>%
      dplyr::arrange(recall.linking) %>%
      dplyr::mutate(method = m)
  }
  agg = dplyr::bind_rows(agg_list)
  agg_lines  = dplyr::filter(agg, method %in% curve_models)
  agg_points = dplyr::filter(agg, method %in% binary_models, threshold == 1)

  # One shared scale for every panel, so the panels can be compared by eye.
  ylim = min(50, max(c(df_plot$enrichment, df_binary$enrichment, agg$enrichment, 1), na.rm = TRUE))
  xlim = min(1, max(c(df_plot$recall.linking, df_binary$recall.linking, agg$recall.linking, 0.01), na.rm = TRUE) * 1.05)

  # Recall counts eVariant-eGene PAIRS (total.variants); enrichment counts unique
  # eVariants. Each axis states its own denominator. The variants file is
  # headerless: column 4 is the variant ID and column 6 the eQTL tissue.
  filtered_variants = fread(snakemake@input$variantsFile, sep = "\t", header = FALSE,
                            select = c(4, 6), col.names = c("varID", "tissue"))
  n_unique_variants = length(unique(filtered_variants[tissue == this_tissue, varID]))
  n_variant_gene_pairs = round(mean(df$total.variants))

  x_label = paste0("Recall (eQTLs overlapping prediction linked to eGene)\n",
                   n_variant_gene_pairs, " ", this_tissue, " eQTLs (eVariant-eGene pairs)")
  y_label = paste0("Enrichment (eVariants versus common SNPs)\n",
                   n_unique_variants, " unique ", this_tissue, " eVariants")

  # Mark each pooled curve at the method's configured threshold. score_thresholds
  # is already sign-flipped for inverse predictors (see Snakefile).
  marker_list = list()
  for (m in curve_models) {
    agg_m = dplyr::filter(agg_lines, method == m)
    if (nrow(agg_m) == 0) next
    signed_threshold = scores$score_threshold[scores$method == m]
    if (is.na(signed_threshold)) next  # no threshold configured: no marker
    marker_list[[m]] = agg_m[which.min(abs(agg_m$threshold - signed_threshold)), ]
  }
  markers = dplyr::bind_rows(marker_list)

  build_plot = function(methods_to_draw, title, subtitle, show_legend = TRUE) {
    df_plot_sub = dplyr::filter(df_plot, method %in% methods_to_draw)
    lines_sub   = dplyr::filter(agg_lines, method %in% methods_to_draw)
    points_sub  = dplyr::filter(agg_points, method %in% methods_to_draw)
    binary_sub  = dplyr::filter(df_binary, method %in% methods_to_draw)
    markers_sub = if (nrow(markers) > 0) dplyr::filter(markers, method %in% methods_to_draw) else markers

    p = ggplot() +
      geom_line(data = df_plot_sub, aes(x = recall.linking, y = enrichment, color = method, group = interaction(method, Biosample)),
                linewidth = 0.35, alpha = 0.25) +
      geom_line(data = lines_sub, aes(x = recall.linking, y = enrichment, color = method),
                linewidth = 1.3, alpha = 1)
    # Binary predictors have one operating point, so they are drawn as points:
    # faint per biosample, solid for the pooled value.
    if (nrow(binary_sub) > 0) {
      p = p + geom_point(data = binary_sub, aes(x = recall.linking, y = enrichment, color = method),
                         size = 1.5, alpha = 0.35, show.legend = FALSE)
    }
    if (nrow(points_sub) > 0) {
      p = p + geom_point(data = points_sub, aes(x = recall.linking, y = enrichment, color = method), size = 3, shape = 16)
    }
    p = p +
      scale_color_manual(values = method_colors[models_here], labels = method_labels[models_here], breaks = models_here) +
      ylab(y_label) + xlab(x_label) +
      labs(col = "Model", title = title, subtitle = subtitle) +
      coord_cartesian(xlim = c(0, xlim), ylim = c(0, ylim)) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 11, color = "#333333"),
        axis.text = element_text(size = 10, color = "#000000"),
        # The rotated y title is two lines long; at this size it still fits
        # inside one panel's height (see panel heights below).
        axis.title = element_text(size = 12),
        axis.ticks = element_line(color = "#000000"),
        panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.position = if (show_legend) "right" else "none",
        legend.direction = "vertical",
        plot.margin = margin(t = 12, r = 12, b = 12, l = 12),
        aspect.ratio = 1
      ) +
      guides(color = guide_legend(override.aes = list(linewidth = 1.5, alpha = 1)))

    if (nrow(markers_sub) > 0) {
      # show.legend = FALSE: the line already gives each model a legend key, and
      # stacking a vertical segment on it renders the key as a "+".
      p = p +
        geom_linerange(data = markers_sub, aes(x = recall.linking, ymin = CI_enr_low, ymax = CI_enr_high, color = method),
                       linewidth = 0.75, show.legend = FALSE) +
        geom_point(data = markers_sub, aes(x = recall.linking, y = enrichment, color = method),
                   size = 3, shape = 16, show.legend = FALSE)
    }
    p
  }

  # ---- pairing table: every biosample the pipeline paired with this tissue ----
  map_files = snakemake@input$biosampleTissueMaps %>% strsplit(" ") %>% unlist()
  pairs = rbindlist(lapply(seq_along(map_files), function(i) {
    mp = fread(map_files[i], sep = "\t", colClasses = "character")
    mp$method = method_names[i]
    mp
  }), fill = TRUE)
  pairs = pairs[tissue == this_tissue]
  paired_biosamples = sort(unique(pairs$biosample))
  preds = fread(snakemake@params$predictionsTable, sep = "\t", colClasses = "character")
  preds$biosample = trimws(preds$biosample)

  # A readable name for a biosample: the predictions table's SampleTermName column if
  # it has one, else the "# SampleTermName:" line of one of its prediction files'
  # headers, else the biosample ID (the table and the headers are both optional).
  sample_term_name = function(b) {
    row = preds[biosample == b]
    if (nrow(row) == 0) return(b)
    if ("SampleTermName" %in% colnames(row) && nzchar(row$SampleTermName[1])) return(row$SampleTermName[1])
    for (m in intersect(method_names, colnames(row))) {
      f = row[[m]][1]
      if (is.na(f) || !nzchar(f) || !file.exists(f)) next
      con = if (grepl("\\.gz$", f)) gzfile(f) else file(f)
      head_lines = tryCatch(readLines(con, n = 40, warn = FALSE), error = function(e) character(0))
      close(con)
      hit = grep("^#\\s*SampleTermName\\s*:", head_lines, value = TRUE)
      if (length(hit) > 0) return(trimws(sub("^#\\s*SampleTermName\\s*:", "", hit[1])))
    }
    b
  }

  # A tissue paired with exactly one biosample gets one panel, subtitled with that
  # biosample's name: per-model panels would only repeat the same single pairing.
  n_paired = length(paired_biosamples)
  subtitle = if (n_paired == 1) {
    paste0("All models -- ", sample_term_name(paired_biosamples))
  } else {
    paste0("All models -- ", n_paired, " biosample pairings")
  }
  g = build_plot(models_here, this_tissue, subtitle)

  # Per-model panels add nothing when there is only one model or only one paired
  # biosample. They carry no legend: the one model they show is named in the title.
  sub_plots = list()
  if (length(models_here) > 1 && n_paired > 1) {
    sub_plots = lapply(models_here, function(m) {
      n_bio_m = length(unique(dplyr::filter(df, method == m)$Biosample))
      build_plot(m, paste0(this_tissue, " -- ", method_labels[[m]]),
                 paste0(n_bio_m, " biosample", ifelse(n_bio_m == 1, "", "s")),
                 show_legend = FALSE)
    })
  }

  pooled_keys = unique(paste(df$method, df$Biosample))
  pairs$pooled = paste(pairs$method, pairs$biosample) %in% pooled_keys

  table_df = pairs %>%
    group_by(biosample) %>%
    summarize(`Models pooled` = paste0(sum(pooled), "/", n()), .groups = "drop") %>%
    dplyr::rename(Biosample = biosample) %>%
    as.data.frame()

  # Optional annotation columns, shown only if the predictions table has them.
  ann_cols = intersect(c("CellAnnotation", "SampleTermID"), colnames(preds))
  if (length(ann_cols) > 0) {
    ann = unique(preds[, c("biosample", ann_cols), with = FALSE])
    ann$biosample = trimws(ann$biosample)
    ann = ann[!duplicated(ann$biosample)]
    table_df = dplyr::left_join(table_df, ann, by = c("Biosample" = "biosample"))
    if ("SampleTermID" %in% ann_cols) table_df = dplyr::rename(table_df, `Ontology ID` = SampleTermID)
  }

  # Candidate-universe enrichment (no score cutoff), one column per distinct
  # candidate universe. A different quantity from the recall-based curves.
  wrap_col = function(x, width) {
    vapply(x, function(s) if (is.na(s)) "" else paste(strwrap(s, width = width), collapse = "\n"), character(1), USE.NAMES = FALSE)
  }
  cand_files = snakemake@input$candidateEnrichmentTables %>% strsplit(" ") %>% unlist()
  cand_names = unlist(snakemake@params$candidateUniverseNames)
  sort_key = NULL
  for (k in seq_along(cand_files)) {
    cand = fread(cand_files[k], sep = "\t", colClasses = list(character = "Biosample"))
    cand = dplyr::filter(cand, GTExTissue == this_tissue)
    hit = cand[match(table_df$Biosample, cand$Biosample), ]
    cell = ifelse(is.na(hit$enrichment), "",
                  paste0(sprintf("%.2f", hit$enrichment), "\np=",
                         ifelse(hit$p_adjust_enr == 0, "0", formatC(hit$p_adjust_enr, format = "e", digits = 1))))
    table_df[[paste0("Cand. enr.\n", wrap_col(cand_names[k], 16))]] = cell
    if (is.null(sort_key)) sort_key = hit$enrichment
  }
  if (!is.null(sort_key)) table_df = table_df[order(-sort_key, na.last = TRUE), , drop = FALSE]

  # Biosample IDs have no spaces, so make "_"/"-" breakable before wrapping.
  wrap_identifier = function(x, width) {
    gsub("([_-]) ", "\\1", wrap_col(gsub("([_-])", "\\1 ", x), width))
  }
  table_df$Biosample = wrap_identifier(table_df$Biosample, 16)
  if ("CellAnnotation" %in% colnames(table_df)) table_df$CellAnnotation = wrap_col(table_df$CellAnnotation, 34)
  table_df[is.na(table_df)] = ""

  count_lines = function(x) vapply(strsplit(as.character(x), "\n", fixed = TRUE), function(v) max(length(v), 1L), integer(1))
  n_lines = do.call(pmax, lapply(table_df, count_lines))
  header_lines = max(count_lines(colnames(table_df)))

  col_width_in = function(name) {
    if (name == "Biosample") return(1.7)
    if (name == "CellAnnotation") return(2.7)
    if (name == "Ontology ID") return(1.1)
    if (name == "Models pooled") return(0.9)
    1.1
  }
  col_widths_in = vapply(colnames(table_df), col_width_in, numeric(1))

  tbl = tableGrob(table_df, rows = NULL,
                  theme = ttheme_minimal(base_size = 7,
                                         core = list(fg_params = list(hjust = 0, x = 0.02)),
                                         colhead = list(fg_params = list(hjust = 0, x = 0.02, fontface = "bold"))))
  tbl$widths = unit(col_widths_in, "in")
  row_line_height_in = 0.14
  header_height_in = header_lines * row_line_height_in + 0.11
  row_heights_in = pmax(n_lines * row_line_height_in + 0.08, 0.22)
  tbl$heights = unit(c(header_height_in, row_heights_in), "in")
  table_height_in = header_height_in + sum(row_heights_in)

  caption = textGrob(
    paste0("\"Cand. enr.\": enrichment over each candidate element-gene universe with no score cutoff (as in the\n",
           "candidate-universe heatmaps), with adjusted p-value. \"Models pooled\": models whose curve includes this\n",
           "biosample, out of models paired with it; pairings with 20 or fewer eQTL pairs are not pooled."),
    x = 0.01, hjust = 0, gp = gpar(fontsize = 6, col = "#666666")
  )

  # Every panel is a fixed 6in tall so the two-line rotated y title fits inside
  # it; the figure grows with the number of panels and table rows instead.
  out_width = max(sum(col_widths_in), 7)
  panel_height_in = 6
  gap_in = 0.4

  panel_grobs = c(list(g), sub_plots)

  # A blank spacer after each panel, so one panel's x title does not run into
  # the next panel's title.
  interleaved_grobs = list()
  interleaved_heights_in = c()
  for (i in seq_along(panel_grobs)) {
    interleaved_grobs = c(interleaved_grobs, list(panel_grobs[[i]], grid::nullGrob()))
    interleaved_heights_in = c(interleaved_heights_in, panel_height_in, gap_in)
  }

  caption_height_in = 0.5
  combined = arrangeGrob(grobs = c(interleaved_grobs, list(tbl, caption)), ncol = 1,
                         widths = unit(out_width, "in"),
                         heights = unit(c(interleaved_heights_in, table_height_in, caption_height_in), "in"))

  out_height = sum(interleaved_heights_in) + table_height_in + caption_height_in
  ggsave(out_plot, combined, width = out_width, height = out_height, limitsize = FALSE)
}

## Interactive (plotly) version of a per-tissue enrichment-recall plot. Reads
## the plot data saved by plot_enrichment_recall.R, so the pooling and CI math
## live in one place. Hovering a faint line names its biosample.
suppressPackageStartupMessages({
  library(dplyr)
  library(plotly)
  library(htmlwidgets)
})

d = readRDS(snakemake@input$plot_data)
out_html = snakemake@output$html
stopifnot(length(d) > 0)  # only built for tissues with per-biosample curves

fmt = function(x) formatC(x, format = "g", digits = 4)
br = function(x) gsub("\n", "<br>", x)
bio_name = function(b) {
  n = unname(d$biosample_names[b])
  if (is.null(n) || is.na(n)) b else n
}
# thresholds are stored negated for inverse predictors; show the real cutoff
cutoff_text = function(m, t) if (isTRUE(d$inverse[[m]])) paste0("score <= ", fmt(-t)) else paste0("score >= ", fmt(t))

p = plot_ly()
for (m in d$models) {
  col = d$colors[[m]]
  lab = d$labels[[m]]

  # one faint trace per biosample, grouped with its model in the legend
  cur = dplyr::filter(d$curves, method == m)
  for (b in unique(cur$Biosample)) {
    cb = dplyr::filter(cur, Biosample == b) %>% dplyr::arrange(recall.linking)
    hover = paste0("<b>", lab, "</b><br>", bio_name(b), "<br><i>", b, "</i>",
                   "<br>", cutoff_text(m, cb$threshold),
                   "<br>recall: ", fmt(cb$recall.linking), "<br>enrichment: ", fmt(cb$enrichment))
    p = add_trace(p, x = cb$recall.linking, y = cb$enrichment, type = "scatter", mode = "lines",
                  line = list(color = col, width = 1.2), opacity = 0.35, name = bio_name(b),
                  legendgroup = m, showlegend = FALSE, text = hover, hoverinfo = "text", meta = "biosample")
  }
  bp = dplyr::filter(d$binary_points, method == m)
  if (nrow(bp) > 0) {
    hover = paste0("<b>", lab, "</b><br>", vapply(bp$Biosample, bio_name, character(1)), "<br><i>", bp$Biosample, "</i>",
                   "<br>recall: ", fmt(bp$recall.linking), "<br>enrichment: ", fmt(bp$enrichment))
    p = add_trace(p, x = bp$recall.linking, y = bp$enrichment, type = "scatter", mode = "markers",
                  marker = list(color = col, size = 7), opacity = 0.35, legendgroup = m, showlegend = FALSE,
                  text = hover, hoverinfo = "text")
  }

  # pooled curve (or point, for a binary model): the one legend entry per model
  n_bio = length(unique(c(cur$Biosample, bp$Biosample)))
  pooled_label = paste0("pooled across ", n_bio, " biosample", ifelse(n_bio == 1, "", "s"))
  pl = dplyr::filter(d$pooled_lines, method == m)
  pp = dplyr::filter(d$pooled_points, method == m)
  if (nrow(pl) > 0) {
    hover = paste0("<b>", lab, "</b><br>", pooled_label, "<br>", cutoff_text(m, pl$threshold),
                   "<br>recall: ", fmt(pl$recall.linking), "<br>enrichment: ", fmt(pl$enrichment))
    p = add_trace(p, x = pl$recall.linking, y = pl$enrichment, type = "scatter", mode = "lines",
                  line = list(color = col, width = 3.5), name = lab, legendgroup = m, showlegend = TRUE,
                  text = hover, hoverinfo = "text")
  } else if (nrow(pp) > 0) {
    hover = paste0("<b>", lab, "</b><br>", pooled_label,
                   "<br>recall: ", fmt(pp$recall.linking), "<br>enrichment: ", fmt(pp$enrichment))
    p = add_trace(p, x = pp$recall.linking, y = pp$enrichment, type = "scatter", mode = "markers",
                  marker = list(color = col, size = 11), name = lab, legendgroup = m, showlegend = TRUE,
                  text = hover, hoverinfo = "text")
  }

  # threshold marker with its 95% CI
  mk = if (nrow(d$markers) > 0) dplyr::filter(d$markers, method == m) else d$markers
  if (nrow(mk) > 0) {
    hover = paste0("<b>", lab, "</b><br>configured threshold (", cutoff_text(m, mk$threshold), ")",
                   "<br>recall: ", fmt(mk$recall.linking), "<br>enrichment: ", fmt(mk$enrichment),
                   " (95% CI ", fmt(mk$CI_enr_low), "-", fmt(mk$CI_enr_high), ")")
    p = add_trace(p, x = mk$recall.linking, y = mk$enrichment, type = "scatter", mode = "markers",
                  marker = list(color = col, size = 11),
                  error_y = list(type = "data", symmetric = FALSE, color = col,
                                 array = mk$CI_enr_high - mk$enrichment, arrayminus = mk$enrichment - mk$CI_enr_low),
                  legendgroup = m, showlegend = FALSE, text = hover, hoverinfo = "text")
  }
}

p = layout(p,
  title = list(text = paste0("<b>", d$tissue, "</b><br><sup>hover a faint line to see its biosample; ",
                             "click a model in the legend to hide it, double-click to show only that model</sup>")),
  xaxis = list(title = br(d$x_label), range = c(0, d$xlim), zeroline = FALSE, gridcolor = "#ebebeb"),
  yaxis = list(title = br(d$y_label), range = c(0, d$ylim), zeroline = FALSE, gridcolor = "#ebebeb"),
  legend = list(title = list(text = "<b>Model</b>")),
  hovermode = "closest", plot_bgcolor = "white", margin = list(t = 90))

# Brighten the hovered biosample line so it stands out from the others.
p = htmlwidgets::onRender(p, "
function(el) {
  var restyle = function(e, opacity, width) {
    var i = e.points[0].curveNumber;
    var meta = el.data[i].meta;  // plotly R recycles meta to one value per point
    if ((Array.isArray(meta) ? meta[0] : meta) !== 'biosample') return;
    Plotly.restyle(el, {opacity: opacity, 'line.width': width}, [i]);
  };
  el.on('plotly_hover', function(e) { restyle(e, 1, 3); });
  el.on('plotly_unhover', function(e) { restyle(e, 0.35, 1.2); });
}")

htmlwidgets::saveWidget(p, out_html, selfcontained = TRUE, title = paste0("Enrichment-recall: ", d$tissue))
# saveWidget can leave its temporary <name>_files folder next to the HTML
unlink(sub("\\.html$", "_files", out_html), recursive = TRUE)

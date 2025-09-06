## Reformat pillar project predictors
## April 24, 2025

library(optparse)

# Process input arguments --------------------------------------------------------------------------

# create arguments list
option_list = list(
  make_option(c("-i", "--input_file"), type = "character", default = NULL,
              help = "Path to transcripts input file", metavar = "character"),
  make_option(c("-o", "--output_file"), type = "character", default = NULL,
              help = "Path to output file", metavar = "character"),
  # make_option(c("-c", "--cell_type"), type = "character", default = NULL,
  #             help = "Cell type", metavar = "character"),
  # make_option(c("-d", "--term_id"), type = "character", default = NULL,
  #             help = "IGVF Sample Term ID", metavar = "character"),
  # make_option(c("-s", "--summary"), type = "character", default = NULL,
  #             help = "Short description of the sample including treatments", metavar = "character"),
  # make_option(c("-m", "--method"), type = "character", default = "scE2G",
  #             help = "E2G method that produced the predictions", metavar = "character"),
  # make_option(c("-v", "--version"), type = "character", default = NULL,
  #             help = "E2G method version", metavar = "character"),
  make_option(c("--threshold"), type = "character", default = NULL,
              help = "Used score threshold if applicable", metavar = "character"), 
  make_option(c("-l", "--portal_link"), type = "character", default = NULL,
              help = "Link to metadata of this file on IGVF data portal", metavar = "character"), 
  make_option(c("-a", "--all_columns"), action = "store_true", default = FALSE,
              help = "Include all columns in element/gene lists")
)

# parse arguments
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# function to check for required arguments
check_required_args <- function(arg, opt, opt_parser) {
  if (is.null(opt[[arg]])) {
    print_help(opt_parser)
    stop(arg, " argument is required!", call. = FALSE)
  }
}

# check that all required parameters are provided
required_args <- c("input_file", "output_file")
for (i in required_args) {
  check_required_args(i, opt = opt, opt_parser = opt_parser)
}

# Process file -------------------------------------------------------------------------------------

# required packages
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# load input file
pred <- fread(opt$input_file)

# get all score columns (all columns except EG-pair defining columns)
message("Reformatting predictions...")

# set summarized sample description if specified
if (!is.null(opt$summary)) {
  pred$SampleSummaryShort <- opt$summary
}

# create header lines
header <- c(
  paste("# Source:", opt$method),
  # paste("# Version:", opt$version),
  "# GenomeBuild: GRCh38",
  "# URL: https://github.com/EngreitzLab/scE2G/tree/main",
  "# Assays: 10x Multiome",
  "# SampleAgnostic: False"

  # paste("# SampleTermName:", unique(pred$CellType)),
  # paste("# SampleTermID:", opt$term_id),
  # paste("# SampleSummaryShort:", opt$summary)
)

# add threshold if applicable
if (!is.null(opt$threshold)) {
  header <- c(header, paste("# ScoreThreshold:", opt$threshold))
}

# add ScoreType if these are predictions (not metadata List)
if (!grepl("_list", opt$input_file)) {
  header <- c(header, "# ScoreType: positive_score")
}

# add link to where the metadata of the file is stored on the IGVF data portal if available
if (!is.null(opt$portal_link)) {  
    header <- c(header, "# Metadata:", opt$portal_link)
}

# add additional columns and extract output columns
if (opt$all_columns) { # include all columns
  pred <- pred %>%
    mutate(
          name = paste0(chr, ":", start, "-", end)) %>%
    select(ElementChr = chr,
          ElementStart = start,
          ElementEnd = end,
          ElementName = name,
          ElementClass = class,
          GeneSymbol = TargetGene,
          GeneEnsemblID = TargetGeneEnsembl_ID,
          GeneTSS = TargetGeneTSS,
          CellType = CellType,
          # SampleSummaryShort,
          isSelfPromoter = isSelfPromoter,
          ABC.Score,
          distance,
          Kendall,
          everything())
} else { # E2G predictions
    pred <- pred %>%
    mutate(
          name = paste0(chr, ":", start, "-", end)) %>%
    select(ElementChr = chr,
          ElementStart = start,
          ElementEnd = end,
          ElementName = name,
          ElementClass = class,
          GeneSymbol = TargetGene,
          GeneEnsemblID = TargetGeneEnsembl_ID,
          GeneTSS = TargetGeneTSS,
          CellType = CellType,
          # SampleSummaryShort,
          isSelfPromoter = isSelfPromoter,
          ABC.Score,
          distance,
          Kendall,
          ARC.E2G.Score,
          E2G.Score.qnorm)
}

# save to output file
message("Writing to output file...")
if (tools::file_ext(opt$output_file) == "gz") {
  
  # save to gzip compressed file
  tmp_file <- tools::file_path_sans_ext(opt$output_file)
  writeLines(header, con = tmp_file)
  fwrite(pred, file = tmp_file, sep = "\t", quote = FALSE, na = "NA", append = TRUE,
         col.names = TRUE)
  system2("gzip", args = c("-f", tmp_file))
  
} else {
  
  # save to uncompressed file
  writeLines(header, con = opt$output_file)
  fwrite(pred, file = opt$output_file, sep = "\t", quote = FALSE, na = "NA", append = TRUE,
         col.names = TRUE)
  
}

message("Done!")

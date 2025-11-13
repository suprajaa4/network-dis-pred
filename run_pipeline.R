#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(igraph)
})

# Load project functions
source("R/string_api.R")
source("R/build_graph.R")
source("R/features.R")
source("R/ranking.R")
source("R/evaluation.R")

option_list <- list(
  make_option("--seed_file", type = "character", default = NA),
  make_option("--outdir", type = "character", default = "results"),
  make_option("--species", type = "integer", default = 9606),
  make_option("--partner_limit", type = "integer", default = 150),
  make_option("--score_threshold", type = "double", default = 0.7),
  make_option("--val_frac", type = "double", default = 0.30)
)

opt <- parse_args(OptionParser(option_list = option_list))

DEFAULT_SEEDS <- c(
  "APP","PSEN1","PSEN2","APOE","MAPT","TREM2",
  "BIN1","CLU","PICALM","ABCA7","SORL1","CR1",
  "CD33","PLCG2"
)

if (!is.na(opt$seed_file) && file.exists(opt$seed_file)) {
  seeds <- read_tsv(opt$seed_file, show_col_types = FALSE)[[1]] |> unique()
  message(sprintf("Loaded %d seed genes", length(seeds)))
} else {
  seeds <- DEFAULT_SEEDS
  message("Using default Alzheimer’s seeds.")
}

ppis <- string_small_network(
  seeds,
  species = opt$species,
  partner_limit = opt$partner_limit,
  score_threshold = opt$score_threshold
)

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
write_tsv(ppis$named, "data/raw/string_edges_named.tsv")
write_tsv(ppis$ids,   "data/raw/string_edges_ids.tsv")

G <- build_ppi_graph(ppis$named)
feat_info <- compute_graph_features(
  G, seeds,
  damping = 0.85,
  val_frac = opt$val_frac
)

ranked <- rank_candidates(feat_info$features)
metrics <- evaluate_ranking(ranked, feat_info$val_seeds)

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

write_csv(head(ranked %>% select(gene, score), 100),
          file.path(opt$outdir, "candidates_top100.csv"))
write_lines(
  sprintf(
    "AUROC=%.4f\nRecall@50=%.4f\nRecall@100=%.4f",
    metrics$auroc, metrics$recall50, metrics$recall100
  ),
  file.path(opt$outdir, "metrics.txt")
)

cat(sprintf(
  "\nSaved → %s/candidates_top100.csv\nSaved → %s/metrics.txt\nAUROC %.3f\nRecall@50 %.3f | Recall@100 %.3f\n\n",
  opt$outdir, opt$outdir, metrics$auroc, metrics$recall50, metrics$recall100
))

# End-to-end R pipeline for Alzheimer’s: seeds → STRING API → PPI (named + IDs) → Graph → RWR → Ranking → Metrics
#
# What you get when you run this once:
# - data/raw/string_edges_named.tsv        (src, dst, score in 0..1)
# - data/raw/string_edges_ids.tsv          (stringId_A, stringId_B, score)
# - results/candidates_top100.csv          (top-ranked genes)
# - results/metrics.txt                    (AUROC, Recall@50/100)
#
# Dependencies (install once):
# install.packages(c("httr","jsonlite","dplyr","readr","purrr","stringr","igraph","pROC"))

suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
  library(dplyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(igraph)
  library(pROC)
})

# =========================
# Parameters (tweak here)
# =========================
SPECIES <- 9606                      # human
PARTNER_LIMIT <- 150                 # partners per seed to fetch (higher = bigger network)
SCORE_THRESHOLD <- 0.7               # keep interactions with score >= this (API scores are 0..1)
DAMPING <- 0.85                      # PageRank damping (1-restart prob)
VAL_FRAC <- 0.30                     # fraction of seeds to hold-out for validation
SEED_FILE <- NULL  # set to NULL to use DEFAULT_SEEDS below

DEFAULT_SEEDS <- c("APP","PSEN1","PSEN2","APOE","MAPT","TREM2","BIN1","CLU","PICALM","ABCA7",
                   "SORL1","CR1","CD33","PLCG2")

# =========================
# Helpers — STRING API
# =========================
STRING_BASE <- "https://string-db.org/api"

string_map_symbols <- function(symbols, species = SPECIES) {
  symbols <- unique(symbols)
  id_param <- paste(symbols, collapse = "%0d")  # %0d = CR between identifiers
  url <- sprintf("%s/json/get_string_ids?identifiers=%s&species=%s",
                 STRING_BASE, URLencode(id_param, reserved = TRUE), species)
  resp <- GET(url); stop_for_status(resp)
  dat <- fromJSON(content(resp, as = "text", encoding = "UTF-8"))
  if (length(dat) == 0) return(tibble())
  as_tibble(dat) %>%
    transmute(queryItem, stringId, preferredName) %>%
    distinct()
}

string_partners_one <- function(string_id, species = SPECIES, limit = PARTNER_LIMIT) {
  url <- sprintf("%s/tsv/interaction_partners?identifiers=%s&species=%s&limit=%d",
                 STRING_BASE, URLencode(string_id, reserved = TRUE), species, limit)
  resp <- GET(url); stop_for_status(resp)
  read_tsv(content(resp, as = "raw"), show_col_types = FALSE)
}

string_small_network <- function(seed_symbols, species = SPECIES,
                                 partner_limit = PARTNER_LIMIT,
                                 score_threshold = SCORE_THRESHOLD) {
  message(sprintf("Mapping %d seed symbols to STRING IDs...", length(seed_symbols)))
  map_df <- string_map_symbols(seed_symbols, species)
  if (nrow(map_df) == 0) stop("No STRING IDs found for provided symbols.")
  
  message(sprintf("Fetching partners (limit=%d, score >= %.2f)...", partner_limit, score_threshold))
  edges_list <- map(map_df$stringId, ~ string_partners_one(.x, species = species, limit = partner_limit))
  edges_raw <- bind_rows(edges_list)
  
  if (!all(c("preferredName_A","preferredName_B","stringId_A","stringId_B","score") %in% names(edges_raw))) {
    stop("Unexpected columns from STRING API. Check API response format.")
  }
  
  # Named version (human-readable gene symbols)
  edges_named <- edges_raw %>%
    filter(score >= score_threshold) %>%
    transmute(src = preferredName_A, dst = preferredName_B, score = score) %>%
    distinct() %>% filter(src != dst)
  
  # ID version (original STRING IDs)
  edges_ids <- edges_raw %>%
    filter(score >= score_threshold) %>%
    transmute(stringId_A, stringId_B, score = score) %>%
    distinct() %>% filter(stringId_A != stringId_B)
  
  list(named = edges_named, ids = edges_ids, map = map_df)
}

# =========================
# Load seeds
# =========================
if (!is.null(SEED_FILE) && file.exists(SEED_FILE)) {
  seeds <- read_tsv(SEED_FILE, show_col_types = FALSE)[[1]]
  seeds <- unique(na.omit(seeds))
  message(sprintf("Loaded %d seed symbols from %s", length(seeds), SEED_FILE))
} else {
  seeds <- DEFAULT_SEEDS
  message(sprintf("Using DEFAULT_SEEDS (%d). Set SEED_FILE to use your DisGeNET export.", length(seeds)))
}
stopifnot(length(seeds) >= 5)

# =========================
# Fetch small PPI and save both named + ID versions
# =========================
ppis <- string_small_network(seeds)

if (!dir.exists("data/raw")) dir.create("data/raw", recursive = TRUE)
write_tsv(ppis$named, "data/raw/string_edges_named.tsv")
write_tsv(ppis$ids,   "data/raw/string_edges_ids.tsv")
message("Saved: data/raw/string_edges_named.tsv and data/raw/string_edges_ids.tsv")

# =========================
# Build graph (named) and compute features
# =========================
edges <- ppis$named %>% rename(combined_score = score)
G <- graph_from_data_frame(edges %>% select(src, dst), directed = FALSE)
G <- simplify(G, remove.multiple = TRUE, remove.loops = TRUE)

# Largest connected component
comp <- components(G)
G <- induced_subgraph(G, which(comp$membership == which.max(comp$csize)))
V(G)$name <- make.names(V(G)$name, unique = TRUE)

# Intersect seeds with present nodes
seed_genes <- intersect(seeds, V(G)$name)
if (length(seed_genes) < 3) stop("<3 seeds present in graph after filtering; try lowering SCORE_THRESHOLD or raising PARTNER_LIMIT.")

# Centralities
feat <- tibble(
  gene = V(G)$name,
  degree = degree(G),
  betweenness = betweenness(G, directed = FALSE, normalized = TRUE),
  closeness = suppressWarnings(closeness(G, normalized = TRUE))
)
pr <- page_rank(G, directed = FALSE, damping = DAMPING)$vector
feat$pagerank <- pr[feat$gene]

# =========================
# Personalized PageRank (RWR) with hold-out validation
# =========================
set.seed(7)
val_n <- max(1, floor(length(seed_genes) * VAL_FRAC))
val_seeds <- sample(seed_genes, val_n)
train_seeds <- setdiff(seed_genes, val_seeds)

pers <- rep(0, vcount(G)); names(pers) <- V(G)$name
pers[train_seeds] <- 1/length(train_seeds)
pr_personal <- page_rank(G, directed = FALSE, damping = DAMPING, personalized = pers)$vector
feat$rwr <- pr_personal[feat$gene]

# =========================
# Aggregate into a single score (rank → z → sum)
# =========================
num_cols <- c("degree","betweenness","closeness","pagerank","rwr")
rank_mat <- apply(as.matrix(feat[, num_cols]), 2, function(x) rank(x, ties.method = "average"))
z_mat <- scale(rank_mat)
feat$score <- rowSums(z_mat)

# =========================
# Validation metrics
# =========================
labels <- ifelse(feat$gene %in% val_seeds, 1, 0)
roc_obj <- tryCatch(roc(response = labels, predictor = feat$score, quiet = TRUE), error = function(e) NULL)
auroc <- if (!is.null(roc_obj)) as.numeric(auc(roc_obj)) else NA_real_

ranked <- feat %>% arrange(desc(score))
recall_at_k <- function(k) {
  topk <- head(ranked$gene, k)
  sum(topk %in% val_seeds) / length(val_seeds)
}
recall50  <- recall_at_k(min(50, nrow(ranked)))
recall100 <- recall_at_k(min(100, nrow(ranked)))

# =========================
# Save outputs
# =========================
if (!dir.exists("results")) dir.create("results", recursive = TRUE)
write_csv(ranked %>% select(gene, score) %>% head(100), "results/candidates_top100.csv")
write_lines(sprintf("AUROC=%.4f\nRecall@50=%.4f\nRecall@100=%.4f", auroc, recall50, recall100), "results/metrics.txt")

cat(sprintf("\nSaved results -> results/candidates_top100.csv and results/metrics.txt\nAUROC: %s\nRecall@50: %.3f | Recall@100: %.3f\n\n",
            ifelse(is.na(auroc), "NA", sprintf("%.3f", auroc)), recall50, recall100))

# =========================
# Tips
# =========================
# - If AUROC is NA or very low, increase PARTNER_LIMIT or lower SCORE_THRESHOLD (e.g., 0.6) to include more edges.
# - If seeds drop out, ensure seed symbols match STRING preferred names; consider adding synonyms to DEFAULT_SEEDS.
# - For larger graphs, consider approximate betweenness or skip it for speed.
# - To add expression features later, join a data.frame with columns gene + logFC (or |logFC|) to `feat` before aggregation.

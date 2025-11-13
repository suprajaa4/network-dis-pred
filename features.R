# R/features.R

suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
})

compute_graph_features <- function(G, seeds,
                                   damping = 0.85,
                                   val_frac = 0.30,
                                   rng_seed = 7) {
  
  seed_genes <- intersect(seeds, V(G)$name)
  if (length(seed_genes) < 3)
    stop("<3 seeds present after filtering — adjust STRING parameters.")
  
  feat <- tibble(
    gene = V(G)$name,
    degree = degree(G),
    betweenness = betweenness(G, normalized = TRUE),
    closeness = suppressWarnings(closeness(G, normalized = TRUE))
  )
  
  pr <- page_rank(G, damping = damping)$vector
  feat$pagerank <- pr[feat$gene]
  
  set.seed(rng_seed)
  val_n <- max(1, floor(length(seed_genes) * val_frac))
  val_seeds <- sample(seed_genes, val_n)
  train_seeds <- setdiff(seed_genes, val_seeds)
  
  pers <- rep(0, vcount(G)); names(pers) <- V(G)$name
  pers[train_seeds] <- 1/length(train_seeds)
  
  pr_personal <- page_rank(G, damping = damping, personalized = pers)$vector
  feat$rwr <- pr_personal[feat$gene]
  
  list(
    features = feat,
    train_seeds = train_seeds,
    val_seeds = val_seeds
  )
}

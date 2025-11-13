# R/ranking.R

suppressPackageStartupMessages(library(dplyr))

rank_candidates <- function(features) {
  num_cols <- c("degree", "betweenness", "closeness", "pagerank", "rwr")
  
  rank_mat <- apply(
    as.matrix(features[, num_cols]),
    2,
    function(x) rank(x, ties.method = "average")
  )
  
  z_mat <- scale(rank_mat)
  features$score <- rowSums(z_mat)
  
  arrange(features, desc(score))
}

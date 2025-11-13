# R/build_graph.R

suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(dplyr))

build_ppi_graph <- function(edges_named) {
  edges <- rename(edges_named, combined_score = score)
  
  G <- graph_from_data_frame(edges %>% select(src, dst), directed = FALSE)
  G <- simplify(G, remove.multiple = TRUE, remove.loops = TRUE)
  
  comp <- components(G)
  G <- induced_subgraph(G, which(comp$membership == which.max(comp$csize)))
  V(G)$name <- make.names(V(G)$name, unique = TRUE)
  
  G
}

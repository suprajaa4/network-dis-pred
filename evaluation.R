# R/evaluation.R

suppressPackageStartupMessages({
  library(dplyr)
  library(pROC)
})

evaluate_ranking <- function(ranked, val_seeds) {
  labels <- ifelse(ranked$gene %in% val_seeds, 1, 0)
  
  roc_obj <- tryCatch(
    roc(response = labels, predictor = ranked$score, quiet = TRUE),
    error = function(e) NULL
  )
  
  auroc <- if (!is.null(roc_obj)) as.numeric(auc(roc_obj)) else NA_real_
  
  recall_at_k <- function(k) {
    topk <- head(ranked$gene, min(k, nrow(ranked)))
    sum(topk %in% val_seeds) / length(val_seeds)
  }
  
  list(
    auroc = auroc,
    recall50 = recall_at_k(50),
    recall100 = recall_at_k(100)
  )
}

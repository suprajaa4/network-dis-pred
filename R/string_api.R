# R/string_api.R

suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
  library(dplyr)
  library(readr)
  library(purrr)
})

STRING_BASE <- "https://string-db.org/api"

string_map_symbols <- function(symbols, species = 9606) {
  symbols <- unique(symbols)
  id_param <- paste(symbols, collapse = "%0d")
  url <- sprintf(
    "%s/json/get_string_ids?identifiers=%s&species=%s",
    STRING_BASE,
    URLencode(id_param, reserved = TRUE),
    species
  )
  resp <- GET(url); stop_for_status(resp)
  dat <- fromJSON(content(resp, as = "text", encoding = "UTF-8"))
  if (length(dat) == 0) return(tibble())
  
  as_tibble(dat) %>%
    transmute(queryItem, stringId, preferredName) %>%
    distinct()
}

string_partners_one <- function(string_id, species = 9606, limit = 150) {
  url <- sprintf(
    "%s/tsv/interaction_partners?identifiers=%s&species=%s&limit=%d",
    STRING_BASE,
    URLencode(string_id, reserved = TRUE),
    species,
    limit
  )
  resp <- GET(url); stop_for_status(resp)
  read_tsv(content(resp, as = "raw"), show_col_types = FALSE)
}

string_small_network <- function(seed_symbols,
                                 species = 9606,
                                 partner_limit = 150,
                                 score_threshold = 0.7) {
  message(sprintf("Mapping %d seed symbols to STRING IDs...", length(seed_symbols)))
  map_df <- string_map_symbols(seed_symbols, species)
  if (nrow(map_df) == 0) stop("No STRING IDs found for provided symbols.")
  
  message(sprintf("Fetching partners (limit=%d, score>=%.2f)...",
                  partner_limit, score_threshold))
  
  edges_list <- purrr::map(
    map_df$stringId,
    ~ string_partners_one(.x, species = species, limit = partner_limit)
  )
  edges_raw <- bind_rows(edges_list)
  
  stopifnot(all(c("preferredName_A","preferredName_B",
                  "stringId_A","stringId_B","score") %in% names(edges_raw)))
  
  edges_named <- edges_raw %>%
    filter(score >= score_threshold) %>%
    transmute(src = preferredName_A, dst = preferredName_B, score) %>%
    distinct() %>%
    filter(src != dst)
  
  edges_ids <- edges_raw %>%
    filter(score >= score_threshold) %>%
    transmute(stringId_A, stringId_B, score) %>%
    distinct() %>%
    filter(stringId_A != stringId_B)
  
  list(named = edges_named, ids = edges_ids, map = map_df)
}

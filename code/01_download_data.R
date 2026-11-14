# =============================================================================
# 01_download_data.R  (v3)
# Download influenza A (H3N2) haemagglutinin (HA) sequences from NCBI GenBank.
#
# Run with the working directory set to the code/ folder.
# Output: ../data/h3n2_ha.fasta
#
# v3: tries several query styles per year and uses the first one that returns
# records. It also prints how NCBI interpreted the query, which makes any
# remaining problem easy to diagnose.
# =============================================================================

# install.packages("rentrez")
library(rentrez)

set.seed(2026)

years    <- 2014:2024   # seasons to sample
per_year <- 25          # max sequences sampled per year
out_dir  <- "../data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Optional: an NCBI API key allows more requests per second
# set_entrez_key("YOUR_KEY_HERE")

safe_search <- function(term) {
  tryCatch(entrez_search(db = "nuccore", term = term, retmax = 300),
           error = function(e) {
             message("  Search error: ", conditionMessage(e))
             NULL
           })
}

# Candidate queries, from strictest to most forgiving.
# Note: no quoted "complete cds" phrase, which NCBI sometimes fails to match.
candidate_terms <- function(yr) {
  c(
    paste0('"Influenza A virus"[Organism] AND H3N2[Title] AND hemagglutinin[Title] ',
           'AND complete[Title] AND cds[Title] AND ', yr, '[Title]'),
    paste0('"Influenza A virus"[Organism] AND H3N2[Title] AND hemagglutinin[Title] ',
           'AND complete[Title] AND cds[Title] AND ', yr, '[PDAT]'),
    paste0('Influenza A virus H3N2 hemagglutinin complete cds AND ', yr, '[PDAT]'),
    paste0('H3N2 hemagglutinin AND ', yr, '[PDAT]')
  )
}

fetch_year <- function(yr) {
  res <- NULL
  for (term in candidate_terms(yr)) {
    r <- safe_search(term)
    if (!is.null(r) && length(r$ids) >= 10) { res <- r; break }
    if (!is.null(r)) {
      message("  No hits for: ", term)
      message("  NCBI read it as: ", r$QueryTranslation)
    }
  }

  if (is.null(res)) {
    message("Year ", yr, ": no sequences found with any query")
    return(character(0))
  }

  ids <- res$ids
  if (length(ids) > per_year) ids <- sample(ids, per_year)
  message("Year ", yr, ": fetching ", length(ids), " of ", res$count, " records")

  Sys.sleep(0.5)  # be polite to NCBI
  fa <- tryCatch(entrez_fetch(db = "nuccore", id = ids,
                              rettype = "fasta", retmode = "text"),
                 error = function(e) {
                   message("  Fetch error: ", conditionMessage(e))
                   NA_character_
                 })
  if (is.na(fa)) return(character(0))
  fa
}

all_fasta <- unlist(lapply(years, fetch_year))

if (length(all_fasta) == 0) {
  stop("No sequences were downloaded. Check the messages above.")
}

writeLines(all_fasta, file.path(out_dir, "h3n2_ha.fasta"))
message("Done. Wrote ", file.path(out_dir, "h3n2_ha.fasta"),
        " (", length(all_fasta), " year batches)")
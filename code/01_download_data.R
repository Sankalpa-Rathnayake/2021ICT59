# =============================================================================
# 01_download_data.R
# Download influenza A (H3N2) haemagglutinin (HA) sequences from NCBI GenBank.
#
# Run this script with the working directory set to the code/ folder.
# Output: ../data/h3n2_ha.fasta
#
# A fixed random seed and a fixed number of sequences per year keep the
# download reproducible. Note that NCBI adds new sequences over time, so the
# exact sampled set can change if you re-run this later; the committed FASTA
# file in data/ is the version used for the report.
# =============================================================================

# install.packages("rentrez")
library(rentrez)

set.seed(2026)

years    <- 2014:2024   # seasons to sample
per_year <- 25          # max sequences sampled per year (keeps runtime low)
out_dir  <- "../data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Optional: an NCBI API key allows more requests per second
# set_entrez_key("YOUR_KEY_HERE")

fetch_year <- function(yr) {
  # Strain names look like A/Texas/50/2012(H3N2), so the year is a title token
  term <- paste0('"Influenza A virus"[Organism] AND H3N2[Title] AND ',
                 'hemagglutinin[Title] AND "complete cds"[Title] AND ',
                 yr, '[Title]')

  res <- tryCatch(entrez_search(db = "nuccore", term = term, retmax = 300),
                  error = function(e) NULL)
  if (is.null(res) || length(res$ids) == 0) {
    message("Year ", yr, ": no sequences found")
    return(character(0))
  }

  ids <- res$ids
  if (length(ids) > per_year) ids <- sample(ids, per_year)
  message("Year ", yr, ": fetching ", length(ids), " of ", res$count, " records")

  # Fetch in one request per year (<= 25 records)
  Sys.sleep(0.5)  # be polite to NCBI
  fa <- tryCatch(entrez_fetch(db = "nuccore", id = ids,
                              rettype = "fasta", retmode = "text"),
                 error = function(e) NA_character_)
  if (is.na(fa)) return(character(0))
  fa
}

all_fasta <- unlist(lapply(years, fetch_year))
writeLines(all_fasta, file.path(out_dir, "h3n2_ha.fasta"))

message("Done. Wrote ", file.path(out_dir, "h3n2_ha.fasta"))

library(readr)
suppressPackageStartupMessages(library(dplyr))
library(stringr)
library(purrr)

input <- c(
  "intermediate/pre-impute_filter/output/nacc-gsa_cluster_pops.tsv",
  "intermediate/pre-impute_filter/output/nacc-HumanOmniExpressExome_cluster_pops.tsv",
  "intermediate/pre-impute_filter/output/mssm-gsa_cluster_pops.tsv")

input <- snakemake@input
output <- snakemake@output[[1]]

read_pops <- function(fname) {
  id <- stringr::str_extract(fname, "(?<=output/).+(?=_cluster_pops.tsv$)")
  fname |>
    readr::read_tsv(col_types = cols(
      .default = "-", superpop_infered = "c", cohort = "c")) |>
    dplyr::filter(cohort != "Reference") |>
    dplyr::count(superpop_infered) |>
    dplyr::mutate(identifier = id) |>
    dplyr::rename(ancestry = superpop_infered)
}

input |>
  purrr::map_dfr(read_pops) |>
  tidyr::unite("identifier_ancestry",
    c("identifier", "ancestry"),
    sep = "_", remove = FALSE) |>
  dplyr::select(identifier, ancestry,
                identifier_ancestry, n) |>
  readr::write_tsv(output)

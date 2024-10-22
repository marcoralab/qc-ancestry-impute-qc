library(dplyr)
library(readr)

in_fam <- snakemake@input[[1]]
out_sexcheck <- snakemake@output[[1]]

in_fam |>
  read_table(col_names = c("FID", "IID", "PEDSEX"), col_types = "cc--i-") |>
  mutate(SNPSEX = 0, STATUS = "OK", Fisher = 0) |>
  select(FID, IID, PEDSEX, SNPSEX, STATUS, F = Fisher) |>
  write_tsv(out_sexcheck)
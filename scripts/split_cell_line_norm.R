library(tidyverse)
library(CRISPRcleanR)
data(AVANA_Library)
library(furrr)

dataset <- read_tsv("results/avana_lfc_rename.tsv") |> 
           dplyr::rename("gene" = "Gene") |>
           left_join(AVANA_Library |> 
                    select(CODE, seq, GENES), by = c("guide" = "seq", "gene" = "GENES")) |>
select(-guide) |>
dplyr::rename("sgRNA" = "CODE") |> 
relocate(sgRNA, .before = gene) |> 
filter(!is.na(sgRNA)) 


sample_names <- names(dataset)[-c(1, 2)]

# 2. Determine the grouping prefix for each sample column.
group_prefix <- substr(sample_names, 1, 10)


# Create a list of column indices for each group (e.g., "ACH0999913", etc.)
group_indices <- split(seq_along(sample_names), group_prefix)

group_indices <- map(group_indices, ~ .x + 2)

split_dfs <- map(group_indices, ~ dataset[, c(1, 2, .x), drop = FALSE])

plan(multisession, workers = 8)

ccr_values <- future_map(split_dfs, ~ccr.logFCs2chromPos(.x,AVANA_Library))
ccr_adjustments <- future_map(ccr_values, ~ccr.adjust(.x,AVANA_Library))


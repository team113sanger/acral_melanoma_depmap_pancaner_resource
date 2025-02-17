library(tidyverse)
library(CRISPRcleanR)
data(AVANA_Library)
library(furrr)

dataset <- read_tsv("results/avana_counts_rename.tsv") |> 
           dplyr::rename("gene" = "Gene") |>
           left_join(AVANA_Library |> 
                    select(CODE, seq, GENES), by = c("guide" = "seq", "gene" = "GENES")) |>
select(-guide) |>
dplyr::rename("sgRNA" = "CODE") |> 
relocate(sgRNA, .before = gene) |> 
filter(!is.na(sgRNA)) |> 
relocate(contains("pDNA_batch_Avana"), .after = gene)


sample_names <- names(dataset)[-c(1,2,3,4,5)]
sample_batches <- str_extract(sample_names, "Avana-\\d+")
pdna_map <- paste0("pDNA_batch_", str_replace(sample_batches, "-", "_"))


positions <- split(seq_along(pdna_map), pdna_map)
screens <- dataset[,-c(1,2,3,4,5)]
avana_datasets <- map(positions, ~screens[,.x])


# 2. Determine the grouping prefix for each sample column.

split_cell_lines <- function(avana_dataset, ref){
    ids <- names(avana_dataset)
    group_prefix <- substr(ids, 1, 10)
    group_indices <- split(seq_along(ids), group_prefix)

    relevel <- ref |> bind_cols(avana_dataset)
    new_indexes <- map(group_indices, ~.x + 3)
    map(new_indexes, ~relevel[, c(1,2,3,.x), drop = FALSE])
}


avana_v2 <- split_cell_lines(avana_datasets[[1]], dataset[,c(1,2,5)])
avana_v3 <- split_cell_lines(avana_datasets[[2]], dataset[,c(1,2,4)])
avana_v4 <- split_cell_lines(avana_datasets[[3]], dataset[,c(1,2,3)])



plan(multisession, workers = 4)
normANDfcs <- future_map(avana_v4, ~ccr.NormfoldChanges(Dframe = .x,
                                  saveToFig = FALSE,
                                  min_reads = 30, 
                                  EXPname = "ccr", 
                                  libraryAnnotation = AVANA_Library, 
                                  display = FALSE))



ccr_values <- future_map(normANDfcs, 
                        ~ccr.logFCs2chromPos(.x[["logFCs"]],
                        AVANA_Library))

ccr_adjustments <- future_map(ccr_values, 
                        ~ccr.GWclean(.x,
                            display = FALSE, 
                            label = "crisprcleanr"))

ccr.GWclean(ccr_adjustments[[1]], display = FALSE, 
                            label = "crisprcleanr")
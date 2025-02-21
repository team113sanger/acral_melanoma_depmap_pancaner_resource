
################################################################################
# Setup
################################################################################

library(tidyverse)
library(CRISPRcleanR)
data(AVANA_Library)
library(furrr)
library(glue)


################################################################################
# Read in the data
################################################################################


# 1. Read in the data and join with the AVANA library so that is is annotated 
# for CRISPRcleanR.
dataset <- read_tsv("results/avana_counts_rename.tsv") |>
  dplyr::rename("gene" = "Gene") |>
  left_join(AVANA_Library |>
    select(CODE, seq, GENES), by = c("guide" = "seq", "gene" = "GENES")) |>
  select(-guide) |>
  dplyr::rename("sgRNA" = "CODE") |>
  relocate(sgRNA, .before = gene) |>
  filter(!is.na(sgRNA)) |>
  relocate(contains("pDNA_batch_Avana"), .after = gene)


# Divide the dataset into separate screens based on the pDNA batch.
sample_names <- names(dataset)[-c(1, 2, 3, 4, 5)]
sample_batches <- str_extract(sample_names, "Avana-\\d+")
pdna_map <- paste0("pDNA_batch_", str_replace(sample_batches, "-", "_"))


positions <- split(seq_along(pdna_map), pdna_map)
screens <- dataset[, -c(1, 2, 3, 4, 5)]
avana_datasets <- map(positions, ~ screens[, .x])


# 2. Determine the grouping prefix for each sample column.

split_cell_lines <- function(avana_dataset, ref) {
  ids <- names(avana_dataset)
  group_prefix <- substr(ids, 1, 10)
  group_indices <- split(seq_along(ids), group_prefix)

  relevel <- ref |> bind_cols(avana_dataset)
  new_indexes <- map(group_indices, ~ .x + 3)
  fetched_cols <- map(new_indexes, ~ relevel[, c(1, 2, 3, .x), drop = FALSE])
  map(fetched_cols, na.omit)
}


# 2. Determine the grouping prefix for each sample column.


avana_v2 <- split_cell_lines(avana_datasets[[1]], 
                             dataset[, c(1, 2, 5)])
avana_v3 <- split_cell_lines(avana_datasets[[2]], 
                             dataset[, c(1, 2, 4)])
avana_v4 <- split_cell_lines(avana_datasets[[3]], 
                             dataset[, c(1, 2, 3)])


################################################################################
# CRISPRcleanR correction
################################################################################

# Paralellise work 
plan(multisession, workers = 4)

dataset <- avana_v4
# calculate_cell_line <- function(dataset, dest){
normANDfcs <- future_map(dataset, ~ ccr.NormfoldChanges(
  Dframe = .x,
  saveToFig = FALSE,
  EXPname = "ccr",
  min_reads = 30,
  libraryAnnotation = AVANA_Library,
  display = FALSE
))

ccr_values <- future_map(normANDfcs, ~ ccr.logFCs2chromPos(
  .x[["logFCs"]],
  AVANA_Library
))

ccr_adjustments <- imap(ccr_values, ~ ccr.GWclean(.x,
  display = FALSE,
  label = .y
))

dest <- "avana_v4_"

fetch_corrected_fc <- function(x) {
  x[["corrected_logFCs"]] |>
    tibble() |>
    mutate(sgRNA = x[["SORTED_sgRNAs"]]) |>
    select(sgRNA, genes, correctedFC)
}


fold_ch <- map_dfr(ccr_adjustments, fetch_corrected_fc, .id = "cell_line") |>
  pivot_wider(names_from = cell_line, values_from = correctedFC)

write_tsv(fold_ch, glue("./results/avana_v4_ccr_fold_changes.tsv"))


corrected_counts <- imap(ccr_adjustments, ~ ccr.correctCounts(.y,
  normANDfcs[[.y]][["norm_counts"]],
  .x,
  AVANA_Library,
  minTargetedGenes = 3,
  OutDir = glue("./results/{dest}")
))

imap(corrected_counts, ~ ccr.PlainTsvFile(.x,
  fprefix = .y,
  path = glue("./results/{dest}")
))
################################################################################
# Apply to Avana 3 samples 
###############################################################################


dataset_3 <- avana_v3
# calculate_cell_line <- function(dataset, dest){
normANDfcs_3 <- future_map(dataset_3, ~ ccr.NormfoldChanges(
  Dframe = .x,
  saveToFig = FALSE,
  EXPname = "ccr",
  min_reads = 30,
  libraryAnnotation = AVANA_Library,
  display = FALSE
))

ccr_values_3 <- future_map(normANDfcs_3, ~ ccr.logFCs2chromPos(
  .x[["logFCs"]],
  AVANA_Library
))

ccr_adjustments_3 <- imap(ccr_values_3, ~ ccr.GWclean(.x,
  display = FALSE,
  label = .y
))



fold_ch_3 <- map_dfr(ccr_adjustments_3, fetch_corrected_fc, .id = "cell_line") |>
  pivot_wider(names_from = cell_line, values_from = correctedFC)

fold_ch_3 |> write_tsv(glue("./results/avana_v3_ccr_fold_changes.tsv"))

dest <- "avana_v3_"
corrected_counts_3 <- imap(ccr_adjustments_3, ~ ccr.correctCounts(.y,
  normANDfcs_3[[.y]][["norm_counts"]],
  .x,
  AVANA_Library,
  minTargetedGenes = 3,
  OutDir = glue("./results/{dest}")
))


imap(corrected_counts_3, ~ ccr.PlainTsvFile(.x,
  fprefix = .y,
  path = glue("./results/{dest}")
))


################################################################################
# Apply to Avana 2 samples 
###############################################################################

dataset_2 <- avana_v2
# calculate_cell_line <- function(dataset, dest){
normANDfcs_2 <- future_map(dataset_2, ~ ccr.NormfoldChanges(
  Dframe = .x,
  saveToFig = FALSE,
  EXPname = "ccr",
  min_reads = 30,
  libraryAnnotation = AVANA_Library,
  display = FALSE
))

ccr_values_2 <- future_map(
  normANDfcs_2,
  ~ ccr.logFCs2chromPos(
    .x[["logFCs"]],
    AVANA_Library
  )
)

ccr_adjustments_2 <- imap(ccr_values_2, ~ ccr.GWclean(.x,
  display = FALSE,
  label = .y
))


fold_ch_2 <- map_dfr(ccr_adjustments_2, fetch_corrected_fc, .id = "cell_line") |>
  pivot_wider(names_from = cell_line, values_from = correctedFC)

fold_ch_2 |> write_tsv(glue("./results/avana_v2_ccr_fold_changes.tsv"))

dest <- "avana_v2_"
corrected_counts_2 <- imap(ccr_adjustments_2, ~ ccr.correctCounts(.y,
  normANDfcs_2[[.y]][["norm_counts"]],
  .x,
  AVANA_Library,
  minTargetedGenes = 3,
  OutDir = glue("./results/{dest}")
))


imap(corrected_counts_2, ~ ccr.PlainTsvFile(.x,
  fprefix = .y,
  path = glue("./results/{dest}")
))

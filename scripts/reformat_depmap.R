

################################################################################
# Setup
################################################################################

library(here)
library(tidyverse)
library(glue)
# renv::install("francescojm/CRISPRcleanR")
library(CRISPRcleanR)
library(data.table)
library(depmap)

################################################################################
# Reading in the data
################################################################################
models <- read_csv(here("data/ScreenSequenceMap.csv")) |> 
  mutate(ModelID = case_when(!is.na(Replicate) ~ glue("{ModelID}_{Replicate}_{pDNABatch}_{ModelConditionID}"),
                             TRUE ~ SequenceID)) |> 
  select(ModelID, SequenceID)


guide_map <- read_csv(here("data/AvanaGuideMap.csv")) |> 
  select(sgRNA, Gene) |> 
  mutate(Gene = str_remove(Gene, pattern = " \\(.*\\)" )) |> 
  distinct()


avana_counts <- fread(here("data/AvanaLogfoldChange.csv"), fill = TRUE) 

rename_map <- setNames(models$ModelID, models$SequenceID)

common_cols <- intersect(names(avana_counts), names(rename_map))


# Rename only the matching columns
names(avana_counts) <- c("guide", rename_map[common_cols])
avana_counts |> 
left_join(guide_map, by = c("guide" = "sgRNA")) |>
relocate(Gene, .after = guide) |>
as_tibble()
    


# Load Avana
data(AVANA_Library)


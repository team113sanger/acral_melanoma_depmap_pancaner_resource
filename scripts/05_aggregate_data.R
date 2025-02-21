
################################################################################
# Setup
################################################################################

library(tidyverse)

# Read in PR data output from BAGEL
pr_files <- fs::dir_ls("results/pr_curves")

cell_line_data <- read_tsv(pr_files, id = "cell_line")

cell_line_data <- cell_line_data |>
  mutate(cell_line = str_remove(basename(cell_line), pattern = "_pr.txt")) |>
  separate(cell_line, into = c("avana_replicate", "model"), sep = "_A") |>
  mutate(model = paste0("A", model))


# Find models that appear in multiple Avana pDNA batches and select only one set of 
# replicates for each model to avoid double counting.
unique_models <- cell_line_data |>
  select(model, avana_replicate) |>
  group_by(model) |>
  slice_head(n = 1) |>
  ungroup()


one_per_model <- unique_models |>
  left_join(cell_line_data, by = c("model", "avana_replicate"))


# Pivot the data to wide format with cell line as cols and Gene as rows
one_per_model |>
  select(model, BF, Gene) |>
  pivot_wider(names_from = model, values_from = BF) |>
  write_tsv("avana_unique_bayes_factors_wide.tsv")


one_per_model |>
write_tsv("avana_unique_bayes_factors.tsv")

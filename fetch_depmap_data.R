
################################################################################
# Setup
################################################################################

library(here)
library(tidyverse)
library(glue)
library(data.table)
library(curl)


################################################################################
# Downloads
################################################################################


# Fetch the LogFoldChange file from figshare
curl::curl_download("https://figshare.com/ndownloader/files/51063578", destfile = here("data", "AvanaLogfoldChange.csv"))
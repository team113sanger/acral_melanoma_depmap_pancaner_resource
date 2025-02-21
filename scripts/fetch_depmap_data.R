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

# Fetch the raw counts file from figshare
curl::curl_download("https://figshare.com/ndownloader/files/51064568", destfile = here("data", "AvanaRawReadcounts.csv"))

# Fetch the LogFoldChange file from figshare
curl::curl_download("https://figshare.com/ndownloader/files/51063578", destfile = here("data", "AvanaLogfoldChange.csv"))

# Fetch the LogFoldChange file from figshare
curl::curl_download("https://figshare.com/ndownloader/files/51063575", destfile = here("data", "AvanaGuideMap.csv"))

# Fetch the ScreenMap file from figshare
curl::curl_download("https://figshare.com/ndownloader/files/51065828", destfile = here("data", "ScreenSequenceMap.csv"))


curl::curl_download("https://figshare.com/ndownloader/files/51065828", destfile = here("data", "ScreenSequenceMap.csv"))

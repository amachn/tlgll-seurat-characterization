# ~~~ scripts/00_config.R ~~~
library(here)

# - project paths
DATA_DIR <- here("data")
RAW_DIR <- here("data", "raw")
PROCESSED_DIR <- here("data", "processed")
RESULTS_DIR <- here("results")

dir.create(DATA_DIR, showWarnings = FALSE)
dir.create(RAW_DIR, showWarnings = FALSE)
dir.create(PROCESSED_DIR, showWarnings = FALSE)
dir.create(RESULTS_DIR, showWarnings = FALSE)

# - dataset
ACCESSION <- "GSE168859"

# - key files

# - reproducibility
SEED <- 1
set.seed(SEED)

# - Seurat defaults

# - param grid
#PARAM_GRID <- expand.grid()

message("Project configuration loaded.")
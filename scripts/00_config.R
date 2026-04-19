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

# - dataset info
GEO_ACCESSION <- "GSE168859"
GROUP_LEVELS <- c("healthy", "pretreatment", "posttreatment")

# - key files
MANIFEST_FILE <- here(PROCESSED_DIR, "gse168859_sample_manifest.csv")

# - reproducibility
SEED <- 1
set.seed(SEED)

# - Seurat defaults
DEFAULT_MIN_FEATURES <- 200 # ?
DEFAULT_MAX_MT <- 15 # ?

# - param grid
#PARAM_GRID <- expand.grid()

message("Project configuration loaded.")
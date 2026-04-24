# ~~~ scripts/00_config.R ~~~
library(here)

get_config <- function() {
  cfg <- list(
    # - project paths
    DATA_DIR = here("data"),
    RAW_DIR = here("data", "raw"),
    PROCESSED_DIR = here("data", "processed"),
    RESULTS_DIR = here("results"),
    
    # - dataset info
    GEO_ACCESSION = "GSE168859",
    GROUP_LEVELS = c("healthy", "pretreatment", "posttreatment"),
    BASELINE_GROUP = c("healthy", "pretreatment"),
    TREATMENT_GROUP = c("pretreatment", "posttreatment"),
  
    # - reproducibility
    SEED = 1,
    
    # - Seurat defaults
    DEFAULT_MIN_FEATURES = 200,
    DEFAULT_MAX_FEATURES = 2500,
    DEFAULT_MAX_MT = 10,
    
    DEFAULT_NPCS = 20, 
    DEFAULT_K_PARAM = 20,
    DEFAULT_RESOLUTION = 0.5,
    DEFAULT_PERPLEXITY = 30,
    
    # - param grid
    PARAM_GRID = expand.grid(
      npcs = c(10, 20, 30),
      # k.param = c(10, 20, 30),
      resolution = c(0.2, 0.5, 0.8)
      # perplexity = c(20, 30, 50)
    )
  )
  
  # - key files
  cfg$MANIFEST_FILE <- here(cfg$PROCESSED_DIR, "gse168859_sample_manifest.csv")
  cfg$MERGED_SEURAT_FILE <- here(cfg$PROCESSED_DIR, "tlgll_merged_seurat.rds")
  cfg$BASELINE_SEURAT_FILE <- here(cfg$PROCESSED_DIR, "tlgll_baseline_seurat.rds")
  cfg$TREATMENT_SEURAT_FILE <- here(cfg$PROCESSED_DIR, "tlgll_treatment_seurat.rds")
  
  # - create directories if missing
  dir.create(cfg$DATA_DIR, showWarnings = FALSE)
  dir.create(cfg$RAW_DIR, showWarnings = FALSE)
  dir.create(cfg$PROCESSED_DIR, showWarnings = FALSE)
  dir.create(cfg$RESULTS_DIR, showWarnings = FALSE)
  
  # - apply seed
  set.seed(cfg$SEED)
  
  cfg
}

# ~~~ scripts/04_seurat_analysis.R ~~~
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(here)
})

source(here("scripts", "00_config.R"), local = TRUE)
cfg <- get_config()

# - paths -

# - parameters - 

# - module score gene sets -
cytotoxic_genes <- c() # placeholder: use cytotoxic pathway set?
apoptosis_genes <- c() # placeholder: use apoptosis pathway set?
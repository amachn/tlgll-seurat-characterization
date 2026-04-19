# ~~~ scripts/02_seurat_prep.R ~~~
suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(dplyr)
  library(here)
  library(purrr)
  library(readr)
  library(stringr)
})

source(here("scripts", "00_config.R"))

# - 1 ~ load manifest
manifest <- read_csv(MANIFEST_FILE, show_col_types = FALSE)

expr_manifest <- manifest |>
  filter(data_type == "geneexpression")

# - 2 ~ find extracted count files
EXTRACTED_DIR <- file.path(RAW_DIR, GEO_ACCESSION, "extracted")

count_files <- list.files(EXTRACTED_DIR, full.names = TRUE)
file_tbl <- tibble(
  file_path = count_files,
  file_name = basename(count_files),
  gsm = str_extract(file_name, "GSM\\d+")
)

expr_manifest <- expr_manifest |>
  left_join(file_tbl, by = "gsm")

# - 3 ~ build one Seurat obj per sample
make_seurat_obj <- function(file_path, sample_row) {
  counts <- as.matrix(read.csv(file_path, row.names = 1)) |> 
    Matrix(sparse = TRUE)
  
  obj <- suppressWarnings(
    CreateSeuratObject(counts = counts, project = "T-LGLL"))

  obj$gsm <- sample_row$gsm
  obj$subject_id <- sample_row$subject_id
  obj$group <- factor(sample_row$group, levels = GROUP_LEVELS)
  obj$title <- sample_row$title
  
  obj
}

# TODO: this is highly inefficient and is not completing in a reasonable amount
# of time, need to come up with an alternative solution.
obj_list <- map(
  seq_len(nrow(expr_manifest)),
  function(i) {
    message(paste0("Building Seurat object ", i, "/", nrow(expr_manifest)))
    make_seurat_obj(
      file_path = expr_manifest$file_path[[i]],
      sample_row = expr_manifest[i, ]
    )
  }
)

names(obj_list) <- expr_manifest$gsm

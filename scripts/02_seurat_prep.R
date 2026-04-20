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
manifest <- read_csv(MANIFEST_FILE, show_col_types = FALSE) |>
  filter(data_type == "geneexpression")

# - 2 ~ find extracted count files
EXTRACTED_DIR <- file.path(RAW_DIR, GEO_ACCESSION, "extracted")

count_files <- list.files(EXTRACTED_DIR, full.names = TRUE)
file_tbl <- tibble(
  file_path = count_files,
  file_name = basename(count_files),
  gsm = str_extract(file_name, "GSM\\d+")
)

manifest <- manifest |> left_join(file_tbl, by = "gsm")

# ~ output directories
SAMPLE_OBJ_DIR <- here(PROCESSED_DIR, "sample_objects")
BATCH_OBJ_DIR <- here(PROCESSED_DIR, "batch_objects")

dir.create(SAMPLE_OBJ_DIR, showWarnings = FALSE)
dir.create(BATCH_OBJ_DIR, showWarnings = FALSE)

# - 3 ~ build one Seurat obj per sample
# * will likely take a long time to run (30-60 minutes)
for (i in seq_len(nrow(manifest))) {
  row <- manifest[i, ]
  
  out_file <- here(SAMPLE_OBJ_DIR, paste0(row$gsm, ".rds"))
  if (file.exists(out_file)) {
    next
  }
  
  message("Building object for ", row$gsm, " (", i, "/", nrow(manifest), ")")
  counts <- as.matrix(read.csv(row$file_path, row.names = 1)) |>
    Matrix(sparse = TRUE)
  
  obj <- suppressWarnings(CreateSeuratObject(
    counts = counts, 
    project = "T-LGLL",
    min.features = DEFAULT_MIN_FEATURES
  )) |>
    RenameCells(add.cell.id = row$gsm)
    
  obj$gsm <- row$gsm
  obj$subject_id <- row$subject_id
  obj$group <- factor(row$group, levels = GROUP_LEVELS)
  obj$title <- row$title
  
  Idents(obj) <- "gsm"
  
  # save to disk to avoid keeping all objects in memory
  saveRDS(obj, out_file)
  rm(counts, obj)
  gc()
}

# - 4 ~ merge Seurat objects in batches
sample_obj_files <- list.files(SAMPLE_OBJ_DIR, full.names = TRUE)

batch_size <- 4
batch_ids <- split(
  seq_along(sample_obj_files), 
  ceiling(seq_along(sample_obj_files) / batch_size)
)

batch_files <- character(length(batch_ids))

for (b in seq_along(batch_ids)) {
  batch_file <- here(BATCH_OBJ_DIR, paste0("batch_", b, ".rds"))
  if (file.exists(batch_file)) {
    batch_files[b] <- batch_file
    next
  }
  
  message(paste0("Building batch ", b, "/", length(batch_ids)))
  idx <- batch_ids[[b]]
  objs <- lapply(sample_obj_files[idx], readRDS)
  
  batch_obj <- Reduce(function(x, y) merge(x, y, merge.data = FALSE), objs)
  
  saveRDS(batch_obj, batch_file)
  batch_files[b] <- batch_file
  
  rm(objs, batch_obj)
  gc()
}

# - 5 ~ final Seurat object merge
message("Performing final merge...")

batch_objs <- lapply(batch_files, readRDS)
merged_obj <- Reduce(function(x, y) merge(x, y, merge.data = FALSE), batch_objs)

rm(batch_objs)
gc()

merged_obj <- JoinLayers(merged_obj, assay = "RNA") # join split count layers
Idents(merged_obj) <- "group" # group as active.ident going forward

saveRDS(merged_obj, MERGED_SEURAT_FILE)

message("Merged Seurat object saved to: ", MERGED_SEURAT_FILE)
message("02_seurat_prep.R complete.")

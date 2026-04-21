# ~~~ scripts/03_merge_and_subset.R ~~~
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(here)
  library(readr)
})

source(here("scripts", "00_config.R"), local = TRUE)
cfg <- get_config()

# - paths -
sample_obj_dir <- here(cfg$PROCESSED_DIR, "sample_objects")
batch_obj_dir <- here(cfg$PROCESSED_DIR, "batch_objects")

dir.create(batch_obj_dir, showWarnings = FALSE)

# - merge sample objects in batches -
sample_obj_files <- list.files(sample_obj_dir, full.names = TRUE)

batch_size <- 4
idx_vec <- seq_len(length(sample_obj_files))
batch_ids <- split(idx_vec, ceiling(idx_vec / batch_size))
batch_files <- character(length(batch_ids))

for (batch in seq_along(batch_ids)) {
  batch_file <- here(batch_obj_dir, paste0("batch_", batch, ".rds"))
  if (file.exists(batch_file)) {
    batch_files[batch] <- batch_file
    next
  }
  
  message(paste0("Building batch ", batch, "/", length(batch_ids)))
  
  objs <- lapply(sample_obj_files[batch_ids[[batch]]], readRDS)
  batch_obj <- Reduce(function(x, y) merge(x, y, merge.data = FALSE), objs)
  
  saveRDS(batch_obj, batch_file)
  batch_files[batch] <- batch_file
  
  rm(batch_file, objs, batch_obj)
  gc()
}

# - final merge - 
# *very memory intensive
message("Performing final merge...")

batch_objs <- lapply(batch_files, readRDS)
merged_obj <- Reduce(function(x, y) merge(x, y, merge.data = FALSE), batch_objs)
rm(batch_objs)
gc()

merged_obj <- JoinLayers(merged_obj, assay = "RNA") # join split count layers
Idents(merged_obj) <- "group" # "group" as active.ident going forward

saveRDS(merged_obj, cfg$MERGED_SEURAT_FILE)
message("Merged Seurat object saved to: ", cfg$MERGED_SEURAT_FILE)

# - create & save subsets -
baseline_obj <- subset(merged_obj, subset = group %in% cfg$BASELINE_GROUP)
treatment_obj <- subset(merged_obj, subset = group %in% cfg$TREATMENT_GROUP)
saveRDS(baseline_obj, cfg$BASELINE_SEURAT_FILE)
saveRDS(treatment_obj, cfg$TREATMENT_SEURAT_FILE)

message("Baseline Seurat object saved to: ", cfg$BASELINE_SEURAT_FILE)
message("Treatment Seurat object saved to: ", cfg$TREATMENT_SEURAT_FILE)
message("03_merge_and_subset.R complete.")
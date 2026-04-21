# ~~~ scripts/02_sample_qc.R ~~~
suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(dplyr)
  library(here)
  library(readr)
  library(stringr)
})

source(here("scripts", "00_config.R"), local = TRUE)
cfg <- get_config()

# - paths -
extracted_dir <- here(cfg$RAW_DIR, cfg$GEO_ACCESSION, "extracted")
sample_obj_dir <- here(cfg$PROCESSED_DIR, "sample_objects")
qc_summary_file <- here(cfg$PROCESSED_DIR, "sample_qc_summary.csv")

dir.create(sample_obj_dir, showWarnings = FALSE)

# - match files to manifest -
manifest <- read_csv(cfg$MANIFEST_FILE, show_col_types = FALSE) |>
  filter(data_type == "geneexpression")

count_files <- list.files(extracted_dir, full.names = TRUE)

file_tbl <- tibble(
  file_path = count_files,
  file_name = basename(count_files),
  gsm = str_extract(file_name, "GSM\\d+")
)

manifest <- manifest |> left_join(file_tbl, by = "gsm")

# - build & QC one sample at a time -
# * will likely take some time to complete (~40 mins on my machine)
qc_results <- list()

for (i in seq_len(nrow(manifest))) {
  # build object
  row <- manifest[i, ]
  out_file <- here(sample_obj_dir, paste0(row$gsm, ".rds"))
  if (file.exists(out_file)){
    next
  }
  
  message("Building object for ", row$gsm, " (", i, "/", nrow(manifest), ")")
  counts <- as.matrix(read.csv(row$file_path, row.names = 1)) |>
    Matrix(sparse = TRUE)
  
  obj <- suppressWarnings(CreateSeuratObject(
    counts = counts,
    project = "T-LGLL"
  )) |> RenameCells(add.cell.id = row$gsm)
  
  # attach sample metadata
  obj$gsm <- row$gsm
  obj$subject_id <- row$subject_id
  obj$group <- factor(row$group, levels = cfg$GROUP_LEVELS)
  obj$title <- row$title
  
  # QC metrics
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  cells_before <- ncol(obj)
  median_features_before <- median(obj$nFeature_RNA)
  median_counts_before <- median(obj$nCount_RNA)
  median_mt_before <- median(obj$percent.mt)
  
  # QC filtering (implement adaptive QC filtering if time allows)
  obj <- subset(
    obj,
    subset = 
      nFeature_RNA >= cfg$DEFAULT_MIN_FEATURES &
      nFeature_RNA <= cfg$DEFAULT_MAX_FEATURES &
      percent.mt <= cfg$DEFAULT_MAX_MT
  )
  
  # post-QC metrics
  cells_after <- ncol(obj)
  median_features_after <- median(obj$nFeature_RNA)
  median_counts_after <- median(obj$nCount_RNA)
  median_mt_after <- median(obj$percent.mt)
  
  # update active.idents -> saveRDS -> memory cleanup
  Idents(obj) <- "gsm"
  saveRDS(obj, out_file)
  rm(counts, obj)
  gc()
  
  # save QC results
  qc_results[[i]] <- tibble(
    gsm = row$gsm,
    subject_id = row$subject_id,
    group = row$group,
    cells_before = cells_before,
    cells_after = cells_after,
    pct_retained = 100 * (cells_after / cells_before),
    median_features_before = median_features_before,
    median_features_after = median_features_after,
    median_counts_before = median_counts_before,
    median_counts_after = median_counts_after,
    median_mt_before = median_mt_before,
    median_mt_after = median_mt_after
  )
}

# - save QC summary -
write_csv(bind_rows(qc_results), qc_summary_file)

message("Cleaned sample objects saved to: ", sample_obj_dir)
message("QC summary saved to: ", qc_summary_file)
message("02_sample_qc.R complete.")
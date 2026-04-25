# ~~~ scripts/04_seurat_analysis.R ~~~
suppressPackageStartupMessages({
  library(Seurat)
  library(cluster)
  library(dplyr)
  library(here)
  library(msigdbr)
  library(purrr)
  library(readr)
})

source(here("scripts", "00_config.R"), local = TRUE)
cfg <- get_config()

# - paths -
analysis_table_dir <- here(cfg$RESULTS_DIR, "analysis_tables")
dir.create(analysis_table_dir, showWarnings = FALSE)

baseline_summary_file <- here(analysis_table_dir, "baseline_summary.csv")
baseline_cluster_file <- here(analysis_table_dir, "baseline_cluster_summary.csv")
baseline_module_file <- here(analysis_table_dir, "baseline_group_module_scores.csv")

treatment_summary_file <- here(analysis_table_dir, "treatment_summary.csv")
treatment_cluster_file <- here(analysis_table_dir, "treatment_cluster_summary.csv")
treatment_module_file <- here(analysis_table_dir, "treatment_group_module_scores.csv")

# - parameters - 
param_grid <- as_tibble(cfg$PARAM_GRID)
max_npcs <- max(param_grid$npcs)

# - module score seed gene sets -
apoptosis_genes_seed <- msigdbr(species = "Homo sapiens", collection = "H") |>
  filter(gs_name == "HALLMARK_APOPTOSIS") |> 
  pull(gene_symbol) |> 
  unique()

cytotoxic_genes_seed <- msigdbr(
  species = "Homo sapiens", collection = "C5", subcollection = "GO:BP"
) |>
  filter(gs_name == "GOBP_T_CELL_MEDIATED_CYTOTOXICITY") |> 
  pull(gene_symbol) |> 
  unique()

# - helper functions -
preprocess <- function(obj) {
  obj |> 
    NormalizeData(
      normalization.method = "LogNormalize",
      scale.factor = 10000,
      verbose = FALSE
    ) |> 
    FindVariableFeatures(
      selection.method = "vst",
      nfeatures = 2000,
      verbose = FALSE
    ) |>
    ScaleData(verbose = FALSE) |>
    RunPCA(npcs = max_npcs, seed.use = cfg$SEED, verbose = FALSE)
}

compute_cluster_entropy <- function(cluster_labels) {
  counts <- table(cluster_labels) # count cells in cluster
  props <- counts / sum(counts) # convert to proportions
  -sum(props * log(props)) # compute Shannon entropy
}

compute_silhouette_mean <- function(obj, npcs, max_cells = 3000) {
  cluster_labels <- Idents(obj) # get cluster labels for each cell
  if (length(unique(cluster_labels)) < 2) { return (NA_real_) }
  pca_emb <- Embeddings(obj, reduction = "pca")[, seq_len(npcs), drop = FALSE]
  # ^ pulls the PCA coordinates for 'npcs' PCs
  
  # subsample to save time
  keep_idx <- seq_len(nrow(pca_emb))
  if (length(keep_idx) > max_cells) {
    set.seed(cfg$SEED)
    keep_idx <- sample(keep_idx, max_cells)
  }
  
  pca_emb <- pca_emb[keep_idx, , drop = FALSE]
  cl_sub <- droplevels(cluster_labels[keep_idx])
  if (length(unique(cl_sub)) < 2) { return (NA_real_) }
  
  distance <- dist(pca_emb) # compute pairwise distances
  sil <- silhouette(as.integer(cl_sub), distance) # calculate silhouette width
  mean(sil[, "sil_width"]) # return mean width across cells
}

add_program_scores <- function(obj) {
  apoptosis_genes <- intersect(apoptosis_genes_seed, rownames(obj))
  cytotoxic_genes <- intersect(cytotoxic_genes_seed, rownames(obj))
  
  if (length(apoptosis_genes) > 0) {
    obj <- AddModuleScore(
      obj, features = list(apoptosis_genes), 
      name = "ApoptosisScore", seed = cfg$SEED
    )
  } else {
    obj$ApoptosisScore1 <- NA_real_
  }
  
  if (length(cytotoxic_genes) > 0) {
    obj <- AddModuleScore(
      obj, features = list(cytotoxic_genes),
      name = "CytotoxicScore", seed = cfg$SEED
    )
  } else {
    obj$CytotoxicScore1 <- NA_real_
  }
  
  obj
}

run_setting <- function(obj, name, npcs, k.param, resolution) {
  obj <- FindNeighbors(
      obj, dims = seq_len(npcs), k.param = k.param, verbose = FALSE
    ) |> FindClusters(
      resolution = resolution, random.seed = cfg$SEED, verbose = FALSE
    ) 
  
  Idents(obj) <- "seurat_clusters"
  
  meta <- obj@meta.data |>
    as_tibble(rownames = "cell") |>
    mutate(cluster = as.character(Idents(obj)))
  
  summary_tbl <- tibble(
    dataset = name,
    npcs = npcs,
    k.param = k.param,
    resolution = resolution,
    n_cells = ncol(obj),
    n_clusters = n_distinct(Idents(obj)),
    largest_cluster_prop = max(table(Idents(obj))) / ncol(obj),
    cluster_entropy = compute_cluster_entropy(Idents(obj)),
    silhouette_mean = compute_silhouette_mean(obj, npcs)
  )
  
  cluster_summary_tbl <- meta |>
    count(group, cluster, name = "n_group_cluster") |>
    group_by(group) |>
    mutate(group_cluster_prop = n_group_cluster / sum(n_group_cluster)) |>
    ungroup() |>
    left_join(
      meta |>
        group_by(cluster) |>
        summarize(
          n_cells = n(),
          cluster_prop = n_cells / nrow(meta),
          mean_apoptosis = mean(ApoptosisScore1, na.rm = TRUE),
          mean_cytotoxic = mean(CytotoxicScore1, na.rm = TRUE),
          median_apoptosis = median(ApoptosisScore1, na.rm = TRUE),
          median_cytotoxic = median(CytotoxicScore1, na.rm = TRUE),
          .groups = "drop"
        ),
      by = "cluster"
    )
  
  list(
    summary = summary_tbl,
    cluster_summary = cluster_summary_tbl
  )
}

analyze_dataset <- function(obj_file, name) {
  message("Loading ", name, " object...")
  obj <- readRDS(obj_file)
  
  message("Preprocessing ", name, " object...")
  obj <- preprocess(obj) |> add_program_scores()
  
  group_module_tbl <- obj@meta.data |>
    as_tibble() |>
    group_by(group) |>
    summarise(
      mean_apoptosis = mean(ApoptosisScore1, na.rm = TRUE),
      mean_cytotoxic = mean(CytotoxicScore1, na.rm = TRUE),
      median_apoptosis = median(ApoptosisScore1, na.rm = TRUE),
      median_cytotoxic = median(CytotoxicScore1, na.rm = TRUE),
      .groups = "drop"
    )
  
  n_param_sets <- nrow(param_grid)
  results <- vector("list", n_param_sets)
  
  for (i in seq_len(n_param_sets)) {
    p <- param_grid[i, ]
    
    message(
      "Running ", name, " setting ", i, "/", n_param_sets,
      " [npcs=", p$npcs,
      ", k=", p$k.param,
      ", res=", p$resolution, "]"
    )
    
    results[[i]] <- run_setting(
      obj, name, p$npcs, p$k.param, p$resolution
    )
    
    gc()
  }
  
  list(
    summary = bind_rows(map(results, "summary")),
    cluster_summary = bind_rows(map(results, "cluster_summary")),
    group_modules = group_module_tbl
  )
}

# - run workflow -
baseline <- analyze_dataset(cfg$BASELINE_SEURAT_FILE, "baseline")
treatment <- analyze_dataset(cfg$TREATMENT_SEURAT_FILE, "treatment")

# - save results -
write_csv(baseline$summary, baseline_summary_file)
write_csv(baseline$cluster_summary, baseline_cluster_file)
write_csv(baseline$group_modules, baseline_module_file)

write_csv(treatment$summary, treatment_summary_file)
write_csv(treatment$cluster_summary, treatment_cluster_file)
write_csv(treatment$group_modules, treatment_module_file)

message("04_seurat_analysis.R complete.")
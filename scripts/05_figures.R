# ~~~ scripts/05_figures.R ~~~
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(here)
  library(msigdbr)
  library(patchwork)
  library(readr)
  library(stringr)
  library(tidyr)
})

source(here("scripts", "00_config.R"), local = TRUE)
cfg <- get_config()

# - paths -
qc_summary_file <- here(cfg$PROCESSED_DIR, "sample_qc_summary.csv")

baseline_summary_file <- here(cfg$TABLES_DIR, "baseline_summary.csv")
baseline_cluster_file <- here(cfg$TABLES_DIR, "baseline_cluster_summary.csv")
baseline_module_file <- here(cfg$TABLES_DIR, "baseline_group_module_scores.csv")

treatment_summary_file <- here(cfg$TABLES_DIR, "treatment_summary.csv")
treatment_cluster_file <- here(cfg$TABLES_DIR, "treatment_cluster_summary.csv")
treatment_module_file <- here(cfg$TABLES_DIR, "treatment_group_module_scores.csv")

# - read summary outputs -
qc_summary <- read_csv(qc_summary_file, show_col_types = FALSE)

baseline_summary <- read_csv(baseline_summary_file, show_col_types = FALSE)
baseline_cluster <- read_csv(baseline_cluster_file, show_col_types = FALSE)
baseline_modules <- read_csv(baseline_module_file, show_col_types = FALSE)

treatment_summary <- read_csv(treatment_summary_file, show_col_types = FALSE)
treatment_cluster <- read_csv(treatment_cluster_file, show_col_types = FALSE)
treatment_modules <- read_csv(treatment_module_file, show_col_types = FALSE)

# - choose settings -
pick_best_setting <- function(summary_tbl) {
  summary_tbl |> arrange(desc(silhouette_mean)) |> slice(1)
}

best_baseline <- pick_best_setting(baseline_summary)
best_treatment <- pick_best_setting(treatment_summary)

# - plot helpers -
plot_metric_lines <- function(summary_tbl, metric, y_label) {
  summary_tbl |>
    mutate(resolution = factor(resolution)) |>
    ggplot(aes(x = npcs, y = .data[[metric]], 
               color = resolution, group = resolution)) +
    geom_line() + geom_point(size = 2) +
    facet_wrap(~ k.param, nrow = 3, 
               labeller = as_labeller(function(x) paste("k.param =", x))) +
    scale_x_continuous(breaks = c(10, 20, 30), limits = c(10, 30)) +
    labs(
      x = "Number of PCs",
      y = y_label,
      color = "Resolution"
    ) +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "grey85", color = "black"),
      strip.text = element_text(face = "bold"),
      panel.grid.major = element_line(linewidth = 0.3, color = "grey80"),
      panel.grid.minor = element_blank()
    )
}

plot_group_modules <- function(group_module_tbl) {
  group_module_tbl |>
    select(group, mean_apoptosis, mean_cytotoxic) |>
    pivot_longer(
      cols = c(mean_apoptosis, mean_cytotoxic),
      names_to = "score_type",
      values_to = "score"
    ) |>
    mutate(
      group = str_to_title(group),
      score_type = recode_values(
        score_type,
        "mean_apoptosis" ~ "Apoptosis",
        "mean_cytotoxic" ~ "Cytotoxicity"
      )
    ) |>
    ggplot(aes(x = group, y = score, fill = score_type)) +
    geom_col(position = "dodge") +
    labs(
      x = "Group",
      y = "Mean Module Score",
      fill = "Pathway"
    ) +
    theme_bw()
}

plot_cluster_modules <- function(cluster_tbl, chosen_setting) {
  cluster_tbl |>
    filter(
      k.param == chosen_setting$k.param,
      npcs == chosen_setting$npcs,
      resolution == chosen_setting$resolution
    ) |>
    distinct(cluster, mean_apoptosis, mean_cytotoxic) |>
    pivot_longer(
      cols = c(mean_apoptosis, mean_cytotoxic),
      names_to = "score_type",
      values_to = "score"
    ) |>
    mutate(
      cluster = as.factor(cluster),
      score_type = recode_values(
        score_type,
        "mean_apoptosis" ~ "Apoptosis",
        "mean_cytotoxic" ~ "Cytotoxicity"
      )
    ) |>
    ggplot(aes(x = cluster, y = score, fill = score_type)) +
    geom_col(position = "dodge") +
    labs(
      caption = paste0(
        "PCs = ", chosen_setting$npcs,
        ", k = ", chosen_setting$k.param,
        ", res = ", chosen_setting$resolution
      ),
      x = "Cluster",
      y = "Mean Module Score",
      fill = "Pathway"
    ) +
    theme_bw() 
}

plot_cluster_group_comp <- function(cluster_tbl, chosen_setting) {
  cluster_tbl |>
    filter(
      k.param == chosen_setting$k.param,
      npcs == chosen_setting$npcs,
      resolution == chosen_setting$resolution
    ) |>
    mutate(cluster = as.factor(cluster), group = str_to_title(group)) |>
    ggplot(aes(x = cluster, y = group_cluster_prop, fill = group)) +
    geom_col(position = "fill") +
    labs(
      caption = paste0(
        "PCs = ", chosen_setting$npcs,
        ", k = ", chosen_setting$k.param,
        ", res = ", chosen_setting$resolution
      ),
      x = "Cluster",
      y = "Within-group Proportion",
      fill = "Group"
    ) +
    theme_bw()
}

# - seurat helpers for TSNE -
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

preprocess <- function(obj, max_npcs) {
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

add_program_scores <- function(obj) {
  apoptosis_genes <- intersect(apoptosis_genes_seed, rownames(obj))
  cytotoxic_genes <- intersect(cytotoxic_genes_seed, rownames(obj))
  
  obj <- AddModuleScore(
    object = obj, features = list(apoptosis_genes),
    name = "ApoptosisScore", seed = cfg$SEED
  )
  
  obj <- AddModuleScore(
    object = obj, features = list(cytotoxic_genes),
    name = "CytotoxicScore", seed = cfg$SEED
  )
  
  obj
}

run_final_setting <- function(obj_file, chosen_setting) {
  readRDS(obj_file) |>
    preprocess(chosen_setting$npcs) |>
    add_program_scores() |>
    FindNeighbors(
      dims = seq_len(chosen_setting$npcs), k.param = chosen_setting$k.param,
      verbose = FALSE
    ) |>
    FindClusters(
      resolution = chosen_setting$resolution, random.seed = cfg$SEED,
      verbose = FALSE
    ) |>
    RunTSNE(
      dims = seq_len(chosen_setting$npcs), perplexity = cfg$DEFAULT_PERPLEXITY,
      seed.use = cfg$SEED, verbose = FALSE
    )
}

# - build QC figures -
qc_violin <- qc_summary |>
  select(
    gsm, subject_id, group,
    median_features_before, median_features_after,
    median_counts_before, median_counts_after,
    median_mt_before, median_mt_after
  ) |>
  pivot_longer(
    cols = -c(gsm, subject_id, group),
    names_to = "metric_stage",
    values_to = "value"
  ) |>
  mutate(
    metric = case_when(
      grepl("^median_features_", metric_stage) ~ "Detected Features",
      grepl("^median_counts_", metric_stage) ~ "UMI Counts",
      grepl("^median_mt_", metric_stage) ~ "Percent Mitochondrial"
    ),
    stage = case_when(
      grepl("_before$", metric_stage) ~ "Before QC",
      grepl("_after$", metric_stage) ~ "After QC"
    ),
    group = factor(
      str_to_title(group), 
      levels = c("Healthy", "Pretreatment", "Posttreatment")
    ),
    metric = factor(metric, levels = unique(metric)),
    stage = factor(stage, levels = unique(stage))
  ) |>
  select(gsm, subject_id, group, metric, stage, value) |>
  ggplot(aes(x = stage, y = value, fill = stage)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_jitter(width = 0.08, size = 1.5, alpha = 0.8) +
  facet_wrap(group ~ metric, scales = "free_y") +
  labs(
    x = NULL,
    y = "Per-sample Median Value",
  ) +
  theme_bw() +
  theme(legend.position = "none")

qc_retention <- qc_summary |>
  mutate(
    gsm = reorder(gsm, pct_retained),
    group = str_to_title(group)
  ) |>
  ggplot(aes(x = gsm, y = pct_retained, fill = group)) +
  geom_col() +
  geom_hline(yintercept = 80, linetype = "dashed") +
  coord_flip() +
  labs(
    x = "Sample",
    y = "Cells Retained after QC (%)",
    fill = "Group"
  ) +
  theme_bw()

ggsave(filename = here(cfg$FIGURES_DIR, "qc_violin.png"),
       plot = qc_violin, width = 8, height = 8, dpi = 300)
ggsave(filename = here(cfg$FIGURES_DIR, "qc_retention.png"),
       plot = qc_retention, width = 8, height = 8, dpi = 300)

# - build silhouette figures -
baseline_silhouetteF <- plot_metric_lines(baseline_summary, "silhouette_mean", 
                                          "Mean Silhouette")
treatment_silhouetteF <- plot_metric_lines(treatment_summary, "silhouette_mean", 
                                           "Mean Silhouette")
ggsave(filename = here(cfg$FIGURES_DIR, "baseline_silhouette.png"),
       plot = baseline_silhouetteF, width = 5, height = 5, dpi = 300)
ggsave(filename = here(cfg$FIGURES_DIR, "treatment_silhouette.png"),
       plot = treatment_silhouetteF, width = 5, height = 5, dpi = 300)

# - build module figures -
baseline_moduleF <- 
  plot_group_modules(baseline_modules) / 
  plot_cluster_modules(baseline_cluster, best_baseline) + 
  plot_layout(guides = "collect")
treatment_moduleF <-
  plot_group_modules(treatment_modules) /
  plot_cluster_modules(treatment_cluster, best_treatment) +
  plot_layout(guides = "collect")
ggsave(filename = here(cfg$FIGURES_DIR, "baseline_module.png"),
       plot = baseline_moduleF, width = 6, height = 8, dpi = 300)
ggsave(filename = here(cfg$FIGURES_DIR, "treatment_module.png"),
       plot = treatment_moduleF, width = 6, height = 8, dpi = 300)

# - build appendix figures: number of clusters -
baseline_clusterF <- plot_metric_lines(baseline_summary, "n_clusters",
                                       "Number of Clusters")
treatment_clusterF <- plot_metric_lines(treatment_summary, "n_clusters",
                                        "Number of Clusters")
ggsave(filename = here(cfg$FIGURES_DIR, "baseline_cluster.png"),
       plot = baseline_clusterF, width = 5, height = 5, dpi = 300)
ggsave(filename = here(cfg$FIGURES_DIR, "treatment_cluster.png"),
       plot = treatment_clusterF, width = 5, height = 5, dpi = 300)

# - build appendix figures: group composition by cluster -
baseline_cluster_compF <- plot_cluster_group_comp(baseline_cluster,
                                                  best_baseline)
treatment_cluster_compF <- plot_cluster_group_comp(treatment_cluster,
                                                   best_treatment)
ggsave(filename = here(cfg$FIGURES_DIR, "baseline_cluster_group_comp.png"),
       plot = baseline_cluster_compF, width = 6, height = 4, dpi = 300)
ggsave(filename = here(cfg$FIGURES_DIR, "treatment_cluster_group_comp.png"),
       plot = treatment_cluster_compF, width = 6, height = 4, dpi = 300)

# - build t-SNE for final visualization -
baseline_obj_final <- run_final_setting(cfg$BASELINE_SEURAT_FILE, 
                                        best_baseline)
baseline_obj_final$group <- str_to_title(baseline_obj_final$group)
treatment_obj_final <- run_final_setting(cfg$TREATMENT_SEURAT_FILE, 
                                         best_treatment)
treatment_obj_final$group <- str_to_title(treatment_obj_final$group)

baseline_cluster_plot <-
  DimPlot(baseline_obj_final, label = TRUE, raster = TRUE) + 
  labs(title = NULL)
baseline_group_plot <-
  DimPlot(baseline_obj_final, group.by = "group", raster = TRUE) + 
  labs(title = NULL)
baseline_apop_plot <-
  FeaturePlot(baseline_obj_final, features = "ApoptosisScore1", raster = TRUE) + 
  labs(title = NULL)
baseline_cyto_plot <-
  FeaturePlot(baseline_obj_final, features = "CytotoxicScore1", raster = TRUE) +
  labs(title = NULL)

treatment_cluster_plot <-
  DimPlot(treatment_obj_final, label = TRUE, raster = TRUE) +
  labs(title = NULL)
treatment_group_plot <-
  DimPlot(treatment_obj_final, group.by = "group", raster = TRUE) +
  labs(title = NULL)
treatment_apop_plot <-
  FeaturePlot(treatment_obj_final, features = "ApoptosisScore1", raster = TRUE) +
  labs(title = NULL)
treatment_cyto_plot <-
  FeaturePlot(treatment_obj_final, features = "CytotoxicScore1", raster = TRUE) +
  labs(title = NULL)
  
baseline_row <-
  baseline_cluster_plot + baseline_group_plot +
  baseline_apop_plot + baseline_cyto_plot +
  plot_layout(widths = c(1, 1, 1, 1))

treatment_row <-
  treatment_cluster_plot + treatment_group_plot +
  treatment_apop_plot + treatment_cyto_plot +
  plot_layout(widths = c(1, 1, 1, 1))

full_tsne <- baseline_row / treatment_row + plot_layout(heights = c(1, 1))
ggsave(filename = here(cfg$FIGURES_DIR, "full_tsne.png"),
       plot = full_tsne, width = 21, height = 12, dpi = 300)

message("05_figures.R complete.")
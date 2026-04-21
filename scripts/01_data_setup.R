# ~~~ scripts/01_data_setup.R ~~~
suppressPackageStartupMessages({
  library(GEOquery)
  library(dplyr)
  library(here)
  library(readr)
  library(stringr)
  library(tibble)
})

source(here("scripts", "00_config.R"), local = TRUE)
cfg <- get_config()

# - download GEO metadata & build manifest -
gsm_list <- suppressMessages(GSMList(
  getGEO(cfg$GEO_ACCESSION, GSEMatrix = FALSE)
))

# extract metadata from each GSM and return a vector of tibbles
manifest <- lapply(names(gsm_list), function(gsm_id) {
  meta <- Meta(gsm_list[[gsm_id]])
  
  tibble(
    gsm = gsm_id,
    title = meta$title,
    characteristics = paste(meta$characteristics_ch1, collapse = " | ")
  )
}) |> bind_rows() |> # bind all tibbles into one tbl_df
  mutate( # create manifest attributes
    subject_id = str_extract(title, "UPN\\d+|Healthydonor\\d+"),
    data_type = str_extract(title, "geneexpression|VDJ"),
    group = str_to_lower(
      str_extract(title, "Healthy|pretreatment|posttreatment")),
    is_healthy = group == "healthy",
    is_patient = group != "healthy",
    use_baseline = data_type == "geneexpression" & 
      group %in% cfg$BASELINE_GROUP,
    use_treatment = data_type == "geneexpression" & 
      group %in% cfg$TREATMENT_GROUP
  ) |>
  select( # rearrange columns 
    gsm,
    subject_id,
    data_type,
    group,
    is_healthy,
    is_patient,
    use_baseline,
    use_treatment,
    title,
    characteristics
  )

write_csv(manifest, cfg$MANIFEST_FILE)
message("Manifest saved to: ", cfg$MANIFEST_FILE)

# - download and unpack supplementary files (raw data) -
suppressMessages(getGEOSuppFiles(cfg$GEO_ACCESSION, baseDir = cfg$RAW_DIR))
# ^ creates supp_dir by default

supp_dir <- here(cfg$RAW_DIR, cfg$GEO_ACCESSION)
supp_tar_file <- here(supp_dir, str_c(cfg$GEO_ACCESSION, "_RAW.tar"))

if (file.exists(supp_tar_file)) {
  untar(supp_tar_file, exdir = here(supp_dir, "extracted"))
}

message("Supplementary files saved to: ", supp_dir)
message("01_data_setup.R complete.")
# ~~~ scripts/01_data_setup.R ~~~
suppressPackageStartupMessages({
  library(GEOquery)
  library(dplyr)
  library(here)
  library(readr)
  library(stringr)
  library(tibble)
})

source(here("scripts", "00_config.R"))

# - 1 ~ download GEO metadata & build manifest
gsm_list <- suppressMessages(GSMList(getGEO(GEO_ACCESSION, GSEMatrix = FALSE)))

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
    use_baseline = data_type == "geneexpression" & group %in% c(
      "healthy", "pretreatment"),
    use_treatment = data_type == "geneexpression" & group %in% c(
      "pretreatment", "posttreatment")
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

write_csv(manifest, MANIFEST_FILE)

# - 2 ~ download supplementary files (raw data)

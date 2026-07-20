
setwd("/Users/bohe/Desktop/Metabolomics/Raw_data")

dir.create(
  "MTBLS12932_LCMS_RAW/positive",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "MTBLS12932_LCMS_RAW/negative",
  recursive = TRUE,
  showWarnings = FALSE
)


options(timeout = 86400)
for (i in seq_len(nrow(lcms_files))) {
  
  output_dir <- file.path(
    "MTBLS12932_LCMS_RAW",
    lcms_files$Polarity[i]
  )
  
  output_file <- file.path(
    output_dir,
    lcms_files$File_name[i]
  )
  
  download.file(
    url = lcms_files$Download_URL[i],
    destfile = output_file,
    mode = "wb"
  )
}
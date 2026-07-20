
library(xcms)
library(MsExperiment)
library(Spectra)
library(BiocParallel)




project_dir <- "/Users/bohe/Desktop/Metabolomics"

raw_data_dir <- file.path(
  project_dir,
  "Raw_data",
  "MTBLS12932_LCMS_RAW"
)

processing_dir <- file.path(
  project_dir,
  "Data_processing",
  "MTBLS12932_XCMS"
)



dir.create(
  file.path(processing_dir, "positive"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(processing_dir, "negative"),
  recursive = TRUE,
  showWarnings = FALSE
)




positive_mzml <- list.files(
  path = file.path(raw_data_dir, "mzML_positive"),
  pattern = "\\.mzML$",
  full.names = TRUE,
  ignore.case = TRUE
)

negative_mzml <- list.files(
  path = file.path(raw_data_dir, "mzML_negative"),
  pattern = "\\.mzML$",
  full.names = TRUE,
  ignore.case = TRUE
)




positive_sample_names <- tools::file_path_sans_ext(
  basename(positive_mzml)
)

positive_sample_data <- data.frame(
  sample_name = positive_sample_names,
  sample_group = rep("positive", length(positive_mzml)),
  row.names = basename(positive_mzml),
  stringsAsFactors = FALSE
)



positive_data <- readMsExperiment(
  spectraFiles = positive_mzml,
  sampleData = positive_sample_data
)



positive_centwave_param <- CentWaveParam(
  ppm = 10,
  peakwidth = c(5, 30),
  snthresh = 10,
  prefilter = c(3, 1000),
  mzCenterFun = "wMean",
  integrate = 1,
  mzdiff = -0.001,
  noise = 0
)



positive_peaks <- findChromPeaks(
  object = positive_data,
  param = positive_centwave_param,
  msLevel = 1L,
  chunkSize = 1L,
  BPPARAM = SerialParam()
)



saveRDS(
  positive_peaks,
  file = file.path(
    processing_dir,
    "positive",
    "MTBLS12932_positive_01_peak_detection.rds"
  )
)




positive_obiwarp_param <- ObiwarpParam(
  binSize = 0.5,
  response = 1,
  distFun = "cor_opt",
  gapInit = 0.3,
  gapExtend = 2.4,
  factorDiag = 2,
  factorGap = 1,
  localAlignment = FALSE
)



positive_rt <- adjustRtime(
  object = positive_peaks,
  param = positive_obiwarp_param,
  chunkSize = 1L,
  BPPARAM = SerialParam()
)


saveRDS(
  positive_rt,
  file = file.path(
    processing_dir,
    "positive",
    "MTBLS12932_positive_02_retention_time_adjusted.rds"
  )
)




positive_group_param <- PeakDensityParam(
  sampleGroups = rep(1, length(positive_mzml)),
  minFraction = 0.25,
  minSamples = 1,
  bw = 5,
  binSize = 0.01
)



positive_grouped <- groupChromPeaks(
  object = positive_rt,
  param = positive_group_param
)




saveRDS(
  positive_grouped,
  file = file.path(
    processing_dir,
    "positive",
    "MTBLS12932_positive_03_peak_grouping.rds"
  )
)



positive_fill_param <- ChromPeakAreaParam()



positive_filled <- fillChromPeaks(
  object = positive_grouped,
  param = positive_fill_param,
  chunkSize = 1L,
  BPPARAM = SerialParam()
)



saveRDS(
  positive_filled,
  file = file.path(
    processing_dir,
    "positive",
    "MTBLS12932_positive_04_filled_peaks.rds"
  )
)



positive_feature_information <- as.data.frame(
  featureDefinitions(positive_filled)
)

positive_feature_information <- positive_feature_information[
  ,
  !vapply(
    positive_feature_information,
    is.list,
    logical(1)
  ),
  drop = FALSE
]

positive_feature_information$Feature_ID <- rownames(
  positive_feature_information
)

positive_feature_information <- positive_feature_information[
  ,
  c(
    "Feature_ID",
    setdiff(
      colnames(positive_feature_information),
      "Feature_ID"
    )
  ),
  drop = FALSE
]


positive_peak_area <- featureValues(
  object = positive_filled,
  value = "into",
  method = "maxint",
  filled = TRUE
)

positive_peak_area <- as.data.frame(
  positive_peak_area,
  check.names = FALSE
)

colnames(positive_peak_area) <- positive_sample_names

positive_peak_area$Feature_ID <- rownames(
  positive_peak_area
)

positive_peak_area <- positive_peak_area[
  ,
  c(
    "Feature_ID",
    positive_sample_names
  ),
  drop = FALSE
]




positive_feature_table <- cbind(
  positive_feature_information,
  positive_peak_area[
    ,
    positive_sample_names,
    drop = FALSE
  ]
)




write.table(
  positive_feature_information,
  file = file.path(
    processing_dir,
    "positive",
    "MTBLS12932_positive_feature_information.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  positive_peak_area,
  file = file.path(
    processing_dir,
    "positive",
    "MTBLS12932_positive_peak_area_matrix.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  positive_feature_table,
  file = file.path(
    processing_dir,
    "positive",
    "MTBLS12932_positive_feature_table.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)






negative_sample_names <- tools::file_path_sans_ext(
  basename(negative_mzml)
)

negative_sample_data <- data.frame(
  sample_name = negative_sample_names,
  sample_group = rep("negative", length(negative_mzml)),
  row.names = basename(negative_mzml),
  stringsAsFactors = FALSE
)



negative_data <- readMsExperiment(
  spectraFiles = negative_mzml,
  sampleData = negative_sample_data
)



negative_centwave_param <- CentWaveParam(
  ppm = 10,
  peakwidth = c(5, 30),
  snthresh = 10,
  prefilter = c(3, 1000),
  mzCenterFun = "wMean",
  integrate = 1,
  mzdiff = -0.001,
  noise = 0
)




negative_peaks <- findChromPeaks(
  object = negative_data,
  param = negative_centwave_param,
  msLevel = 1L,
  chunkSize = 1L,
  BPPARAM = SerialParam()
)


saveRDS(
  negative_peaks,
  file = file.path(
    processing_dir,
    "negative",
    "MTBLS12932_negative_01_peak_detection.rds"
  )
)



negative_obiwarp_param <- ObiwarpParam(
  binSize = 0.5,
  response = 1,
  distFun = "cor_opt",
  gapInit = 0.3,
  gapExtend = 2.4,
  factorDiag = 2,
  factorGap = 1,
  localAlignment = FALSE
)



negative_rt <- adjustRtime(
  object = negative_peaks,
  param = negative_obiwarp_param,
  chunkSize = 1L,
  BPPARAM = SerialParam()
)


saveRDS(
  negative_rt,
  file = file.path(
    processing_dir,
    "negative",
    "MTBLS12932_negative_02_retention_time_adjusted.rds"
  )
)



negative_group_param <- PeakDensityParam(
  sampleGroups = rep(1, length(negative_mzml)),
  minFraction = 0.25,
  minSamples = 1,
  bw = 5,
  binSize = 0.01
)



negative_grouped <- groupChromPeaks(
  object = negative_rt,
  param = negative_group_param
)


saveRDS(
  negative_grouped,
  file = file.path(
    processing_dir,
    "negative",
    "MTBLS12932_negative_03_peak_grouping.rds"
  )
)


negative_fill_param <- ChromPeakAreaParam()



negative_filled <- fillChromPeaks(
  object = negative_grouped,
  param = negative_fill_param,
  chunkSize = 1L,
  BPPARAM = SerialParam()
)


saveRDS(
  negative_filled,
  file = file.path(
    processing_dir,
    "negative",
    "MTBLS12932_negative_04_filled_peaks.rds"
  )
)



negative_feature_information <- as.data.frame(
  featureDefinitions(negative_filled)
)



negative_feature_information <- negative_feature_information[
  ,
  !vapply(
    negative_feature_information,
    is.list,
    logical(1)
  ),
  drop = FALSE
]

negative_feature_information$Feature_ID <- rownames(
  negative_feature_information
)

negative_feature_information <- negative_feature_information[
  ,
  c(
    "Feature_ID",
    setdiff(
      colnames(negative_feature_information),
      "Feature_ID"
    )
  ),
  drop = FALSE
]




negative_peak_area <- featureValues(
  object = negative_filled,
  value = "into",
  method = "maxint",
  filled = TRUE
)

negative_peak_area <- as.data.frame(
  negative_peak_area,
  check.names = FALSE
)

colnames(negative_peak_area) <- negative_sample_names

negative_peak_area$Feature_ID <- rownames(
  negative_peak_area
)

negative_peak_area <- negative_peak_area[
  ,
  c(
    "Feature_ID",
    negative_sample_names
  ),
  drop = FALSE
]





negative_feature_table <- cbind(
  negative_feature_information,
  negative_peak_area[
    ,
    negative_sample_names,
    drop = FALSE
  ]
)





write.table(
  negative_feature_information,
  file = file.path(
    processing_dir,
    "negative",
    "MTBLS12932_negative_feature_information.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  negative_peak_area,
  file = file.path(
    processing_dir,
    "negative",
    "MTBLS12932_negative_peak_area_matrix.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  negative_feature_table,
  file = file.path(
    processing_dir,
    "negative",
    "MTBLS12932_negative_feature_table.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)








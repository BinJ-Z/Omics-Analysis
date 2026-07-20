# ============================================================
# LC-MS Data Processing
#
# Workflow
# 1. Feature filtering
# 2. Missing value imputation
# 3. TIC normalization
# 4. Log2 transformation
#
# Output:
# Final log2 feature tables
# ============================================================


project_dir <- "/Users/bohe/Desktop/Metabolomics"

xcms_dir <- file.path(
  project_dir,
  "Data_processing",
  "MTBLS12932_XCMS"
)

output_dir <- file.path(
  project_dir,
  "Data_processing",
  "Final_feature_tables"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


process_feature_table <- function(
    input_file,
    sample_prefix,
    output_file
){
  
  feature_data <- read.delim(
    input_file,
    check.names = FALSE
  )
  
  sample_columns <- grep(
    paste0("^", sample_prefix),
    colnames(feature_data),
    value = TRUE
  )
  
  feature_matrix <- as.matrix(
    feature_data[
      ,
      sample_columns,
      drop = FALSE
    ]
  )
  
  storage.mode(feature_matrix) <- "numeric"
  
  feature_matrix[feature_matrix == 0] <- NA
  
  
  keep <- rowSums(
    !is.na(feature_matrix)
  ) >= ceiling(
    ncol(feature_matrix) * 0.5
  )
  
  feature_data <- feature_data[
    keep,
    ,
    drop = FALSE
  ]
  
  feature_matrix <- feature_matrix[
    keep,
    ,
    drop = FALSE
  ]
  
  
  for(i in seq_len(nrow(feature_matrix))){
    
    minimum_value <- min(
      feature_matrix[i, ],
      na.rm = TRUE
    )
    
    feature_matrix[
      i,
      is.na(feature_matrix[i, ])
    ] <- minimum_value / 2
    
  }
  
  
  sample_tic <- colSums(
    feature_matrix
  )
  
  feature_matrix <- sweep(
    feature_matrix,
    2,
    sample_tic,
    "/"
  )
  
  feature_matrix <- feature_matrix *
    median(sample_tic)
  
  
  feature_matrix <- log2(
    feature_matrix
  )
  
  
  feature_data[
    ,
    sample_columns
  ] <- feature_matrix
  
  
  write.table(
    feature_data,
    file = output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
}


process_feature_table(
  
  input_file = file.path(
    xcms_dir,
    "positive",
    "MTBLS12932_positive_feature_table.tsv"
  ),
  
  sample_prefix = "POS_",
  
  output_file = file.path(
    output_dir,
    "MTBLS12932_positive_feature_table_log2.tsv"
  )
  
)


process_feature_table(
  
  input_file = file.path(
    xcms_dir,
    "negative",
    "MTBLS12932_negative_feature_table.tsv"
  ),
  
  sample_prefix = "NEG_",
  
  output_file = file.path(
    output_dir,
    "MTBLS12932_negative_feature_table_log2.tsv"
  )
  
)
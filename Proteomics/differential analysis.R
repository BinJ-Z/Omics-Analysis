
# ============================================================
# Proteomics differential protein analysis
# ============================================================


project_dir <- "/Users/bohe/Desktop/Proteomics"

abundance_file <- file.path(
  project_dir,
  "Quantification",
  "protein_abundance.tsv"
)

output_dir <- file.path(
  project_dir,
  "Differential_Proteins"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Read protein abundance table

protein_data <- read.delim(
  abundance_file,
  comment.char = "#",
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# Find OpenMS abundance columns

abundance_columns <- grep(
  "^abundance_sample",
  colnames(protein_data),
  value = TRUE
)


# Define sample IDs according to OpenMS sample order

sample_names <- c(
  "221005_SR_CRC_23",
  "221005_SR_CRC_24",
  "221005_SR_CRC_25",
  "221005_SR_CRC_26",
  "221005_SR_HC_23",
  "221005_SR_HC_24",
  "221005_SR_HC_25",
  "221005_SR_HC_26"
)


# Create protein abundance matrix

abundance_matrix <- as.matrix(
  protein_data[
    ,
    abundance_columns,
    drop = FALSE
  ]
)

storage.mode(abundance_matrix) <- "numeric"

colnames(abundance_matrix) <- sample_names

rownames(abundance_matrix) <- make.unique(
  as.character(protein_data$protein)
)


# Define groups

crc_samples <- c(
  "221005_SR_CRC_23",
  "221005_SR_CRC_24",
  "221005_SR_CRC_25",
  "221005_SR_CRC_26"
)

hc_samples <- c(
  "221005_SR_HC_23",
  "221005_SR_HC_24",
  "221005_SR_HC_25",
  "221005_SR_HC_26"
)


# Remove decoy proteins

keep_target <- !grepl(
  "DECOY_",
  rownames(abundance_matrix),
  ignore.case = TRUE
)

abundance_matrix <- abundance_matrix[
  keep_target,
  ,
  drop = FALSE
]


# Convert zero and negative values to NA

abundance_matrix[
  abundance_matrix <= 0
] <- NA


# Keep proteins detected in at least 50% of both groups

crc_detected <- rowSums(
  !is.na(
    abundance_matrix[
      ,
      crc_samples,
      drop = FALSE
    ]
  )
)

hc_detected <- rowSums(
  !is.na(
    abundance_matrix[
      ,
      hc_samples,
      drop = FALSE
    ]
  )
)

keep_proteins <- (
  crc_detected >= ceiling(length(crc_samples) / 2)
) & (
  hc_detected >= ceiling(length(hc_samples) / 2)
)

abundance_matrix <- abundance_matrix[
  keep_proteins,
  ,
  drop = FALSE
]


# Log2 transformation

log2_matrix <- log2(
  abundance_matrix
)


# Column-wise median normalization

sample_medians <- apply(
  log2_matrix,
  2,
  median,
  na.rm = TRUE
)

normalized_matrix <- sweep(
  log2_matrix,
  2,
  sample_medians,
  FUN = "-"
)


# Differential analysis

results <- data.frame(
  Protein = rownames(normalized_matrix),
  mean_CRC = NA_real_,
  mean_HC = NA_real_,
  log2FC = NA_real_,
  FC = NA_real_,
  p_value = NA_real_,
  stringsAsFactors = FALSE
)


for (i in seq_len(nrow(normalized_matrix))) {
  
  crc_values <- as.numeric(
    normalized_matrix[
      i,
      crc_samples
    ]
  )
  
  hc_values <- as.numeric(
    normalized_matrix[
      i,
      hc_samples
    ]
  )
  
  crc_values <- crc_values[
    !is.na(crc_values)
  ]
  
  hc_values <- hc_values[
    !is.na(hc_values)
  ]
  
  results$mean_CRC[i] <- mean(
    crc_values
  )
  
  results$mean_HC[i] <- mean(
    hc_values
  )
  
  results$log2FC[i] <- (
    results$mean_CRC[i] -
      results$mean_HC[i]
  )
  
  results$FC[i] <- 2^results$log2FC[i]
  
  if (
    length(crc_values) >= 2 &&
    length(hc_values) >= 2
  ) {
    
    results$p_value[i] <- tryCatch(
      t.test(
        crc_values,
        hc_values,
        var.equal = FALSE
      )$p.value,
      error = function(e) NA_real_
    )
  }
}


# BH correction

results$FDR <- p.adjust(
  results$p_value,
  method = "BH"
)


# Sort all results by FDR

results <- results[
  order(
    results$FDR,
    na.last = TRUE
  ),
  ,
  drop = FALSE
]


# Significant proteins: FDR < 0.05 and FC >= 2 or FC <= 0.5

significant_results <- results[
  !is.na(results$FDR) &
    results$FDR < 0.05 &
    (
      results$FC >= 2 |
        results$FC <= 0.5
    ),
  ,
  drop = FALSE
]


# Save normalized protein abundance matrix

write.table(
  data.frame(
    Protein = rownames(normalized_matrix),
    normalized_matrix,
    check.names = FALSE
  ),
  file = file.path(
    output_dir,
    "protein_abundance_log2_normalized.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# Save all differential analysis results

write.table(
  results,
  file = file.path(
    output_dir,
    "differential_proteins_all.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# Save significant differential proteins

write.table(
  significant_results,
  file = file.path(
    output_dir,
    "differential_proteins_significant.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


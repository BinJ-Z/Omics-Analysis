# Proteomics differential protein analysis

project_dir <- "/Users/bohe/Desktop/Proteomics"
protein_data <- read.delim(file.path(project_dir, "Quantification", "protein_abundance.tsv"), comment.char = "#", check.names = FALSE)
output_dir <- file.path(project_dir, "Differential_Proteins")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sample_names <- c(
  "221005_SR_CRC_23", "221005_SR_CRC_24",
  "221005_SR_CRC_25", "221005_SR_CRC_26",
  "221005_SR_HC_23", "221005_SR_HC_24",
  "221005_SR_HC_25", "221005_SR_HC_26"
)

crc_samples <- sample_names[grepl("_CRC_", sample_names)]
hc_samples <- sample_names[grepl("_HC_", sample_names)]

abundance_columns <- grep("^abundance_sample", colnames(protein_data), value = TRUE)
abundance_matrix <- as.matrix(protein_data[, abundance_columns, drop = FALSE])
storage.mode(abundance_matrix) <- "numeric"
colnames(abundance_matrix) <- sample_names
rownames(abundance_matrix) <- make.unique(as.character(protein_data$protein))

abundance_matrix <- abundance_matrix[
  !grepl("DECOY_", rownames(abundance_matrix), ignore.case = TRUE),
  ,
  drop = FALSE
]

abundance_matrix[abundance_matrix <= 0] <- NA

keep <- rowSums(!is.na(abundance_matrix[, crc_samples, drop = FALSE])) >= 2 &
  rowSums(!is.na(abundance_matrix[, hc_samples, drop = FALSE])) >= 2

abundance_matrix <- abundance_matrix[keep, , drop = FALSE]

log2_matrix <- log2(abundance_matrix)
normalized_matrix <- sweep(
  log2_matrix,
  2,
  apply(log2_matrix, 2, median, na.rm = TRUE),
  FUN = "-"
)

results <- do.call(
  rbind,
  lapply(seq_len(nrow(normalized_matrix)), function(i) {
    
    crc <- na.omit(as.numeric(normalized_matrix[i, crc_samples]))
    hc <- na.omit(as.numeric(normalized_matrix[i, hc_samples]))
    
    mean_crc <- mean(crc)
    mean_hc <- mean(hc)
    log2fc <- mean_crc - mean_hc
    
    p_value <- if (length(crc) >= 2 && length(hc) >= 2) {
      tryCatch(t.test(crc, hc, var.equal = FALSE)$p.value, error = function(e) NA_real_)
    } else {
      NA_real_
    }
    
    data.frame(
      Protein = rownames(normalized_matrix)[i],
      mean_CRC = mean_crc,
      mean_HC = mean_hc,
      log2FC = log2fc,
      FC = 2^log2fc,
      p_value = p_value
    )
  })
)

results$FDR <- p.adjust(results$p_value, method = "BH")
results <- results[order(results$FDR, na.last = TRUE), ]

significant_results <- results[
  !is.na(results$FDR) &
    results$FDR < 0.05 &
    (results$FC >= 2 | results$FC <= 0.5),
]

write.table(
  results,
  file.path(output_dir, "differential_proteins_all.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  significant_results,
  file.path(output_dir, "differential_proteins_significant.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
library(tximport)
library(DESeq2)


setwd("/Users/bohe/Desktop/CRC ssRNA")


dir.create("results/13_DESeq2", recursive = TRUE, showWarnings = FALSE)


metadata <- read.delim("results/07_Differential_Taxa/metadata.tsv", row.names = 1, check.names = FALSE)

metadata$Group <- factor(metadata$Group, levels = c("healthy", "CRC"))

files <- file.path("results/11_Salmon_Quant", rownames(metadata), "quant.sf")

names(files) <- rownames(metadata)

txip <- tximport(files, type = "salmon", txOut = TRUE)


write.table(
  txip$counts,
  file = "results/12_Transcript_Counts/transcript_counts.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)


write.table(
  txip$abundance,
  file = "results/12_Transcript_Counts/transcript_TPM.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)


dds <- DESeqDataSetFromTximport(txip, colData = metadata, design = ~ Group)

dds <- dds[rowSums(counts(dds)) >= 10, ]

dds <- DESeq(dds)

res <- results(dds, contrast = c("Group", "CRC", "healthy"), alpha = 0.05)

res <- res[order(res$padj), ]

res_1 <- as.data.frame(res)
res_1$Transcript <- rownames(res_1)

res_1 <- res_1[, c("Transcript", "baseMean", "log2FoldChange",
                   "lfcSE", "stat", "pvalue", "padj")]


write.table(res_1, file = "results/13_DESeq2/DETs_all.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


significant <- res_1[
  !is.na(res_1$padj) &
    res_1$padj < 0.05 &
    abs(res_1$log2FoldChange) >= 1,
]


write.table(significant, file = "results/13_DESeq2/DETs_significant.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


normalized_counts <- counts(dds, normalized = TRUE)


write.table(normalized_counts, file = "results/12_Transcript_Counts/transcript_normalized_counts.tsv", sep = "\t", quote = FALSE, col.names = NA)


saveRDS(dds, file = "results/13_DESeq2/DESeq2_dds.rds")
cat("DESeq2 analysis completed\n")
cat("Tested transcripts:", nrow(res_1), "\n")
cat("Significant DETs:", nrow(significant), "\n")
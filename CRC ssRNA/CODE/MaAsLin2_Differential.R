library(Maaslin2)

table(metadata$Group)


abundance <- read.delim("/Users/bohe/Desktop/CRC ssRNA/results/07_Differential_Taxa/metaphlan_merged.tsv")
metadata <- read.delim(
 "/Users/bohe/Desktop/CRC ssRNA/metadata.tsv",
  row.names = 1,
)

taxon_column <- colnames(abundance)[1]
species <- abundance[
  grepl("\\|s__", abundance[[taxon_column]]) &
    !grepl("\\|t__", abundance[[taxon_column]]),
  ,
  drop = FALSE
]

rownames(species) <- make.unique(
  sub(".*\\|s__", "", species[[taxon_column]])
)
colnames(species) <- sub("_profile$", "", colnames(species))


species[[taxon_column]] <- NULL
species[] <- lapply(
  species,
  function(x) as.numeric(as.character(x))
)

species <- species[
  rowSums(species > 0, na.rm = TRUE) >= 2,
  ,
  drop = FALSE
]



species <- as.data.frame(t(species))

common_samples <- intersect(rownames(metadata), rownames(species))

metadata <- metadata[common_samples, , drop = FALSE]
species <- species[common_samples, , drop = FALSE]

metadata$Group <- factor(
  metadata$Group,
  levels = c("Healthy", "CRC")
)

Maaslin2(
  input_data = species,
  input_metadata = metadata,
  output = "/Users/bohe/Desktop/CRC ssRNA/results/07_Differential_Taxa/MaAsLin2_results",
  fixed_effects = "Group",
  reference = "Group,Healthy",
  normalization = "NONE",
  transform = "LOG",
  analysis_method = "LM",
  standardize = FALSE,
  plot_heatmap = TRUE
)
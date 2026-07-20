LC-MS Untargeted Metabolomics Analysis


1.Data source
SOURCE:Metabolights
STUDY ID:MTBLS12932
URL:https://www.ebi.ac.uk/metabolights/editor/MTBLS12932/overview

2.Software
R 4.5
ThermoRawFileParser
xcms
MSnbase
ggplot2

3.Workflow
*Convert Data RAW files to mzML using( ThermoRawFileParser.)
*Perform chromatographic peak detection  (XCMS).
*Fill missing peaks.
*Filter features detected in fewer than 50% of samples.
*Impute missing values.
*Perform TIC normalization.
*Apply log2 transformation.
*Perform differential feature analysis （FDR < 0.05 and |log2FC| ≥ 1）



4.Output Files

For both positive and negative ionization modes, the following files are generated:

* `*_all_features.tsv` — Differential analysis results for all features.
* `*_significant_features.tsv` — Significant differential features.
* `*_differential_summary.tsv` — Summary of differential analysis results.


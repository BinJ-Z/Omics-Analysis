#set environment first
#conda activate rachis-qiime2-2026.4
#bash 16s_data.sh
#cd "/Users/bohe/Desktop/CRC 16S"

printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > manifest.tsv

for f in raw_data/*_1.fastq.gz
do
    sample=$(basename "$f" _1.fastq.gz)

    printf "%s\t%s\t%s\n" \
    "$sample" \
    "$PWD/raw_data/${sample}_1.fastq.gz" \
    "$PWD/raw_data/${sample}_2.fastq.gz" \
    >> manifest.tsv
done


qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.tsv \
  --output-path paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2


#below code seperate to run 

qiime cutadapt trim-paired \
  --i-demultiplexed-sequences paired-end-demux.qza \
  --p-front-f CCTACGGGNGGCWGCAG \
  --p-front-r GGACTACHVGGGTWTCTAAT \
  --p-discard-untrimmed \
  --p-cores 0 \
  --o-trimmed-sequences paired-end-trimmed.qza \
  --verbose



qiime dada2 denoise-paired \
  --i-demultiplexed-seqs paired-end-trimmed.qza \
  --p-trunc-len-f 220 \
  --p-trunc-len-r 220 \
  --p-max-ee-f 2 \
  --p-max-ee-r 2 \
  --p-n-threads 4 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats stats.qza \
  --o-base-transition-stats base-transition-stats.qza



qiime metadata tabulate \
  --m-input-file stats.qza \
  --o-visualization stats.qzv



qiime rescript get-silva-data \
  --p-version 138.2 \
  --p-target SSURef_NR99 \
  --o-silva-sequences silva-seqs.qza \
  --o-silva-taxonomy silva-tax.qza


qiime rescript reverse-transcribe \
  --i-rna-sequences silva-seqs.qza \
  --o-dna-sequences silva-dna-seqs.qza


qiime feature-classifier extract-reads \
  --i-sequences silva-dna-seqs.qza \
  --p-f-primer CCTACGGGNGGCWGCAG \
  --p-r-primer GGACTACHVGGGTWTCTAAT \
  --p-min-length 100 \
  --p-max-length 500 \
  --o-reads silva-341f-806r-seqs.qza


qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-341f-806r-seqs.qza \
  --i-reference-taxonomy silva-tax.qza \
  --o-classifier silva-341f-806r-classifier.qza


qiime feature-classifier classify-sklearn \
  --i-classifier silva-341f-806r-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza


qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv



qiime metadata tabulate \
  --m-input-file sample-metadata.tsv \
  --o-visualization sample-metadata.qzv


qiime taxa filter-table \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast,eukaryota \
  --o-filtered-table table-filtered.qza


qiime taxa filter-seqs \
  --i-sequences rep-seqs.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast,eukaryota \
  --o-filtered-sequences rep-seqs-filtered.qza



qiime taxa barplot \
  --i-table table-filtered.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization taxa-bar-plots-filtered.qzv


qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs-filtered.qza \
  --o-alignment aligned-rep-seqs-filtered.qza \
  --o-masked-alignment masked-aligned-rep-seqs-filtered.qza \
  --o-tree unrooted-tree-filtered.qza \
  --o-rooted-tree rooted-tree-filtered.qza


qiime feature-table summarize \
  --i-table table-filtered.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-feature-frequencies feature-frequencies-filtered.qza \
  --o-sample-frequencies sample-frequencies-filtered.qza \
  --o-summary table-filtered.qzv

# below code can run together

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree-filtered.qza \
  --i-table table-filtered.qza \
  --p-sampling-depth 45000 \
  --m-metadata-file sample-metadata.tsv \
  --output-dir core-metrics-results-filtered


qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results-filtered/shannon_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization core-metrics-results-filtered/shannon-group-significance.qzv


qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results-filtered/observed_features_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization core-metrics-results-filtered/observed-features-group-significance.qzv


qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results-filtered/faith_pd_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization core-metrics-results-filtered/faith-pd-group-significance.qzv


qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results-filtered/evenness_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization core-metrics-results-filtered/evenness-group-significance.qzv


qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results-filtered/bray_curtis_distance_matrix.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column Group \
  --p-method permanova \
  --p-pairwise \
  --o-visualization core-metrics-results-filtered/bray-curtis-permanova.qzv


qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results-filtered/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column Group \
  --p-method permanova \
  --p-pairwise \
  --o-visualization core-metrics-results-filtered/weighted-unifrac-permanova.qzv


qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results-filtered/bray_curtis_distance_matrix.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column Group \
  --p-method permdisp \
  --o-visualization core-metrics-results-filtered/bray-curtis-permdisp.qzv


qiime taxa collapse \
  --i-table table-filtered.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table table-filtered-genus.qza

qiime feature-table summarize \
  --i-table table-filtered-genus.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-feature-frequencies feature-frequencies-genus.qza \
  --o-sample-frequencies sample-frequencies-genus.qza \
  --o-summary table-filtered-genus.qzv

qiime composition ancombc \
  --i-table table-filtered-genus.qza \
  --m-metadata-file sample-metadata.tsv \
  --p-formula Group \
  --o-differentials ancombc-genus.qza

qiime composition tabulate \
  --i-data ancombc-genus.qza \
  --o-visualization ancombc-genus.qzv


qiime composition da-barplot \
  --i-data ancombc-genus.qza \
  --p-significance-threshold 0.05 \
  --o-visualization ancombc-genus-barplot.qzv

qiime composition da-barplot \
  --i-data ancombc-genus.qza \
  --p-significance-threshold 0.05 \
  --o-visualization ancombc-genus-barplot.qzv

mkdir -p result

mv *.qza *.qzv result/

mv core-metrics-results-filtered result/

cd result

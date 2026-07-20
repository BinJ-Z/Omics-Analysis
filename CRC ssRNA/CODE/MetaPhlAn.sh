
conda activate metaRNA
mkdir -p results/06_MetaPhlAn/profiles
mkdir -p results/06_MetaPhlAn/bowtie2


for r1 in results/05_ncRNA_Removal/clean_reads/*_non_rRNA_tRNA_fwd.fq
do
  sample=$(basename "$r1" _non_rRNA_tRNA_fwd.fq)

  metaphlan \
    "results/05_ncRNA_Removal/clean_reads/${sample}_non_rRNA_tRNA_fwd.fq,results/05_ncRNA_Removal/clean_reads/${sample}_non_rRNA_tRNA_rev.fq" \
    --input_type fastq \
    --bowtie2db database/MetaPhlAn \
    --nproc 8 \
    --bowtie2out "results/06_MetaPhlAn/bowtie2/${sample}.bowtie2.bz2" \
    -o "results/06_MetaPhlAn/profiles/${sample}_profile.tsv"
done



mkdir -p results/07_Differential_Taxa

merge_metaphlan_tables.py \
  results/06_MetaPhlAn/profiles/*_profile.tsv \
  > results/07_Differential_Taxa/metaphlan_merged.tsv



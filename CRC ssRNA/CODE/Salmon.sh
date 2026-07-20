conda activate metaRNA
cd "/Users/bohe/Desktop/CRC ssRNA"

mkdir -p results/10_Salmon/index

salmon index \
  -t results/08_rnaSPAdes/assembly/transcripts.fasta \
  -i results/10_Salmon/index/common_transcriptome_index \
  -p 8 



mkdir -p results/11_Salmon_Quant


for r1 in results/05_ncRNA_Removal/clean_reads/*_non_rRNA_tRNA_fwd.fq
do
  sample=$(basename "$r1" _non_rRNA_tRNA_fwd.fq)
  mkdir -p "results/11_Salmon_Quant/${sample}"

  salmon quant \
    -i results/10_Salmon/index/common_transcriptome_index \
    -l A \
    -1 "results/05_ncRNA_Removal/clean_reads/${sample}_non_rRNA_tRNA_fwd.fq" \
    -2 "results/05_ncRNA_Removal/clean_reads/${sample}_non_rRNA_tRNA_rev.fq" \
    -p 8 \
    --meta \
    --validateMappings \
    --seqBias \
    --gcBias \
    -o "results/11_Salmon_Quant/${sample}" 
done



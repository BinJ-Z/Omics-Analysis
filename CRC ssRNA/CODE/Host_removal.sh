
mkdir -p results/04_HostRemoval/clean_reads


for f in results/03_fastp/clean_reads/*_1.clean.fastq.gz
do
  sample=$(basename "$f" _1.clean.fastq.gz)

  bowtie2 \
    -x database/GRCh38/GRCh38_noalt_as \
    -1 "results/03_fastp/clean_reads/${sample}_1.clean.fastq.gz" \
    -2 "results/03_fastp/clean_reads/${sample}_2.clean.fastq.gz" \
    --sensitive \
    --threads 8 \
    --un-conc-gz "results/04_HostRemoval/clean_reads/${sample}_nonhost.fastq.gz" \
    -S /dev/null
done


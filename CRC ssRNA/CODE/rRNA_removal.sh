mkdir -p results/05_ncRNA_Removal/clean_reads
mkdir -p results/05_ncRNA_Removal/removed_reads
mkdir -p results/05_ncRNA_Removal/logs
mkdir -p results/05_ncRNA_Removal/work

for r1 in results/04_HostRemoval/clean_reads/*_nonhost.fastq.1.gz
do
  sample=$(basename "$r1" _nonhost.fastq.1.gz)

  sortmerna \
    --ref "database/SortMeRNA/rRNA_databases/smr_v4.3_default_db.fasta" \
    --ref "database/SortMeRNA/tRNA_databases/microbial_tRNA.fasta" \
    --reads "results/04_HostRemoval/clean_reads/${sample}_nonhost.fastq.1.gz" \
    --reads "results/04_HostRemoval/clean_reads/${sample}_nonhost.fastq.2.gz" \
    --paired_in \
    --out2 \
    --fastx \
    --aligned "results/05_ncRNA_Removal/removed_reads/${sample}_rRNA_tRNA" \
    --other "results/05_ncRNA_Removal/clean_reads/${sample}_non_rRNA_tRNA" \
    --threads 8 \
    --workdir "results/05_ncRNA_Removal/work/${sample}" \
    > "results/05_ncRNA_Removal/logs/${sample}.sortmerna.log" 2>&1
done




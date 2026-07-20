conda activate metaRNA
cd "/Users/bohe/Desktop/CRC ssRNA"


mkdir -p results/01_FastQC
fastqc \
  ssRNA_data/*.fastq.gz \
  --threads 8 \
  --outdir results/01_FastQC


mkdir -p results/02_MultiQC
multiqc \
  results/01_FastQC \
  --outdir results/02_MultiQC


mkdir -p results/03_fastp/clean_reads
mkdir -p results/03_fastp/reports

for f in ssRNA_data/*_1.fastq.gz
do
  sample=$(basename "$f" _1.fastq.gz)

  fastp \
    --in1 "ssRNA_data/${sample}_1.fastq.gz" \
    --in2 "ssRNA_data/${sample}_2.fastq.gz" \
    --out1 "results/03_fastp/clean_reads/${sample}_1.clean.fastq.gz" \
    --out2 "results/03_fastp/clean_reads/${sample}_2.clean.fastq.gz" \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --unqualified_percent_limit 40 \
    --length_required 50 \
    --thread 8 \
    --html "results/03_fastp/reports/${sample}.fastp.html" \
    --json "results/03_fastp/reports/${sample}.fastp.json"
done


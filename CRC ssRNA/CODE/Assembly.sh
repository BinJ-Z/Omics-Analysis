
conda activate metaRNA
cd "/Users/bohe/Desktop/CRC ssRNA"


conda activate metaRNA
cd "/Users/bohe/Desktop/CRC ssRNA"
mkdir -p results/08_rnaSPAdes/input
mkdir -p results/08_rnaSPAdes/assembly

cat results/05_ncRNA_Removal/clean_reads/*_non_rRNA_tRNA_fwd.fq \
  > results/08_rnaSPAdes/input/all_samples_fwd.fq

cat results/05_ncRNA_Removal/clean_reads/*_non_rRNA_tRNA_rev.fq \
  > results/08_rnaSPAdes/input/all_samples_rev.fq

rnaspades.py \
  -1 results/08_rnaSPAdes/input/all_samples_fwd.fq \
  -2 results/08_rnaSPAdes/input/all_samples_rev.fq \
  -t 8 \
  -m 18 \
  -o results/08_rnaSPAdes/assembly



mkdir -p results/09_rnaQUAST


rnaQUAST.py \
  --transcripts results/08_rnaSPAdes/assembly/transcripts.fasta \
  --labels common_transcriptome \
  --threads 8 \
  --output_dir results/09_rnaQUAST/report


1. Software and Database Requirements

1.1 conda
FastQC 0.12.1
MultiQC 1.35
fastp 1.1.0
Bowtie2 2.5.5
SortMeRNA 4.3.7
MetaPhlAn 4.0.6
SPAdes (RNA mode) 4.2.0
rnaQUAST 2.3.1
Salmon 2.3.3


1.2 R package
R 4.3.3
tximport 1.38.1
DESeq2 1.42.1
MaAsLin2 1.16.0

1.3 Database
Human reference genome GRCh38_noalt_as
SortMeRNA rRNA database smr_v4.3_default_db.fasta (v4.3)
Microbial tRNA database microbial_tRNA.fasta
MetaPhlAn marker database mpa_vJan25_CHOCOPhlAnSGB_202503



2. Workflow
Raw data prepare（Raw_data.sh）
Quality assessment and Read preprocessing (Quality_control.sh)
Host sequence removal (Host_removal.sh)
rRNA and tRNA removal (rRNA_removal.sh)
Microbial taxonomic profiling (MetaPhlAn.sh)
Differential microbial abundance analysis (MaAsLin2_Differential)
Transcriptome assembly(Assembly.sh)
Transcriptome indexing and quantification (Salmon.sh)
transcript expression analysis (Expression_DESeq.sh)



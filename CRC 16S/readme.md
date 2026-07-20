16S rRNA Gene Sequencing Analysis Pipeline

1.Software
 
- QIIME 2 (2026.4)
- RESCRIPt
- Cutadapt
- DADA2
- SILVA database v138.2 

2..Workflow

- Prepare manifest file and import sequencing data
- Trim primer sequences and perform DADA2 denoising
- Train the SILVA classifier and assign taxonomy
- Filter non-bacterial sequences and generate taxonomic profiles
- Construct the phylogenetic tree
- Perform alpha and beta diversity analyses
- Perform differential abundance analysis using ANCOM-BC
- Organize all results

1.Software
Database：PRIDE Archive - PRoteomics IDEntifications Database
Study ID: PXD078610
URL: https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/07/PXD078610/
openMS download：https://openms.readthedocs.io/en/latest/about/installation/installation-on-macos.html
Proteomes Homo sapiens (Human) database：https://www.uniprot.org/proteomes/UP000005640  （Download only reviewed (Swiss-Prot) canonical proteins (20,416)）
Proteomics contaminant databases： https://zenodo.org/records/15115102?utm_source=chatgpt.com （crap_gpm.fasta:）
Java: https://www.oracle.com/java/technologies/downloads/#jdk26-mac
MSFragger: http://msfragger-upgrader.nesvilab.org/upgrader/


2. Workflow

- Generate target–decoy protein database
- Perform peptide identification using MSFragger
- Annotate peptide-to-protein relationships
- Apply 1% false discovery rate (FDR) filtering
- Detect MS1 peptide features
- Align retention times across samples
- Link corresponding peptide features
- Perform label-free protein quantification
- Normalize protein abundances
- Perform differential protein analysis (Welch's t-test, BH correction; FDR < 0.05 and FC ≥ 2 or ≤ 0.5)

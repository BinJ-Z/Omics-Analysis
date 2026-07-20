conda activate openms

# 1. Paths


PROJECT="/Users/bohe/Desktop/Proteomics"

MZML_DIR="${PROJECT}/mzML"
DATABASE_DIR="${PROJECT}/Database"
IDENTIFICATION_DIR="${PROJECT}/Identification"
QUANTIFICATION_DIR="${PROJECT}/Quantification"
DIFFERENTIAL_DIR="${PROJECT}/Differential_Proteins"

MSFRAGGER_JAR="${PROJECT}/MSFragger-4.4.1/MSFragger-4.4.1.jar"

TARGET_FASTA="${DATABASE_DIR}/human_swissprot_cRAP.fasta"

TARGET_DECOY_FASTA="${DATABASE_DIR}/human_swissprot_cRAP_target_decoy.fasta"


MSFRAGGER_OUT="${IDENTIFICATION_DIR}/MSFragger"
MZBIN_OUT="${IDENTIFICATION_DIR}/mzBIN"
INDEXED_OUT="${IDENTIFICATION_DIR}/PeptideIndexer"
FDR_OUT="${IDENTIFICATION_DIR}/FDR_1percent"

FEATURE_OUT="${QUANTIFICATION_DIR}/FeatureFinderIdentification"
ALIGNED_OUT="${QUANTIFICATION_DIR}/Aligned_features"

CONSENSUS_FILE="${QUANTIFICATION_DIR}/proteomics_consensus.consensusXML"

PROTEIN_CSV="${QUANTIFICATION_DIR}/protein_abundance.csv"
PROTEIN_TSV="${QUANTIFICATION_DIR}/protein_abundance.tsv"

PEPTIDE_CSV="${QUANTIFICATION_DIR}/peptide_abundance.csv"
PEPTIDE_TSV="${QUANTIFICATION_DIR}/peptide_abundance.tsv"


mkdir -p "${MSFRAGGER_OUT}"
mkdir -p "${MZBIN_OUT}"
mkdir -p "${INDEXED_OUT}"
mkdir -p "${FDR_OUT}"
mkdir -p "${FEATURE_OUT}"
mkdir -p "${ALIGNED_OUT}"
mkdir -p "${QUANTIFICATION_DIR}"
mkdir -p "${DIFFERENTIAL_DIR}"


# 2. Generate target-decoy FASTA database


if [ ! -f "${TARGET_DECOY_FASTA}" ]; then

    echo "Generating target-decoy FASTA database"

    DecoyDatabase \
        -in "${TARGET_FASTA}" \
        -out "${TARGET_DECOY_FASTA}" \
        -method reverse \
        -enzyme Trypsin \
        -decoy_string DECOY_ \
        -decoy_string_position prefix

fi



# 3. MSFragger database search


for mzml in "${MZML_DIR}"/*.mzML
do

    sample=$(basename "${mzml}" .mzML)

    idxml="${MSFRAGGER_OUT}/${sample}.idXML"

    if [ ! -f "${idxml}" ]; then

        echo "MSFraggerAdapter: ${sample}"

        MSFraggerAdapter \
            -license yes \
            -java_executable java \
            -java_heapmemory 6500 \
            -executable "${MSFRAGGER_JAR}" \
            -in "${mzml}" \
            -out "${idxml}" \
            -database "${TARGET_DECOY_FASTA}" \
            -tolerance:precursor_mass_tolerance_lower 10 \
            -tolerance:precursor_mass_tolerance_upper 10 \
            -tolerance:precursor_mass_unit ppm \
            -tolerance:fragment_mass_tolerance 20 \
            -tolerance:fragment_mass_unit ppm \
            -digest:search_enzyme_name Trypsin \
            -digest:allowed_missed_cleavage 2 \
            -statmod:unimod "Carbamidomethyl (C)" \
            -varmod:unimod "Oxidation (M)" \
            -threads 4

    fi

done


# Move MSFragger mzBIN files

mv "${MZML_DIR}"/*_uncalibrated.mzBIN \
   "${MZBIN_OUT}/" 2>/dev/null || true


# 4. PeptideIndexer

for idxml in "${MSFRAGGER_OUT}"/*.idXML
do

    sample=$(basename "${idxml}" .idXML)

    indexed_idxml="${INDEXED_OUT}/${sample}_indexed.idXML"

    if [ ! -f "${indexed_idxml}" ]; then

        echo "PeptideIndexer: ${sample}"

        PeptideIndexer \
            -in "${idxml}" \
            -out "${indexed_idxml}" \
            -fasta "${TARGET_DECOY_FASTA}" \
            -enzyme:name Trypsin \
            -decoy_string DECOY_ \
            -decoy_string_position prefix \
            -missing_decoy_action warn \
            -threads 4

    fi

done



# 5. PSM-level 1% FDR


for idxml in "${INDEXED_OUT}"/*_indexed.idXML
do

    sample=$(basename "${idxml}" _indexed.idXML)

    fdr_idxml="${FDR_OUT}/${sample}_fdr_1percent.idXML"

    if [ ! -f "${fdr_idxml}" ]; then

        echo "FalseDiscoveryRate: ${sample}"

        FalseDiscoveryRate \
            -in "${idxml}" \
            -out "${fdr_idxml}" \
            -PSM true \
            -peptide false \
            -protein false \
            -FDR:PSM 0.01 \
            -threads 4

    fi

done


# 6. MS1 feature detection and peptide quantification
#

for mzml in "${MZML_DIR}"/*.mzML
do

    sample=$(basename "${mzml}" .mzML)

    fdr_idxml="${FDR_OUT}/${sample}_fdr_1percent.idXML"

    featurexml="${FEATURE_OUT}/${sample}.featureXML"

    if [ ! -f "${featurexml}" ]; then

        echo "FeatureFinderIdentification: ${sample}"

        FeatureFinderIdentification \
            -in "${mzml}" \
            -id "${fdr_idxml}" \
            -out "${featurexml}" \
            -threads 4

    fi

done


# 7. Retention-time alignment

FEATURE_FILES=()
ALIGNED_FILES=()

for featurexml in "${FEATURE_OUT}"/*.featureXML
do

    sample=$(basename "${featurexml}" .featureXML)

    FEATURE_FILES+=("${featurexml}")

    ALIGNED_FILES+=(
        "${ALIGNED_OUT}/${sample}_aligned.featureXML"
    )

done


ALIGNMENT_REQUIRED=false

for aligned_file in "${ALIGNED_FILES[@]}"
do

    if [ ! -f "${aligned_file}" ]; then
        ALIGNMENT_REQUIRED=true
        break
    fi

done


if [ "${ALIGNMENT_REQUIRED}" = true ]; then

    echo "MapAlignerIdentification"

    MapAlignerIdentification \
        -in "${FEATURE_FILES[@]}" \
        -out "${ALIGNED_FILES[@]}" \
        -threads 4

fi


# 8. Link corresponding features across samples

if [ ! -f "${CONSENSUS_FILE}" ]; then

    echo "FeatureLinkerUnlabeledKD"

    FeatureLinkerUnlabeledKD \
        -in "${ALIGNED_FILES[@]}" \
        -out "${CONSENSUS_FILE}" \
        -threads 4

fi



# 9. Protein and peptide quantification


if [ ! -f "${PROTEIN_CSV}" ]; then

    echo "ProteinQuantifier"

    ProteinQuantifier \
        -in "${CONSENSUS_FILE}" \
        -out "${PROTEIN_CSV}" \
        -peptide_out "${PEPTIDE_CSV}" \
        -method top \
        -top:N 3 \
        -top:aggregate median \
        -threads 4

fi

cp "${PROTEIN_CSV}" "${PROTEIN_TSV}"
cp "${PEPTIDE_CSV}" "${PEPTIDE_TSV}"


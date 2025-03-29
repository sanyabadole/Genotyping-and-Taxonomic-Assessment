#!/bin/bash

# ===============================================================================
# Bacterial Genome Analysis Pipeline v2.0
# 
# A comprehensive pipeline for bacterial genome analysis including:
# - Genus and species-level genotyping
# - Strain typing using MLST
# - Genome quality assessment
# - Taxonomic classification and contamination detection
# ===============================================================================

set -e  # Exit on error

# ========================== CONFIGURATION OPTIONS ==============================

# Default input and output paths
ASSEMBLY_DIR="${1:-../Data/assemblies}"
GFF_DIR="${2:-../Data/gffs}"
OUTPUT_DIR="${3:-$HOME/Bacterial_Genome_Analysis}"
THREADS="${4:-8}"           # Default to 8 threads
KRAKEN_DB_DIR="${5}"        # Optional Kraken database directory
MIN_ANI_THRESHOLD=95        # Minimum ANI value for species assignment
CONDA_ENV_NAME="bacterial_analysis_env"

# ========================== INITIALIZATION ====================================

# Function for displaying usage information
usage() {
    echo "Usage: $0 [ASSEMBLY_DIR] [GFF_DIR] [OUTPUT_DIR] [THREADS] [KRAKEN_DB_DIR]"
    echo ""
    echo "Arguments:"
    echo "  ASSEMBLY_DIR   Directory containing genome assembly FASTA files (default: ../Data/assemblies)"
    echo "  GFF_DIR        Directory containing GFF annotation files (default: ../Data/gffs)"
    echo "  OUTPUT_DIR     Directory for output files (default: $HOME/Bacterial_Genome_Analysis)"
    echo "  THREADS        Number of CPU threads to use (default: 8)"
    echo "  KRAKEN_DB_DIR  Directory containing Kraken2 database (optional)"
    echo ""
    echo "Example: $0 ./my_assemblies ./my_annotations ./results 16 /path/to/kraken_db"
    exit 1
}

# Print banner
echo "======================================================================"
echo "           Bacterial Genome Analysis Pipeline v2.0                    "
echo "======================================================================"
echo ""

# Create output directory structure
mkdir -p "${OUTPUT_DIR}/logs"
mkdir -p "${OUTPUT_DIR}/genotyping/genus_level"
mkdir -p "${OUTPUT_DIR}/genotyping/species_level/mlst"
mkdir -p "${OUTPUT_DIR}/genotyping/species_level/fastani"
mkdir -p "${OUTPUT_DIR}/quality_assessment"
mkdir -p "${OUTPUT_DIR}/taxonomy"

LOG_DIR="${OUTPUT_DIR}/logs"

# Set up logging
log_file="${LOG_DIR}/pipeline_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$log_file") 2>&1

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting pipeline execution"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Output directory: ${OUTPUT_DIR}"

# Verify input files exist
if [ -z "$(ls -A ${ASSEMBLY_DIR}/*.fasta 2>/dev/null)" ]; then
    echo "[ERROR] No FASTA files found in ${ASSEMBLY_DIR}"
    echo "Please check that your assembly files exist and have .fasta extension"
    usage
fi

# ========================== CONDA ENVIRONMENT SETUP ==========================

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Setting up conda environment"

if ! command -v conda &> /dev/null; then
    echo "[ERROR] Conda is not installed or not in your PATH"
    echo "Please install Conda or add it to your PATH before running this script."
    exit 1
fi

# Create conda environment if it doesn't exist
if ! conda env list | grep -q "${CONDA_ENV_NAME}"; then
    echo "Creating conda environment: ${CONDA_ENV_NAME}"
    conda create -n ${CONDA_ENV_NAME} -c conda-forge -c bioconda mlst checkm-genome fastani kraken2 taxonkit -y
fi

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate ${CONDA_ENV_NAME}

# Verify environment activation
if [[ ! "$CONDA_DEFAULT_ENV" == "${CONDA_ENV_NAME}" ]]; then
    echo "[ERROR] Failed to activate conda environment: ${CONDA_ENV_NAME}"
    exit 1
fi

echo "Successfully activated conda environment: ${CONDA_DEFAULT_ENV}"

# ========================== UTILITY FUNCTIONS ================================

# Function to track execution time and memory usage
run_with_metrics() {
    local command="$1"
    local description="$2"
    local log_file="$3"
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting: $description"
    
    # Run the command with time measurement
    start_time=$(date +%s)
    bash -c "gtime -f '\nRuntime: %E\nMax Memory (KB): %M' bash -c \"$command\"" 2>> "$log_file"
    exit_status=$?
    end_time=$(date +%s)
    
    # Calculate duration
    duration=$((end_time - start_time))
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed: $description (Duration: ${duration}s, Exit status: $exit_status)"
    
    return $exit_status
}

# Function to handle errors
handle_error() {
    local error_message="$1"
    echo "[ERROR] $error_message" >&2
    echo "[ERROR] $error_message" >> "${LOG_DIR}/errors.log"
}

# Function to get taxonomic hierarchy for common bacterial genera
get_taxonomy_hierarchy() {
    local genus=$1
    # Convert to lowercase for case-insensitive matching using tr
    genus=$(echo "$genus" | tr '[:upper:]' '[:lower:]')
    
    case "${genus}" in
        "campylobacter")
            echo "Epsilonproteobacteria|Campylobacterales|Campylobacteraceae|Campylobacter"
            ;;
        "salmonella")
            echo "Gammaproteobacteria|Enterobacterales|Enterobacteriaceae|Salmonella"
            ;;
        "escherichia")
            echo "Gammaproteobacteria|Enterobacterales|Enterobacteriaceae|Escherichia"
            ;;
        "listeria")
            echo "Bacilli|Bacillales|Listeriaceae|Listeria"
            ;;
        "vibrio")
            echo "Gammaproteobacteria|Vibrionales|Vibrionaceae|Vibrio"
            ;;
        "staphylococcus")
            echo "Bacilli|Bacillales|Staphylococcaceae|Staphylococcus"
            ;;
        "bacillus")
            echo "Bacilli|Bacillales|Bacillaceae|Bacillus"
            ;;
        *)
            # For unknown genera, just return the genus name
            echo "$genus"
            ;;
    esac
}

# ========================== PART 1: GENUS-LEVEL GENOTYPING ==================

echo ""
echo "===================== PART 1: GENUS-LEVEL GENOTYPING ====================="
echo "Using GFF files for taxonomic classification at genus level"

TAXONKIT_DIR="$(pwd)/databases/taxonkit"
export TAXONKIT_DB="${TAXONKIT_DIR}"

TAXON_KIT_OUTDIR="${OUTPUT_DIR}/genotyping/genus_level"

if [ ! -d "$GFF_DIR" ] || [ -z "$(ls -A ${GFF_DIR}/*.gff 2>/dev/null)" ]; then
    echo "[WARNING] GFF directory not found or empty: ${GFF_DIR}"
    echo "Skipping genus-level genotyping with GFF files"
else
    for GFF_FILE in ${GFF_DIR}/*.gff; do
        ACCESSION=$(basename "${GFF_FILE}" .gff)
        
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing genus-level classification for ${ACCESSION}"
        
        # Record start of analysis in log
        echo -e "\n\n======== TAXONKIT ANALYSIS: $ACCESSION ========" >> "${LOG_DIR}/taxonkit.log"
        
        run_with_metrics "grep -o 'em_target=[0-9]\+' ${GFF_FILE} | \
            sed 's/em_target=//' | \
            taxonkit lineage | \
            taxonkit reformat -f '{g}' | \
            sort | \
            uniq -c > ${TAXON_KIT_OUTDIR}/${ACCESSION}_genus_classification.txt" \
            "TaxonKit analysis for ${ACCESSION}" "${LOG_DIR}/taxonkit.log"
        
        # Determine the most abundant genus
        MOST_ABUNDANT_GENUS=$(sort -nr "${TAXON_KIT_OUTDIR}/${ACCESSION}_genus_classification.txt" | head -n 1 | awk '{print $2}')
        echo "[INFO] Most abundant genus for ${ACCESSION}: ${MOST_ABUNDANT_GENUS:-Unknown}"
        echo "${ACCESSION}${MOST_ABUNDANT_GENUS:+\t$MOST_ABUNDANT_GENUS}" >> "${TAXON_KIT_OUTDIR}/summary_genus_classification.txt"
    done
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed genus-level genotyping"
fi

# ========================== PART 2: SPECIES-LEVEL GENOTYPING ================

echo ""
echo "=================== PART 2: SPECIES-LEVEL GENOTYPING ==================="
echo "Using MLST and FastANI for strain typing and species identification"

# Set directories
GENOTYPING_OUTDIR="${OUTPUT_DIR}/genotyping/species_level"
MLST_DIR="${GENOTYPING_OUTDIR}/mlst"
FASTANI_DIR="${GENOTYPING_OUTDIR}/fastani"
REFERENCE_GENOMES_DIR="${GENOTYPING_OUTDIR}/reference_genomes"

mkdir -p "${MLST_DIR}" "${FASTANI_DIR}" "${REFERENCE_GENOMES_DIR}"

# ------------------- Part 2a: MLST Analysis -----------------------------

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting MLST analysis"

# Run MLST on all assemblies
run_with_metrics "mlst ${ASSEMBLY_DIR}/*.fasta" "MLST analysis" "${LOG_DIR}/mlst.log" > "${MLST_DIR}/mlst_raw.tsv"

# Process MLST results
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing MLST results"
# Check if file exists and has content
if [ -s "${MLST_DIR}/mlst_raw.tsv" ]; then
  # Process only lines that look like MLST results (contain .fasta)
  grep "\.fasta" "${MLST_DIR}/mlst_raw.tsv" | awk '{
    # Get the basename without path or extension
    accession = $1;
    gsub(".*/", "", accession);  # Remove path
    gsub("\\.fasta$", "", accession);  # Remove extension
    
    # Store the original line
    original = $0;
    
    # Replace the first column with just the accession ID
    sub($1, accession, original);
    
    # Print the modified line
    print original;
  }' > "${MLST_DIR}/mlst_processed.tsv"
else
  echo "Error: MLST raw results file is empty or missing"
  touch "${MLST_DIR}/mlst_processed.tsv"  # Create empty file to prevent further errors
fi

# Create the header and final MLST file
echo -e "file\tscheme\tST\tallele1\tallele2\tallele3\tallele4\tallele5\tallele6\tallele7" > "${MLST_DIR}/mlst_results.tsv"
cat "${MLST_DIR}/mlst_processed.tsv" >> "${MLST_DIR}/mlst_results.tsv"

# Extract genus information from MLST results for downstream analysis
DETECTED_GENUS=$(head -n 2 "${MLST_DIR}/mlst_results.tsv" 2>/dev/null | tail -n 1 | cut -f2 2>/dev/null || echo "unknown")
echo "[INFO] Primary detected genus from MLST: ${DETECTED_GENUS}"

# ------------------- Part 2b: FastANI Analysis --------------------------

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting FastANI analysis"

# Function to download reference genomes based on detected genus
download_reference_genomes() {
    local genus="$1"
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Downloading reference genomes for ${genus}"
    cd "${REFERENCE_GENOMES_DIR}" || return 1
    
    case "$(echo "${genus}" | tr '[:upper:]' '[:lower:]')" in
        "campylobacter")
            # Download C. jejuni reference (NCTC11168)
            if [ ! -f "c_jejuni_reference.fna" ]; then
                echo "Downloading Campylobacter jejuni reference genome..."
                wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/085/GCF_000009085.1_ASM908v1/GCF_000009085.1_ASM908v1_genomic.fna.gz
                gunzip GCF_000009085.1_ASM908v1_genomic.fna.gz
                mv GCF_000009085.1_ASM908v1_genomic.fna c_jejuni_reference.fna
            fi
            
            # Download C. coli reference (RM2228)
            if [ ! -f "c_coli_reference.fna" ]; then
                echo "Downloading Campylobacter coli reference genome..."
                wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/018/045/GCF_000018045.1_ASM1804v1/GCF_000018045.1_ASM1804v1_genomic.fna.gz
                gunzip GCF_000018045.1_ASM1804v1_genomic.fna.gz
                mv GCF_000018045.1_ASM1804v1_genomic.fna c_coli_reference.fna
            fi
            ;;
            
        "salmonella")
            # Download S. enterica reference
            if [ ! -f "s_enterica_reference.fna" ]; then
                echo "Downloading Salmonella enterica reference genome..."
                wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.fna.gz
                gunzip GCF_000006945.2_ASM694v2_genomic.fna.gz
                mv GCF_000006945.2_ASM694v2_genomic.fna s_enterica_reference.fna
            fi
            ;;
            
        "escherichia")
            # Download E. coli reference
            if [ ! -f "e_coli_reference.fna" ]; then
                echo "Downloading Escherichia coli reference genome..."
                wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
                gunzip GCF_000005845.2_ASM584v2_genomic.fna.gz
                mv GCF_000005845.2_ASM584v2_genomic.fna e_coli_reference.fna
            fi
            ;;
            
        *)
            echo "[WARNING] No pre-configured reference genomes for genus: ${genus}"
            echo "Please manually add reference genomes to: ${REFERENCE_GENOMES_DIR}"
            return 1
            ;;
    esac
    
    cd - > /dev/null || return 1
    return 0
}

# Download reference genomes based on detected genus
download_reference_genomes "${DETECTED_GENUS}"

# Check if reference genomes were downloaded or exist
if [ -z "$(ls -A ${REFERENCE_GENOMES_DIR}/*.fna 2>/dev/null)" ]; then
    echo "[WARNING] No reference genomes found for FastANI analysis"
    echo "Consider manually adding reference genomes to: ${REFERENCE_GENOMES_DIR}"
else
    # Create a file listing all reference genomes
    find "${REFERENCE_GENOMES_DIR}" -name "*.fna" > "${REFERENCE_GENOMES_DIR}/reference_list.txt"
    echo "[INFO] Found $(wc -l < "${REFERENCE_GENOMES_DIR}/reference_list.txt") reference genome(s) for FastANI analysis"
    
    # Run FastANI for each assembly against all references
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running FastANI analysis"
    
    # Create header for species predictions file
    echo -e "Accession\tBest_Match\tANI_Value\tSpecies_Prediction" > "${FASTANI_DIR}/species_predictions.tsv"
    
    for ASSEMBLY in ${ASSEMBLY_DIR}/*.fasta; do
        ACCESSION=$(basename "$ASSEMBLY" .fasta)
        
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing ${ACCESSION} with FastANI"
        echo -e "\n\n======== FASTANI ANALYSIS: $ACCESSION ========" >> "${LOG_DIR}/fastani.log"
        
        # Run FastANI comparing query genome to all references
        run_with_metrics "fastANI -q $ASSEMBLY --rl ${REFERENCE_GENOMES_DIR}/reference_list.txt -o ${FASTANI_DIR}/${ACCESSION}_fastani.out" \
                        "FastANI for ${ACCESSION}" "${LOG_DIR}/fastani.log"
        
        # Extract the best hit (highest ANI value)
        if [[ -s "${FASTANI_DIR}/${ACCESSION}_fastani.out" ]]; then
            BEST_HIT=$(sort -k3,3nr "${FASTANI_DIR}/${ACCESSION}_fastani.out" | head -n 1)
            
            # Parse the best hit
            ANI_VALUE=$(echo "$BEST_HIT" | awk '{print $3}')
            REF_GENOME=$(echo "$BEST_HIT" | awk '{print $2}')
            REF_SPECIES=$(basename "$REF_GENOME" | sed 's/_reference.fna//')
            
            # Determine species based on ANI threshold
            if (( $(echo "$ANI_VALUE >= $MIN_ANI_THRESHOLD" | bc -l) )); then
                SPECIES_PREDICTION="$REF_SPECIES"
                echo "[INFO] ${ACCESSION}: Identified as ${SPECIES_PREDICTION} (ANI: ${ANI_VALUE}%)"
            else
                SPECIES_PREDICTION="Unknown_${DETECTED_GENUS}_sp"
                echo "[INFO] ${ACCESSION}: Unknown ${DETECTED_GENUS} species (Best match: ${REF_SPECIES}, ANI: ${ANI_VALUE}%)"
            fi
            
            # Add to species predictions file
            echo -e "${ACCESSION}\t${REF_SPECIES}\t${ANI_VALUE}\t${SPECIES_PREDICTION}" >> "${FASTANI_DIR}/species_predictions.tsv"
        else
            SPECIES_PREDICTION="No_match"
            ANI_VALUE="NA"
            echo "[WARNING] ${ACCESSION}: No FastANI match found"
            echo -e "${ACCESSION}\tNo_match\tNA\tNo_match" >> "${FASTANI_DIR}/species_predictions.tsv"
        fi
    done
fi

# ------------------- Part 2c: Combine MLST and FastANI Results ----------

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Combining MLST and FastANI results"

# Create header for combined results
echo -e "Accession\tMLST_Genus\tMLST_ST\tFastANI_Species\tANI_Value\tAllele1\tAllele2\tAllele3\tAllele4\tAllele5\tAllele6\tAllele7" > "${GENOTYPING_OUTDIR}/combined_results.tsv"

# Process each genome
for ASSEMBLY in ${ASSEMBLY_DIR}/*.fasta; do
    ACCESSION=$(basename "$ASSEMBLY" .fasta)
    
    # Get MLST information
    MLST_INFO=$(grep "^$ACCESSION" "${MLST_DIR}/mlst_results.tsv" || echo -e "$ACCESSION\tUnknown\tUnknown\t-\t-\t-\t-\t-\t-\t-")
    MLST_GENUS=$(echo "$MLST_INFO" | cut -f2)
    MLST_ST=$(echo "$MLST_INFO" | cut -f3)
    ALLELE1=$(echo "$MLST_INFO" | cut -f4)
    ALLELE2=$(echo "$MLST_INFO" | cut -f5)
    ALLELE3=$(echo "$MLST_INFO" | cut -f6)
    ALLELE4=$(echo "$MLST_INFO" | cut -f7)
    ALLELE5=$(echo "$MLST_INFO" | cut -f8)
    ALLELE6=$(echo "$MLST_INFO" | cut -f9)
    ALLELE7=$(echo "$MLST_INFO" | cut -f10)
    
    # Get FastANI information if available
    FASTANI_INFO=$(grep "^$ACCESSION" "${FASTANI_DIR}/species_predictions.tsv" 2>/dev/null)
    if [[ -n "$FASTANI_INFO" ]]; then
        SPECIES_PREDICTION=$(echo "$FASTANI_INFO" | cut -f4)
        ANI_VALUE=$(echo "$FASTANI_INFO" | cut -f3)
    else
        SPECIES_PREDICTION="No_data"
        ANI_VALUE="NA"
    fi
    
    # Write to combined results file
    echo -e "$ACCESSION\t$MLST_GENUS\t$MLST_ST\t$SPECIES_PREDICTION\t$ANI_VALUE\t$ALLELE1\t$ALLELE2\t$ALLELE3\t$ALLELE4\t$ALLELE5\t$ALLELE6\t$ALLELE7" >> "${GENOTYPING_OUTDIR}/combined_results.tsv"
done

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Genotyping analysis complete"

# ========================== PART 3: QUALITY ASSESSMENT =======================

echo ""
echo "===================== PART 3: QUALITY ASSESSMENT ======================"
echo "Using CheckM to assess genome completeness and contamination"

QA_OUTDIR="${OUTPUT_DIR}/quality_assessment"
CHECKM_TEMP="${QA_OUTDIR}/checkm_temp"
TEMP_QUALITY_DIR="${CHECKM_TEMP}/quality_files"

mkdir -p "${QA_OUTDIR}" "${CHECKM_TEMP}" "${TEMP_QUALITY_DIR}"

# Set CheckM database location
checkm data setRoot "$(realpath databases/checkm)"

# Read the combined results from previous step
COMBINED_RESULTS="${GENOTYPING_OUTDIR}/combined_results.tsv"

if [ ! -f "$COMBINED_RESULTS" ]; then
    echo "[WARNING] Combined results file not found: ${COMBINED_RESULTS}"
    echo "Using assembly files directly for CheckM analysis"
    
    # Process each assembly directly
    for ASSEMBLY in ${ASSEMBLY_DIR}/*.fasta; do
        ACCESSION=$(basename "$ASSEMBLY" .fasta)
        
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running CheckM analysis for ${ACCESSION}"
        
        # Create temporary directory for this assembly
        ACCESSION_DIR="${CHECKM_TEMP}/${ACCESSION}"
        mkdir -p "${ACCESSION_DIR}"/{bin,output}
        
        # Create a symbolic link to the assembly
        ln -sf "$(realpath "$ASSEMBLY")" "${ACCESSION_DIR}/bin/${ACCESSION}.fasta"
        
        # Run CheckM with the general bacterial marker set
        run_with_metrics "checkm lineage_wf --threads ${THREADS} --extension fasta --file ${TEMP_QUALITY_DIR}/${ACCESSION}_quality.tsv --tab_table ${ACCESSION_DIR}/bin ${ACCESSION_DIR}/output" \
                        "CheckM for ${ACCESSION}" "${LOG_DIR}/checkm.log"
        
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed CheckM analysis for ${ACCESSION}"
    done
else
    # Iterate through each line of the combined results file, skipping the header
    tail -n +2 "$COMBINED_RESULTS" | while read -r line; do
        # Extract the accession and taxonomic information
        ACCESSION=$(echo "$line" | cut -f1)
        GENUS=$(echo "$line" | cut -f2)
        SPECIES=$(echo "$line" | cut -f4 | cut -d "_" -f 2 2>/dev/null)
        
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running CheckM analysis for ${ACCESSION}"
        
        # Create temporary directory for this assembly
        ACCESSION_DIR="${CHECKM_TEMP}/${ACCESSION}"
        mkdir -p "${ACCESSION_DIR}"/{bin,output}
        
        # Create a symbolic link to the assembly
        ln -sf "$(realpath "${ASSEMBLY_DIR}/${ACCESSION}.fasta")" "${ACCESSION_DIR}/bin/${ACCESSION}.fasta"
        
        # Try to find appropriate CheckM taxonomic level
        echo -e "\n\n======== CHECKM ANALYSIS: $ACCESSION ========" >> "${LOG_DIR}/checkm.log"
        
        # Convert genus and species to lowercase for case-insensitive matching
        genus_species_lower=$(echo "${GENUS} ${SPECIES}" | tr '[:upper:]' '[:lower:]')
        
        # Check if lineage workflow should be used (safer, use this by default)
        use_lineage_workflow=true
        
        # Only attempt taxon-specific analysis if taxon_list command works
        if checkm taxon_list > /dev/null 2>&1; then
            # Check if there's a specific marker set for this species
            checkm_taxon=""
            checkm_species=$(checkm taxon_list 2>/dev/null | grep -i "$genus_species_lower" | grep -i "species" | head -n 1 || echo "")
            
            if [ -n "$checkm_species" ]; then
                # Extract the taxon name for the species
                checkm_taxon=$(echo "$checkm_species" | awk '{print $2" "$3}')
                checkm_level="species"
                echo "[INFO] Using species-specific marker set: $checkm_taxon"
                use_lineage_workflow=false
            else
                # Try to find a genus-level marker set
                genus_lower=$(echo "${GENUS}" | tr '[:upper:]' '[:lower:]')
                checkm_genus=$(checkm taxon_list 2>/dev/null | grep -i "$genus_lower" | grep -i "genus" | head -n 1 || echo "")
                
                if [ -n "$checkm_genus" ]; then
                    # Extract the taxon name for the genus
                    checkm_taxon=$(echo "$checkm_genus" | awk '{print $2}')
                    checkm_level="genus"
                    echo "[INFO] Using genus-specific marker set: $checkm_taxon"
                    use_lineage_workflow=false
                else
                    echo "[INFO] No specific marker set found. Using domain bacteria."
                fi
            fi
        else
            echo "[INFO] CheckM taxon_list command failed. Using lineage workflow."
        fi
        
        # Run CheckM with the appropriate marker set
        if [ "$use_lineage_workflow" = true ]; then
            # Use lineage workflow for domain level or if taxon-specific failed
            run_with_metrics "checkm lineage_wf --threads ${THREADS} --extension fasta --file ${TEMP_QUALITY_DIR}/${ACCESSION}_quality.tsv --tab_table ${ACCESSION_DIR}/bin ${ACCESSION_DIR}/output" \
                            "CheckM domain-level for ${ACCESSION}" "${LOG_DIR}/checkm.log"
        else
            # Attempt to create marker file for the specific taxon
            if ! checkm taxon_set "$checkm_level" "$checkm_taxon" "${ACCESSION_DIR}/${ACCESSION}.markers" 2>/dev/null; then
                echo "[WARNING] Taxon-specific approach failed. Falling back to lineage workflow."
                run_with_metrics "checkm lineage_wf --threads ${THREADS} --extension fasta --file ${TEMP_QUALITY_DIR}/${ACCESSION}_quality.tsv --tab_table ${ACCESSION_DIR}/bin ${ACCESSION_DIR}/output" \
                                "CheckM domain-level for ${ACCESSION}" "${LOG_DIR}/checkm.log"
            else
                # Run CheckM analyze and qa
                run_with_metrics "checkm analyze --threads ${THREADS} --extension fasta ${ACCESSION_DIR}/${ACCESSION}.markers ${ACCESSION_DIR}/bin ${ACCESSION_DIR}/output" \
                                "CheckM analyze for ${ACCESSION}" "${LOG_DIR}/checkm.log"
                
                run_with_metrics "checkm qa --file ${TEMP_QUALITY_DIR}/${ACCESSION}_quality.tsv --out_format 2 --threads ${THREADS} ${ACCESSION_DIR}/${ACCESSION}.markers ${ACCESSION_DIR}/output" \
                                "CheckM qa for ${ACCESSION}" "${LOG_DIR}/checkm.log"
            fi
        fi
        
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed CheckM analysis for ${ACCESSION}"
    done
fi

# Process all quality files to create a properly formatted final TSV
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Combining CheckM results"

# Check if any quality files were generated
if [ -z "$(ls -A ${TEMP_QUALITY_DIR}/*_quality.tsv 2>/dev/null)" ]; then
    echo "[WARNING] No CheckM quality files were generated"
else
    # Create a temporary file for all data rows
    DATA_FILE="${CHECKM_TEMP}/data_rows.tsv"
    > "${DATA_FILE}"  # Empty the file
    
    # Extract data from each quality file, removing leading spaces and fixing spacing
    for file in ${TEMP_QUALITY_DIR}/*_quality.tsv; do
        # Extract accession ID from filename
        accession=$(basename "$file" _quality.tsv)
        
        # Process the file to extract data rows and fix formatting
        CONTENT=$(cat "$file" | grep -v "^----" | awk '{gsub(/[ ]{2,}/,"\t"); print}')
        
        # Get the data row and add accession ID as the first column
        DATA_ROW=$(echo "$CONTENT" | tail -n 1 | sed 's/^[[:space:]]*//' | sed 's/[[:space:]][[:space:]]\+/\t/g')
        echo -e "${accession}\t${DATA_ROW}" >> "${DATA_FILE}"
        
        # Extract header if not already done
        if [ ! -f "${CHECKM_TEMP}/header.tsv" ]; then
            HEADER=$(echo "$CONTENT" | head -n 1 | sed 's/^[[:space:]]*//' | sed 's/[[:space:]][[:space:]]\+/\t/g')
            echo -e "Accession\t${HEADER}" > "${CHECKM_TEMP}/header.tsv"
        fi
    done
    
    # Combine header and data rows into final output
    cat "${CHECKM_TEMP}/header.tsv" "${DATA_FILE}" > "${QA_OUTDIR}/checkm_quality.tsv"
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] CheckM quality assessment complete"
fi

# ========================== PART 4: TAXONOMIC CLASSIFICATION =================

echo ""
echo "================= PART 4: TAXONOMIC CLASSIFICATION ================="
echo "Contig-by-contig taxonomy assessment with Kraken2"

# Skip this section if Kraken database is not provided
if [ -z "${KRAKEN_DB_DIR}" ]; then
    echo "[WARNING] Kraken database directory not provided, skipping taxonomy assessment"
else
    # Check if Kraken database exists
    if [ ! -d "${KRAKEN_DB_DIR}" ]; then
        echo "[WARNING] Kraken database directory not found: ${KRAKEN_DB_DIR}"
        echo "Skipping taxonomy assessment"
    else
        TAXONOMY_OUTDIR="${OUTPUT_DIR}/taxonomy"
        TEMP_KRAKEN="${TAXONOMY_OUTDIR}/temp"
        SUMMARY_DIR="${TAXONOMY_OUTDIR}/summaries"
        
        mkdir -p "${TAXONOMY_OUTDIR}" "${TEMP_KRAKEN}" "${SUMMARY_DIR}"
        
        for ASSEMBLY in ${ASSEMBLY_DIR}/*.fasta; do
            ACCESSION=$(basename "$ASSEMBLY" .fasta)
            
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running Kraken2 analysis for ${ACCESSION}"
            echo -e "\n\n======== KRAKEN ANALYSIS: $ACCESSION ========" >> "${LOG_DIR}/kraken.log"
            
            # Create output directories for this accession
            mkdir -p "${TEMP_KRAKEN}/${ACCESSION}/"{contigs,reports,kraken}
            
            # Set path variables
            CONTIGS_DIR="${TEMP_KRAKEN}/${ACCESSION}/contigs"
            REPORTS_DIR="${TEMP_KRAKEN}/${ACCESSION}/reports"
            KRAKEN_OUT_DIR="${TEMP_KRAKEN}/${ACCESSION}/kraken"
            
            echo "Step 1: Splitting assembly into contigs..."
            # Split the assembly into individual contigs

            cat > /tmp/split_contigs.awk << 'EOF'
                BEGIN {count=0; seq=""}
                /^>/ {
                if (count > 0 && seq != "") {
                    print seq > (outdir "/contig_" count ".fasta");
                }
                count++;
                header=$0;
                print header > (outdir "/contig_" count ".fasta");
                seq="";
                next;
                }
                {
                seq = seq $0 "\n";
                }
                END {
                if (seq != "") {
                    print seq > (outdir "/contig_" count ".fasta");
                }
                }
EOF

            # Then run it
            run_with_metrics "awk -v outdir='${CONTIGS_DIR}' -f /tmp/split_contigs.awk '${ASSEMBLY}'" \
            "Splitting contigs for ${ACCESSION}" "${LOG_DIR}/kraken.log"

            # Clean up
            rm /tmp/split_contigs.awk
            
            # Verify files were created
            echo "Verifying contig files in ${CONTIGS_DIR}..."
            ls -la "${CONTIGS_DIR}" >> "${LOG_DIR}/kraken.log"

            echo "Step 2: Running Kraken2 on each contig..."
            # Create output file for kraken results
            KRAKEN_OUT="${KRAKEN_OUT_DIR}/${ACCESSION}_kraken.out"
            KRAKEN_REPORT="${REPORTS_DIR}/${ACCESSION}_kraken.report"
            
            # Run Kraken2 on the entire assembly first (for overall classification)
            run_with_metrics "kraken2 --db ${KRAKEN_DB_DIR} --threads ${THREADS} \
                              --quick \
                              --memory-mapping \
                              --output ${KRAKEN_OUT} \
                              --report ${KRAKEN_REPORT} \
                              ${ASSEMBLY}" \
                             "Kraken2 on ${ACCESSION} (full assembly)" "${LOG_DIR}/kraken.log"
            
            # Run Kraken2 on each contig
            for CONTIG in ${CONTIGS_DIR}/*.fasta; do
                CONTIG_NAME=$(basename "${CONTIG}" .fasta)
                
                CONTIG_OUT="${KRAKEN_OUT_DIR}/${CONTIG_NAME}_kraken.out"
                CONTIG_REPORT="${REPORTS_DIR}/${CONTIG_NAME}_kraken.report"
                
                run_with_metrics "kraken2 --db ${KRAKEN_DB_DIR} --threads 1 \
                                 --quick \
                                 --memory-mapping \
                                 --output ${CONTIG_OUT} \
                                 --report ${CONTIG_REPORT} \
                                 ${CONTIG}" \
                                "Kraken2 on ${CONTIG_NAME}" "${LOG_DIR}/kraken.log"
            done
            
            echo "Step 3: Summarizing Kraken2 results..."
            # Create a summary of contig-by-contig taxonomy
            SUMMARY_FILE="${SUMMARY_DIR}/${ACCESSION}_taxonomy_summary.tsv"
            
            # Create header for summary file
            echo -e "Contig\tLength\tTaxID\tRank\tTaxName\tPercentage" > "${SUMMARY_FILE}"
            
            # Process each contig report
            for REPORT in ${REPORTS_DIR}/contig_*_kraken.report; do
                CONTIG_NAME=$(basename "${REPORT}" _kraken.report)
                
                # Extract the length of the contig
                CONTIG_LENGTH=$(grep -v "^>" "${CONTIGS_DIR}/${CONTIG_NAME}.fasta" | tr -d '\n' | wc -c)
                
                # Get the top hit from the Kraken report
                TOP_HIT=$(grep -v "unclassified" "${REPORT}" | sort -k1,1nr | head -n 1 || echo "0\t0\t0\tunclassified\tNo hits")
                
                # Extract relevant information
                PERCENTAGE=$(echo "${TOP_HIT}" | cut -f1)
                TAXID=$(echo "${TOP_HIT}" | cut -f5)
                RANK=$(echo "${TOP_HIT}" | cut -f4)
                TAXNAME=$(echo "${TOP_HIT}" | cut -f6 | sed 's/^ *//')
                
                # Add to summary file
                echo -e "${CONTIG_NAME}\t${CONTIG_LENGTH}\t${TAXID}\t${RANK}\t${TAXNAME}\t${PERCENTAGE}" >> "${SUMMARY_FILE}"
            done
            
            # Create an assembly-level summary
            ASSEMBLY_SUMMARY="${SUMMARY_DIR}/${ACCESSION}_assembly_summary.tsv"
            
            # Extract top hits at each taxonomic level
            echo -e "Level\tTaxID\tName\tPercentage" > "${ASSEMBLY_SUMMARY}"
            for LEVEL in "D" "P" "C" "O" "F" "G" "S"; do
                TOP_HIT=$(grep "${LEVEL}" "${KRAKEN_REPORT}" | sort -k1,1nr | head -n 1 || echo "0\t0\t${LEVEL}\tNo hit\t0")
                PERCENTAGE=$(echo "${TOP_HIT}" | cut -f1)
                TAXID=$(echo "${TOP_HIT}" | cut -f5)
                NAME=$(echo "${TOP_HIT}" | cut -f6 | sed 's/^ *//')
                
                # Map level code to full name
                case "${LEVEL}" in
                    "D") LEVEL_NAME="Domain" ;;
                    "P") LEVEL_NAME="Phylum" ;;
                    "C") LEVEL_NAME="Class" ;;
                    "O") LEVEL_NAME="Order" ;;
                    "F") LEVEL_NAME="Family" ;;
                    "G") LEVEL_NAME="Genus" ;;
                    "S") LEVEL_NAME="Species" ;;
                    *) LEVEL_NAME="Unknown" ;;
                esac
                
                echo -e "${LEVEL_NAME}\t${TAXID}\t${NAME}\t${PERCENTAGE}" >> "${ASSEMBLY_SUMMARY}"
            done
            
            # Calculate contamination metrics
            echo "Step 4: Assessing potential contamination..."
            CONTAM_REPORT="${SUMMARY_DIR}/${ACCESSION}_contamination_report.txt"
            
            # Get expected genus from previous analysis
            EXPECTED_GENUS=$(grep "^${ACCESSION}" "${GENOTYPING_OUTDIR}/combined_results.tsv" 2>/dev/null | cut -f2 || echo "Unknown")
            
            # Start the contamination report
            echo "Contamination Assessment for ${ACCESSION}" > "${CONTAM_REPORT}"
            echo "Expected genus: ${EXPECTED_GENUS}" >> "${CONTAM_REPORT}"
            echo "-------------------------------------------" >> "${CONTAM_REPORT}"
            
            # Calculate percent of contigs matching expected genus
            TOTAL_CONTIGS=$(wc -l < "${SUMMARY_FILE}")
            # Adjust for header line
            TOTAL_CONTIGS=$((TOTAL_CONTIGS - 1))
            
            # Get taxonomic pattern to search for based on the expected genus
            TAXONOMY_PATTERN=$(get_taxonomy_hierarchy "${EXPECTED_GENUS}")

            # Count matching genus contigs
            MATCHING_GENUS=$(grep -i -E "${TAXONOMY_PATTERN}" "${SUMMARY_FILE}" | wc -l || echo "0")
            
            # Calculate percentage
            if [ ${TOTAL_CONTIGS} -gt 0 ]; then
                MATCHING_PERCENT=$(echo "scale=2; ${MATCHING_GENUS} * 100 / ${TOTAL_CONTIGS}" | bc)
            else
                MATCHING_PERCENT="0.00"
            fi
            
            echo "Total contigs analyzed: ${TOTAL_CONTIGS}" >> "${CONTAM_REPORT}"
            echo "Contigs matching expected genus: ${MATCHING_GENUS}" >> "${CONTAM_REPORT}"
            echo "Percent matching expected genus: ${MATCHING_PERCENT}%" >> "${CONTAM_REPORT}"
            echo "-------------------------------------------" >> "${CONTAM_REPORT}"
            
            # Identify potential contaminants (contigs with different genus)
            echo "Potential contamination detected in the following contigs:" >> "${CONTAM_REPORT}"
            grep -v "^Contig" "${SUMMARY_FILE}" | grep -v -i "${EXPECTED_GENUS}" >> "${CONTAM_REPORT}" || echo "None detected" >> "${CONTAM_REPORT}"

            # First remove the header line, then find non-matching lines
            if [ ${MATCHING_GENUS} -lt ${TOTAL_CONTIGS} ]; then
                sed '1d' "${SUMMARY_FILE}" | grep -v -i -E "${TAXONOMY_PATTERN}" >> "${CONTAM_REPORT}" || true
            else
                echo "None detected" >> "${CONTAM_REPORT}"
            fi
            
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed taxonomy assessment for ${ACCESSION}"
        done
        
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Taxonomy assessment complete"
    fi
fi

# ========================== PART 5: FINAL SUMMARY ===========================

echo ""
echo "====================== PART 5: GENERATING SUMMARY ====================="
echo "Creating final summary reports combining all analysis results"

# Create final summary directory
FINAL_DIR="${OUTPUT_DIR}/final_summary"
mkdir -p "${FINAL_DIR}"

# Create a master summary file
MASTER_SUMMARY="${FINAL_DIR}/master_summary.tsv"

# Create header
echo -e "Accession\tGenus\tSpecies\tMLST_ST\tFastANI_Value\tCompleteness\tContamination\tN50\tGenome_Size\tGC_Content\tContigs\tPotential_Contamination" > "${MASTER_SUMMARY}"

# Process each genome and combine all results
for ASSEMBLY in ${ASSEMBLY_DIR}/*.fasta; do
    ACCESSION=$(basename "$ASSEMBLY" .fasta)
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Generating summary for ${ACCESSION}"
    
    # Get genotyping data
    GENOTYPING_DATA=$(grep "^${ACCESSION}" "${GENOTYPING_OUTDIR}/combined_results.tsv" 2>/dev/null || echo "${ACCESSION}\tUnknown\tUnknown\tUnknown\tNA")
    GENUS=$(echo "${GENOTYPING_DATA}" | cut -f2)
    SPECIES=$(echo "${GENOTYPING_DATA}" | cut -f4 | cut -d "_" -f 2 2>/dev/null || echo "Unknown")
    MLST_ST=$(echo "${GENOTYPING_DATA}" | cut -f3)
    ANI_VALUE=$(echo "${GENOTYPING_DATA}" | cut -f5)
    
    # Get quality assessment data
    QA_DATA=$(grep "^${ACCESSION}" "${QA_OUTDIR}/checkm_quality.tsv" 2>/dev/null || echo "${ACCESSION}\tNA\tNA")
    COMPLETENESS=$(echo "${QA_DATA}" | cut -f12 2>/dev/null || echo "NA")
    CONTAMINATION=$(echo "${QA_DATA}" | cut -f13 2>/dev/null || echo "NA")
    
    # Calculate genome metrics
    # Count contigs
    NUM_CONTIGS=$(grep -c "^>" "${ASSEMBLY}" || echo "NA")
    
    # Calculate genome size
    GENOME_SIZE=$(grep -v "^>" "${ASSEMBLY}" | tr -d '\n' | wc -c || echo "NA")
    
    # Calculate GC content
    GC_COUNT=$(grep -v "^>" "${ASSEMBLY}" | tr -d '\n' | tr -cd 'GCgc' | wc -c || echo "0")
    TOTAL_BASES=${GENOME_SIZE}
    if [ "${TOTAL_BASES}" != "NA" ] && [ ${TOTAL_BASES} -gt 0 ]; then
        GC_CONTENT=$(echo "scale=2; ${GC_COUNT} * 100 / ${TOTAL_BASES}" | bc)
    else
        GC_CONTENT="NA"
    fi
    
    # Calculate N50
    TMP_LENGTHS="${FINAL_DIR}/${ACCESSION}_contig_lengths.tmp"
    grep -A 1 "^>" "${ASSEMBLY}" | grep -v "^>" | grep -v "^--" | awk '{print length($0)}' | sort -nr > "${TMP_LENGTHS}"
    
    TOTAL_LENGTH=$(awk '{sum+=$1} END {print sum}' "${TMP_LENGTHS}")
    HALF_LENGTH=$((TOTAL_LENGTH / 2))
    
    current_sum=0
    N50="NA"
    
    while read -r length; do
        current_sum=$((current_sum + length))
        if [ ${current_sum} -ge ${HALF_LENGTH} ]; then
            N50=${length}
            break
        fi
    done < "${TMP_LENGTHS}"
    
    rm "${TMP_LENGTHS}"
    
    # Check for potential contamination
    CONTAM_REPORT="${OUTPUT_DIR}/taxonomy/summaries/${ACCESSION}_contamination_report.txt"
    if [ -f "${CONTAM_REPORT}" ]; then
        POTENTIAL_CONTAM=$(grep "Percent matching expected genus" "${CONTAM_REPORT}" | cut -d ":" -f2 | sed 's/%//' | sed 's/ //g')
        if [ -n "${POTENTIAL_CONTAM}" ]; then
            # If less than 50% match, flag as potential contamination
            if (( $(echo "${POTENTIAL_CONTAM} < 50" | bc -l) )); then
                CONTAM_FLAG="Yes (${POTENTIAL_CONTAM}%)"
            else
                CONTAM_FLAG="No"
            fi
        else
            CONTAM_FLAG="Unknown"
        fi
    else
        CONTAM_FLAG="Not analyzed"
    fi
    
    # Add to master summary
    echo -e "${ACCESSION}\t${GENUS}\t${SPECIES}\t${MLST_ST}\t${ANI_VALUE}\t${COMPLETENESS}\t${CONTAMINATION}\t${N50}\t${GENOME_SIZE}\t${GC_CONTENT}\t${NUM_CONTIGS}\t${CONTAM_FLAG}" >> "${MASTER_SUMMARY}"
done

# Create a quality category field
awk 'BEGIN {FS="\t"; OFS="\t"} 
     NR==1 {print $0, "Quality_Category"}
     NR>1 {
         # Default to "Unknown"
         quality="Unknown";
         
         # Only assign if we have completeness and contamination values
         if ($6 != "NA" && $7 != "NA") {
             if ($6 >= 95 && $7 <= 5) quality="High";
             else if ($6 >= 90 && $7 <= 10) quality="Medium";
             else quality="Low";
         }
         
         print $0, quality;
     }' "${MASTER_SUMMARY}" > "${FINAL_DIR}/master_summary_with_quality.tsv"

# Generate a simple HTML report
HTML_REPORT="${FINAL_DIR}/report.html"

cat > "${HTML_REPORT}" << EOL
<!DOCTYPE html>
<html>
<head>
    <title>Bacterial Genome Analysis Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        h1, h2 { color: #333366; }
        table { border-collapse: collapse; width: 100%; margin-bottom: 20px; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
        tr:nth-child(even) { background-color: #f9f9f9; }
        .high-quality { background-color: #d6f5d6; }
        .medium-quality { background-color: #fff4cc; }
        .low-quality { background-color: #ffcccc; }
        .summary { margin-bottom: 30px; }
    </style>
</head>
<body>
    <h1>Bacterial Genome Analysis Report</h1>
    <p>Report generated on $(date '+%Y-%m-%d %H:%M:%S')</p>
    
    <div class="summary">
        <h2>Analysis Summary</h2>
        <p>Total genomes analyzed: $(tail -n +2 "${FINAL_DIR}/master_summary_with_quality.tsv" | wc -l)</p>
EOL

# Add genome quality statistics to HTML
HIGH_QUALITY=$(grep "High" "${FINAL_DIR}/master_summary_with_quality.tsv" | wc -l || echo "0")
MEDIUM_QUALITY=$(grep "Medium" "${FINAL_DIR}/master_summary_with_quality.tsv" | wc -l || echo "0")
LOW_QUALITY=$(grep "Low" "${FINAL_DIR}/master_summary_with_quality.tsv" | wc -l || echo "0")

cat >> "${HTML_REPORT}" << EOL
        <ul>
            <li>High-quality genomes: ${HIGH_QUALITY}</li>
            <li>Medium-quality genomes: ${MEDIUM_QUALITY}</li>
            <li>Low-quality genomes: ${LOW_QUALITY}</li>
        </ul>
    </div>
    
    <h2>Genome Analysis Results</h2>
    <table>
        <tr>
            <th>Accession</th>
            <th>Genus</th>
            <th>Species</th>
            <th>MLST ST</th>
            <th>Completeness</th>
            <th>Contamination</th>
            <th>Genome Size</th>
            <th>Quality</th>
        </tr>
EOL

# Add rows for each genome
awk 'BEGIN {FS="\t"} 
     NR>1 {
         quality_class = "";
         if ($13 == "High") quality_class = "high-quality";
         else if ($13 == "Medium") quality_class = "medium-quality";
         else if ($13 == "Low") quality_class = "low-quality";
         
         print "<tr class=\"" quality_class "\">";
         print "    <td>" $1 "</td>";
         print "    <td>" $2 "</td>";
         print "    <td>" $3 "</td>";
         print "    <td>" $4 "</td>";
         print "    <td>" $6 "</td>";
         print "    <td>" $7 "</td>";
         print "    <td>" $9 "</td>";
         print "    <td>" $13 "</td>";
         print "</tr>";
     }' "${FINAL_DIR}/master_summary_with_quality.tsv" >> "${HTML_REPORT}"

cat >> "${HTML_REPORT}" << EOL
    </table>
    
    <h2>Pipeline Information</h2>
    <ul>
        <li>Pipeline version: 2.0</li>
        <li>Execution date: $(date '+%Y-%m-%d')</li>
        <li>Assembly directory: ${ASSEMBLY_DIR}</li>
        <li>Output directory: ${OUTPUT_DIR}</li>
    </ul>
</body>
</html>
EOL

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Generated HTML report: ${HTML_REPORT}"

# Create a simple CSV version for easy import
sed 's/\t/,/g' "${FINAL_DIR}/master_summary_with_quality.tsv" > "${FINAL_DIR}/master_summary.csv"

# ========================== CLEANUP =========================================

echo ""
echo "========================= CLEANUP ========================="
echo "Cleaning up temporary files and finalizing output"

# Ask if user wants to remove temporary files
read -p "Would you like to remove temporary files to save disk space? (y/n) " -n 1 -r
echo    # Move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "Removing temporary files..."
    rm -rf "${CHECKM_TEMP}" "${TEMP_KRAKEN}" "${QA_OUTDIR}/tmp"
    echo "Temporary files removed"
else
    echo "Keeping all temporary files"
fi

# Deactivate conda environment
conda deactivate

# ========================== CONCLUSION ======================================

echo ""
echo "======================= PIPELINE COMPLETED ========================"
echo "Analysis completed successfully!"
echo ""
echo "Summary reports are available in: ${FINAL_DIR}"
echo "  - Master summary (TSV): ${FINAL_DIR}/master_summary_with_quality.tsv"
echo "  - Master summary (CSV): ${FINAL_DIR}/master_summary.csv"
echo "  - HTML report: ${HTML_REPORT}"
echo ""
echo "Pipeline execution log: ${log_file}"
echo ""
echo "Thank you for using the Bacterial Genome Analysis Pipeline v2.0"
echo "======================================================================"

exit 0
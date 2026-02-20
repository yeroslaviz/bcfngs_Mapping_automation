#!/bin/bash
#SBATCH --job-name=fastq_conversion
#SBATCH --output=/fs/pool/pool-bcfngs/scripts/fastq_conversion_%j.log
#SBATCH --error=/fs/pool/pool-bcfngs/scripts/fastq_conversion_%A_%a.err
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --partition=p.hpcl93

# Does it have to be on HPCL93*??

##########################
### User Configuration ###
##########################

####################################################
####### Parameter to change before running ######### 
# Instrument type: "AVITTI" or "ILLUMINA"
INSTRUMENT="ILLUMINA"  # Change this to "ILLUMINA" for Illumina processing
# Project information
PROJECT_NUMBER="P935"  # e.g., P947, P1024, etc.
# Folder name of project
RAW_DATA="P935_250513_A01878_0164_BHNC27DRX5"
####################################################

# Input and Output paths
INPUT_FOLDER="/fs/pool/pool-bcfngs/Raw_data/${RAW_DATA}"  # Replace with your input folder path
OUTPUT_BASE="/fs/pool/pool-bcfngs/fastq_files"    # Replace with your output base path

#SAMPLE_SHEET="/path/to/SampleSheet.csv"  # Only needed for Illumina, this is always inside the input folder of Raw_data

#####  Processing parameters
## Aviti
SPLIT_LANES_ELEMENT=TRUE  # TRUE or FALSE for Element's --split-lanes
# Element Bioscience specific settings
ELEMENT_THREADS=200
ELEMENT_SETTINGS="'I1MismatchThreshold,0,1' 'I2MismatchThreshold,0,1'"
ELEMENT_EXCLUDE_TILE="'L2R..C..S.'"

## Illumina specific settings
NO_LANE_SPLITTING_ILLUMINA=TRUE  # TRUE or FALSE for Illumina's --no-lane-splitting
ILLUMINA_THREADS=65
ILLUMINA_LOADER_THREADS=24
ILLUMINA_WRITER_THREADS=24
ILLUMINA_BARCODE_MISMATCHES=0

###########################
### Script Starts Here  ###
###########################

# Create output directory with project number and timestamp
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
OUTPUT_FOLDER="${OUTPUT_BASE}/${PROJECT_NUMBER}_${TIMESTAMP}"
# mkdir -p "${OUTPUT_FOLDER}"

# Load appropriate modules based on instrument type
if [ "${INSTRUMENT}" = "AVITTI" ]; then
    echo "Loading Element Biosciences bases2fastq environment"
    module purge
    module load bases2fastq  # Adjust version as needed
    
    # Build the bases2fastq command
    CMD="bases2fastq ${INPUT_FOLDER} ${OUTPUT_FOLDER} -p ${ELEMENT_THREADS}"
    
    # Add --split-lanes if TRUE
    if [ "${SPLIT_LANES_ELEMENT}" = "TRUE" ]; then
        CMD="${CMD} --split-lanes"
    fi
    
    # Add settings
    for setting in ${ELEMENT_SETTINGS}; do
        CMD="${CMD} --settings ${setting}"
    done
    
    # Add exclude tile if specified
    if [ -n "${ELEMENT_EXCLUDE_TILE}" ]; then
        CMD="${CMD} --exclude-tile ${ELEMENT_EXCLUDE_TILE}"
    fi
    
elif [ "${INSTRUMENT}" = "ILLUMINA" ]; then
    echo "Loading Illumina bcl2fastq environment"
    module purge
    module load bclconvert  # Adjust version as needed
    
    # Build the bcl2fastq command
#    CMD="bcl2fastq --runfolder-dir ${INPUT_FOLDER} --sample-sheet ${INPUT_FOLDER}/SampleSheet.csv --output-dir ${OUTPUT_FOLDER}"
#    CMD="${CMD} -p ${ILLUMINA_THREADS} -r ${ILLUMINA_LOADER_THREADS} -w ${ILLUMINA_WRITER_THREADS}"
#    CMD="${CMD} --barcode-mismatches ${ILLUMINA_BARCODE_MISMATCHES}"
    
    # build bcl-convert
    CMD="bcl-convert --bcl-input-directory ${INPUT_FOLDER} --output-directory ${OUTPUT_FOLDER}"
    CMD="${CMD} --bcl-num-compression-threads ${ILLUMINA_THREADS} --bcl-num-conversion-threads ${ILLUMINA_THREADS}"
    CMD="${CMD} --sample-name-column-enabled true"
    
    # Add --no-lane-splitting if TRUE
    if [ "${NO_LANE_SPLITTING_ILLUMINA}" = "TRUE" ]; then
        CMD="${CMD} --no-lane-splitting true"
    fi
else
    echo "ERROR: Unknown instrument type '${INSTRUMENT}'. Must be 'AVITTI' or 'Novaseq'."
    exit 1
fi

# Log the command and execute
echo "Running command:"
echo "${CMD}"
echo "Start time: $(date)"

# Execute the command
eval "${CMD}"

# Check exit status
if [ $? -eq 0 ]; then
    echo "Successfully completed at: $(date)"
    
    # Only consolidate for Illumina runs
#    if [ "${INSTRUMENT}" = "ILLUMINA" ]; then
        # Extract the subfolder name (e.g., P935_20250423_JM_DB_AG)
        # This assumes the Illumina output structure has a subfolder like /1/, /2/, etc.
        # We take the first existing subfolder's parent name
#        FIRST_SAMPLE_DIR=$(find "${OUTPUT_FOLDER}" -mindepth 1 -maxdepth 1 -type d -name "${PROJECT_NUMBER}_*" | head -n 1)
#        SUBFOLDER_NAME=$(basename "${FIRST_SAMPLE_DIR%/*}" 2>/dev/null || basename "${OUTPUT_FOLDER}")

        # Define the conc.fastq path (e.g., OUTPUT_FOLDER/P935_20250423_JM_DB_AG/conc.fastq)
#        CONC_FOLDER="${FIRST_SAMPLE_DIR}/conc.fastq"
#        mkdir -p "${CONC_FOLDER}"

        # Move all FASTQ files to conc.fastq directory
#        echo "Moving all FASTQ files to ${CONC_FOLDER}"
#        find "${OUTPUT_FOLDER}" -type f \( -name "*.fastq.gz" -o -name "*.fastq" \) ! -name "Undetermined*" -exec mv {} "${CONC_FOLDER}" \;

        # Remove empty sample directories (optional)
#        echo "Cleaning up empty sample directories"
#        find "${OUTPUT_FOLDER}" -type d -empty -delete

#        echo "All FASTQ files have been consolidated in: ${CONC_FOLDER}"
#    fi
else
    echo "Error encountered during processing. Check log files."
    exit 1
fi

# Print final output location
echo "Output files created in: ${OUTPUT_FOLDER}"

#!/bin/bash
#SBATCH --export=ALL
#SBATCH -p mrcq 
#SBATCH --time=03:00:00 
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
#SBATCH --mem=4G 
#SBATCH --mail-type=END 
#SBATCH --output=oxBS_5mc%j.log
#SBATCH --error=oxBS_5mc%j.err
#SBATCH --job-name=oxBS_5mc

usage() {
cat <<EOF
================================================================================
2_binarize_oxBS_5mC_data.sh
================================================================================
Purpose: Create binary files for dense and sparse regions of 5mC from 
oxidative bisulphite sequencing bed files. Bed files are expected to be in the 
format:
  chr \t start \t end \t 5mC_reads \t total_reads  
Author: Sam Fletcher
Contact: s.o.fletcher@exeter.ac.uk
Dependencies: R, awk, bedtools
================================================================================
EOF
    exit 0
}

if [ "$#" -eq 0 ]; then usage; fi 

config_file_location=$1
source "${config_file_location}" || { echo "could not find config file at:
${config_file_location}"; exit 1; }

source "${REPO_DIR}/parameters.txt" || exit 1

for file in "${FUNCTIONS_DIR}"/*; do source "$file" || exit 1; done

move_log_files oxBS_5mc

## =================================== ##
##   EXTRACT HYDROXYMETHYLATED SITES   ##
## =================================== ##

processing_directory="${BASE_DIR}/oxBS_5mc"

rm -rf "${processing_directory}"
mkdir -p "${processing_directory}"

purification_convertBSBedToMethylBedFormat \
  "m" \
  "${oxBS_bed_file_location}" \
  "${processing_directory}/formatted.bed"

purification_extractSitesWithLowMethylation \
  "m" \
  "${processing_directory}/formatted.bed" \
  "${processing_directory}/unmethylated_reads.bed"
purification_filterOutLowReadDepthSites \
  "m" \
  "${processing_directory}/formatted.bed" \
  "${processing_directory}/filtered_reads.bed"
purification_calculateSiteMethylationProbability \
  "${processing_directory}" \
  "unmethylated_reads.bed" \
  "filtered_reads.bed" \
  "processed_reads.bed"
purification_removeDeterminedUnmethylatedSites \
  "${processing_directory}/processed_reads.bed" \
  "${processing_directory}/purified_reads.bed"

## ======================== ##
##   BINARIZATION PROCESS   ##
## ======================== ##

binarization_createDirectories \
  "${processing_directory}"
binarization_splitIntoChromosomes \
  "${processing_directory}" \
  "purified_reads.bed"
binarization_createBlankBins \
  "${processing_directory}"
binarization_countSignalIntersectionWithBins \
  "${processing_directory}"
binarization_createChromhmmBinaryFiles \
  "${processing_directory}" \
  "${BINARY_DIR}/oxBS_5mC" \
  "oxBS_5mC"

if [[ ! "${DEBUG_MODE:=0}" -eq 1 ]]; then
  rm -rf "${processing_directory}"
fi

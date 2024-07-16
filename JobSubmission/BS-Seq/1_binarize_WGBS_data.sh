#!/bin/bash
#SBATCH --export=ALL
#SBATCH -p mrcq 
#SBATCH --time=03:00:00 
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
#SBATCH --mem=4G 
#SBATCH --mail-type=END 
#SBATCH --output=WGBS_5mc_5hmc%j.log
#SBATCH --error=WGBS_5mc_5hmc%j.err
#SBATCH --job-name=WGBS_5mc_5hmc

usage() {
cat <<EOF
================================================================================
1_binarize_WGBS_data.sh
================================================================================
Purpose: Create binary files for dense and sparse regions of methylation 
from whole genome bisulphite sequencing bed files. Bed files are expected to be 
in the format:
  chr \t start \t end \t methylated_reads \t total_reads  
Author: Sam Fletcher
Contact: s.o.fletcher@exeter.ac.uk
Dependencies: R, awk, bedtools
================================================================================
EOF
    exit 0
}

if [ "$#" -eq 0 ]; then usage; fi 

config_file_location=$1
source "${config_file_location}" || { usage; exit 1; }

source "${REPO_DIR}/parameters.txt" || exit 1

source "${FUNCTIONS_DIR}/move_log_files.sh" || exit 1
move_log_files WGBS_5mc_5hmc

processing_directory="${BASE_DIR}/WGBS_5mc_5hmc"

## =================================== ##
##   EXTRACT HYDROXYMETHYLATED SITES   ##
## =================================== ##

rm -rf "${processing_directory}"
mkdir -p "${processing_directory}"
source "${FUNCTIONS_DIR}/purification.sh" || exit 1

purification_convertBSBedToMethylBedFormat "${processing_directory}/formatted.bed" \
  "${WGBS_bed_file_location}" \
  "mh"

purification_extractSitesWithHighMethylation "${processing_directory}" \
  "${processing_directory}/formatted.bed" \
  "mh"
purification_extractSitesWithLowMethylation "${processing_directory}" \
  "${processing_directory}/formatted.bed" \
  "mh"
purification_filterOutLowReadDepthSites "${processing_directory}" "${processing_directory}/formatted.bed" "h"
purification_calculateSiteMethylationProbability "${processing_directory}" 
purification_removeDeterminedUnmethylatedSites "${processing_directory}"

## ======================== ##
##   BINARIZATION PROCESS   ##
## ======================== ##

source "${FUNCTIONS_DIR}/binarization.sh" || exit 1

binarization_createDirectories "${processing_directory}"
binarization_splitIntoChromosomes "${processing_directory}"
binarization_createBlankBins "${processing_directory}"
binarization_countSignalIntersectionWithBins "${processing_directory}"

binarization_createChromhmmBinaryFiles "${processing_directory}" \
  "${BINARY_DIR}/WGBS_5mC_5hmC" \
  "WGBS_5mC_5hmC"

if [[ ! "${debug_mode:='false'}" == "true" ]]; then
  rm -rf "${processing_directory}"
fi
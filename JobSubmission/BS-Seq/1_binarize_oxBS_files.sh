#!/bin/bash
#SBATCH --export=ALL
#SBATCH -p mrcq 
#SBATCH --time=03:00:00 
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
#SBATCH --mem=4G 
#SBATCH --mail-type=END 
#SBATCH --output=hydroxy%j.log
#SBATCH --error=hydroxy%j.err
#SBATCH --job-name=hydroxy

usage() {
cat <<EOF
================================================================================
1_binarize_oxBS_files.sh
================================================================================
Purpose: Create binary files for dense and sparse regions of hydroxymethylation 
from oxidative bisulphite sequencing bed files. Bed files are expected to be 
in the format:
  chr \t start \t end \t hydroxymethylated_reads \t total_reads  
Author: Sam Fletcher
Contact: s.o.fletcher@exeter.ac.uk
Dependencies: R, awk, bedtools
================================================================================
EOF
    exit 0
}

if [ "$#" -eq 0 ]; then usage; fi 

config_file_location=$1
source "${config_file_location}" || exit 1

source "${ROOT_DIR}/parameters.txt" || exit 1

source "${FUNCTIONS_DIR}/move_log_files.sh" || exit 1
move_log_files hydroxy

output_directory="${BASE_DIR}/5hmc"

## =================================== ##
##   EXTRACT HYDROXYMETHYLATED SITES   ##
## =================================== ##

mkdir -p "${output_directory}"
source "${FUNCTIONS_DIR}/purification.sh" || exit 1

purification_convertBSBedToMethylBedFormat "${output_directory}/formatted.bed"

purification_extractSitesWithHighMethylation "${output_directory}" "${output_directory}/formatted.bed" "h"
purification_extractSitesWithLowMethylation "${output_directory}" "${output_directory}/formatted.bed" "h"
purification_filterOutLowReadDepthSites "${output_directory}" "${output_directory}/formatted.bed" "h"
purification_calculateSiteMethylationProbability "${output_directory}" 
purification_removeDeterminedUnmethylatedSites "${output_directory}"

## ======================== ##
##   BINARIZATION PROCESS   ##
## ======================== ##

source "${FUNCTIONS_DIR}/binarization.sh" || exit 1

binarization_createDirectories "${output_directory}"
binarization_splitIntoChromosomes "${output_directory}"
binarization_createBlankBins "${output_directory}"
binarization_countSignalIntersectionWithBins "${output_directory}"
binarization_createChromhmmBinaryFiles "${output_directory}"


#!/bin/bash
#SBATCH --export=ALL
#SBATCH -p mrcq 
#SBATCH --time=03:00:00 
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
#SBATCH --mem=4G 
#SBATCH --mail-type=END 
#SBATCH --output=binomial%j.log
#SBATCH --error=binomial%j.err
#SBATCH --job-name=Binomial

usage() {
cat <<EOF
================================================================================
1_purify_reads.sh
================================================================================
Purpose: Filters input bed file on sites that are significantly (un)methylated
Author: Sam Fletcher
Contact: s.o.fletcher@exeter.ac.uk
Dependencies: R, awk
================================================================================
EOF
    exit 0
}

if [ "$#" -eq 0 ]; then usage; fi 

## ======== ##
##   MAIN   ##
## ======== ##

config_file_location=$1
source "${config_file_location}" || exit 1

source "${ROOT_DIR}/parameters.txt" || exit 1

source "${FUNCTIONS_DIR}/move_log_files.sh" || exit 1
move_log_files binomial

rm -rf "${BASE_DIR}/5mc" "${BASE_DIR}/5hmc"
mkdir -p "${BASE_DIR}/5mc" "${BASE_DIR}/5hmc"

source "${FUNCTIONS_DIR}/purification.sh" || exit 1

## ==== ##
##  5mC ##
## ==== ##

purification_extractSitesWithHighMethylation "${BASE_DIR}/5mc" "${bed_file_location}" "m"
purification_extractSitesWithLowMethylation "${BASE_DIR}/5mc" "${bed_file_location}" "m"
purification_filterOutLowReadDepthSites "${BASE_DIR}/5mc" "${bed_file_location}" "m"
purification_calculateSiteMethylationProbability "${BASE_DIR}/5mc"
purification_removeDeterminedUnmethylatedSites "${BASE_DIR}/5mc"

## ===== ##
##  5hmC ##
## ===== ##

purification_extractSitesWithHighMethylation "${BASE_DIR}/5hmc" "${bed_file_location}" "h"
purification_extractSitesWithLowMethylation "${BASE_DIR}/5hmc" "${bed_file_location}" "h"
purification_filterOutLowReadDepthSites "${BASE_DIR}/5hmc" "${bed_file_location}" "h"
purification_calculateSiteMethylationProbability "${BASE_DIR}/5hmc"
purification_removeDeterminedUnmethylatedSites "${BASE_DIR}/5hmc"


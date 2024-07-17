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
Purpose: Filters input ONT file on sites that are significantly methylated
Author: Sam Fletcher
Contact: s.o.fletcher@exeter.ac.uk
Dependencies: R, awk
================================================================================
EOF
    exit 0
}

if [ "$#" -eq 0 ]; then usage; fi 

config_file_location=$1
source "${config_file_location}" || usage

source "${REPO_DIR}/parameters.txt" || exit 1

for file in "${FUNCTIONS_DIR}"/*; do source "$file" || exit 1; done

move_log_files binomial

rm -rf "${BASE_DIR}/5mc" "${BASE_DIR}/5hmc"
mkdir -p "${BASE_DIR}/5mc" "${BASE_DIR}/5hmc"

## ==== ##
##  5mC ##
## ==== ##

processing_directory="${BASE_DIR}/5mc"

purification_extractSitesWithLowMethylation \
  "m" \
  "${ONT_bed_file_location}" \
  "${processing_directory}/unmethylated_reads.bed"
purification_filterOutLowReadDepthSites \
  "m" \
  "${ONT_bed_file_location}" \
  "${processing_directory}/filtered_reads.bed" 
purification_calculateSiteMethylationProbability \
  "${processing_directory}" \
  "unmethylated_reads.bed" \
  "filtered_reads.bed" \
  "processed_reads.bed"
purification_removeDeterminedUnmethylatedSites \
  "${processing_directory}/processed_reads.bed" \
  "${processing_directory}/purified_reads.bed"

## ===== ##
##  5hmC ##
## ===== ##

processing_directory="${BASE_DIR}/5hmc"

purification_extractSitesWithLowMethylation \
  "h" \
  "${ONT_bed_file_location}" \
  "${processing_directory}/unmethylated_reads.bed"
purification_filterOutLowReadDepthSites \
  "h" \
  "${ONT_bed_file_location}" \
  "${processing_directory}/filtered_reads.bed" 
purification_calculateSiteMethylationProbability \
  "${processing_directory}" \
  "unmethylated_reads.bed" \
  "filtered_reads.bed" \
  "processed_reads.bed"
purification_removeDeterminedUnmethylatedSites \
  "${processing_directory}/processed_reads.bed" \
  "${processing_directory}/purified_reads.bed"

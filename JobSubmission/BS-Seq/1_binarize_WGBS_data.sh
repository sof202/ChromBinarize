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
Optional argument: -c -> train erroneous reads probability on CpGs in CGIs only
Author: Sam Fletcher
Contact: s.o.fletcher@exeter.ac.uk
Dependencies: R, awk, bedtools
================================================================================
EOF
    exit 0
}

if [ "$#" -eq 0 ]; then usage; fi 

while getopts c OPT; do
    case "$OPT" in
        c )       use_cpg_islands="TRUE" ;;
        * )       usage ;;
    esac
done
shift $((OPTIND-1))

config_file_location=$1
source "${config_file_location}" || { echo "could not find config file at:
${config_file_location}"; exit 1; }

for file in "${FUNCTIONS_DIR}"/*; do source "$file" || exit 1; done

move_log_files WGBS_5mc_5hmc

## =================================== ##
##   EXTRACT HYDROXYMETHYLATED SITES   ##
## =================================== ##

processing_directory="${BASE_DIR}/WGBS_5mc_5hmc"

rm -rf "${processing_directory}"
mkdir -p "${processing_directory}"

purification_convertBSBedToMethylBedFormat \
  "m" \
  "${WGBS_bed_file_location}" \
  "${processing_directory}/formatted.bed"

if [[ -n "${use_cpg_islands}" ]]; then
  mv "${processing_directory}/formatted.bed" \
    "${processing_directory}/formatted_all_cpgs.bed"

  purification_extractSitesInCpGIslands \
    "${processing_directory}/formatted_all_cpgs.bed" \
    "${processing_directory}/formatted.bed"
fi

purification_calculateSiteMethylationProbability \
  "m" \
  "${processing_directory}/formatted.bed" \
  "processed_reads.bed"
purification_filterOnReadDepth \
  "m" \
  "${processing_directory}/processed_reads.bed" \
  "${processing_directory}/filtered_reads.bed" 
purification_removeDeterminedUnmethylatedSites \
  "${processing_directory}/filtered_reads.bed" \
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
  "${processing_directory}" \
  "${BIN_SIZE}"
binarization_countSignalIntersectionWithBins \
  "${processing_directory}"
binarization_createChromhmmBinaryFiles \
  "${processing_directory}" \
  "${BINARY_DIR}/WGBS_5mC_5hmC" \
  "WGBS_5mC_5hmC"

if [[ "${DEBUG_MODE:=0}" -eq 0 ]]; then
  rm -rf "${processing_directory}"
fi

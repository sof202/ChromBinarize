#!/bin/bash
#SBATCH --export=ALL
#SBATCH -p mrcq 
#SBATCH --time=03:00:00 
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
#SBATCH --mem=4G 
#SBATCH --mail-type=END 
#SBATCH --output=methylation%j.log
#SBATCH --error=methylation%j.err
#SBATCH --job-name=methylation

usage() {
cat <<EOF
================================================================================
2_binarize_WGBS_files.sh
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

# config will source all of the variables seen below
config_file_location=$1
source "${config_file_location}" || { usage; exit 1; } 

source "${REPO_DIR}/parameters.txt" || exit 1

source "${FUNCTIONS_DIR}/move_log_files.sh" || exit 1
move_log_files methylation

processing_directory="${BASE_DIR}/WGBS_5mc"

## ============================ ##
##   EXTRACT METHYLATED SITES   ##
## ============================ ##

## ------------------------------------ ##
##   REMOVE HYDROXYMETHYLATION SIGNAL   ##
## ------------------------------------ ##

# This requires oxBS and WGBS files to work. WGBS captures 5mC AND 5hmC signal.
# We only want to capture the 5hmC signal, so we need to remove the 5mc signal
# using oxBS data.

rm -rf "${processing_directory}"
mkdir -p "${processing_directory}"

if [[ -n ${oxBS_bed_file_location} ]]; then
  module purge
  module load BEDTools

  bedtools intersect -wo \
    -a "${WGBS_bed_file_location}" \
    -b "${oxBS_bed_file_location}" |
    awk \
      -v read_threshold="${reference_read_depth_threshold_h}" \
      'function convert_to_percent(reads, total_reads) {
         return (int(reads / total_reads * 10000) / 100)
       }
       {OFS="\t"} 
       $10 >= read_threshold {print $1,$2,$3,$5,convert_to_percent($4,$5),convert_to_percent($9,$10)}' > \
    "${processing_directory}/combined.bed"

  awk '
      function ReLU_distance(i,j) {
        return ((i - j) > 0 ? (i - j) : 0)
      }
      {OFS="\t"}
      {print $1,$2,$3,"m",$4,"+",ReLU_distance($5,$6)}
      ' "${processing_directory}/combined.bed" > \
        "${processing_directory}/WGBS_5mc_removed.bed"
else
  cp "${WGBS_bed_file_location}" "${processing_directory}/WGBS_5mc_removed.bed"
fi

source "${FUNCTIONS_DIR}/purification.sh" || exit 1

purification_extractSitesWithHighMethylation "${processing_directory}" "${processing_directory}/WGBS_5mc_removed.bed" "m"
purification_extractSitesWithLowMethylation "${processing_directory}" "${processing_directory}/WGBS_5mc_removed.bed" "m"
purification_filterOutLowReadDepthSites "${processing_directory}" "${processing_directory}/WGBS_5mc_removed.bed" "m"
purification_calculateSiteMethylationProbability "${processing_directory}"
purification_removeDeterminedUnmethylatedSites "${processing_directory}"

## ======================== ##
##   BINARIZATION PROCESS   ##
## ======================== ##
#
source "${FUNCTIONS_DIR}/binarization.sh" || exit 1

binarization_createDirectories "${processing_directory}"
binarization_splitIntoChromosomes "${processing_directory}"
binarization_createBlankBins "${processing_directory}"
binarization_countSignalIntersectionWithBins "${processing_directory}"
binarization_createChromhmmBinaryFiles "${processing_directory}" "${BINARY_DIR}/WGBS_5mC" "WGBS_5mC"

if [[ ! "${debug_mode:='false'}" == "true" ]]; then
  rm -rf "${processing_directory}"
fi

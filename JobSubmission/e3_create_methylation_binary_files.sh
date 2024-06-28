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
e2_create_methylation_binary_files.sh
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
source "${config_file_location}" || exit 1

source "${ROOT_DIR}/parameters.txt" || exit 1

source "${FUNCTIONS_DIR}/move_log_files.sh" || exit 1
move_log_files methylation

## ============================ ##
##   EXTRACT METHYLATED SITES   ##
## ============================ ##

## ------------------------------------ ##
##   REMOVE HYDROXYMETHYLATION SIGNAL   ##
## ------------------------------------ ##

# This requires oxBS and WGBS files to work. WGBS captures 5mC AND 5hmC signal.
# We only want to capture the 5mC signal, so we need to remove the 5hmC signal
# using oxBS data.

mkdir -p "${BASE_DIR}/5mc"

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
    "${BASE_DIR}/5mc/combined.bed"

  awk '
      function convert_to_reads(percentage, total_reads) {
        return (int(percentage * total_reads / 100))
      }
      function ReLU_distance(i,j) {
        return ((i - j) > 0 ? (i - j) : 0)
      }
      {OFS="\t"}
      {print $1,$2,$3,convert_to_reads(ReLU_distance($5,$10), $4), $4}
      ' "${BASE_DIR}/5mc/combined.bed" > \
        "${BASE_DIR}/5mc/WGBS_5hmc_removed.bed"
else
  cp "${WGBS_bed_file_location}" "${BASE_DIR}/5mc/WGBS_5hmc_removed.bed"
fi

## ------------------------------- ##
##   EXTRACT CONFIDENT POSITIONS   ##
## --------------------------------##

awk -v percent_threshold="${reference_percentage_threshold_h}" \
  -v read_threshold="${reference_read_depth_threshold_h}" \
  'function convert_to_percent(reads, total_reads) {
     return (int(reads / total_reads * 10000) / 100)
   }
   $5 >= read_threshold && convert_to_percent($4,$5) >= percent_threshold {print $5","convert_to_percent($4,$5)}' \
  "${oxBS_bed_file_location}" > \
    "${BASE_DIR}/5mc/methylated.csv"

awk -v percent_threshold=$((100 - ${reference_percentage_threshold_h})) \
  -v read_threshold="${reference_read_depth_threshold_h}" \
  'function convert_to_percent(reads, total_reads) {
     return (int(reads / total_reads * 10000) / 100)
   }
   $5 >= read_threshold && convert_to_percent($4,$5) >= percent_threshold {print $5","convert_to_percent($4,$5)}' \
  "${oxBS_bed_file_location}" > \
    "${BASE_DIR}/5mc/unmethylated.csv"

awk -v read_threshold="${minimum_read_depth}" \
  'function convert_to_percent(reads, total_reads) {
     return (int(reads / total_reads * 10000) / 100)
   }
  {OFS="\t"} 
  $5 >= read_threshold {print $1,$2,$3,"h",$5,"+",convert_to_percent($4,$5)}' \
  "${oxBS_bed_file_location}" > \
    "${BASE_DIR}/5mc/filtered_reads.bed"


## ------------------------- ##
##   RUN BINOMIAL ANALYSIS   ##
## ------------------------- ##

module purge
module load R/4.2.1-foss-2022a

Rscript "${RSCRIPT_DIR}/binom.R" "${BASE_DIR}/5mc"

module purge

## ----------------------------- ##
##   REMOVE UNMETHYLATED SITES   ##
## ----------------------------- ##

awk -v threshold="${binomial_threshold}" \
  '$9 < threshold' \
  "${BASE_DIR}/5mc/processed_reads.bed" > \
  "${BASE_DIR}/5mc/purified_reads.bed"

## ======================== ##
##   BINARIZATION PROCESS   ##
## ======================== ##
#
source "${FUNCTIONS_DIR}/binarization.sh" || exit 1

binarization_createDirectories "${BASE_DIR}/5mc"
binarization_splitIntoChromosomes "${BASE_DIR}/5mc"
binarization_createBlankBins "${BASE_DIR}/5mc"
binarization_countSignalIntersectionWithBins "${BASE_DIR}/5mc"
binarization_createChromhmmBinaryFiles "${BASE_DIR}/5mc"


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
e2_create_hydroxy_binary_files.sh
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

## =================================== ##
##   EXTRACT HYDROXYMETHYLATED SITES   ##
## =================================== ##

## ------------------------------- ##
##   EXTRACT CONFIDENT POSITIONS   ##
## --------------------------------##

mkdir -p "${BASE_DIR}/5hmc"

awk -v percent_threshold="${reference_percentage_threshold_h}" \
  -v read_threshold="${reference_read_depth_threshold_h}" \
  'function convert_to_percent(reads, total_reads) {
     return (int(reads / total_reads * 10000) / 100)
   }
   $5 >= read_threshold && convert_to_percent($4,$5) >= percent_threshold {print $5","convert_to_percent($4,$5)}' \
  "${oxBS_bed_file_location}" > \
    "${BASE_DIR}/5hmc/methylated.csv"

awk -v percent_threshold=$((100 - ${reference_percentage_threshold_h})) \
  -v read_threshold="${reference_read_depth_threshold_h}" \
  'function convert_to_percent(reads, total_reads) {
     return (int(reads / total_reads * 10000) / 100)
   }
   $5 >= read_threshold && convert_to_percent($4,$5) >= percent_threshold {print $5","convert_to_percent($4,$5)}' \
  "${oxBS_bed_file_location}" > \
    "${BASE_DIR}/5hmc/unmethylated.csv"

awk -v read_threshold="${minimum_read_depth}" \
  'function convert_to_percent(reads, total_reads) {
     return (int(reads / total_reads * 10000) / 100)
   }
  {OFS="\t"} 
  $5 >= read_threshold {print $1,$2,$3,"h",$5,"+",convert_to_percent($4,$5)}' \
  "${oxBS_bed_file_location}" > \
    "${BASE_DIR}/5hmc/filtered_reads.bed"

## ------------------------- ##
##   RUN BINOMIAL ANALYSIS   ##
## ------------------------- ##

module purge
module load R/4.2.1-foss-2022a

Rscript "${RSCRIPT_DIR}/binom.R" "${BASE_DIR}/5hmc"

module purge

## -------------------------==== ##
##   REMOVE UNMETHYLATED SITES   ##
## -------------------------==== ##

awk -v threshold="${binomial_threshold}" \
  '$9 < threshold' \
  "${BASE_DIR}/5hmc/processed_reads.bed" > \
  "${BASE_DIR}/5hmc/purified_reads.bed"

## ======================== ##
##   BINARIZATION PROCESS   ##
## ======================== ##

source "${FUNCTIONS_DIR}/binarization.sh" || exit 1

binarization_createDirectories "${BASE_DIR}/5hmc"
binarization_splitIntoChromosomes "${BASE_DIR}/5hmc"
binarization_createBlankBins "${BASE_DIR}/5hmc"
binarization_countSignalIntersectionWithBins "${BASE_DIR}/5hmc"
binarization_createChromhmmBinaryFiles "${BASE_DIR}/5hmc"


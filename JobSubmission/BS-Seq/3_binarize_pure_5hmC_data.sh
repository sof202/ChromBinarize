#!/bin/bash
#SBATCH --export=ALL
#SBATCH -p mrcq 
#SBATCH --time=03:00:00 
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
#SBATCH --mem=4G 
#SBATCH --mail-type=END 
#SBATCH --output=5mC%j.log
#SBATCH --error=5mC%j.err
#SBATCH --job-name=5mC

usage() {
cat <<EOF
================================================================================
3_binarize_pure_5hmC_data.
================================================================================
Purpose: Create binary files for dense and sparse regions of 5hmc using a
combination of whole genome bisulphite sequencing and oxidative bisulphite
sequencing bed files. Bed files are expected to be in the format:
  chr \t start \t end \t number_of_methylated_reads \t total_reads  
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

for file in "${FUNCTIONS_DIR}"/*; do source "$file" || exit 1; done

move_log_files 5mC

## =================================== ##
##   EXTRACT HYDROXYMETHYLATED SITES   ##
## =================================== ##

## ------------- ##
##   FILTERING   ##
## ------------- ##

processing_directory="${BASE_DIR}/WGBS_5hmc"

rm -rf "${processing_directory}"
mkdir -p "${processing_directory}"

# We remove WGBS reads with 0% methylation as there cannot be hydroxymethylation
# here (saves computational time later)
awk -v min_read_depth="${minimum_read_depth}" \
  -v max_read_depth="${maximum_read_depth}" \
  '$5 >= min_read_depth && $5 <= max_read_depth && $4 > 0' \
  "${WGBS_bed_file_location}" > \
  "${processing_directory}/WGBS_filtered.bed"

if [[ ! -s "${processing_directory}/WGBS_filtered.bed" ]]; then
errors "${processing_directory}/WGBS_filtered.bed is empty. 
The most likely cause of this is your read threshold is too strict."
fi

# We remove oxBS reads with 100% methylation as there cannot be
# hydroxymethylation here (it is shown to be 100% methylation). This saves
# computational time later.
awk -v min_read_depth="${minimum_read_depth}" \
  -v max_read_depth="${maximum_read_depth}" \
  '$5 >= min_read_depth && $5 <= max_read_depth && $4 < $5' \
  "${oxBS_bed_file_location}" > \
  "${processing_directory}/oxBS_filtered.bed"

if [[ ! -s "${processing_directory}/oxBS_filtered.bed" ]]; then
errors "${processing_directory}/oxBS_filtered.bed is empty. 
The most likely cause of this is your read threshold is too strict."
fi

## -------------------------------- ##
##   GENERATE CONFIDENCE INTERVAL   ##
## -------------------------------- ##

module purge
module load R/4.2.1-foss-2022a

za=$(Rscript -e "cat(qnorm(1 - ${confidence_interval_alpha:=0.05}/2))")

module purge

# We use Agresti-Coull interval for obtaining the binomial confidence interval 
# as it is inexpensive whilst still providing good coverage (i.e in pratical 
# cases, it has been shown this interval contains the true value roughly as 
# common as alpha suggests).
awk \
  -v za="$za" \
  -v is_wgbs=1 \
  -f "${AWK_DIR}/Agresti_Coull.awk" \
  "${processing_directory}/WGBS_filtered.bed" > \
  "${processing_directory}/WGBS_Agresti_Coull.bed"

awk \
  -v za="$za" \
  -v is_wgbs=0 \
  -f "${AWK_DIR}/Agresti_Coull.awk" \
  "${processing_directory}/oxBS_filtered.bed" > \
  "${processing_directory}/oxBS_Agresti_Coull.bed"

## --------------------------- ##
##   INTERSECT OXBS AND WGBS   ##
## --------------------------- ##

module load BEDTools

bedtools intersect -wo \
  -a "${processing_directory}/WGBS_Agresti_Coull.bed" \
  -b "${processing_directory}/oxBS_Agresti_Coull.bed" > \
  "${processing_directory}/WGBS_oxBS_combined.bed"

module purge

## ----------------------------------- ##
##   EXTRACT HYDROXYMETHYLATED SITES   ##
## ----------------------------------- ##

# For a significant hydroxymethylation signal, the lower confidence interval
# bound for WGBS signal must be greater than the upper confidence interval 
# bound for the oxBS signal [see README.md]
awk \
  '{OFS="\t"}
  $6 > $12 {print $1,$2,$3,$4,$5}' \
  "${processing_directory}/WGBS_oxBS_combined.bed" > \
  "${processing_directory}/WGBS_5mC_removed.bed"

purification_convertBSBedToMethylBedFormat \
  "h" \
  "${processing_directory}/WGBS_5mC_removed.bed" \
  "${processing_directory}/purified_reads.bed"

if [[ ! -s "${processing_directory}/purified_reads.bed" ]]; then
errors "${processing_directory}/purified_reads.bed is empty. 
The most likely cause of this is your significance level is too small."
fi

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
  "${BINARY_DIR}/WGBS_5hmC" \
  "WGBS_5hmC"

if [[ "${DEBUG_MODE:=0}" -eq 0 ]]; then
  rm -rf "${processing_directory}"
fi

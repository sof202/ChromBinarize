#!/bin/bash
#SBATCH --export=ALL
#SBATCH -p mrcq 
#SBATCH --time=03:00:00 
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
#SBATCH --mem=20G 
#SBATCH --mail-type=END 
#SBATCH --output=robustness%j.log
#SBATCH --error=robustness%j.err
#SBATCH --job-name=robustness

usage() {
cat <<EOF
================================================================================
b_CpG_robustness.sh
================================================================================
Purpose: Calculates the correlation of percent methylation seen in local 
clusters of CpGs. High correlation suggests robustness in basecalling.
Author: Sam Fletcher
Contact: s.o.fletcher@exeter.ac.uk
Dependencies: R
================================================================================
EOF
    exit 0
}

if [ "$#" -eq 0 ]; then usage; fi 


# config will source all of the variables seen below
config_file_location=$1
source "${config_file_location}" || usage 

source "${REPO_DIR}/parameters.txt" || exit 1

source "${FUNCTIONS_DIR}/move_log_files.sh" || exit 1
move_log_files robustness

filtered_reads_directory="${BASE_DIR}/CpG_robustness"
rm -rf "${filtered_reads_directory}"
mkdir -p "${filtered_reads_directory}"

source "${FUNCTIONS_DIR}/purification.sh" || exit 1

number_of_columns=$(awk '{print NF; exit}' "${bed_file_location}")
if [[ "${number_of_columns}" -eq 5 ]]; then
    purification_convertBSBedToMethylBedFormat \
      "${mark}" \
      "${bed_file_location}" \
      "${filtered_reads_directory}/converted.bed"

    bed_file_location="${filtered_reads_directory}/converted.bed"
fi

purification_filterOutLowReadDepthSites \
  "${mark}" \
  "${bed_file_location}" \
  "${filtered_reads_directory}"

## ============ ##
##   PLOTTING   ##
## ============ ##

mkdir -p "${BASE_DIR}/plots"

module purge
module load R/4.2.1-foss-2022a

Rscript "${RSCRIPT_DIR}/CpG_robustness.R" \
  "${filtered_reads_directory}/filtered_reads.bed" \
  "${min_distance}" \
  "${max_distance}" \
  "${BASE_DIR}/plots/cpg_robustness_min_${min_distance}_max_${max_distance}.png"

rm -rf "${filtered_reads_directory}"

#!/bin/bash
#SBATCH --export=ALL
#SBATCH -p mrcq 
#SBATCH --time=03:00:00 
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
#SBATCH --mem=4G 
#SBATCH --mail-type=END 
#SBATCH --output=pvalues%j.log
#SBATCH --error=pvalues%j.err
#SBATCH --job-name=pvalues

usage() {
cat <<EOF
================================================================================
a_erroneous_rate_plot.sh
================================================================================
Purpose: Outputs plots that show how approximated erroneous methylation call 
probability changes with read depth/methylation percent thresholds.
Author: Sam Fletcher
Contact: s.o.fletcher@exeter.ac.uk
Dependencies: R
================================================================================
EOF
    exit 0
}

if [ "$#" -eq 0 ]; then usage; fi 

config_file_location=$1
source "${config_file_location}" || { usage; exit 1; }

source "${REPO_DIR}/parameters.txt" || exit 1

source "${FUNCTIONS_DIR}/move_log_files.sh" || exit 1
move_log_files pvalues

number_of_columns=$(awk '{print NF; exit}' "${bed_file_location}")

if [[ "${number_of_columns}" -eq 5 ]]; then
    source "${FUNCTIONS_DIR}/purification.sh" || exit 1

    purification_convertBSBedToMethylBedFormat "${BASE_DIR}/converted.bed" "${bed_file_location}" "${mark}"

    bed_file_location="${BASE_DIR}/converted.bed"
fi

## ======== ##
##   MAIN   ##
## ======== ##

mkdir -p "${BASE_DIR}/plots/"

module purge
module load R/4.2.1-foss-2022a

if [[ "${run_type}" == "N" ]]; then
  Rscript "${RSCRIPT_DIR}/erroneous_rate_plot_N.R" \
    "${bed_file_location}" \
    "${mark}" \
    "${max_N_value}" \
    "${plot_type}" \
    "${BASE_DIR}/plots/erroneous_rate_plot_${plot_type}_${run_type}_x_axis.png"
elif [[ "${run_type}" == "percent" ]]; then
  Rscript "${RSCRIPT_DIR}/erroneous_rate_plot_percent.R" \
    "${bed_file_location}" \
    "${mark}" \
    "${plot_type}" \
    "${BASE_DIR}/plots/erroneous_rate_plot_${plot_type}_${run_type}_x_axis.png"
fi

if [[ "${number_of_columns}" -eq 5 ]]; then
    rm "${BASE_DIR}/converted.bed" 
fi

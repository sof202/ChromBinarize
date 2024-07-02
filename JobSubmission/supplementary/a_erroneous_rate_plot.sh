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
Purpose: Outputs plots that show how p1/p2 change with read depth/methylation
percent thresholds.
Author: Sam Fletcher
Contact: s.o.fletcher@exeter.ac.uk
Dependencies: R
================================================================================
EOF
    exit 0
}

if [ "$#" -eq 0 ]; then usage; fi 

## ======== ##
##   MAIN   ##
## ======== ##

# config will source all of the variables seen below
config_file_location=$1
source "${config_file_location}" || exit 1

source "${ROOT_DIR}/parameters.txt" || exit 1

source "${FUNCTIONS_DIR}/move_log_files.sh" || exit 1
move_log_files pvalues

mkdir -p "$ROOT_DIR/plots/"

module purge
module load R/4.2.1-foss-2022a

if [ "$run_type" == "N" ]; then
  Rscript "$RSCRIPT_DIR/erroneous_rate_plot_N.R" \
    "$bed_file_location" \
    "$mark" \
    "$max_N_value" \
    "$plot_type" \
    "$ROOT_DIR/plots/erroneous_rate_plot_${plot_type}_${run_type}_x_axis.png"
else
  Rscript "$RSCRIPT_DIR/erroneous_rate_plot_percent.R" \
    "$bed_file_location" \
    "$mark" \
    "$plot_type" \
    "$ROOT_DIR/plots/erroneous_rate_plot_${plot_type}_${run_type}_x_axis.png"
fi
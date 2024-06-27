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

SCRIPT_PATH=$(scontrol show job "$SLURM_JOBID" | \
  awk '/Command=/{print $1}' | \
  cut -d= -f1)
SCRIPT_DIR=$(realpath "$(dirname "$SCRIPT_PATH")")

ROOT_DIR="${SCRIPT_DIR}/.."
RSCRIPT_DIR="${ROOT_DIR}/Rscripts"

usage() {
cat <<EOF
================================================================================
a1_erroneous_rate_plot.sh
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

mkdir -p "${LOG_DIR}/"
mv "${SLURM_SUBMIT_DIR}/pvalues${SLURM_JOB_ID}.log" \
  "${LOG_DIR}/pvalues${SLURM_JOB_ID}.log"
mv "${SLURM_SUBMIT_DIR}/pvalues${SLURM_JOB_ID}.err" \
  "${LOG_DIR}/pvalues${SLURM_JOB_ID}.err"

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

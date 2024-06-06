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

mkdir -p "${ROOT_DIR}/logs/"
mv "${SLURM_SUBMIT_DIR}/pvalues${SLURM_JOB_ID}.log" \
  "${ROOT_DIR}/logs/pvalues${SLURM_JOB_ID}.log"
mv "${SLURM_SUBMIT_DIR}/pvalues${SLURM_JOB_ID}.err" \
  "${ROOT_DIR}/logs/pvalues${SLURM_JOB_ID}.err"

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
Inputs:
\$1 -> input bed file
\$2 -> mark (m for 5mc, h for 5hmc)
\$3 -> maximum read depth to consider (in plots)
\$4 -> plot type ("p1" or "p2")
\$5 -> run type ("N" -> Read depth x axis, else -> percent methylation x axis)
================================================================================
EOF
    exit 0
}

if [ -z "$5" ]; then usage; fi 

## ======== ##
##   MAIN   ##
## ======== ##

bed_file=$1
mark=$2 # m for methylation, h for hydroxymethylation
max_read_depth=$3
plot_type=$4
run_type=$5 # N for read_number, otherwise looks at percentage

mkdir -p "$ROOT_DIR/plots/"

module purge
module load R/4.2.1-foss-2022a

if [ "$run_type" == "N" ]; then
  Rscript "$RSCRIPT_DIR/erroneous_rate_plot_N.R" \
    "$bed_file" \
    "$mark" \
    "$max_read_depth" \
    "$plot_type" \
    "$ROOT_DIR/plots/erroneous_rate_plot_${plot_type}_${run_type}_x_axis.png"
else
  Rscript "$RSCRIPT_DIR/erroneous_rate_plot_percent.R" \
    "$bed_file" \
    "$mark" \
    "$plot_type" \
    "$ROOT_DIR/plots/erroneous_rate_plot_${plot_type}_${run_type}_x_axis.png"
fi

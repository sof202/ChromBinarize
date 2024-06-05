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

mkdir -p "${SCRIPT_DIR}/logs/"
mv "${SLURM_SUBMIT_DIR}/pvalues${SLURM_JOB_ID}.log" \
  "${SCRIPT_DIR}/logs/pvalues${SLURM_JOB_ID}.log"
mv "${SLURM_SUBMIT_DIR}/pvalues${SLURM_JOB_ID}.err" \
  "${SCRIPT_DIR}/logs/pvalues${SLURM_JOB_ID}.err"

bed_file=$1
mark=$2 # m for methylation, h for hydroxymethylation
max_read_depth=$3
plot_type=$4
run_type=$5 # N for read_number, otherwise looks at percentage

mkdir -p "$SCRIPT_DIR/plots/"

module purge
module load R/4.2.1-foss-2022a

if [ "$run_type" == "N" ]; then
  Rscript "$SCRIPT_DIR/erroneous_rate_plot_N.R" \
    "$bed_file" \
    "$mark" \
    "$max_read_depth" \
    "$plot_type" \
    "$SCRIPT_DIR/plots/erroneous_rate_plot_N_x_axis.png"
else
  Rscript "$SCRIPT_DIR/erroneous_rate_plot_percent.R" \
    "$bed_file" \
    "$mark" \
    "$plot_type" \
    "$SCRIPT_DIR/plots/erroneous_rate_plot_percent_x_axis.png"
fi

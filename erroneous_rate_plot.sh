#!/bin/bash
#SBATCH --export=ALL
#SBATCH -p mrcq 
#SBATCH --time=03:00:00 
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
#SBATCH --mem=4G 
#SBATCH --mail-type=END 
#SBATCH --output=binom%j.log
#SBATCH --error=binom%j.err
#SBATCH --job-name=Binom

SCRIPT_PATH=$(scontrol show job "$SLURM_JOBID" | \
  awk '/Command=/{print $1}' | \
  cut -d= -f1)
SCRIPT_DIR=$(realpath "$(dirname "$SCRIPT_PATH")")


bed_file=$1
mark=$2
max_read_depth=$3
methylation_threshold=$4


module purge
module load R

Rscript "$SCRIPT_DIR/erroneous_rate_plot.R" \
  "$bed_file" \
  "$mark" \
  "$max_read_depth" \
  "$methylation_threshold"

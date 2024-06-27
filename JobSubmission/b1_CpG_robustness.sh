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

SCRIPT_PATH=$(scontrol show job "$SLURM_JOBID" | \
  awk '/Command=/{print $1}' | \
  cut -d= -f1)
SCRIPT_DIR=$(realpath "$(dirname "$SCRIPT_PATH")")

ROOT_DIR="${SCRIPT_DIR}/.."
RSCRIPT_DIR="${ROOT_DIR}/Rscripts"

usage() {
cat <<EOF
================================================================================
b1_CpG_robustness.sh
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

## ======== ##
##   MAIN   ##
## ======== ##

# config will source all of the variables seen below
config_file_location=$1
source "${config_file_location}" || exit 1

mkdir -p "${LOG_DIR}/"
mv "${SLURM_SUBMIT_DIR}/robustness${SLURM_JOB_ID}.log" \
  "${LOG_DIR}/robustness${SLURM_JOB_ID}.log"
mv "${SLURM_SUBMIT_DIR}/robustness${SLURM_JOB_ID}.err" \
  "${LOG_DIR}/robustness${SLURM_JOB_ID}.err"

mkdir -p "${base_folder}/5mc"
awk -v min_read_threshold="${minimum_read_depth}" \
  -v max_read_threshold="${maximum_read_depth}" \
  -v mark="${mark}" \
  '$4 == mark && $5 >= min_read_threshold && $5 <= max_read_threshold' \
  "${bed_file_location}" > "${base_folder}/5mc/filtered_reads.bed"

mkdir -p "${ROOT_DIR}/plots"

module purge
module load R/4.2.1-foss-2022a

Rscript "${RSCRIPT_DIR}/CpG_robustness.R" \
  "${base_folder}/5mc/filtered_reads.bed" \
  "${min_distance}" \
  "${max_distance}" \
  "${ROOT_DIR}/plots/cpg_robustness_min_${min_distance}_max_${max_distance}.png"



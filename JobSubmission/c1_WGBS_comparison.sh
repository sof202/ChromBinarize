#!/bin/bash
#SBATCH --export=ALL
#SBATCH -p mrcq 
#SBATCH --time=03:00:00 
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
#SBATCH --mem=4G 
#SBATCH --mail-type=END 
#SBATCH --output=compare%j.log
#SBATCH --error=compare%j.err
#SBATCH --job-name=compare

SCRIPT_PATH=$(scontrol show job "$SLURM_JOBID" | \
  awk '/Command=/{print $1}' | \
  cut -d= -f1)
SCRIPT_DIR=$(realpath "$(dirname "$SCRIPT_PATH")")

ROOT_DIR="${SCRIPT_DIR}/.."
RSCRIPT_DIR="${ROOT_DIR}/Rscripts"

usage() {
cat <<EOF
================================================================================
c1_WGBS_comparison.sh
================================================================================
Purpose: Creates several plots comparing a two bed files. One for ONT, the 
other for WGBS data.
Author: Sam Fletcher
Contact: s.o.fletcher@exeter.ac.uk
Dependencies: R, awk, bedtools
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
mv "${SLURM_SUBMIT_DIR}/compare${SLURM_JOB_ID}.log" \
  "${LOG_DIR}/compare${SLURM_JOB_ID}.log"
mv "${SLURM_SUBMIT_DIR}/compare${SLURM_JOB_ID}.err" \
  "${LOG_DIR}/compare${SLURM_JOB_ID}.err"


module purge
module load BEDTools

WGBS_folder=$(dirname "${WGBS_bed_file}")

bedtools intersect \
  -wo \
  -a "${ONT_bed_file}" \
  -b "${WGBS_bed_file}" | \
    awk \
    '{OFS="\t"} {print $1,$2,$3,$4,$5,$7,$12,int($11 / $12 * 10000) / 100}' > \
    "${WGBS_folder}/ONT_WGBS_intersect.bed"

module purge
module load R/4.2.1-foss-2022a

mkdir -p "${ROOT_DIR}/plots"

Rscript "${RSCRIPT_DIR}/WGBS_comparison.R" \
  "${WGBS_folder}/ONT_WGBS_intersect.bed" \
  "${ROOT_DIR}/plots"

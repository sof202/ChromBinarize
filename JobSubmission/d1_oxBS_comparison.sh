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

mkdir -p "${ROOT_DIR}/logs/"
mv "${SLURM_SUBMIT_DIR}/compare${SLURM_JOB_ID}.log" \
  "${ROOT_DIR}/logs/compare${SLURM_JOB_ID}.log"
mv "${SLURM_SUBMIT_DIR}/compare${SLURM_JOB_ID}.err" \
  "${ROOT_DIR}/logs/compare${SLURM_JOB_ID}.err"

usage() {
cat <<EOF
================================================================================
d1_oxBS_comparison.sh
================================================================================
Purpose: Creates several plots comparing a two bed files. One for ONT, the 
other for oxBS data.
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


module purge
module load BEDTools

oxBS_folder=$(dirname "${oxBS_bed_file}")

bedtools intersect \
  -wo \
  -a "${ONT_bed_file}" \
  -b "${oxBS_bed_file}" | \
    awk \
    '{OFS="\t"} {print $1,$2,$3,$4,$5,$7,$12,int($11 / $12 * 10000) / 100}' > \
    "${oxBS_folder}/ONT_oxBS_intersect.bed"

module purge
module load R/4.2.1-foss-2022a

Rscript "${RSCRIPT_DIR}/oxBS_comparison.R" \
  "${oxBS_folder}/ONT_oxBS_intersect.bed" \
  "${oxBS_folder}"

#!/bin/bash
#SBATCH --export=ALL
#SBATCH -p mrcq 
#SBATCH --time=03:00:00 
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
#SBATCH --mem=4G 
#SBATCH --mail-type=END 
#SBATCH --output=binomial%j.log
#SBATCH --error=binomial%j.err
#SBATCH --job-name=Binomial


SCRIPT_PATH=$(scontrol show job "$SLURM_JOBID" | \
  awk '/Command=/{print $1}' | \
  cut -d= -f1)
SCRIPT_DIR=$(realpath "$(dirname "$SCRIPT_PATH")")

ROOT_DIR="${SCRIPT_DIR}/.."
RSCRIPT_DIR="${ROOT_DIR}/Rscripts"

mkdir -p "${ROOT_DIR}/logs/"
mv "${SLURM_SUBMIT_DIR}/binomial${SLURM_JOB_ID}.log" \
  "${ROOT_DIR}/logs/binomial${SLURM_JOB_ID}.log"
mv "${SLURM_SUBMIT_DIR}/binomial${SLURM_JOB_ID}.err" \
  "${ROOT_DIR}/logs/binomial${SLURM_JOB_ID}.err"

usage() {
cat <<EOF
================================================================================
1_binomial.sh
================================================================================
Purpose: Filters input bed file on sites that are significantly (un)methylated
Author: Sam Fletcher
Contact: s.o.fletcher@exeter.ac.uk
Dependencies: R, awk
Inputs:
\$1 -> base folder (for data)
\$2 -> bed file location
================================================================================
EOF
    exit 0
}

if [ -z "$2" ]; then usage; fi 

## ======== ##
##   MAIN   ##
## ======== ##

base_folder=$1
bed_file_location=$2

rm -rf "${base_folder}/5mc" "${base_folder}5hmc"
mkdir -p "${base_folder}/5mc" "${base_folder}5hmc"

awk '$4 == "m" && $5 >= 500 && $7 >= 95 {print $5","$7}' "${bed_file_location}" > "${base_folder}/5mc/methylated.csv"
awk '$4 == "m" && $5 >= 500 && $7 <= 5 {print $5","$7}' "${bed_file_location}" > "${base_folder}/5mc/unmethylated.csv"
awk '$4 == "m" && $5 >= 30' "${bed_file_location}" > "${base_folder}/5mc/filtered_reads.bed"

awk '$4 == "h" && $5 >= 50 && $7 >= 95 {print $5","$7}' "${bed_file_location}" > "${base_folder}/5hmc/methylated.csv"
awk '$4 == "h" && $5 >= 50 && $7 <= 5 {print $5","$7}' "${bed_file_location}" > "${base_folder}/5hmc/unmethylated.csv"
awk '$4 == "h" && $5 >= 30' "${bed_file_location}" > "${base_folder}/5hmc/filtered_reads.bed"

module purge
module load R/4.2.1-foss-2022a

Rscript "${RSCRIPT_DIR}/binom.R" "${base_folder}/5mc"
Rscript "${RSCRIPT_DIR}/binom.R" "${base_folder}/5hmc"

module purge

awk '$8 > 1.7e-8 && $9 < 1.7e-8' "${base_folder}/5mc/processed_reads.bed" > "${base_folder}/5mc/purified_reads.bed"
awk '$8 > 1.7e-8 && $9 < 1.7e-8' "${base_folder}/5hmc/processed_reads.bed" > "${base_folder}/5hmc/purified_reads.bed"

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

usage() {
cat <<EOF
================================================================================
1_binomial.sh
================================================================================
Purpose: Filters input bed file on sites that are significantly (un)methylated
Author: Sam Fletcher
Contact: s.o.fletcher@exeter.ac.uk
Dependencies: R, awk
================================================================================
EOF
    exit 0
}

if [ "$#" -eq 0 ]; then usage; fi 

## ======== ##
##   MAIN   ##
## ======== ##

config_file_location=$1
source "${config_file_location}" || exit 1

source "${ROOT_DIR}/parameters.txt" || exit 1

mkdir -p "${LOG_DIR}/"
mv "${SLURM_SUBMIT_DIR}/binomial${SLURM_JOB_ID}.log" \
  "${LOG_DIR}/binomial${SLURM_JOB_ID}.log"
mv "${SLURM_SUBMIT_DIR}/binomial${SLURM_JOB_ID}.err" \
  "${LOG_DIR}/binomial${SLURM_JOB_ID}.err"

rm -rf "${BASE_DIR}/5mc" "${BASE_DIR}5hmc"
mkdir -p "${BASE_DIR}/5mc" "${BASE_DIR}5hmc"

## ================ ##
##  5mC filtering   ##
## ================ ##

awk -v percent_threshold="${reference_percentage_threshold_m}" \
  -v read_threshold="${reference_read_depth_threshold_m}" \
  '$4 == "m" && $5 >= read_threshold && $7 >= percent_threshold {print $5","$7}' \
  "${bed_file_location}" > "${BASE_DIR}/5mc/methylated.csv"

awk -v percent_threshold=$((100 - ${reference_percentage_threshold_m})) \
  -v read_threshold="${reference_read_depth_threshold_m}" \
  '$4 == "m" && $5 >= read_threshold && $7 <= percent_threshold {print $5","$7}' \
  "${bed_file_location}" > "${BASE_DIR}/5mc/unmethylated.csv"

awk -v read_threshold="${minimum_read_depth}" \
  '$4 == "m" && $5 >= read_threshold' \
  "${bed_file_location}" > "${BASE_DIR}/5mc/filtered_reads.bed"

## ================= ##
##  5hmC filtering   ##
## ================= ##

awk -v percent_threshold="${reference_percentage_threshold_h}" \
  -v read_threshold="${reference_read_depth_threshold_h}" \
  '$4 == "h" && $5 >= read_threshold && $7 >= percent_threshold {print $5","$7}' \
  "${bed_file_location}" > "${BASE_DIR}/5hmc/methylated.csv"

awk -v percent_threshold=$((100 - ${reference_percentage_threshold_h})) \
  -v read_threshold="${reference_read_depth_threshold_h}" \
  '$4 == "h" && $5 >= read_threshold && $7 <= percent_threshold {print $5","$7}' \
  "${bed_file_location}" > "${BASE_DIR}/5hmc/unmethylated.csv"

awk -v read_threshold="${minimum_read_depth}" \
  '$4 == "h" && $5 >= read_threshold' \
  "${bed_file_location}" > "${BASE_DIR}/5hmc/filtered_reads.bed"

## ========================= ##
##   RUN BINOMIAL ANALYSIS   ##
## ========================= ##

module purge
module load R/4.2.1-foss-2022a

Rscript "${RSCRIPT_DIR}/binom.R" "${BASE_DIR}/5mc"
Rscript "${RSCRIPT_DIR}/binom.R" "${BASE_DIR}/5hmc"

module purge

## ============================= ##
##   REMOVE UNMETHYLATED SITES   ##
## ============================= ##

awk -v threshold="${binomial_threshold}" \
  '$9 < threshold' \
  "${BASE_DIR}/5mc/processed_reads.bed" > \
  "${BASE_DIR}/5mc/purified_reads.bed"

awk -v threshold="${binomial_threshold}" \
  '$9 < threshold' \
  "${BASE_DIR}/5hmc/processed_reads.bed" > \
  "${BASE_DIR}/5hmc/purified_reads.bed"

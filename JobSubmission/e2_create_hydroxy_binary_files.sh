#!/bin/bash
#SBATCH --export=ALL
#SBATCH -p mrcq 
#SBATCH --time=03:00:00 
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
#SBATCH --mem=4G 
#SBATCH --mail-type=END 
#SBATCH --output=hydroxy%j.log
#SBATCH --error=hydroxy%j.err
#SBATCH --job-name=hydroxy

SCRIPT_PATH=$(scontrol show job "$SLURM_JOBID" | \
  awk '/Command=/{print $1}' | \
  cut -d= -f1)
SCRIPT_DIR=$(realpath "$(dirname "$SCRIPT_PATH")")

ROOT_DIR="${SCRIPT_DIR}/.."
RSCRIPT_DIR="${ROOT_DIR}/Rscripts"

usage() {
cat <<EOF
================================================================================
e2_create_hydroxy_binary_files.sh
================================================================================
Purpose: Create binary files for dense and sparse regions of hydroxymethylation 
from oxidative bisulphite sequencing bed files. Bed files are expected to be 
in the format:
  chr \t start \t end \t hydroxymethylated_reads \t total_reads  
Author: Sam Fletcher
Contact: s.o.fletcher@exeter.ac.uk
Dependencies: R, awk, bedtools
================================================================================
EOF
    exit 0
}

if [ "$#" -eq 0 ]; then usage; fi 

## =================================== ##
##   EXTRACT HYDROXYMETHYLATED SITES   ##
## =================================== ##

# config will source all of the variables seen below
config_file_location=$1
source "${config_file_location}" || exit 1

mkdir -p "${LOG_DIR}/"
mv "${SLURM_SUBMIT_DIR}/hydroxy${SLURM_JOB_ID}.log" \
  "${LOG_DIR}/hydroxy${SLURM_JOB_ID}.log"
mv "${SLURM_SUBMIT_DIR}/hydroxy${SLURM_JOB_ID}.err" \
  "${LOG_DIR}/hydroxy${SLURM_JOB_ID}.err"

## ------------------------------- ##
##   EXTRACT CONFIDENT POSITIONS   ##
## --------------------------------##

mkdir -p "${base_folder}/5hmc"

awk -v percent_threshold="${reference_percentage_threshold_h}" \
  -v read_threshold="${reference_read_depth_threshold_h}" \
  'function convert_to_percent(reads, total_reads) {
     return (int(reads / total_reads * 10000) / 100)
   }
   $5 >= read_threshold && convert_to_percent($4,$5) >= percent_threshold {print $5","convert_to_percent($4,$5)}' \
  "${oxBS_bed_file_location}" > \
    "${base_folder}/5hmc/methylated.csv"

awk -v percent_threshold=$((100 - ${reference_percentage_threshold_h})) \
  -v read_threshold="${reference_read_depth_threshold_h}" \
  'function convert_to_percent(reads, total_reads) {
     return (int(reads / total_reads * 10000) / 100)
   }
   $5 >= read_threshold && convert_to_percent($4,$5) >= percent_threshold {print $5","convert_to_percent($4,$5)}' \
  "${oxBS_bed_file_location}" > \
    "${base_folder}/5hmc/unmethylated.csv"

awk -v read_threshold="${minimum_read_depth}" \
  'function convert_to_percent(reads, total_reads) {
     return (int(reads / total_reads * 10000) / 100)
   }
  {OFS="\t"} 
  $5 >= read_threshold {print $1,$2,$3,"h",$5,"+",convert_to_percent($4,$5)}' \
  "${oxBS_bed_file_location}" > \
    "${base_folder}/5hmc/filtered_reads.bed"

## ------------------------- ##
##   RUN BINOMIAL ANALYSIS   ##
## ------------------------- ##

module purge
module load R/4.2.1-foss-2022a

Rscript "${RSCRIPT_DIR}/binom.R" "${base_folder}/5hmc"

module purge

## -------------------------==== ##
##   REMOVE UNMETHYLATED SITES   ##
## -------------------------==== ##

awk -v threshold="${binomial_threshold}" \
  '$9 < threshold' \
  "${base_folder}/5hmc/processed_reads.bed" > \
  "${base_folder}/5hmc/purified_reads.bed"

## ======================== ##
##   BINARIZATION PROCESS   ##
## ======================== ##

rm -rf "${base_folder}/5hmc/split" "${base_folder}/5hmc/blanks" "${base_folder}/5hmc/bin_counts" "${base_folder}/5hmc/binarized"
mkdir -p "${base_folder}/5hmc/split" "${base_folder}/5hmc/blanks" "${base_folder}/5hmc/bin_counts" "${base_folder}/5hmc/binarized"

for chromosome in {1..22} X; do
  awk \
    -v chromosome="$chromosome" \
    '$1 == "chr"chromosome' \
    "purified_reads.bed" > \
    "${base_folder}/5hmc/split/purified_chr${chromosome}.bed"
done

## ---------------- ##
##   BIN CREATION   ##
## ---------------- ##

module purge
module load R/4.2.1-foss-2022a

Rscript "$RSCRIPT_DIR/create_blank_bed_files.R" \
  "$chromosome_sizes" \
  "$bin_size" \
  "${base_folder}/5hmc/blanks"

## ---------------- ##
##   INTERSECTION   ##
## ---------------- ##

# We want to only use bins that have methylated sites within them when 
# binarizing. This is so that we can more noticably discern bins with little
# methylation and those with an actual peak in methylation. This is required
# as the majority of the genome is unmethylated, using a global 'baseline'
# signal will not be 'powerful' enough.

module purge
module load BEDTools

for chromosome in {1..22} X; do
  bedtools intersect \
    -wa \
    -c \
    -a "${base_folder}/5hmc/blanks/chromosome${chromosome}.bed" \
    -b "${base_folder}/5hmc/split/purified_chr${chromosome}.bed" > \
    "${base_folder}/5hmc/bin_counts/chromosome${chromosome}.bed"
done

## ------------ ##
##   BINARIZE   ##
## ------------ ##

module purge
module load R/4.2.1-foss-2022a

for chromosome in {1..22} X; do
  dense_file="${base_folder}/5hmc/binarized/dense/${cell_type}_chr${chromosome}_binary.txt"
  echo -e "${cell_type}\tchr${chromosome}" > "$dense_file" 
  echo "5hmC_dense" >> "$dense_file"

  sparse_file="${base_folder}/5hmc/binarized/sparse/${cell_type}_chr${chromosome}_binary.txt"
  echo -e "${cell_type}\tchr${chromosome}" > "$sparse_file" 
  echo "5hmC_sparse" >> "$sparse_file"

  Rscript "$RSCRIPT_DIR/binarize.R" \
    "${base_folder}/5hmc/bin_counts/chromosome${chromosome}.bed" \
    "$dense_file" \
    "$sparse_file"

  gzip "$dense_file" "$sparse_file"
done

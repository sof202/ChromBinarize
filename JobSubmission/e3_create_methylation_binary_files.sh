#!/bin/bash
#SBATCH --export=ALL
#SBATCH -p mrcq 
#SBATCH --time=03:00:00 
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
#SBATCH --mem=4G 
#SBATCH --mail-type=END 
#SBATCH --output=methylation%j.log
#SBATCH --error=methylation%j.err
#SBATCH --job-name=methylation

SCRIPT_PATH=$(scontrol show job "$SLURM_JOBID" | \
  awk '/Command=/{print $1}' | \
  cut -d= -f1)
SCRIPT_DIR=$(realpath "$(dirname "$SCRIPT_PATH")")

ROOT_DIR="${SCRIPT_DIR}/.."
RSCRIPT_DIR="${ROOT_DIR}/Rscripts"

usage() {
cat <<EOF
================================================================================
e2_create_methylation_binary_files.sh
================================================================================
Purpose: Create binary files for dense and sparse regions of methylation
from whole genome bisulphite sequencing bed files. Bed files are expected to be 
in the format:
  chr \t start \t end \t methylationmethylated_reads \t total_reads  
Author: Sam Fletcher
Contact: s.o.fletcher@exeter.ac.uk
Dependencies: R, awk, bedtools
================================================================================
EOF
    exit 0
}

if [ "$#" -eq 0 ]; then usage; fi 

## ============================ ##
##   EXTRACT METHYLATED SITES   ##
## ============================ ##

# config will source all of the variables seen below
config_file_location=$1
source "${config_file_location}" || exit 1

mkdir -p "${LOG_DIR}/"
mv "${SLURM_SUBMIT_DIR}/methylation${SLURM_JOB_ID}.log" \
  "${LOG_DIR}/methylation${SLURM_JOB_ID}.log"
mv "${SLURM_SUBMIT_DIR}/methylation${SLURM_JOB_ID}.err" \
  "${LOG_DIR}/methylation${SLURM_JOB_ID}.err"

## ------------------------------------ ##
##   REMOVE HYDROXYMETHYLATION SIGNAL   ##
## ------------------------------------ ##

# This requires oxBS and WGBS files to work. WGBS captures 5mC AND 5hmC signal.
# We only want to capture the 5mC signal, so we need to remove the 5hmC signal
# using oxBS data.

module purge
module load BEDTools

bedtools intersect -wo \
  -a "${WGBS_bed_file_location}" \
  -b "${oxBS_bed_file_location}" |
  awk \
    -v read_threshold="${reference_read_depth_threshold_h}" \
    'function convert_to_percent(reads, total_reads) {
       return (int(reads / total_reads * 10000) / 100)
     }
     {OFS="\t"} 
     $10 >= read_threshold {print $1,$2,$3,$5,convert_to_percent($4,$5),convert_to_percent($9,$10)}' > \
  "${base_folder}/5mc/combined.bed"

awk '
    function convert_to_reads(percentage, total_reads) {
      return (int(percentage * total_reads / 100))
    }
    function ReLU_distance(i,j) {
      return ((j - i) > 0 ? (j - i) : 0)
    }
    {OFS="\t"}
    {print $1,$2,$3,convert_to_reads(ReLU_distance($5,$10), $4), $4}
    ' "${base_folder}/5mc/combined.bed" > \
      "${base_folder}/5mc/WGBS_5hmc_removed.bed"


## ------------------------------- ##
##   EXTRACT CONFIDENT POSITIONS   ##
## --------------------------------##

awk -v percent_threshold="${reference_percentage_threshold_h}" \
  -v read_threshold="${reference_read_depth_threshold_h}" \
  '$5 >= read_threshold && (int($4/$5 * 10000)/100) >= percent_threshold {print $5","$7}' \
  "${base_folder}/5mc/WGBS_5hmc_removed.bed" > "${base_folder}/5mc/methylated.csv"

awk -v percent_threshold=$((100 - ${reference_percentage_threshold_h})) \
  -v read_threshold="${reference_read_depth_threshold_h}" \
  '$5 >= read_threshold && (int($4/$5 * 10000)/100) <= percent_threshold {print $5","$7}' \
  "${base_folder}/5mc/WGBS_5hmc_removed.bed" > "${base_folder}/5mc/unmethylated.csv"

awk -v read_threshold="${minimum_read_depth}" \
  '{OFS="\t"} $5 >= read_threshold {print $1,$2,$3,"m",$5,"+",int($4/$5 * 10000) / 100}' \
  "${base_folder}/5mc/WGBS_5hmc_removed.bed" > "${base_folder}/5mc/filtered_reads.bed"

## ------------------------- ##
##   RUN BINOMIAL ANALYSIS   ##
## ------------------------- ##

module purge
module load R/4.2.1-foss-2022a

Rscript "${RSCRIPT_DIR}/binom.R" "${base_folder}/5mc"

module purge

## -------------------------==== ##
##   REMOVE UNMETHYLATED SITES   ##
## -------------------------==== ##

awk -v threshold="${binomial_threshold}" \
  '$9 < threshold' \
  "${base_folder}/5mc/processed_reads.bed" > \
  "${base_folder}/5mc/purified_reads.bed"

## ======================== ##
##   BINARIZATION PROCESS   ##
## ======================== ##

rm -rf "${base_folder}/5mc/split" "${base_folder}/5mc/blanks" "${base_folder}/5mc/bin_counts" "${base_folder}/5mc/binarized"
mkdir -p "${base_folder}/5mc/split" "${base_folder}/5mc/blanks" "${base_folder}/5mc/bin_counts" "${base_folder}/5mc/binarized"

for chromosome in {1..22} X; do
  awk \
    -v chromosome="$chromosome" \
    '$1 == "chr"chromosome' \
    "purified_reads.bed" > \
    "${base_folder}/5mc/split/purified_chr${chromosome}.bed"
done

## ---------------- ##
##   BIN CREATION   ##
## ---------------- ##

module purge
module load R/4.2.1-foss-2022a

Rscript "$RSCRIPT_DIR/create_blank_bed_files.R" \
  "$chromosome_sizes" \
  "$bin_size" \
  "${base_folder}/5mc/blanks"

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
    -a "${base_folder}/5mc/blanks/chromosome${chromosome}.bed" \
    -b "${base_folder}/5mc/split/purified_chr${chromosome}.bed" > \
    "${base_folder}/5mc/bin_counts/chromosome${chromosome}.bed"
done

## ------------ ##
##   BINARIZE   ##
## ------------ ##

module purge
module load R/4.2.1-foss-2022a

for chromosome in {1..22} X; do
  dense_file="${base_folder}/5mc/binarized/dense/${cell_type}_chr${chromosome}_binary.txt"
  echo -e "${cell_type}\tchr${chromosome}" > "$dense_file" 
  echo "5mc_dense" >> "$dense_file"

  sparse_file="${base_folder}/5mc/binarized/sparse/${cell_type}_chr${chromosome}_binary.txt"
  echo -e "${cell_type}\tchr${chromosome}" > "$sparse_file" 
  echo "5mc_sparse" >> "$sparse_file"

  Rscript "$RSCRIPT_DIR/binarize.R" \
    "${base_folder}/5mc/bin_counts/chromosome${chromosome}.bed" \
    "$dense_file" \
    "$sparse_file"

  gzip "$dense_file" "$sparse_file"
done

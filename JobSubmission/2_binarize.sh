#!/bin/bash
#SBATCH --export=ALL
#SBATCH -p mrcq 
#SBATCH --time=03:00:00 
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
#SBATCH --mem=4G 
#SBATCH --mail-type=END 
#SBATCH --output=binarize%j.log
#SBATCH --error=binarize%j.err
#SBATCH --job-name=Binarize


SCRIPT_PATH=$(scontrol show job "$SLURM_JOBID" | \
  awk '/Command=/{print $1}' | \
  cut -d= -f1)
SCRIPT_DIR=$(realpath "$(dirname "$SCRIPT_PATH")")
ROOT_DIR="${SCRIPT_DIR}/.."
RSCRIPT_DIR="${ROOT_DIR}/Rscripts"

mkdir -p "${ROOT_DIR}/logs/"
mv "${SLURM_SUBMIT_DIR}/binarize${SLURM_JOB_ID}.log" \
  "${ROOT_DIR}/logs/binarize${SLURM_JOB_ID}.log"
mv "${SLURM_SUBMIT_DIR}/binarize${SLURM_JOB_ID}.err" \
  "${ROOT_DIR}/logs/binarize${SLURM_JOB_ID}.err"

usage() {
cat <<EOF
================================================================================
2_binarize.sh
================================================================================
Purpose: Turns processed data into binarized data
Author: Sam Fletcher
Contact: s.o.fletcher@exeter.ac.uk
Dependencies: R, bedtools
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

if [ "${mark}" == "m" ]; then
  cd "${base_folder}/5mc" || exit 1
else
  cd "${base_folder}/5hmc" || exit 1
fi

## ======================== ##
##   SPLITING CHROMOSOMES   ##
## ======================== ##

rm -rf split blanks bin_counts binarized
mkdir -p split blanks bin_counts binarized

chromosomes=$(seq 1 22)
chromosomes=$(echo -e "$chromosomes\nX")

for chromosome in $(echo "$chromosomes"); do
  awk \
    -v chromosome="$chromosome" \
    '$1 == "chr"chromosome' \
    "purified_reads.bed" > \
    "split/purified_chr${chromosome}.bed"
done

## ======== ##
##   BINS   ##
## ======== ##

module purge
module load R/4.2.1-foss-2022a

Rscript "$RSCRIPT_DIR/create_blank_bed_files.R" \
  "$chromosome_sizes" \
  "$bin_size" \
  "$(pwd)/blanks"

## ================ ##
##   INTERSECTION   ##
## ================ ##

# We want to only use bins that have methylated sites within them when 
# binarizing. This is so that we can more noticably discern bins with little
# methylation and those with an actual peak in methylation. This is required
# as the majority of the genome is unmethylated, using a global 'baseline'
# signal will not be 'powerful' enough.

module purge
module load BEDTools

for chromosome in $(echo "$chromosomes"); do
  bedtools intersect \
    -wa \
    -c \
    -a "blanks/chromosome${chromosome}.bed" \
    -b "split/purified_chr${chromosome}.bed" > \
    "bin_counts/chromosome${chromosome}.bed"
done

## ============ ##
##   BINARIZE   ##
## ============ ##

module purge
module load R/4.2.1-foss-2022a

if [ "${mark}" == "m" ]; then
  mark="5-methylcytosine"
else
  mark="5-hydroxymethylcytosine"
fi


for chromosome in $(echo "$chromosomes"); do
  dense_file="binarized/dense/${cell_type}_chr${chromosome}_binary.txt"
  echo -e "${cell_type}\tchr${chromosome}" > "$dense_file" 
  echo "$mark" >> "$dense_file"

  sparse_file="binarized/sparse/${cell_type}_chr${chromosome}_binary.txt"
  echo -e "${cell_type}\tchr${chromosome}" > "$sparse_file" 
  echo "$mark" >> "$sparse_file"

  Rscript "$RSCRIPT_DIR/binarize.R" \
    "bin_counts/chromosome${chromosome}.bed" \
    "$dense_file" \
    "$sparse_file"

  gzip "$dense_file" "$sparse_file"
done

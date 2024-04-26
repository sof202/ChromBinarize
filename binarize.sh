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
              awk -F= '/Command=/{print $2}')
SCRIPT_DIR=$(dirname "$SCRIPT_PATH")
cd "$SCRIPT_DIR" || exit 1
cd ..

mkdir "${SCRIPT_DIR}/logs/"
mv "${SLURM_SUBMIT_DIR}/binarize${SLURM_JOB_ID}.log" \
  "${SCRIPT_DIR}/logs/binarize${SLURM_JOB_ID}.log"
mv "${SLURM_SUBMIT_DIR}/binarize${SLURM_JOB_ID}.err" \
  "${SCRIPT_DIR}/logs/binarize${SLURM_JOB_ID}.err"

mark=$1             #5mc or 5hmc
chromosome_sizes=$2 #file of chromosome sizes
bin_size=$3         #bin size expected to be used with ChromHMM
cell_type=$4        #celltype to be used in ChromHMM

cd "$mark" || exit 1

## ======================== ##
##   SPLITING CHROMOSOMES   ##
## ======================== ##

mkdir split blanks bin_counts binarized

chromosomes=$(seq 1 22)
chromosomes=$(echo -e "$chromosomes\nX")

for chromosome in $(echo "$chromosomes"); do
  awk \
    -v chromosome="$chromosome" \
    '$1 == "chr"chromsome' \
    "purified_reads.bed" > \
    "split/purified_chr${chromosome}.bed"
done

## ======== ##
##   BINS   ##
## ======== ##

module purge
module load R

Rscript "$SCRIPT_DIR/create_blank_bed_files.R" \
  "$chromosome_sizes" \
  "$bin_size" \
  "blanks"

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
    -a "split/purified_chr${chromosome}.bed" \
    -b "blanks/chromosome${chromosome}.bed" > \
    "bin_counts/chromosome${chromosome}.bed"
done

## ============ ##
##   BINARIZE   ##
## ============ ##

module purge
module load R

for chromosome in $(echo "$chromosomes"); do
  file="binarized/${cell_type}_chr${chromosome}_binary.txt"
  echo -e "${cell_type}\tchr${chromosome}" > "$file" 
  echo "$mark"

  Rscript binarize.R \
    "bin_counts/chromosome${chromosome}.bed" \
    "binarized/${mark}.binarized.bed"

  gzip "$file"
done

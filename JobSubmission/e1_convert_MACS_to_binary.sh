#!/bin/bash
#SBATCH --export=ALL
#SBATCH -p mrcq 
#SBATCH --time=03:00:00 
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
#SBATCH --mem=4G 
#SBATCH --mail-type=END 
#SBATCH --output=convert%j.log
#SBATCH --error=convert%j.err
#SBATCH --job-name=convert

SCRIPT_PATH=$(scontrol show job "$SLURM_JOBID" | \
  awk '/Command=/{print $1}' | \
  cut -d= -f1)
SCRIPT_DIR=$(realpath "$(dirname "$SCRIPT_PATH")")

ROOT_DIR="${SCRIPT_DIR}/.."
RSCRIPT_DIR="${ROOT_DIR}/Rscripts"

mkdir -p "${LOG_DIR}/"
mv "${SLURM_SUBMIT_DIR}/convert${SLURM_JOB_ID}.log" \
  "${LOG_DIR}/convert${SLURM_JOB_ID}.log"
mv "${SLURM_SUBMIT_DIR}/convert${SLURM_JOB_ID}.err" \
  "${LOG_DIR}/convert${SLURM_JOB_ID}.err"

usage() {
cat <<EOF
================================================================================
e1_convert_MACS_to_binary.sh
================================================================================
Purpose: Converts narrow peak or broad peaks to binary format for chromHMM
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
module load R/4.2.1-foss-2022a

mkdir -p "${base_folder}/4_BinarizedFiles/${epigenetic_mark}"

Rscript ${RSCRIPT_DIR}/create_blank_bed_files.R \
  "${chromosome_sizes}" \
  "${bin_size}" \
  "${base_folder}/4_BinarizedFiles/${epigenetic_mark}"

module purge
module load BEDTools

for chr in {1..22} X; do
  output_binary_file="${base_folder}/4_BinarizedFiles/${epigenetic_mark}/${cell_type}_chr${chr}_binary.txt.gz"
  echo -e "${cell_type}\tchr${chr}" > "${output_binary_file}"
  echo "${epigenetic_mark}" >> "${output_binary_file}"

  bedtools intersect \
    -wa \
    -c \
    -a "${base_folder}/4_BinarizedFiles/${epigenetic_mark}/chromosome${chr}.bed" \
    -b "${input_MACS_file}" | \
    awk '{OFS="\t"} {print ($4 > 0 ? 1 : 0)}' | \
    gzip >> \
    "${output_binary_file}"

  rm "${base_folder}/4_BinarizedFiles/${epigenetic_mark}/chromosome${chr}.bed"
done





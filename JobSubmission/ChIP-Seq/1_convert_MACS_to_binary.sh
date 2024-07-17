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

usage() {
cat <<EOF
================================================================================
1_convert_MACS_to_binary.sh
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

config_file_location=$1
source "${config_file_location}" || { echo "could not find config file at:
${config_file_location}"; exit 1; }

source "${REPO_DIR}/parameters.txt" || exit 1

for file in "${FUNCTIONS_DIR}"/*; do source "$file" || exit 1; done

move_log_files convert

module purge
module load R/4.2.1-foss-2022a

## ================ ##
##   BINARIZATION   ##
## ================ ##

rm -rf "${BINARY_DIR}/${epigenetic_mark_name:?}"
mkdir -p "${BINARY_DIR}/${epigenetic_mark_name}"

Rscript ${RSCRIPT_DIR}/create_blank_bed_files.R \
  "${chromosome_sizes}" \
  "${bin_size}" \
  "${BINARY_DIR}/${epigenetic_mark_name}"

module purge
module load BEDTools

for chr in {1..22} X; do
  output_binary_file="${BINARY_DIR}/${epigenetic_mark_name}/${cell_type}_chr${chr}_binary.txt"
  echo -e "${cell_type}\tchr${chr}" > "${output_binary_file}"
  echo "${epigenetic_mark_name}" >> "${output_binary_file}"

  bedtools intersect \
    -wa \
    -c \
    -a "${BINARY_DIR}/${epigenetic_mark_name}/chromosome${chr}.bed" \
    -b "${input_MACS_file}" | \
    awk '{OFS="\t"} {print ($4 > 0 ? 1 : 0)}' >> \
    "${output_binary_file}"

  gzip "${output_binary_file}"

  rm "${BINARY_DIR}/${epigenetic_mark_name}/chromosome${chr}.bed"
done


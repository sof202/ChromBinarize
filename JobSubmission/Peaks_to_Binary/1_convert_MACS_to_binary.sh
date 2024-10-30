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
$(basename "$0")
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

for file in "${FUNCTIONS_DIR}"/*; do source "$file" || exit 1; done

move_log_files convert

## ================ ##
##   BINARIZATION   ##
## ================ ##

rm -rf "${BINARY_DIR}/${epigenetic_mark_name:?}"
mkdir -p "${BINARY_DIR}/${epigenetic_mark_name}"

conda activate ChromBinarize-R
Rscript ${RSCRIPT_DIR}/create_blank_bed_files.R \
  "${REPO_DIR}" \
  "${chromosome_sizes}" \
  "${BIN_SIZE}" \
  "${BINARY_DIR}/${epigenetic_mark_name}"
conda deactivate

logs "${DEBUG_MODE:0}" \
"Generating ChromHMM binary files..."

for chromosome in {1..22} X; do
  output_binary_file="${BINARY_DIR}/${epigenetic_mark_name}/${cell_type}_chr${chromosome}_binary.txt"
  echo -e "${cell_type}\tchr${chromosome}" > "${output_binary_file}"
  echo "${epigenetic_mark_name}" >> "${output_binary_file}"

  conda activate ChromBinarize-bedtools
  bedtools intersect \
    -wa \
    -c \
    -a "${BINARY_DIR}/${epigenetic_mark_name}/chromosome${chromosome}.bed" \
    -b "${input_MACS_file}" | \
    awk '{OFS="\t"} {print ($4 > 0 ? 1 : 0)}' >> \
    "${output_binary_file}"
  conda deactivate

  number_of_signatures=$(awk 'NR>2 && $1>0' "${output_binary_file}" | wc -l)

logs 1 \
"chromosome ${chromosome} has:
${number_of_signatures} signatures."

    if [[ "${number_of_signatures}" -eq 0 ]]; then
errors "${chromosome}'s binary file has no true/1 entries. 
Please check to see if your input MACS file is empty."
    fi

  gzip "${output_binary_file}"

  rm "${BINARY_DIR}/${epigenetic_mark_name}/chromosome${chromosome}.bed"
done


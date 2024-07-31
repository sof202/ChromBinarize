#!/bin/bash
#SBATCH --export=ALL
#SBATCH -p mrcq 
#SBATCH --time=03:00:00 
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
#SBATCH --mem=4G 
#SBATCH --mail-type=END 
#SBATCH --output=binchange%j.log
#SBATCH --error=binchange%j.err
#SBATCH --job-name=binchange

usage() {
cat <<EOF
================================================================================
e_change_bin_size.sh
================================================================================
Purpose: Converts a binary file from one bin size to another.
WARNING: For best results, the new bin size should be smaller than the original
bin size.
Inputs:
  \$1 -> Config file location
  \$2 -> Directory with binary files
  \$3 -> Directory for new binary files
  \$4 -> Original bin size
  \$5 -> New bin size
Optional: -r -> remove the original binary files
Author: Sam Fletcher
Contact: s.o.fletcher@exeter.ac.uk
Dependencies: R, bedtools
================================================================================
EOF
    exit 0
}

while getopts r OPT; do
    case "$OPT" in
        r )       remove_original_files="TRUE" ;;
        * )       usage ;;
    esac
done
shift $((OPTIND-1))

if [ "$#" -eq 0 ]; then usage; fi 

config_file_location=$1
source "${config_file_location}" || { echo "could not find config file at:
${config_file_location}"; exit 1; } 

for file in "${FUNCTIONS_DIR}"/*; do source "$file" || exit 1; done

move_log_files binchange

## ============= ##
##   ARGUMENTS   ##
## ============= ##

old_binary_directory=$2
new_binary_directory=$3
original_bin_size=$4
new_bin_size=$5

if [[ "${old_binary_directory}" == "${new_binary_directory}" ]]; then
errors \
"Setting the new binary directory to the same position as the old directory
will result in errors. Exiting..."
  exit 1
fi

if [[ "${original_bin_size}" -lt "${new_bin_size}" ]]; then
errors \
"WARNING: Setting the new bin size to be larger than the original bin size
is ill advised. See README.md."
fi

if [[ "$#" -ne 5 ]]; then usage; fi

## ======== ##
##   MAIN   ##
## ======== ##

rm -rf "${new_binary_directory}/blanks" \
  "${new_binary_directory}/old_binary_files" \
  "${old_binary_directory}/blanks"
mkdir -p "${new_binary_directory}/blanks" \
  "${new_binary_directory}/old_binary_files" \
  "${old_binary_directory}/blanks"

## ------------------------------------- ##
##   GENERATE BED FILE FROM OLD BINARY   ##
## ------------------------------------- ##

binarization_createBlankBins "${new_binary_directory}" "${new_bin_size}"
binarization_createBlankBins "${old_binary_directory}" "${original_bin_size}"

for chromosome in {1..22}X; do
  binary_file=$(find "${old_binary_directory}" -name "*chr${chromosome}*binary.txt.gz")

logs "${DEBUG_MODE:0}" \
"Generating bed file from old binary file: ${binary_file}..."

  gzip -dc "${binary_file}" | awk 'NR>2' > \
    "${new_binary_directory}/old_binary_files/chr${chromosome}.bed"

  mark_name=$(head -2 "${new_binary_directory}/old_binary_files/chr${chromosome}.bed" | tail -1)

  paste \
    "${old_binary_directory}/blanks/chromosome${chromosome}.bed"
    "${new_binary_directory}/old_binary_files/chr${chromosome}.bed" > \
      "${new_binary_directory}/old_binary_files/combined_chr${chromosome}.bed"

  gzip "${new_binary_directory}/old_binary_files/combined_chr${chromosome}.bed"
done

## --------------------------- ##
##   CREATE NEW BINARY FILES   ##
## --------------------------- ##

for chromosome in {1..22}X; do
logs "${DEBUG_MODE:0}" \
"Generating new binary file for chromosome ${chromosome}..."

  new_binary_file="${new_binary_directory}/${cell_type}_chr${chromosome}_binary.txt"
  echo -e "${cell_type}\tchr${chromosome}" > "${new_binary_file}" 
  echo "${mark_name}" >> "${new_binary_file}"

  bedtools intersect -wo \
    -a "${new_binary_directory}/blanks/chromosome${chromosome}.bed" \
    -b "${new_binary_directory}/old_binary_files/combined_chr${chromosome}.bed.gz" | \
    awk '{OFS="\t"} {print $7}' >> \
      "${new_binary_file}"
done 

## ============ ##
##   CLEAN UP   ##
## ============ ##

rm -rf "${new_binary_directory}/blanks" \
  "${new_binary_directory}/old_binary_files" \
  "${old_binary_directory}/blanks"

if [[ ${remove_original_files} == "TRUE" ]]; then
  rm -rf "${old_binary_directory}"
fi

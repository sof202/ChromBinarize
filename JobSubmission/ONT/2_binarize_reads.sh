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

usage() {
cat <<EOF
================================================================================
$(basename "$0")
================================================================================
Purpose: Turns processed ONT data into dense and sparsed ChromHMM compliable
binarized data.
Author: Sam Fletcher
Contact: s.o.fletcher@exeter.ac.uk
Dependencies: R, bedtools
================================================================================
EOF
    exit 0
}

if [ "$#" -eq 0 ]; then usage; fi 

config_file_location=$1
source "${config_file_location}" || { echo "could not find config file at:
${config_file_location}"; exit 1; } 

for file in "${FUNCTIONS_DIR}"/*; do source "$file" || exit 1; done

move_log_files binarize

## ================ ##
##   BINARIZATION   ##
## ================ ##

if [ "${mark}" == "m" ]; then
  processing_directory="${BASE_DIR}/ONT_5mc"
  mark_name="ONT_5mC"
elif [ "${mark}" == "h" ]; then
  processing_directory="${BASE_DIR}/ONT_5hmc"
  mark_name="ONT_5hmC"
else
  >&2 echo "config file needs 'm' or 'h' in the 'mark' field"
  exit 1
fi


binarization_createDirectories \
  "${processing_directory}"
binarization_splitIntoChromosomes \
  "${processing_directory}" \
  "purified_reads.bed"
binarization_createBlankBins \
  "${processing_directory}" \
  "${BIN_SIZE}"
binarization_countSignalIntersectionWithBins \
  "${processing_directory}"
binarization_createChromhmmBinaryFiles \
  "${processing_directory}" \
  "${BINARY_DIR}/${mark_name}" \
  "${mark_name}"

if [[ "${DEBUG_MODE:=0}" -eq 0 ]]; then
  rm -rf "${processing_directory}"
fi

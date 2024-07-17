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
2_binarize_reads.sh
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

config_file_location=$1
source "${config_file_location}" || { usage; exit 1; } 

source "${REPO_DIR}/parameters.txt" || exit 1

source "${FUNCTIONS_DIR}/move_log_files.sh" || exit 1
move_log_files binarize

## ======== ##
##   MAIN   ##
## ======== ##

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

source "${FUNCTIONS_DIR}/binarization.sh" || exit 1

binarization_createDirectories \
  "${processing_directory}"
binarization_splitIntoChromosomes \
  "${processing_directory}" \
  "purified_reads.bed"
binarization_createBlankBins \
  "${processing_directory}"
binarization_countSignalIntersectionWithBins \
  "${processing_directory}"
binarization_createChromhmmBinaryFiles \
  "${processing_directory}" \
  "${BINARY_DIR}/${mark_name}" \
  "${mark_name}"

if [[ ! "${debug_mode:='false'}" == "true" ]]; then
  rm -rf "${processing_directory}"
fi

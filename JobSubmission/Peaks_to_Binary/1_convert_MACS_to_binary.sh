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

all_marks=$(\
  find "${macs_directory}" -type f -name "*.bed" | \
  xargs -n 1 basename | \
  cut -d. -f1\
)
for mark in ${all_marks}; do
  binarization_convertToChromHMM "$mark"
done

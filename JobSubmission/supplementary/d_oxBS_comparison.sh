#!/bin/bash
#SBATCH --export=ALL
#SBATCH -p mrcq 
#SBATCH --time=03:00:00 
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
#SBATCH --mem=4G 
#SBATCH --mail-type=END 
#SBATCH --output=compare%j.log
#SBATCH --error=compare%j.err
#SBATCH --job-name=compare

usage() {
cat <<EOF
================================================================================
d_oxBS_comparison.sh
================================================================================
Purpose: Creates several plots comparing a two bed files. One for ONT, the 
other for oxBS data.
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
source "${config_file_location}" || usage

source "${REPO_DIR}/parameters.txt" || exit 1

source "${FUNCTIONS_DIR}/move_log_files.sh" || exit 1
move_log_files compare

## =========================== ##
##   INTERSECT ONT WITH OXBS   ##
## =========================== ##

oxBS_folder=$(dirname "${oxBS_bed_file_location}")

source "${FUNCTIONS_DIR}/intersect.sh" || exit 1

intersect_intersectBSWithONT \
  "${ONT_bed_file_location}" \
  "${oxBS_bed_file_location}" \
  "${oxBS_folder}/ONT_oxBS_intersect.bed"  

## ============== ##
##   COMPARISON   ##
## ============== ##

module purge
module load R/4.2.1-foss-2022a

mkdir -p "${BASE_DIR}/plots"

Rscript "${RSCRIPT_DIR}/oxBS_comparison.R" \
  "${oxBS_folder}/ONT_oxBS_intersect.bed" \
  "${BASE_DIR}/plots"

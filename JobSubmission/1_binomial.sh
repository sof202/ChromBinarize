#!/bin/bash
#SBATCH --export=ALL
#SBATCH -p mrcq 
#SBATCH --time=03:00:00 
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
#SBATCH --mem=4G 
#SBATCH --mail-type=END 
#SBATCH --output=binom%j.log
#SBATCH --error=binom%j.err
#SBATCH --job-name=Binom

SCRIPT_PATH=$(scontrol show job "$SLURM_JOBID" | \
  awk '/Command=/{print $1}' | \
  cut -d= -f1)
SCRIPT_DIR=$(realpath "$(dirname "$SCRIPT_PATH")")
cd "$SCRIPT_DIR" || exit 1
cd ..

mkdir -p "${SCRIPT_DIR}/logs/"
mv "${SLURM_SUBMIT_DIR}/binom${SLURM_JOB_ID}.log" \
  "${SCRIPT_DIR}/logs/binom${SLURM_JOB_ID}.log"
mv "${SLURM_SUBMIT_DIR}/binom${SLURM_JOB_ID}.err" \
  "${SCRIPT_DIR}/logs/binom${SLURM_JOB_ID}.err"


bed_file_location=$1
base_folder=$(pwd)

rm -rf 5mc
rm -rf 5hmc
mkdir -p 5mc 5hmc

awk '$4 == "m" && $5 >= 500 && $7 >= 95 {print $5","$7}' "${bed_file_location}" > "${base_folder}/5mc/methylated.csv"
awk '$4 == "m" && $5 >= 500 && $7 <= 5 {print $5","$7}' "${bed_file_location}" > "${base_folder}/5mc/unmethylated.csv"
awk '$4 == "m" && $5 >= 30' "${bed_file_location}" > "${base_folder}/5mc/filtered_reads.bed"

awk '$4 == "h" && $5 >= 50 && $7 >= 95 {print $5","$7}' "${bed_file_location}" > "${base_folder}/5hmc/methylated.csv"
awk '$4 == "h" && $5 >= 50 && $7 <= 5 {print $5","$7}' "${bed_file_location}" > "${base_folder}/5hmc/unmethylated.csv"
awk '$4 == "h" && $5 >= 30' "${bed_file_location}" > "${base_folder}/5hmc/filtered_reads.bed"

module purge
module load R/4.2.1-foss-2022a

Rscript "$SCRIPT_DIR/binom.R" "${base_folder}/5mc"
Rscript "$SCRIPT_DIR/binom.R" "${base_folder}/5hmc"

module purge

awk '$8 > 1.7e-8 && $9 < 1.7e-8' "${base_folder}/5mc/processed_reads.bed" > "${base_folder}/5mc/purified_reads.bed"
awk '$8 > 1.7e-8 && $9 < 1.7e-8' "${base_folder}/5hmc/processed_reads.bed" > "${base_folder}/5hmc/purified_reads.bed"

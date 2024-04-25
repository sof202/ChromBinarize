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

SCRIPT_PATH=$(scontrol show job "$SLURM_JOBID" | awk -F= '/Command=/{print $2}')
SCRIPT_DIR=$(dirname "$SCRIPT_PATH")
cd "$SCRIPT_DIR" || exit 1
cd ..

mv "${SLURM_SUBMIT_DIR}/binom${SLURM_JOB_ID}.log" \
  "${SCRIPT_DIR}/logs/binom${SLURM_JOB_ID}.log"
mv "${SLURM_SUBMIT_DIR}/binom${SLURM_JOB_ID}.err" \
  "${SCRIPT_DIR}/logs/binom${SLURM_JOB_ID}.err"

rm -rf 5mc
rm -rf 5hmc
mkdir 5mc 5hmc

awk '$4 == "m" && $5 >= 500 && $7 >= 95 {print $5","$7}' NeuN_3_sample_sam.bed > 5mc/methylated.csv
awk '$4 == "m" && $5 >= 500 && $7 <= 5 {print $5","$7}' NeuN_3_sample_sam.bed > 5mc/unmethylated.csv
awk '$4 == "m" && $5 >= 30' NeuN_3_sample_sam.bed > 5mc/filtered_reads.bed

awk '$4 == "h" && $5 >= 50 && $7 >= 95 {print $5","$7}' NeuN_3_sample_sam.bed > 5hmc/methylated.csv
awk '$4 == "h" && $5 >= 50 && $7 <= 5 {print $5","$7}' NeuN_3_sample_sam.bed > 5hmc/unmethylated.csv
awk '$4 == "h" && $5 >= 30' NeuN_3_sample_sam.bed > 5hmc/filtered_reads.bed

module purge
module load R

Rscript "scripts/binom.R" 5mc
Rscript "scripts/binom.R" 5hmc

module purge

awk '$8 > 1.7e-8 && $9 < 1.7e-8' 5mc/processed_reads.bed > 5mc/Subsampled.100.5mc.bed
awk '$8 > 1.7e-8 && $9 < 1.7e-8' 5hmc/processed_reads.bed > 5hmc/Subsampled.100.5hmc.bed

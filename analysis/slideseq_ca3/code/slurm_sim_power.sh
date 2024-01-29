#!/bin/bash
#SBATCH --job-name=cs_sim_power
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=philnicol740@gmail.com
#SBATCH --mem=32gb
#SBATCH --output=../logs/sim_power.log

module load R/4.3.2b
module load gcc/9.2.0

Rscript sim_power.R

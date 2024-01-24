#!/bin/bash
#SBATCH --job-name=cs_sim_power
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=philnicol740@gmail.com
#SBATCH --mem=32gb
#SBATCH --output=../logs/sim_power.log

module load R/4.1
Rscript run.R

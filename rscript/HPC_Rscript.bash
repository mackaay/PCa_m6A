#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --mem=50gb
#SBATCH --ntasks-per-node=2

module load R/3.6.3


cd ${SLURM_SUBMIT_DIR}
Rscript rscript_sbatch.R
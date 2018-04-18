#!/bin/bash
#SBATCH -J dada2
#SBATCH --mem=128000
#SBATCH --output=dada2.out
#SBATCH --error=dada2.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=justin.silverman@duke.edu

module load gcc/5.3.0-fasrc01
module load R/3.4.0-gcb01
module load lapack/3.6.1-gcb01

R CMD BATCH --no-save --no-restore dada2.R dada2.Rout

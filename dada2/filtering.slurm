#!/bin/bash
#SBATCH -J filtering
#SBATCH --mem=64000
#SBATCH --output=filtering.out
#SBATCH --error=filtering.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=justin.silverman@duke.edu

module load gcc/5.3.0-fasrc01
module load R/3.4.0-gcb01
module load lapack/3.6.1-gcb01

R CMD BATCH --no-save --no-restore filtering.R filtering.Rout


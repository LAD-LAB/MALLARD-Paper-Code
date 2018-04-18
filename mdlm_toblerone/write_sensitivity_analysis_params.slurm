#!/bin/bash
#SBATCH -J write_params
#SBATCH --mem=2000
#SBATCH --output=write_params.out
#SBATCH --error=write_params.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=justin.silverman@duke.edu

module load gcc/5.3.0-fasrc01
module load R/3.4.0-gcb01
module load lapack/3.6.1-gcb01

R CMD BATCH --no-save --no-restore write_sensitivity_analysis_params.R write_sensitivity_analysis_params.Rout


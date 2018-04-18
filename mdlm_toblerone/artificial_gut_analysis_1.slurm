#!/bin/bash
#SBATCH -J art_gut_analysis
#SBATCH --mem=32000
#SBATCH --output=art_gut_analysis.out
#SBATCH --error=art_gut_analysis.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=justin.silverman@duke.edu
#SBATCH --array=1-15

module load gcc/5.3.0-fasrc01
module load R/3.4.0-gcb01
module load lapack/3.6.1-gcb01

R CMD BATCH --no-save --no-restore artificial_gut_analysis.R artificial_gut_analysis.${SLURM_ARRAY_TASK_ID}.Rout


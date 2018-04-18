# Run Simulation Analysis
simaJID=$(sbatch simulation.slurm)

# Write Sensitivity Analysis Params
senaJID=$(sbatch write_sensitivity_analysis_params.slurm)

# Run MCMC
rstanJID=$(sbatch -n 5 -N 1 --dependency=afterok:${senaJID##* }  artificial_gut_analysis_1.slurm)

# Run Sensitivity Analysis
sbatch -n 1 -N 1 --dependency=afterok:${rstanJID##* } artificial_gut_analysis_2.slurm

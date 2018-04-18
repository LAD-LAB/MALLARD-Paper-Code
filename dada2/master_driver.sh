# Driver Script for dada2

# Remove Primers 
rpJID=$(sbatch remove_primers.slurm)

# Sync Barcodes
sbJID=$(sbatch --dependency=afterok:${rpJID##* } sync_barcodes.slurm)

# Demultiplex
dmJID=$(sbatch --dependency=afterok:${sbJID##* } demultiplex.slurm)

# Filtering
filtJID=$(sbatch --dependency=afterok:${dmJID##* } filtering.slurm)

# Dada2
sbatch -n 8 -N 1 --dependency=afterok:${filtJID##* } dada2.slurm

This directory provides the code used in the MALLARD paper by Silverman et al.

--------
OVERVIEW
--------

The primary analysis is split into two parts:
(1) Producing sequence variant table
(2) Running the mdlm_toblerone model (aka the MALLARD model using in our paper;
    we call this the "toblerone" model due to the similarity of the latent
    state space to the container a toblerone chocolate bar comes in)

For both of these steps, the primary script calls are given in files called
master_driver.sh. Most of these scripts have been designed to be run on a cluster using the Slurm
job scheduler. With slight modeification, all scripts can be run on a local machine as well. 
Numerous scripts include a "machine" variable that allows, local evaluation 
with appropriately updated path variables. 

There is also a small directory called creating_better_family_tree that
just writes a tree file manually created to describe phylogenetic relationships
between bacterial families included in analysis (see Methods). 

-------------------------------------------
Notes on downloading data and running dada2
-------------------------------------------

The data has been uploaded to NCBI SRA server and indexed by BioProject ID: PRJNA450809
The data consists of two sequencing runs of paired forward and reverse reads. 

The scripts in the dada2 directory were designed to operate on the multiplexed sequencing files. 
However, to accord with the requirements of the NCBI SRA database we uploaded the demultiplexed
sequencing files that had also been trimmed to remove primers.  One complication is that SRA also 
required that forward and reverse reads files for a given sample be distinguished by file name 
rather than directory location as these scripts use. Therefore the data downloaded from 
the SRA database corresponds to the sequencing files that would be found in the dada2/3_demultiplex/
directory however these files must be sorted into the subdirectories

s1r1_split_samples/
s1r2_split_samples/
s2r1_split_samples/
s2r2_split_samples/

based on whether they come from the forward (r1) or reverse (r2) reads of the first (s1) or second (s2) 
sequencing run. 
The variable seq_run in the BioSample metadata describes whether a sample came from s1 or s2. 
Once these files are sorted the file extensions _F and _R should be removed from all file names. 

Files silva_nr_v123_train_set.fa.gz and silva_species_assignment_v123.fa.gz that were used for
assignment can be downloaded from https://zenodo.org/record/158958#.WteIPdPwZDw 
these should be placed into the subdirectory 0_training/

The original mapping files  are provided as well under 0_mapping/ for ease of 
use  although all the information contained there is avaiable in a slightly cleaner 
format on SRA under the associated BioSample metadata. Similarly, for ease of 
use the artificial gut run notes are provided under runnotes/ although this information 
is also provided on SRA. 

The dada2 code given in filtering.slurm and dada2.slurm should then be run once this initial file 
manipulation is completed. 

The key results from the scripts in the dada2 directory are also provided as both a compressed
R data file (phyloseq.rds) and as tab separated text files: taxtab.nochim.tsv, seqtab.nochim.tsv, 
and refseqs.nochim.tsv. The R data file (phyloseq.rds) and the associated mapping files
is all that is required to ruth the MALLARD scripts. 

---------------------------------
Running MALLARD (mdlm_toblerone/)
---------------------------------

All code for both simulations and longitudinal analysis of artificial gut data can be run through the executible master_driver.sh
in the mdlm_toblerone/ directory. Note, a number of R scripts in this folder contain paths that must be manually specified prior to running. All path variables are located at the top of the associated R scripts. 

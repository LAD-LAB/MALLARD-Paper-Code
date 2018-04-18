#!/bin/bash
#SBATCH -J remove_primers
#SBATCH --mem=16000
#SBATCH --output=remove_primers.out
#SBATCH --error=remove_primers.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=justin.silverman@duke.edu

PATH1=0_data_raw/run1
PATH2=0_data_raw/run2
MAP1=0_mapping/2016.03.16MappingFilePool1_ordered_corrected1.txt
MAP2=0_mapping/2016.03.22MappingFile_Pool2_corrected1.txt

# Removing primers see evernote note "Removing Primers from Illumina Data (Paired-End Fastq)
# Using trimmomatic 
# This step is recommended by dada2 in their FAQ on chimera checking

# For Sequencing Run 1
java -jar trimmomatic-0.36.jar PE $PATH1/SAMPLE1_S1_L001_R1_001.fastq.gz $PATH1/SAMPLE1_S1_L001_R2_001.fastq.gz 1_remove_primers/run1/R1.paired.fastq.gz 1_remove_primers/run1/R1.unpaired.fastq.gz 1_remove_primers/run1/R2.paired.fastq.gz 1_remove_primers/run1/R2.unpaired.fastq.gz ILLUMINACLIP:techseqs.fa:2:30:10

# For Sequencing Run 2
java -jar trimmomatic-0.36.jar PE $PATH2/SAMPLE2_S1_L001_R1_001.fastq.gz $PATH2/SAMPLE2_S1_L001_R2_001.fastq.gz 1_remove_primers/run2/R1.paired.fastq.gz 1_remove_primers/run2/R1.unpaired.fastq.gz 1_remove_primers/run2/R2.paired.fastq.gz 1_remove_primers/run2/R2.unpaired.fastq.gz ILLUMINACLIP:techseqs.fa:2:30:10


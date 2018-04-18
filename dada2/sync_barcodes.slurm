#!/bin/bash
#SBATCH -J sync_barcodes
#SBATCH --mem=16000
#SBATCH --output=sync_barcodes.out
#SBATCH --error=sync_barcodes.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=justin.silverman@duke.edu


BARCODE1=0_data_raw/run1/SAMPLE1_S1_L001_I1_001.fastq.gz 
BARCODE2=0_data_raw/run2/SAMPLE2_S1_L001_I1_001.fastq.gz 

./sync_paired_end_reads.py $BARCODE1 $BARCODE1 1_remove_primers/run1/R1.paired.fastq.gz 2_sync_barcodes/run1/I1.synced.fastq.gz 2_sync_barcodes/run1/R1.paired.fastq.gz 
cp 1_remove_primers/run1/R2.paired.fastq.gz 2_sync_barcodes/run1/R2.paired.fastq.gz

./sync_paired_end_reads.py $BARCODE2 $BARCODE2 1_remove_primers/run2/R1.paired.fastq.gz 2_sync_barcodes/run2/I1.synced.fastq.gz 2_sync_barcodes/run2/R1.paired.fastq.gz 
cp 1_remove_primers/run2/R2.paired.fastq.gz 2_sync_barcodes/run2/R2.paired.fastq.gz 


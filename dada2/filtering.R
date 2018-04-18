# Applying dada2 pipeline to bioreactor time-series
## Following tutorial http://benjjneb.github.io/dada2/tutorial.html
## and here http://benjjneb.github.io/dada2/bigdata.html

setwd('/data/davidlab/users/jds/mdlm_publication/dada2/')

library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(ggplot2); packageVersion("ggplot2")

path <- "/data/davidlab/users/jds/mdlm_publication/dada2/"

# Filtering and Trimming --------------------------------------------------

# Forward and Reverse Filenames
fnFs.s1 <-list.files(paste(path, "3_demultiplex/s1r1_split_samples/", sep=""))
fnRs.s1 <-list.files(paste(path, "3_demultiplex/s1r2_split_samples/", sep=""))
fnFs.s2 <-list.files(paste(path, "3_demultiplex/s2r1_split_samples/", sep=""))
fnRs.s2 <-list.files(paste(path, "3_demultiplex/s2r2_split_samples/", sep=""))
# Sort to ensure filenames are in the same order
fnFs.s1 <- sort(fnFs.s1)
fnRs.s1 <- sort(fnRs.s1)
fnFs.s2 <- sort(fnFs.s2)
fnRs.s2 <- sort(fnRs.s2)

sample.names.1 <- sapply(strsplit(fnFs.s1,".fastq", fixed=TRUE), `[`, 1)
sample.names.2 <- sapply(strsplit(fnFs.s2,".fastq", fixed=TRUE), `[`, 1)

# Fully Specify the path for the fnFs and fnRs
fnFs.s1 <- file.path(path, '3_demultiplex/s1r1_split_samples/', fnFs.s1)
fnRs.s1 <- file.path(path, '3_demultiplex/s1r2_split_samples/', fnRs.s1)
fnFs.s2 <- file.path(path, '3_demultiplex/s2r1_split_samples/', fnFs.s2)
fnRs.s2 <- file.path(path, '3_demultiplex/s2r2_split_samples/', fnRs.s2)

# Examine qulaity profiles of the forward and reverse reads
p <- plotQualityProfile(fnFs.s1[[1]]) 
ggsave('Forward_quality_profile_s1.png', plot=p)
p <- plotQualityProfile(fnRs.s1[[1]]) 
ggsave('Reverse_quality_profile_s1.png', plot=p)
p <- plotQualityProfile(fnFs.s2[[1]]) 
ggsave('Forward_quality_profile_s2.png', plot=p)
p <- plotQualityProfile(fnRs.s2[[1]]) 
ggsave('Reverse_quality_profile_s2.png', plot=p)

# Perform filtering and trimming
filtpath <- file.path(path, "4_filter")
dir.create(filtpath)

# For the first sequencing run
dir.create(file.path(filtpath,'s1'))
filtFs.s1 <- file.path(filtpath, 's1', paste0(sample.names.1,"_F_filt.fastq.gz")) 
filtRs.s1 <- file.path(filtpath, 's1', paste0(sample.names.1,"_R_filt.fastq.gz")) 
for (i in seq_along(fnFs.s1)){
    fastqPairedFilter(c(fnFs.s1[i], fnRs.s1[i]), c(filtFs.s1[i], filtRs.s1[i]), 
                      trimLeft=c(10,10), truncLen=c(180, 140),
                      maxN=0, maxEE=2, truncQ=2,
                      compress=TRUE, verbose=TRUE, 
                      rm.phix=TRUE)
}


# For the second sequencing run # Use the same trimming as before so results are comparable
dir.create(file.path(filtpath,'s2'))
filtFs.s2 <- file.path(filtpath, 's2', paste0(sample.names.2,"_F_filt.fastq.gz")) 
filtRs.s2 <- file.path(filtpath, 's2', paste0(sample.names.2,"_R_filt.fastq.gz")) 
for (i in seq_along(fnFs.s2)){
    fastqPairedFilter(c(fnFs.s2[i], fnRs.s2[i]), c(filtFs.s2[i], filtRs.s2[i]), 
                      trimLeft=c(10,10), truncLen=c(180, 140),
                      maxN=0, maxEE=2, truncQ=2, 
                      compress=TRUE, verbose=TRUE,
                      rm.phix=TRUE)
}

# Applying dada2 pipeline to bioreactor time-series
## Following tutorial http://benjjneb.github.io/dada2/tutorial.html
## and here http://benjjneb.github.io/dada2/tutorial.html

setwd('/data/davidlab/users/jds/mdlm_publication/dada2/')

library(dada2); packageVersion("dada2")
library(ggplot2); packageVersion("ggplot2")
library(dplyr); packageVersion("dplyr")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
set.seed(4)

filtpath <- "/data/davidlab/users/jds/mdlm_publication/dada2/4_filter/"

# Find filenames ----------------------------------------------------------

# Forward and reverse filenames
filts.s1 <- list.files(file.path(filtpath,'s1'), full.names=TRUE)
filts.s2 <- list.files(file.path(filtpath,'s2'), full.names=TRUE)

# Only look at Bioreactor Time-series samples
filts.s1  <- filts.s1[grepl("T..", basename(filts.s1), fixed=TRUE)]
filts.s2  <- filts.s2[grepl("T..", basename(filts.s2), fixed=TRUE)]

# Sort to ensure fileneames are in the same order
filts.s1 <- sort(filts.s1)
sample.names.1 <- sapply(strsplit(basename(filts.s1),"_"), `[`, 1)
names(filts.s1) <- sample.names.1

filts.s2 <- sort(filts.s2)
sample.names.2 <- sapply(strsplit(basename(filts.s2),"_"), `[`, 1)
names(filts.s2) <- sample.names.2

# Separate forward and reverse samples
filtFs.s1 <- filts.s1[grepl("_F_filt",filts.s1)]
filtRs.s1 <- filts.s1[grepl("_R_filt",filts.s1)]
filtFs.s2 <- filts.s2[grepl("_F_filt",filts.s2)]
filtRs.s2 <- filts.s2[grepl("_R_filt",filts.s2)]

sample.names.1 <- sapply(strsplit(basename(filtFs.s1), "_"), `[`, 1)
sample.names.2 <- sapply(strsplit(basename(filtFs.s2),"_"), `[`, 1)

# Dereplication -----------------------------------------------------------

# Learn Error Rates
## aim to learn from about 1M total reads - so just need subset of samples 
## source: http://benjjneb.github.io/dada2/bigdata.html
filts.learn.s1 <- sample(sample.names.1, 40)
filts.learn.s2 <- sample(sample.names.2, 40)

derepFs.s1.learn <- derepFastq(filtFs.s1[filts.learn.s1], verbose=TRUE)
derepRs.s1.learn <- derepFastq(filtRs.s1[filts.learn.s1], verbose=TRUE)
derepFs.s2.learn <- derepFastq(filtFs.s2[filts.learn.s2], verbose=TRUE)
derepRs.s2.learn <- derepFastq(filtRs.s2[filts.learn.s2], verbose=TRUE)


# Sample Inference --------------------------------------------------------

dadaFs.s1.learn <- dada(derepFs.s1.learn, err=NULL, selfConsist=TRUE, multithread=8)
dadaRs.s1.learn <- dada(derepRs.s1.learn, err=NULL, selfConsist=TRUE, multithread=8)
dadaFs.s2.learn <- dada(derepFs.s2.learn, err=NULL, selfConsist=TRUE, multithread=8)
dadaRs.s2.learn <- dada(derepRs.s2.learn, err=NULL, selfConsist=TRUE, multithread=8)
rm(derepFs.s1.learn, derepRs.s1.learn, derepFs.s2.learn, derepRs.s2.learn)

# # Visualize estimated error rates
p<- plotErrors(dadaFs.s1.learn[[1]], nominalQ=TRUE)
ggsave("dada_errors_F_s1.png", plot=p)
p<- plotErrors(dadaRs.s1.learn[[1]], nominalQ=TRUE)
ggsave("dada_errors_R_s1.png", plot=p)
p<- plotErrors(dadaFs.s2.learn[[1]], nominalQ=TRUE)
ggsave("dada_errors_F_s2.png", plot=p)
p<- plotErrors(dadaRs.s2.learn[[1]], nominalQ=TRUE)
ggsave("dada_errors_R_s2.png", plot=p)

# Just keep the error profiles
errFs.s1 <- dadaFs.s1.learn[[1]]$err_out
errRs.s1 <- dadaRs.s1.learn[[1]]$err_out
errFs.s2 <- dadaFs.s2.learn[[1]]$err_out
errRs.s2 <- dadaRs.s2.learn[[1]]$err_out
rm(dadaFs.s1.learn, dadaRs.s1.learn, dadaFs.s2.learn, dadaRs.s2.learn)

# Now sample inference for entire dataset
# Run 1
derepFs.s1 <- vector("list", length(sample.names.1))
derepRs.s1 <- vector("list", length(sample.names.1))
dadaFs.s1 <- vector("list", length(sample.names.1))
dadaRs.s1 <- vector("list", length(sample.names.1))
names(dadaFs.s1) <- sample.names.1
names(dadaRs.s1) <- sample.names.1
names(derepFs.s1) <- sample.names.1
names(derepRs.s1) <- sample.names.1
for (sam in sample.names.1){
    cat("Processing:", sam, "\n")
    derepFs.s1[[sam]] <- derepFastq(filtFs.s1[[sam]])
    derepRs.s1[[sam]] <- derepFastq(filtRs.s1[[sam]])
    dadaFs.s1[[sam]] <- dada(derepFs.s1[[sam]], err=errFs.s1, multithread=8)
    dadaRs.s1[[sam]] <- dada(derepRs.s1[[sam]], err=errRs.s1, multithread=8)
}
# Run 1: Merge Paired Reads
mergers.s1 <- mergePairs(dadaFs.s1, derepFs.s1, dadaRs.s1, derepRs.s1, verbose=TRUE)
head(mergers.s1)
# Run 1: Clear up space
rm(derepFs.s1, derepRs.s1, dadaFs.s1, dadaRs.s1)


# Run 2
derepFs.s2 <- vector("list", length(sample.names.2))
derepRs.s2 <- vector("list", length(sample.names.2))
dadaFs.s2 <- vector("list", length(sample.names.2))
dadaRs.s2 <- vector("list", length(sample.names.2))
names(dadaFs.s2) <- sample.names.2
names(dadaRs.s2) <- sample.names.2
names(derepFs.s2) <- sample.names.2
names(derepRs.s2) <- sample.names.2
for (sam in sample.names.2){
    cat("Processing:", sam, "\n")
    derepFs.s2[[sam]] <- derepFastq(filtFs.s2[[sam]])
    derepRs.s2[[sam]] <- derepFastq(filtRs.s2[[sam]])
    dadaFs.s2[[sam]] <- dada(derepFs.s2[[sam]], err=errFs.s2, multithread=8)
    dadaRs.s2[[sam]] <- dada(derepRs.s2[[sam]], err=errRs.s2, multithread=8)
}
# Run 1: Merge Paired Reads
mergers.s2 <- mergePairs(dadaFs.s2, derepFs.s2, dadaRs.s2, derepRs.s2, verbose=TRUE)
head(mergers.s2)
# Run 1: Clear up space
rm(derepFs.s2, derepRs.s2, dadaFs.s2, dadaRs.s2)


# Construct Sequence Table ------------------------------------------------

seqtab.s1 <- makeSequenceTable(mergers.s1)
saveRDS(seqtab.s1, "seqtab.s1.rds")
dim(seqtab.s1)
# Inspect the distributioh of sequence lengths
table(nchar(colnames(seqtab.s1)))

seqtab.s2 <- makeSequenceTable(mergers.s2)
saveRDS(seqtab.s2, "seqtab.s2.rds")
dim(seqtab.s2)
# Inspect the distributioh of sequence lengths
table(nchar(colnames(seqtab.s2)))


# Remove Chimeras ---------------------------------------------------------

seqtab.s1.nochim <- removeBimeraDenovo(seqtab.s1, tableMethod='consensus',
                                       verbose=TRUE)
dim(seqtab.s1.nochim)
sum(seqtab.s1.nochim)/sum(seqtab.s1)
saveRDS(seqtab.s1, "seqtab.s1.nochim.rds")

seqtab.s2.nochim <- removeBimeraDenovo(seqtab.s2, tableMethod='consensus', 
                                       verbose=TRUE)
dim(seqtab.s2.nochim)
sum(seqtab.s2.nochim)/sum(seqtab.s2)
saveRDS(seqtab.s2, "seqtab.s2.nochim.rds")

# Merge Sequence Tables Together ------------------------------------------

seqtab.nochim <- mergeSequenceTables(seqtab.s1.nochim, seqtab.s2.nochim)
saveRDS(seqtab.nochim, "seqtab.nochim.rds")


# Simplify naming ---------------------------------------------------------

seqtab <- seqtab.nochim

# Assign Taxonomy ---------------------------------------------------------
# Following: http://benjjneb.github.io/dada2/species.html

# Assign using Naive Bayes RDP
taxtab <- assignTaxonomy(colnames(seqtab), 
                         '0_training/silva_nr_v123_train_set.fa.gz')

# improve with exact genus-species matches - note: Not allowing multiple
taxtab <- addSpecies(taxtab, '0_training/silva_species_assignment_v123.fa.gz', 
                     verbose=TRUE)

# How many sequences are classified at different levels? (percent)
colSums(!is.na(taxtab))/nrow(taxtab)


# Make phyloseq object ----------------------------------------------------

# Import mapping 
map1 <- read.delim('0_mapping/2016.03.16MappingFilePool1_ordered_corrected1.txt', 
                   stringsAsFactors = F)
map2 <- read.delim('0_mapping/2016.03.22MappingFile_Pool2_corrected1.txt', 
                   stringsAsFactors = F)
map1 <- map1[grepl(pattern = 'T..', x = map1$X.SampleID, fixed = TRUE),]
map2 <- map2[grepl(pattern = 'T..', x = map2$X.SampleID, fixed = TRUE),]
map2$Vessel <- as.integer(map2$Vessel)
# make sure rownames match up (map1 has some we don't care about)
map1 <- map1[,colnames(map2)] 
map <- bind_rows(map1, map2)
map <- as.data.frame(map) # without this line get sam_data slot empty error from phyloseq
rownames(map) <- map$X.SampleID

# Make refseq object and extract sequences from tables
#refseq <- DNAStringSet(colnames(seqtab))
refseq <- colnames(seqtab)
names(refseq) <- paste0('seq_', 1:length(refseq))
colnames(seqtab) <- names(refseq[match(colnames(seqtab), refseq)])
rownames(taxtab) <- names(refseq[match(rownames(taxtab), refseq)])

# Combine into phyloseq object
ps <- phyloseq(otu_table(seqtab, taxa_are_rows = FALSE),
               sample_data(map), 
               tax_table(taxtab), 
               DNAStringSet(refseq))
saveRDS(ps, 'phyloseq.rds')

# Write the taxtable, seqtable, and refseq to ascii ------------------------
write.table(seqtab, file='seqtab.nochim.tsv', quote=FALSE, sep='\t')
write.table(taxtab, file='taxtab.nochim.tsv', quote=FALSE, sep='\t')
write.table(refseq, file='refseqs.nochim.tsv', quote=FALSE, sep='\t', 
            col.names = F)



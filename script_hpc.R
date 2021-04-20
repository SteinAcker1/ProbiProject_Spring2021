### This script performs initial data handling using FASTQ data. It is VERY computationally intensive and should be run on a very powerful computer. ###

library(dada2)

path <- "~/probiData/ProGastro17"
list.files(path)

#Separating forward and reverse reads
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_1_[12]"), `[`, 1)
"Separated forward and reverse reads"

#Getting ready for trimming
forwardPrimer <- "TACGGGAGGCAGCAG"
reversePrimer <- "CCAGGGTATCTAATCC"
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names
"Got ready to trim reads"

#Trimming and filtering
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, trimLeft=c(length(forwardPrimer),length(reversePrimer)))
"Trimmed reads"

#Getting error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
"Got error rates"

#Sample inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
"Inferred sample variants"

#Merging paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
"Merged pairs"

#Making a sequence table
seqtab <- makeSequenceTable(mergers)
"Generated sequence table"

#Removing chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
seqPath <- file.path(path, "output", "probiSeqs.tsv")
write.table(seqtab.nochim, file = seqPath, sep = "\t")
"Removed chimeras"

#Assigning taxonomy with Silva database
taxa <- assignTaxonomy(seqtab.nochim, "~/taxa_dada2/silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread=TRUE)
taxaPath <- file.path(path, "output", "probiTaxa.tsv")
write.table(taxa, file = taxaPath, sep = "\t")
"Assigned taxonomy"


  
  
  

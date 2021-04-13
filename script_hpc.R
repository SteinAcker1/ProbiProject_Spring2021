library(dada2)

path <- "~/probiData_local"
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
write.table(seqtab.nochim, file = "probiSeqs.tsv", sep = "\t")
"Removed chimeras"

#Assigning taxonomy with Silva database
taxa <- assignTaxonomy(seqtab.nochim, "~/taxa_dada2/silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread=TRUE)
write.table(taxa, file = "probiTaxa.tsv", sep = "\t")
#taxa <- read.csv("~/Downloads/dada2_practice_taxa.tsv", sep = "\t")
"Assigned taxonomy"


  
  
  

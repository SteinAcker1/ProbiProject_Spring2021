echo '# Installing from Bioconductor
if(!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos="http://cran.us.r-project.org")
}
BiocManager::install(version = "3.13")
BiocManager::install("dada2")

# Installing from CRAN
if(!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", repos="http://cran.us.r-project.org")
}
devtools::install_version("tidyverse", version = "1.3.1", repos="http://cran.us.r-project.org", quiet = TRUE)
devtools::install_version("vegan", version = "2.5.7", repos="http://cran.us.r-project.org", quiet = TRUE)
devtools::install_version("ggrepel", version = "0.9.1", repos="http://cran.us.r-project.org", quiet = TRUE)
devtools::install_version("viridis", version = "0.6.1", repos="http://cran.us.r-project.org", quiet = TRUE)
devtools::install_version("tidyMicro", version = "1.47", repos="http://cran.us.r-project.org", quiet = TRUE)

### This script performs initial data handling using FASTQ data. It is VERY computationally intensive and should be run on a very powerful computer. ###

set.seed(123)
library(dada2)

path <- "~/probiData/ProGastro17"

#Separating forward and reverse reads
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_1_[12]"), `[`, 1)

#Getting ready for trimming
forwardPrimerLength <- 15
reversePrimerLength <- 16
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

#Trimming and filtering
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, trimLeft=c(forwardPrimerLength,reversePrimerLength))

#Getting error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#Sample inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#Merging paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

#Making a sequence table
seqtab <- makeSequenceTable(mergers)

#Removing chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
seqPath <- file.path(path, "output", "probiSeqs.tsv")
write.table(seqtab.nochim, file = seqPath, sep = "\t")

#Assigning taxonomy with Silva database
taxa <- assignTaxonomy(seqtab.nochim, "~/taxa_dada2/silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread=TRUE)
taxaPath <- file.path(path, "output", "probiTaxa.tsv")
write.table(taxa, file = taxaPath, sep = "\t")

' | R --no-save
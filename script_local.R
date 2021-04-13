library(dada2)
library(tidyverse)
library(ggfortify)
library(factoextra)

#Loading data (this takes a little while)
taxa.df <- read.csv("~/probiData/ProGastro17/output/probiTaxa.tsv", sep = "\t")
seqs.df <- read.csv("~/probiData/ProGastro17/output/probiSeqs.tsv", sep = "\t")

#Defining function to generate count matrix
countTaxa <- function(taxa, seqs, level) {
  count.df <- data.frame(row.names = rownames(seqs))
  sequences <- colnames(seqs)
  N <- length(sequences)
  for(i in sequences) { 
    taxon <- taxa[level][i,]
    if(is.na(taxon) == FALSE) {
      if(taxon %in% colnames(count.df)) {
        count.df[taxon] <- count.df[taxon] + seqs[,i]
      } else {
        count.df[taxon] <- seqs[,i]
      }
    }
  }
  return(count.df)
}

phylumCount.df <- countTaxa(taxa = taxa.df, seqs = seqs.df, level = "Phylum")

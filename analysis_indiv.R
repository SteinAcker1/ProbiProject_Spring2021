source("script_local.R")
library(DESeq2)
library(tidyverse)
#Use this script to do analysis for individuals

counts_indiv.df <- genusCount.df %>%
  t()

ddsFullCountTable <- DESeqDataSetFromMatrix( countData = counts_indiv.df,
                                             colData = clinical.df,
                                             design = ~ Participant + Treatment)

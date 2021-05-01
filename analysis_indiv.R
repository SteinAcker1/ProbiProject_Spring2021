source("script_local.R")
#Use this script to do analysis for individuals

counts_indiv.df <- genusCount.df %>%
  t() %>%
  filter(isFALSE(endswith("NA", row.names())))

ddsFullCountTable <- DESeqDataSetFromMatrix( countData = counts_indiv.df,
                                             colData = clinical.df,
                                             design = ~ Participant + Treatment)

output <- DESeq(ddsFullCountTable)
result <- results(output)

# Coding from scratch

# Generating a sequence matrix with NAs scrubbed out at each taxonomic level
seqsPhylum.df <- filterNA(quo(Phylum))
seqsClass.df <- filterNA(quo(Class))
seqsOrder.df <- filterNA(quo(Order))
seqsFamily.df <- filterNA(quo(Family))
seqsGenus.df <- filterNA(quo(Genus))

# Creating new count matrices using the new scrubbed sequence matrices
phylumCount.df <- countTaxa(taxa = taxa.df, seqs = seqsPhylum.df, level = "Phylum")
classCount.df <- countTaxa(taxa = taxa.df, seqs = seqsClass.df, level = "Class")
orderCount.df <- countTaxa(taxa = taxa.df, seqs = seqsOrder.df, level = "Order")
familyCount.df <- countTaxa(taxa = taxa.df, seqs = seqsFamily.df, level = "Family")
genusCount.df <- countTaxa(taxa = taxa.df, seqs = seqsGenus.df, level = "Genus")

# Converting count matrices to proportion matrices
phylumProp.df <- countToProp(phylumCount.df)
classProp.df <- countToProp(classCount.df)
orderProp.df <- countToProp(orderCount.df)
familyProp.df <- countToProp(familyCount.df)
genusProp.df <- countToProp(genusCount.df)

# Appending clinical data to each proportion matrix
phylumProp_appended.df <- appendData(phylumProp.df)
classProp_appended.df <- appendData(classProp.df)
orderProp_appended.df <- appendData(orderProp.df)
familyProp_appended.df <- appendData(familyProp.df)
genusProp_appended.df <- appendData(genusProp.df)

# Performing statistical analysis on proportion matrices
phylumChanges.df <- phylumProp_appended.df %>%
  getChange() %>%
  evalChange()
classChanges.df <- classProp_appended.df %>%
  getChange() %>%
  evalChange()
orderChanges.df <- orderProp_appended.df %>%
  getChange() %>%
  evalChange()
familyChanges.df <- familyProp_appended.df %>%
  getChange() %>%
  evalChange()
genusChanges.df <- genusProp_appended.df %>%
  getChange() %>%
  evalChange()


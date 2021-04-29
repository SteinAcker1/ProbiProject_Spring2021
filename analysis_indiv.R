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

phylumProp.df <- countToProp(phylumCount.df)
genusProp.df <- countToProp(genusCount.df)

phylumProp_appended.df <- appendData(phylumProp.df)
genusProp_appended.df <- appendData(genusProp.df)



source("script_local.R")
#Use this script to do analysis for individuals

### Preliminary DESeq2 code that may or may not be further developed
# counts_indiv.df <- genusCount.df %>%
#   t() %>%
#   filter(isFALSE(endswith("NA", row.names())))
# 
# ddsFullCountTable <- DESeqDataSetFromMatrix( countData = counts_indiv.df,
#                                              colData = clinical.df,
#                                              design = ~ Participant + Treatment)
# 
# output <- DESeq(ddsFullCountTable)
# result <- results(output)

### Getting changes in Shannon's H at genus and phylum level

genus_shannonH <- analyzeDiversityChange(tidymicro.df, "Genus")
phylum_shannonH <- analyzeDiversityChange(tidymicro.df, "Phylum")

boxplot(genus_shannonH[[1]],
        main = "Genus-level diversity change",
        sub = paste("p = ", round(genus_shannonH[[2]]$p.value, 4)),
        xlab = "Phase",
        ylab = "Change in Shannon's H")

boxplot(phylum_shannonH[[1]],
        main = "Phylum-level diversity change",
        sub = paste("p = ", round(phylum_shannonH[[2]]$p.value, 4)),
        xlab = "Phase",
        ylab = "Change in Shannon's H")

### F/B Ratio

fb <- getFBchange(phylumProp_appended.df)
boxplot(fb[[1]],
        main = "F/B ratio change",
        sub = paste("p = ", round(fb[[2]]$p.value, 4)),
        xlab = "Phase",
        ylab = "Change in F/B ratio")

### Getting taxon-level changes

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
  getTaxonChange() %>%
  evalTaxonChange()
classChanges.df <- classProp_appended.df %>%
  getTaxonChange() %>%
  evalTaxonChange()
orderChanges.df <- orderProp_appended.df %>%
  getTaxonChange() %>%
  evalTaxonChange()
familyChanges.df <- familyProp_appended.df %>%
  getTaxonChange() %>%
  evalTaxonChange()
genusChanges.df <- genusProp_appended.df %>%
  getTaxonChange() %>%
  evalTaxonChange()

# Analysis on L. plantarum

Lplantarum.mat <- genusCount.df %>%
  appendData() %>%
  filterIncompletes() %>%
  findDiffs(measure = "Firmicutes/Bacilli/Lactobacillales/Lactobacillaceae/Lactiplantibacillus",
            fold = F)

Lp_present.df <- getLPstatus()
with(Lp_present.df, chisq.test(Treatment, Lp_present))

Lp_present_active.df <- filter(Lp_present.df, Treatment == "Lp299v")
with(Lp_present_active.df, chisq.test(Overweight, Lp_present))
with(Lp_present_active.df, chisq.test(Site, Lp_present))
with(Lp_present_active.df, chisq.test(Gender, Lp_present))
t.test(BSF ~ Lp_present, data = Lp_present_active.df)


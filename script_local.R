### This script performs data fine-tuning in preparation for analysis. It can be run on a typical personal computer. ###
source("funcLibrary.R")

### Initial data handling ###

#Loading data and formatting it properly (this takes a little while)
taxa.df <- read.csv("~/probiData/ProGastro17/output/probiTaxa.tsv", sep = "\t") %>%
  mutate(Class = paste(Phylum, Class, sep = "/")) %>%
  mutate(Order = paste(Class, Order, sep = "/")) %>%
  mutate(Family = paste(Order, Family, sep = "/")) %>%
  mutate(Genus = paste(Family, Genus, sep = "/"))
seqsRaw.df <- read.csv("~/probiData/ProGastro17/output/probiSeqs.tsv", sep = "\t")

#Getting a list of IDs
ids <- c()
for(name in rownames(seqsRaw.df)) {
  idx <- strsplit(name, split = "erfext_|_lib")[[1]][2]
  idx <- paste("p", idx, sep = "_")
  ids <- c(ids, idx)
}
rownames(seqsRaw.df) <- ids

#Remove the reads where Order was not identified to clean data up
seqs.df <- filterNA()

#Cleaning up the demographic data so that it can be used in analysis
demographics.df <- read.csv("~/probiData/ProGastro17/otherInfo/probiDemographics.csv")
demographics.df$Screening.number <- as.character(demographics.df$Screening.number)
demographics.df$Overweight <- with(demographics.df, ifelse(BMI > 25, TRUE, FALSE))
BSF.df <- read.csv("~/probiData/ProGastro17/otherInfo/probiBSF.csv")
BSF.df$Screening.number <- as.character(BSF.df$Screening.number)

### Creating dataframes to be used in analysis ###

#Producing count and classic OTU dataframes
phylumCount.df <- countTaxa(taxa = taxa.df, seqs = seqs.df, level = "Phylum")
classCount.df <- countTaxa(taxa = taxa.df, seqs = seqs.df, level = "Class")
orderCount.df <- countTaxa(taxa = taxa.df, seqs = seqs.df, level = "Order")
familyCount.df <- countTaxa(taxa = taxa.df, seqs = seqs.df, level = "Family")
genusCount.df <- countTaxa(taxa = taxa.df, seqs = seqs.df, level = "Genus")

phylumOTU.df <- getOTUtable(phylumCount.df)
classOTU.df <- getOTUtable(classCount.df)
orderOTU.df <- getOTUtable(orderCount.df)
familyOTU.df <- getOTUtable(familyCount.df)
genusOTU.df <- getOTUtable(genusCount.df)

#Producing clinical data dataframe
clinical.df <- appendData(data.frame(ids = ids, row.names = ids))

#Producing tidyMicro dataframe
tidymicro.df <- tidy_micro(otu_tabs = list(Phylum = phylumOTU.df,
                                        Class = classOTU.df,
                                        Order = orderOTU.df,
                                        Family = familyOTU.df,
                                        Genus = genusOTU.df),
                        clinical = clinical.df,
                        complete_clinical = TRUE,
                        library_name = "ids")

#Splitting tidyMicro dataframe into testing periods
tidymicro1.df <- tidymicro.df %>%
  filter(TestingPeriod %in% c("1","2")) %>%
  mutate(Treatment = factor(Treatment, levels = c("PreTrial_1", "Placebo", "Lp299v")))

tidymicro2.df <- tidymicro.df %>%
  filter(TestingPeriod %in% c("3","4")) %>%
  mutate(Treatment = factor(Treatment, levels = c("PreTrial_2", "Placebo", "Lp299v")))

### Old code that may be recycled in the future ###

#Converting count matrices to proportion matrices
# phylumProp.df <- countToProp(phylumCount.df)
# genusProp.df <- countToProp(genusCount.df)
# 
# phylumProp_appended.df <- appendData(phylumProp.df)
# genusProp_appended.df <- appendData(genusProp.df)
# 
# phylumProp_appended.df$FBratio <- log2(phylumProp_appended.df$Firmicutes / phylumProp_appended.df$Bacteroidota)


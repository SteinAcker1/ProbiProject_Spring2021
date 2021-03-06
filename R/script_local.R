### This script performs data fine-tuning in preparation for analysis. It can be run on a typical personal computer. ###
source("script_hpc.R")
theme_set(theme_bw())

### Initial data handling ###

#Loading data and formatting it properly (this takes a little while)
taxa.df <- read.csv("../data/processed/probiTaxa.tsv", sep = "\t")
taxa_rownames <- rownames(taxa.df)
taxa.df <- taxa.df %>%
#taxa.df <- read.csv("~/probiData/ProGastro17/output_silva_nospecies/probiTaxa.tsv", sep = "\t") %>%
  lapply(function(x) {
    str_replace_all(x, "[:punct:]| ", "")
    }) %>% #Removing all punctuation to avoid problems with shiftRockyMtn()
  data.frame(row.names = taxa_rownames) %>%
  mutate(Class = paste(Phylum, Class, sep = "/")) %>%
  mutate(Order = paste(Class, Order, sep = "/")) %>%
  mutate(Family = paste(Order, Family, sep = "/")) %>%
  mutate(Genus = paste(Family, Genus, sep = "/"))
#seqsRaw.df <- read.csv("~/probiData/ProGastro17/output_silva_nospecies/probiSeqs.tsv", sep = "\t")
seqs.df <- read.csv("../data/processed/probiSeqs.tsv", sep = "\t")

#Getting a list of IDs
ids <- c()
for(name in rownames(seqs.df)) {
  idx <- strsplit(name, split = "erfext_|_lib")[[1]][2]
  idx <- paste("p", idx, sep = "_")
  ids <- c(ids, idx)
}
rownames(seqs.df) <- ids

#Cleaning up the demographic data so that it can be used in analysis
demographics.df <- read.csv("../data/probiDemographics.csv")
demographics.df$Screening.number <- as.character(demographics.df$Screening.number)
demographics.df$Overweight <- with(demographics.df, ifelse(BMI > 25, TRUE, FALSE))
BSF.df <- read.csv("../data/probiBSF.csv")
BSF.df$Screening.number <- as.character(BSF.df$Screening.number)

### Creating dataframes to be used in analysis ###

#Producing count and classic OTU dataframes
phylumCount.df <- countTaxa(taxa = taxa.df, seqs = seqs.df, level = "Phylum") %>%
  select(-ends_with("NA"))
classCount.df <- countTaxa(taxa = taxa.df, seqs = seqs.df, level = "Class") %>%
  select(-ends_with("NA"))
orderCount.df <- countTaxa(taxa = taxa.df, seqs = seqs.df, level = "Order") %>%
  select(-ends_with("NA"))
familyCount.df <- countTaxa(taxa = taxa.df, seqs = seqs.df, level = "Family") %>%
  select(-ends_with("NA"))
genusCount.df <- countTaxa(taxa = taxa.df, seqs = seqs.df, level = "Genus") %>%
  select(-ends_with("NA"))

phylumOTU.df <- getOTUtable(phylumCount.df)
classOTU.df <- getOTUtable(classCount.df)
orderOTU.df <- getOTUtable(orderCount.df)
familyOTU.df <- getOTUtable(familyCount.df)
genusOTU.df <- getOTUtable(genusCount.df)

#Producing proportion matrices with appended data
phylumProp.df <- countToProp(phylumCount.df)
classProp.df <- countToProp(classCount.df)
orderProp.df <- countToProp(orderCount.df)
familyProp.df <- countToProp(familyCount.df)
genusProp.df <- countToProp(genusCount.df)

phylumProp_appended.df <- appendData(phylumProp.df)
classProp_appended.df <- appendData(classProp.df)
orderProp_appended.df <- appendData(orderProp.df)
familyProp_appended.df <- appendData(familyProp.df)
genusProp_appended.df <- appendData(genusProp.df)

#Producing clinical data dataframe (including whether or not Lp299v was detected)
genusCount_appended.df <- appendData(genusCount.df)
clinical.df <- genusCount_appended.df %>%
  select(-colnames(genusCount.df)) %>%
  rownames_to_column(var = "ids")
# genusCount.df <- select(genusCount.df, str_detect("Lactiplantibacillus", negate = TRUE))

#Producing tidyMicro dataframe
tidymicro.df <- tidy_micro(otu_tabs = list(Phylum = phylumOTU.df,
                                        Class = classOTU.df,
                                        Order = orderOTU.df,
                                        Family = familyOTU.df,
                                        Genus = genusOTU.df),
                        clinical = clinical.df,
                        complete_clinical = TRUE,
                        library_name = "ids")

# Splitting tidyMicro dataframe into testing periods
tidymicro1.df <- tidymicro.df %>%
  filter(TestingPeriod %in% c("1","2")) %>%
  mutate(Treatment = factor(Treatment, levels = c("PreTrial_1", "Placebo", "Lp299v")))

tidymicro2.df <- tidymicro.df %>%
  filter(TestingPeriod %in% c("3","4")) %>%
  mutate(Treatment = factor(Treatment, levels = c("PreTrial_2", "Placebo", "Lp299v")))

# Splitting tidyMicro dataframe into relevant treatments
tidymicro_Lp299v.df <- filter(tidymicro.df, Treatment == "Lp299v")%>%
  mutate(Overweight = factor(ifelse(Overweight, "yes", "no"))) %>%
  mutate(Lp_present = factor(ifelse(Lp_present, "yes", "no"))) %>%
  mutate(Gender = factor(Gender))

tidymicro_Lp_baseline.df <- filter(tidymicro.df, (Treatment == "PreTrial_1" & TestingOrder == "Lp299v - Placebo") |
                                     (Treatment == "PreTrial_2" & TestingOrder == "Placebo - Lp299v")) %>%
  mutate(Overweight = factor(ifelse(Overweight, "yes", "no"))) %>%
  mutate(Lp_present = factor(ifelse(Lp_present, "yes", "no"))) %>%
  mutate(Gender = factor(Gender))

### This script performs data fine-tuning in preparation for analysis. It can be run on a typical personal computer. ###
library(tidyverse)
library(tidyMicro) #This version of tidyMicro was downloaded directly from the CharlieCarpenter/tidyMicro GitHub repository on 19 April 2021, rather than the CRAN repository

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
filterNA <- function(taxlevel = quo(Order)) {
  NAs <- taxa.df %>%
    filter(endsWith(!!taxlevel, "NA")) %>%
    rownames()
  filteredSeqs.df <- seqsRaw.df[setdiff(colnames(seqsRaw.df),NAs)]
  return(filteredSeqs.df)
}
seqs.df <- filterNA()

#Cleaning up the demographic data so that it can be used in analysis
demographics.df <- read.csv("~/probiData/ProGastro17/otherInfo/probiDemographics.csv")
demographics.df$Screening.number <- as.character(demographics.df$Screening.number)
demographics.df$Overweight <- with(demographics.df, ifelse(BMI > 25, TRUE, FALSE))
BSF.df <- read.csv("~/probiData/ProGastro17/otherInfo/probiBSF.csv")
BSF.df$Screening.number <- as.character(BSF.df$Screening.number)



### Defining some important functions for data cleaning ###

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

#Defining function to generate a classic OTU table so tidyMicro can parse the data
getOTUtable <- function(df) {
  otus.df <- t(df) %>%
    data.frame() %>%
    rownames_to_column("OTU")
  return(otus.df)
}

#Defining function to append demographic data
appendData <- function(df) {
  #Create empty vectors
  ScreenNums <- c()
  TestingPeriod <- c()
  TestingOrder <- c()
  Treatment <- c()
  Gender <- c()
  Age <- c()
  BMI <- c()
  Site <- c()
  BSF <- c()
  Overweight <- c()
  rows <- df[,1]
  #Looping through rows
  for(i in rows) {
    #Getting info from names
    id <- strsplit(i, split = "_")
    screeningNum <- id[[1]][2]
    period <- id[[1]][3]
    #Grabbing demographic info
    info <- subset(demographics.df, Screening.number == screeningNum)
    group <- info$TestingOrder
    if(period == "1") {
      treat <- "PreTrial_1"
    } else if(period == "3") {
      treat <- "PreTrial_2"
    } else if(isTRUE(period == "2" & group == "Placebo - Lp299v") | 
              isTRUE(period == "4" & group == "Lp299v - Placebo")) {
      treat <- "Placebo"
    } else {
      treat <- "Lp299v"
    }
    indiv_gender <- info$Gender
    indiv_age <- info$Age
    indiv_bmi <- info$BMI
    indiv_site <- info$Site
    indiv_overweight <- info$Overweight
    #Grabbing BSF info
    bsf_info <- subset(BSF.df, Screening.number == screeningNum)
    if(period == "1") {
      indiv_bsf <- bsf_info$BSF1
    } else if(period == "2") {
      indiv_bsf <- bsf_info$BSF2
    } else if(period == "3") {
      indiv_bsf <- bsf_info$BSF3
    } else {
      indiv_bsf <- bsf_info$BSF4
    }
    #Appending grabbed info to vectors
    ScreenNums[i] <- paste("p", screeningNum, sep = "_")
    TestingPeriod[i] <- period
    TestingOrder[i] <- group
    Treatment[i] <- treat
    Gender[i] <- indiv_gender
    Age[i] <- indiv_age
    BMI[i] <- indiv_bmi
    Site[i] <- indiv_site
    BSF[i] <- indiv_bsf
    Overweight[i] <- indiv_overweight
  }
  #Appending the new columns to the dataframe
  df$Participant <- ScreenNums
  df$TestingPeriod <- TestingPeriod
  df$TestingOrder <- TestingOrder
  df$Treatment <- Treatment
  df$Gender <- Gender
  df$Age <- Age
  df$BMI <- BMI
  df$Site <- Site
  df$BSF <- BSF
  df$Overweight <- Overweight
  #Return the altered dataframe
  return(df)
}



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
clinical.df <- appendData(data.frame(ids))

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

#Defining function to get proportions dataframe rather than count dataframe
#NOTE: as of right now, this is not used in analysis
# countToProp <- function(count) {
#   prop.matrix <- t(apply(count, 1, FUN = function(vec) vec / sum(vec)))
#   return(data.frame(prop.matrix, row.names = rownames(count)))
# }

#Converting count matrices to proportion matrices
# phylumProp.df <- countToProp(phylumCount.df)
# genusProp.df <- countToProp(genusCount.df)
# 
# phylumProp_appended.df <- appendData(phylumProp.df)
# genusProp_appended.df <- appendData(genusProp.df)
# 
# phylumProp_appended.df$FBratio <- log2(phylumProp_appended.df$Firmicutes / phylumProp_appended.df$Bacteroidota)


library(dada2)
library(tidyverse)
library(ggfortify)
library(factoextra)
library(tidyMicro)

#Loading data (this takes a little while)
taxa.df <- read.csv("~/probiData/ProGastro17/output/probiTaxa.tsv", sep = "\t") %>%
  mutate(Class = paste(Phylum, "/", Class, sep = "")) %>%
  mutate(Order = paste(Class, "/", Order, sep = "")) %>%
  mutate(Family = paste(Order, "/", Family, sep = "")) %>%
  mutate(Genus = paste(Family, "/", Genus, sep = ""))
seqs.df <- read.csv("~/probiData/ProGastro17/output/probiSeqs.tsv", sep = "\t")
ids <- c()
for(name in rownames(seqs.df)) {
  idx <- strsplit(name, split = "erfext_|_lib")[[1]][2]
  idx <- paste("p", idx, sep = "_")
  ids <- c(ids, idx)
}
rownames(seqs.df) <- ids
getOTUtable <- function(df) {
  otus.df <- t(df) %>%
    data.frame() %>%
    rownames_to_column("OTU")
}
  

demographics.df <- read.csv("~/probiData/ProGastro17/otherInfo/probiDemographics.csv")
demographics.df$Screening.number <- as.character(demographics.df$Screening.number)
BSF.df <- read.csv("~/probiData/ProGastro17/otherInfo/probiBSF.csv")
BSF.df$Screening.number <- as.character(BSF.df$Screening.number)

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

#Defining function to get proportions dataframe rather than count dataframe
countToProp <- function(count) {
  prop.matrix <- t(apply(count, 1, FUN = function(vec) vec / sum(vec)))
  return(data.frame(prop.matrix, row.names = rownames(count)))
}

#Defining function to generate a classic OTU table
getOTUtable <- function(df) {
  otus.df <- t(df) %>%
    data.frame() %>%
    rownames_to_column("OTU")
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
    }
     else if(isTRUE(period == "2" & group == "Placebo - Lp299v") | isTRUE(period == "4" & group == "Lp299v - Placebo")) {
      treat <- "Placebo"
    } else {
      treat <- "Lp299v"
    }
    indiv_gender <- info$Gender
    indiv_age <- info$Age
    indiv_bmi <- info$BMI
    indiv_site <- info$Site
    #Grabbing BSF info
    bsf_info <- subset(BSF.df, Screening.number == screeningNum)
    if(period == "1") indiv_bsf <- bsf_info$BSF1
    else if(period == "2") indiv_bsf <- bsf_info$BSF2
    else if(period == "3") indiv_bsf <- bsf_info$BSF3
    else indiv_bsf <- bsf_info$BSF4
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
  #Return the altered dataframe
  return(df)
}

#Producing dataframes
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

clinical.df <- appendData(data.frame(ids))

tidymicro.df <- tidy_micro(otu_tabs = list(Phylum = phylumOTU.df,
                                        Class = classOTU.df,
                                        Order = orderOTU.df,
                                        Family = familyOTU.df,
                                        Genus = genusOTU.df),
                        clinical = clinical.df,
                        complete_clinical = TRUE,
                        library_name = "ids")

tidymicro1.df <- tidymicro.df %>%
  filter(TestingPeriod %in% c("1","2")) %>%
  mutate(Treatment = factor(Treatment, levels = c("PreTrial_1", "Placebo", "Lp299v")))

tidymicro2.df <- tidymicro.df %>%
  filter(TestingPeriod %in% c("3","4")) %>%
  mutate(Treatment = factor(Treatment, levels = c("PreTrial_2", "Placebo", "Lp299v")))


# phylumProp.df <- countToProp(phylumCount.df)
# genusProp.df <- countToProp(genusCount.df)
# 
# phylumProp_appended.df <- appendData(phylumProp.df)
# genusProp_appended.df <- appendData(genusProp.df)
# 
# phylumProp_appended.df$FBratio <- log2(phylumProp_appended.df$Firmicutes / phylumProp_appended.df$Bacteroidota)


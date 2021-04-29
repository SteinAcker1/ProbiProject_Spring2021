# This script contains essential functions for other scripts to function
library(tidyverse)
library(tidyMicro) #This version of tidyMicro was downloaded directly from the CharlieCarpenter/tidyMicro GitHub repository on 19 April 2021, rather than the CRAN repository
library(DESeq2)

# Eliminates reads which DADA2 failed to identify at the Order level or below
filterNA <- function(taxlevel = quo(Order)) {
  NAs <- taxa.df %>%
    filter(endsWith(!!taxlevel, "NA")) %>%
    rownames()
  filteredSeqs.df <- seqsRaw.df[setdiff(colnames(seqsRaw.df),NAs)]
  return(filteredSeqs.df)
}

# Converts DADA2 output to a count table
countTaxa <- function(taxa, seqs, level) {
  count.df <- data.frame(row.names = rownames(seqs))
  sequences <- colnames(seqs)
  N <- length(sequences)
  for(i in sequences) { 
    taxon <- taxa[level][i,]
    if(!is.na(taxon)) {
      if(taxon %in% colnames(count.df)) {
        count.df[taxon] <- count.df[taxon] + seqs[,i]
      } else {
        count.df[taxon] <- seqs[,i]
      }
    }
  }
  return(count.df)
}

# Converts count table to a traditional OTU table
getOTUtable <- function(df) {
  otus.df <- t(df) %>%
    data.frame() %>%
    rownames_to_column("OTU")
  return(otus.df)
}

# Converts count table to proportion table
countToProp <- function(count) {
  prop.df <- count %>%
    apply(1, FUN = function(vec) vec / sum(vec)) %>%
    t() %>%
    data.frame(row.names = rownames(count))
  return(prop.df)
}

# Appends clinical data to a given dataframe
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
  rows <- rownames(df)
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

# Calculates the change in abundance of each taxon in each individual in response to both placebo and Lp299v
# STILL IN PROGRESS
getChange <- function(prop.df, exclude = TRUE) {
  if(exclude) {
    completeness <- prop.df$Participant %>%
      factor() %>%
      summary()
    participants <- names(completeness)
    incomplete_rows <- c()
    for(i in participants) {
      if(completeness[i] != 4) {
        incomplete_rows <- c(incomplete_rows, i)
      }
    }
    prop.df <- filter(prop.df, !(Participant %in% incomplete_rows))
  }
  changes <- list()
  N <- ncol(prop.df)
  for(i in 1:N) {
    mat <- matrix(ncol = 2, nrow = nrow(prop.df) / 4)
    colnames(mat) <- c("Placebo", "Treatment")
    
  }
}





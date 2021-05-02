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

filterIncompletes <- function(data.df) {
  raw_participants <- unique(data.df$Participant)
  # filter out participants who did not provide all 4 samples
  completeness <- data.df$Participant %>%
    factor() %>%
    summary()
  incomplete_rows <- c()
  for(i in raw_participants) {
    if(completeness[i] != 4) {
      incomplete_rows <- c(incomplete_rows, i)
    }
  }
  clean.df <- filter(data.df, !(Participant %in% incomplete_rows))
  return(clean.df)
}

findDiffs <- function(data.df, measure, particips) {
  mat <- matrix(ncol = 2, nrow = nrow(data.df) / 4)
  colnames(mat) <- c("Placebo", "Treatment")
  for(i in 1:length(particips)) {
    participantInfo.df <- data.df %>%
      filter(Participant == particips[i])
    if(participantInfo.df[1,]$TestingOrder == "Placebo - Lp299v") {
      placebo_change <- filter(participantInfo.df, Treatment == "Placebo")[,measure] -
        filter(participantInfo.df, Treatment == "PreTrial_1")[,measure]
      treatment_change <- filter(participantInfo.df, Treatment == "Lp299v")[,measure] -
        filter(participantInfo.df, Treatment == "PreTrial_2")[,measure]
    } else {
      treatment_change <- filter(participantInfo.df, Treatment == "Lp299v")[,measure] -
        filter(participantInfo.df, Treatment == "PreTrial_1")[,measure]
      placebo_change <- filter(participantInfo.df, Treatment == "Placebo")[,measure] -
        filter(participantInfo.df, Treatment == "PreTrial_2")[,measure]
    }
    mat[i,1] <- placebo_change
    mat[i,2] <- treatment_change
  }
  return(mat)
}

# Calculates the change in abundance of each taxon in each individual in response to both placebo and Lp299v
# This function returns a list with matrices for each taxon showing the change in abundance in response to both the treatment and placebo
getTaxonChange <- function(prop.df, exclude = TRUE) {
  prop.df <- filterIncompletes(prop.df)
  participants <- unique(prop.df$Participant)
  changes <- list()
  taxa <- setdiff(colnames(prop.df), colnames(clinical.df))
  for(taxon in taxa) {
    changes[[taxon]] <- findDiffs(prop.df, taxon, participants)
  }
  return(changes)
}

# Calculates essential statistics about changes in each taxon's abundance
# This function runs a paired Student's t-test for each matrix generated by getChange() and
#   performs a Benjamini-Hochberg correction on the p-values to account for multiple testing
evalTaxonChange <- function(change_list) {
  N_taxa <- length(change_list)
  measures <- c("taxon", "diff", "CI95_lower", "CI95_upper", "p_val")
  results_mat <- matrix(nrow = N_taxa, ncol = length(measures))
  colnames(results_mat) <- measures
  for(i in 1:N_taxa) {
    matrix <- change_list[[i]]
    taxon <- names(change_list[i])
    t_results <- t.test(matrix[,"Treatment"], matrix[,"Placebo"], paired = TRUE) #outputs mean change from using treatment over placebo
    diff <- t_results$estimate[[1]]
    CI95_lower <- t_results$conf.int[1]
    CI95_upper <- t_results$conf.int[2]
    p_val <- t_results$p.value
    results_mat[i,] <- c(taxon, diff, CI95_lower, CI95_upper, p_val)
  }
  p_adj <- p.adjust(results_mat[,"p_val"], method = "fdr")
  results_mat <- cbind(results_mat, p_adj)
  results.df <- data.frame(results_mat)
  return(results.df)
}

analyzeDiversityChange <- function(tidy.df,
                                   taxlevel,
                                   mindepth = 5000,
                                   statistic = "ShannonH") {
  tidy.df <- filter(tidy.df, Table == taxlevel & str_detect(Taxa, "NA$", negate = TRUE))
  alpha.df <- alpha_div(tidy.df, table = taxlevel, min_depth = mindepth)
  diversity.df <- distinct(data.frame(Participant = alpha.df$Participant,
                                      Diversity = alpha.df[,statistic],
                                      Treatment = alpha.df$Treatment,
                                      TestingOrder = alpha.df$TestingOrder))
  diversity.df <- filterIncompletes(diversity.df)
  participants <- unique(diversity.df$Participant)
  diversityChange.mat <- findDiffs(data.df = diversity.df,
                                   measure = "Diversity",
                                   particips = participants)
  t_results <- t.test(diversityChange.mat[,"Treatment"],
                      diversityChange.mat[,"Placebo"],
                      paired = TRUE)
  output <- list(diversityChange.mat, t_results)
  return(output)
}



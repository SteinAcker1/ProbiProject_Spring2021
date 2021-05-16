# This script contains essential functions for other scripts to function
library(tidyverse)
library(tidyMicro) #This version of tidyMicro was downloaded directly from the CharlieCarpenter/tidyMicro GitHub repository on 19 April 2021, rather than the CRAN repository
library(DESeq2)
library(vegan)
library(ggvegan) #This version of tidyMicro was downloaded directly from the gavinsimpson/ggvegan GitHub repository on 10 May 2021
library(viridis)
library(ggfortify)

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

# Generate a list of individuals who had L. plantarum present in gut in testing phase

getLPstatusList <- function(df_appended) {
  LpList <- df_appended %>%
    filter(Treatment == "Lp299v") %>%
    select(Participant, Lp_present) %>%
    remove_rownames() %>%
    column_to_rownames("Participant") %>%
    t()
  return(LpList)
}

# Appends clinical data to a given dataframe
appendData <- function(df) {
  #Save Lp299v as a variable and see if it should be evaluated
  Lp_names <- c("Firmicutes/Bacilli/Lactobacillales/Lactobacillaceae/Lactiplantibacillus",
                "Firmicutes.Bacilli.Lactobacillales.Lactobacillaceae.Lactiplantibacillus")
  if(Lp_names[1] %in% colnames(df)) {
    Lp_name <- Lp_names[1]
    eval_Lp <- TRUE
  } else if(Lp_names[2] %in% colnames(df)) {
    Lp_name <- Lp_names[2]
    eval_Lp <- TRUE
  } else {
    eval_Lp <- FALSE
  }
  if(eval_Lp) Lp_present <- c()
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
    #Appending value to see if L. plantarum is present
    if(eval_Lp) {
      if(df[i,Lp_name] != 0) {
        Lplantarum = TRUE
      } else {
        Lplantarum = FALSE
      }
      Lp_present[i] <- Lplantarum
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
  if(eval_Lp) {
    df$Lp_present <- Lp_present
    LpList.df <- getLPstatusList(df)
    Lp_present_treat <- c()
    for(i in df$Participant) {
      print(i)
      if(i %in% colnames(LpList.df)) {
        Lp_present_treat <- c(Lp_present_treat, LpList.df[,i])
      } else {
        Lp_present_treat <- c(Lp_present_treat, NA)
      }
    }
    df$Lp_present_treat <- Lp_present_treat
  }
  #Return the altered dataframe
  return(df)
}

# filter out participants who did not provide all 4 samples
filterIncompletes <- function(data.df) {
  raw_participants <- unique(data.df$Participant)
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

# Finds the change in a certain measure for each individual
findDiffs <- function(data.df, measure, particips = unique(data.df$Participant), fold = FALSE) {
  if(fold) warning("You have selected the \"fold\" option. This adds a pseudocount of +1 to both the numerator and denominator, so it should not be used with proportion tables.")
  mat <- matrix(ncol = 2, nrow = nrow(data.df) / 4)
  colnames(mat) <- c("Placebo", "Treatment")
  for(i in 1:length(particips)) {
    participantInfo.df <- data.df %>%
      filter(Participant == particips[i])
    if(fold) {
      if(participantInfo.df[1,]$TestingOrder == "Placebo - Lp299v") {
        placebo_change <- (filter(participantInfo.df, Treatment == "Placebo")[,measure] + 1) /
          (filter(participantInfo.df, Treatment == "PreTrial_1")[,measure] + 1)
        treatment_change <- (filter(participantInfo.df, Treatment == "Lp299v")[,measure] + 1) /
          (filter(participantInfo.df, Treatment == "PreTrial_2")[,measure] + 1)
      } else {
        treatment_change <- (filter(participantInfo.df, Treatment == "Lp299v")[,measure] + 1) /
          (filter(participantInfo.df, Treatment == "PreTrial_1")[,measure] + 1)
        placebo_change <- (filter(participantInfo.df, Treatment == "Placebo")[,measure] + 1) /
          (filter(participantInfo.df, Treatment == "PreTrial_2")[,measure] + 1)
      }
    } else {
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
    }
    mat[i,1] <- placebo_change
    mat[i,2] <- treatment_change
  }
  return(mat)
}

# Calculates the change in abundance of each taxon in each individual in response to both placebo and Lp299v
# This function returns a list with matrices for each taxon showing the change in abundance in response to both the treatment and placebo
getTaxonChange <- function(prop.df) {
  prop.df <- filterIncompletes(prop.df)
  participants <- unique(prop.df$Participant)
  changes <- list()
  taxa <- setdiff(colnames(prop.df), c(colnames(clinical.df), "Lp_present"))
  for(taxon in taxa) {
    changes[[taxon]] <- findDiffs(prop.df, taxon, participants)
  }
  class(changes) <- "TaxonShifts"
  return(changes)
}

# Calculates essential statistics about changes in each taxon's abundance
# This function runs a paired Student's t-test for each matrix generated by getChange() and
#   performs a Benjamini-Hochberg correction on the p-values to account for multiple testing
evalTaxonChange <- function(change_list) {
  if(class(change_list) != "TaxonShifts") {
    stop("Only the output from getTaxonChange() may be used with this function")
  }
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

getFBchange <- function(phylum.df) {
  phylum.df <- filterIncompletes(phylum.df)
  participants <- unique(phylum.df$Participant)
  phylum.df$FBratio <- phylum.df$Firmicutes / phylum.df$Bacteroidota
  FBchange.mat <- findDiffs(data.df = phylum.df,
                            measure = "FBratio",
                            particips = participants)
  t_results <- t.test(FBchange.mat[,"Treatment"],
                      FBchange.mat[,"Placebo"],
                      paired = TRUE)
  output <- list(FBchange.mat, t_results)
  return(output)
}

getLPstatus <- function(df = genusProp_appended.df) {
  treatmentNew <- c()
  for(i in 1:nrow(df)) {
    row <- df[i,]
    if(isTRUE(row["TestingPeriod"] == "1" & row["TestingOrder"] == "Placebo - Lp299v") | 
       isTRUE(row["TestingPeriod"] == "3" & row["TestingOrder"] == "Lp299v - Placebo")) {
      treatmentNew[i] <- "PrePlacebo"
    } else if(isTRUE(row["TestingPeriod"] == "3" & row["TestingOrder"] == "Placebo - Lp299v") | 
              isTRUE(row["TestingPeriod"] == "1" & row["TestingOrder"] == "Lp299v - Placebo")) {
      treatmentNew[i] <- "PreLp299v"
    } else {
      treatmentNew[i] <- row["Treatment"][[1]]
    }
  }
  lp.df <- df %>%
    select(Participant,
           TestingPeriod,
           TestingOrder,
           Treatment,
           Gender,
           Age,
           BMI,
           Site,
           BSF,
           Overweight,
           Lp_present) %>%
    mutate(Treatment = treatmentNew)
  return(lp.df)
}

getDemographicDiff <- function(data.df, var, treat = "PreTrial_1") {
  treat.df <- data.df %>%
    appendData() %>%
    filter(Treatment %in% treat) %>%
    select_if(function(col) max(col) != 0)
  taxa <- setdiff(colnames(treat.df), c(colnames(clinical.df), "Lp_present"))
  tax_name <- c()
  group1_name <- c()
  group1_mean <- c()
  group2_name <- c()
  group2_mean <- c()
  p_value <- c()
  for(taxon in taxa) {
    t_results <- t.test(eval(parse(text = taxon)) ~ eval(parse(text = var)), data = treat.df)
    tax_name <- c(tax_name, taxon)
    group1_name <- c(group1_name, names(t_results$estimate[1]))
    group2_name <- c(group2_name, names(t_results$estimate[2]))
    group1_mean <- c(group1_mean, t_results$estimate[1])
    group2_mean <- c(group2_mean, t_results$estimate[2])
    p_value <- c(p_value, t_results$p.value)
  }
  results.df <- data.frame(tax_name,
                           group1_name,
                           group1_mean,
                           group2_name,
                           group2_mean,
                           p_value) %>%
    mutate(p_adj = p.adjust(p_value, method = "fdr")) 
  return(results.df)
}


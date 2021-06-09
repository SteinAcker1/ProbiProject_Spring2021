:"
TITLE: Probi ProGastro17 Analysis Executable Readme
AUTHOR: Stein Acker
DATE: 9 June 2021
REQUIRED R VERSION: 4.1.0
INSTRUCTIONS FOR USE:
1) Ensure you have the right R version on your computer (4.1.0)
2) Ensure you have access to the internet during analysis
3) Ensure your file tree matches the one below (reference taxa may be downloaded with 'wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_wSpecies_train_set.fa.gz', participant data must be accessed from Probi)
4) Navigate to the directory you are using for analysis
5) Enter ./executableReadme.sh
6) Done!

GITHUB LINK: https://github.com/SteinAcker1/ProbiProject_Spring2021

REQUIRED FILE TREE:
.
├── data
│   ├── *participant FASTQ data*
│   ├── probiBSF.csv
│   └── probiDemographics.csv
├── executableReadme.sh
└── refTaxa
    └── silva_nr99_v138.1_wSpecies_train_set.fa.gz
    
RESULTING FILE TREE:
.
├── data
│   ├── *participant FASTQ data*
│   ├── filtered
│   │   └── *filtered data*
│   ├── processed
│   │   ├── probiTaxa.tsv
│   │   └── probiSeqs.tsv
│   ├── probiBSF.csv
│   └── probiDemographics.csv
├── output
│   ├── Fig2.jpeg
│   ├── Fig3A.jpeg
│   ├── Fig3B.jpeg
│   ├── Fig4.jpeg
│   ├── statistics.txt
│   ├── suppFig1A.jpeg
│   ├── suppFig1B.jpeg
│   └── suppFig2.jpeg
├── executableReadme.sh
└── refTaxa
    └── silva_nr99_v138.1_wSpecies_train_set.fa.gz"

echo '# Installing from Bioconductor
if(!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos="http://cran.us.r-project.org", quiet = TRUE)
}
BiocManager::install(version = "3.13")
BiocManager::install("dada2")

# Installing from CRAN
if(!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", repos="http://cran.us.r-project.org", quiet = TRUE)
}
devtools::install_version("tidyverse", version = "1.3.1", repos="http://cran.us.r-project.org", quiet = TRUE)
devtools::install_version("vegan", version = "2.5.7", repos="http://cran.us.r-project.org", quiet = TRUE)
devtools::install_version("ggrepel", version = "0.9.1", repos="http://cran.us.r-project.org", quiet = TRUE)
devtools::install_version("viridis", version = "0.6.1", repos="http://cran.us.r-project.org", quiet = TRUE)
devtools::install_version("tidyMicro", version = "1.47", repos="http://cran.us.r-project.org", quiet = TRUE)
devtools::install_version("latex2exp", version = "0.5.0", repos="http://cran.us.r-project.org", quiet = TRUE)

# Loading packages
library(dada2, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(tidyMicro, quietly = TRUE)
library(vegan, quietly = TRUE)
library(viridis, quietly = TRUE)
library(ggrepel, quietly = TRUE)
library(latex2exp, quietly = TRUE)

### This script performs initial data handling using FASTQ data. It is VERY computationally intensive and should be run on a very powerful computer. ###

set.seed(456)

path <- "data/"

#Separating forward and reverse reads
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_1_[12]"), `[`, 1)

#Getting ready for trimming
forwardPrimerLength <- 15
reversePrimerLength <- 16
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

#Trimming and filtering
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, trimLeft=c(forwardPrimerLength,reversePrimerLength))

#Getting error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#Sample inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#Merging paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

#Making a sequence table
seqtab <- makeSequenceTable(mergers)

#Removing chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
seqPath <- file.path(path, "processed", "probiSeqs.tsv")
write.table(seqtab.nochim, file = seqPath, sep = "\t")

#Assigning taxonomy with Silva database
taxa <- assignTaxonomy(seqtab.nochim, "refTaxa/silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread=TRUE)
taxaPath <- file.path(path, "processed", "probiTaxa.tsv")
write.table(taxa, file = taxaPath, sep = "\t")

# This script contains essential functions for other scripts to function


# Converts DADA2 output to a count table
countTaxa <- function(taxa, seqs, level) {
  count.df <- data.frame(row.names = rownames(seqs))
  sequences <- colnames(seqs)
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
  } else if(exists("clinical.df")) {
    if("Lp_present" %in% colnames(clinical.df)) {
    df$Lp_present <- clinical.df$Lp_present
    df$Lp_present_treat <- clinical.df$Lp_present_treat
    }
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

# Calculates essential statistics about changes in each taxons abundance
# This function runs a paired Students t-test for each matrix generated by getChange() and
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
  var_sym <- rlang::sym(var)
  treat.df <- data.df %>%
    filter(Treatment %in% treat)%>%
    filter(!is.na(!!var_sym)) %>%
    select_if(function(col) max(col, na.rm = TRUE) != 0)
  taxa <- setdiff(colnames(treat.df), c(colnames(clinical.df)))
  taxon <- c()
  group1_name <- c()
  group1_mean <- c()
  group2_name <- c()
  group2_mean <- c()
  p_val <- c()
  for(i in taxa) {
    t_results <- t.test(eval(parse(text = i)) ~ eval(parse(text = var)), data = treat.df)
    taxon <- c(taxon, i)
    group1_name <- c(group1_name, names(t_results$estimate[1]))
    group2_name <- c(group2_name, names(t_results$estimate[2]))
    group1_mean <- c(group1_mean, t_results$estimate[1])
    group2_mean <- c(group2_mean, t_results$estimate[2])
    p_val <- c(p_val, t_results$p.value)
  }
  results.df <- data.frame(taxon,
                           group1_name,
                           group1_mean,
                           group2_name,
                           group2_mean,
                           p_val) %>%
    mutate(p_adj = p.adjust(p_val, method = "fdr")) %>%
    mutate(diff = group2_mean - group1_mean)
  return(results.df)
}

shiftRockyMtn <- function(data.changes) {
  data.changes <- data.changes %>%
    mutate(Phylum = taxon %>%
             str_split("[:punct:]") %>%
             lapply("[[", 1) %>%
             unlist()
    ) %>%
    mutate(Lowest_level = taxon %>%
             str_split("[:punct:]") %>%
             lapply(tail, 1) %>%
             unlist()
    ) %>%
    mutate(p_val = as.numeric(p_val)) %>%
    mutate(diff = as.numeric(diff)) %>%
    mutate(p_inv = 1 / p_val) %>%
    mutate(p_log_inv = log10(p_inv) * sign(diff))
    print(data.changes$Lowest_level)
  p <- ggplot(data = data.changes, mapping = aes(x = taxon, y = p_log_inv)) +
    theme(axis.text.x = element_blank()) +
    geom_label_repel(data = filter(data.changes, p_val < 0.01),
                     mapping = aes(label = Lowest_level)) +
    geom_segment(mapping = aes(xend = taxon, yend = 0, color = Phylum)) +
    geom_point(data = filter(data.changes, p_val < 0.01))
  output <- list(data.changes, p)
  return(output)
}

theme_set(theme_bw())

### Initial data handling ###

#Loading data and formatting it properly (this takes a little while)
taxa.df <- read.csv("data/processed/probiTaxa.tsv", sep = "\t")
taxa_rownames <- rownames(taxa.df)
taxa.df <- taxa.df %>%
  lapply(function(x) {
    str_replace_all(x, "[:punct:]| ", "")
    }) %>% #Removing all punctuation to avoid problems with shiftRockyMtn()
  data.frame(row.names = taxa_rownames) %>%
  mutate(Class = paste(Phylum, Class, sep = "/")) %>%
  mutate(Order = paste(Class, Order, sep = "/")) %>%
  mutate(Family = paste(Order, Family, sep = "/")) %>%
  mutate(Genus = paste(Family, Genus, sep = "/"))
seqs.df <- read.csv("data/processed/probiSeqs.tsv", sep = "\t")

#Getting a list of IDs
ids <- c()
for(name in rownames(seqs.df)) {
  idx <- strsplit(name, split = "erfext_|_lib")[[1]][2]
  idx <- paste("p", idx, sep = "_")
  ids <- c(ids, idx)
}
rownames(seqs.df) <- ids

#Cleaning up the demographic data so that it can be used in analysis
demographics.df <- read.csv("data/probiDemographics.csv")
demographics.df$Screening.number <- as.character(demographics.df$Screening.number)
demographics.df$Overweight <- with(demographics.df, ifelse(BMI > 25, TRUE, FALSE))
BSF.df <- read.csv("data/probiBSF.csv")
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
  
### Figures ###

# Fig. 2

Lp_present.df <- getLPstatus()

Lp_present_count.df <- Lp_present.df %>%
  with(table(Treatment, Lp_present, TestingOrder)) %>%
  data.frame()

Fig2.df <- Lp_present_count.df %>%
  group_by(Treatment, TestingOrder) %>%
  mutate(pct = 100 * Freq / sum(Freq)) %>%
  filter(Lp_present == TRUE) %>%
  select(-Lp_present)

jpeg(filename = "output/Fig2.jpeg", width = 667, height = 453)
ggplot(data = Fig2.df, mapping = aes(fill = TestingOrder, y = pct, x = Treatment)) +
  geom_bar(stat = "identity", position = "dodge")
dev.off()

# Fig. 3A

rmplot_Lp_present <- genusProp_appended.df %>%
  getDemographicDiff(var = "Lp_present", treat = "Lp299v") %>%
  filter(str_detect(taxon, "Lactiplantibacillus$", negate = TRUE)) %>%
  shiftRockyMtn()

jpeg(filename = "output/Fig3A.jpeg", width = 667, height = 453)
rmplot_Lp_present[[2]] +
  ylab(TeX(r"($log(1/p) * sgn(\Delta)$)")) +
  labs(tag = "A")
dev.off()

# Fig. 3B

rmplot_Lp_present_baseline <- genusProp_appended.df %>%
  mutate(Treatment = ifelse((Treatment == "PreTrial_1" & TestingOrder == "Lp299v - Placebo") |
                              (Treatment == "PreTrial_2" & TestingOrder == "Placebo - Lp299v"),
                            "PreLp299v",
                            Treatment)) %>%
  getDemographicDiff(var = "Lp_present_treat", treat = "PreLp299v") %>%
  filter(str_detect(taxon, "Lactiplantibacillus$", negate = TRUE)) %>%
  shiftRockyMtn()

jpeg(filename = "output/Fig3B.jpeg", width = 667, height = 453)
rmplot_Lp_present_baseline[[2]] +
  ylab(TeX(r"($log(1/p) * sgn(\Delta)$)")) +
  labs(tag = "B")
dev.off()

# Fig. 4

genus_alpha.df <- tidymicro.df %>%
  alpha_div(table = "Genus", min_depth = 5000) %>%
  select(c("Treatment", "TestingOrder", "ShannonH", "Participant", "Lp_present_treat")) %>%
  distinct()

genus_alpha_boxplot.df <- genus_alpha.df %>%
  filter((Treatment == "PreTrial_1" & TestingOrder == "Lp299v - Placebo") |
           (Treatment == "PreTrial_2" & TestingOrder == "Placebo - Lp299v") |
           Treatment == "Lp299v") %>%
  mutate(Treatment = ifelse(Treatment == "Lp299v", Treatment, "Baseline")) %>%
  filter(!is.na(Lp_present_treat))

jpeg(filename = "output/Fig4.jpeg", width = 667, height = 453)
ggplot(data = genus_alpha_boxplot.df, mapping = aes(x = Treatment, y = ShannonH, fill = Lp_present_treat)) +
  geom_boxplot()
dev.off()

# Supp. fig. 1A

jpeg(filename = "output/suppFig1A.jpeg", width = 667, height = 453)
ra_bars(tidymicro_Lp299v.df,
        table = "Genus",
        Lp_present_treat,
        top_taxa = 15) +
  labs(tag = "A")
dev.off()

# Supp. fig. 1B

jpeg(filename = "output/suppFig1B.jpeg", width = 667, height = 453)
ra_bars(tidymicro_Lp_baseline.df,
        table = "Genus",
        Lp_present_treat,
        top_taxa = 15) +
  labs(tag = "B")
dev.off()

# Supp. fig. 2

jpeg(filename = "output/suppFig2.jpeg", width = 667, height = 453)
Lp299v_beta <- beta_div(tidymicro_Lp299v.df, table = "Genus", method = "bray")
Lp299v_beta %>%
  beta_heatmap(micro_set = tidymicro_Lp299v.df, Lp_present) +
  scale_fill_viridis(option = "magma")
dev.off()' | R --save

echo '
library(tidyverse, quietly = TRUE)
library(tidyMicro, quietly = TRUE)
library(vegan, quietly = TRUE)

### Basic demographic info (Table 1 data) ###

# Age by group

"Age"
clinical.df %>%
  filter(Treatment == "PreTrial_1") %>%
  with(tapply(Age, TestingOrder, summary))

# BMI by group

"BMI"
clinical.df %>%
  filter(Treatment == "PreTrial_1") %>%
  with(tapply(BMI, TestingOrder, summary))

# BSF by group and treatment

"BSF at Lp299v baseline"
clinical.df %>%
  filter((Treatment == "PreTrial_1" & TestingOrder == "Lp299v - Placebo") |
           (Treatment == "PreTrial_2" & TestingOrder == "Placebo - Lp299v")) %>%
  with(tapply(BSF, TestingOrder, summary)) #Pre-Lp299v

"BSF at placebo baseline"
clinical.df %>%
  filter((Treatment == "PreTrial_2" & TestingOrder == "Lp299v - Placebo") |
           (Treatment == "PreTrial_1" & TestingOrder == "Placebo - Lp299v")) %>%
  with(tapply(BSF, TestingOrder, summary)) #Pre-Placebo

"BSF after Lp299v treatment"
clinical.df %>%
  filter(Treatment == "Lp299v") %>%
  with(tapply(BSF, TestingOrder, summary)) #Lp299v

"BSF after placebo treatment"
clinical.df %>%
  filter(Treatment == "Placebo") %>%
  with(tapply(BSF, TestingOrder, summary)) #Placebo

# Diversity

"Shannons H at Lp299v baseline"
genus_alpha.df %>%
  filter((Treatment == "PreTrial_1" & TestingOrder == "Lp299v - Placebo") |
           (Treatment == "PreTrial_2" & TestingOrder == "Placebo - Lp299v")) %>%
  with(tapply(ShannonH, TestingOrder, summary)) #Pre-Lp299v

"Shannons H at placebo baseline"
genus_alpha.df %>%
  filter((Treatment == "PreTrial_2" & TestingOrder == "Lp299v - Placebo") |
           (Treatment == "PreTrial_1" & TestingOrder == "Placebo - Lp299v")) %>%
  with(tapply(ShannonH, TestingOrder, summary)) #Pre-Placebo

"Shannons H after Lp299v treatment"
genus_alpha.df %>%
  filter(Treatment == "Lp299v") %>%
  with(tapply(ShannonH, TestingOrder, summary)) #Lp299v

"Shannons H after placebo treatment"
genus_alpha.df %>%
  filter(Treatment == "Placebo") %>%
  with(tapply(ShannonH, TestingOrder, summary)) #Placebo

"Lactiplantibacillus in gut during treatment - Pre-Lp299v Shannons H"
genus_alpha.df %>%
  filter((Treatment == "PreTrial_1" & TestingOrder == "Lp299v - Placebo") |
           (Treatment == "PreTrial_2" & TestingOrder == "Placebo - Lp299v")) %>%
  with(tapply(ShannonH, Lp_present_treat, summary)) #Pre-Lp299v, based on active status

"Lactiplantibacillus in gut during treatment - Post-Lp299v Shannons H"
genus_alpha.df %>%
  filter(Treatment == "Lp299v") %>%
  with(tapply(ShannonH, Lp_present_treat, summary)) #Lp299v, based on active status

# F/B Ratio

fb.df <- phylumProp_appended.df %>%
  mutate(FBratio = Firmicutes / Bacteroidota)

"F/B ratio at Lp299v baseline"
fb.df %>%
  filter((Treatment == "PreTrial_1" & TestingOrder == "Lp299v - Placebo") |
           (Treatment == "PreTrial_2" & TestingOrder == "Placebo - Lp299v")) %>%
  with(tapply(FBratio, TestingOrder, summary)) #Pre-Lp299v

"F/B ratio at placebo baseline"
fb.df %>%
  filter((Treatment == "PreTrial_2" & TestingOrder == "Lp299v - Placebo") |
           (Treatment == "PreTrial_1" & TestingOrder == "Placebo - Lp299v")) %>%
  with(tapply(FBratio, TestingOrder, summary)) #Pre-Placebo

"F/B ratio after Lp299v treatment"
fb.df %>%
  filter(Treatment == "Lp299v") %>%
  with(tapply(FBratio, TestingOrder, summary)) #Lp299v

"F/B ratio after placebo treatment"
fb.df %>%
  filter(Treatment == "Placebo") %>%
  with(tapply(FBratio, TestingOrder, summary)) #Placebo

### P-values: t-tests ###
"Difference between placebo and Lp299v - t-test, Shannons H"
genus_alpha.df %>%
  filter(Treatment %in% c("Placebo", "Lp299v")) %>%
  with(t.test(ShannonH ~ Treatment)) # based on treatment

"Difference between baselines - t-test, Shannons H"
genus_alpha.df %>%
  filter(Treatment %in% c("PreTrial_1", "PreTrial_2")) %>%
  with(t.test(ShannonH ~ Treatment)) # based on baselines

"Lactiplantibacillus in gut during treatment - Pre-Lp299v t-test, Shannons H"
genus_alpha.df %>%
  filter((Treatment == "PreTrial_1" & TestingOrder == "Lp299v - Placebo") |
           (Treatment == "PreTrial_2" & TestingOrder == "Placebo - Lp299v")) %>%
  with(t.test(ShannonH ~ Lp_present_treat)) # based on baselines to Lp299v status

"Lactiplantibacillus in gut during treatment - Post-Lp299v t-test, Shannons H"
genus_alpha.df %>%
  filter(Treatment == "Lp299v") %>%
  with(t.test(ShannonH ~ Lp_present_treat)) # based on Lp299v status

### P-values: RDA ###

"Lactiplantibacillus in gut during treatment - Pre-Lp299v RDA"

clinical_Lp_baseline.df <- filter(clinical.df, (Treatment == "PreTrial_1" & TestingOrder == "Lp299v - Placebo") |
                                    (Treatment == "PreTrial_2" & TestingOrder == "Placebo - Lp299v"))

genusProp_appended_Lp_baseline_pca.df <- genusCount_appended.df %>%
  filter((Treatment == "PreTrial_1" & TestingOrder == "Lp299v - Placebo") |
           (Treatment == "PreTrial_2" & TestingOrder == "Placebo - Lp299v")) %>%
  select(-c(colnames(clinical.df[,-1]),
            Lp_present)) %>%
  select_if(function(col) max(col) != 0) %>%
  select(!ends_with("Lactiplantibacillus")) %>%
  countToProp() #converting back to prop table after eliminating L. plantarum just to control for ANY extraneous impact of including Lp in calculations

lp_baseline.rda <- rda(genusProp_appended_Lp_baseline_pca.df ~ Lp_present_treat,
                       scale = TRUE,
                       data = clinical_Lp_baseline.df,
                       na.rm = TRUE)

anova(lp_baseline.rda, permutations = 9999)

"Lactiplantibacillus in gut during treatment - Post-Lp299v RDA"

clinical_Lp299v.df <- filter(clinical.df, Treatment == "Lp299v")

genusProp_appended_Lp299v_pca.df <- genusCount_appended.df %>%
  filter(Treatment == "Lp299v") %>%
  select(-c(colnames(clinical.df[,-1]),
            Lp_present)) %>%
  select_if(function(col) max(col) != 0) %>%
  select(!ends_with("Lactiplantibacillus")) %>%
  countToProp() #converting back to prop table after eliminating L. plantarum just to control for ANY extraneous impact of including Lp in calculations

lp299v.rda <- rda(genusProp_appended_Lp299v_pca.df ~ Lp_present,
                  scale = TRUE,
                  data = clinical_Lp299v.df,
                  na.rm = TRUE)

anova(lp299v.rda, permutations = 9999)' | R --save > output/statistics.txt



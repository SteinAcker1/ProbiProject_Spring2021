# Probi ProGastro17 Analysis

## Requirements
R v4.1.0

## Instructions for use
1) Ensure you have the right R version on your computer (4.1.0)
2) Ensure you have access to the internet during analysis
3) Ensure your file tree matches the one below (reference taxa may be downloaded with 'wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_wSpecies_train_set.fa.gz', participant data must be accessed from Probi)
4) Navigate to the directory you are using for analysis
5) Enter ./executableReadme.sh
6) Done!

## GitHub Link

https://github.com/SteinAcker1/ProbiProject_Spring2021

## Required File Tree

        .
        ├── data
        │   ├── *participant FASTQ data*
        │   ├── probiBSF.csv
        │   └── probiDemographics.csv
        ├── executableReadme.sh
        └── refTaxa
            └── silva_nr99_v138.1_wSpecies_train_set.fa.gz
    
## Resulting File Tree
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
            └── silva_nr99_v138.1_wSpecies_train_set.fa.gz

## Scripts Guide

packageSetup.R: Downloads all required packages

funcLibrary.R: Defines a variety of functions to be used by the program

script_hpc.R: DADA2 pipeline for taxon assignment (should be run on a very powerful computer

script_local.R: Data wrangling and cleanup (may be run on local computer)

paperScript.R: Performs all analyses used in paper

analysis_aggregate.R and analysis_indiv.R: Scratchpad scripts

# Probi ProGastro17 Analysis

## Introduction
This repository contains all the code required to replicate the analyses performed in my paper, *Assessing the persistence and impact of* Lactiplantibacillus
plantarum *probiotic strain Lp299v in the human digestive tract*. The following packages are used:

- BiocManager (v. 3.13): a helper package for installing other Bioconductor packages.

- DADA2 (v. 1.20.0): a pipeline for cleaning up and assigning taxonomy to 16S rRNA FASTQ amplicon data.

- tidyverse (v. 1.3.1): a suite of packages for data handling and plotting in R.

- vegan (v. 2.5.7): a community ecology package containing a wide variety of helpful statistical tools.

- ggrepel (v. 0.9.1): a package which allows the user to add attractive labels to ggplot2 objects.

- viridis (v. 0.6.1): a set of ggplot2 color palettes which minimize visual bias and are accessible to individuals with a variety of forms of colorblindness while maximizing contrast.

- tidyMicro (v. 1.47): a suite of microbiome assessment tools which work with tidyverse packages.

- latex2exp (v. 0.5.0): a package which allows the user to add LaTeX text to ggplot2 objects.

## Workflow

1) packageSetup.R: Sets up all required packages.

2) funcLibrary.R: Defines a variety of functions to be used by the program.

3) script_hpc.R: DADA2 pipeline for taxon assignment (should be run on a very powerful computer).

4) script_local.R: Data wrangling and cleanup (may be run on local computer).

5) paperScript.R: Performs all analyses used in paper.

## Instructions for use
1) Ensure you have the right R version on your computer (4.1.0). If you have conda installed, this can be easily done with:

        conda create -n myEnv -c conda-forge r-base=4.1.0
        conda activate myEnv

2) Ensure you have access to the internet during analysis, as some packages will be installed from the Internet.
3) Ensure your file tree matches the one below (reference taxa may be downloaded with the command below, participant data must be accessed from Probi)

        wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_wSpecies_train_set.fa.gz
        
4) Navigate to the directory you are using for analysis.
5) Execute packageSetup.R with the following command to set up the required packages:
        
        Rscript R/packageSetup.R
        
6) Execute the analyses using the following command:

        nohup Rscript R/paperScript.R &

7) Done!

*Or*, download executableReadme.sh and follow the instructions inside if you wish to only use one script.

## Required File Tree

        .
        ├── data
        │   ├── *participant FASTQ data*
        │   ├── probiBSF.csv
        │   └── probiDemographics.csv
        ├── R
        │   ├── packageSetup.R
        │   ├── funcLibrary.R
        │   ├── script_hpc.R
        │   ├── script_local.R
        │   └── paperScript.R
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
        ├── R
        │   ├── packageSetup.R
        │   ├── funcLibrary.R
        │   ├── script_hpc.R
        │   ├── script_local.R
        │   └── paperScript.R
        └── refTaxa
            └── silva_nr99_v138.1_wSpecies_train_set.fa.gz

## GitHub Link

https://github.com/SteinAcker1/ProbiProject_Spring2021

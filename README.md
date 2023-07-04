[![DOI](https://zenodo.org/badge/659820253.svg)](https://zenodo.org/badge/latestdoi/659820253)
# Biotic interactions outweigh abiotic factors as drivers of bark microbial communities in Central European forests 

![Important effects on beta diversity](https://github.com/LukDrey/bark_microbiome_drivers/blob/main/Figure4_effects_beta.png)

This repository contains the data and pre-processing pipeline, as well as the R script for the analysis accompanying our paper: 

> Dreyling L, Penone C, Schenk N, Schmitt I, Dal Grande F 2023. Biotic interactions outweigh abiotic factors as drivers of bark microbial communities in Central European forests.
> You can cite the dataset and code with the following [![DOI](https://zenodo.org/badge/659820253.svg)](https://zenodo.org/badge/latestdoi/659820253).

## Contacts

**Lukas Dreyling**  
[E-Mail](mailto:lukas.dreyling@senckenberg.de)  

**Caterina Penone**  
[E-Mail](mailto:caterina.penone@unibe.ch)

**Noëlle Schenk**  
[E-Mail](mailto:noelle.schenk@unibe.ch)

**Imke Schmitt**  
[E-Mail](mailto:imke.schmitt@senckenberg.de)  

**Francesco Dal Grande**  
[E-Mail](mailto:francesco.dalgrande@unipd.it)  

## Contents

1. [Data](Data)
2. [Pre-Processing Pipeline](00_Read_Processing_pipeline.txt)
3. [ASV Curation](01_ASV_Curation.R)
4. [Data Cleaning](02_Data_Cleaning.R)
5. [Alpha Diversity Analysis](03_Alpha_Diversity.R)
6. [Beta Diversity Analysis](04_Beta_Diversity.R)

Additionally you can download the raw reads [here](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA932736).  

If you lack the computing power to process the raw reads, the resulting ASV tables, FASTA files, taxonomy table and metadata are located [here](Data).  

## Before starting

### You will need to have the following software installed.

#### Pre-Processing 
* fastQC http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
* Cutadapt https://cutadapt.readthedocs.io/
* R https://www.r-project.org/
    - dada2 https://benjjneb.github.io/dada2/index.html
    - ShortRead https://kasperdanielhansen.github.io/genbioconductor/html/ShortRead.html
    - Biostrings https://bioconductor.org/packages/release/bioc/html/Biostrings.html
* BLASTn https://www.ncbi.nlm.nih.gov/books/NBK279690/

#### Analysis
* R https://www.r-project.org/
* Rstudio https://www.rstudio.com/
  - here https://here.r-lib.org/
  - tidyverse https://www.tidyverse.org/
  - decontam https://benjjneb.github.io/decontam/
  - phyloseq https://joey711.github.io/phyloseq/
  - LULU https://github.com/tobiasgf/lulu
  - Biostrings https://bioconductor.org/packages/release/bioc/html/Biostrings.html
  - microbiome https://github.com/microbiome/microbiome
  - fantaxtic https://github.com/gmteunisse/Fantaxtic
  - vegan https://rdrr.io/cran/vegan/man/vegan-package.html
  - hillR https://github.com/daijiang/hillR
  - gdm https://github.com/fitzLab-AL/gdm
  - effects https://cran.r-project.org/web/packages/effects/index.html
  - ggeffects https://strengejacke.github.io/ggeffects/index.html
  - https://ggobi.github.io/ggally/

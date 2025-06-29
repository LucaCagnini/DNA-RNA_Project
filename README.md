# DNA-RNA_Project

## Table of contents
- [Overview](#overview)
- [Requirements](#requirements)
- [Workflow](#workflow)
  - [1. Data Preparation and Import](#1-data-preparation-and-import)
  - [2. Signal Extraction](#2-signal-extraction)
  - [3. Quality Control](#3-quality-control)
  - [4. Beta and M Value Calculation](#4-beta-and-m-value-calculation)
  - [5. Normalization and Batch Effect Analysis](#5-normalization-and-batch-effect-analysis)
  - [6. Statistical Analysis](#6-statistical-analysis)
  - [7. Visualization](#7-visualization)
- [Contacts](#contacts)

## Overview

This repository contains the final project of Group 5 for the DNA/RNA Dynamics course (MSc in Bioinformatics, University of Bologna, a.y. 2024/2025). It features a complete pipeline for analyzing Illumina 450K methylation data using R. The pipeline includes preprocessing (with the preprocessNoob method), quality control, normalization, statistical analysis (using  Mann-Whitney U test and P-value threshold of 0.05), principal component analysis (PCA), and the identification of differentially methylated positions (DMPs) between control (CTRL) and disease (DIS) samples.

## Requirements
The project was performed using R and Rstudio, analysing data from the platform Illumina HumanMethylation450k (input_data). Necessary packages were installed in our R enviroment:

```r
install.packages(minfi)
install.packages(ggplot2)
install.packages(nitr)
install.packages(BiocManager)
install.packages(factoextra)
install.packages(cluster)
install.packages(qqman)
install.packages(gplots)
```


# Project Workflow

This document outlines the step-by-step workflow of our DNA/RNA methylation project, comparing CTRL and DIS sample groups using Illumina microarray data.

---

## 1. Data Preparation and Import
- Install and load required R packages: `minfi`, `BiocManager`, `knitr`.
- Clean the R environment and set the working directory.
- Load the sample sheet using `read.metharray.sheet()` to import metadata.
- Read raw data using `read.metharray.exp()` to generate the `RGset` object.

---

## 2. Signal Extraction
- Extract fluorescence intensity data for Red (Cy5) and Green (Cy3) channels from `RGset` using `getRed()` and `getGreen()`.
- Store the data into two separate dataframes: `Red` and `Green`.

---

## 3. Quality Control
- Classify sample quality based on the percentage of failed probes:
  - **High quality**: < 0.01%
  - **Good quality**: < 0.2%
  - **Low quality**: > 0.2%
  - **Critical quality**: around 1% (may require exclusion)

---

## 4. Beta and M Value Calculation
- Split samples into CTRL and DIS groups using metadata.
- Create `MSet.raw` objects for each group.
- Compute **Beta values** with `getBeta()` and **M values** with  `getM()`.
- Calculate and plot mean methylation values for both groups.

---
## 5. Normalization and Batch Effect Analysis
- Perform Principal Component Analysis (PCA) to assess batch effects (e.g., by `Sentrix_ID`).
- PCA suggested that normalization did not fully correct batch effects.

---

## 6. Statistical Analysis
- Perform t-tests for each probe comparing CTRL vs DIS.
- Create a dataframe containing p-values.
- Filter probes with p ≤ 0.05 for significance.
- Plot the distribution of p-values.

---

## 7. Visualization
- Generate a **heatmap** using the top 100 most significant probes.
- Apply **hierarchical clustering** with various linkage methods:
  - Single
  - Complete
  - Average (found to be most biologically interpretable)
  
  ---

  
## Contacts

Project members:

-Marco Cuscunà (marco.cuscuna@studio.unibo.it)
-Marco Centenaro (marco.centenaro@studio.unibo.it)
-Michele Carbonieri (michele.carbonieri@studio.unibo.it)
-Marina Mariano (marina.mariano@studio.unibo.it)
-Luca Cagnini (luca.cagnini@studio.unibo.it)
-Massimo Lanari (massimo.lanari@studio.unibo.it)




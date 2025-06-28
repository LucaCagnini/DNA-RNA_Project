# DNA-RNA_Project

## Table of contents
- [Overview](#Overview)
- [Requirements](#requirements)
- [Workflow](#Workflow)
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


## Workflow




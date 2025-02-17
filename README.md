# Overview

This repository contains the source c code for the paper Detection of Cell-type-specific Differentially Methylated Regions
in Epigenome-Wide Association Studies. The code will later be packed into a R function

# Introduction

DNA methylation at cytosine-phosphate-guanine (CpG) sites is one of the most important epigenetic markers. Therefore, epidemiologists are interested in investigating DNA methylation in large cohorts through epigenome-wide association studies (EWAS). However, the observed EWAS data are bulk data with signals aggregated from distinct cell types. Deconvolution of cell-type-specific signals from EWAS data is challenging because phenotypes can affect both cell-type proportions and cell-type-specific methylation levels. Recently, there has been active research on detecting cell-type-specific risk CpG sites for EWAS data. However, since existing methods all assume that the methylation levels of different CpG sites are independent and perform association detection for each CpG site separately, although they significantly improve the detection at the aggregated-level---identifying a CpG site as a risk CpG site as long as it is associated with the phenotype in any cell type, they have low power in detecting cell-type-specific associations for EWAS with typical sample sizes. Here, we develop a new method, Fine-scale inference for Differentially Methylated Regions (FineDMR), to borrow strengths of nearby CpG sites to improve the cell-type-specific association detection. Via a Bayesian hierarchical model built upon Gaussian process functional regression, FineDMR takes advantage of the spatial dependencies between CpG sites. FineDMR can provide cell-type-specific association detection as well as output subject-specific and cell-type-specific methylation profiles for each subject. 

# System Requirements

Before using the FineDMR function, users should install the [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/). 

The following packages are also needed in order to use FineDMR: splines, quadprog

# Input

The input should include:

$X$: A $n$ by $Q$ matrix where $n$ is the number of sample size and Q is the number of covariates

$O$: A $n$ by $G$ matrix where $G$ is the number of CpG sites

$t$: A vector of length $G$ containing the genomic locations of the CpG sites

$K$: The number of cell types

$FWER$: the FWER 

# Output

The output includes:

$beta$: A $K$ by $Q$ by $G$ array containing the estimated effect size

$pvalue$: A $K$ by $Q$ by $G$ array containing the p-values estimated from the model

$prop$: The estimated cellular proportions

$M$: An array with dim $K$ by $G$ by $n$ containing the cell-type-specific profiles 

$DMRs$: A list with length $K$ containing the DMRs detected for each cell type


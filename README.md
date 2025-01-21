# Overview

This repository contains the source c code for the paper Cell-type-specific DMR Detection]{Detection of Cell-type-specific Differentially Methylated Regions
in Epigenome-Wide Association Studies. The code will later be packed into a R function

# Introduction

Epigenome-wide association studies (EWAS) aim to identify the
associations between the DNA methylation levels of (CpG) sites
and phenotypes of interest. EWAS data is aggregated from mixed cell types. It is of importance to do cell type deconvolution and detect the associations at the cell-type-specific levels. This paper aims to find the cell-type-specific risk CpG sites by leverage the information across nearby CpG sites. As a result, this proposed method can detect the Differentially Methylated Regions (DMRs).

# Input

The input should include:

$X$: A $n$ by $Q$ matrix where $n$ is the number of sample size and Q is the number of covariates

$O$: A $n$ by $G$ matrix where $G$ is the number of CpG sites

$t$: A vector of length $G$ containing the genomic locations of the CpG sites

$nb$: The number of blocks

$block_{vec}$: A vector of length $nb$ containing the number of CpG sites in each block

$df$: A vector of length $nb$ containing the degree of freedom of B-splines in each block

$K$: The number of cell types

# Output

The output includes:

$pvalue$: A $K$ by $Q$ by $G$ array containing the p-values estimated from the model

$prop$: The estimated cellular proportions



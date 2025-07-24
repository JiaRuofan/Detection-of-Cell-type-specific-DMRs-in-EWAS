# Overview

This repository contains the source c code for the paper Detection of Cell-type-specific Differentially Methylated Regions
in Epigenome-Wide Association Studies. The code will later be packed into a R function

# Introduction

DNA methylation at cytosine-phosphate-guanine (CpG) sites is one of the most important epigenetic markers. Therefore, epidemiologists are interested in investigating DNA methylation in large cohorts through epigenome-wide association studies (EWAS). However, the observed EWAS data are bulk data with signals aggregated from distinct cell types. Deconvolution of cell-type-specific signals from EWAS data is challenging because phenotypes can affect both cell-type proportions and cell-type-specific methylation levels. Recently, there has been active research on detecting cell-type-specific risk CpG sites for EWAS data. However, since existing methods all assume that the methylation levels of different CpG sites are independent and perform association detection for each CpG site separately, although they significantly improve the detection at the aggregated-level---identifying a CpG site as a risk CpG site as long as it is associated with the phenotype in any cell type, they have low power in detecting cell-type-specific associations for EWAS with typical sample sizes. Here, we develop a new method, Fine-scale inference for Differentially Methylated Regions (FineDMR), to borrow strengths of nearby CpG sites to improve the cell-type-specific association detection. Via a Bayesian hierarchical model built upon Gaussian process functional regression, FineDMR takes advantage of the spatial dependencies between CpG sites. FineDMR can provide cell-type-specific association detection as well as output subject-specific and cell-type-specific methylation profiles for each subject. 

# Maintainer

Ruofan Jia 1155165249@link.cuhk.edu.hk

# System Requirements

Before using the FineDMR function, users should install the [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/). 

The following R packages are also needed in order to use FineDMR: splines, quadprog

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

$profile$: An array with dim $K$ by $G$ by $n$ containing the cell-type-specific profiles 

$DMRs$: A list with length $K$ containing the DMRs detected for each cell type

# Data availability
One can access the simulation data via the [Onedrive](https://mycuhk-my.sharepoint.com/:f:/g/personal/1155165249_link_cuhk_edu_hk/EpoAC32GLP1Hm9lQNvxaCe8BzHwd5Si9N64tHmTOhdvQmA?e=gA6nrT) with PW FineDMR2025.

# Example

Below is an example showing how to use FineDMR to analyze the data in the alternative case of the simulation study when $n=300$.

```
#Input: sample size n=300; number of CpG sites G=10000; number of covariates Q=1
load('O.Rdata')#The observed methylation data
load('X.Rdata')#The phenotypes
load('t.Rdata')#Genomic locations

# Load required dependencies
library(quadprog)
library(splines)

source('FineDMR.R')

#Run the FineDMR algorithm
FD_res<-FineDMR(O,X,t,K=6) 

#Adjust the p values for the covariates X with Bonferroni method
p_array_adjust<-padjust(FD_res$`p values`[,2,],method='bonferroni')

#Detect the DMRs in cell type 1-6
K<-6
FWER<-0.01
for (k in 1:K){

  print(DMR_detect(p_array_adjust[k,]<=FWER))
}


```
In order to get the heatmap shown as following, one can run the following codes with ggplots and ggsci required:

 ![heatmap](https://github.com/JiaRuofan/Detection-of-Cell-type-specific-DMRs-in-EWAS/blob/main/simu_heatmap.png?raw=true)

```
library(ggplot2)
library(ggsci)

G<-dim(O)[2]

FineDMR_heatmap(G,K,FD_res$`beta`[,2,])

```
To get the manhatten plot shown as following, one can run the following codes:

 ![manhatten](https://github.com/JiaRuofan/Detection-of-Cell-type-specific-DMRs-in-EWAS/blob/main/simu_manhatten.png?raw=true)

```
library(ggplot2)
library(ggsci)
#Draw a manhatten plot for cell type 1
#P is a sequence of pvalue, Location is the genomic location, clas corresponds to the significance of P
Q<-2
P<-FD_res$`p values`[1,2,]
Location<-t/0.01#de-normalization
clas<-p_array_adjust[1,]
FineDMR_manhatten(P,Location,clas,FWER,G,K,Q)

```
Below shows the scatter plots of subject-specific and cell-type-specific methylation profiles of CpG sites from genomic location 3148005 to 3195414 in cell type 1.

 ![scatter](https://github.com/JiaRuofan/Detection-of-Cell-type-specific-DMRs-in-EWAS/blob/main/simu_scatterplot_smooth.png?raw=true)

You can run the following code to derive the plot:
```
n<-300
startp<-4030
endp<-4120
bl<-endp-startp+1

rt<-rep(t[startp:endp],n)#71
for (i in 1:n){
  rProf<-c(rProf,as.vector(FD_res$`profile`[1,starp:endp,i]))
}

rX<-c()

for (i in 1:n){
  rX<-c(rX,rep(X[i,2],bl))
}


dat<-data.frame(rt,X=as.factor(rX),rProf)

ggplot(dat,aes(x=rt,y=rProf,group=X,color=X))+scale_color_npg()+ 
  geom_point(size=1)+scale_color_manual(values = c('#4DBBD5FF','#E64B35FF'))+
  geom_smooth(se=FALSE)+ 
  theme_bw() + 
  theme(legend.position='top', 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )+theme(axis.text=element_text(size=20,face = "bold"),axis.title=element_text(size=20,face = "bold"),
          legend.text=element_text(size=20),
          legend.title=element_text(size=20))+labs(fill='',x = 'Genomic Location', y = 'Cell-type-specific Profile')
```
And you can also plot the corresponding effect size of $X$.

 ![effect](https://github.com/JiaRuofan/Detection-of-Cell-type-specific-DMRs-in-EWAS/blob/main/simu_effect%20size.png?raw=true)

You can run the following code to derive the plot:

```
rt<-BP[startp:endp]
rmean=FD_res$`beta`[1,2,startp:endp]
clas<-c()

for (i in 1:bl)
if (p_array_adjust[1,(startp:endp)[i]]<=0.01){
  clas[i]<-'risk'
}else{
  clas[i]<-'non-risk'
}

dat<-data.frame(t=rt,rmean,clas)


ggplot(data = dat, mapping = aes(x = t)) +scale_color_npg()+scale_color_manual(values = c('#F39B7FFF','#91D1C2FF'))+
  geom_line(mapping = aes(y = rmean, color = clas, group =1),linewidth=1.5)+ 
  theme_bw() + 
  theme(legend.position='top', 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )+theme(axis.text=element_text(size=20,face = "bold"),axis.title=element_text(size=20,face = "bold"),
          legend.text=element_text(size=20),
          legend.title=element_text(size=20))+labs(x = 'Genomic Location', y = 'Effect Size')+ guides(color=guide_legend(title = ""))

```

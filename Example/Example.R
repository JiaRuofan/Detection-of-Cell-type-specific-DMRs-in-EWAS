library(quadprog)
library(splines)

source('FineDMR.R')

#sample size n=100; number of CpG sites G=1000; number of covariates Q=1

O<-read.table('O.txt')#A n by G matrix with G observed DNA methylation levels for n samples
X<-read.table('X.txt')#A n by Q matrix with covariates for n samples 
t<-scan('t.txt') #genomic locations with length of G

#########Run the FineDMR algorithm
FD_res<-FineDMR(O,X,t,K=6) 

#########Adjust the p values for the covariates X with Bonferroni method
p_array_adjust<-padjust(FD_res$`p values`[,2,],method='bonferroni')

#########Detect the DMRs in cell type 1-6
K<-6

for (k in 1:K){

  print(DMR_detect(p_array_adjust[k,]<=0.01))
}

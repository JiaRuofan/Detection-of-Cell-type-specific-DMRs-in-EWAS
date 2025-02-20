library(quadprog)
library(splines)

source('FineDMR.R')

O<-read.table('O.txt')
X<-read.table('X.txt')


#########Run the FineDMR algorithm
FD_res<-FineDMR(O,X,tj,K=6) 

#########Adjust the p values for the covariates X with Bonferroni method
p_array_adjust<-padjust(FD_res$`p values`[,2,],method='bonferroni')

#########Detect the DMRs for cell type 1-6
K<-6

for (k in 1:K){

  print(DMR_detect(p_array_adjust[k,]<=0.01))
}

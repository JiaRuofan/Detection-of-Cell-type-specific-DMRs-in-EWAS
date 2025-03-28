n<-1000
G<-10000
Q<-2
K<-6
x_vec<-scan('D:\\HP_simu\\new simu\\semi\\throu\\set23\\n1000\\X.txt')


nblock<-134
block_list<-scan('D:\\HP_simu\\new simu\\semi\\throu\\set23\\n1000\\block.txt')

O<-matrix(0,nrow=n,ncol=G)
X<-matrix(0,nrow=n,ncol=Q)

nt<-0
for (b in 1:nblock){
  O_vec<-scan(paste0('D:\\HP_simu\\new simu\\semi\\throu\\set23\\n1000\\o\\O',b,'.txt'))
  
  for (j in 1:block_list[b]){
    nt<-nt+1
    for (i in 1:n){
      O[i,nt]<-O_vec[(i-1)*block_list[b]+j]
    }
  }
}


for (i in 1:n){
  for (j in 1:Q){
    X[i,j]<-x_vec[(i-1)*Q+j]
  }
}





re_semi_set23_n1000<-HIRE(t(O), t(X[,-1]), 6, tol=10^(-5), num_iter=2000, alpha=0.01)






K <- 6
outT <- myRefFreeCellMix(t(O), mu0=myRefFreeCellMixInitialize(t(O), K = K))
estProp_RF <- outT$Omega
colnames(estProp_RF)<-c('cell1','cell2','cell3','cell4','cell5','cell6')

design<-data.frame(disease=X[,2])

Design_out <- makeDesign(design, estProp_RF)






fitted_model <- fitModel(Design_out, t(O))

res_table1 <- csTest(fitted_model, coef = "disease", 
                     cell_type = "cell1", contrast_matrix = NULL)

res_table2 <- csTest(fitted_model, coef = "disease", 
                     cell_type = "cell2", contrast_matrix = NULL)

res_table3 <- csTest(fitted_model, coef = "disease", 
                     cell_type = "cell3", contrast_matrix = NULL)

res_table4 <- csTest(fitted_model, coef = "disease", 
                     cell_type = "cell4", contrast_matrix = NULL)

res_table5 <- csTest(fitted_model, coef = "disease", 
                     cell_type = "cell5", contrast_matrix = NULL)

res_table6 <- csTest(fitted_model, coef = "disease", 
                     cell_type = "cell6", contrast_matrix = NULL)






#####################




#6-1;5-6;4-2;3-5;2-4;1-3

design<-data.frame(disease=X[,2])


ox<-t(O)
rownames(ox)<-1:G
colnames(ox)<-1:n
rownames(estProp_RF)<-1:n
rownames(design)<-1:n

library(TCA)
tca.mdl <- tca(X = ox,
               W = estProp_RF,
               C1 = design)



####################

celldmc.o <- CellDMC(ox, X[,2], estProp_RF)


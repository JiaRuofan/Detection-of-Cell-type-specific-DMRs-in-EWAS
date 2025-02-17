library(quadprog)


#generate samples from a Dirichlet distribution
rDirichlet <- function(alpha_vec){
  num <- length(alpha_vec)
  temp <- rgamma(num, shape = alpha_vec, rate = 1)
  return(temp / sum(temp))
}


#CorDescent
CorDescent <- function(MethMatr, num_celltype, tol = 0.01, showIter = FALSE){
  err0 <- 0
  err1 <- 1000
  m <- nrow(MethMatr) #number of CpG sites
  n <- ncol(MethMatr) #number of samples
  K <- num_celltype
  if(m < n){
    stop("The CpG site number must be larger than the sample number!")
  }
  #initialize P_matr
  P_matr_t <- vapply(seq_len(n), function(i){ rDirichlet(rep(2,K)) }, FUN.VALUE=rep(-1,K))
  while(abs(err1 - err0) >= tol){
    err0 <- err1
    #update U_matr
    Dmat <- 2*P_matr_t%*%t(P_matr_t)
    Amat <- cbind(diag(rep(1,K)), diag(rep(-1,K)))
    bvec <- c(rep(0,K), rep(-1,K))
    U_matr_t <- t( vapply(seq_len(m), function(j){
      dvec <- 2*P_matr_t %*% as.numeric(MethMatr[j,])
      
      solu <- solve.QP(Dmat, dvec, Amat, bvec)
      solu$solution
    }, FUN.VALUE=rep(-1,K)) )
    
    #update P_matr
    Dmat <- 2*t(U_matr_t) %*% U_matr_t
    Amat <- cbind(matrix(1, K, K), diag(rep(1,K)))
    bvec <- c(rep(1, K), rep(0, K))
    P_matr_t <- vapply(seq_len(n), function(i){
      dvec <- 2 * t(U_matr_t) %*% as.numeric(MethMatr[ ,i])
      solu <- solve.QP(Dmat, dvec, Amat, bvec, meq = K)
      solu$solution 
    }, FUN.VALUE=rep(-1,K))
    
    #calculate errors
    err1 <- sum((MethMatr - U_matr_t %*% P_matr_t)^2)
    if(showIter == TRUE){
      message("  ", err1, "\n")
    }
  }
  
  return(list(U=U_matr_t, P=P_matr_t))	
}

#split the CpG sites into blocks
psearch<-function(cop_list,min.block=50,max.block=200,methy_loc){
  
  if ((length(methy_loc)-cop_list[length(cop_list)])<=max.block){
    return(cop_list)
  }else{
    
    nextp<-which.max(diff(methy_loc[(cop_list[length(cop_list)]+min.block):(cop_list[length(cop_list)]+max.block)]))
    return(psearch(c(cop_list,cop_list[length(cop_list)]+min.block+nextp),min.block,max.block,methy_loc))
    
  }
  
}



FineDMR<-function(O_mat, X, t, df_prop=1/5, K){
  
  
  block_list<-diff(c(psearch(1,50,200,as.numeric(t)),length(t)+1))
  nblock<-length(block_list)
  
  nbasis_sub_list<-list()
  B_sub_list<-list()
  BTBIBT_sub_list<-list()
  for (j in 1:(nblock)){
    nbasis_sub<-floor(block_list[j]*df_prop)+3
    nbasis_sub_list[[j]]<-nbasis_sub
    B_sub_list[[j]]<-bs(t[(sum(block_list[1:(j-1)])+1):(sum(block_list[1:(j-1)])+block_list[j])], df=nbasis_sub,degree=3)
    BTB<-t(B_sub_list[[j]])%*%B_sub_list[[j]]
    BTBI<-solve(BTB)
    BTBIBT_sub_list[[j]]<-BTBI%*%t(B_sub_list[[j]])
  }
  
  nbasis<-c()
  for (i in 1:nblock){
    nbasis<-c(nbasis,nbasis_sub_list[[i]])
  }
  
  B_r<-c()
  
  
  for (j in 1:nblock){
    for (l in 1:block_list[j]){
      B_r<-c(B_r,B_sub_list[[j]][l,])
    }
  }
  
  
  BTBIBT_r<-c()
  
  
  for (j in 1:nblock){
    for (l in 1:nbasis_sub_list[[j]]){
      BTBIBT_r<-c(BTBIBT_r,BTBIBT_sub_list[[j]][l,])
    }
  }
  
  P_t<-t(CorDescent(t(O_mat),K)[[2]])
  
  X<-cbind(rep(1,dim(as.data.frame(X))[1]),X)
  
  system("R CMD SHLIB FineDMR.c -lgsl -lgslcblas")
  dyn.load("FineDMR.so")
  
  #B_r;BTBIBT_r
  FineDMRRcallC <- function(O_mat, X, P_t, t, nblock, block_list, nbasis, K, n_batch=50, num_iter=1, sigma2_init=0.01, w_init=20, v_init=0.0001){
    args <- list("P_init"=as.numeric(P_t), "Omat_r"=as.numeric(O_mat), "B_r"=as.numeric(B_r), 'BTBIBT'=as.numeric(BTBIBT_r),
                 "Omat_r_dim"=as.integer(dim(O_mat)), "X_r"=as.numeric(X), "X_r_dim"=as.integer(dim(X)), "num_iter" = as.integer(num_iter),
                 'nblock_r'=as.integer(nblock),'t_r'=as.numeric(t),'block_r'=as.integer(block_list),'nbasis_r'=as.integer(nbasis),'K_r'=as.integer(K),
                 'n_batch_r'=as.integer(n_batch),'sigma2_init'=as.numeric(sigma2_init),'w_init'=as.numeric(w_init),'v_init'=as.numeric(v_init))
    res_list <- .Call("FineDMRRcallC", args)		
    return(res_list)
  }
  
  
  res_list <- FineDMRRcallC(O_mat, X, P_t, t, nblock, block_list, nbasis, K)
  
  return(res_list)
  
}


##################Detect the DMRs

DMR_detect<-function(sequ){
  
  res<-list()
  i<-1
  DMR<-sequ[i]
  while(i<length(sequ)){
    if (sequ[i+1]-sequ[i]>1){
      res<-c(res,list(DMR))
      DMR<-sequ[i+1]
    }else{DMR<-c(DMR,sequ[i+1])}
    i<-i+1
  }
  return(res)
}

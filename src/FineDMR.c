#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_splinalg.h>
#include <gsl/gsl_spblas.h>
#include <time.h>

#include "FineDMR.h"


extern SEXP FineDMRRcallC(SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"FineDMRRcallC",          (DL_FUNC) &FineDMRRcallC,          1},
  {NULL, NULL, 0}
};

void R_init_FineDMR(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}




SEXP getListElement (SEXP list, char *str) {
  
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i;
  
  for(i = 0; i < length(list); i++) {
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0){
      elmt = VECTOR_ELT(list, i);
      break;
    }
  }
  
  return elmt;
}

SEXP FineDMRRcallC(SEXP args) {
  
  //prior
  int nProtected = 0;
  
  
  int nblock;
  
  SEXP list, P_init, Omat_r, Omat_r_dim, X_r, X_r_dim, num_iter, nblock_r, t_r, block_r, nbasis_r, K_r, B_r, BTBIBT, \
  n_batch_r, sigma2_init, w_init, v_init;
  
  list = args;
  P_init = getListElement(list, "P_init");
  Omat_r = getListElement(list, "Omat_r");
  Omat_r_dim = getListElement(list, "Omat_r_dim");
  X_r = getListElement(list, "X_r");
  X_r_dim = getListElement(list, "X_r_dim");
  num_iter = getListElement(list, "num_iter");
  nblock_r = getListElement(list, "nblock_r");
  t_r = getListElement(list, "t_r");
  block_r = getListElement(list, "block_r");
  nbasis_r = getListElement(list, "nbasis_r");
  K_r = getListElement(list, "K_r");
  B_r = getListElement(list, "B_r");
  n_batch_r = getListElement(list, "n_batch_r");
  sigma2_init = getListElement(list, "sigma2_init");
  w_init = getListElement(list, "w_init");
  v_init = getListElement(list, "v_init");
  
  nblock=INTEGER(nblock_r)[0];
  
  
  int *block=(int *) malloc(nblock*sizeof(int));
  int *nbasis=(int *) malloc(nblock*sizeof(int));
  
  for(int j=0; j<nblock; j++){
    block[j] = INTEGER(block_r)[j];
  }
  
  
  for(int j=0; j<nblock; j++){
    nbasis[j] = INTEGER(nbasis_r)[j];
  }
  
  int Q;
  Q=INTEGER(X_r_dim)[1];
  
  int G;
  G=INTEGER(Omat_r_dim)[1];
  
  int n;
  n=INTEGER(Omat_r_dim)[0];
  
  int type;
  type=INTEGER(K_r)[0];
  
  
  int n_batch;
  n_batch=INTEGER(n_batch_r)[0];
  
  
  int maxit;
  maxit=INTEGER(num_iter)[0];
  
  
  double** t=C2darray_diy2(nblock,block);
  double*** O=C3darray_diy3(nblock,n,block);
  double*** B=C3darray_diy23(nblock,block,nbasis);
  double** X=make2Darray(n,Q);
  double** Pest=make2Darray(n,type);
  double*** ThetaHat=C3darray_diy23(nblock,nbasis,block);//C
  
  double Msigma;
  Msigma=REAL(sigma2_init)[0];
  
  
  double Mw;
  Mw=REAL(w_init)[0];
  double Mv;
  Mv=REAL(v_init)[0];
  
  
   int coun=0;
  
  for(int i=0;i<nblock;i++)
  {
    
    for(int j=0;j<block[i];j++)
    {
      
      t[i][j]=REAL(t_r)[coun];
      coun=coun+1;
    }
    
  }
  
  coun=0;
  for (int j=0;j<(nblock);j++){
    
    for(int l=0;l<block[j];l++)
    {
      for(int i=0;i<n;i++)
      {
        O[j][i][l]=REAL(Omat_r)[coun];
        coun=coun+1;
      }
    }
    
    
  }
  
  
  coun=0;
  for (int l=0;l<nblock;l++){
    
    for(int i=0;i<block[l];i++)
    {
      for(int j=0;j<nbasis[l];j++)
      {
        
        B[l][i][j]=REAL(B_r)[coun];
        coun=coun+1;
      }
    }}
  
  
  
  for(int i=0;i<n;i++)
  {
    for(int j=0;j<Q;j++)
    {
      X[i][j]=REAL(X_r)[i*Q+j];
    }
  }
  
  coun=0;
  for(int l=0;l<nblock;l++){
    for(int i=0;i<nbasis[l];i++)
    {
      for(int j=0;j<block[l];j++)
      {
        ThetaHat[l][i][j]=REAL(BTBIBT)[coun];
        coun=coun+1;
      }
      
    }}
  
  
  for(int i=0;i<n;i++)
  {
    for(int j=0;j<type;j++)
    {
      Pest[i][j]=REAL(P_init)[i*type+j];
    }
  }
  
  
  double ***Emu;
  Emu = (double ***)malloc(nblock*sizeof(double **));
  for(int i=0; i<nblock; i++){
    Emu[i] = (double **) malloc(n*sizeof(double *));
    for(int j=0; j<n; j++){
      Emu[i][j] = (double *) malloc(type*block[i]*sizeof(double));
    }
  }
  
  
  double ****beta_est;
  beta_est = (double ****)malloc(nblock*sizeof(double ***));
  for(int i=0; i<nblock; i++){
    beta_est[i] = (double ***) malloc(type*sizeof(double **));
    for(int j=0; j<type; j++){
      beta_est[i][j] = (double **) malloc(block[i]*sizeof(double *));
      for(int l=0; l<block[i]; l++){
        beta_est[i][j][l] = (double *) malloc(Q*sizeof(double));
      }
    }
  }
  
  
  double ****beta_z;
  beta_z = (double ****)malloc(nblock*sizeof(double ***));
  for(int i=0; i<nblock; i++){
    beta_z[i] = (double ***) malloc(type*sizeof(double **));
    for(int j=0; j<type; j++){
      beta_z[i][j] = (double **) malloc(block[i]*sizeof(double *));
      for(int l=0; l<block[i]; l++){
        beta_z[i][j][l] = (double *) malloc(Q*sizeof(double));
      }
    }
  }
  
  
  FineDMR(nblock, block, nbasis, Q, G, n, type, n_batch, maxit, t, O,
          B, X, Msigma, Mw, Mv, Emu, beta_est, beta_z, Pest, ThetaHat);
  
  
  //return values
  SEXP res_P_c, res_P_c_dim;
  PROTECT(res_P_c = allocVector(REALSXP, n*type));
  ++nProtected;
  for(int k=0; k < type; k++){
    for(int i=0; i < n; i++){
      REAL(res_P_c)[k + type*i] = Pest[i][k];
    }
  }
  PROTECT(res_P_c_dim = allocVector(INTSXP, 2));
  ++nProtected;	
  INTEGER(res_P_c_dim)[0] = n;
  INTEGER(res_P_c_dim)[1] = type;
  setAttrib(res_P_c, R_DimSymbol, res_P_c_dim);
  
  
  //return profile
  SEXP res_profile_c, res_profile_c_dim;
  PROTECT(res_profile_c = allocVector(REALSXP, G*n*type));
  ++nProtected;
  
  coun=0;
  for (int k=0;k<type;k++){
    for (int j = 0; j < nblock; j++){
      
      for (int i = 0; i < n; i++){
        for (int l = 0; l < (block[j]);l++){
          REAL(res_profile_c)[coun] = Emu[j][i][k*block[j]+l];
          coun=coun+1;
        }
      }
    }
  }
  
  PROTECT(res_profile_c_dim = allocVector(INTSXP, 3));
  ++nProtected;	
  INTEGER(res_profile_c_dim)[0] = type;
  INTEGER(res_profile_c_dim)[1] = n;
  INTEGER(res_profile_c_dim)[2] = G;
  setAttrib(res_profile_c, R_DimSymbol, res_profile_c_dim);
  
  
  //return beta
  SEXP res_beta_c, res_beta_c_dim;
  PROTECT(res_beta_c = allocVector(REALSXP, G*Q*type));
  ++nProtected;
  
  coun=0;
  for (int ty=0;ty<type;ty++){
    
    for (int q=1;q<Q;q++){
      
      for (int i=0;i<nblock;i++){
        
        for (int j=0;j<block[i];j++){
          
          REAL(res_beta_c)[coun] = beta_est[i][ty][j][q];
          coun=coun+1;
          
        }}}}
  
  
  
  PROTECT(res_beta_c_dim = allocVector(INTSXP, 3));
  ++nProtected;	
  INTEGER(res_beta_c_dim)[0] = type;
  INTEGER(res_beta_c_dim)[1] = n;
  INTEGER(res_beta_c_dim)[2] = G;
  setAttrib(res_beta_c, R_DimSymbol, res_beta_c_dim);
  
  //return z
  SEXP res_z_c, res_z_c_dim;
  PROTECT(res_z_c = allocVector(REALSXP, G*Q*type));
  ++nProtected;
  
  coun=0;
  for (int ty=0;ty<type;ty++){
    
    for (int q=1;q<Q;q++){
      
      for (int i=0;i<nblock;i++){
        
        for (int j=0;j<block[i];j++){
          
          REAL(res_z_c)[coun] = beta_z[i][ty][j][q];
          coun=coun+1;
          
        }}}}
  
  
  
  PROTECT(res_z_c_dim = allocVector(INTSXP, 3));
  ++nProtected;	
  INTEGER(res_z_c_dim)[0] = type;
  INTEGER(res_z_c_dim)[1] = n;
  INTEGER(res_z_c_dim)[2] = G;
  setAttrib(res_z_c, R_DimSymbol, res_z_c_dim);
  
  
  //return list
  SEXP resList; 
  PROTECT(resList = allocVector(VECSXP, 4));
  ++nProtected;
  
  SEXP names; // components names in the return list
  PROTECT(names = allocVector(STRSXP, 4));
  ++nProtected;
  
  SET_STRING_ELT(names, 0, mkChar("Proportions"));
  SET_STRING_ELT(names, 1, mkChar("Profile"));
  SET_STRING_ELT(names, 2, mkChar("beta"));
  SET_STRING_ELT(names, 3, mkChar("z value"));
  
  
  //add elements to the return list
  SET_VECTOR_ELT(resList, 0, res_P_c);
  SET_VECTOR_ELT(resList, 1, res_profile_c);
  SET_VECTOR_ELT(resList, 2, res_beta_c);
  SET_VECTOR_ELT(resList, 3, res_z_c);
  
  
  UNPROTECT(nProtected);
  
  
  return resList;		     							
}

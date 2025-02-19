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

#include "FineDMR2.h"


extern SEXP FineDMRRcallC(SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"FineDMRRcallC", (DL_FUNC) &FineDMRRcallC, 1},
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
  
  SEXP list, P_init, Omat_r, B_r, BTBIBT, Omat_r_dim, X_r, X_r_dim, num_iter, nblock_r, t_r, block_r, nbasis_r, K_r,  \
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
  BTBIBT= getListElement(list, "BTBIBT");
  
  
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
  for(int l=0;l<nblock;l++){
    
    for(int j=0;j<block[l];j++)
    {
      for(int i=0;i<nbasis[l];i++)
      {
        B[l][j][i]=REAL(B_r)[coun];
        coun=coun+1;
      }
    }}
  
  
  
  for(int j=0;j<Q;j++)
  {
    for(int i=0;i<n;i++)
    {
      X[i][j]=REAL(X_r)[j*n+i];
    }
  }
  
  coun=0;
  for(int l=0;l<nblock;l++){
    
    for(int j=0;j<block[l];j++)
    {
      for(int i=0;i<nbasis[l];i++)
      {
        ThetaHat[l][i][j]=REAL(BTBIBT)[coun];
        coun=coun+1;
      }
      
    }}
  
  
  
  for(int j=0;j<type;j++)
  {
    for(int i=0;i<n;i++)
    {
      Pest[i][j]=REAL(P_init)[j*n+i];
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
  
  
  double Ll=0;
  double Ll_old;
  double Ll_pre=100000;
  double alpha_w=0.001;
  double mu_w=1;
  double mu_v=0.01;
  double sigma_v=100;
  
  int batch;
  
  double BLOCK_del; 
  
  double res_old=0;
  double res_new=0;
  int it=0;
  
  double prior_alpha_sigma=0.0001;
  double prior_beta_sigma=1;
  
  double Msigma_old;
  
  double w_old=0;
  double v_old=0;
  
  double var_trans;
  
  double ****PSP_mat;
  PSP_mat = (double ****)malloc(nblock*sizeof(double ***));
  for(int i=0; i<nblock; i++){
    PSP_mat[i] = (double ***) malloc(n*sizeof(double **));
    for(int j=0; j<n; j++){
      PSP_mat[i][j] = (double **) malloc(block[i]*sizeof(double *));
      for(int l=0; l<(block[i]); l++){
        PSP_mat[i][j][l] = (double *) malloc(block[i]*sizeof(double));
      }
    }
  }
  
  
  double ***L;
  L = (double ***)malloc(nblock*sizeof(double **));
  for(int i=0; i<nblock; i++){
    L[i] = (double **) malloc(n*sizeof(double *));
    for(int j=0; j<n; j++){
      L[i][j] = (double *) malloc(block[i]*sizeof(double));
    }
  }
  
  double** Pest_old=make2Darray(n,type);
  double** Pest_init=make2Darray(n,type);
  
  double ****Best;
  Best = (double ****)malloc(nblock*sizeof(double ***));
  for(int i=0; i<nblock; i++){
    Best[i] = (double ***) malloc(type*sizeof(double **));
    for(int j=0; j<type; j++){
      Best[i][j] = (double **) malloc(nbasis[i]*sizeof(double *));
      for(int l=0; l<nbasis[i]; l++){
        Best[i][j][l] = (double *) malloc(Q*sizeof(double));
      }
    }
  }
  
  
  double value;//Esigma
  int row;
  int col;
  
  
  double*** y=C3darray_diy3(nblock,n,block);
  double** residual=make2Darray(n,G);
  double*** A=C3darray_diy3(nblock,n,nbasis);
  double** UHat=make2Darray(Q,n);
  double** UHatT=make2Darray(Q*type,n);
  double** U=make2Darray(n,Q*type);
  double*** BT=C3darray_diy3(nblock,Q,nbasis);
  double*** BTT=C3darray_diy3(nblock,Q*type,nbasis);
  double thre=0.000001;
  
  
  double ****XCB;
  XCB = (double ****)malloc(nblock*sizeof(double ***));
  for(int i=0; i<nblock; i++){
    XCB[i] = (double ***) malloc(n*sizeof(double **));
    for(int j=0; j<n; j++){
      XCB[i][j] = (double **) malloc(type*sizeof(double *));
      for(int l=0; l<type; l++){
        XCB[i][j][l] = (double *) malloc(block[i]*sizeof(double));
      }
    }
  }
  
  double*** CB=C3darray_diy2(nblock,block,Q);
  double *sigma_vec=(double *) malloc(n*sizeof(double));
  double*** Dmat=make3Darray(n,type,type);
  double** d=make2Darray(n,type);
  
  double dldv=10;
  double dldw=10;
  double *dldv_vec=(double *) malloc(nblock*sizeof(double));
  double *dldw_vec=(double *) malloc(nblock*sizeof(double));
  double *dldv1_vec=(double *) malloc(nblock*sizeof(double));
  double *dldw1_vec=(double *) malloc(nblock*sizeof(double));
  double *dldv2_vec=(double *) malloc(nblock*sizeof(double));
  double *dldw2_vec=(double *) malloc(nblock*sizeof(double));
  double *dldv3_vec=(double *) malloc(nblock*sizeof(double));
  double *dldw3_vec=(double *) malloc(nblock*sizeof(double));
  
  double Mv_old=0;
  double Mw_old=0;
  double sigma_old=1000;
  
  
  double dldv_old=10000000;
  double dldw_old=100000;
  
  double stepsize_w=10;
  double stepsize_v=10;
  
  double *dld2=(double *) malloc(nblock*sizeof(double));
  
  double diff1=1;
  double diff2=1;
  int maxiter=4;
  int iter=1;
  double thr=0.00001;
  double step=0.001;
  
  double *vv=(double *) malloc(n*sizeof(double));
  
  hat1(X,n,Q,UHat);
  
  
  for(int i=0;i<nblock;i++)
  {
    for(int j=0;j<n;j++)
    {
      for(int l=0;l<block[i];l++)
      {
        y[i][j][l]=O[i][j][l];
      }
    }
  }
  
  
  
  for (int i=0;i<n;i++){
    for (int ty=0;ty<type;ty++){
      for (int j=0;j<Q;j++){
        U[i][ty*Q+j]=Pest[i][ty]*X[i][j];
      }
    }
  }
  
  
  hat1(U,n,type*Q,UHatT);
  
  
  for (int b=0;b<nblock;b++){
    
#pragma omp parallel for
    for (int i=0;i<n;i++){
      MatrixVector(ThetaHat[b],O[b][i],A[b][i],nbasis[b],block[b]);
    }
    
    
    MultiplyMatrix(UHatT,A[b],BTT[b],Q*type,nbasis[b],n);
    
    
    for (int ty=0;ty<type;ty++){
      for (int i=0;i<nbasis[b];i++){
        for (int j=0;j<Q;j++){
          
          Best[b][ty][i][j]=BTT[b][Q*ty+j][i];
        }
      }
    }
    
    
  } 
  
  for (int j=0;j<n;j++){
    for (int k=0;k<type;k++){
      Pest_init[j][k]=Pest[j][k];
    }
  }
  
  
  double ***Emup;
  Emup = (double ***)malloc(nblock*sizeof(double **));
  for(int i=0; i<nblock; i++){
    Emup[i] = (double **) malloc(n*sizeof(double *));
    for(int j=0; j<n; j++){
      Emup[i][j] = (double *) malloc(type*block[i]*sizeof(double));
    }
  }
  
  
  double ***BLOCK_mat;
  BLOCK_mat = (double ***)malloc(nblock*sizeof(double **));
  for(int i=0; i<nblock; i++){
    BLOCK_mat[i] = (double **) malloc(block[i]*sizeof(double *));
    for(int j=0; j<block[i]; j++){
      BLOCK_mat[i][j] = (double *) malloc(block[i]*sizeof(double));
    }
  }
  
  
  double ***XWX;
  XWX = (double ***)malloc(nblock*sizeof(double **));
  for(int i=0; i<nblock; i++){
    XWX[i] = (double **) malloc(Q*type*nbasis[i]*sizeof(double *));
    for(int j=0; j<(Q*type*nbasis[i]); j++){
      XWX[i][j] = (double *) malloc(Q*type*nbasis[i]*sizeof(double));
    }
  }
  
  
  double ****beta_sd;
  beta_sd = (double ****)malloc(nblock*sizeof(double ***));
  for(int i=0; i<nblock; i++){
    beta_sd[i] = (double ***) malloc(type*sizeof(double **));
    for(int j=0; j<type; j++){
      beta_sd[i][j] = (double **) malloc(block[i]*sizeof(double *));
      for(int l=0; l<block[i]; l++){
        beta_sd[i][j][l] = (double *) malloc(Q*sizeof(double));
      }
    }
  }
  
  
  
  double*** coef_mat=make3Darray(n,Q*type,Q*type);
  
  
  while(1){
    
    it=it+1;
    if (it>maxit){break;}
    //if (Msigma>0.04){break;}
    
    
    Ll_pre=Ll;
    
    w_old=Mw;
    v_old=Mv;
    
    Msigma_old=Msigma;
    
    for (int j=0;j<nblock;j++){
      
      for (int ty=0;ty<type;ty++){
        MultiplyMatrix(B[j],Best[j][ty],CB[j],block[j],Q,nbasis[j]);
        for (int i=0;i<n;i++){
          MatrixVector(CB[j],X[i],XCB[j][i][ty],block[j],Q);
        }
      }
      
    }
    
    
    
    
    
    for (int j=0;j<n;j++){
      for (int k=0;k<type;k++){
        Pest_old[j][k]=Pest[j][k];
      }
    }
    
    
    
    ///M step
    
    for (int i=0;i<n;i++){
      
      for (int m1=0;m1<type;m1++){
        for (int m2=0;m2<(m1+1);m2++){
          Dmat[i][m1][m2]=0;
          
        }}
      
      
      for (int m1=0;m1<type;m1++){
        d[i][m1]=0;
      }
      
    } 
    
    
    //////M step :p
    
    
    for (int b=0;b<nblock;b++){
      
      gsl_matrix *BLOCK = gsl_matrix_alloc(block[b],block[b]);
      gsl_matrix_set_zero(BLOCK);
      
      for (int ii=0; ii<block[b]; ii++){
        for (int j=0; j<(ii+1); j++){
          value=0;
          if (ii==j){
            value=Mv+0.0000001;
          }
          else{value=Mv*exp((-0.5)*Mw*pow((t[b][ii]-t[b][j]),2));}
          gsl_matrix_set(BLOCK, j, ii, value);
          gsl_matrix_set(BLOCK, ii, j, value);
        }
      }
      gsl_matrix_inv(BLOCK);
      
      
      
      for (int i=0;i<n;i++){
        
        sigma_vec[i]=0; 
        
      }
      
#pragma omp parallel for
      for (int i=0;i<n;i++){
        
        
        gsl_matrix *Sigma_tilde = gsl_matrix_alloc(block[b],block[b]);
        gsl_matrix *BS = gsl_matrix_alloc(block[b],block[b]);
        gsl_vector *mu_k1 = gsl_vector_alloc(block[b]);
        gsl_vector *mu_k2 = gsl_vector_alloc(block[b]);
        gsl_vector *Bmu_k2 = gsl_vector_alloc(block[b]);
        gsl_matrix *Pi = gsl_matrix_alloc(block[b],type*block[b]);
        gsl_matrix *Pi_T = gsl_matrix_alloc(type*block[b],block[b]);
        gsl_matrix *Pi_TB = gsl_matrix_alloc(type*block[b],block[b]);
        gsl_matrix *PBP = gsl_matrix_alloc(type*block[b],type*block[b]);
        gsl_vector *O_sub = gsl_vector_alloc(block[b]);
        gsl_vector *Emu_sub = gsl_vector_alloc(type*block[b]);
        gsl_vector *Emu_sub_nor = gsl_vector_alloc(type*block[b]);
        gsl_vector *Emu_sub1 = gsl_vector_alloc(type*block[b]);
        
        gsl_matrix_set_zero(Pi_T);
        gsl_matrix_set_zero(Pi);
        gsl_matrix_set_zero(Sigma_tilde);
        gsl_matrix_set_zero(BS);
        gsl_matrix_set_zero(Pi_TB);
        gsl_matrix_set_zero(PBP);
        gsl_vector_set_zero(mu_k1);
        gsl_vector_set_zero(mu_k2);
        gsl_vector_set_zero(Bmu_k2);
        gsl_vector_set_zero(O_sub);
        gsl_vector_set_zero(Emu_sub);
        gsl_vector_set_zero(Emu_sub_nor);
        gsl_vector_set_zero(Emu_sub1);
        
        
        
        
        for(int k=0;k<type;k++){
          for(int j=0;j<block[b];j++){
            
            gsl_matrix_set(Pi, j, j*type+k, Pest[i][k]);
            gsl_matrix_set(Pi_T,  j*type+k,j, Pest[i][k]);
            
          }
        }
        
        
        
        gsl_matrix_mul(Pi_T, BLOCK, Pi_TB);
        
        
        
        
        gsl_matrix_mul(Pi_TB, Pi, PBP);
        
        
        
        for(int j=0;j<(type*block[b]);j++){
          
          vv[i]=gsl_matrix_get(PBP, j, j);
          gsl_matrix_set(PBP, j, j, vv[i]+1/Msigma);
          
        }
        
        gsl_matrix_inv(PBP);
        
        
        
        
        
        sigma_vec[i]=sigma_vec[i]+gsl_matrix_trace(PBP);
        
        for(int j=0;j<block[b];j++){
          
          gsl_vector_set(O_sub, j, O[b][i][j]);
          
        }
        
        
        //gsl_matrixvector_mul(Pi_TB,O_sub,Emu_sub,type*block[b],block[b]);
        
        for (int j1=0;j1<(type*block[b]);j1++){
          
          vv[i]=0;
          
          for (int j2=0;j2<(block[b]);j2++){
            vv[i]=vv[i]+gsl_matrix_get(Pi_TB, j1, j2)*gsl_vector_get(O_sub, j2);
          }
          
          gsl_vector_set(Emu_sub, j1, vv[i]);
        }
        
        for (int j1=0;j1<(type);j1++){
          for (int j2=0;j2<block[b];j2++){
            gsl_vector_set(Emu_sub1,type*j2+j1,gsl_vector_get(Emu_sub,type*j2+j1)+XCB[b][i][j1][j2]/Msigma);
          }}
        
        
        gsl_matrixvector_mul(PBP,Emu_sub1,Emu_sub_nor,type*block[b],type*block[b]);
        
        
        
        
        for (int j1=0;j1<(type);j1++){
          for (int j2=0;j2<block[b];j2++){
            Emu[b][i][type*j2+j1]=gsl_vector_get(Emu_sub_nor,type*j2+j1);
            
          }}
        
        
        for (int j1=0;j1<(type);j1++){
          for (int j2=0;j2<block[b];j2++){
            Emup[b][i][type*j2+j1]=gsl_vector_get(Emu_sub,type*j2+j1);
          }}
        
        for (int m1=0;m1<type;m1++){
          for (int m2=0;m2<(m1+1);m2++){
            
            
            
            
            for (int i1=0; i1<block[b]; i1++){
              
              for (int i2=0; i2<block[b]; i2++){
                gsl_matrix_set(Sigma_tilde, i1, i2, gsl_matrix_get(PBP,i1*type+m1,i2*type+m2));
                
              }
            }
            
            
            gsl_matrix_mul(BLOCK, Sigma_tilde, BS);
            Dmat[i][m1][m2]= Dmat[i][m1][m2]+gsl_matrix_trace(BS);
            for (int i1=0; i1<block[b]; i1++){
              
              gsl_vector_set(mu_k1, i1, Emu[b][i][i1*type+m1]);
              
            }
            
            
            for (int i1=0; i1<block[b]; i1++){
              
              gsl_vector_set(mu_k2, i1, Emu[b][i][i1*type+m2]);
              
            }
            
            
            for (int j1=0;j1<(block[b]);j1++){
              
              vv[i]=0;
              
              for (int j2=0;j2<(block[b]);j2++){
                vv[i]=vv[i]+gsl_matrix_get(BLOCK, j1, j2)*gsl_vector_get(mu_k2, j2);
              }
              
              gsl_vector_set(Bmu_k2, j1, vv[i]);
            }
            
            
            //gsl_matrixvector_mul(BLOCK,mu_k2,Bmu_k2,block[b],block[b]);
            
            for (int i1=0; i1<block[b]; i1++){
              
              Dmat[i][m1][m2]=Dmat[i][m1][m2]+gsl_vector_get(mu_k1, i1)*gsl_vector_get(Bmu_k2, i1);
              
            }
            Dmat[i][m2][m1]= Dmat[i][m1][m2];
            
            
            
            
            
          }
          
          
        }
        
        
        for (int m=0;m<type;m++){
          
          
          for (int i1=0; i1<block[b]; i1++){
            
            gsl_vector_set(mu_k1, i1, Emu[b][i][i1*type+m]);
            
          }
          
          for (int j1=0;j1<(block[b]);j1++){
            
            vv[i]=0;
            
            for (int j2=0;j2<(block[b]);j2++){
              vv[i]=vv[i]+gsl_matrix_get(BLOCK, j1, j2)*gsl_vector_get(mu_k1, j2);
            }
            
            gsl_vector_set(Bmu_k2, j1, vv[i]);
          }
          
          
          for (int i1=0; i1<block[b]; i1++){
            
            d[i][m]=d[i][m]+gsl_vector_get(Bmu_k2, i1)*O[b][i][i1];
            
          }
          
          
        }
        
        
        gsl_matrix_free(BS);
        gsl_matrix_free(Sigma_tilde);
        gsl_vector_free(mu_k1);
        gsl_vector_free(mu_k2);
        gsl_vector_free(Bmu_k2);
        gsl_matrix_free(Pi);
        gsl_matrix_free(Pi_T);
        gsl_matrix_free(Pi_TB);
        gsl_matrix_free(PBP);
        gsl_vector_free(O_sub);
        gsl_vector_free(Emu_sub);
        gsl_vector_free(Emu_sub1);
        gsl_vector_free(Emu_sub_nor);
      }
      
      gsl_matrix_free(BLOCK);
      
      
    }
    
    
    for (int i=0;i<n;i++){
      quadprog(Dmat[i],d[i],Pest[i],type);
    }
    
   
    
    
    
    Ll=0;
    for (int b=0;b<nblock;b++){
      gsl_matrix *BLOCK = gsl_matrix_alloc(block[b],block[b]);
      gsl_vector *x_sub = gsl_vector_alloc(block[b]);
      
      gsl_matrix_set_zero(BLOCK);
      gsl_vector_set_zero(x_sub);
      for (int i=0;i<n;i++){
        
        for (int ii=0; ii<block[b]; ii++){
          
          for (int j=0; j<(ii+1); j++){
            
            if (ii==j){
              value=Mv;
              for (int ty=0;ty<type;ty++){
                value=value+Pest[i][ty]*Pest[i][ty]*Msigma;
              }
              
            }
            else{value=Mv*exp((-0.5)*Mw*pow((t[b][ii]-t[b][j]),2));}
            gsl_matrix_set(BLOCK, j, ii, value);
            gsl_matrix_set(BLOCK, ii, j, value);
          }
        }
        BLOCK_del=get_det(BLOCK);
        gsl_matrix_inv(BLOCK);
        Ll=Ll+0.5*BLOCK_del;
        
        for (int j=0;j<block[b];j++){
          L[b][i][j]=O[b][i][j];
          for (int ty=0;ty<type;ty++){
            L[b][i][j]=L[b][i][j]-Pest[i][ty]*XCB[b][i][ty][j];
          }
          gsl_vector_set(x_sub,  j, L[b][i][j]);
        }
        Ll=Ll-0.5*gsl_VMV(x_sub,BLOCK,x_sub,block[b],block[b]);
      }
      gsl_matrix_free(BLOCK);
      gsl_vector_free(x_sub);
    }
    
    
    
    
    Ll_old=0;
    for (int b=0;b<nblock;b++){
      gsl_matrix *BLOCK = gsl_matrix_alloc(block[b],block[b]);
      gsl_vector *x_sub = gsl_vector_alloc(block[b]);
      gsl_matrix_set_zero(BLOCK);
      gsl_vector_set_zero(x_sub);
      
      for (int i=0;i<n;i++){
        
        for (int ii=0; ii<block[b]; ii++){
          
          for (int j=0; j<(ii+1); j++){
            
            if (ii==j){
              value=Mv;
              for (int ty=0;ty<type;ty++){
                value=value+Pest_old[i][ty]*Pest_old[i][ty]*Msigma;
              }
              
            }
            else{value=Mv*exp((-0.5)*Mw*pow((t[b][ii]-t[b][j]),2));}
            gsl_matrix_set(BLOCK, j, ii, value);
            gsl_matrix_set(BLOCK, ii, j, value);
          }
        }
        BLOCK_del=get_det(BLOCK);
        gsl_matrix_inv(BLOCK);
        Ll_old=Ll_old+0.5*BLOCK_del;
        
        for (int j=0;j<block[b];j++){
          L[b][i][j]=O[b][i][j];
          for (int ty=0;ty<type;ty++){
            L[b][i][j]=L[b][i][j]-Pest_old[i][ty]*XCB[b][i][ty][j];
          }
          gsl_vector_set(x_sub,  j, L[b][i][j]);
        }
        Ll_old=Ll_old-0.5*gsl_VMV(x_sub,BLOCK,x_sub,block[b],block[b]);
      }
      gsl_matrix_free(BLOCK);
      gsl_vector_free(x_sub);
    }
    
    
    
    //Ck
    
    for (int k=0;k<type;k++){
      
#pragma omp parallel for
      for (int j=0;j<nblock;j++){
        for (int i=0;i<n;i++){
          
          for (int l=0;l<block[j];l++){
            
            y[j][i][l]=Emu[j][i][type*l+k];
            
          }
          
          MatrixVector(ThetaHat[j],y[j][i],A[j][i],nbasis[j],block[j]);
        }
        
        MultiplyMatrix(UHat,A[j],BT[j],Q,nbasis[j],n);
        
        
        for (int i=0;i<nbasis[j];i++){
          for (int l=0;l<Q;l++){
            Best[j][k][i][l]=BT[j][l][i];
          }
        }
        
        
        
      }
      
      
    }
    
    for (int j=0;j<nblock;j++){
      for (int i=0;i<type;i++){
        MultiplyMatrix(B[j],Best[j][i],beta_est[j][i],block[j],Q,nbasis[j]);
      }}
    
    
    //epsilon
    
    Msigma=0;
    for (int i=0;i<n;i++){
      
      Msigma=Msigma+sigma_vec[i];
    }
    
    
#pragma omp parallel for
    for (int j=0;j<nblock;j++){
      
      for (int ty=0;ty<type;ty++){
        MultiplyMatrix(B[j],Best[j][ty],CB[j],block[j],Q,nbasis[j]);
        for (int i=0;i<n;i++){
          MatrixVector(CB[j],X[i],XCB[j][i][ty],block[j],Q);
        }
      }
      
    }
    
    
    
    
    
    
    
    
    
    
    
    
    for (int l=0;l<nblock;l++){
      
      for (int j=0;j<block[l];j++){
        
        for (int i=0;i<n;i++){
          
          for (int k=0;k<type;k++){
            Msigma=Msigma+XCB[l][i][k][j]*XCB[l][i][k][j]-2*Emu[l][i][j*type+k]*XCB[l][i][k][j]+Emu[l][i][j*type+k]*Emu[l][i][j*type+k];
            
          }
        }}}
    
    Msigma=(Msigma+2*prior_beta_sigma)/(n*G*type+2*(prior_alpha_sigma+1));    
    
    
    //v,w
    
    
    value=0;
    diff1=1;
    diff2=1;
    iter=0;
    
    
    batch=it%10;
    
    //#pragma omp parallel for
    for (int b=0;b<nblock;b++){    
      
      gsl_matrix *BLOCK = gsl_matrix_alloc(block[b],block[b]);
      gsl_matrix_set_zero(BLOCK);
      
      for (int ii=0; ii<block[b]; ii++){
        for (int j=0; j<(ii+1); j++){
          
          if (ii==j){
            value=Mv+0.00000001;
          }
          else{value=Mv*exp((-0.5)*Mw*pow((t[b][ii]-t[b][j]),2));}
          gsl_matrix_set(BLOCK, j, ii, value);
          gsl_matrix_set(BLOCK, ii, j, value);
        }
      }
      gsl_matrix_inv(BLOCK);
      
    
      
      for (int i=(batch*n_batch+1);i<(batch*n_batch+n_batch);i++){
        
        gsl_matrix *Pi_old = gsl_matrix_alloc(block[b],type*block[b]);
        gsl_matrix *Pi_T_old = gsl_matrix_alloc(type*block[b],block[b]);
        gsl_matrix *Pi_TB_old = gsl_matrix_alloc(type*block[b],block[b]);
        gsl_matrix *Pi = gsl_matrix_alloc(block[b],type*block[b]);
        gsl_matrix *Pi_T = gsl_matrix_alloc(type*block[b],block[b]);
        
        gsl_matrix *PBP = gsl_matrix_alloc(type*block[b],type*block[b]);
        
        gsl_matrix *PSP = gsl_matrix_alloc(block[b],block[b]);
        gsl_matrix *PS_T = gsl_matrix_alloc(block[b],type*block[b]);
        
        gsl_matrix_set_zero(Pi_old);
        gsl_matrix_set_zero(Pi_T_old);
        gsl_matrix_set_zero(Pi_TB_old);
        gsl_matrix_set_zero(Pi);
        gsl_matrix_set_zero(Pi_T);
        gsl_matrix_set_zero(PBP);
        gsl_matrix_set_zero(PSP);
        gsl_matrix_set_zero(PS_T);
        
        
        for(int k=0;k<type;k++){
          for(int j=0;j<block[b];j++){
            
            gsl_matrix_set(Pi_old, j, j*type+k, Pest_old[i][k]);
            gsl_matrix_set(Pi_T_old, j*type+k,j , Pest_old[i][k]);
          }
        }
        
        //gsl_matrix_trans(Pi_old,Pi_T_old);
        gsl_matrix_mul(Pi_T_old, BLOCK, Pi_TB_old);
        gsl_matrix_mul(Pi_TB_old, Pi_old, PBP);
        
        
        for(int j=0;j<(type*block[b]);j++){
          
          gsl_matrix_set(PBP, j, j, gsl_matrix_get(PBP, j, j)+1/Msigma_old);
          
        }
        
        gsl_matrix_inv(PBP);
        
        
        for(int k=0;k<type;k++){
          for(int j=0;j<block[b];j++){
            
            gsl_matrix_set(Pi, j, j*type+k, Pest[i][k]);
            
          }
        }
        
        
        gsl_matrix_mul_pa(Pi, PBP, PS_T);
        gsl_matrix_mul_pa(PS_T, Pi_T, PSP);
        
        
        for (int j1=0;j1<block[b];j1++){
          for (int j2=0;j2<block[b];j2++){
            PSP_mat[b][i][j1][j2]=gsl_matrix_get(PSP,j1,j2);
          }}
        
        
        gsl_matrix_free(Pi_old);
        gsl_matrix_free(Pi_T_old);
        gsl_matrix_free(Pi_TB_old);
        gsl_matrix_free(Pi);
        gsl_matrix_free(Pi_T);
        gsl_matrix_free(PBP);
        
        gsl_matrix_free(PSP);
        gsl_matrix_free(PS_T);
        
        
      }
      
      gsl_matrix_free(BLOCK);
      
    }
    
    
    while (1){
      iter=iter+1;
      
      
#pragma omp parallel for
      for (int i=0;i<nblock;i++){
        
        
        gsl_matrix *Phitov = gsl_matrix_alloc(block[i],block[i]);
        gsl_matrix *Phitow = gsl_matrix_alloc(block[i],block[i]);
        gsl_matrix *Phi = gsl_matrix_alloc(block[i],block[i]); 
        gsl_matrix *IP = gsl_matrix_alloc(block[i],block[i]);//v,w
        gsl_matrix *IPI = gsl_matrix_alloc(block[i],block[i]); 
        gsl_matrix *IPIO = gsl_matrix_alloc(block[i],block[i]); 
        gsl_vector *u = gsl_vector_alloc(block[i]);
        gsl_vector *Su = gsl_vector_alloc(block[i]);
        gsl_matrix *PSP = gsl_matrix_alloc(block[i],block[i]);
        //double *uSu = &dld2;
        
        
        gsl_matrix_set_zero(Phitov);
        gsl_matrix_set_zero(Phitow);
        gsl_matrix_set_zero(Phi);
        gsl_matrix_set_zero(IP);
        gsl_matrix_set_zero(IPI);
        gsl_matrix_set_zero(IPIO);
        gsl_matrix_set_zero(PSP);
        gsl_vector_set_zero(u);
        gsl_vector_set_zero(Su);
        
        PhitoV(Mw, t[i],block[i],Phitov);
        PhitoW(Mw, Mv, t[i],block[i],Phitow);
        PhiMatrix(Mw, Mv, t[i],block[i],Phi);
        
        gsl_matrix_inv(Phi);
        gsl_matrix_mul(Phi,Phitov,IP);
        gsl_matrix_mul(IP,Phi,IPI);
        
        dldv2_vec[i]=0;
        
        for (int r=(batch*n_batch+1);r<(batch*n_batch+n_batch);r++){
          
          for (int m=0;m<block[i];m++){
            value=0;
            for (int ty=0;ty<type;ty++){
              value=value+Pest[r][ty]*Emu[i][r][m*type+ty];
            }
            gsl_vector_set(u, m, O[i][r][m]-value);
          }
          
          gsl_blas_dgemv(CblasNoTrans,1,IPI,u,0,Su);
          
          dld2[i]=0;
          for (int l=0;l<block[i];l++){
            dld2[i]=dld2[i]+gsl_vector_get(u,l)*gsl_vector_get(Su,l);
          }
          
          //gsl_blas_ddot(u, Su, uSu);
          
          dldv2_vec[i]=dldv2_vec[i]+dld2[i];
        }
        
        dldv3_vec[i]=0;
        for (int r=(batch*n_batch+1);r<(batch*n_batch+n_batch);r++){
          
          
          
          for(int j1=0;j1<(block[i]);j1++){
            for(int j2=0;j2<(block[i]);j2++){
              
              gsl_matrix_set(PSP, j1, j2, PSP_mat[i][r][j1][j2]);
              
            }
          }
          
          
          gsl_matrix_mul(IPI,PSP,IPIO);
          
          dldv3_vec[i]=dldv3_vec[i]+0.5*gsl_matrix_trace(IPIO);
        }
        
        dldv1_vec[i]=-0.5*gsl_matrix_trace(IP);
        dldv2_vec[i]=0.5*dldv2_vec[i]/n_batch;
        dldv3_vec[i]=dldv3_vec[i]/n_batch;
        
        gsl_matrix_mul(Phi,Phitow,IP);
        gsl_matrix_mul(IP,Phi,IPI);
        
        gsl_vector_set_zero(Su);
        
        dldw2_vec[i]=0;
        
        for (int r=(batch*n_batch+1);r<(batch*n_batch+n_batch);r++){
          
          for (int m=0;m<block[i];m++){
            value=0;
            for (int ty=0;ty<type;ty++){
              value=value+Pest[r][ty]*Emu[i][r][m*type+ty];
            }
            gsl_vector_set(u, m, O[i][r][m]-value);
          }
          
          gsl_blas_dgemv(CblasNoTrans,1,IPI,u,0,Su);
          
          dld2[i]=0;
          for (int l=0;l<block[i];l++){
            dld2[i]=dld2[i]+gsl_vector_get(u,l)*gsl_vector_get(Su,l);
          }
          
          //gsl_blas_ddot(u, Su, uSu);
          
          dldw2_vec[i]=dldw2_vec[i]+dld2[i];
        }
        
        dldw3_vec[i]=0;
        
        for (int r=(batch*n_batch+1);r<(batch*n_batch+n_batch);r++){
          
          
          
          for(int j1=0;j1<(block[i]);j1++){
            for(int j2=0;j2<(block[i]);j2++){
              
              gsl_matrix_set(PSP, j1, j2, PSP_mat[i][r][j1][j2]);
              
            }
          }
          
          gsl_matrix_mul(IPI,PSP,IPIO);
          
          dldw3_vec[i]=dldw3_vec[i]+0.5*gsl_matrix_trace(IPIO);
        }
        
        dldw1_vec[i]=-0.5*gsl_matrix_trace(IP);
        dldw2_vec[i]=0.5*dldw2_vec[i]/n_batch;
        dldw3_vec[i]=dldw3_vec[i]/n_batch;
        
        
        
        
        dldw_vec[i]=dldw1_vec[i]+dldw2_vec[i]+dldw3_vec[i];
        dldv_vec[i]=dldv1_vec[i]+dldv2_vec[i]+dldv3_vec[i];
        
        
        
        gsl_matrix_free(Phitov);
        gsl_matrix_free(Phitow);
        gsl_matrix_free(Phi);
        gsl_matrix_free(IP);
        gsl_matrix_free(IPI);
        gsl_matrix_free(IPIO);
        gsl_vector_free(u);
        gsl_vector_free(Su);
        gsl_matrix_free(PSP);
        
      }
      
      
      dldv=0;
      dldw=0;
      for (int i=0;i<nblock;i++){
        dldw=dldw+dldw_vec[i];
        dldv=dldv+dldv_vec[i];
      }
      dldw=dldw-1/n_batch*(alpha_w+1)/Mw+n_batch*alpha_w/mu_w/Mw/Mw;
      dldv=dldv-1/n_batch*1/Mv*(1+(log(Mv)-mu_v)/sigma_v);
      
      stepsize_w=(Mw-Mw_old)/(dldw-dldw_old);
      stepsize_v=(Mv-Mv_old)/(dldv-dldv_old);
      
      Mv_old=Mv;
      Mw_old=Mw;
      Mw=Mw-stepsize_w*dldw;
      Mv=Mv-stepsize_v*dldv;
      
      dldw_old=dldw;
      dldv_old=dldv;
      if (Mw>100){Mw=Mw_old;}
      if (Mw<0){Mw=Mw_old;}
      
      if (Mv<0){Mv=Mv_old;}
      if (Mv>1){Mv=Mv_old;break;}
      
      if(fabs(dldw)<thr && fabs(dldv)<thr){
        break;
      }
      if(fabs(Mw-Mw_old)<0.0001 && fabs(Mv-Mv_old)<0.0001){
        break;
      }
      if(fabs(dldw_old) <fabs(dldw)&& fabs(dldv_old) <fabs(dldv)){
        Mv=Mv_old;
        Mw=Mw_old;
        break;
      }
      //if(it>=3&&iter>maxitertry){break;}
      if(iter>maxiter){break;}
      
      
      
    }
    
    
  }
  
  
  for (int j=0;j<nblock;j++){
    for (int i=0;i<type;i++){
      MultiplyMatrix(B[j],Best[j][i],beta_est[j][i],block[j],Q,nbasis[j]);
    }}
  
  
  for (int i1=0;i1<type;i1++){
    for (int i3=0;i3<type;i3++){
      for (int i2=0;i2<Q;i2++){
        for (int i4=0;i4<Q;i4++){
          
          for (int j=0;j<n;j++){
            row=(i1*Q+i2);
            col=(i3*Q+i4);
            coef_mat[j][row][col]=X[j][i2]*Pest[j][i1]*X[j][i4]*Pest[j][i3];
          }
          
        }
      }
    }
  }
  
  for (int b=0;b<nblock;b++){
    gsl_matrix *BLOCK = gsl_matrix_alloc(block[b],block[b]);
    gsl_matrix *B_block = gsl_matrix_alloc(block[b],nbasis[b]); 
    gsl_matrix *B_block_T = gsl_matrix_alloc(nbasis[b],block[b]); 
    gsl_matrix *BG = gsl_matrix_alloc(nbasis[b],block[b]); 
    gsl_matrix *BGB = gsl_matrix_alloc(nbasis[b],nbasis[b]); 
    double** XWX_sub=make2Darray(Q*type*nbasis[b],Q*type*nbasis[b]);//C
    gsl_matrix *XWX_sub_gsl = gsl_matrix_alloc(Q*type*nbasis[b],Q*type*nbasis[b]); 
    
    
    
    gsl_matrix_set_zero(BLOCK);
    gsl_matrix_set_zero(B_block);
    gsl_matrix_set_zero(B_block_T);
    gsl_matrix_set_zero(BG);
    gsl_matrix_set_zero(BGB);
    gsl_matrix_set_zero(XWX_sub_gsl);
    
    for (int i1=0;i1<block[b];i1++){
      
      for (int i2=0;i2<nbasis[b];i2++){
        
        gsl_matrix_set(B_block, i1, i2, B[b][i1][i2]);
        gsl_matrix_set(B_block_T, i2, i1, B[b][i1][i2]);
      }
      
    }
    
    for (int i1=0;i1<(Q*type*nbasis[b]);i1++){
      
      for (int i2=0;i2<(Q*type*nbasis[b]);i2++){
        
        XWX_sub[i1][i2]=0;
      }
      
    }
    
    
    //gsl_matrix_trans(B_block,B_block_T);
    
    
    for (int i=0;i<n;i++){
      
      for (int ii=0; ii<block[b]; ii++){
        
        for (int j=0; j<(ii+1); j++){
          
          if (ii==j){
            value=0;
            for (int ty=0;ty<type;ty++){
              value=Mv+Pest[i][ty]*Pest[i][ty]*Msigma;
            }
            
          }
          else{value=Mv*exp((-0.5)*Mw*pow((t[b][ii]-t[b][j]),2));}
          gsl_matrix_set(BLOCK, j, ii, value);
          gsl_matrix_set(BLOCK, ii, j, value);
        }
      }
      
      gsl_matrix_inv(BLOCK);
      gsl_matrix_mul(B_block_T, BLOCK, BG);
      gsl_matrix_mul(BG, B_block, BGB);
      
      
      for (int i1=0;i1<(type*Q);i1++){
        for (int i2=0;i2<(type*Q);i2++){
          
          for (int j1=0;j1<nbasis[b];j1++){
            for (int j2=0;j2<nbasis[b];j2++){
              XWX_sub[i1*nbasis[b]+j1][i2*nbasis[b]+j2]=XWX_sub[i1*nbasis[b]+j1][i2*nbasis[b]+j2]+coef_mat[i][i1][i2]*gsl_matrix_get(BGB, j1, j2);
            }
          }
          
        }
      }
      
      
      
    }
    
    
    for (int i1=0;i1<(type*Q*nbasis[b]);i1++){
      for (int i2=0;i2<(type*Q*nbasis[b]);i2++){
        gsl_matrix_set(XWX_sub_gsl, i1, i2, XWX_sub[i1][i2]);
      }}
    
    
    gsl_matrix_inv(XWX_sub_gsl);
    
    for (int i1=0;i1<(type*Q*nbasis[b]);i1++){
      for (int i2=0;i2<(type*Q*nbasis[b]);i2++){
        XWX[b][i1][i2]=gsl_matrix_get(XWX_sub_gsl, i1, i2);
      }}
    
    gsl_matrix_free(BLOCK);
    gsl_matrix_free(B_block);
    gsl_matrix_free(B_block_T);
    gsl_matrix_free(BG);
    gsl_matrix_free(BGB);
    gsl_matrix_free(XWX_sub_gsl);
    delet2Darray(XWX_sub,Q*type*nbasis[b],Q*type*nbasis[b]);
    
  }
  
  
  for (int ty=0;ty<type;ty++){
    
    for (int q=1;q<Q;q++){
      
      for (int i=0;i<nblock;i++){
        
        for (int j=0;j<block[i];j++){
          
          
          var_trans=0;
          for (int j1=0;j1<(nbasis[i]-1);j1++){
            row=j1+nbasis[i]*q+ty*nbasis[i]*Q;
            var_trans=var_trans+pow(B[i][j][j1],2)*XWX[i][row][row];
            for (int j2=(j1+1);j2<nbasis[i];j2++){
              col=j2+nbasis[i]*q+ty*nbasis[i]*Q;
              var_trans=var_trans+2*B[i][j][j1]*B[i][j][j2]*XWX[i][row][col];
            }
          }
          row=nbasis[i]*(q+1)+ty*nbasis[i]*Q-1;
          var_trans=var_trans+pow(B[i][j][nbasis[i]-1],2)*XWX[i][row][row];
          beta_sd[i][ty][j][q]=sqrt(var_trans);
          beta_z[i][ty][j][q]=beta_est[i][ty][j][q]/beta_sd[i][ty][j][q];
        }
      }
    }
    
  }
  
  
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
  for (int j = 0; j < nblock; j++){
    for (int l = 0; l < (block[j]);l++){
      
      for (int i = 0; i < n; i++){
        for (int k=0;k<type;k++){
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
  
  for (int i=0;i<nblock;i++){
    
    for (int j=0;j<block[i];j++){    
      for (int q=1;q<Q;q++){
        for (int ty=0;ty<type;ty++){
          
          REAL(res_beta_c)[coun] = beta_est[i][ty][j][q];
          coun=coun+1;
          
        }}}}
  
  
  
  PROTECT(res_beta_c_dim = allocVector(INTSXP, 3));
  ++nProtected;	
  INTEGER(res_beta_c_dim)[0] = type;
  INTEGER(res_beta_c_dim)[1] = Q;
  INTEGER(res_beta_c_dim)[2] = G;
  setAttrib(res_beta_c, R_DimSymbol, res_beta_c_dim);
  
  //return z
  SEXP res_z_c, res_z_c_dim;
  PROTECT(res_z_c = allocVector(REALSXP, G*Q*type));
  ++nProtected;
  
  coun=0;
  for (int i=0;i<nblock;i++){
    
    for (int j=0;j<block[i];j++){
      for (int q=0;q<Q;q++){
        for (int ty=0;ty<type;ty++){
          
          REAL(res_z_c)[coun] = beta_z[i][ty][j][q];
          coun=coun+1;
          
        }}}}
  
  
  
  PROTECT(res_z_c_dim = allocVector(INTSXP, 3));
  ++nProtected;	
  INTEGER(res_z_c_dim)[0] = type;
  INTEGER(res_z_c_dim)[1] = Q;
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
  setAttrib(resList, R_NamesSymbol, names);
  
  UNPROTECT(nProtected);
  
  
  return resList;		     				     							
}

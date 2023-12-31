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

double** make2Darray(int a1, int a2){
  double **tmp;
  tmp = (double **)malloc(a1*sizeof(double *));
  for(int i=0; i<a1; i++){
    tmp[i] = (double *) malloc(a2*sizeof(double));
  }
  return tmp;
}

void delet2Darray(double **tmp, int a1, int a2){
  for(int i=0; i<a1; i++){
    free(tmp[i]);
  }
  free(tmp);
}

double*** make3Darray(int a1, int a2, int a3){
  double ***tmp;
  tmp = (double ***)malloc(a1*sizeof(double **));
  for(int i=0; i<a1; i++){
    tmp[i] = (double **) malloc(a2*sizeof(double *));
    for(int j=0; j<a2; j++){
      tmp[i][j] = (double *) malloc(a3*sizeof(double));
    }
  }
  return tmp;
}

void delet3Darray(double ***tmp, int a1, int a2, int a3){
  for(int i=0; i<a1; i++){
    for(int j=0; j<a2; j++){
      free(tmp[i][j]);
    }
    free(tmp[i]);
  }
  free(tmp);
}

void gsl_matrix_inv(gsl_matrix *a)
{
  size_t n=a->size1;
  size_t m=a->size2;
  
  gsl_matrix *temp1=gsl_matrix_calloc(n,n);
  gsl_matrix_memcpy(temp1,a);
  
  gsl_permutation *p=gsl_permutation_calloc(n);
  int sign=0;
  gsl_linalg_LU_decomp(temp1,p,&sign);
  gsl_matrix *inverse=gsl_matrix_calloc(n,n);
  
  gsl_linalg_LU_invert(temp1,p,inverse);
  gsl_matrix_memcpy(a,inverse);
  
  gsl_permutation_free(p);
  gsl_matrix_free(temp1);
  gsl_matrix_free(inverse);
  
}


void MultiplyMatrix(double **a,double **b,double **c,int ROW, int COL,int RC){
  int i,j,k;
  
  for (i=0; i<ROW; i++){
    for(j=0; j<COL; j++){
      c[i][j]= 0 ;
      for(k=0; k<RC; k++){
        c[i][j]= c[i][j] + a[i][k]*b[k][j] ;
        //c[i][j] += a[i][k]*b[k][j] ;
      }
    }
  }
}

double Determinant(double **a,int n){
  int i,j,j1,j2;
  double det = 0;
  double **m=NULL;
  
  if (n < 1) { /* Error */
    
  } else if (n == 1) { /* Shouldn't get used */
    det = a[0][0];
  } else if (n == 2) {
    det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
  } else {
    det = 0;
    for (j1=0;j1<n;j1++) {
      m = (double **) malloc((n-1)*sizeof(double *));
      for (i=0;i<n-1;i++)
        m[i] = (double *) malloc((n-1)*sizeof(double));
      for (i=1;i<n;i++) {
        j2 = 0;
        for (j=0;j<n;j++) {
          if (j == j1)
            continue;
          m[i-1][j2] = a[i][j];
          j2++;
        }
      }
      det += pow(-1.0,j1+2.0) * a[0][j1] * Determinant(m,n-1);
      for (i=0;i<n-1;i++)
        free(m[i]);
      free(m);
    }
  }
  return(det);
}

/*Find the cofactor matrix of a square matrix*/

void CoFactor(double **a, int n, double **b){
  int i,j,ii,jj,i1,j1;
  double det;
  double **c;
  
  c = (double **) malloc((n-1)*sizeof(double *));
  for (i=0;i<n-1;i++)
    c[i] = (double *) malloc((n-1)*sizeof(double));
  
  for (j=0;j<n;j++) {
    for (i=0;i<n;i++) {
      
      /* Form the adjoint a_ij */
      i1 = 0;
      for (ii=0;ii<n;ii++) {
        if (ii == i)
          continue;
        j1 = 0;
        for (jj=0;jj<n;jj++) {
          if (jj == j)
            continue;
          c[i1][j1] = a[ii][jj];
          j1++;
        }
        i1++;
      }
      
      /* Calculate the determinate */
      det = Determinant(c,n-1);
      
      /* Fill in the elements of the cofactor */
      b[i][j] = pow(-1.0,i+j+2.0) * det;
    }
  }
  for (i=0;i<n-1;i++)
    free(c[i]);
  free(c);
}

/*Transpose of a square matrix, do it in place*/
void Transpose(double **a,int n){
  int i,j;
  double tmp;
  
  for (i=1;i<n;i++) {
    for (j=0;j<i;j++) {
      tmp = a[i][j];
      a[i][j] = a[j][i];
      a[j][i] = tmp;
    }
  }
}

void transpose(double **x, double **y,int p,int q)//Transpose
{
  
  for (int i = 0; i < p; i++)
  {
    for (int j = 0; j < q; j++)
    {
      y[j][i] = x[i][j];
    }
  }
}

/*calculate the inverse*/
void inverse(double **a, int n, double **a_inv){
  double det;
  double  **cofac_a;
  cofac_a = (double **)malloc(n*sizeof(double*));
  for(int i =0; i <n;i++){
    cofac_a[i] = (double *)malloc(n*sizeof(double));
  }
  CoFactor(a, n, cofac_a);
  Transpose(cofac_a, n); //turn the cofacotor matrix into the adjoint matrix
  det = Determinant(a, n); 
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      a_inv[i][j] = cofac_a[i][j] / det;
    }
  }
  
  for(int i=0; i<n; i++){
    free(cofac_a[i]);
  }
  free(cofac_a);
}

void hat(double **B, int row, int col, double **Theta){
  double** BT=make2Darray(col,row);
  transpose(B,BT,row,col);
  double** BTB=make2Darray(col,col);
  MultiplyMatrix(BT,B,BTB,col,col,row);
  double** BTBI=make2Darray(col,col);
  inverse(BTB,col,BTBI);
  MultiplyMatrix(BTBI,BT,Theta,col,row,col);
  
  delet2Darray(BT,col,row);
  delet2Darray(BTB,col,col);
  delet2Darray(BTBI,col,col);
}

void hat1(double **B, int row, int col, double **Theta){
  gsl_matrix *BTB = gsl_matrix_alloc(col,col); 
  double** BT=make2Darray(col,row);
  transpose(B,BT,row,col);
  double value=0;
  double** BTBI=make2Darray(col,col);
  for (int i=0;i<col;i++){
    for (int j=0;j<(i+1);j++){
      value=0;
      for (int l=0;l<row;l++){
        value=value+B[l][i]*B[l][j];
      }
      gsl_matrix_set(BTB, i, j, value);
      gsl_matrix_set(BTB, j, i, value);
    }
  }
  gsl_matrix_inv(BTB);
  for (int i=0;i<col;i++){
    for (int j=0;j<col;j++){
      BTBI[i][j]=gsl_matrix_get(BTB,i,j);
    }
  }
  MultiplyMatrix(BTBI,BT,Theta,col,row,col);
  
  delet2Darray(BT,col,row);
  gsl_matrix_free(BTB);
  delet2Darray(BTBI,col,col);
}

void MatrixVector(double **X,double *y,double *b,int row, int col){
  
  for (int i=0; i<row; i++){
    b[i]=0;
    for(int j=0; j<col; j++){
      b[i]=b[i]+X[i][j]*y[j];
    }
  }
}

//========================================================================================================================
//quadratic programming
//========================================================================================================================
void quadprog(double **Dmat, double *dvec, double *xvec, int K){ //K is the dimension of xvec
  if(K==1){
    xvec[0] = 1;
  }else{
    double **Dmat_inv, *x_star, s1, s2;
    double Dmat_inv_sum=0, x_star_sum=0, *Dmat_inv_rowsum;
    int num_negatives=0;
    double *ind_negative;
    Dmat_inv = (double **)malloc(K*sizeof(double*));
    x_star = (double *)malloc(K*sizeof(double));
    Dmat_inv_rowsum = (double *)malloc(K*sizeof(double));
    ind_negative = (double *)malloc(K*sizeof(double));
    for(int k=0; k<K; k++){
      Dmat_inv[k] = (double *)malloc(K*sizeof(double));
    }
    inverse(Dmat, K, Dmat_inv);
    for(int k=0; k<K; k++){
      s1 = 0;
      s2 = 0;
      for(int k1=0; k1<K; k1++){
        s1 += Dmat_inv[k][k1]*dvec[k1];
        Dmat_inv_sum += Dmat_inv[k][k1];
        s2 += Dmat_inv[k][k1];
      }
      x_star[k] = s1;
      Dmat_inv_rowsum[k] = s2;
      x_star_sum += s1;
    }
    for(int k=0; k<K; k++){
      xvec[k] = x_star[k] + (1-x_star_sum)/Dmat_inv_sum*Dmat_inv_rowsum[k];
      if(xvec[k]<0){
        num_negatives++;
        ind_negative[k] = 1; 
      }else{
        ind_negative[k] = 0;
      }
    }
    free(x_star);
    free(Dmat_inv_rowsum);
    for(int k=0; k<K; k++){
      free(Dmat_inv[k]);
    }
    free(Dmat_inv);
    
    if(num_negatives == 0){
      free(ind_negative);
    }else{
      int Knew = K-num_negatives, i, j;
      double ** Dmat_new, *dvec_new, *xvec_sub;
      Dmat_new = (double **)malloc((Knew)*sizeof(double*));
      for(int k=0; k<Knew; k++){
        Dmat_new[k] = (double *) malloc((Knew)*sizeof(double));
      }
      dvec_new = (double *) malloc((Knew)*sizeof(double));
      xvec_sub = (double *) malloc((Knew)*sizeof(double));
      i = 0;
      
      for(int k1=0; k1<K; k1++){
        if(ind_negative[k1]==0){
          dvec_new[i] = dvec[k1];
          j = 0;
          for(int k2=0; k2<K; k2++){
            if(ind_negative[k2]==0){
              Dmat_new[i][j] = Dmat[k1][k2];
              j++;
            }
          }
          i++;
        }else{
          xvec[k1] = 0;
        }
      }
      
      quadprog(Dmat_new, dvec_new, xvec_sub, Knew);
      i=0;
      for(int k=0; k<K; k++){
        if(ind_negative[k]==0){
          xvec[k] = xvec_sub[i];
          i++;
        }
      }
      free(dvec_new);
      free(xvec_sub);
      for(int k=0; k<Knew;k++){
        free(Dmat_new[k]);
      }
      free(Dmat_new);
    }
    
  }
  
} 

double gsl_matrix_trace(gsl_matrix *a)
{
  double tr=0;
  for (size_t i=0;i<a->size1;i++)
  {
    tr=tr+gsl_matrix_get(a,i,i);
  }
  return(tr);
}


void gsl_matrix_mul(gsl_matrix *a,gsl_matrix *b,gsl_matrix *c)
{
  for (size_t i=0;i<a->size1;i++)
  {
    for (size_t j=0;j<b->size2;j++)
    {
      double sum=0.0;
      for (size_t k=0;k<b->size1;k++)
      {				   sum+=gsl_matrix_get(a,i,k)*gsl_matrix_get(b,k,j);
      }
      gsl_matrix_set(c,i,j,sum);
    }
  }
}

void gsl_matrixvector_mul(gsl_matrix *a,gsl_vector *b,gsl_vector *c,int row,int col)
{
  double sum=0.0;
  for (int i=0;i<row;i++)
  {
    sum=0.0;
    for (int j=0;j<col;j++)
    {
      sum+=gsl_matrix_get(a,i,j)*gsl_vector_get(b,j);
      
    }
    gsl_vector_set(c,i,sum);
  }
}

void PhitoV(double w, double *t,int G,gsl_matrix *Phitov){
  double value;
  for (int i=0; i<G; i++){
    for (int j=0; j<(i+1); j++){
      value=exp((-0.5)*w*pow((t[i]-t[j]),2));
      gsl_matrix_set(Phitov, j, i, value);
      gsl_matrix_set(Phitov, i, j, value);
    }
  }
  
}

void PhitoW(double w, double v, double *t,int G,gsl_matrix *Phitow){
  double value;
  for (int i=0; i<G; i++){
    for (int j=0; j<(i+1); j++){
      value=v*exp((-0.5)*w*pow((t[i]-t[j]),2))*(-0.5)*pow(t[i]-t[j],2);
      gsl_matrix_set(Phitow, j, i, value);
      gsl_matrix_set(Phitow, i, j, value);
    }
  }
  
}

void PhiMatrix(double w, double v, double *t,int G,gsl_matrix *Phi){
  double value;
  for (int i=0; i<G; i++){
    for (int j=0; j<(i+1); j++){
      if (i==j){
        value=v+0.000001;
        gsl_matrix_set(Phi, i, i, value);
      }
      else{value=v*exp((-0.5)*w*pow((t[i]-t[j]),2));
        gsl_matrix_set(Phi, j, i, value);
        gsl_matrix_set(Phi, i, j, value);}
      
    }
  }
  
}


double get_det(gsl_matrix * A)
{
  double det=0.0; 
  int n = A->size1;
  gsl_permutation *p = gsl_permutation_calloc(n);
  gsl_matrix *tmpA = gsl_matrix_calloc(n, n);
  int signum;
  gsl_matrix_memcpy(tmpA, A);
  gsl_linalg_LU_decomp(tmpA, p, &signum);
  det = gsl_linalg_LU_det(tmpA, signum);
  gsl_permutation_free(p);
  gsl_matrix_free(tmpA);
  return det;
}

double VMV(double *t1,double *t2,double **Dmat,int a, int b)
{
  double v=0;
  for (int i=0;i<a;i++){
    for (int j=0;j<b;j++){
      v=v+t1[i]*Dmat[i][j]*t2[j];
    }
  }
  return(v);
}

void gsl_matrix_trans(gsl_matrix *a,gsl_matrix *b)
{
  for (size_t i=0;i<a->size1;i++)
  {
    for (size_t j=0;j<a->size2;j++)
    {
      gsl_matrix_set(b,j,i,gsl_matrix_get(a,i,j));
    }
  }
}

void IterativeSolver(gsl_spmatrix *A,gsl_vector *f,gsl_vector *u,int n){
  /* convert to compressed column format */
  gsl_spmatrix *C;
  C = gsl_spmatrix_ccs(A);
  
  /* now solve the system with the GMRES iterative solver */
  {
    const double tol = 1.0e-6;  /* solution relative tolerance */
  const size_t max_iter = 10; /* maximum iterations */
  const gsl_splinalg_itersolve_type *T = gsl_splinalg_itersolve_gmres;
  gsl_splinalg_itersolve *work =
    gsl_splinalg_itersolve_alloc(T, n, 0);
  size_t iter = 0;
  double residual;
  int status;
  
  /* initial guess u = 0 */
  gsl_vector_set_zero(u);
  
  /* solve the system A u = f */
  do
  {
    status = gsl_splinalg_itersolve_iterate(C, f, tol, u, work);
    
    /* print out residual norm ||A*u - f|| */
    residual = gsl_splinalg_itersolve_normr(work);
    //fprintf(stderr, "iter %zu residual = %.12e\n", iter, residual);
    
    //if (status == GSL_SUCCESS)
    // fprintf(stderr, "Converged\n");
  }
  while (status == GSL_CONTINUE && ++iter < max_iter);
  
  /* output solution */
  //for (i = 0; i < n; ++i)
  //  {
  //    double xi = (i + 1) * h;
  //    double u_exact = sin(M_PI * xi);
  //   double u_gsl = gsl_vector_get(u, i);
  
  //   printf("%f %.12e %.12e\n", xi, u_gsl, u_exact);
  //  }
  
  gsl_splinalg_itersolve_free(work);
  }
}

int main(void){
  
  //prior
  
  double alpha_w=0.001;
  double mu_w=1;
  double mu_v=0.01;
  double sigma_v=100;
  
  //int maxitertry=20;
  int nbasis=1003;
  int Q=7;
  int G=10000;
  int n=687;
  int type=6;
  int nblock=20;
  int block=500; 
  double res_old=0;
  double res_new=0;
  int it=0;
  int maxit=2000;
  
  double prior_alpha_sigma=0.0001;
  double prior_beta_sigma=0.0001;
  
  double** t=make2Darray(nblock,block);
  double** O=make2Darray(n,G);
  double** B=make2Darray(G,nbasis);
  double** X=make2Darray(n,Q);
  
  double Msigma=0.2;
  double Mw=0.5;
  double Mv=0.2;
  //double Mw_e=log(20);
  //double Mv_e=log(0.5);
  double w_old=0;
  double v_old=0;
  
  double** ESigma=make2Darray(n,G);
  double** ESigma1=make2Darray(G,G);
  //double** ESigma1_inv=make2Darray(G,G);
  double** Emu=make2Darray(n,G);
  double** L=make2Darray(n,G);
  
  double** Pest=make2Darray(n,type);
  double*** Best=make3Darray(type,nbasis,Q);
  double*** beta_est=make3Darray(type,G,Q);
  
  double value;//Esigma
  int row;
  int col;
  gsl_matrix *BLOCK = gsl_matrix_alloc(block,block);
  
  double** ThetaHat=make2Darray(nbasis,G);//C
  double** BTB=make2Darray(nbasis,nbasis);//C
  // double *yp=(double *) malloc(G*sizeof(double));
  double** y=make2Darray(n,G);
  double** residual=make2Darray(n,G);
  double** A=make2Darray(n,nbasis);
  double** UHat=make2Darray(Q*type,n);
  double** U=make2Darray(n,Q*type);
  double** BT=make2Darray(Q*type,nbasis);
  double thre=0.000001;
  
  
  double*** XCB=make3Darray(n,type,G);//p
  double** CB=make2Darray(G,Q);
  double** Dmat=make2Darray(type,type);
  double *d=(double *) malloc(type*sizeof(double));
  
  gsl_matrix *IP = gsl_matrix_alloc(block,block);//v,w
  gsl_matrix *IPI = gsl_matrix_alloc(block,block); 
  gsl_matrix *IPIO = gsl_matrix_alloc(block,block); 
  gsl_matrix *Phitov = gsl_matrix_alloc(block,block);
  gsl_matrix *Phitow = gsl_matrix_alloc(block,block);
  gsl_matrix *Phi = gsl_matrix_alloc(block,block);  
  
  gsl_matrix *omega = gsl_matrix_alloc(block,block);  
  
  gsl_vector *u = gsl_vector_alloc(block);
  gsl_vector *Su = gsl_vector_alloc(block);
  gsl_matrix *SigmaB = gsl_matrix_alloc(block,block); 
  
  double dldv1=0;
  double dldv3=0;
  double dldw1=0;
  double dldw3=0;
  double dldv2=0;
  double dldw2=0;
  
  double dldv=10;
  double dldw=10;
  
  double Mv_old=1000;
  double Mw_old=1000;
  
  double dldv_old=1000;
  double dldw_old=1000;
  
  double stepsize_w=10;
  double stepsize_v=10;
  
  double dld2;
  double *uSu = &dld2;
  
  double diff1=1;
  double diff2=1;
  int maxiter=100;
  int iter=1;
  double thr=0.0001;
  double step=0.001;
  
  /// Input
  
  FILE* fp1=fopen("t10000.txt","r"); //Ã¦â€°â€œÃ¥Â¼â‚¬Ã¦â€“â€¡Ã¤Â»Â¶
  if(fp1==NULL)
  {
    printf("no file");
    return -1;
  }
  for(int i=0;i<nblock;i++)
  {
    for(int j=0;j<block;j++)
    {
      fscanf(fp1,"%lf",&t[i][j]);/*Ã¦Â¯ÂÃ¦Â¬Â¡Ã¨Â¯Â»Ã¥Ââ€“Ã¤Â¸â‚¬Ã¤Â¸ÂªÃ¦â€¢Â°Ã¯Â¼Å’fscanfÃ¥â€¡Â½Ã¦â€¢Â°Ã©Ââ€¡Ã¥Ë†Â°Ã§Â©ÂºÃ¦Â Â¼Ã¦Ë†â€“Ã¨â‚¬â€¦Ã¦ÂÂ¢Ã¨Â¡Å’Ã§Â»â€œÃ¦ÂÅ¸*/
    }
    fscanf(fp1,"\n");
  }
  fclose(fp1);
  
  FILE* fp2=fopen("O10000.txt","r"); //Ã¦â€°â€œÃ¥Â¼â‚¬Ã¦â€“â€¡Ã¤Â»Â¶
  if(fp2==NULL)
  {
    printf("no file");
    return -1;
  }
  for(int i=0;i<n;i++)
  {
    for(int j=0;j<G;j++)
    {
      fscanf(fp2,"%lf",&O[i][j]);/*Ã¦Â¯ÂÃ¦Â¬Â¡Ã¨Â¯Â»Ã¥Ââ€“Ã¤Â¸â‚¬Ã¤Â¸ÂªÃ¦â€¢Â°Ã¯Â¼Å’fscanfÃ¥â€¡Â½Ã¦â€¢Â°Ã©Ââ€¡Ã¥Ë†Â°Ã§Â©ÂºÃ¦Â Â¼Ã¦Ë†â€“Ã¨â‚¬â€¦Ã¦ÂÂ¢Ã¨Â¡Å’Ã§Â»â€œÃ¦ÂÅ¸*/
    }
  }
  fclose(fp2);
  
  
  FILE* fp3=fopen("B10000.txt","r"); //Ã¦â€°â€œÃ¥Â¼â‚¬Ã¦â€“â€¡Ã¤Â»Â¶
  if(fp3==NULL)
  {
    printf("no file");
    return -1;
  }
  for(int i=0;i<G;i++)
  {
    for(int j=0;j<nbasis;j++)
    {
      fscanf(fp3,"%lf",&B[i][j]);/*Ã¦Â¯ÂÃ¦Â¬Â¡Ã¨Â¯Â»Ã¥Ââ€“Ã¤Â¸â‚¬Ã¤Â¸ÂªÃ¦â€¢Â°Ã¯Â¼Å’fscanfÃ¥â€¡Â½Ã¦â€¢Â°Ã©Ââ€¡Ã¥Ë†Â°Ã§Â©ÂºÃ¦Â Â¼Ã¦Ë†â€“Ã¨â‚¬â€¦Ã¦ÂÂ¢Ã¨Â¡Å’Ã§Â»â€œÃ¦ÂÅ¸*/
    }
  }
  fclose(fp3);
  
  FILE* fp4=fopen("X_vec.txt","r"); //Ã¦â€°â€œÃ¥Â¼â‚¬Ã¦â€“â€¡Ã¤Â»Â¶
  if(fp4==NULL)
  {
    printf("no file");
    return -1;
  }
  for(int i=0;i<n;i++)
  {
    for(int j=0;j<Q;j++)
    {
      fscanf(fp4,"%lf",&X[i][j]);/*Ã¦Â¯ÂÃ¦Â¬Â¡Ã¨Â¯Â»Ã¥Ââ€“Ã¤Â¸â‚¬Ã¤Â¸ÂªÃ¦â€¢Â°Ã¯Â¼Å’fscanfÃ¥â€¡Â½Ã¦â€¢Â°Ã©Ââ€¡Ã¥Ë†Â°Ã§Â©ÂºÃ¦Â Â¼Ã¦Ë†â€“Ã¨â‚¬â€¦Ã¦ÂÂ¢Ã¨Â¡Å’Ã§Â»â€œÃ¦ÂÅ¸*/
    }
  }
  fclose(fp4);
  
  FILE* fp5=fopen("BTBIBT10000.txt","r"); //Ã¦â€°â€œÃ¥Â¼â‚¬Ã¦â€“â€¡Ã¤Â»Â¶
  if(fp5==NULL)
  {
    printf("no file");
    return -1;
  }
  for(int i=0;i<nbasis;i++)
  {
    for(int j=0;j<G;j++)
    {
      fscanf(fp5,"%lf",&ThetaHat[i][j]);/*Ã¦Â¯ÂÃ¦Â¬Â¡Ã¨Â¯Â»Ã¥Ââ€“Ã¤Â¸â‚¬Ã¤Â¸ÂªÃ¦â€¢Â°Ã¯Â¼Å’fscanfÃ¥â€¡Â½Ã¦â€¢Â°Ã©Ââ€¡Ã¥Ë†Â°Ã§Â©ÂºÃ¦Â Â¼Ã¦Ë†â€“Ã¨â‚¬â€¦Ã¦ÂÂ¢Ã¨Â¡Å’Ã§Â»â€œÃ¦ÂÅ¸*/
    }
  }
  fclose(fp5);
  
  FILE* fp6=fopen("pinit_20000_hire.txt","r"); //Ã¦â€°â€œÃ¥Â¼â‚¬Ã¦â€“â€¡Ã¤Â»Â¶
  if(fp6==NULL)
  {
    printf("no file");
    return -1;
  }
  for(int i=0;i<n;i++)
  {
    for(int j=0;j<type;j++)
    {
      fscanf(fp6,"%lf",&Pest[i][j]);/*Ã¦Â¯ÂÃ¦Â¬Â¡Ã¨Â¯Â»Ã¥Ââ€“Ã¤Â¸â‚¬Ã¤Â¸ÂªÃ¦â€¢Â°Ã¯Â¼Å’fscanfÃ¥â€¡Â½Ã¦â€¢Â°Ã©Ââ€¡Ã¥Ë†Â°Ã§Â©ÂºÃ¦Â Â¼Ã¦Ë†â€“Ã¨â‚¬â€¦Ã¦ÂÂ¢Ã¨Â¡Å’Ã§Â»â€œÃ¦ÂÅ¸*/
    }
  }
  fclose(fp6);
  
  FILE* fp7=fopen("BTB10000.txt","r"); //Ã¦â€°â€œÃ¥Â¼â‚¬Ã¦â€“â€¡Ã¤Â»Â¶
  if(fp7==NULL)
  {
    printf("no file");
    return -1;
  }
  for(int i=0;i<nbasis;i++)
  {
    for(int j=0;j<nbasis;j++)
    {
      fscanf(fp7,"%lf",&BTB[i][j]);/*Ã¦Â¯ÂÃ¦Â¬Â¡Ã¨Â¯Â»Ã¥Ââ€“Ã¤Â¸â‚¬Ã¤Â¸ÂªÃ¦â€¢Â°Ã¯Â¼Å’fscanfÃ¥â€¡Â½Ã¦â€¢Â°Ã©Ââ€¡Ã¥Ë†Â°Ã§Â©ÂºÃ¦Â Â¼Ã¦Ë†â€“Ã¨â‚¬â€¦Ã¦ÂÂ¢Ã¨Â¡Å’Ã§Â»â€œÃ¦ÂÅ¸*/
    }
  }
  fclose(fp7);
  
  
  for(int i=0;i<n;i++)
  {
    for(int j=0;j<G;j++)
    {
      y[i][j]=O[i][j];
    }
  }
  
  
#pragma omp parallel for
  for (int i=0;i<n;i++){
    MatrixVector(ThetaHat,y[i],A[i],nbasis,G);
  }
  
  
  for (int i=0;i<n;i++){
    for (int ty=0;ty<type;ty++){
      for (int j=0;j<Q;j++){
        U[i][ty*Q+j]=Pest[i][ty]*X[i][j];
      }
    }
  }
  //printf("%d\n",100);
  
  hat1(U,n,type*Q,UHat);
  MultiplyMatrix(UHat,A,BT,Q*type,nbasis,n);
  
  
  for (int ty=0;ty<type;ty++){
    for (int i=0;i<nbasis;i++){
      for (int j=0;j<Q;j++){
        //if (fabs(BT[Q*ty+j][i])>0.01){
        //  Best[ty][i][j]=BT[Q*ty+j][i];}else{
        //  Best[ty][i][j]=0;
        //}
        Best[ty][i][j]=BT[Q*ty+j][i];
      }
    }
  }
  
  
  
  for (int ty=0;ty<type;ty++){
    MultiplyMatrix(B,Best[ty],CB,G,Q,nbasis);
    for (int i=0;i<n;i++){
      MatrixVector(CB,X[i],XCB[i][ty],G,Q);
    }
  }
  
  
  ////E step
  while(1){
    
    it=it+1;
    if (it>maxit){break;}
    
    res_old=res_new;
    res_new=0;
    
    
    for (int i=0;i<n;i++){
      for (int j=0;j<G;j++){
        L[i][j]=O[i][j];
        for (int ty=0;ty<type;ty++){
          L[i][j]=L[i][j]-Pest[i][ty]*XCB[i][ty][j];
        }
        res_new=res_new+L[i][j]*L[i][j];
        L[i][j]=L[i][j]/Msigma;
      }
    }
    //if (it>maxit){break;}
    res_new=res_new/n/G;
    //printf("%lf\n",res_new);
    //printf("%lf\n",res_old);
    //if ((fabs(res_old-res_new)<thre)&&(fabs(Mv-v_old)<0.001)&&(fabs(Mw-w_old)<0.001)){break;}
    if (fabs(res_old-res_new)<thre){break;}
    // if (fabs(res_new-res_old)<thre){break;}
    //Mw_old=w_old;
    //Mv_old=v_old;
    w_old=Mw;
    v_old=Mv;
    
    for (int b=0;b<nblock;b++){
      
      for (int i=0; i<block; i++){
        for (int j=0; j<(i+1); j++){
          if (i==j){
            value=Mv+0.00001;
          }
          else{value=Mv*exp((-0.5)*Mw*pow((t[b][i]-t[b][j]),2));}
          gsl_matrix_set(BLOCK, j, i, value);
          gsl_matrix_set(BLOCK, i, j, value);
        }
      }
      gsl_matrix_inv(BLOCK);
      for (int i=0; i<block; i++){
        gsl_matrix_set(BLOCK, i, i, gsl_matrix_get(BLOCK, i, i)+1/Msigma);
      }
      // for (int j1=0;j1<block;j1++){
      //    for (int j2=0;j2<block;j2++){
      //       row=b*block+j1;
      //      col=b*block+j2;
      //       ESigma1_inv[row][col]=gsl_matrix_get(BLOCK,j1,j2);
      //     }}
      gsl_matrix_inv(BLOCK);
      
      // for (int i=0;i<n;i++){
      //   for (int j=0;j<block;j++){
      //     col=b*block+j;
      //     ESigma[i][col]=gsl_matrix_get(BLOCK,j,j)+Emu[i][col]*Emu[i][col];
      //   }
      // }
      
      
      for (int j1=0;j1<block;j1++){
        for (int j2=0;j2<block;j2++){
          row=b*block+j1;
          col=b*block+j2;
          ESigma1[row][col]=gsl_matrix_get(BLOCK,j1,j2);
        }}
      
      
    }
    
    
    
#pragma omp parallel for
    for (int i=0;i<n;i++){
      MatrixVector(ESigma1,L[i],Emu[i],G,G);
    }
    
    for (int i=0;i<n;i++){
      for (int j=0;j<G;j++){
        ESigma[i][j]=ESigma1[j][j]+Emu[i][j]*Emu[i][j];
      }
    }
    
    
    
    for(int i=0;i<n;i++)
    {
      for(int j=0;j<G;j++)
      {
        y[i][j]=O[i][j]-Emu[i][j];
      }
    }
    
    
    //  for(int i=0;i<n;i++)
    //  {
    //    for(int j=0;j<G;j++)
    //    {
    //      y[i][j]=O[i][j]-Emu[i][j];
    //    }
    //  }
#pragma omp parallel for
    for (int i=0;i<n;i++){
      MatrixVector(ThetaHat,y[i],A[i],nbasis,G);
    }
    
    
    for (int i=0;i<n;i++){
      for (int ty=0;ty<type;ty++){
        for (int j=0;j<Q;j++){
          U[i][ty*Q+j]=Pest[i][ty]*X[i][j];
        }
      }
    }
    // printf("%d\n",100);
    
    hat1(U,n,type*Q,UHat);
    MultiplyMatrix(UHat,A,BT,Q*type,nbasis,n);
    
    // printf("%d\n",2);
    
    //  for (int ty=0;ty<type;ty++){
    //   for (int i=0;i<nbasis;i++){
    //     for (int j=0;j<Q;j++){
    //         Best[ty][i][j]=BT[Q*ty+j][i];
    //     }
    //   }
    // }
    
    for (int ty=0;ty<type;ty++){
      for (int i=0;i<nbasis;i++){
        for (int j=0;j<Q;j++){
          //if (fabs(BT[Q*ty+j][i])>0.01){
          //  Best[ty][i][j]=BT[Q*ty+j][i];}else{
          //  Best[ty][i][j]=0;
          //}
          Best[ty][i][j]=BT[Q*ty+j][i];
        }
      }
    }
    
    
    //    for (int i=0;i<type;i++){
    //    for (int j=0;j<nbasis;j++){
    //       for (int q=0;q<Q;q++){
    //         Best[i][j][q]=0;
    //        }
    //     }
    //  }
    
    
    // for (int i=0;i<type;i++){
    //   for (int j=0;j<nbasis;j++){
    //     for (int m=0;m<Q;m++){
    //       if (fabs(Best[i][j][m])<0.01){
    //         Best[i][j][m]=0;
    //        }
    //     }
    //    }
    // }
    
    //////M step :p
    
    
    
    for (int ty=0;ty<type;ty++){
      MultiplyMatrix(B,Best[ty],CB,G,Q,nbasis);
      for (int i=0;i<n;i++){
        MatrixVector(CB,X[i],XCB[i][ty],G,Q);
      }
    }
    
    //  for (int i=0;i<5;i++){
    //    for (int j=0;j<3;j++){
    //      printf("%lf\n",Pest[i][j]);
    //    }
    //  }
    
    //  for (int i=0;i<5;i++){
    //    for (int j=0;j<3;j++){
    //      printf("%lf\n",XCB[i][j][1]);
    //    }
    //  }
    for (int i=0;i<n;i++){
      
      for (int m1=0;m1<type;m1++){
        for (int m2=0;m2<(m1+1);m2++){
          Dmat[m1][m2]=0;
          if (m1==m2){
            for (int j=0;j<G;j++){
              Dmat[m1][m2]=Dmat[m1][m2]+XCB[i][m1][j]*XCB[i][m1][j];
            }
            // Dmat[m1][m2]=Dmat[m1][m2]/G;
            Dmat[m1][m2]=Dmat[m1][m2]+0.00001; 
          }else{
            for (int j=0;j<G;j++){
              Dmat[m1][m2]=Dmat[m1][m2]+ XCB[i][m1][j]*XCB[i][m2][j];
            }
            //  Dmat[m1][m2]=Dmat[m1][m2]/G;
            Dmat[m2][m1]=Dmat[m1][m2];
          }
        }}
      
      for (int m=0;m<type;m++){
        d[m]=0;
        for (int j=0;j<G;j++){
          d[m]=d[m]+XCB[i][m][j]*(O[i][j]-Emu[i][j]);
        }
        //d[m]=d[m]/G;
      }
      
      quadprog(Dmat,d,Pest[i],type);
      
    }
    
    
    
    // for (int j=0;j<3;j++){
    //      printf("%lf\n",d[j]);
    //   }
    
    //  for (int i=0;i<5;i++){
    //   for (int j=0;j<3;j++){
    //     printf("%lf\n",Pest[i][j]);
    //   }
    // }
    
    //  for (int i=0;i<3;i++){
    //    for (int j=0;j<3;j++){
    //      printf("%lf\n",Dmat[i][j]);
    //    }
    //  }
    
    /////M sigma
    
    for(int i=0;i<n;i++)
    {
      for(int j=0;j<G;j++)
      {
        y[i][j]=O[i][j];
      }
    }
    
    for (int ty=0;ty<type;ty++){
      for (int i=0;i<n;i++){
        for (int j=0;j<G;j++){
          y[i][j]=y[i][j]-Pest[i][ty]*XCB[i][ty][j];
        }
      }}
    Msigma=0;
    for (int i=0;i<n;i++){
      for (int j=0;j<G;j++){
        Msigma=Msigma+ESigma[i][j]-2*Emu[i][j]*y[i][j]+y[i][j]*y[i][j];
      }
    }
    
    Msigma=(Msigma+2*prior_beta_sigma)/(n*G+2*(prior_alpha_sigma+1));

    dldv1=0;
    dldv3=0;
    dldw1=0;
    dldw3=0;
    dldv2=0;
    dldw2=0;
    value=0;
    diff1=1;
    diff2=1;
    iter=0;
    Mw_old=0;
    Mv_old=0;
    dldw_old=1000;
    dldv_old=1000;

    
    while (1){
      iter=iter+1;
      dldv=0;
      dldw=0;
      
      for (int i=0;i<nblock;i++){
        
        for (int j1=0;j1<block;j1++){
          for (int j2=0;j2<block;j2++){
            row=block*i+j1;
            col=block*i+j2;
            gsl_matrix_set(SigmaB, j1, j2, ESigma1[row][col]);
          }
        }
        
        PhitoV(Mw, t[i],block,Phitov);
        PhitoW(Mw, Mv, t[i],block,Phitow);
        PhiMatrix(Mw, Mv, t[i],block,Phi);
        
        gsl_matrix_inv(Phi);
        gsl_matrix_mul(Phi,Phitov,IP);
        gsl_matrix_mul(IP,Phi,IPI);
        
        dldv2=0;
        
        for (int r=0;r<n;r++){
          
          for (int m=0;m<block;m++){
            col=i*block+m;
            gsl_vector_set(u, m, Emu[r][col]);
          }
          
          gsl_blas_dgemv(CblasNoTrans,1,IPI,u,0,Su);
          
          gsl_blas_ddot(u, Su, uSu);
          
          dldv2=dldv2+dld2;
        }
        
        gsl_matrix_mul(IPI,SigmaB,IPIO);
        
        dldv1=-0.5*gsl_matrix_trace(IP);
        dldv2=0.5*dldv2/n;
        dldv3=0.5*gsl_matrix_trace(IPIO);
        
        gsl_matrix_mul(Phi,Phitow,IP);
        gsl_matrix_mul(IP,Phi,IPI);
        gsl_matrix_mul(IPI,SigmaB,IPIO);
        
        gsl_vector_set_zero(Su);
        
        dldw2=0;
        
        for (int row=0;row<n;row++){
          
          for (int m=0;m<block;m++){
            col=i*block+m;
            gsl_vector_set(u, m, Emu[row][col]);
          }
          
          gsl_blas_dgemv(CblasNoTrans,1,IPI,u,0,Su);
          
          gsl_blas_ddot(u, Su, uSu);
          
          dldw2=dldw2+dld2;
        }
        
        
        dldw1=-0.5*gsl_matrix_trace(IP);
        dldw2=0.5*dldw2/n;
        dldw3=0.5*gsl_matrix_trace(IPIO);
        
        dldw=dldw+dldw1+dldw2+dldw3;
        dldv=dldv+dldv1+dldv2+dldv3;
      }
      dldw=dldw-1/n*(alpha_w+1)/Mw+n*alpha_w/mu_w/Mw/Mw;
      dldv=dldv-1/n*1/Mv*(1+(log(Mv)-mu_v)/sigma_v);
      stepsize_w=(Mw-Mw_old)/(dldw-dldw_old);
      stepsize_v=(Mv-Mv_old)/(dldv-dldv_old);
      //    printf("%d\n",iter);
      //  printf("%lf\n",Mw);
      //   printf("%lf\n",Mv);
      //printf("%lf\n",dldw);
      //printf("%lf\n",dldv);
      Mv_old=Mv;
      Mw_old=Mw;
      Mw=Mw-stepsize_w*dldw;
      Mv=Mv-stepsize_v*dldv;
      // Mw=Mw+stepsize_w*dldw;
      // Mv=Mv+stepsize_v*dldv;
      dldw_old=dldw;
      dldv_old=dldv;
      
      if (Mw<0.00001){Mw=0.1;}
      //if (Mw>20){Mw=2;}
      //if (Mv>0.1){Mv=Mv_old;}
      if (Mv<0.00001){Mv=0.00001;}
      if (Mv>1){Mv=0.00001;Mw=100;break;}
      if (Mw>1000){Mv=0.00001;Mw=100;break;}
      //if (Mv>0.1){Mv=0.00001;break;}
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
    printf("%d\n",it);
    printf("%lf\n",Msigma);
    printf("%lf\n",Mw);
    printf("%lf\n",Mv);
    
  }
  
  
  //printf("%d\n",100);
  
  for (int l = 0; l < n; l++)
  {//Ã¥Â¾ÂªÃ§Å½Â¯Ã¨Â¾â€œÃ¥â€¡ÂºÃ¦â€¢Â°Ã¦ÂÂ®Ã¥Â¹Â¶Ã¥â€ â„¢Ã¥â€¦Â¥
    // printf("%.3f\n", 3.1415926);
    for (int i = 0; i < G; i++){
      L[l][i]=Msigma*L[l][i];
    }}
  
  
  for (int i=0;i<type;i++){
    MultiplyMatrix(B,Best[i],beta_est[i],G,Q,nbasis);
  }
  
  FILE * pFile1;
  const char *pFileName1="XCB.txt";
  pFile1 = fopen(pFileName1, "w");
  if (NULL == pFile1)
  {
    printf("error");
    return 0;
  }
  for (int l = 0; l < n; l++){
    for (int j = 0; j < type; j++)
    {
      // printf("%.3f\n", 3.1415926);
      for (int i = 0; i < G; i++){
        fprintf(pFile1,"%.3f\n", XCB[l][j][i]);/
      }}}
  fclose(pFile1);
  
  FILE * pFile2;
  const char *pFileName2="L.txt";
  pFile2 = fopen(pFileName2, "w");
  if (NULL == pFile2)
  {
    printf("error");
    return 0;
  }
  for (int l = 0; l < n; l++)
  {
    // printf("%.3f\n", 3.1415926);
    for (int i = 0; i < G; i++){
      fprintf(pFile2,"%.3f\n", L[l][i]);
    }}
  fclose(pFile2);
  
  
  FILE * pFile3;
  const char *pFileName3="Emu.txt";
  pFile3 = fopen(pFileName3, "w");
  if (NULL == pFile3)
  {
    printf("error");
    return 0;
  }
  for (int l = 0; l < n; l++)
  {
    // printf("%.3f\n", 3.1415926);
    for (int i = 0; i < G; i++){
      fprintf(pFile3,"%.3f\n", Emu[l][i]);
    }}
  fclose(pFile3);
  
  
  FILE * pFile4;
  const char *pFileName4="Esigma.txt";
  pFile4 = fopen(pFileName4, "w");
  if (NULL == pFile4)
  {
    printf("error");
    return 0;
  }
  for (int l = 0; l < n; l++)
  {
    // printf("%.3f\n", 3.1415926);
    for (int i = 0; i < G; i++){
      fprintf(pFile4,"%.3f\n", ESigma[l][i]);
    }}
  fclose(pFile4);
  
  
  FILE * pFile5;
  const char *pFileName5="Pest.txt";
  pFile5 = fopen(pFileName5, "w");
  if (NULL == pFile5)
  {
    printf("error");
    return 0;
  }
  for (int l = 0; l < n; l++)
  {
    // printf("%.3f\n", 3.1415926);
    for (int i = 0; i < type; i++){
      fprintf(pFile5,"%.3f\n", Pest[l][i]);
    }}
  fclose(pFile5);
  
  
  FILE * pFile01;
  const char *pFileName01="beta1.txt";
  pFile01 = fopen(pFileName01, "w");
  if (NULL == pFile01)
  {
    printf("error");
    return 0;
  }
  for (int i = 0; i < Q; i++)
  {
    // printf("%.3f\n", 3.1415926);
    for (int ty = 0; ty < G; ty++){
      
      fprintf(pFile01,"%.3f\n", beta_est[0][ty][i]);
    }}
  fclose(pFile01);
  
  
  
  FILE * pFile02;
  const char *pFileName02="beta2.txt";
  pFile02 = fopen(pFileName02, "w");
  if (NULL == pFile02)
  {
    printf("error");
    return 0;
  }
  for (int i = 0; i < Q; i++)
  {
    // printf("%.3f\n", 3.1415926);
    for (int ty = 0; ty < G; ty++){
      
      fprintf(pFile02,"%.3f\n", beta_est[1][ty][i]);
    }}
  fclose(pFile02);
  
  
  
  FILE * pFile03;
  const char *pFileName03="beta3.txt";
  pFile03 = fopen(pFileName03, "w");
  if (NULL == pFile03)
  {
    printf("error");
    return 0;
  }
  for (int i = 0; i < Q; i++)
  {
    
    for (int ty = 0; ty < G; ty++){
      
      fprintf(pFile03,"%.3f\n", beta_est[2][ty][i]);
    }}
  fclose(pFile03);
  
  
  
  FILE * pFile04;
  const char *pFileName04="beta4.txt";
  pFile04 = fopen(pFileName04, "w");
  if (NULL == pFile04)
  {
    printf("error");
    return 0;
  }
  for (int i = 0; i < Q; i++)
  {
   
    for (int ty = 0; ty < G; ty++){
      
      fprintf(pFile04,"%.3f\n", beta_est[3][ty][i]);
    }}
  fclose(pFile04);
  
  
  
  FILE * pFile05;
  const char *pFileName05="beta5.txt";
  pFile05 = fopen(pFileName05, "w");
  if (NULL == pFile05)
  {
    printf("error");
    return 0;
  }
  for (int i = 0; i < Q; i++)
  {
  
    for (int ty = 0; ty < G; ty++){
      
      fprintf(pFile05,"%.3f\n", beta_est[4][ty][i]);
    }}
  fclose(pFile05);
  
  FILE * pFile06;
  const char *pFileName06="beta6.txt";
  pFile06 = fopen(pFileName06, "w");
  if (NULL == pFile06)
  {
    printf("error");
    return 0;
  }
  for (int i = 0; i < Q; i++)
  {
   
    for (int ty = 0; ty < G; ty++){
      
      fprintf(pFile06,"%.3f\n", beta_est[5][ty][i]);
    }}
  fclose(pFile06);
  
  delet2Darray(Dmat,type,type);
  delet2Darray(CB,G,Q);
  delet3Darray(XCB,n,type,G);
  delet2Darray(ThetaHat,nbasis,G);//C
  delet2Darray(y,n,G);//C
  delet2Darray(A,n,nbasis);//C
  delet2Darray(UHat,Q*type,n);//C
  delet2Darray(U,n,Q*type);//C
  delet2Darray(BT,Q*type,nbasis);//C
  delet2Darray(ESigma,n,G);//C
  delet2Darray(ESigma1,G,G);//C
  delet2Darray(Emu,n,G);//C
  
  
  gsl_matrix *X_lm = gsl_matrix_alloc(n,type*Q);
  gsl_matrix *X_lm_t = gsl_matrix_alloc(type*Q,n);
  gsl_vector *y_lm = gsl_vector_alloc(n);
  gsl_matrix *XTX_lm = gsl_matrix_alloc(type*Q,type*Q);
  gsl_matrix *XTXIXT_lm = gsl_matrix_alloc(type*Q,n);
  gsl_vector *beta_lm = gsl_vector_alloc(type*Q);
  gsl_vector *yest_lm = gsl_vector_alloc(n);
  double*** beta_sd=make3Darray(type,G,Q);
  double*** beta_z=make3Darray(type,G,Q);
  
  for (int ty=0;ty<type;ty++){
    for (int j=0;j<Q;j++){
      for (int i=0;i<n;i++){
      col=ty*Q+j;
      gsl_matrix_set(X_lm, i, col, X[i][j]*Pest[i][ty]);
      }
    }
  }
  
  gsl_matrix_trans(X_lm,X_lm_t);
  gsl_matrix_mul(X_lm_t, X_lm, XTX_lm);
  gsl_matrix_inv(XTX_lm);
  gsl_matrix_mul(XTX_lm, X_lm_t, XTXIXT_lm);
  
  for (int g=0;g<G;g++){
    for (int i=0;i<n;i++){
    gsl_vector_set(y_lm, i, O[i][g]);
      
    }
    gsl_matrixvector_mul(XTXIXT_lm,y_lm,beta_lm,type*Q,n);
    gsl_matrixvector_mul(X_lm,beta_lm,yest_lm,n,type*Q);
    value=0;
    for (int i=0;i<n;i++){
      value=value+pow(gsl_vector_get(y_lm, i)-gsl_vector_get(yest_lm, i),2);
      
    }
    value=value/(n-1);
    for (int ty=0;ty<type;ty++){
      for (int j=0;j<Q;j++){
        col=ty*Q+j;
        beta_sd[ty][g][j]=sqrt(gsl_matrix_get(XTX_lm, col, col)*value);
        beta_z[ty][g][j]=beta_est[ty][g][j]/beta_sd[ty][g][j];
      }}
  }
  
  
  
  
  FILE * pFile13;
  const char *pFileName13="beta_z.txt";
  pFile13 = fopen(pFileName13, "w");
  if (NULL == pFile13)
  {
    printf("error");
    return 0;
  }
  for (int i = 0; i < G; i++){
    for (int ty = 0; ty < type; ty++){
      for (int q = 1; q < Q; q++){
        fprintf(pFile13,"%.3f\n", beta_z[ty][i][q]);
      }
    }
  }
  fclose(pFile13);
  
  
  FILE * pFile14;
  const char *pFileName14="beta_sd.txt";
  pFile14 = fopen(pFileName14, "w");
  if (NULL == pFile14)
  {
    printf("error");
    return 0;
  }
  for (int i = 0; i < G; i++){
    for (int ty = 0; ty < type; ty++){
      for (int q = 1; q < Q; q++){
        fprintf(pFile14,"%.3f\n", beta_sd[ty][i][q]);
      }
    }
  }
  fclose(pFile14);
  
}


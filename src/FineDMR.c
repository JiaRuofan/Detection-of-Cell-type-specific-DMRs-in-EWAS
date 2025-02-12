
//NR method
//Current
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

double** C2darray_diy2(int a1, int *a2){
  double **tmp;
  tmp = (double **)malloc(a1*sizeof(double *));
  for(int i=0; i<a1; i++){
    tmp[i] = (double *) malloc(a2[i]*sizeof(double));
  }
  return tmp;
}

void delet2darray_diy2(double **tmp, int a1, int a2){
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

double*** C3darray_diy2(int a1, int *a2, int a3){
  double ***tmp;
  tmp = (double ***)malloc(a1*sizeof(double **));
  for(int i=0; i<a1; i++){
    tmp[i] = (double **) malloc(a2[i]*sizeof(double *));
    for(int j=0; j<a2[i]; j++){
      tmp[i][j] = (double *) malloc(a3*sizeof(double));
    }
  }
  return tmp;
}

void delet3darray_diy2(double ***tmp, int a1, int *a2, int a3){
  for(int i=0; i<a1; i++){
    for(int j=0; j<a2[i]; j++){
      free(tmp[i][j]);
    }
    free(tmp[i]);
  }
  free(tmp);
}

double*** C3darray_diy3(int a1, int a2, int *a3){
  double ***tmp;
  tmp = (double ***)malloc(a1*sizeof(double **));
  for(int i=0; i<a1; i++){
    tmp[i] = (double **) malloc(a2*sizeof(double *));
    for(int j=0; j<a2; j++){
      tmp[i][j] = (double *) malloc(a3[i]*sizeof(double));
    }
  }
  return tmp;
}

void delet3darray_diy3(double ***tmp, int a1, int a2, int *a3){
  for(int i=0; i<a1; i++){
    for(int j=0; j<a2; j++){
      free(tmp[i][j]);
    }
    free(tmp[i]);
  }
  free(tmp);
}

double*** C3darray_diy23(int a1, int *a2, int *a3){
  double ***tmp;
  tmp = (double ***)malloc(a1*sizeof(double **));
  for(int i=0; i<a1; i++){
    tmp[i] = (double **) malloc(a2[i]*sizeof(double *));
    for(int j=0; j<a2[i]; j++){
      tmp[i][j] = (double *) malloc(a3[i]*sizeof(double));
    }
  }
  return tmp;
}

void delet3darray_diy23(double ***tmp, int a1, int *a2, int *a3){
  for(int i=0; i<a1; i++){
    for(int j=0; j<a2[i]; j++){
      free(tmp[i][j]);
    }
    free(tmp[i]);
  }
  free(tmp);
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
      double su=0.0;
      for (size_t k=0;k<b->size1;k++)
      {				   su+=gsl_matrix_get(a,i,k)*gsl_matrix_get(b,k,j);
      }
      gsl_matrix_set(c,i,j,su);
    }
  }
}

void gsl_matrix_mul_pa(gsl_matrix *a,gsl_matrix *b,gsl_matrix *c)
{
  
  double** ab=make2Darray(a->size1,b->size2);
  //#pragma omp parallel for
  for (size_t i=0;i<a->size1;i++)
  {
    for (size_t j=0;j<b->size2;j++)
    {
      ab[i][j]=0;
      for (size_t k=0;k<b->size1;k++)
      {				   ab[i][j]+=gsl_matrix_get(a,i,k)*gsl_matrix_get(b,k,j);
      }
      gsl_matrix_set(c,i,j,ab[i][j]);
    }
  }
  delet2Darray(ab,a->size1,b->size2);
  
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

double gsl_VMV(gsl_vector *b,gsl_matrix *a,gsl_vector *c,int row,int col)
{
  double v=0;
  for (int i=0;i<row;i++){
    for (int j=0;j<col;j++){
      v=v+gsl_vector_get(b,i)*gsl_matrix_get(a,i,j)*gsl_vector_get(c,j);
    }
  }
  return(v);
}


int main(void){
  
  //prior
  double Ll=0;
  double Ll_old;
  double Ll_pre=100000;
  double alpha_w=0.001;
  double mu_w=1;
  double mu_v=0.01;
  double sigma_v=100;
  
  //int maxitertry=20;
  int nblock=50;
  int *block=(int *) malloc(nblock*sizeof(int));
  int *nbasis=(int *) malloc(nblock*sizeof(int));
  
  FILE* fp8=fopen("block.txt","r"); //æ‰“å¼€æ–‡ä»¶
  if(fp8==NULL)
  {
    printf("no file");
    return -1;
  }
  for(int j=0;j<nblock;j++)
  {
    fscanf(fp8,"%d",&block[j]);/*æ¯æ¬¡è¯»å–ä¸€ä¸ªæ•°ï¼Œfscanfå‡½æ•°é‡åˆ°ç©ºæ ¼æˆ–è€…æ¢è¡Œç»“æŸ*/
  }
  
  fclose(fp8);
  
  FILE* fp9=fopen("nbasis.txt","r"); //æ‰“å¼€æ–‡ä»¶
  if(fp9==NULL)
  {
    printf("no file");
    return -1;
  }
  for(int j=0;j<nblock;j++)
  {
    fscanf(fp9,"%d",&nbasis[j]);/*æ¯æ¬¡è¯»å–ä¸€ä¸ªæ•°ï¼Œfscanfå‡½æ•°é‡åˆ°ç©ºæ ¼æˆ–è€…æ¢è¡Œç»“æŸ*/
  }
  
  fclose(fp9);
  
  
  
  
  int Q=2;
  int G=10000;
  int n=600;
  int type=6;
  
  int n_batch=n/10;
  int batch;
  
  double BLOCK_del; 
  
  double res_old=0;
  double res_new=0;
  int it=0;
  int maxit=10;
  
  double prior_alpha_sigma=0.0001;
  //double prior_beta_sigma=10000;
  double prior_beta_sigma=1;
  
  double** t=C2darray_diy2(nblock,block);
  double*** O=C3darray_diy3(nblock,n,block);
  double*** B=C3darray_diy23(nblock,block,nbasis);
  double** X=make2Darray(n,Q);
  
  double Msigma=0.0001;
  double Msigma_old;
  double Mw=2;
  double Mv=0.01;
  //double Mw_e=log(20);
  //double Mv_e=log(0.5);
  double w_old=0;
  double v_old=0;
  
  //double*** ESigma=C3darray_diy3(nblock,n,block);
  //double*** ESigma1=C3darray_diy23(nblock,block,block);
  
  
  
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
  
  
  double ***Emu;
  Emu = (double ***)malloc(nblock*sizeof(double **));
  for(int i=0; i<nblock; i++){
    Emu[i] = (double **) malloc(n*sizeof(double *));
    for(int j=0; j<n; j++){
      Emu[i][j] = (double *) malloc(type*block[i]*sizeof(double));
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
  
  double** Pest=make2Darray(n,type);
  double** Pest_old=make2Darray(n,type);
  double** Pest_init=make2Darray(n,type);
  //double*** Best=make3Darray(type,nbasis,Q);
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
  
  //double*** beta_est=make3Darray(type,G,Q);
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
  
  double value;//Esigma
  int row;
  int col;
  
  
  double*** ThetaHat=C3darray_diy23(nblock,nbasis,block);//C
  // double *yp=(double *) malloc(G*sizeof(double));
  double*** y=C3darray_diy3(nblock,n,block);
  double** residual=make2Darray(n,G);
  double*** A=C3darray_diy3(nblock,n,nbasis);
  double** UHat=make2Darray(Q,n);
  double** UHatT=make2Darray(Q*type,n);
  double** U=make2Darray(n,Q*type);
  double*** BT=C3darray_diy3(nblock,Q,nbasis);
  double*** BTT=C3darray_diy3(nblock,Q*type,nbasis);
  double thre=0.000001;
  
  
  //double*** XCB=make3Darray(n,type,G);//p
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
  //double** Dmat=make2Darray(type,type);
  //double *d=(double *) malloc(type*sizeof(double));
  double *sigma_vec=(double *) malloc(n*sizeof(double));
  double*** Dmat=make3Darray(n,type,type);
  double** d=make2Darray(n,type);
  
  
  
  // double dldv1=0;
  // double dldv3=0;
  // double dldw1=0;
  // double dldw3=0;
  // double dldv2=0;
  // double dldw2=0;
  
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
  
  
  double dldv_old=100000;
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
  
  FILE* fp1=fopen("t.txt","r"); //æ‰“å¼€æ–‡ä»¶
  if(fp1==NULL)
  {
    printf("no file");
    return -1;
  }
  for(int i=0;i<nblock;i++)
  {
    for(int j=0;j<block[i];j++)
    {
      fscanf(fp1,"%lf",&t[i][j]);/*æ¯æ¬¡è¯»å–ä¸€ä¸ªæ•°ï¼Œfscanfå‡½æ•°é‡åˆ°ç©ºæ ¼æˆ–è€…æ¢è¡Œç»“æŸ*/
    }
    
  }
  fclose(fp1);
  
  char Ofilename[10];
  FILE *fpSet[nblock];
  for (int j=0;j<(nblock);j++){
    sprintf(Ofilename,"O%d.txt",(j+1));
    fpSet[j]=fopen(Ofilename,"r");  
    if(!fpSet[j])
    {
      printf("no file");
      return -1;
    }
    
    for(int i=0;i<n;i++)
    {
      for(int l=0;l<block[j];l++)
      {
        fscanf(fpSet[j],"%lf",&O[j][i][l]);/*æ¯æ¬¡è¯»å–ä¸€ä¸ªæ•°ï¼Œfscanfå‡½æ•°é‡åˆ°ç©ºæ ¼æˆ–è€…æ¢è¡Œç»“æŸ*/
      }
    }
    fclose(fpSet[j]);
  }
  
  
  FILE* fp3=fopen("B.txt","r"); //æ‰“å¼€æ–‡ä»¶
  if(fp3==NULL)
  {
    printf("no file");
    return -1;
  }
  for (int l=0;l<nblock;l++){
    for(int i=0;i<block[l];i++)
    {
      for(int j=0;j<nbasis[l];j++)
      {
        fscanf(fp3,"%lf",&B[l][i][j]);/*æ¯æ¬¡è¯»å–ä¸€ä¸ªæ•°ï¼Œfscanfå‡½æ•°é‡åˆ°ç©ºæ ¼æˆ–è€…æ¢è¡Œç»“æŸ*/
      }
    }}
  fclose(fp3);
  
  
  FILE* fp4=fopen("X.txt","r"); //æ‰“å¼€æ–‡ä»¶
  if(fp4==NULL)
  {
    printf("no file");
    return -1;
  }
  for(int i=0;i<n;i++)
  {
    for(int j=0;j<Q;j++)
    {
      fscanf(fp4,"%lf",&X[i][j]);/*æ¯æ¬¡è¯»å–ä¸€ä¸ªæ•°ï¼Œfscanfå‡½æ•°é‡åˆ°ç©ºæ ¼æˆ–è€…æ¢è¡Œç»“æŸ*/
    }
  }
  fclose(fp4);
  
  FILE* fp5=fopen("BTBIBT.txt","r"); //æ‰“å¼€æ–‡ä»¶
  if(fp5==NULL)
  {
    printf("no file");
    return -1;
  }
  for(int l=0;l<nblock;l++){
    for(int i=0;i<nbasis[l];i++)
    {
      for(int j=0;j<block[l];j++)
      {
        fscanf(fp5,"%lf",&ThetaHat[l][i][j]);/*æ¯æ¬¡è¯»å–ä¸€ä¸ªæ•°ï¼Œfscanfå‡½æ•°é‡åˆ°ç©ºæ ¼æˆ–è€…æ¢è¡Œç»“æŸ*/
      }
    }}
  fclose(fp5);
  
  FILE* fp6=fopen("pinit.txt","r"); //æ‰“å¼€æ–‡ä»¶
  if(fp6==NULL)
  {
    printf("no file");
    return -1;
  }
  for(int i=0;i<n;i++)
  {
    for(int j=0;j<type;j++)
    {
      fscanf(fp6,"%lf",&Pest[i][j]);/*æ¯æ¬¡è¯»å–ä¸€ä¸ªæ•°ï¼Œfscanfå‡½æ•°é‡åˆ°ç©ºæ ¼æˆ–è€…æ¢è¡Œç»“æŸ*/
    }
  }
  fclose(fp6);
  
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
  
  
  //for (int l=0;l<nblock;l++){
  
  //   for (int i=0;i<nbasis[l];i++){
  //   for (int j=0;j<Q;j++){
  //      //if (fabs(BT[Q*ty+j][i])>0.01){
  //  Best[ty][i][j]=BT[Q*ty+j][i];}else{
  //  Best[ty][i][j]=0;
  //}
  //    Best[l][0][i][j]=0.1;
  //   Best[l][1][i][j]=0.9;
  //   Best[l][2][i][j]=0.8;
  //   Best[l][3][i][j]=0.2;
  //   Best[l][4][i][j]=0.3;
  //   Best[l][5][i][j]=0.7;
  // }
  // }
  //  }
  
  
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
  
  
  
  //  for (int b=0;b<nblock;b++){
  
  
  
  
  //   MultiplyMatrix(UHatT,A[b],BTT[b],Q*type,nbasis[b],n);
  
  
  
  //   for (int i=0;i<nbasis[b];i++){
  //     for (int j=0;j<Q;j++){
  
  //    Best[b][0][i][0]=0.1;
  //     Best[b][1][i][0]=0.9;
  //    Best[b][2][i][0]=0.2;
  //    Best[b][3][i][0]=0.8;
  //    Best[b][4][i][0]=0.3;
  //    Best[b][5][i][0]=0.7;
  //  }
  //    }
  //
  
  
  // }  
  
  
  
  // for (int j=0;j<nblock;j++){
  //  for (int i=0;i<type;i++){
  //    MultiplyMatrix(B[j],Best[j][i],beta_est[j][i],block[j],Q,nbasis[j]);
  // }}
  
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
  
  
  double ***PiTB;
  PiTB = (double ***)malloc(n*sizeof(double **));
  for(int i=0; i<n; i++){
    PiTB[i] = (double **) malloc(1200*sizeof(double *));
    for(int j=0; j<1200; j++){
      PiTB[i][j] = (double *) malloc(200*sizeof(double));
    }
  }
  
  
  //E step
  
  while(1){
    
    it=it+1;
    if (it>maxit){break;}
    //if (abs(Ll_pre-Ll)<0.01){
    //  break;
    // }
    
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
    
    
    FILE * pFile1;
    const char *pFileName1="XCB.txt";
    pFile1 = fopen(pFileName1, "w");
    for (int ty = 0; ty < nblock; ty++){
      for (int l = 0; l < n; l++){
        for (int j = 0; j < type; j++)
        {
          
          for (int i = 0; i < block[ty]; i++){
            fprintf(pFile1,"%.3f\n", XCB[ty][l][j][i]);
          }}}}
    fclose(pFile1);
    
    
    
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
    
    //Msigma=0;
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
        //sigma_vec[i]=0;
        
        
        
        for(int k=0;k<type;k++){
          for(int j=0;j<block[b];j++){
            
            gsl_matrix_set(Pi, j, j*type+k, Pest[i][k]);
            gsl_matrix_set(Pi_T,  j*type+k,j, Pest[i][k]);
            // if (abs(Pest[i][k])>1){
            //   printf("%lf\n",Pest[i][k]);
            // }
            //printf("%lf\n",gsl_matrix_get(Pi_T,j*type+k,j));
          }
        }
        
        //gsl_matrix_trans(Pi,Pi_T);
        
        for (int j1=0;j1<(type*block[b]);j1++){
          for (int j2=0;j2<block[b];j2++){
            PiTB[i][j1][j2]=gsl_matrix_get(Pi_T, j1, j2);
            if (abs(PiTB[i][j1][j2])>1){
              printf("%lf\n",PiTB[i][j1][j2]);
            }
          }
        }
        
        
        gsl_matrix_mul(Pi_T, BLOCK, Pi_TB);
        
        
        
        
        gsl_matrix_mul(Pi_TB, Pi, PBP);
        
        //  for (int j1=0;j1<(type*block[b]);j1++){
        //    for (int j2=0;j2<block[b];j2++){
        //    PiTB[i][j1][j2]=gsl_matrix_get(Pi_TB, j1, j2);
        //    }
        //  }
        // for (int j1=0;j1<(type*block[b]);j1++){
        //   for (int j2=0;j2<(type*block[b]);j2++){
        //     vv[i]=0;
        //     for (int j3=0;j3<(block[b]);j3++){
        //      vv[i]=vv[i]+gsl_matrix_get(Pi_TB, j1, j3)*gsl_matrix_get(Pi, j3, j2);
        
        //    }
        //    gsl_matrix_set(PBP, j1, j2, vv[i]);
        //  }
        
        // }
        
        
        for(int j=0;j<(type*block[b]);j++){
          
          vv[i]=gsl_matrix_get(PBP, j, j);
          gsl_matrix_set(PBP, j, j, vv[i]+1/Msigma);
          
        }
        
        gsl_matrix_inv(PBP);
        
        
        
        //  FILE * pFile0100;
        //   const char *pFileName0100="PBP.txt";
        //   pFile0100 = fopen(pFileName0100, "w");
        
        //   for (int j1 = 0; j1 < (type*block[b]); j1++)
        //   {
        //    for (int j2 = 0; j2 < (type*block[b]); j2++)
        //    {  
        
        
        //        fprintf(pFile0100,"%.3f\n", gsl_matrix_get(PBP, j1, j2));
        //       }}
        //   fclose(pFile0100);
        
        
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
        
        
        
        
        //   for (int j2=1100;j2<1200;j2++){
        //    printf("%lf\n",gsl_vector_get(Emu_sub,j2));
        
        //  }
        
        for (int j1=0;j1<(type);j1++){
          for (int j2=0;j2<block[b];j2++){
            Emu[b][i][type*j2+j1]=gsl_vector_get(Emu_sub_nor,type*j2+j1);
            //printf("%lf\n",gsl_vector_get(Emu_sub_nor,type*j2+j1));
            //printf("%lf\n",gsl_vector_get(Emu_sub_nor,type*j2+j1));
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
            //printf("%lf\n",gsl_matrix_trace(BS));
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
          
          //gsl_matrixvector_mul(BLOCK,mu_k1,Bmu_k2,block[b],block[b]);
          
          for (int i1=0; i1<block[b]; i1++){
            
            d[i][m]=d[i][m]+gsl_vector_get(Bmu_k2, i1)*O[b][i][i1];
            
          }
          
          //printf("%lf\n",d[i][m]);
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
      //printf("%lf\n",Ll);
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
        Ll=Ll-0.5*BLOCK_del;
        
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
      //printf("%lf\n",Ll);
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
        Ll_old=Ll_old-0.5*BLOCK_del;
        
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
    
    
    
    // if (Ll_old>Ll){
    
    //   for (int j=0;j<n;j++){
    // for (int k=0;k<type;k++){
    //    Pest[j][k]=Pest_old[j][k];
    //  }
    //  }
    
    //   Ll=Ll_old;
    // }
    
    
    //Ll=Ll/n/G;
    
    
    
    //   for (int j = 0; j < nblock; j++)
    //  {
    //   for (int l = 0; l < n; l++)
    //    {
    
    //    for (int i = 0; i < (type*block[j]); i++){
    //     if (Emu[j][l][i]<0){
    //      Emu[j][l][i]=0;
    //     }
    //   }}}
    
    
    //Ck
    
    for (int k=0;k<type;k++){
      
      //#pragma omp parallel for
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
      for (int j = 0; j < nblock; j++)
      {  
        for (int ty = 0; ty < block[j]; ty++){
          
          fprintf(pFile01,"%.3f\n", beta_est[j][0][ty][i]);
        }}}
    fclose(pFile01);
    
    
    
    FILE * pFile02;
    const char *pFileName02="beta2.txt";
    pFile02 = fopen(pFileName02, "w");//è¿™ä¸ªç”¨â€œwâ€æ˜¯å†™æ–‡ä»¶ï¼Œè¦†ç›–åŽŸå†…å®¹ï¼Œè‹¥ä¸æƒ³è¦†ç›–åˆ™ç”¨â€œaâ€
    if (NULL == pFile02)
    {//æ–‡ä»¶æ‰“å¼€é”™è¯¯
      printf("error");
      return 0;
    }
    for (int i = 0; i < Q; i++)
    {
      for (int j = 0; j < nblock; j++)
      {  
        for (int ty = 0; ty < block[j]; ty++){
          
          fprintf(pFile01,"%.3f\n", beta_est[j][1][ty][i]);
        }}}
    fclose(pFile02);//æœ€åŽä¸€å®šè®°å¾—å…³é—­æ–‡ä»¶
    
    
    
    FILE * pFile03;
    const char *pFileName03="beta3.txt";
    pFile03 = fopen(pFileName03, "w");//è¿™ä¸ªç”¨â€œwâ€æ˜¯å†™æ–‡ä»¶ï¼Œè¦†ç›–åŽŸå†…å®¹ï¼Œè‹¥ä¸æƒ³è¦†ç›–åˆ™ç”¨â€œaâ€
    if (NULL == pFile03)
    {//æ–‡ä»¶æ‰“å¼€é”™è¯¯
      printf("error");
      return 0;
    }
    for (int i = 0; i < Q; i++)
    {
      for (int j = 0; j < nblock; j++)
      {  
        for (int ty = 0; ty < block[j]; ty++){
          
          fprintf(pFile01,"%.3f\n", beta_est[j][2][ty][i]);
        }}}
    fclose(pFile03);//æœ€åŽä¸€å®šè®°å¾—å…³é—­æ–‡ä»¶
    
    
    
    FILE * pFile04;
    const char *pFileName04="beta4.txt";
    pFile04 = fopen(pFileName04, "w");//è¿™ä¸ªç”¨â€œwâ€æ˜¯å†™æ–‡ä»¶ï¼Œè¦†ç›–åŽŸå†…å®¹ï¼Œè‹¥ä¸æƒ³è¦†ç›–åˆ™ç”¨â€œaâ€
    if (NULL == pFile04)
    {//æ–‡ä»¶æ‰“å¼€é”™è¯¯
      printf("error");
      return 0;
    }
    for (int i = 0; i < Q; i++)
    {
      for (int j = 0; j < nblock; j++)
      {  
        for (int ty = 0; ty < block[j]; ty++){
          
          fprintf(pFile01,"%.3f\n", beta_est[j][3][ty][i]);
        }}}
    fclose(pFile04);//æœ€åŽä¸€å®šè®°å¾—å…³é—­æ–‡ä»¶
    
    
    
    FILE * pFile05;
    const char *pFileName05="beta5.txt";
    pFile05 = fopen(pFileName05, "w");//è¿™ä¸ªç”¨â€œwâ€æ˜¯å†™æ–‡ä»¶ï¼Œè¦†ç›–åŽŸå†…å®¹ï¼Œè‹¥ä¸æƒ³è¦†ç›–åˆ™ç”¨â€œaâ€
    if (NULL == pFile05)
    {//æ–‡ä»¶æ‰“å¼€é”™è¯¯
      printf("error");
      return 0;
    }
    for (int i = 0; i < Q; i++)
    {
      for (int j = 0; j < nblock; j++)
      {  
        for (int ty = 0; ty < block[j]; ty++){
          
          fprintf(pFile01,"%.3f\n", beta_est[j][4][ty][i]);
        }}}
    fclose(pFile05);//æœ€åŽä¸€å®šè®°å¾—å…³é—­æ–‡ä»¶
    
    FILE * pFile06;
    const char *pFileName06="beta6.txt";
    pFile06 = fopen(pFileName06, "w");//è¿™ä¸ªç”¨â€œwâ€æ˜¯å†™æ–‡ä»¶ï¼Œè¦†ç›–åŽŸå†…å®¹ï¼Œè‹¥ä¸æƒ³è¦†ç›–åˆ™ç”¨â€œaâ€
    if (NULL == pFile06)
    {//æ–‡ä»¶æ‰“å¼€é”™è¯¯
      printf("error");
      return 0;
    }
    for (int i = 0; i < Q; i++)
    {
      for (int j = 0; j < nblock; j++)
      {  
        for (int ty = 0; ty < block[j]; ty++){
          
          fprintf(pFile01,"%.3f\n", beta_est[j][5][ty][i]);
        }}}
    fclose(pFile06);
    
    
    
    //epsilon
    
    Msigma=0;
    for (int i=0;i<n;i++){
      
      Msigma=Msigma+sigma_vec[i];
    }
    
    
    printf("%lf\n",Msigma);
    
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
            //printf("%lf\n",XCB[l][i][k][j]*XCB[l][i][k][j]-2*Emu[l][i][j*type+k]*XCB[l][i][k][j]);
          }
        }}}
    
    Msigma=(Msigma+2*prior_beta_sigma)/(n*G*type+2*(prior_alpha_sigma+1));    
    
    
    //v,w
    
    
    value=0;
    diff1=1;
    diff2=1;
    iter=0;
    
    //Mw=5;
    //Mv=0.03;
    // for (int i=0;i<3;i++){
    //  for (int j=0;j<3;j++){
    //      printf("%lf\n",ESigma1[i][j]);
    //    }
    //  }
    
    batch=it%10;
    
    
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
      
      
      //#pragma omp parallel for
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
        //printf("%d\n",i);
      }
      
      
      dldv=0;
      dldw=0;
      for (int i=0;i<nblock;i++){
        dldw=dldw+dldw_vec[i];
        dldv=dldv+dldv_vec[i];
      }
      dldw=dldw-1/n_batch*(alpha_w+1)/Mw+n_batch*alpha_w/mu_w/Mw/Mw;
      dldv=dldv-1/n_batch*1/Mv*(1+(log(Mv)-mu_v)/sigma_v);
      printf("%d\n",iter);
      printf("%lf\n",dldw);
      printf("%lf\n",dldv);
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
      
      if (Mw<0.00001){Mw=Mw_old;}
      //if (Mw>20){Mw=2;}
      //if (Mv>0.1){Mv=Mv_old;}
      if (Mv<0.00001){Mv=Mv_old;}
      if (Mv>1){Mv=0.1;Mw=Mw_old;break;}
      if (Mw>1000){Mv=Mv_old;Mw=100;break;}
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
    
    
    
    
    
    
    //   Ll_init=0;
    //   for (int b=0;b<nblock;b++){
    //printf("%lf\n",Ll);
    //     gsl_matrix *BLOCK = gsl_matrix_alloc(block[b],block[b]);
    //     gsl_vector *x_sub = gsl_vector_alloc(block[b]);
    
    
    //     for (int i=0;i<n;i++){
    
    //       for (int ii=0; ii<block[b]; ii++){
    
    //        for (int j=0; j<(ii+1); j++){
    
    //        if (ii==j){
    //          value=0;
    //          for (int ty=0;ty<type;ty++){
    //            value=Mv+Pest_init[i][ty]*Pest_init[i][ty]*Msigma;
    //         }
    
    //       }
    //       else{value=Mv*exp((-0.5)*Mw*pow((t[b][ii]-t[b][j]),2));}
    //       gsl_matrix_set(BLOCK, j, ii, value);
    //       gsl_matrix_set(BLOCK, ii, j, value);
    //      }
    //   }
    //    BLOCK_del=get_det(BLOCK);
    //    gsl_matrix_inv(BLOCK);
    //    Ll_init=Ll_init-0.5*BLOCK_del;
    
    //    for (int j=0;j<block[b];j++){
    //     L[b][i][j]=O[b][i][j];
    //    for (int ty=0;ty<type;ty++){
    //      L[b][i][j]=L[b][i][j]-Pest_init[i][ty]*XCB[b][i][ty][j];
    //    }
    //     gsl_vector_set(x_sub,  j, L[b][i][j]);
    //    }
    //    Ll_init=Ll_init-0.5*gsl_VMV(x_sub,BLOCK,x_sub,block[b],block[b]);
    //  }
    // gsl_matrix_free(BLOCK);
    //  gsl_vector_free(x_sub);
    //  }
    
    //   if (Ll_old>Ll){
    //    if (Ll_old>Ll_init){
    //      for (int j=0;j<n;j++){
    //       for (int k=0;k<type;k++){
    //        Pest[j][k]=Pest_old[j][k];
    //       }
    //     }
    //    }else{
    //     for (int j=0;j<n;j++){
    //     for (int k=0;k<type;k++){
    //       Pest[j][k]=Pest_init[j][k];
    //     }
    //     }
    //   }
    //  }else{
    //    if (Ll<Ll_init){
    //    for (int j=0;j<n;j++){
    //      for (int k=0;k<type;k++){
    //       Pest[j][k]=Pest_init[j][k];
    //       }
    //    }
    
    //   }
    
    //   }
    
    printf("%lf\n",Ll);
    printf("%lf\n",Ll_old);
    //   printf("%lf\n",Ll_init);
    
    printf("%d\n",it);
    printf("%lf\n",Msigma);
    printf("%lf\n",Mw);
    printf("%lf\n",Mv);
    
    if (it%10<1){
      for (int j=0;j<nblock;j++){
        for (int i=0;i<type;i++){
          MultiplyMatrix(B[j],Best[j][i],beta_est[j][i],block[j],Q,nbasis[j]);
        }}
      
      
      
      //FILE * pFile2;
      //const char *pFileName2="L.txt";
      //pFile2 = fopen(pFileName2, "w");//è¿™ä¸ªç”¨â€œwâ€æ˜¯å†™æ–‡ä»¶ï¼Œè¦†ç›–åŽŸå†…å®¹ï¼Œè‹¥ä¸æƒ³è¦†ç›–åˆ™ç”¨â€œaâ€
      //if (NULL == pFile2)
      //{//æ–‡ä»¶æ‰“å¼€é”™è¯¯
      //  printf("error");
      //  return 0;
      //}
      //for (int l = 0; l < n; l++)
      //{//å¾ªçŽ¯è¾“å‡ºæ•°æ®å¹¶å†™å…¥
      // printf("%.3f\n", 3.1415926);
      //  for (int i = 0; i < G; i++){
      //   fprintf(pFile2,"%.3f\n", L[l][i]);//è¿™é‡Œå¾ªçŽ¯å†™å…¥æ–‡ä»¶ 3ä¸ª  3.14
      //  }}
      //fclose(pFile2);//æœ€åŽä¸€å®šè®°å¾—å…³é—­æ–‡ä»¶
      
      
      
      //       FILE * pFile4;
      //      const char *pFileName4="Esigma.txt";
      //     pFile4 = fopen(pFileName4, "w");
      //      if (NULL == pFile4)
      //     {
      //      printf("error");
      //      return 0;
      //    }
      //    for (int j = 0; j < nblock; j++)
      //    {
      //      for (int l = 0; l < n; l++)
      //     {
      
      //      for (int i = 0; i < block[j]; i++){
      //        fprintf(pFile4,"%.3f\n", ESigma[j][l][i]);
      //      }}}
      //  fclose(pFile4);
      
      
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
        for (int j = 0; j < nblock; j++)
        {  
          for (int ty = 0; ty < block[j]; ty++){
            
            fprintf(pFile01,"%.3f\n", beta_est[j][0][ty][i]);
          }}}
      fclose(pFile01);
      
      
      
      FILE * pFile02;
      const char *pFileName02="beta2.txt";
      pFile02 = fopen(pFileName02, "w");//è¿™ä¸ªç”¨â€œwâ€æ˜¯å†™æ–‡ä»¶ï¼Œè¦†ç›–åŽŸå†…å®¹ï¼Œè‹¥ä¸æƒ³è¦†ç›–åˆ™ç”¨â€œaâ€
      if (NULL == pFile02)
      {//æ–‡ä»¶æ‰“å¼€é”™è¯¯
        printf("error");
        return 0;
      }
      for (int i = 0; i < Q; i++)
      {
        for (int j = 0; j < nblock; j++)
        {  
          for (int ty = 0; ty < block[j]; ty++){
            
            fprintf(pFile01,"%.3f\n", beta_est[j][1][ty][i]);
          }}}
      fclose(pFile02);//æœ€åŽä¸€å®šè®°å¾—å…³é—­æ–‡ä»¶
      
      
      
      FILE * pFile03;
      const char *pFileName03="beta3.txt";
      pFile03 = fopen(pFileName03, "w");//è¿™ä¸ªç”¨â€œwâ€æ˜¯å†™æ–‡ä»¶ï¼Œè¦†ç›–åŽŸå†…å®¹ï¼Œè‹¥ä¸æƒ³è¦†ç›–åˆ™ç”¨â€œaâ€
      if (NULL == pFile03)
      {//æ–‡ä»¶æ‰“å¼€é”™è¯¯
        printf("error");
        return 0;
      }
      for (int i = 0; i < Q; i++)
      {
        for (int j = 0; j < nblock; j++)
        {  
          for (int ty = 0; ty < block[j]; ty++){
            
            fprintf(pFile01,"%.3f\n", beta_est[j][2][ty][i]);
          }}}
      fclose(pFile03);//æœ€åŽä¸€å®šè®°å¾—å…³é—­æ–‡ä»¶
      
      
      
      FILE * pFile04;
      const char *pFileName04="beta4.txt";
      pFile04 = fopen(pFileName04, "w");//è¿™ä¸ªç”¨â€œwâ€æ˜¯å†™æ–‡ä»¶ï¼Œè¦†ç›–åŽŸå†…å®¹ï¼Œè‹¥ä¸æƒ³è¦†ç›–åˆ™ç”¨â€œaâ€
      if (NULL == pFile04)
      {//æ–‡ä»¶æ‰“å¼€é”™è¯¯
        printf("error");
        return 0;
      }
      for (int i = 0; i < Q; i++)
      {
        for (int j = 0; j < nblock; j++)
        {  
          for (int ty = 0; ty < block[j]; ty++){
            
            fprintf(pFile01,"%.3f\n", beta_est[j][3][ty][i]);
          }}}
      fclose(pFile04);//æœ€åŽä¸€å®šè®°å¾—å…³é—­æ–‡ä»¶
      
      
      
      FILE * pFile05;
      const char *pFileName05="beta5.txt";
      pFile05 = fopen(pFileName05, "w");//è¿™ä¸ªç”¨â€œwâ€æ˜¯å†™æ–‡ä»¶ï¼Œè¦†ç›–åŽŸå†…å®¹ï¼Œè‹¥ä¸æƒ³è¦†ç›–åˆ™ç”¨â€œaâ€
      if (NULL == pFile05)
      {//æ–‡ä»¶æ‰“å¼€é”™è¯¯
        printf("error");
        return 0;
      }
      for (int i = 0; i < Q; i++)
      {
        for (int j = 0; j < nblock; j++)
        {  
          for (int ty = 0; ty < block[j]; ty++){
            
            fprintf(pFile01,"%.3f\n", beta_est[j][4][ty][i]);
          }}}
      fclose(pFile05);//æœ€åŽä¸€å®šè®°å¾—å…³é—­æ–‡ä»¶
      
      FILE * pFile06;
      const char *pFileName06="beta6.txt";
      pFile06 = fopen(pFileName06, "w");//è¿™ä¸ªç”¨â€œwâ€æ˜¯å†™æ–‡ä»¶ï¼Œè¦†ç›–åŽŸå†…å®¹ï¼Œè‹¥ä¸æƒ³è¦†ç›–åˆ™ç”¨â€œaâ€
      if (NULL == pFile06)
      {//æ–‡ä»¶æ‰“å¼€é”™è¯¯
        printf("error");
        return 0;
      }
      for (int i = 0; i < Q; i++)
      {
        for (int j = 0; j < nblock; j++)
        {  
          for (int ty = 0; ty < block[j]; ty++){
            
            fprintf(pFile01,"%.3f\n", beta_est[j][5][ty][i]);
          }}}
      fclose(pFile06);
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
  
  double*** coef_mat=make3Darray(n,Q*type,Q*type);
  
  
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
  
  
  double ***XWX;
  XWX = (double ***)malloc(nblock*sizeof(double **));
  for(int i=0; i<nblock; i++){
    XWX[i] = (double **) malloc(Q*type*nbasis[i]*sizeof(double *));
    for(int j=0; j<(Q*type*nbasis[i]); j++){
      XWX[i][j] = (double *) malloc(Q*type*nbasis[i]*sizeof(double));
    }
  }
  
  //#pragma omp parallel for
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
  
  
  double var_trans;
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
          //printf("%lf\n",var_trans);
          beta_sd[i][ty][j][q]=sqrt(var_trans);
          beta_z[i][ty][j][q]=beta_est[i][ty][j][q]/beta_sd[i][ty][j][q];
        }
      }
    }
    
  }
  
  
  FILE * pFile13;
  const char *pFileName13="beta_z.txt";
  pFile13 = fopen(pFileName13, "w");//è¿™ä¸ªç”¨â€œwâ€æ˜¯å†™æ–‡ä»¶ï¼Œè¦†ç›–åŽŸå†…å®¹ï¼Œè‹¥ä¸æƒ³è¦†ç›–åˆ™ç”¨â€œaâ€
  if (NULL == pFile13)
  {
    printf("error");
    return 0;
  }
  
  for (int ty = 0; ty < type; ty++){
    
    for (int q = 1; q < Q; q++){
      for (int l = 0; l < nblock;l++){
        for (int i = 0; i < block[l]; i++){
          fprintf(pFile13,"%.3f\n", beta_z[l][ty][i][q]);//è¿™é‡Œå¾ªçŽ¯å†™å…¥æ–‡ä»¶ 3ä¸ª 
        }
      }
    }
  }
  fclose(pFile13);
  
  
  FILE * pFile14;
  const char *pFileName14="beta_sd1.txt";
  pFile14 = fopen(pFileName14, "w");//è¿™ä¸ªç”¨â€œwâ€æ˜¯å†™æ–‡ä»¶ï¼Œè¦†ç›–åŽŸå†…å®¹ï¼Œè‹¥ä¸æƒ³è¦†ç›–åˆ™ç”¨â€œaâ€
  if (NULL == pFile14)
  {
    printf("error");
    return 0;
  }
  for (int ty = 0; ty < type; ty++){
    
    for (int q = 1; q < Q; q++){
      for (int l = 0; l < nblock;l++){
        for (int i = 0; i < block[l]; i++){
          fprintf(pFile14,"%.3f\n", beta_sd[l][ty][i][q]);//è¿™é‡Œå¾ªçŽ¯å†™å…¥æ–‡ä»¶ 3ä¸ª 
        }
      }
    }
  }
  fclose(pFile14);
  
  
  
  
  time_t timep;        
  
  time (&timep); 
  
  printf("Time is: %s", asctime( gmtime(&timep) ) );
  
}

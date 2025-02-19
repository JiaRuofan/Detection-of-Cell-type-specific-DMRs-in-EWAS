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

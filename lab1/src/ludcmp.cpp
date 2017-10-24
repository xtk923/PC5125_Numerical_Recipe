#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double **initMat(int n){
  double **a = (double ** )malloc(n * sizeof(double *));
  for(int i = 0; i <n; i++)a[i] = (double *)malloc(n * sizeof(double));
  return a;
};

void printMat(double **a, int n){
  printf("\n");
  for(int i = 0; i < n ; i++){
    for(int j = 0; j < n ; j++){
      printf("  %f  ", a[i][j]);
    }
    printf("\n");
  }
};

void printVec(double *b, int n){
  printf("\n");
  for(int i = 0; i < n; i++){
    printf("  %f  ", b[i]);
  }
  printf("\n");
};

void ludcmp(double **a, int n, int *indx, float *d, int noPivot){

  const double TINY = 1.0e-40;

  int i, imax, j, k;
  double big, temp;
  double *vv;

  vv = (double *)malloc(n * sizeof(double));
  // d stores the information about row interchanges, for future
  // calculation of det(A)
  *d = 1.0;                     /* No row interchanges yet */


  for(i = 0; i<n; i++){        // loop over rows to get the implicit
    big = 0.0;                  // scalling information
    for(j = 0; j < n; j++)
      if((temp=fabs(a[i][j])) > big){ // fabs: returns absolute value
        big = temp;
      }
    if(big == 0.0) {
      printf("Singular matrix in routine ludcmp"); // throw error
    }
                                // No nonezero largest element.
    vv[i] = 1.0/big;            //save the scalling
  }
  for (k = 0; k< n; k++){       // This is the outermost k_ij loop
    big = 0.0;
    for (i = k; i < n; i++){
      temp = vv[i] * fabs(a[i][k]);
      if (temp > big){          // Is the figure of merit for the
        big = temp;             // pivot better than the best so far?
        imax = i;
      }
    }
    if (k != imax && !noPivot){
      for (j = 0; j < n; j++){
        temp = a[imax][j];
        a[imax][j] = a[k][j];
        a[k][j] = temp;
      }
      *d = -(*d);
      vv[imax] = vv[k];
    }
    indx[k] = imax;
    if (a[k][k] == 0.0){
      a[k][k] = TINY;
      // If the pivot element is zero, the matrix is singular
      // (at least to the precision of the algorithm).
      // For some applications on singular matricies, it is
      // desirable to substitute TINY to zero.
    }
    for (i = k +1; i < n; i++){
      temp = a[i][k] /= a[k][k]; // Divide by the pivot element
      for (j = k + 1; j < n; j++){ // Innermost loop: reduce remaining
                                // submatirx
        a[i][j] -= temp * a[k][j];
      }
    }
  }

  free(vv);

}

void lubksb(double **lu, double *b, double *x, int n, int *indx, int noPivot){

  int i, ii=0, ip, j;

  double sum;

  if(noPivot){
    for(int k = 0; k < n; k++){
      indx[k] = k;
    }
  }

  // make sure that the dimensions match

  for (i=0; i<n; i++) x[i] = b[i];
  for (i=0; i<n; i++){          // When ii is set to a positive value
      ip = indx[i];               // it will become the index of the first
      sum = x[ip];                // nonvanishing element of b. We now do
      x[ip] = x[i];               // the forward substitution.

    // The only new wrinkle is to unscramble the permutation as we go.
    if (ii != 0){
      for (j = ii-1; j<i; j++) sum -= lu[i][j]*x[j];
    }else if(sum != 0.0){       // A nonzero element was encountered,
      ii = i+1;                 // so from now on we will have to do the
    }                           // sums in the loop above
    x[i] = sum;
  }
  for (i = n-1; i>=0; i--){     // Now we do the bksb
    sum=x[i];
    for(j = i+1; j<n; j++) sum -= lu[i][j] * x[j];
    x[i] = sum/lu[i][i];        // Store a component of the solution vector X
  }

}

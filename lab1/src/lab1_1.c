#include <stdio.h>
#include <stdlib.h>
#include "ludcmp.cpp"


int main(){

  int dim = 4;

  int *indx = (int *)malloc(dim * sizeof(int));
  float *d = malloc(sizeof(float));
  /* d needs to be a pointer to be refered to in det(A) */

  double **a = initMat(dim);

  a[0][0] = 1.0; a[0][1] = 3.0; a[0][2] = 3.0; a[0][3] = -5.0;
  a[1][0] = 2.0; a[1][1] = -4.0; a[1][2] = 7.0; a[1][3] = -1.0;
  a[2][0] = 7.0; a[2][1] = 1.0/2.0; a[2][2] = 3.0; a[2][3] = -6.0;
  a[3][0] = 9.0; a[3][1] = -2.0; a[3][2] = 3.0; a[3][3] = 8.0;

  printf("Matrix A recorded as:\n");
  printMat(a, dim);

  double *x = (double *)malloc(dim * sizeof(double));
  double *b = (double *)malloc(dim * sizeof(double));
  b[0] = 0.0; b[1] = 2.0; b[2] = 3.0; b[3] = -10.0;

  ludcmp(a, dim, indx, d, 1);

  printf("Matrix A after LU decomposition, stored in one single matrix:\n");
  printMat(a, dim);



  printf("For linear system Ax = b\n");
  printf("b = \n");
  printVec(b, dim);

  lubksb(a, b, x, dim, indx, 1);

  printf("Solving thought backward substitution, x = \n");
  printVec(x, dim);
}

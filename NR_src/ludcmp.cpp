#include <stdio.h>
#include <math.h>
#include <stdlib.h>

struct ludcmp{
  int n;                        // dimension of the matrix
  double **lu;                   // pointer to the matrix
  int *indx;                    // index of permutations
  float *d;                     // for det
  ludcmp(double **a);           // Constructor. Argument is the matrix A
  void solve(double *b, double *x); // Solve a single right-hand side
  void solve(double **b, double **x); // Solve for multiple right-hand sides
  void inverse(double **ainv);        // Calculate matrix inverse A^-1
  double det();                       // return det(A)

};



ludcmp::ludcmp(double **a, int noPivot){

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

#include <stdio.h>
#include <math.h>
#include <time.h>
#include "ludcmp.cpp"



int main(){

  time_t start, end;


  // Solve by setting the voltage of each node as a unknown. We have
  // L^2 unknowns and L^2 equations thus a unique solution can be
  // obtained. Notice that we already have V_00 = 0 and V_LL = 1
  int L;
  printf("Insert the dimension of the grid:\n");
  scanf("%d", &L);

  start = clock();              /* Start the timer */
  int dim = (L+1)*(L+1);        /* we have (L+1)^2 nodes in total */

  double **a = initMat(dim);

  for(int i = 0; i < dim; i++){
    for(int j = 0; j < dim; j++){
      a[i][j] = 0;              /* Since most values will be zero */
    }
  }

  double *x = (double *)malloc(dim * sizeof(double));
  double *b = (double *)malloc(dim * sizeof(double));
  for (int i = 0; i < dim; i++){b[i] = 0;} /* Summation of current is
                                              zero at each node */
  b[dim-1] = 1;                    /* This makes sure that V_LL = 1 */

  for(int i = 0; i<= L; i++){
    for(int j = 0; j<=L; j++){
      /* Assign equations for each node on the grid: */
      /* The row number (L * i + j) denote [i][j], the node we are
         concerning in this equation. The column number denotes the
         source of current from that node to node [i][j]*/
      if((i == 0 && j == 0) || (i == L && j == L)){
        /* This corresponds to point A or B where the voltage is 0 or
           1 hence do not need to be determined */
        a[(L+1) * i + j ][(L+1) * i + j] = 1; /* This will give us V_00 = 0
                                         and V_LL = 1*/

      }else if(i == 0 && j == L){ /* deal with the right bottom
                                       corner */
        /* Take L = 4 as an example */
        /* sqrt(L) [(V_03 - V_04) + (V_14 - V_04)] = 0 */
        a[(L+1) * i + j][(L+1) * i + j - 1]   = sqrt(L);
        a[(L+1) * i + j][(L+1) * (i + 1) + j] = sqrt(L);
        a[(L+1) * i + j][(L+1) * i + j]       = -2 * sqrt(L);

      }else if(i == L && j == 0){ /* Left top corner */
        a[(L+1) * i + j][(L+1) * (i - 1) + j] = sqrt(L);
        a[(L+1) * i + j][(L+1) * i + j +1]    = sqrt(L);
        a[(L+1) * i + j][(L+1) * i + j]       = -2 * sqrt(L);
      }else if(i == 0){
        a[(L+1) * i + j][(L+1) * (i + 1) + j] = sqrt(L);
        a[(L+1) * i + j][(L+1) * i + j + 1]   = sqrt(L);
        a[(L+1) * i + j][(L+1) * i + j - 1]   = sqrt(L);
        a[(L+1) * i + j][(L+1) * i + j]       = -3 * sqrt(L);
      }else if(i == L){
        a[(L+1) * i + j][(L+1) * (i - 1) + j] = sqrt(L);
        a[(L+1) * i + j][(L+1) * i + j + 1]   = sqrt(L);
        a[(L+1) * i + j][(L+1) * i + j - 1]   = sqrt(L);
        a[(L+1) * i + j][(L+1) * i + j]       = -3 * sqrt(L);
      }else if(j == 0){
        a[(L+1) * i + j][(L+1) * (i - 1) + j] = sqrt(L);
        a[(L+1) * i + j][(L+1) * (i + 1) + j] = sqrt(L);
        a[(L+1) * i + j][(L+1) * i + j + 1]   = sqrt(L);
        a[(L+1) * i + j][(L+1) * i + j]       = -3 * sqrt(L);
      }else if(j == L){
        a[(L+1) * i + j][(L+1) * (i - 1) + j] = sqrt(L);
        a[(L+1) * i + j][(L+1) * (i + 1) + j] = sqrt(L);
        a[(L+1) * i + j][(L+1) * i + j - 1]   = sqrt(L);
        a[(L+1) * i + j][(L+1) * i + j]       = -3 * sqrt(L);
      }else{
        a[(L+1) * i + j][(L+1) * (i - 1) + j] = sqrt(L);
        a[(L+1) * i + j][(L+1) * (i + 1) + j] = sqrt(L);
        a[(L+1) * i + j][(L+1) * i + j - 1]   = sqrt(L);
        a[(L+1) * i + j][(L+1) * i + j + 1]   = sqrt(L);
        a[(L+1) * i + j][(L+1) * i + j]       = -4 * sqrt(L);
      }
    }
  }

  int *indx = (int *)malloc(dim * sizeof(int));
  float *d = malloc(sizeof(float));

  ludcmp(a, dim, indx, d, 1);
  lubksb(a, b, x, dim, indx, 1);

  /* We need the two voltages next to point A to solve for the total
     current */

  double V_01, V_10;

  V_01 = x[1];
  V_10 = x[(L+1) * 1];

  double totalCurrent = sqrt(L)*(V_01 + V_10);
  double totalR = 1/totalCurrent;


  printf("Total current is %f\n", totalCurrent);
  printf("Total resistance is %f\n", totalR);

  end = clock();                /* End the timer */

  int t = (end-start) * 1000/CLOCKS_PER_SEC; /* Convert to ms */
  printf("Total running time to determine size of L = %d is %d ms\n", L, t);
}

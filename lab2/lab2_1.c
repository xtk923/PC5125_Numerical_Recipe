#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lab2.h"


int main(){

  // The function polint() is written in the fashion that the index
  // starts from 0
  double xa[4] = {-1.0, 1.0, 2.0, 4.0};
  double ya[4] = {1.25, 2.0, 3.0, 0.0};

  int n = 4;

  double *y = (double *)malloc(sizeof(double));
  double *dy = (double *)malloc(sizeof(double));

  FILE *fp;

  fp = fopen("Q1.dat", "w+");
  for(double x = -1.0; x<=4.0; x+=0.01){

    polint(xa, ya, n, x, y, dy);
    fprintf(fp, "%f  %f\n", x, *y);

  }
  fclose(fp);

  return(0);


}

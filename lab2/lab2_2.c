#include <stdio.h>
#include "lab2.h"


double _function(double x){
  return (pow(x, 4) * log(x + sqrt(pow(x, 2) +1)));
};

int main(){

  double integral = qromb(_function, 0.0, 2.0);
  printf("K = %d\n", K);
  printf("Using double precision.\n");
  printf("Integral result is \n%.30f\n", integral);
  printf("Result from Mathematica: \n8.1533641198111650205\n");
  printf("Difference: \n%.30f\n", integral - 8.1533641198111650205);

}

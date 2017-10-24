#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define EPS 1.0e-20
#define JMAX 10
#define JMAXP (JMAX+1)
#define K 1

#define FUNC(x) ((*func)(x))


void polint(double xa[], double ya[], int n, double x, double *y, double *dy){

  int i, m, ns = 0;
  double den, dif, dift, ho, hp, w;
  double *c, *d;

  dif = fabs(x - xa[0]);
  c = (double *)malloc(n * sizeof(double));
  d = (double *)malloc(n * sizeof(double));

  for(i = 0; i < n; i++){       /* Here we find the index ns of the
                                   closest table entry */
    if((dift = fabs(x-xa[i])) < dif){
      ns = i;
      dif = dift;
    }
    c[i] = ya[i];               /* and initialize the tableau of c's
                                   and d's */
    d[i] = ya[i];
    /* This should return c and d as xa and ya, if n = dim(xa) */
  }
  *y = ya[ns--];                /* This is the initial approximation
                                   to y */
  for(m =1; m<n; m++){
    for(i = 0; i< n-m; i++){
      ho = xa[i]-x;
      hp = xa[i+m] -x;
      w = c[i+1] - d[i];
      if((den = ho-hp) == 0.0) printf("Error in routine polint\n");
      /* This error can occur only if two input xa'a are(to within
         roundoff) identical */
      den = w/den;
      d[i] = hp*den;            /* Here the c's and d's are updated */
      c[i] = ho*den;
    }
    *y += (*dy = (2*(ns + 1 ) < (n-m) ? c[ns + 1] : d[ns--]));
    /* After each column in the tableau is completed, we decide which
       correction, c or d, we want to add to our accumulating value of
       y. i.e. which path to take through the tableau--forking up or
       down. We do this in such a way as to take the most "straight
       line" route through the tableau to its apex, updating ns
       accordingly to keep track of where we are. This route keeps the
       partial approximations centered(insofar as possible) on the
       target x. The last dy added is thus the error indication. */
  }
  free(d);
  free(c);


};

double trapzd(double (*func)(double), double a, double b, int n){
  double x, tnm, sum, del;
  static double s;
  int it, j;

  if(n == 1){
    return (s = 0.5*(b-a)*(FUNC(a)+FUNC(b)));
  }else{
    for (it = 1, j = 1; j < n -1; j++){
      it <<= 1;
    }
    tnm = it;
    del = (b-a)/tnm;
    x = a + 0.5*del;
    for(sum = 0.0, j = 0; j<it; j++, x+=del){
      sum += FUNC(x);
    }
    s = 0.5*(s+(b-a)*sum/tnm);
    return s;

  }

};




double qromb(double (*func)(double), double a, double b){


  void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
  double trapzd(double (*func)(double), double a, double b, int n);
  double ss, dss;
  double s[JMAX], h[JMAXP];
  int j;

  h[0] = 1.0;
  for (j = 1; j<=JMAX; j++){
    s[j-1] = trapzd(func, a, b, j-1);
    if (j >= K){
      polint(&h[j-K], &s[j-K], K, 0.0, &ss, &dss);
      if(fabs(dss) <= EPS*fabs(ss)){
        return ss;
      }
    }
    h[j] = 0.25*h[j-1];
  }
  printf("Too many steps in routine qromb\n");
  return 0.0;
};

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

#define EPS 1.0e-10
#define STEP 1.0e-2

void randomR(double *r){

  double dx = ((double)rand()/RAND_MAX - 0.5)*2;
  double dy = ((double)rand()/RAND_MAX - 0.5)*2;
  double dz = ((double)rand()/RAND_MAX - 0.5)*2;



  // in the rare case where r = (0, 0, 0) will be obtained, redo
  if(r[0] == -dx && r[1] == -dy && r[2] == -dz){
    randomR(r);
    printf("This actually happend!\n");
  }

  r[0] += dx;
  r[1] += dy;
  r[2] += dz;


};

void updateR(double *r){
  r[3] = r[0];
  r[4] = r[1];
  r[5] = r[2];

};

void resetR(double *r){
  r[0] = r[3];
  r[1] = r[4];
  r[2] = r[5];
}

double norm(double *u, double *v, int n){

  double result = 0.0;
  for(int i = 0; i< n; i++){
    result += pow((u[i]-v[i]), 2);
  }
  return sqrt(result);
};


double psi(double *r1, double *r2, double *rA, double *rB){

  double psiA1 = exp(-norm(rA, r1, 3));
  double psiB2 = exp(-norm(rB, r2, 3));
  double psiA2 = exp(-norm(rA, r2, 3));
  double psiB1 = exp(-norm(rB, r1, 3));

  double psi = psiA1*psiB2 + psiA2*psiB1;

  return psi;

};

double php(double *r1, double *r2, double *rA, double *rB){
  double PHP = 0.0;
  double P = psi(r1, r2, rA, rB);
  double K, V;                  /* Kinetic and potential */

  // result copied from mathematica
  K = exp(-norm(r2, rA, 3) - norm(r1, rB, 3))*(-1 + 1/(norm(r2, rA, 3)) + 1/norm(r1, rB, 3))+
    exp(-norm(r1, rA, 3) - norm(r2, rB, 3))*(-1 + 1/norm(r1, rA, 3) + 1/norm(r2, rB, 3)    );



  V = (-1/norm(rA, r1, 3) - 1/norm(rB, r1, 3) - 1/norm(rA, r2, 3) - 1/norm(rB, r2, 3)
       +1/norm(r1, r2, 3) + 1/norm(rA, rB, 3))*P;


  double HP = K + V;

  PHP = P*HP;
  return PHP;
};


int main(){

  srand(time(NULL));            /* seed the RNG */
  // r1 and r2 will store the current position and there previous positions
  double r1[6];
  double r2[6];


  double rA[3] = {0.0, 0.0, 0.0}; /* fix A at origin */
  double rB[3] = {0.0, 0.0, 0.0}; /* w/o loss of generality fix B
                                         on the x-axis */

  double p1, p2;                /* previous P(r1, r2) and current P(r1, r2) */

  double PHP;                /* \psi H \psi */
  double PP;                 /* \psi \psi */



  double Energy1,  Energy2;

  FILE *OUTPUT_FILE;
  OUTPUT_FILE = fopen("Q3.dat", "w+");
  int counter;


  for(int i = 5; i*STEP < 10; ++i){
    PHP = 0.0;
    PP = 0.0;
    r1[0] = i*STEP/2; r1[1] = 0.0; r1[2] = 0.0; /* initiate the electrons between the two nuclei */
    r2[0] = i*STEP/2; r2[1] = 0.0; r2[2] = 0.0;
    randomR(r1); randomR(r2);     /* generate random r1, r2 */
    updateR(r1); updateR(r2);     /* update data */
    rB[0] = i*STEP;
    counter = 0;
    Energy1 = -0.1;             /* reset energies */
    Energy2 = 0.0;

    p1 = pow(psi(r1, r2, rA, rB), 2);    /* initiate p1 */


    FILE *RECORD_FILE;
    if(i*STEP == 0.5){
      RECORD_FILE = fopen("0.5.dat", "w+");
    }else if(i*STEP == 1.0){
      RECORD_FILE = fopen("1.0.dat", "w+");
    }else if(i*STEP == 8.0){
      RECORD_FILE = fopen("8.0.dat", "w+");
    }

    do{


      randomR(r1);                /* generate random r1 */
      randomR(r2);                /* and r2 */
      p2 = pow(psi(r1, r2, rA, rB), 2); /* calculate new proba */

      if(p2 > p1 || p2/p1 > (double)rand()/RAND_MAX){
        Energy1 = Energy2;
        /* if this case if more likely, or it is taken with a pabability p2/p1 */
        updateR(r1); updateR(r2); /* update the positions */
        p1 = p2;                  /* and probability */
        PHP += php(r1, r2, rA, rB);
        PP += p2;                 /* add the value to PP */
        counter++;
        Energy2 = PHP/PP;

        if(i*STEP ==0.5||i*STEP == 1.0 || i*STEP == 8.0){
          fprintf(RECORD_FILE, "%f %f %f\n", r1[0], r1[1], r1[2]);
          fprintf(RECORD_FILE, "%f %f %f\n", r2[0], r2[1], r2[2]);

        }


      }else{
        resetR(r1); resetR(r2);
      }




    }while(fabs(Energy1-Energy2) >= EPS);
    if(RECORD_FILE) fclose(RECORD_FILE);
    printf("%d steps for i= %d\n", counter, i);

    fprintf(OUTPUT_FILE, "%f  %f\n", i*STEP, Energy2);
  }
  fclose(OUTPUT_FILE);




}

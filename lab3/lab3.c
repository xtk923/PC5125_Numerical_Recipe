#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define EPS 0.1

void create_lattice(int ***lattice, int N){
  /*
    The lattice is created by grouping two C atoms in the same unit
    cell. By doing so, we have a 2D array of unit cells, with each
    unit cell equivalent to the others. Consequently, when we want to
    join two atoms from two boundaries, it will always to atom 1
    connected to atom 2.
   */
  for(int i = 0; i< N; i++){
    lattice[i] = malloc(N * sizeof(int*));
    for(int j = 0; j<N; j++){
      lattice[i][j] = malloc(2 * sizeof(int));
      for(int k = 0; k<2; k++){
        lattice[i][j][k] = 1;
        /* This is better for the lattice at low temperature because
           it starts at equilibrium state. For high temperature, the
           kinetic energy is sufficient to alter the state to
           equilibrium */
      }
    }
  }
};

void flip_lattice(int ***lattice, int i, int j, int k){
  lattice[i][j][k] *= -1;
};

float calculate_energy(int ***lattice, int N){
  float total_energy = 0;
  for(int i = 0; i< N; i++){
    for(int j = 0; j< N; j++){
      for(int k = 0; k < 2; k++){
        if(k == 0){
          /* for the atom on the left hand side of the unit cell */
          /* its nearest neighboors are: */
          /* the other atom in the same unit cell */
          total_energy -= lattice[i][j][0] * lattice[i][j][1];
          /* the right hand side atom in the unit cell on its left */
          total_energy -= lattice[i][j][0] * lattice[i][(j-1 + N) % N][1];
          /* the right hand side atom in the unit cell below it */
          total_energy -= lattice[i][j][0] * lattice[(i-1 + N)%N][j][1];

          /* (i-1 +N)%N is used so that it joins boundaries */
        }
        else if(k == 1){
          total_energy -= lattice[i][j][1] * lattice[i][j][0];

          total_energy -= lattice[i][j][1] * lattice[i][(j+1 + N) % N][0];

          total_energy -= lattice[i][j][1] * lattice[(i+1 + N)%N][j][0];

        }else{
          printf("error!\n");
        }
      }
    }
  }
  total_energy /= 2;            /* each pair has been counted twice */
  return total_energy;
};

/* If only one site is flipped, we don't have to iterate over the
   entire lattice */
float calculate_difference(int ***lattice, int i, int j, int k, int N){
  float energy_diff = 0;
  if(k == 0){
    energy_diff += 2* lattice[i][j][0] * lattice[i][j][1];
    energy_diff += 2* lattice[i][j][0] * lattice[i][(j-1 + N) % N][1];
    energy_diff += 2* lattice[i][j][0] * lattice[(i-1 + N)%N][j][1];

  }else if(k == 1){
    energy_diff += 2* lattice[i][j][1] * lattice[i][j][0];

    energy_diff += 2* lattice[i][j][1] * lattice[i][(j+1 + N) % N][0];

    energy_diff += 2* lattice[i][j][1] * lattice[(i+1 + N)%N][j][0];

  }else{
    printf("error! k = %d\n", k);
  }

  return energy_diff;
};

void random_site(int *i, int *j, int *k, int N){
  *i = rand() % N;
  *j = rand() % N;
  *k = rand() % 2;
};




void print_lattice(int ***lattice, int N){
  for(int i =0; i< N; i++){
    printf("\n");
    for(int j = 0; j< N; j++){
      printf(" [ ");
      for(int k = 0; k< 2; k++){
        printf(" %d ", lattice[i][j][k]);
      }
      printf("] ");
    }
    printf("\n");
  }
};


float heat_capacity(float *energies, int instances, int N, float T){
  float c = 0;
  c = 1/(T*T*N*N*2);            /* remark: total number of spin is N*N*2 */
  float E_square_mean = 0, E_mean_square = 0;
  for(int i = 0; i< instances; i++){
    E_square_mean += pow(energies[i], 2);
    E_mean_square += energies[i];
  }
  E_square_mean /= instances;
  E_mean_square /= instances;
  E_mean_square = pow(E_mean_square, 2);
  c *= (E_square_mean - E_mean_square);
  return c;

};

/* the average spin of the lattice */
float average_spin(int ***lattice, int N){
  float spin = 0.0;
  for(int i = 0; i < N; i++){
    for(int j = 0; j < N; j++){
      for(int k = 0; k < 2; k++){
        spin += lattice[i][j][k];
      }
    }
  }
  spin /= (pow(N, 2) * 2);
  spin = fabs(spin);

  return spin;
};

/* the mean average spin of all iterations */
float mean_average_spin(float *spins, int instances){
  float result = 0;
  for(int i = 0; i < instances; i++){
    result += spins[i];
  }

  result /= instances;
  return result;
};

float _mean(float *sample, int instances){
  float sum = 0;
  for(int i = 0; i< instances; i++){
    sum += sample[i];
  }
  float mean = sum/(float)instances;
  return mean;
};


/* this might exist in stdlib.h, but I just redefine it here */
float corrected_sample_error(float *sample, int instances){
  float mean = _mean(sample, instances);

  float result = 0;
  for(int i = 0; i < instances; i++){
    result += pow((sample[i] - mean), 2);
  }
  return sqrt(result/(instances*(instances - 1)));
};





int main(){


  printf("What's the dimension you want to use? (total 2*N^2 sits)\n");
  char strN[10];
  scanf("%s", strN);

  printf("How many samples would you like to take?\n");
  char strSample_size[10];
  scanf("%s", strSample_size);

  printf("How many iterations would you like to run?\n");
  char strIterations[10];
  scanf("%s", strIterations);

  FILE *c_FILE;
  c_FILE = fopen("heat_capacity.dat", "w+");

  FILE *s_FILE;
  s_FILE = fopen("spin.dat", "w+");
  float total_energy;

  int N = atoi(strN);                   /* this gives us N*N*2 sites */


  int ***lattice = malloc(N * sizeof(int**));
  int i, j, k;

  int sample_size = atoi(strSample_size);
  int total_iteration = atoi(strIterations);  /* total number of sweeps */

  for(float T = 0.1; T < 5; T+=EPS){

    float *spin_record = malloc(sample_size * sizeof(float));
    /* This is to record average spin so that we can take its standard error */
    float *heat_capacity_record = malloc(sample_size * sizeof(float));
    /* The same for heat capacity */

    for(int r = 0; r < sample_size; r++){
      /* Take 30 samples such that it has statistical significance */
      create_lattice(lattice, N);


      float *energies = malloc(total_iteration * sizeof(float));
      /* an array of total energies */
      float *spins = malloc(total_iteration * sizeof(float));
      /* an array of average spings */
      int instances = 0;      /* a counter for the above two arrays */
      for(int iteration = 0; iteration < total_iteration; iteration++){
        /* Sweep over the entire lattice */
        float prob;               /* the probability of such flip */
        for(int counter = 0; counter < N*N*2; counter++){
          random_site(&i, &j, &k, N);    /* choose a random site */
          /* decide whether to accept the flip */
          /* P(configuration1)/P(configuration2) = e^{-(H2 - H1)} */
          /* k_b = 1 */
          float E_diff = calculate_difference(lattice, i, j, k, N);
          prob = exp(-E_diff/T);

          if(prob > rand()/((float)RAND_MAX+1)){

            /* if the ratio > 1, it must be bigger than the random
               value, hence we accept */
            /* else, we check it against the random float in [0, 1] */
            flip_lattice(lattice, i, j, k);
          }
        }
        /* calculate the total energy */
        total_energy = calculate_energy(lattice, N);
        if(iteration > total_iteration/10){
          energies[instances] = total_energy;
          spins[instances] = average_spin(lattice, N);
        instances++;
        }

      }
      heat_capacity_record[r] = heat_capacity(energies, instances, N, T);
      spin_record[r] = mean_average_spin(spins, instances);

    }
    fprintf(c_FILE, "%f    %f    %f\n",
            T, _mean(heat_capacity_record, sample_size),
            corrected_sample_error(heat_capacity_record, sample_size));
    fprintf(s_FILE, "%f    %f    %f\n",
            T, _mean(spin_record, sample_size),
            corrected_sample_error(spin_record, sample_size));
  }

}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define RAND01 ((double)random() / (double)RAND_MAX)

void random_fill_LR(int nU, int nI, int nF, double *L, double *R) {
  srandom(0);
  for(int i = 0; i < nU; i++)
    for(int j = 0; j < nF; j++)
      L[i*nF + j] = RAND01 / (double) nF;
  for(int i = 0; i < nF; i++)
    for(int j = 0; j < nI; j++)
      R[j*nF + i] = RAND01 / (double) nF;
}

void estimate_B(int num_non_zeros, int num_feats, int num_columns, double *B,
                double *L, double *R, int *A_r, int *A_c) {
  for (int i = 0; i < num_non_zeros; i++) {
    double sum = 0;
    for (int j = 0; j < num_feats; j++)
      sum += L[A_r[i] * num_feats + j] * R[A_c[i] * num_feats + j];
    B[A_r[i] * num_columns + A_c[i]] = sum;
  }
}

void calculate_B(int num_rows, int num_columns, int num_feats, double *B, double *L, double *R){
  for(int i = 0; i < num_rows; i++){
    for(int j = 0; j < num_columns; j++){
      double sum = 0;	
      for(int k = 0; k < num_feats; k++){
        sum += L[i*num_feats + k]*R[j*num_feats + k];
      }
      B[i*num_columns + j] = sum;
    }	
  }
}

void calculate_L_and_R(int num_rows, int num_columns, int num_feats, int num_non_zeros,
                       double alpha, double *A_val, int *A_r, int *A_c,
                       double *B, double *L, double *R){

  double *L_copy = (double *)malloc(sizeof(double) * num_rows * num_feats);
  memcpy(L_copy, L, sizeof(double) * num_rows * num_feats);
  double *R_copy = (double *)malloc(sizeof(double) * num_feats * num_columns);
  memcpy(R_copy, R, sizeof(double) * num_feats * num_columns);

  for (int n = 0; n < num_non_zeros; n++) {
    double delta = A_val[n] - B[A_r[n]*num_columns + A_c[n]];
    for (int f = 0; f < num_feats; f++) {
      L[A_r[n]*num_feats + f] += alpha*2*delta*R_copy[A_c[n]*num_feats + f];
      R[A_c[n]*num_feats + f] += alpha*2*delta*L_copy[A_r[n]*num_feats + f];
    }
  }

  free(L_copy);
  free(R_copy);
}

int main(int argc, char **argv){
  if (argc != 2) {
    printf("Wrong usage: matFact [input file]\n");
    return 1;
  }
  
  FILE *f = fopen(argv[1], "r");
  int num_iters, num_feats, num_rows, num_columns, num_non_zeros;
  double alpha; 

  fscanf (f, "%d", &num_iters);
  fscanf (f, "%lf", &alpha);
  fscanf (f, "%d", &num_feats);
  fscanf (f, "%d %d %d", &num_rows, &num_columns, &num_non_zeros);
    
  double* L = (double *)malloc(sizeof(double) * num_rows * num_feats);
    
  double *R = (double *)malloc(sizeof(double) * num_feats * num_columns);

  double *B = (double *)malloc(sizeof(double) * num_rows * num_columns);
  
  double *A_val = (double *)malloc(sizeof(double) * num_non_zeros);
  int *A_r = (int *)malloc(sizeof(int) * num_non_zeros);
  int *A_c = (int *)malloc(sizeof(int) * num_non_zeros);
    
  for (int i = 0; i < num_non_zeros; i++) {
    int row, col;
    double value;
    fscanf(f, "%d %d %lf", &row, &col, &value);
    if (value < 1 || value > 5){
      printf("Rating is not between 1 and 5\n");
      return 1;
    }
    A_r[i] = row;
    A_c[i] = col;
    A_val[i] = value;
  }
  fclose(f);

  random_fill_LR(num_rows, num_columns, num_feats, L, R);
  estimate_B(num_non_zeros, num_feats, num_columns, B, L, R, A_r, A_c);

  for (int i = 0; i < num_iters; i++) {
    calculate_L_and_R(num_rows, num_columns, num_feats, num_non_zeros,
                      alpha, A_val, A_r, A_c, B, L, R);
    estimate_B(num_non_zeros, num_feats, num_columns, B, L, R, A_r, A_c);
  }

  calculate_B(num_rows, num_columns, num_feats, B, L, R);
  
  int nz = 0;
  for(int i = 0; i < num_rows; i++){
    double max = 0;
    int index;
    for(int j = 0; j < num_columns; j++) {
      if(A_r[nz] != i || A_c[nz] != j) {
        if(B[i*num_columns + j] > max) {
          max = B[i*num_columns + j];
          index = j;
        }
      } else
        nz++;
    }
    printf("%d\n", index);
  }

  free(L);
  free(R);
  free(B);
  free(A_val);
  free(A_r);
  free(A_c);

  return 0;
}	

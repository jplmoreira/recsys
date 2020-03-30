#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

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
    
int is_non_zero(int row, int col, int num_non_zeros, int *A_r, int *A_c) {
  for (int i = 0; i < num_non_zeros; i++) {
    if (row == A_r[i] && col == A_c[i])
      return 1;
  }
  return 0;
}

int main(int argc, char **argv){
  if (argc != 2) {
    printf("Wrong usage: matFact [input file]\n");
    return 1;
  }
  
  time_t start_t, end_t;
  double diff_t;
  
  time(&start_t);

  
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

  time(&end_t);
  printf("----> Inicializacao das matrizes ----> %f\n", difftime(end_t, start_t));
  time(&start_t);

  random_fill_LR(num_rows, num_columns, num_feats, L, R);
  calculate_B(num_rows, num_columns, num_feats, B, L, R);
  
  time(&end_t);
  printf("----> Random fill dos Ls r Rs e primeiro B ----> %f\n", difftime(end_t, start_t));
  time(&start_t);
  
  time_t start_LR;
  time_t start_B;
  time_t end_LR;
  time_t end_B;
  double time_LR = 0; 
  double time_B = 0;
  for (int i = 0; i < num_iters; i++) {
	time(&start_LR);
    calculate_L_and_R(num_rows, num_columns, num_feats, num_non_zeros,
                      alpha, A_val, A_r, A_c, B, L, R);
	time(&end_LR);
	time_LR += difftime(end_LR, start_LR);
	time(&start_B);
    calculate_B(num_rows, num_columns, num_feats, B, L, R);
	time(&end_B);
	time_B += difftime(end_B, start_B);
  }
  
  
  time(&end_t);
  printf("---->Itereçoes essenciais----> %f\n", difftime(end_t, start_t));
  time(&start_t);
  
  printf("        ---->Calculo dos Ls e Rs-----> %f\n", time_LR);
  printf("        ---->Calculo dos B-----> %f\n", time_B);

  FILE *out = fopen("output.txt", "w");
  int nz= 0;
  for(int i = 0; i < num_rows; i++){
    double max = 0;
    int index;
    for(int j = 0; j < num_columns; j++){
      if(A_r[nz] == i && A_c[nz] == j){
	nz++;
      } else {
        if(B[i*num_columns + j] > max){
          max = B[i*num_columns + j];
          index = j;
        }
      }
    }
	fprintf(out, "%d", index);
	fputs("\n", out);
    //printf("%d\n", index);
  }
  time(&end_t);
  printf("----> Output dos resultados ----> %f\n", difftime(end_t, start_t));
  fclose(out);
  free(L);
  free(R);
  free(B);
  free(A_val);
  free(A_r);
  free(A_c);

  return 0;
}	

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define RAND01 ((double)random() / (double)RAND_MAX)

double*** L;
double*** R;
double** B;
double** A;
int cur = 0, prev = 1;

void random_fill_LR(int nU, int nI, int nF)
{
  srandom(0);
  for(int i = 0; i < nU; i++)
    for(int j = 0; j < nF; j++)
      L[cur][i][j] = RAND01 / (double) nF;
  for(int i = 0; i < nF; i++)
    for(int j = 0; j < nI; j++)
      R[cur][i][j] = RAND01 / (double) nF;
  /* memcpy(L[cur][0], (double[]){ 0.420094, 0.197191 }, sizeof(double)*2); */
  /* memcpy(L[cur][1], (double[]){ 0.391550, 0.399220 }, sizeof(double)*2); */
  /* memcpy(L[cur][2], (double[]){ 0.455824, 0.098776 }, sizeof(double)*2); */
  /* memcpy(L[cur][3], (double[]){ 0.167611, 0.384115 }, sizeof(double)*2); */
  /* memcpy(R[cur][0], (double[]){ 0.138887, 0.276985, 0.238699, 0.314435, 0.182392 }, sizeof(double)*5); */
  /* memcpy(R[cur][1], (double[]){ 0.256700, 0.476115, 0.458098, 0.317856, 0.358648 }, sizeof(double)*5); */
}


void calculate_B(int num_rows, int num_colums, int num_feats){
  int i, j, k;
  double sum = 0;	
  for(i = 0; i < num_rows; i++){
    for(j = 0; j < num_colums; j++){
      for(k = 0; k < num_feats; k++){
        sum = sum + L[cur][i][k]*R[cur][k][j];
      }
      B[i][j] = sum;
      sum = 0;
    }	
  }
}

void estimate_A(){
    
}

void calculate_L_and_R(int num_rows, int num_colums, int num_feats, double alpha){
  int l, c;
  int k, i;
  double delta;
  double sum_L = 0, sum_R = 0;
  prev = cur;
  cur = ++cur % 2;
  for(l = 0; l < num_rows; l++){
    for(c = 0; c < num_colums; c++){
      if(A[l][c] != 0){
        delta = A[l][c] - B[l][c];
        for(k = 0; k < num_feats; k++){
          for(i = 0; i < num_colums; i++){
            if (A[l][i] != 0)
              sum_L += 2*delta*(-R[prev][k][i]);
          }
          L[cur][l][k] = L[prev][l][k] - alpha*sum_L;
          for(i = 0; i < num_rows; i++){
            if (A[i][c] != 0)
              sum_R += 2*delta*(-L[prev][i][k]);
          }
          R[cur][k][c] = R[prev][k][c] - alpha*sum_R;
          sum_L = sum_R = 0;
        }
      } 
    }
  }
    
}
    


int main(int argc, char **argv){
  FILE *f = fopen(argv[1], "r");
  int num_iters;
  double alpha;
  int num_feats;
  int num_rows, num_colums, num_non_zeros;
  int row, col;
  double value;
  fscanf (f, "%d", &num_iters);
  fscanf (f, "%lf", &alpha);
  fscanf (f, "%d", &num_feats);
  fscanf (f, "%d %d %d", &num_rows, &num_colums, &num_non_zeros);
    
  
  L = (double ***)malloc(2 * sizeof(double**));
  L[0] = (double **)malloc(num_rows * sizeof(double*));
  L[1] = (double **)malloc(num_rows * sizeof(double*));
  for(int d = 0; d < num_rows; d++) {	
    L[0][d] = (double *)malloc(num_feats * sizeof(double));
    L[1][d] = (double *)malloc(num_feats * sizeof(double));
  }
    
  R = (double ***)malloc(2 * sizeof(double**));
  R[0] = (double **)malloc(num_feats * sizeof(double*));
  R[1] = (double **)malloc(num_feats * sizeof(double*));
  for(int d = 0; d < num_feats; d++) {
    R[0][d] = (double *)malloc(num_colums * sizeof(double));
    R[1][d] = (double *)malloc(num_colums * sizeof(double));
  }
    
  B = (double **)malloc(num_rows * sizeof(double));
  for(int d = 0; d < num_rows; d++){
    B[d] = (double *)malloc(num_colums * sizeof(double));
  }
    
  A = (double **)malloc(num_rows * sizeof(double));
  for(int d = 0; d < num_rows; d++){
    A[d] = (double *)malloc(num_colums * sizeof(double));
  }
    
  //double A[num_rows][num_colums];
  for(int i = 0; i < num_rows; i++){
    for(int j = 0; j < num_colums; j++){
      A[i][j] = 0;
    }
  }
  while(fscanf(f, "%d %d %lf", &row, &col, &value) != EOF){
    if (value < 1 || value > 5){
      printf("Rating is not between 1 and 5");
      return 1;
    }
    A[row][col] = value;
    //printf("----->%f\n", A[row][col]);
  }
  fclose(f);
  random_fill_LR(num_rows, num_colums, num_feats);
  calculate_B(num_rows, num_colums, num_feats);
  int i = 0;
  estimate_A();
  while(i < num_iters){
    calculate_L_and_R(num_rows, num_colums, num_feats, alpha);
    calculate_B(num_rows, num_colums, num_feats);
    i++;
  }
  for(int i = 0; i < num_rows; i++){
    double max = 0;
    int index;
    for(int j = 0; j < num_colums; j++){
      if(A[i][j] == 0){
        if(B[i][j] > max){
          /* printf("----->%f\n", max); */
          max = B[i][j];
          index = j;
        }
      }
    }
    max = 0;
    printf("%d\n", index);
  }
  /*int flot;
	fscanf (f, "%d", &flot);
	fclose(f);
	printf("----> %d", flot);*/
  /* printf("\nFinal Matrix L\n"); */
  /* for (int i = 0; i < num_rows; i++) { */
  /*   for (int j = 0; j < num_feats; j++) */
  /*     printf("%lf ", L[cur][i][j]); */
  /*   printf("\n"); */
  /* } */
  /* printf("\nFinal Matrix R\n"); */
  /* for (int i = 0; i < num_feats; i++) { */
  /*   for (int j = 0; j < num_colums; j++) */
  /*     printf("%lf ", R[cur][i][j]); */
  /*   printf("\n"); */
  /* } */
  /* printf("\nFinal Matrix B\n"); */
  /* for (int l = 0; l < num_rows; l++) { */
  /*   for (int c = 0; c < num_colums; c++) */
  /*     printf("%lf ", B[l][c]); */
  /*   printf("\n"); */
  /* } */
  return 1;
}	

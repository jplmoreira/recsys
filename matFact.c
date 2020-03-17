#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define RAND01 ((double)random() / (double)RAND_MAX)

double **L;
double **R;

double **A;
double **B;

void random_fill_LR(int nU, int nI, int nF) {
  srandom(0);

  for (int i = 0; i < nU; i++)
    for (int j = 0; j < nF; j++)
      L[i][j] = RAND01 / (double) nF;

  for (int i = 0; i < nF; i++)
    for (int j = 0; j < nI; j++)
      R[i][j] = RAND01 / (double) nF;
}

int main(int argc, const char* argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage: matFact [input file]\n");
    exit(-1);
  }

  int iter, feats, rows, cols, nZ;
  double alfa;

  char b1[10], b2[10], b3[10];

  FILE *f;
  f = fopen(argv[1], "r");
  if (f == NULL) {
    fprintf(stderr, "Could not open file");
    exit(-1);
  }

  fscanf(f, "%d\n", &iter);
  fscanf(f, "%lf\n", &alfa);
  fscanf(f, "%d\n", &feats);
  fscanf(f, "%d %d %d\n", &rows, &cols, &nZ);

  A = (double **) malloc (rows * sizeof(double*));
  B = (double **) malloc (rows * sizeof(double*));
  for (int i = 0; i < rows; i++) {
    A[i] = (double *) malloc (cols * sizeof(double));
    B[i] = (double *) malloc (cols * sizeof(double));
    for (int j = 0; j < cols; j++)
      A[i][j] = 0;
  }

  for (int i = 0; i < nZ; i++) {
    int r, c;
    double val;
    fscanf(f, "%d %d %lf\n", &r, &c, &val);
    A[r][c] = val;
  }

  return 0;
}

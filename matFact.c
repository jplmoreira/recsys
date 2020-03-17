#include <stdio.h>
#include <stdlib.h>


#define RAND01 ((double)random() / (double)RAND_MAX)

double** L;
double** R;
double** B;
double** A;

void random_fill_LR(int nU, int nI, int nF)
{
	srandom(0);
	for(int i = 0; i < nU; i++)
	for(int j = 0; j < nF; j++)
	L[i][j] = RAND01 / (double) nF;
	for(int i = 0; i < nF; i++)
	for(int j = 0; j < nI; j++)
	R[i][j] = RAND01 / (double) nF;
}


void calculate_B(int num_rows, int num_colums, int num_feats){
	int i, j, k;
	double sum = 0;	
	for(i = 0; i < num_rows; i++){
		for(j = 0; j < num_colums; j++){
			for(k = 0; k < num_feats; k++){
				sum = sum + L[i][k]*R[k][j];
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
	double prev;
	double delta;
	double sum = 0;
	for(l = 0; l < num_rows; l++){
		for(c = 0; c < num_colums; c++){
			if(A[l][c] != 0){
				for(k = 0; k < num_feats; k++){
					prev = L[l][k];
					for(i = 0; i < num_colums; i++){
						sum += 2*(A[l][i] - B[l][i])*(-R[k][i]);
					}
					L[l][k] = prev - alpha*sum;
					sum = 0;
				}
				for(k = 0; k < num_feats; k++){
					prev = R[k][c];
					for(i = 0; i < num_rows; i++){
						sum += 2*(A[l][i] - B[l][i])*(-L[i][k]);
					}
					R[k][c] = prev - alpha*sum;
					sum = 0;
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
	
	L = (double **)malloc(num_rows * sizeof(double));
	for(int d = 0; d < num_rows; d++)	
		L[d] = (double *)malloc(num_feats * sizeof(double));
	
	R = (double **)malloc(num_feats * sizeof(double));
	for(int d = 0; d < num_feats; d++)	
		R[d] = (double *)malloc(num_colums * sizeof(double));
	
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
			if(A[i][j] != 0){		
				if(B[i][j] > max){
					//printf("----->%f\n", max);
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
	return 1;
	
}
all: serial parallel

serial: matFact.c
	gcc -o matFact matFact.c

parallel: matFact-omp.c
	gcc -fopenmp -o matFact-omp matFact-omp.c

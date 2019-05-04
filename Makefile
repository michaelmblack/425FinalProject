make:
	mpicc -o project ParallelJacobi.c -g -fopenmp -lm
debug:
	mpicc -o debug Jacobi.c -g
serial:
	gcc -o serial Jacobi.c -fopenmp
gen:
	gcc -o gen Generator.c


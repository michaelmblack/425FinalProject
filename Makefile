make:
	gcc -o project Jacobi.c -lm
debug:
	gcc -o debug Jacobi.c -g -lm

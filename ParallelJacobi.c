#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#define threadCount 5

void readInMatrix(char *fileName, long double *Arr[], long double Ans[]);
long double RelDif(long double a, long double b);
void sendRecieveDebug(long double **Arr, int size, int rank, MPI_Comm comm, int np);
int checkTolerance(long double tolerance, long double *A, long double *B, int size);
void ParallelJacobian(long double *A[], long double Ans[], int n, long double tolerance, int np, int rank, MPI_Comm comm);
void print1DArray(long double *A, int size);
void print2DArray(long double **A, int size);
void reArrangeArray(long double **A, long double *B, int size);
int getSize(char *fileName);


int main(){
   MPI_Init(NULL, NULL);
   int rank;
   int comm_sz;
   int i;
   int size;
   char *fileName = "array.txt";
   MPI_Comm comm;

   long double tolerance = 0.0001;

   comm = MPI_COMM_WORLD;
   MPI_Comm_size(comm, &comm_sz);
   MPI_Comm_rank(comm, &rank);

   if(rank == 0){
      /* Get all of the data and allocate arrays */
      size = getSize(fileName);
   }
   /* Send and recieve size to every other process */
   MPI_Bcast(&size, 1, MPI_INT, 0, comm);

   /* Allocate the arrays */
   long double Ans[size];
   long double **Arr = (long double **) malloc(size * sizeof(long double *));
   for(i = 0; i < size; i++){
      Arr[i] = (long double *) malloc(size * sizeof(long double));
   }

   if(rank == 0){
      readInMatrix(fileName, Arr, Ans);
      reArrangeArray(Arr, Ans, size);
   }


   MPI_Barrier(MPI_COMM_WORLD);
   for(i = 0; i < size; i++){
      MPI_Bcast(Arr[i], size, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
   }

   for(i = 0; i < size; i++){
      Ans[i] = 0;
   }
   
   /* Print answers */
   ParallelJacobian(Arr, Ans, size, tolerance, comm_sz, rank, comm);

   MPI_Barrier(MPI_COMM_WORLD);

   /* Print answers */
   if(rank == 0){
      print1DArray(Ans, size);
   }

   MPI_Barrier(MPI_COMM_WORLD);
   /* Free Memory */
   for(i = 0; i < size; i++){
      free(Arr[i]);
   }
   free(Arr);

   /* Exit program */
   MPI_Finalize();
   return 0;
}

void sendRecieveDebug(long double **Arr, int size, int rank, MPI_Comm comm, int np){
   MPI_Barrier(MPI_COMM_WORLD);
   int i;
   if(rank != 0){
      MPI_Send(&size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
   } else{
      for(i = 1; i < np; i++){
         MPI_Recv(&size, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         printf("%d\n", size);
      }
   }
}

void ParallelJacobian(long double *A[], long double Ans[], int n, long double tolerance, int np, int rank, MPI_Comm comm){
   int start, end;
   int average = n / np;
   int numRuns = 0;
   int toleranceMet;
   int i, j;
   long double gap;

   start = average * rank + (n % np > rank ? rank : n % np);
   end = start + average + (n % np > rank ? 1 : 0);

   int count = end - start;

   long double allAns[n];

   MPI_Barrier(MPI_COMM_WORLD);
   do {
      #pragma omp parallel for num_threads(threadCount) default(none) shared(A, Ans, tolerance, n, np, rank, comm, start, end, average, numRuns, count, allAns, toleranceMet, gap) private(i, j)
      for(i = start; i < end; i++){
         allAns[i] = 0;
         for(j = 0; j < n; j++){
            if(j == i){
               allAns[i] += A[i][j];
            } else {
               allAns[i] += A[i][j] * Ans[j];
            }
         }
      }


      /* Send all your answers to the root node and recieve their answers */
      /* Not positive this is risk free */
      // #pragma omp parallel for num_threads(threadCount) default(none) shared(A, Ans, tolerance, n, np, rank, comm, start, end, average, numRuns, count, allAns, toleranceMet, gap, j) private(i) schedule(static)
      for(i = 0; i < np; i++){
         MPI_Bcast((allAns + (average * i + (n % np > i ? i : n % np))),
          average + (n % np > i ? 1 : 0),
          MPI_LONG_DOUBLE, i, MPI_COMM_WORLD);
      }
      MPI_Barrier(MPI_COMM_WORLD);


      toleranceMet = 1;
      #pragma omp parallel for num_threads(threadCount) default(none) shared(A, Ans, tolerance, n, np, rank, comm, start, end, average, numRuns, count, allAns, toleranceMet, j) private(i, gap)
      for(i = 0; i < n; i++){
         gap = fabs(allAns[i] - Ans[i]);
         if(!RelDif(allAns[i], Ans[i]) <= tolerance){
            #pragma omp critical
            toleranceMet = 0;
         }
         Ans[i] = allAns[i];
      }


      {
      numRuns++;
      }
   } while(!toleranceMet && numRuns < 400);
   if(rank == 0){
      printf("%d\n", numRuns);
   }
   MPI_Barrier(MPI_COMM_WORLD);

   return;
}

long double RelDif(long double a, long double b){
   long double c = ((a) < 0 ? -(a) : (a));
   long double d = ((a) < 0 ? -(a) : (a));
   
   d = ((c) > (d) ? (c) : (d));

   return d == 0.0 ? 0.0 : (((a-b) < 0 ? -1 * (a-b) : a-b)) / d;
}


/* A = old Answers, B = new Answers */
int checkTolerance(long double tolerance, long double *A, long double *B, int size){
   int i;
   int isDone = 1;
   long double gap;
   for(i = 0; i < size; i++){
      gap = A[i] - B[i];
      gap = (gap > 0 ? gap : gap * -1);
      if(gap > tolerance){
         isDone = 0;
      }
      A[i] = B[i];
   }
   return isDone;
}

int getSize(char *fileName){
   FILE *fp;
   int size;
   fp = fopen(fileName, "r");
   if(fp == NULL){
      exit(1);
   }
   fscanf(fp, "%d", &size);
   fclose(fp);
   return size;
}

void readInMatrix(char *fileName, long double *Arr[], long double Ans[]){
   FILE *fp;
   int size, i, j;
   fp = fopen(fileName, "r");
   if(fp == NULL){
      printf("ERROR OPENING FILE!");
      exit(1);
   }

   fscanf(fp, "%d", &size);

   for(i = 0; i < size; i++){
      for(j = 0; j < size; j++){
         fscanf(fp, "%llf", &Arr[i][j]);
      }
   }
   for(i = 0; i < size; i++){
      fscanf(fp, "%llf", &Ans[i]);
   }

   fclose(fp);
}

void reArrangeArray(long double **A, long double *B, int size){
   int i, j;
   long double divideBy;
   for(i = 0; i < size; i++){
      divideBy = A[i][i];
      A[i][i] = B[i] * -1;
      for(j = 0; j < size; j++){
         A[i][j] = (-1 * A[i][j]) / divideBy;
      }
   }
}

void print1DArray(long double *A, int size){
   int j; 
   printf("\n");
   for(j = 0; j < size; j++){
      printf("%llf ", A[j]);
   }
   printf("\n");
}

void print2DArray(long double **A, int size){
   int j;
   int i; 
   for(i = 0; i < size; i++){
      for(j = 0; j < size; j++){
         printf("%llf ", A[i][j]);
      }
      printf("\n");
   }
   printf("\n");
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#define threadCount 5

void readInMatrix(char *fileName, double *Arr[], double Ans[]);
void sendRecieveDebug(double **Arr, int size, int rank, MPI_Comm comm, int np);
int checkTolerance(double tolerance, double *A, double *B, int size);
void ParallelJacobian(double *A[], double Ans[], int n, double tolerance, int np, int rank, MPI_Comm comm);
void print1DArray(double *A, int size);
void print2DArray(double **A, int size);
void reArrangeArray(double **A, double *B, int size);
int getSize(char *fileName);


int main(){
   MPI_Init(NULL, NULL);
   int rank;
   int comm_sz;
   int i;
   int size;
   char *fileName = "array.txt";
   MPI_Comm comm;

   int tolerance = 0.001;

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
   double Ans[size];
   double **Arr = (double **) malloc(size * sizeof(double *));
   for(i = 0; i < size; i++){
      Arr[i] = (double *) malloc(size * sizeof(double));
   }

   if(rank == 0){
      readInMatrix(fileName, Arr, Ans);
      reArrangeArray(Arr, Ans, size);
   }


   MPI_Barrier(MPI_COMM_WORLD);
   for(i = 0; i < size; i++){
      MPI_Bcast(Arr[i], size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
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

void sendRecieveDebug(double **Arr, int size, int rank, MPI_Comm comm, int np){
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

void ParallelJacobian(double *A[], double Ans[], int n, double tolerance, int np, int rank, MPI_Comm comm){
   int i, j;
   int toleranceMet;
   double gap;
   int start, end;
   int average = n / np;
   int numRuns = 0;

   start = average * rank + (n % np > rank ? rank : n % np);
   end = start + average + (n % np > rank ? 1 : 0);

   int count = end - start;

   double allAns[n];

   MPI_Barrier(MPI_COMM_WORLD);
   do {
      if(rank == 0){
   //      print1DArray(Ans, n);
      }
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
      for(i = 0; i < np; i++){
         MPI_Bcast((allAns + (average * i + (n % np > i ? i : n % np))),
          average + (n % np > i ? 1 : 0),
          MPI_DOUBLE, i, MPI_COMM_WORLD);
      }
      MPI_Barrier(MPI_COMM_WORLD);

      toleranceMet = 1;
      for(i = 0; i < n; i++){
         gap = fabs(allAns[i] - Ans[i]);
         if(! (gap <= tolerance)){
            if(rank == 0){
     //          printf("TRUE: %lf %lf %lf\n", allAns[i], Ans[i], gap);
            }
            toleranceMet = 0;
         }
         Ans[i] = allAns[i];
      }

      numRuns++;
   } while(!toleranceMet && numRuns < 1000);
   MPI_Barrier(MPI_COMM_WORLD);

   return;
}


/* A = old Answers, B = new Answers */
int checkTolerance(double tolerance, double *A, double *B, int size){
   int i;
   int isDone = 1;
   double gap;
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

void readInMatrix(char *fileName, double *Arr[], double Ans[]){
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
         fscanf(fp, "%lf", &Arr[i][j]);
      }
   }
   for(i = 0; i < size; i++){
      fscanf(fp, "%lf", &Ans[i]);
   }

   fclose(fp);
}

void reArrangeArray(double **A, double *B, int size){
   int i, j;
   double divideBy;
   for(i = 0; i < size; i++){
      divideBy = A[i][i];
      A[i][i] = B[i] * -1;
      for(j = 0; j < size; j++){
         A[i][j] = (-1 * A[i][j]) / divideBy;
      }
   }
}

void print1DArray(double *A, int size){
   int j; 
   printf("\n");
   for(j = 0; j < size; j++){
      printf("%lf ", A[j]);
   }
   printf("\n");
}

void print2DArray(double **A, int size){
   int j;
   int i; 
   for(i = 0; i < size; i++){
      for(j = 0; j < size; j++){
         printf("%lf ", A[i][j]);
      }
      printf("\n");
   }
   printf("\n");
}

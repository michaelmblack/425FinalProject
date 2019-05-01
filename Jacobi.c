#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int readInMatrix();
void SerialJacobian(double **A, double *Ans, int n, double tolerance);
void print1DArray(double *A, int size);
void reArrangeArray(double **A, double *B, int size);
void print2DArray(double **A, int size);

double **Arr;
double *Ans;
char *fileName = "array.txt";
int main(){
   int i;
   double err = 0.001;
   int size = readInMatrix(); 
   double *Answers = (double *) malloc (size * sizeof(double));
   for(i = 0; i < size; i++){
      Answers[i] = 0;
   }
   reArrangeArray(Arr, Ans, size);
   SerialJacobian(Arr, Answers, size, err);
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

/*
 * A: co-efficient matrix
 * Ans: answer matrix
 * n: size of matrix
 * err: acceptable error distance
 */
void SerialJacobian(double **A, double *Ans, int n, double tolerance){
   int i, j;
   int toleranceMet;
   double gap;
   double v;
   double *newAns = (double *) malloc( n * sizeof(double));
   /* For every value in the matrix do the following */
   do {
      for(i = 0; i < n; i++){
         newAns[i] = 0;
         for(j = 0; j < n; j++){
            if(j == i){
               newAns[i] += A[i][j];
            }
            else{
               newAns[i] += A[i][j] * Ans[j];
            }
         }
      }
      print1DArray(newAns, n);
      toleranceMet = 1;
      for(i = 0; i < n; i++){
         gap = Ans[i] - newAns[i];
         gap = (gap > 0 ? gap : gap * -1);
         if(gap > tolerance){
            toleranceMet = 0;
         }
         Ans[i] = newAns[i];
      }
   } while(!toleranceMet);

}



void print2DArray(double **A, int size){
   int i, j;
   for(i = 0; i < size; i++){
      for(j = 0; j < size; j++){
         printf("%lf ", A[i][j]);
      }
      printf("\n");
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

/*
 * Form:
 * size
 * array
 * answer array
 */
int readInMatrix(){
   FILE *fp;
   int size, i, j;
   fp = fopen(fileName, "r");
   if(fp == NULL){
      printf("ERROR OPENING FILE!");
      exit(1);
   }

   /* Get the size of the array and allocate the arrays */
   fscanf(fp, "%d", &size);
   Ans = (double *) malloc(size * sizeof(double));
   Arr = (double **) malloc(size * sizeof(double *));
   for(i = 0; i < size; i++)
      Arr[i] = (double *) malloc(size * sizeof(double));

   for(i = 0; i < size; i++){
      for(j = 0; j < size; j++){
         fscanf(fp, "%lf", &Arr[i][j]);
      }
   }
   for(i = 0; i < size; i++){
      fscanf(fp, "%lf", &Ans[i]);
   }
   return size;
}

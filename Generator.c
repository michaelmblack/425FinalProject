#include <stdio.h>
#include <stdlib.h>

#define size 12
int main(){
   int i, j, sum;
   int **A = (int **) malloc (size * sizeof(int *));
   for(i = 0; i < size; i++)
      A[i] = (int *) malloc(size * sizeof(int));
   
   
   srand(21);
   for(i = 0; i < size; i++){
      for(j = 0; j < size; j++){
         A[i][j] = rand() % 10;
      }
   }
   
   for(i = 0; i < size; i++){
      sum = 0;
      for(j = 0; j < size; j++){
         if(i != j)
            sum += A[i][j];
      }
      A[i][i] = sum + 1;
   }

   printf("%d\n", size); 
   for(i = 0; i < size; i++){
      for(j = 0; j < size; j++){
         printf("%d ", A[i][j]);
      }
      printf("\n");
   }
   for(i = 0; i < size; i++){
      printf("%d ", rand() % 10);
   }
   printf("\n");

}

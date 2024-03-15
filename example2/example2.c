#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "Timming.h"

int main(int argc, char **argv){
  int n_proc; // # total de procesos
  int my_proc; // El proceso actual
  long unsigned int suma;
  long unsigned int streaming;
  int i,j;
  float suma0;
  int a,b;
  int dx;
  MPI_Status status;
  double utime0, stime0, wtime0,
         utime1, stime1, wtime1;


  uswtime(&utime0, &stime0, &wtime0); //tomando el tiempo  

  MPI_Init (&argc, &argv); /* Inicializar MPI */
  MPI_Comm_rank(MPI_COMM_WORLD,&my_proc); /* Determinar el rango del proceso invocado*/
  MPI_Comm_size(MPI_COMM_WORLD,&n_proc); /* Determinar el numero de procesos */
  MPI_Barrier (MPI_COMM_WORLD);
  
  dx = 1000000/n_proc;
  suma = 0;
  a = my_proc*dx;
  b = (my_proc+1)*dx;
  //printf("%i\t%i\t%i\n",my_proc,a,b);

  for (i=a;i<b;i++){

    suma += i;


  }
 
  printf("Process %i: %lu\n",my_proc,suma);

  MPI_Barrier (MPI_COMM_WORLD);
 
  if (my_proc == 0){ //Master
    // Recolectar las sumas
    for(i=1; i<n_proc;i++){
      MPI_Recv(&streaming, 1, MPI_UNSIGNED_LONG, i, 0, MPI_COMM_WORLD,&status);
      suma += streaming;
    }
    printf("Total %lu\n",suma);
			
  }else{ //Slaves
    //Enviar el resultado
    MPI_Send(&suma, 1, MPI_UNSIGNED_LONG, 0, 0,MPI_COMM_WORLD);
  }



  // code block
 

  //  suma = suma0+suma1+suma2+suma3;
  // printf("total %i\n",suma);

  MPI_Barrier (MPI_COMM_WORLD);
  MPI_Finalize();

  if (my_proc==0){
    uswtime(&utime1, &stime1, &wtime1);
    printf("\nBenchmarks (sec):\n"); 
    printf("real %.3f\n", wtime1 - wtime0); 
    printf("user %.3f\n", utime1 - utime0); 
    printf("sys %.3f\n", stime1 - stime0); 
    printf("\n"); 
    printf("CPU/Wall %.3f %% \n",
	   100.0 * (utime1 - utime0 + stime1 - stime0) / (wtime1 - wtime0));
    printf("\n");
  }
  

  return 0;
}

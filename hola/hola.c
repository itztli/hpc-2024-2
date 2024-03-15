#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

int main(int argc, char **argv){
  int n_proc; // # total de procesos
  int my_proc; // El proceso actual
  int suma;
  int i,j;
  
  MPI_Init (&argc, &argv); /* Inicializar MPI */
  MPI_Comm_rank(MPI_COMM_WORLD,&my_proc); /* Determinar el rango del proceso invocado*/
  MPI_Comm_size(MPI_COMM_WORLD,&n_proc); /* Determinar el numero de procesos */
  MPI_Barrier (MPI_COMM_WORLD);
  switch (my_proc) {
  case 0:
    float suma0;
    suma0=0;
    for (j=0;j<102400;j++){
      for (i=0;i<10000;i++){
	suma0 += sqrt((float)i);
      }
    }
    printf("Process 0: %f\n",suma0);
    // code block
    break;
  case 1:
    float suma1;
    suma1=0;
    for (j=0;j<102400;j++){
      for (i=10000;i<20000;i++){
	suma1 += sqrt((float)i);
      }
    }
    printf("Process 1: %f\n",suma1);
    
    // code block
    break;
  case 2:
    // code block
    float suma2;
    suma2=0;
    for (j=0;j<102400;j++){
      for (i=20000;i<30000;i++){
	suma2 += sqrt((float)i);
      }
    }
    printf("Process 2: %f\n",suma2);

    break;
  case 3:
    // code block
    float suma3;
    suma3=0;
    for (j=0;j<102400;j++){
      for (i=30000;i<40000;i++){
	suma3 += sqrt((float)i);
      }
    }
    printf("Process 3: %f\n",suma3);

    break;
  default:
    // code block
}

  //  suma = suma0+suma1+suma2+suma3;
  // printf("total %i\n",suma);

  MPI_Barrier (MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}

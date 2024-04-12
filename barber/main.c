// Barbero
// by vdelaluz@enesmorelia.unam.mx
// GNU/GPL License
// 20230512

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "myvar.h"

double f(double x){
  return x*x;
}


int main(int argn, char **argc){
  int miproc, numproc;    
  MPI_Status status;
  int data;
  double a,b,dx,sum,F;
  MPI_myvar range;
  
  MPI_Init(&argn, &argc); /* Inicializar MPI */
  MPI_Comm_rank(MPI_COMM_WORLD,&miproc); /* Determinar el rango del proceso invocado*/
  MPI_Comm_size(MPI_COMM_WORLD,&numproc); /* Determinar el numero de procesos */
  MPI_Barrier (MPI_COMM_WORLD);


  if (miproc == 0){ //master    
    if (argn < 4){
      printf("Faltan parametros [a] [b] [dt]\n");
      return 0;
    }
    if (sscanf(argc[1],"%lf",&a) != 1){
      printf("Error al convertir a.\n");
      return 0;
    }
    if(sscanf(argc[2],"%lf",&b) != 1){
      printf("Error al convertir b.\n");
      return 0;
    }
    if(sscanf(argc[3],"%lf",&dx)!=1){
      printf("Error al convertir dx\n");
      return 0;
    }
    printf("[%lf, %lf] dx=%lf\n",a,b,dx);
  } //master reading command line

  MPI_Barrier (MPI_COMM_WORLD);
    

  
  if (miproc != 0) { // slaves
    //int flag_start = 1;
    double F = 0.0; //result
    data = 1;
    range.F = 0.0;
    
    while(1){
    MPI_Send(&range, sizeof(range), MPI_CHARACTER, 0, 0, MPI_COMM_WORLD); 
      MPI_Recv(&range, sizeof(range), MPI_CHARACTER, 0, 0, MPI_COMM_WORLD, &status);
      //integral de Riemman
      F = f(range.a)*range.dx;
      range.F = F;
      //printf("%i:[%lf,%lf] dx=%lf F=%lf\n",miproc,range.a,range.b,range.dx,range.F);
      //Parallel processing
    }

  }else{ //Master
    int flag = -1;
    MPI_Request request;
    double sum = 0;
    int n=0;
    while (1) {
      if(flag != 0){
	MPI_Irecv(&range, sizeof(range), MPI_CHARACTER, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &request); 
        flag = 0;
	
      }
      MPI_Test(&request, &flag, &status);
      if (flag != 0) {
	if (status.MPI_SOURCE != -1){
	  // sending information
	  // segmentar la informacion para enviarla al nodo disponible

	  //printf("%i: %lf\n",status.MPI_SOURCE,range.F);
	  sum += range.F;

	  range.a = a + n*dx;
	  range.b = a + (n+1)*dx;
	  range.dx = dx;
	  range.F = 0.0;
	  n++;
	  MPI_Send(&range, sizeof(range), MPI_CHARACTER, status.MPI_SOURCE, 0, MPI_COMM_WORLD); 	
	}	
      flag = -1;
      }

      //stop condition
      //if ((a+n*dx) >= b){
      if (n == numproc + (int)((b-a)/dx) ){
	printf("F=%le\n",sum);
	break;
      }
    } //while(1)      
  }//Master

   MPI_Abort(MPI_COMM_WORLD,MPI_SUCCESS);

  //MPI_Finalize ();
  return 0;
}

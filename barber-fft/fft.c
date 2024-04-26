#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "myvar.h"
#include "Timming.h"

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#ifndef PI
#define PI 3.14159265358979323846
#endif


/*Replaces data[1..2*nn] by its discrete Fourier transform, if isign is input as 1; or replaces data[1..2*nn] by nn times its inverse discrete Fourier transform, if isign is input as −1. data is a complex array of length nn or, equivalently, a real array of length 2*nn. nn MUST be an integer power of 2 (this is not checked for!).*/
void four1(float data[], unsigned long nn, int isign){
  unsigned long n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta;
  float tempr,tempi;

  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2) {
    if (j > i) {
      SWAP(data[j],data[i]);
      SWAP(data[j+1],data[i+1]);
    }
    m=nn;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }

  mmax=2;
  while (n > mmax) {
    istep=mmax << 1;
    theta=isign*(6.28318530717959/mmax);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;

    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) {
      for (i=m;i<=n;i+=istep) {
	j=i+mmax;
	tempr=wr*data[j]-wi*data[j+1];
	tempi=wr*data[j+1]+wi*data[j];
	data[j]=data[i]-tempr;
	data[j+1]=data[i+1]-tempi;
	data[i] += tempr;
	data[i+1] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}

/*
  Calculates the Fourier transform of a set of n real-valued data points. Replaces this data (which is stored in array data[1..n]) by the positive frequency half of its complex Fourier transform.
  The real-valued ﬁrst and last components of the complex transform are returned as elements data[1] and data[2], respectively. n must be a power of 2. This routine also calculates the inverse transform of a complex data array if it is the transform of real data. (Result in this case must be multiplied by 2/n.)*/
void realft(float data[], unsigned long n, int isign){
  //void four1(float data[], unsigned long nn, int isign);
  unsigned long i,i1,i2,i3,i4,np3;
  float c1=0.5,c2,h1r,h1i,h2r,h2i;
  double wr,wi,wpr,wpi,wtemp,theta;
  //Double precision for the trigonometric recurrences.
  theta=3.141592653589793/(double) (n>>1);
  //Initialize the recurrence.
  if (isign == 1) {
    c2 = -0.5;
    four1(data,n>>1,1);
    //The forward transform is here.
  } else {
    c2=0.5;
    //Otherwise set up for an inverse transform
    theta = -theta;
  }
  wtemp=sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi=sin(theta);
  wr=1.0+wpr;
  wi=wpi;
  np3=n+3;
  for (i=2;i<=(n>>2);i++) {
    //Case i=1 done separately below.
    i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
    h1r=c1*(data[i1]+data[i3]);
    //The two separate transforms are separated out of data.
    h1i=c1*(data[i2]-data[i4]);
    h2r = -c2*(data[i2]+data[i4]);
    h2i=c2*(data[i1]-data[i3]);
    data[i1]=h1r+wr*h2r-wi*h2i;
    //Here they are recombined to form the true transform of the original real data.
    data[i2]=h1i+wr*h2i+wi*h2r;
    data[i3]=h1r-wr*h2r+wi*h2i;
    data[i4] = -h1i+wr*h2i+wi*h2r;
    wr=(wtemp=wr)*wpr-wi*wpi+wr;
    //The recurrence.
    wi=wi*wpr+wtemp*wpi+wi;
  }
  if (isign == 1) {
    data[1] = (h1r=data[1])+data[2];
    data[2] = h1r-data[2];
  } else {
    data[1]=c1*((h1r=data[1])+data[2]);
    data[2]=c1*(h1r-data[2]);
    four1(data,n>>1,-1);
  }
}


char* Clean(char *linea){
  int n = strlen(linea);
  int i;
  char *rstr;
  rstr = malloc(sizeof(char)*64);
  for (i=0; i<n; i++){
    if (linea[i] == '#'){
      if (i==0){
	return "";
      }else{
	strncpy(rstr,linea,i);
	return rstr;
      }
    }
  }
  return linea;
}




//2^n = 1024  ->   n = log_2(1024) = 10
int main(int argc, char **argv){
  unsigned long nn = 10;
  int isign = 1;
  float theta = 0.0;
  float frequency = 120.0;
  float sampling = 4096.0; // # muestras x segundo
  float dtheta = frequency*2.0*PI/sampling;
  //float signal[((int)sampling)*2];
  float *signal;
  float *subset;
  int i,itera;
  char comando[256];
  char inputfile[256];
  char outputfile[256];
  char window_str[256];
  int signal_flag=0;
  int spectrum_flag=0;
  int synthetic_flag=0;
  int inputfile_flag=0;
  FILE *file;
  char *line;
  char *clean;
  int N;
  int window = 1024;
  MPI_Status status;
  int n_proc; // # total de procesos
  int my_proc; // El proceso actual
  //MPI_myvar data_proc;
  int index;
  double utime0, stime0, wtime0,
    utime1, stime1, wtime1,
    utime2, stime2, wtime2;
  
  MPI_Init (&argc, &argv); /* Inicializar MPI */
  MPI_Comm_rank(MPI_COMM_WORLD,&my_proc); /* Determinar el rango del proceso invocado*/
  MPI_Comm_size(MPI_COMM_WORLD,&n_proc); /* Determinar el numero de procesos */
  MPI_Barrier (MPI_COMM_WORLD);
  
  if (my_proc==0){

    
    for (i=1; i<argc;i++){

      sprintf(comando,"%s",argv[i]);
      
      if (strcmp(comando,"-signal") == 0){
	signal_flag = 1;
      }
      
      if (strcmp(comando,"-spectrum") == 0){
	spectrum_flag = 1;
      }
      if (strcmp(comando,"-synthetic") == 0){
	synthetic_flag = 1;
      }
      if (strcmp(comando,"-help") == 0){
	printf("Flags: -signal, -spectrum, -synthetic, -file <archive>\n");
	return 0;
      }
      
      if (strcmp(comando,"-file") == 0){
        if (sprintf(inputfile,"%s",argv[++i]) > 0){
          //printf(".");
	  inputfile_flag=1;
	  synthetic_flag = 0;
        }else{
	  printf("Error: -file <archive>\n");
          return 0;
        }
      }
      
      if (strcmp(comando,"-window") == 0){
        if (sprintf(window_str,"%s",argv[++i]) > 0){
          //printf(".");
	  if (sscanf(window_str,"%i", &window) != 1){
	    printf("Error: -window <size window (2^n)>\n");
	    return 0;
	  }
	  // falta checar que window = 2^n
        }else{
	  printf("Error: -window <size window (2^n)>\n");
          return 0;
        }
      }
      
    }
    
    if (synthetic_flag){
      for (int i=0;i<((int)sampling);i++){
	signal[2*i] = cos(theta)+cos(0.5*theta)+cos(10.0*theta);
	signal[(2*i)+1] = sin(theta)+sin(0.5*theta)+sin(10.0*theta);
	theta += dtheta;
      }
    }
    
    if (inputfile_flag){
      file = fopen(inputfile,"r");
      if (file == NULL){
	printf("Error 1: File %s not found.\n",inputfile);
	exit(0);
      }

      line = malloc(sizeof(char)*64);
      N=0;
      uswtime(&utime2, &stime2, &wtime2); //tomando el tiempo  
      printf("Calculating memory allocation.\n");
      while (fgets(line, 64, file ) != NULL){
	//clean = Clean(line);
	//if (strlen(clean)>0)
	N++;
	//free(clean);
      }
      printf("%fMB\n",sizeof(float)*N*2/1024.0/1024.0);
      signal = malloc(sizeof(float)*N*2);
      
      rewind(file);

      printf("Reading data.\n");
      i=0;
      while (fgets(line, 64, file ) != NULL){
	//clean = Clean(line);
	//if (strlen(clean)>0){
	sscanf(line,"%f", &signal[2*i]);
	signal[(2*i)+1] = 0.0;
	i++;
	//}
	//free(clean);
      }
      printf("Ready!\n");
      free(line);
      fclose(file);

      uswtime(&utime0, &stime0, &wtime0); //tomando el tiempo  
	printf("\nBenchmarks (sec):\n"); 
	printf("real %.3f\n", wtime0 - wtime2); 
	printf("user %.3f\n", utime0 - utime2); 
	printf("sys %.3f\n", stime0 - stime2); 
	printf("\n"); 
	printf("CPU/Wall %.3f %% \n",
	       100.0 * (utime0 - utime2 + stime0 - stime2) / (wtime0 - wtime2));
	printf("\n");

      
      /*  
      for(i=0;i<N;i++){
	printf("%f,",signal[i]);
      }
      printf("\n");
      */
  
  
  
    if (signal_flag){
      for (int i=0;i < ((int)sampling) ;i++){
	printf("%f\n",signal[2*i]);
      }
    }
  }

  } // master MPI


  MPI_Barrier (MPI_COMM_WORLD);


  //send window

  if (my_proc==0){
    if (spectrum_flag){
      subset = malloc(sizeof(float)*window*2);
      //itera =  N / window;
      //if (itera > n_proc){
      //printf("Warning: iterations greather [%i] than numproc [%i]\n",itera, n_proc);
      //}
      for (i=0; i < n_proc-1;i++){
	//memcpy(subset, &signal[2*i*window], window*2*sizeof(float));
	//printf("Sending window %i\n",window);
	MPI_Send(&window, 1, MPI_INT, i+1, 0,MPI_COMM_WORLD);
	//MPI_Barrier (MPI_COMM_WORLD);	
	//printf("0 > %i subset[0]=%f\n",i+1, subset[2048]);
	//MPI_Send(subset, 2*window+1, MPI_FLOAT, i+1, 0,MPI_COMM_WORLD);
      }
    }
    //free(subset);
    //free(signal);
  }else{ // master MPI
    MPI_Recv(&window, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,&status);
    //printf("Receving window %i\n",window);
    //MPI_Barrier (MPI_COMM_WORLD);
    subset = malloc(sizeof(float)*window*2);
    //MPI_Recv(subset, 2*window+1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD,&status);
    //printf("%i < 0 subset[0]=%f\n",my_proc, subset[2048]);
    //isign=1;
    //realft(subset, 2*window, isign);
    //printf("%i < 0 fft[10]=%f\n",my_proc, subset[2048]);
    
    //sprintf(outputfile,"spectrum-mpi-%i.dat",my_proc-1);
    //file = fopen(outputfile,"w");
    //for (int i=0;i< (int)window; i++){
    //  fprintf(file,"%f\n",subset[2*i]);
    //}
    //free(subset);
    //fclose(file);
  }

  MPI_Barrier (MPI_COMM_WORLD);
  //MPI_Finalize();






  
  if (my_proc != 0) { // slaves
    //int flag_start = 1;
    //double F = 0.0; //result
    //data = 1;
    //range.F = 0.0;
    index=0;
    //data_proc.subset = subset;
    while(1){
      //printf("%i:CODE1\n",my_proc);
      MPI_Send(&index, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD);
      //printf("%i:CODE2\n",my_proc);
      MPI_Recv(&index, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(subset, 2*window+1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status);
      //printf("%i:Processing FFT\n",my_proc);
      isign=1;

      realft(subset, 2*window, isign);

      sprintf(outputfile,"output/spectrum-mpi-%i.dat",index);
      file = fopen(outputfile,"w");
      for (int i=0;i< (int)window; i++){
	fprintf(file,"%f\n",subset[2*i]);
      }
      //fflush(file);
      fclose(file);


      //integral de Riemman
      //F = f(range.a)*range.dx;
      //range.F = F;
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
	//printf("A:receiving to...\n");
	MPI_Irecv(&index, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
	//printf("B:ok\n");
        flag = 0;
	
      }
      MPI_Test(&request, &flag, &status);
      if (flag != 0) {
	if (status.MPI_SOURCE != -1){
	  // sending information
	  // segmentar la informacion para enviarla al nodo disponible
	  //printf("C:copying\n");

	  memcpy(subset, &signal[2*n*window], window*2*sizeof(float));
	  
	  //printf("C:ok %f\n",subset[0]);
	  index = n;
	  //data_proc.subset = subset;
	  //printf("B:sending to...\n");
	  MPI_Send(&index, 1, MPI_INTEGER, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
	  MPI_Send(subset, 2*window+1, MPI_FLOAT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
	  //printf("B:Ready!\n");
	  n++;
	}	
      flag = -1;
      }

      //stop condition
      //if ((a+n*dx) >= b){
      if ((n)*window >= N ){
	printf("Finish: data computed:%i total:%i!\n",(n)*window,N);
	printf("Spectrum computed: %i\n",n);

	uswtime(&utime1, &stime1, &wtime1);
	printf("\nBenchmarks (sec):\n"); 
	printf("real %.3f\n", wtime1 - wtime0); 
	printf("user %.3f\n", utime1 - utime0); 
	printf("sys %.3f\n", stime1 - stime0); 
	printf("\n"); 
	printf("CPU/Wall %.3f %% \n",
	       100.0 * (utime1 - utime0 + stime1 - stime0) / (wtime1 - wtime0));
	printf("\n");
	break;
      }
    } //while(1)      
  }//Master

   MPI_Abort(MPI_COMM_WORLD,MPI_SUCCESS);















  
  

  return 0;
}

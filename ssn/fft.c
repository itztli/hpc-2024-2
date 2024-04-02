#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
  int i;
  char comando[256];
  char inputfile[256];
  int signal_flag=0;
  int spectrum_flag=0;
  int synthetic_flag=0;
  int inputfile_flag=0;
  FILE *file;
  char *line;
  char *clean;
  int N;

  
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
    while (fgets(line, 64, file ) != NULL){
      clean = Clean(line);
      if (strlen(clean)>0)
	N++;
    }
    signal = malloc(sizeof(float)*N);

    rewind(file);

    i=0;
    while (fgets(line, 64, file ) != NULL){
      clean = Clean(line);
      if (strlen(clean)>0){
	sscanf(clean,"%f", &signal[i]);
	i++;
      }
    }
    
    free(line);
    fclose(file);

    for(i=0;i<N;i++){
      printf("%f,",signal[i]);
    }
    printf("\n");

  }
  
  
  if (signal_flag){
    for (int i=0;i < ((int)sampling) ;i++){
      printf("%f\n",signal[2*i]);
    }
  }

  //printf("Applying FFT\n");
  //four1(signal, nn, isign);
  if (spectrum_flag){
    realft(signal, (int)sampling * 2, isign);
    
    for (int i=0;i< (int)sampling; i++){
      printf("%f\n",signal[2*i]);
    }
  }

  free(signal);
  return 0;
}

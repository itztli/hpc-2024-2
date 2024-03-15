#include <math.h>
#include <stdio.h>
#include <stdlib.h>
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






//2^n = 1024  ->   n = log_2(1024) = 10
int main(){
  float signal[1024];
  unsigned long nn = 10;
  int isign = 1;
  float theta = 0.0;
  float dtheta = 60.0*PI/1024.0;

  for (int i=0;i<512;i++){
    signal[2*i] = cos(theta);
    signal[(2*i)+1] = sin(theta);
    theta += dtheta;
  }
  
  //for (int i=0;i<1024;i++){
  //  printf("%f\n",signal[i]);
  //}

  //printf("Applying FFT\n");
  //four1(signal, nn, isign);
  realft(signal, 1024, isign);

    for (int i=0;i<1024;i++){
    printf("%f\n",signal[i]);
  }
  
  return 0;
}

#include <stdio.h>
#include <stdlib.h>

int main(){
  FILE *f;
  int index;
  float data;
  
  f = fopen("bigdata.dat","w");

  for(int k=0;k<1024;k++){ // ~Gb
    for(int j=0;j<1024;j++){ // ~MB
      data= 0.0;  
      for(int i=0;i< 1024;i++){ // ~kB
	fprintf(f,"%i\t%f\n",index,data);
	data = data + 0.1;
      }
    }
  }

  
  fclose(f);
  return 0;
}

#include <stdio.h>
#include <stdlib.h>

int main(){
  FILE *f;
  char *line = NULL;
  size_t len = 0;
  ssize_t read;
  long unsigned int count = 0;
  int *index;
  float *data;
  float result;

  printf("size of count %zu\n",sizeof(count));

  printf("Reading data\n");
  f = fopen("bigdata.dat","r");
  while ((read = getline(&line, &len, f)) != -1) {
    //printf("Retrieved line of length %zu:\n", read);
    //printf("%s", line);
    if (line[0] != '#'){
      count++;
    }
  }
  fclose(f);
  
  printf("Ready: %lu\n",count);
  index = malloc(count*sizeof(int));
  data = malloc(count*sizeof(float));

  count=0;

  f = fopen("bigdata.dat","r");
  while ((read = getline(&line, &len, f)) != -1) {
    //printf("Retrieved line of length %zu:\n", read);
    //printf("%s", line);
    if (line[0] != '#'){
      sscanf(line, "%i\t%f",&index[count],&data[count]);
      count++;
    }
  }
  fclose(f);

  printf("%f\n",data[3]);

  for(int i=0; i<count; i++){
    data[i] = data[i]*data[i];
    printf("%f\n",data[i]);
  }

  result = 0.0;
  for(int i=0; i<count; i++){
    result += data[i];
  }

  printf("%f\n",result);
    
  free(index);
  free(data);
  
  return 0;
}

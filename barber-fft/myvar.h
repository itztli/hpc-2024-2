#ifndef MYVAR
#define MYVAR

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


typedef struct{
  int index;
  float *subset;
}MPI_myvar;
  

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* myvar */

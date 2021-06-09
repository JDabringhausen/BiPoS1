#include "Ptr.h"

float *fvector(unsigned int elmnts)
{
  float *ptr;

  ptr=(float *)malloc(elmnts*sizeof(float));

  return(ptr);
}

double *dvector(unsigned int elmnts)
{
  double *ptr;

  ptr=(double *)malloc(elmnts*sizeof(double));

  return(ptr);
}

float **fmatrix(unsigned int zeilen,unsigned int spalten)
{
  unsigned int i;
  float **ptr;

  ptr=(float **)malloc(zeilen*sizeof(float *));
  if(!ptr) fprintf(stderr,"Allocation failure (rows)!");
  for(i=0;i<zeilen;i++) {
    ptr[i]=(float *)malloc(spalten*sizeof(float));
    if(!ptr[i]) fprintf(stderr,"Allocation failure (column %i)!",i);
  }

  return(ptr);
}

double **dmatrix(unsigned int zeilen,unsigned int spalten)
{
  unsigned int i;
  double **ptr;

  ptr=(double **)calloc(zeilen,sizeof(double *));
  if(!ptr) fprintf(stderr,"Allocation failure (rows)!");
  for(i=0;i<zeilen;i++) {
    ptr[i]=(double *)calloc(spalten,sizeof(double));
    if(!ptr[i]) fprintf(stderr,"Allocation failure (column %i)!",i);
  }

  return(ptr);
}

long double **ldmatrix(unsigned int zeilen,unsigned int spalten)
{
  unsigned int i;
  long double **ptr;

  ptr=(long double **)malloc(zeilen*sizeof(long double *));
  if(!ptr) fprintf(stderr,"Allocation failure (rows)!");
  for(i=0;i<zeilen;i++) {
    ptr[i]=(long double *)malloc(spalten*sizeof(long double));
    if(!ptr[i]) fprintf(stderr,"Allocation failure (column %i)!",i);
  }

  return(ptr);
}

void fmatrix_free(float **ptr,unsigned int zeilen)
{
  unsigned int i;

  for(i=0;i<zeilen;i++)
    free(ptr[i]);
  free(ptr);
}

void dmatrix_free(double **ptr,unsigned int zeilen)
{
  unsigned int i;

  for(i=0;i<zeilen;i++)
    free(ptr[i]);
  free(ptr);
}

void ldmatrix_free(long double **ptr,unsigned int zeilen)
{
  unsigned int i;

  for(i=0;i<zeilen;i++)
    free(ptr[i]);
  free(ptr);
}

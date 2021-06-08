#include <stdlib.h>
#include <stdio.h>

float *fvector(unsigned int elmnts);
/***********************************************
 * returns a float vector with elmnts entries  *
 ***********************************************/
double *dvector(unsigned int elmnts);
/***********************************************
 * same as *fvector but returns double pointer *
 ***********************************************/
float **fmatrix(unsigned int zeilen,unsigned int spalten);
/***********************************************
 * returns a float matrix with 'zeilen' rows   *
 * and 'spalten' columns: ptr[zeilen][spalten] *
 ***********************************************/
double **dmatrix(unsigned int zeilen,unsigned int spalten);
/***********************************************
 * same as **fmatrix but returns double matrix *
 ***********************************************/
long double **ldmatrix(unsigned int zeilen,unsigned int spalten);
/***********************************************
 * same as **dmatrix but for long double matrix*
 ***********************************************/
void fmatrix_free(float **ptr,unsigned int zeilen);
/***********************************************
 * frees the float matrix assigned in routine  *
 * **fmatrix                                   *
 ***********************************************/
void dmatrix_free(double **ptr,unsigned int zeilen);
/***********************************************
 * frees the double matrix assigned in routine *
 * **dmatrix                                   *
 ***********************************************/
void ldmatrix_free(long double **ptr,unsigned int zeilen);
/***********************************************
 * frees the long double matrix assigned in    *
 * routine **ldmatrix                          *
 ***********************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int sort_desc(float *x1, float *x2);
int sort_asc(double *x1, double *x2);
/**********************************************
 * ascending and descending sorting for qsort *
 *********************************************/
float fabsolute(float x);
double dabsolute(double x);
/**********************************************
 * returnes absolute value of x               *
 *********************************************/
double dmax(double x,double y);
double dmin(double x,double y);
/**********************************************
 * returnes max or min  of x and y            *
 *********************************************/
float fround(float Zahl, unsigned int Stellen);
double dround(double Zahl, unsigned int Stellen);
/**********************************************
 * rounds 'Zahl' to 'Stellen' digits          *
 *********************************************/
/* double *convol(double *y1,double *y2,unsigned int elmnts); */
/******************************************
* convolves y1 with y2 and returns result *
* all pointers have elmnts Elements       *
*******************************************/
double dsum(double *ptr,unsigned int elmnts);
double dsum2(double *ptr,unsigned int elmnts);

#include "Small_routines.h"
#include "Ptr.h"

int sort_desc(float *x1, float *x2) {
    if(*x1>*x2) return -1;
    return 1;
}

int sort_asc(double *x1, double *x2) {
    if(*x1<*x2) return -1;
    return 1;
}

float fabsolute(float x) {
  return (x >= 0) ? x : (-x);
}

double dabsolute(double x) {
  return (x >= 0) ? x : (-x);
}

double dmax(double x,double y) {
  return (x >= y) ?  x : y;
}

double dmin(double x,double y) {
  return (x >= y) ?  y : x;
}

float fround(float Zahl, unsigned int Stellen)
{
    Zahl *= pow(10, Stellen);
    if (Zahl >= 0) {
      Zahl=floor(Zahl+0.5);
    } else {
      Zahl=ceil(Zahl-0.5);
    }
    Zahl /= pow(10, Stellen);
    return Zahl;
}

double dround(double Zahl, unsigned int Stellen)
{
    Zahl *= pow(10, Stellen);
    if (Zahl >= 0) {
      Zahl=floor(Zahl+0.5);
    } else {
      Zahl=ceil(Zahl-0.5);
    }
    Zahl /= pow(10, Stellen);
    return Zahl;
}

double dsum(double *ptr,unsigned int elmnts) {
  unsigned int i;
  double sum=0.;

  for(i=0;i<elmnts;i++) {
    sum+=ptr[2*i+1];
  }

  return sum;
}

double dsum2(double *ptr,unsigned int elmnts) {
  unsigned int i;
  double sum=0.;

  for(i=0;i<elmnts;i++) {
    sum+=ptr[i];
  }

  return sum;
}

/*
double *convol(double *y1,double *y2,unsigned int elmnts) {
  unsigned int i,j,zdim;
  double *z;

  zdim=2*elmnts-1;
  z=dvector(zdim);
  for(i=0;i<zdim;i++)
    z[i] = 0.;

  for(i=0;i<zdim;i++)
    for(j=0;j<elmnts;j++)
      if((i-j)>=0 && (i-j)<elmnts)
	z[i] += y1[j]*y2[i-j];

  return(z);
}
*/

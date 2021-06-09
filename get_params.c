#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "Ptr.h"
#include "Small_routines.h"

/* const_bin_36_half_s */
void get_params(double rho,double *Ecut,double *slope,double *asymptote,double t) {
  double e11,e12,s11,s12,s13,a11,a12;
  double e31,e32,s31,s32,s33,a31,a32;
  double e51,e52,s51,s52,s53,a51,a52;
  double e1,e2,s1,s2,s3,a1,a2;
  double logrho,rholim=378000.;

  e11=-5.00;e12=0.65;
  s11=2.36;s12=2.34;s13=0.89;
  a11=2.15;a12=-0.13;

  e31=-4.64;e32=0.59;
  s31=1.24;s32=1.70;s33=0.79;
  a31=2.04;a32=-0.11;

  e51=-4.40;e52=0.54;
  s51=1.47;s52=1.35;s53=0.82;
  a51=1.97;a52=-0.10;

  if(t==1) { 
    e1=e11;e2=e12;
    s1=s11;s2=s12;s3=s13;
    a1=a11;a2=a12;
  }
  else if(t==3) { 
    e1=e31;e2=e32;
    s1=s31;s2=s32;s3=s33;
    a1=a31;a2=a32;
  }
  else if(t==5 || rho>rholim) {
    e1=e51;e2=e52;
    s1=s51;s2=s52;s3=s53;
    a1=a51;a2=a52;
  } else {
    e1=-5.15+0.15*t;e2=0.6775-0.0275*t;
    s1=2.58-0.445*t;s2=2.5875-0.2475*t;s3=0.903-0.035*t;
    a1=2.195-0.045*t;a2=-0.1375+0.0075*t;
  }

  logrho=log10(rho);

  *Ecut=e1+e2*logrho;
  if(*Ecut<-3.2) /* check for nan */
    *Ecut=-3.2;

  *slope=-1.0/(exp(s1*(logrho-s2)))-s3;

  *asymptote=a1+a2*logrho;
  if(*asymptote>2.)
    *asymptote=2.;

}

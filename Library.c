#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<string.h>
#include<stdbool.h>
#include<sys/stat.h>
#include "Constants.h"

#define max(a,b) (a < b) ?  (b) : (a)

double min_allowed_mr(double m,bool ordered);
double eccentric_anomaly(double t,double T,double e);
double true_anomaly(double u,double e);
double usinu(double u,double e);
int make_binary(double m1, double r1[3], double v1[3],
		double m2, double r2[3], double v2[3],
		double a, double e, double phase);
int rotate(double r[3],double v[3],double alpha,double beta,double gamma);
double apparent_separation(double r1[3],double r2[3]);
void hammer_aitov_projection(double phi,double theta,double *x,double *y);

void Library(FILE *dat,long int Nlib,double mmin,double mmax,bool eigen,bool ordered) {
  clock_t t1, t2;
  t1 = clock();

  long int i,j,count;
  double xxx,x; //y
  double c1, c2, k1, k2;
  double mstar1, mstar2, mtot;
  float *marray=malloc(2*Nlib*sizeof(float));
  double lP,P,lE,E,lL,L,a,q,e1,e2,Rperi,rho,rhostar;
  double t,u,phase,ls,r1[3],r2[3],v1[3],v2[3];
  double phi,theta,gamma;
  srand48((unsigned) time(NULL));

  fprintf(stdout,"Generating library...");fflush(stdout);

  c1 = 1.0-ALPHA1;
  c2 = 1.0-ALPHA2;

  /*
    IMF parametrized as \xi = k times a_i times m_i,
    choose a_1 = 2, then a_2 = 1 from continuity at m = 0.5 Msun
  */

  k1 = 2.0/c1*(pow(0.5,c1)-pow(MLOW,c1)); /* 1/k for 1st power-law segment [MLOW,0.5] with a_1 = 2.0 */

  if (MLOW>0.5) { /* single power-law IMF above 0.5 Msun */
    k1 = 0.;
    k2 = 1.0/c2*(pow(MHIGH,c2)-pow(MLOW,c2)); /* 1/k */
  } else
    k2 = k1 + 1.0/c2*(pow(MHIGH,c2)-pow(0.5,c2)); /* if 2-part power-law add power-law segment [0.5,MHIGH] to coefficient */

  if (MHIGH<0.5) { /* single power law below 0.5 Msun */
    k1 = 2.0/c1*(pow(MHIGH,c1)-pow(MLOW,c1));
    k2 = k1;
  }

  fprintf(dat,"#m1\tm2\te\ta\t\tlP\tlE\tlL\ts\tt/T\tphase");

  /* erzeuge ein Massen-array */
  for(i=0;i<2*Nlib;i++) {

    do {
      x = drand48();
      if (x<k1/k2) /* always true if IMF is single power-law below 0.5 Msun; always false if IMF is single power-law above 0.5 */
	marray[i] = pow(0.5*c1*x*k2+pow(MLOW,c1),1.0/c1);
      else
	marray[i] = pow(c2*(x*k2-k1)+pow(max(0.5,MLOW),c2),1.0/c2);
    } while(marray[i]>mmax || marray[i]<mmin);

  }

  for(i=0;i<2*Nlib-1;i++) {

    if(i==2*Nlib-2) fprintf(stderr,"%li...",i);fflush(stderr);

    mstar1=0.;
    mstar2=0.;
    j=0;

    do {

      if(marray[i+j]!=0) {

	mstar1=marray[i+j];
	marray[i+j]=0.0; /* setze Eintrag auf Null, wenn vergeben */

      } else {

	j++;

      }

    } while(mstar1==0.);

    i += j; /* <- eintrag in marray von mstar1, davor sind nur Nullen! */
    j = 1;

    do {

      if(i+j<2*Nlib-1) {

	  q = marray[i+j] / mstar1;
	  if(q > 1.0) q = 1.0 / q;

	  if( marray[i+j] != 0 && q > min_allowed_mr(mstar1,ordered) && q > min_allowed_mr(marray[i+j],ordered) ) {

	    mstar2=marray[i+j];
	    marray[i+j]=0.0;

	  } else {

	    j++;

	  }

      } else { /* if none of the remaining masses match the q-constraint, choose next non-zero mass */

	j=1;
	do {

	  if(marray[i+j]!=0) {

	    mstar2=marray[i+j];
	    marray[i+j]=0.0;

	  } else {

	    j++;

	  }

	} while(mstar2==0.);

    }

  } while(mstar2==0.);

    // if(i==2*Nlib-2) fprintf(stderr,"3...");fflush(stderr);

    mtot=mstar1+mstar2;

    /* select from birth period BDF */
    do {
      xxx = drand48();
      lP=log10(PMIN)+sqrt(DELTA*(exp(2.0*xxx/ETA)-1.0));
    } while ((lP>LPMAX) || (lP<log10(PMIN)));

    P = pow(10,lP); /* d */

    /* birth eccentricity from thermalized distribution */
    xxx = drand48();
    e1=sqrt(xxx);
    e2=e1;

    /* Kroupa 1995b eigenevolution */
    if(eigen) {

      /* eigenvolved eccentricity */
      Rperi=(1-e1)*pow(P/365.25,2.0/3.0)*pow(mtot,1.0/3.0);
      rho=pow(LAMBDA*RSUN/Rperi,CHI);
      if(e1>0)
	e2=pow(exp(1.0),(-rho+log(e1))); /* natural log! */
      else
	e2=e1;
      /* birth mass ratio */
      q=mstar1/mstar2;
      if(q>1.0) q=1.0/q;
      /* evolved mass-ratio */
      if(rho<=1.0) rhostar=rho;
      else rhostar=1.0;
      q=q+(1.0-q)*rhostar;
      /* evolved mass of secondary, primary mass does not change */
      if(mstar1>=mstar2) mstar2=q*mstar1;
      else mstar1=q*mstar2;
      P = P*sqrt(mtot/(mstar1+mstar2))*pow((1-e1)/(1-e2),3.0/2.0);

    }

    lP = log10(P); /* log P [d] */
    P = P*86400; /* s */
    a = pow(GN*P*P*(mstar1+mstar2)/(4*PI*PI),1./3.); /* km */
    E = GN*mstar1*mstar2/(2.*a); /* Msun,km,s */
    lE = log10(E);
    L = sqrt(GN*a*(1.0-e2*e2)/(mstar1+mstar2))*mstar1*mstar2;
    lL = log10(1e10*L/(mstar1+mstar2)); /* cm**2 s**-1 */

    a /= 149597870.691; /* AU */
    P /= 86400.0; /* days */

    /*****************/
    /* orient binary randomly to obtain apparent separation */
    t = drand48() * P; /* orbital time in days */
    u = eccentric_anomaly(t,P,e2);
    phase = true_anomaly(u,e2);

    make_binary(mstar1,r1,v1,mstar2,r2,v2,a,e2,phase);

    phi = drand48() * (2*PI); /* equally distributed in [0:2*PI] */
    theta = drand48() * 2.0 - 1.0; /* equally distributed in [-1:1] */
    theta = 1.0 * acos( theta ); /* cos ist gleichverteilt zwischen [-1:1] */
    gamma = drand48() * (2*PI); /* equally distributed in [0:2*PI] */

    //hammer_aitov_projection(phi,theta,&x,&y);

    rotate(r1,v1,phi,theta,gamma);
    rotate(r2,v2,phi,theta,gamma);

    ls = apparent_separation(r1,r2);
    /*****************/

    fprintf(dat,"\n% .4lf\t% .4lf\t% .4lf\t% .4lf\t% .4lf\t% .3lf\t% .3lf\t% .4lf\t% .4lf\t% .4lf",mstar1,mstar2,e2,log10(a),lP,lE,lL,ls, (t/P) ,phase);

  }

  fprintf(stdout,"done.");fflush(stdout);

  count=0;
  for(i=0;i<2*Nlib;i++) {

    if(marray[i]!=0.0) count++;

  }

  fprintf(stderr,"\nConsistency check: %li non-zero entries in marray\n",count);fflush(stderr);

  free(marray);

  t2 = clock();
  fprintf(stdout,"\nGenerating time: %g sec\n",(double)(t2-t1)/CLOCKS_PER_SEC);
  fflush(stdout);

}

double eccentric_anomaly(double t,double T,double e) {

  double tdivT = t/T;
  double u_low = 0.0,u_up = 2.0*3.1415926;
  double u = (u_up + u_low) / 2.0;
  double diff,eps = 0.001;

  do {

    if( (usinu(u,e) - tdivT) < 0.0)
      u_low = u;
    else
      u_up = u;

    u = (u_up + u_low) / 2.0;
    diff = usinu(u,e) - tdivT;
    if(diff < 0.0) diff = -1.0 * diff;

  } while( diff > eps );

  return u;

}

double min_allowed_mr(double m,bool ordered) {

  double q;

  if(ordered) {
    if(m < 5.0) q = 0.0;
    if(m >= 5.0) q = 0.9;
  } else {
    q = 0.0;
  }
  /*
  if(m <= 1.0) q = 0.0;
  if(m > 1.0 && m < 2.0) q = 0.2*(m-1.0);
  if(m >= 2.0) q = 0.2;
  */

  return q;

}

double true_anomaly(double u,double e) {

  double phi;

  phi = acos( (cos(u) - e) / (1.0-e*cos(u)) );
  if(u > PI && u < (2.0*PI) ) phi = 2*PI - phi;

  return phi;

}

double usinu(double u,double e) {

  return ( ( u - e * sin(u)) / (2.0 * PI) );

}

int make_binary(double m1, double r1[3], double v1[3],
		double m2, double r2[3], double v2[3],
		double a, double e, double phase) {

  double M,mu,k,R;
  double r[3],v[3];
  double phi_dot,R_dot;
  int d;

  double G = 39.48423493; // AU, Msun, yr

  M = m1 + m2;
  mu = m1 * m2 / M;

  k = a*(1-e*e);

  R = k/(1+e*cos(phase));

  phi_dot = sqrt( k * mu * m1 * m2 * G )/( mu * R * R );

  R_dot = k * e * sin(phase) * phi_dot / ( (1+e*cos(phase)) * (1+e*cos(phase)) );

  r[0] = R*cos(phase);
  r[1] = R*sin(phase);
  r[2] = 0;
  v[0] = R_dot*cos(phase)-R*sin(phase)*phi_dot;
  v[1] = R_dot*sin(phase)+R*cos(phase)*phi_dot;
  v[2] = 0;

  for(d=0;d<3;d++){
    r2[d] = r[d]*m1/M;
    r1[d] = -r[d]*m2/M;
    v2[d] = v[d]*m1/M;
    v1[d] = -v[d]*m2/M;
  }

  return 0;

}

int rotate(double r[3],double v[3],double alpha,double beta,double gamma) {
  int d;

  double ra[3],rb[3],va[3],vb[3];

  double Da[3][3],Db[3][3],Dg[3][3];

  /*
  fprintf(stderr,"alpha=%g beta=%g gamma=%g\n",
	  alpha/DEGREE,beta/DEGREE,gamma/DEGREE);
  */

  Da[0][0] = cos(alpha);
  Da[0][1] = -sin(alpha);
  Da[0][2] = 0.;
  Da[1][0] = -Da[0][1];
  Da[1][1] =  Da[0][0];
  Da[1][2] = 0.;
  Da[2][0] = 0.;
  Da[2][1] = 0.;
  Da[2][2] = 1.;

  //  print_matix("Da",3,3,Da);

  Db[0][0] = 1.;
  Db[0][1] = 0.;
  Db[0][2] = 0.;
  Db[1][0] = 0.;
  Db[1][1] = cos(beta);
  Db[1][2] = -sin(beta);
  Db[2][0] = 0;
  Db[2][1] = -Db[1][2];
  Db[2][2] = Db[1][1];

  Dg[0][0] = cos(gamma);
  Dg[0][1] = -sin(gamma);
  Dg[0][2] = 0.;
  Dg[1][0] = -Dg[0][1];
  Dg[1][1] =  Dg[0][0];
  Dg[1][2] = 0.;
  Dg[2][0] = 0.;
  Dg[2][1] = 0.;
  Dg[2][2] = 1.;

  for(d=0;d<3;d++){
    ra[d] = Da[d][0]*r[0]+Da[d][1]*r[1]+Da[d][2]*r[2];
    va[d] = Da[d][0]*v[0]+Da[d][1]*v[1]+Da[d][2]*v[2];
  }

  for(d=0;d<3;d++){
    rb[d] = Db[d][0]*ra[0]+Db[d][1]*ra[1]+Db[d][2]*ra[2];
    vb[d] = Db[d][0]*va[0]+Db[d][1]*va[1]+Db[d][2]*va[2];
  }

  /*  
  for(d=0;d<3;d++){
    r[d] = ra[d];
    v[d] = va[d];
  }
  
  return 0;
  */

  for(d=0;d<3;d++){
    r[d] = Dg[d][0]*rb[0]+Dg[d][1]*rb[1]+Dg[d][2]*rb[2];
    v[d] = Dg[d][0]*vb[0]+Dg[d][1]*vb[1]+Dg[d][2]*vb[2];
  }

  return 0;

}

double apparent_separation(double r1[3],double r2[3]) {

  double s,ls,r[3];

  r[0] = r1[0] - r2[0];
  r[1] = r1[1] - r2[1];

  s = pow( r[0] , 2.0 ) + pow( r[0] , 2.0 );
  s = sqrt(s);

  ls = log10(s);

  return ls;

}

void hammer_aitov_projection(double phi,double theta,double *x,double *y) {

  double b,l;

  l = 0.0;
  b = (PI / 2.0) - theta;

  if(phi <= PI) l = phi;
  if(phi > PI) l = -1.0*PI + (phi - PI);

  *x = 2.0*sqrt(2.0)*cos(b)*sin(l/2.0) / ( sqrt(1.0+cos(b)*cos(l/2.0)) );
  *y = sqrt(2.0) * sin(b) / ( sqrt(1.0+cos(b)*cos(l/2.0)) );

}

#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "Ptr.h"
#include "Small_routines.h"
#include "Settings.h"
#include "Constants.h"
#include "InitFinDistr.h"
#include "get_params.h"

bool Synthesize(char *libname,double rho,double mecl,int t,double beta,double sfr,double meclmin,double rh,double mmin,double mmax,int num_sys_obs,bool use_mk12,bool field,bool clust,bool slidebinfrac,bool init,bool fin,bool SpT[],bool bdf[],bool constrain[],double bin_min[],double bin_max[],int smplpnts[],double binwidth[]) {

  unsigned int i;
  FILE *dat;
  double slope,Ecut,asymptote;
  double *ledf,*sigm,*res,c,s,k_ML,lib_entries;
  double num_stars,step,upper,lower,num_stars_tot,num_bin_in,num_bin_fin,num_sys_fin,num_clust;
  double *fE_in;
  double meclmax,mscs,*nrgspect_in,*nrgspect_fin,kecl;//f_bin
  double width=10.;
  char filename[100];
  bool error = false;

  s=1.34;k_ML=0.0144;
  meclmax=pow(pow(k_ML,s)*sfr/pow(10.,-9.07),1./s);
  if(meclmax>M_ECL_MAX_star) meclmax=M_ECL_MAX_star;

  /* sample initial log-energy distribution function and Sigmoidal */
  ledf = dvector(2*SMPLPNTS_lE);
  sigm = dvector(2*SMPLPNTS_lE);
  nrgspect_in = dvector(2*SMPLPNTS_lE);
  nrgspect_fin = dvector(2*SMPLPNTS_lE);

  /* equidistant x-values (same sampling points for all) */
  for(i=0;i<SMPLPNTS_lE;i++) {

    ledf[2*i] = lE_MIN+BINWIDTH_lE / 2. + i*BINWIDTH_lE;
    sigm[2*i] = ledf[2*i];
    nrgspect_in[2*i] = ledf[2*i];
    nrgspect_fin[2*i] = ledf[2*i];

  }

  /* y-values = 0 */
  for(i=0;i<SMPLPNTS_lE;i++) {

    nrgspect_in[2*i+1] = 0.;
    nrgspect_fin[2*i+1] = 0.;

  }

  /* read iEDF */
  strcpy(filename,"flE_eigenevolved.dat");

  if((dat = fopen(filename,"r")) == NULL) {
    fprintf(stderr,"Reading from %s failed\n",filename);
    exit(1);
  }

  for(i=0;i<SMPLPNTS_lE;i++) {

    rewind(dat);
    ledf[2*i+1] = edf(dat,ledf[2*i]);

  }

  fclose(dat);
  /****************************/

  /* fE_in contains the number of binaries per bin for the number of binaries in the library */
  fE_in=(double *)calloc(SMPLPNTS_lE,sizeof(double));
  error = populate_initial_energy_distribution(libname,&lib_entries,fE_in);
  if(error) return true;

  if(init) {
    fprintf(stderr,"--------------------\n");
    fprintf(stderr,"Synthesizing initial population...\n");
    error = orbital_parameter_distributions( libname , lib_entries , nrgspect_in , nrgspect_fin , fE_in , SpT , bdf , constrain , bin_min , bin_max , smplpnts , binwidth, mmin , mmax , num_sys_obs , false , init , false);
    if(error) return true;
  }

  if(field) {

    mscs = DELTAT * sfr; /* mass of Star Cluster System in Msun */

    /* determine Cluster IMF coeff (integral over ECMF = Mtot) */
    if(beta!=2.0) {

      c = 2.0-beta;
      kecl = mscs * c / ( pow(meclmax,c) - pow(meclmin,c) );

    } else {

      kecl = mscs / ( log(meclmax) - log(meclmin) );

    }

    //    if(fin) {

      /* synthesize galaxy population */
      fprintf(stderr,"Synthesizing Galaxy-wide Population...");

      num_stars_tot = 0.; /* total cumulative number of stars (S+2*B) */
      num_sys_fin = 0.; /* total number of systems (S+B) after dynamical processing */
      num_bin_in = 0.; /* cumulative number of initial binaries */
      num_bin_fin = 0.; /* cumulative number of binaries after dynamical processing */

      mecl = 0.;
      lower = meclmin;
      if(meclmin==5.)
	step=5.;
      else
	step=width;

      c=1.-beta;
      while(mecl<meclmax) {

	/* calculate numb. of clusters per bin */
	upper = lower+step;
	mecl = (upper+lower) / 2.0;

	if(beta!=1.0)
	  num_clust = kecl*(pow(upper,c)-pow(lower,c))/c;
	else
	  num_clust = kecl*(log(upper)-log(lower));

	lower=upper;
	step=width;

	/********************************/

	num_stars = mecl / IMF_AV_MASS; /* number of stars per cluster */
	num_stars_tot += num_stars * num_clust; /* cumulative total number of stars in all clusters */
	num_bin_in += num_clust*num_stars / 2.0; /* cumulative number of initial binaries */

	/* populate initial energy spectrum for binaries with stars */
	for(i=0;i<SMPLPNTS_lE;i++)
	  nrgspect_in[2*i+1] += ledf[2*i+1] * BINWIDTH_lE * ( num_clust * num_stars / 2.0);

	if(use_mk12) rh=0.1*pow(mecl,0.13);
	rho= 3.0 * mecl / ( 8.0 * PI * rh * rh * rh ); /* stellar mass density within rh */

	get_params( rho , &Ecut , &slope , &asymptote , 3 );

	/* Set y-values for logistic function */
	for(i=0;i<SMPLPNTS_lE;i++)
	  sigm[2*i+1] = sigmoid( sigm[2*i] , slope , Ecut , asymptote );

	/* multiply logistic function with initial energy distribution to obtain energy distribution after dynamical processing */
	res = mult( ledf , sigm , SMPLPNTS_lE );

	/* fill resulting (cumulative) energy distribution */
	distribute_nrgs( res , nrgspect_fin , num_stars , &num_sys_fin , &num_bin_fin , num_clust );

	free(res);

      }

      fprintf(stderr,"...done.\n");

      //f_bin = num_bin_fin / num_sys_fin;

      /* construct and output orbital-parameter distributions after dynamical processing and summing up all star clusters */
      error = orbital_parameter_distributions( libname , lib_entries , nrgspect_in , nrgspect_fin , fE_in , SpT , bdf , constrain , bin_min , bin_max , smplpnts , binwidth, mmin , mmax , num_sys_obs , slidebinfrac , false , fin);
      if(error) return true;

      //    }

  }

  if(clust) {

    get_params( rho , &Ecut , &slope , &asymptote , t );

    /* Set y-values for logistic function */
    for(i=0;i<SMPLPNTS_lE;i++)
      sigm[i*2+1] = sigmoid( sigm[i*2] , slope , Ecut , asymptote);

    /* multiply logistic function with initial energy distribution to obtain energy distribution after dynamical processing */
    res = mult( ledf , sigm , SMPLPNTS_lE );

    num_sys_fin = 0.; /* number of systems (S+B) after dynamical processing */
    num_bin_fin = 0.; /* number of binaries (B) after dynamical processing */
    num_clust = 1.; /* number of considered clusters */
    num_stars = mecl / IMF_AV_MASS; /* number of stars (S+2*B) in star cluster */

    for(i=0;i<SMPLPNTS_lE;i++)
      nrgspect_in[2*i+1] += ledf[2*i+1] * BINWIDTH_lE * ( num_clust*num_stars / 2.0 );

    distribute_nrgs( res , nrgspect_fin , num_stars , &num_sys_fin , &num_bin_fin , num_clust );

    //f_bin = num_bin_fin / num_sys_fin; /* total binary fraction */

    /* construct and output orbital-parameter distributions after dynamical processing */
    error = orbital_parameter_distributions( libname , lib_entries , nrgspect_in , nrgspect_fin , fE_in , SpT , bdf , constrain , bin_min , bin_max , smplpnts , binwidth, mmin , mmax , num_sys_obs , slidebinfrac , false , fin);
    if(error) return true;

  }

  /* free memory space allocated for pointers */
  free(fE_in);

  free(nrgspect_fin);

  free(sigm);
  free(ledf);

  return false;

}

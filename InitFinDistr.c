/* generates orbital parameter distributions */

#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "InitFinDistr.h"
#include "Small_routines.h"
#include "Ptr.h"
#include "Settings.h"

double edf(FILE *dat,double x)  {
  double xlow1,xlow2=0.,xhigh1,xhigh2=0.,z1,z2,res;
  long pos1,pos2;
  int dummy=0.;

  if(x>=-6 && x<=6) {

    pos1=ftell(dat);

    do {

      pos2=pos1;
      pos1=ftell(dat);
      dummy=fscanf(dat,"%lf %lf %lf %lf\n",&xlow2,&xhigh2,&z2,&res);

    } while(x > (xlow2+xhigh2) / 2.);

    fseek(dat,pos2,SEEK_SET);
    dummy=fscanf(dat,"%lf %lf %lf %lf\n",&xlow1,&xhigh1,&z1,&res);

    /* interpolate */
    res = z1 + ( z2 - z1 ) * ( x - (xlow1+xhigh1) / 2. ) / ( (xlow2+xhigh2) / 2. - (xlow1+xhigh1) / 2. );

  } else
    res=0.;

  dummy*=2.; /* prevent compiler warning */
  return(res);

}

double sigmoid(double x,double slope,double Ecut,double asymptote) {
  double high,res;//low

  //low=-7.;
  high=7.;

  if(x <= high && x > Ecut)
    res = asymptote / ( 1 + exp( slope * (x-Ecut) ) ) - asymptote / 2.0;
  else
    res = 0.;

  return(res);
}

double *mult(double *x,double *y,unsigned int smplpnts) {

  unsigned int i;
  double *z;

  z = dvector( 2 * smplpnts );

  for(i=0;i<smplpnts;i++) {

    z[i*2] = x[i*2];
    z[i*2+1] = x[i*2+1] * y[i*2+1];

  }

  return(z);

}

void distribute_nrgs(double *res,double *nrgspect_fin,double num_stars,double *num_sys_fin,double *num_bin_fin,double num_clust) {

  unsigned int i;
  double num_bin,num_sys,f_bin;

  /* binary fraction = area below distribution */
  f_bin = 0.0;
  for(i=0;i<SMPLPNTS_lE;i++)
    f_bin += res[2*i+1];
  f_bin *= BINWIDTH_lE;

  num_sys = num_stars / (1.+f_bin); /* number of systems (S+B) follows from number of stars in one cluster and binary fraction (eq. 13 in MK11) */
  num_bin = f_bin * num_sys; /* number of binaries follows from binary fraction and  number of systems */

  num_sys *= num_clust; /* systems in one cluster contribute num_clust times */
  *num_sys_fin += num_sys; /* cumulative number of systems */

  num_bin *= num_clust; /* binaries in one cluster contribute num_clust times */
  *num_bin_fin += num_bin; /* cumulative number of binaries  */

  /* add number of binaries per bin to cumulative spectrum */
  for(i=0;i<SMPLPNTS_lE;i++)
    nrgspect_fin[2*i+1] += res[2*i+1] * BINWIDTH_lE * num_sys;

}

bool populate_initial_energy_distribution(char *libname,double *lib_entries,double *fE_in) {

  FILE *lib;
  unsigned int i;
  double m1,m2,e,a,lP,lE,lL,ls,t,f;
  char tmp[10];
  int dummy;

  *lib_entries = 0.;

  fprintf(stdout,"Reading library...");fflush(stdout);
  if(!(lib=fopen(libname,"r"))) {
    fprintf(stderr,"ERROR: The library you specified does not exist.\n");
    return true;
  }
  dummy=fscanf(lib,"%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s\t%s\t%s\n",tmp,tmp,tmp,tmp,tmp,tmp,tmp,tmp,tmp,tmp); /* library header */
  while(!feof(lib)) {
    dummy=fscanf(lib,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",&m1,&m2,&e,&a,&lP,&lE,&lL,&ls,&t,&f);

    for(i=0;i<SMPLPNTS_lE;i++) {

      if(lE>(lE_MIN+i*BINWIDTH_lE) && lE<=(lE_MIN+(i+1)*BINWIDTH_lE)) {

	fE_in[i] += 1.0;
	(*lib_entries)++;
	break;

      }

    }

  }
  fprintf(stdout,"%g entries...done.\n",(*lib_entries));fflush(stdout);

  fclose(lib);

  dummy*=2.;

  return false;

}

bool orbital_parameter_distributions(char *libname,double lib_entries,double *nrgspect_in,double *nrgspect_fin,double *fE_in,bool SpT[],bool bdf[],bool constrain[],double min[],double max[],int smplpnts[],double binwidth[],double mmin, double mmax,int num_sys_obs,bool slidebinfrac,bool init,bool fin) {

  FILE *lib;
  unsigned int i,j,nrgbin;
  int ST,n_ST=9;//n_BDF=10;
  double *xval,**f[NBDF];/* **fE_ST,**fP_ST,**fq_ST,**fe_ST,**fa_ST */
  double row;//**E_tmp,
  double m1,m2,e,la,lP,lE,lL,q,mu,ls,t,h;
  double num_sys_ST[n_ST],mtmp,lib_size,binfrac[9]={0,0,0,0,0,0,0,0,0},binfrac_alt[9]={0,0,0,0,0,0,0,0,0};
  char tmp[10];
  int dummy;
  double param[NBDF];
  char SpT_names[][9]={"M","K","G","F","A","B","O","c","u"},BDF_names[][12]={"e","q","la","lP","lE","mu","lL","ls","t","phi"},filename[100];
  double num_sing_ST[n_ST],num_bin_ST[n_ST],num_bin_moving[32],num_sys_moving[32]; /* from 0.1Msun to 150Msun (log m=-1 > 2.2) */

  if(init || fin) {

    for(i=0;i<NBDF;i++)
      f[i]=dmatrix(n_ST,smplpnts[i]);

  }

  /*
     scale final number of binaries in each bin to library size and round; this is why:

     fE_in contains the number of binaries per bin for the number of binaries in the library for the initial distribution;

     nrgspect_fin and nrgspect_in contain the *calculated* final and initial numbers of binaries per bin based on embedded cluster masses; has nothing to do with the size of the library, therefore these numbers need to be scaled according to the library size.

     it may occur that nrgspect_fin / nrgspect_in > 1, i.e. more binaries in one bin of the final than of initial distribution, meaning there are hardened binaries!

     LIB_SCALE_FACTOR < 1 ensures there will be sufficient hard binaries from the library that fall into the high energy bin if hardened binaries occur.
  */

  lib_size = 0.0; /* stores the scaled size of the library */
  for(i=0;i<SMPLPNTS_lE;i++) {

    lib_size += fE_in[i]*LIB_SCALE_FACTOR;

    if(nrgspect_in[2*i+1] != 0.) {

      nrgspect_fin[2*i+1] = ( nrgspect_fin[2*i+1] / nrgspect_in[2*i+1] ) * fE_in[i] * LIB_SCALE_FACTOR;

    } else {

      if(init)
	nrgspect_fin[2*i+1] = fE_in[i] * LIB_SCALE_FACTOR;
      //else
	//nrgspect_fin[2*i+1] = 0.;

    }

    nrgspect_fin[2*i+1] = dround(nrgspect_fin[2*i+1],0);

  }

  /* open previously generated library of initial binaries */
  if(!(lib=fopen(libname,"r"))) {
    fprintf(stderr,"ERROR: The library you specified does not exist.\n");
    return true;
  }

  /* zero final number of systems, binaries and singles */
  /* 0 = M-type, 6 = O-type, 7=all SpT, 8=user-defined */
  for(i=0;i<n_ST;i++) {

    num_sys_ST[i] = 0.0;
    num_bin_ST[i] = 0.0;
    num_sing_ST[i] = 0.0;

  }

  if(fin && slidebinfrac) {

    for(i=0;i<32;i++) {

      num_bin_moving[i] = 0.0;
      num_sys_moving[i] = 0.0;

    }

  }

  /* read library and construct orbital-parameter distributions after dynamical processing */
  fprintf(stdout,"--------------------\n");fflush(stdout);
  fprintf(stdout,"Calculating binary-fractions and preparing OPDs...");fflush(stdout);

  row=-1.0; /* row number in library of intial binaries */
  dummy=fscanf(lib,"%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s\t%s\t%s\n",tmp,tmp,tmp,tmp,tmp,tmp,tmp,tmp,tmp,tmp); /* read library header */
  while(!feof(lib)) {

    dummy=fscanf(lib,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",&m1,&m2,&e,&la,&lP,&lE,&lL,&ls,&t,&h);
    row++;

    /* swap masses if necessary */
    if(m2>m1) {
      mtmp=m1;
      m1=m2;
      m2=mtmp;
    }

    if(init || fin) {

      q = m2/m1; /* mass-ratio */
      mu = m1*m2 / (m1+m2); /* reduced-mass */

      param[0] = e;
      param[1] = q;
      param[2] = la;
      param[3] = lP;
      param[4] = lE;
      param[5] = mu;
      param[6] = lL;
      param[7] = ls;
      param[8] = t;
      param[9] = h;

    }

    nrgbin = 1e6; /* stores the number of the energy bin the binary selected from library goes into */

    /* construct resulting orbital-parameter distributions according to the resulting binding energy distribution */
    /* loop until the correct energy bin for the binary selected from the library is found */
    for(i=0;i<SMPLPNTS_lE;i++) {

      if(lE>(lE_MIN+i*BINWIDTH_lE) && lE<=(lE_MIN+(i+1)*BINWIDTH_lE)) {

	if(nrgspect_fin[2*i+1]>0.01) { /* if a binary in spectrum left > leave as binary */

	  nrgspect_fin[2*i+1] -= 1.0; /* reduce number of binaries in the respective final energy bin by 1 */
	  nrgbin = i; /* store bin number for later use */

	  ST=check_spectral_type(m1);

	  /* lib_size < number of binaries in library since multiplied by LIB_SCALE_FACTOR */
	  /* if row in library > scaled lib_size it means that there is a binary in the library although all should already have been distributed > it has to be a hardened binary */
	  //if(row < lib_size) { /* standard binary > add 1 binary */

	    if(SpT[8] && m1 >= mmin && m1 <= mmax) {
	      /* add one system and one binary to the user-defined population */
	      num_sys_ST[8]++;
	      num_bin_ST[8]++;
	    }
	    /* add one system and one binary to the SpT-integrated population */
	    num_sys_ST[7]++;
	    num_bin_ST[7]++;
	    /* add one system and one binary to the SpT-dependent population  */
	    num_sys_ST[ST]++;
	    num_bin_ST[ST]++;

	    if(fin && slidebinfrac) {
		j = 0;
		while(j<32) {
		    if( log10(m1) > ((-1.0+j*0.1)-0.1) && log10(m1) < ((-1.0+j*0.1)+0.1) )
			break;
		    j++;
		}
		num_sys_moving[j]++;
		num_bin_moving[j]++;
	    }

	  /* } else { // hardened binary

	    // one single star with SpT of the primary is now a binary > add a binary
	    num_bin_ST[7]++;
	    num_bin_ST[ST]++;
	    // two stars become one binary > subtract 1 system and 1 single star; number of systems with SpT of the primary remains the same
	    num_sys_ST[7]--;
	    num_sing_ST[7]--;
	    // one single star with SpT becomes secondary of the new binary; one system and one single star with SpT disappear
	    ST=check_spectral_type(m2);
	    num_sys_ST[ST]--;
	    num_sing_ST[ST]--;

	    if(SpT[8] && m1 >= mmin && m1 <= mmax) {
	      num_bin_ST[8]++;
	    }
	    if(SpT[8] && m2 >= mmin && m2 <= mmax) {
	      num_sys_ST[8]--;
	      num_sing_ST[8]--;
	    }

	  } */

	} else { /* if no binary left in final spectrum > binary > 2 single stars */

	  if(row < lib_size) { 

	    if(SpT[8] && m1 >= mmin && m1 <= mmax) {
	      /* add one system and one single star to the user-defined population */
	      num_sys_ST[8]++;
	      num_sing_ST[8]++;
	    }
	    if(SpT[8] && m2 >= mmin && m2 <= mmax) {
	      /* add one system and one single star to the user-defined population */
	      num_sys_ST[8]++;
	      num_sing_ST[8]++;
	    }

	    /* add 2 single stars to SpT integrated spectrum */
	    num_sys_ST[7]+=2.0;
	    num_sing_ST[7]+=2.0;

	    if(fin && slidebinfrac) {
		j = 0;
		while(j<32) {
		    if( log10(m1) > ((-1.0+j*0.1)-0.1) && log10(m1) < ((-1.0+j*0.1)+0.1) )
			break;
		    j++;
		}
		num_sys_moving[j]++;
	    }

	    /* add one single star to SpT dependent spectrum */
	    /* for primary and secondary component each */
	    ST=check_spectral_type(m1);
	    num_sys_ST[ST]++;
	    num_sing_ST[ST]++;

	    ST=check_spectral_type(m2);
	    num_sys_ST[ST]++;
	    num_sing_ST[ST]++;

	    if(fin && slidebinfrac) {
		j = 0;
		while(j<32) {
		    if( log10(m2) > ((-1.0+j*0.1)-0.1) && log10(m2) < ((-1.0+j*0.1)+0.1) )
			break;
		    j++;
		}
		num_sys_moving[j]++;
	    }

	  }

	}

	break; /* when correct energy bin is found, no need to loop through remaining bins */

      }

    }

    /* if(output of OPDs requested) */
    if(init || fin) {

      /* if(is a binary - not ionized) -> put binary into BDFs*/
      if(nrgbin != 1e6) {

	if(passes_constraints(param,constrain,min,max)) {

	  ST=check_spectral_type(m1);

	  for(j=0;j<NBDF;j++) {

	    if(bdf[j]) {
	      if(SpT[8] && m1 >= mmin && m1 <= mmax)
		bin_data(f[j][8],param[j],min[j],binwidth[j],smplpnts[j]);//E_tmp,nrgbin
	      if(SpT[7])
		bin_data(f[j][7],param[j],min[j],binwidth[j],smplpnts[j]);
	      if(SpT[ST])
		bin_data(f[j][ST],param[j],min[j],binwidth[j],smplpnts[j]);
	    }

	  }

	}

      }

    }

  }

  fclose(lib);

  fprintf(stdout,"done.\n");fflush(stdout);

  if(init || fin) {

    /* write orbital-parameter distributions to files */

    fprintf(stdout,"Writing %s OPDs...\n",(init?"initial":"final"));fflush(stdout);

    for(i=0;i<NBDF;i++) {

      if(bdf[i]) { /* if orbital-parameter requested for output */

	/* prepare x-values for writing to file */
	xval = (double *)malloc(smplpnts[i]*sizeof(double));
	for(j=0;j<smplpnts[i];j++)
	  xval[j] = min[i]+((2*j+1)*binwidth[i]) / 2.0;

	for(j=0;j<n_ST;j++) {

	  if(SpT[j]) { /* if additionally spectral type is requested for output > write to file */

	    if(j==8) sprintf(filename,"output/%s_distr_u_%.2lf-%.2lf_%s.dat",BDF_names[i],mmin,mmax,(init?"in":"fin"));
	    else sprintf(filename,"output/%s_distr_%s_%s.dat",BDF_names[i],SpT_names[j],(init?"in":"fin"));
	    binfrac_alt[j] = write_output(filename,xval,f[i][j],num_sys_ST[j],num_sys_obs,binwidth[i],smplpnts[i]);
	    fprintf(stdout,"\t%s %s-OPD\t>\t%s\n",(init?"inital":"final"),BDF_names[i],filename);

	  }

	}

	free(xval);

      }

    }

  }

  if(fin && slidebinfrac) write_moving(num_bin_moving,num_sys_moving,32);

  fprintf(stdout,"done.\n");fflush(stdout);

  fprintf(stdout,"\n\t%s binary fractions\ntype\tintegrated\t(constrained)\n",(init?"initial":"final"));
  for(i=0;i<n_ST;i++) {
    binfrac[8-i] = num_bin_ST[8-i] / num_sys_ST[8-i];

    fprintf(stdout,"%s\t%.1lf%%",SpT_names[8-i],SpT[8-i]?(binfrac[8-i]*100.0):999);
    fprintf(stdout,"\t\t(%.1lf%%)\n",SpT[8-i]?(binfrac_alt[8-i]*100.0):999);

  }

  fprintf(stdout,"--------------------\n");fflush(stdout);

  /* free allocated space for pointers */

  //dmatrix_free(E_tmp,SMPLPNTS_lE);

  if(init || fin) {

    for(i=0;i<NBDF;i++)
      dmatrix_free(f[i],n_ST);

  }

  dummy *= 2.0; /* to prevent 'set but unused' compilation warning */

  return false;

}

int check_spectral_type(double m) {

  if(m > LOWMASS_M && m <= LOWMASS_K) return 0; /* M */

  if(m > LOWMASS_K && m <= LOWMASS_G) return 1; /* K */

  if(m > LOWMASS_G && m <= LOWMASS_F) return 2; /* G */

  if(m > LOWMASS_F && m <= LOWMASS_A) return 3; /* F */

  if(m > LOWMASS_A && m <= LOWMASS_B) return 4; /* A */

  if(m > LOWMASS_B  && m <= LOWMASS_O) return 5; /* B */

  if(m > LOWMASS_O) return 6; /* O */

  return -1;

}

bool passes_constraints(double param[],bool constrain[],double min[],double max[]) {

  unsigned int i;

  for(i=0;i<NBDF;i++)
    if(constrain[i] && ( param[i] > max[i] || param[i] < min[i] ))
      return false;

  return true;

}

void bin_data(double *f,double sort_param,double min,double binwidth,double smplpnts) {//,double **fx_P,unsigned int bin
  unsigned int i;

  for(i=0;i<smplpnts;i++) {

    if(sort_param>min+i*binwidth && sort_param<=min+(i+1)*binwidth) {

      f[i]+=1.0;
      //      fx_P[bin][i]+=1.0;

      break;

    }

  }

}

double write_output(char filename[],double *f_x,double *f_y,double num_sys,int num_sys_obs,double binwidth,unsigned int smplpnts) {

  FILE *dat;
  unsigned int i,num_bin;
  double binfrac;

  dat=fopen(filename,"w");

  num_bin=0;
  for(i=0;i<smplpnts;i++) {
    num_bin+=f_y[i];
  }
  for(i=0;i<smplpnts;i++) {
    /* x-value center-of-bin (e.g. lP,lE,lL,a,p,e,q,mu); normalized binary fraction; binary fraction per binwidth */
    fprintf(dat,"%lf\t%lf\t%lf\n",f_x[i],f_y[i]/num_sys/binwidth,f_y[i]/num_sys*num_sys_obs);
    fflush(dat);
  }

  fclose(dat);

  binfrac=(double)num_bin/num_sys;

  return(binfrac);

}

void write_moving(double num_bin[],double num_sys[],int smplpnts) {

  FILE *dat;
  unsigned int i;
  double mass;

  dat=fopen("output/fb_SpT.dat","w");

  for(i=0;i<smplpnts;i++) {
    mass = -1.0 + i*0.1;
    fprintf(dat,"%lf\t%lf\n",mass,num_bin[i]/num_sys[i]);
  }
  fflush(dat);

  fprintf(stdout,"\tBinary fraction as a continous function of primary mass > output/fb_SpT.dat\n");
  fflush(stdout);

  fclose(dat);

}

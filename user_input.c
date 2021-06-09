/* checks user-provided command line arguments */

#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "Constants.h"
#include "Std_values.h"
#include "Settings.h"

void show_help(bool lib,bool clust,bool field,bool init,bool slidebinfrac,bool fin,bool eigen,bool ordered,bool spt);
void set_std_values(long int *Nlib,char *libname,double *mmin,double *mmax,bool *eigen,bool *ordered,double *rho,double *mecl,int *num_sys_obs,double *t,double *beta,double *sfr,double *meclmin,double *rh,bool *init,bool *fin,double *bin_min,double *bin_max,int *smplpnts,double *binwidth);
bool extract_constraints(char argv[],double *bin_min,double *bin_max,int *smplpnts,double *binwidth);

bool check_user_input(int argc,char *argv[],char *libname,long int *Nlib,double *mmin,double*mmax,int *num_sys_obs,bool *eigen,bool *ordered,double *rho,double *mecl,double *t,double *beta,double *sfr,double *meclmin,double *rh,bool *lib,bool *use_mk12,bool *field,bool *clust,bool *slidebinfrac,bool *init,bool *fin,bool *SpT,bool *bdf,bool *constrain,double *bin_min,double *bin_max,int *smplpnts,double *binwidth) {

  unsigned int i,j;
  double s,k_ML,meclmax;
  bool help=false,spt=false,error=false;
  char *spectral_names[9]={"M","K","G","F","A","B","O","integrated","user"},*bdf_names[10]={"eccentricity","mass-ratio","semi-major-axis","period","energy","reduced mass","angular momentum","apparent separation","orbital-time","phase"};
  char string[100];

  set_std_values(Nlib,libname,mmin,mmax,eigen,ordered,rho,mecl,num_sys_obs,t,beta,sfr,meclmin,rh,init,fin,bin_min,bin_max,smplpnts,binwidth);

  for(i=1;i<argc;i++) {

    if(strncmp(argv[i],"genlib",3)==0) *lib=true;
    if(strncmp(argv[i],"clust",5)==0) *clust=true;
    if(strncmp(argv[i],"field",5)==0) *field=true;
    if(strncmp(argv[i],"help",4)==0) help=true;
    if(strncmp(argv[i],"SpT",3)==0) spt=true;

    if(strncmp(argv[i],"libname=",8)==0) sprintf(libname,"Lib/%s",argv[i]+8);
    if(strncmp(argv[i],"Nlib=",5)==0) *Nlib=atoi(argv[i]+5);
    if(strncmp(argv[i],"mmin=",5)==0) *mmin=atof(argv[i]+5);
    if(strncmp(argv[i],"mmax=",5)==0) *mmax=atof(argv[i]+5);
    if(strncmp(argv[i],"-eigen",6)==0) *eigen=false;
    if(strncmp(argv[i],"-op",3)==0) *ordered=false;

    if(strncmp(argv[i],"scale=",6)==0) *num_sys_obs=atoi(argv[i]+6);
    if(strncmp(argv[i],"+fbspt",6)==0) *slidebinfrac=true;

    if(strncmp(argv[i],"mecl=",5)==0) *mecl=atof(argv[i]+5);
    if(strncmp(argv[i],"-evolve",7)==0) *fin=false;
    if(strncmp(argv[i],"t=",2)==0) *t=atof(argv[i]+2);

    if(strncmp(argv[i],"beta=",5)==0) *beta=atof(argv[i]+5);
    if(strncmp(argv[i],"sfr=",4)==0) *sfr=atof(argv[i]+4);
    if(strncmp(argv[i],"meclmin=",8)==0) *meclmin=atof(argv[i]+8);
    if(strncmp(argv[i],"rh=",3)==0) {
      *rh=atof(argv[i]+3);
      if(*rh!=0.) *use_mk12=false;
    }
    if(strncmp(argv[i],"+init",5)==0) *init=true;
    if(strncmp(argv[i],"-field",6)==0) *fin=false;
    if(strncmp(argv[i],"SpT=",4)==0) {
      for(j=0;j<(int)strlen(argv[i]+4);j++) {
	switch(argv[i][4+j]) {
	case 'u':
	  SpT[8]=true;
	  break;
	case 'c':
	  SpT[7]=true;
	  break;
	case 'O':
	  SpT[6]=true;
	  break;
	case 'B':
	  SpT[5]=true;
	  break;
	case 'A':
	  SpT[4]=true;
	  break;
	case 'F':
	  SpT[3]=true;
	  break;
	case 'G':
	  SpT[2]=true;
	  break;
	case 'K':
	  SpT[1]=true;
	  break;
	case 'M':
	  SpT[0]=true;
	  break;
	default:
	  fprintf(stderr,"\n'%c' is not a valid spectral type!\n\n",argv[i][4+j]);
	  return(true);
	}
      }
    }
    if(strncmp(argv[i],"OPD=",4)==0) {
      for(j=0;j<(int)strlen(argv[i]+4);j++) {
	switch(argv[i][4+j]) {
	case 'f':
	  bdf[9]=true;
	  break;
	case 't':
	  bdf[8]=true;
	  break;
	case 's':
	  bdf[7]=true;
	  break;
	case 'L':
	  bdf[6]=true;
	  break;
	case 'm':
	  bdf[5]=true;
	  break;
	case 'E':
	  bdf[4]=true;
	  break;
	case 'P':
	  bdf[3]=true;
	  break;
	case 'a':
	  bdf[2]=true;
	  break;
	case 'q':
	  bdf[1]=true;
	  break;
	case 'e':
	  bdf[0]=true;
	  break;
	default:
	  fprintf(stderr,"\n'%c' is not a valid OPD!\n\n",argv[i][4+j]);
	  return(true);
	}
      }
    }

    if(strncmp(argv[i],"constrain=",10)==0) {

      switch(argv[i][10]) {
      case 'f':
	constrain[9]=true;
	error = extract_constraints(argv[i]+12,&(bin_min[9]),&(bin_max[9]),&(smplpnts[9]),&(binwidth[9]));
	break;
      case 't':
	constrain[8]=true;
	error = extract_constraints(argv[i]+12,&(bin_min[8]),&(bin_max[8]),&(smplpnts[8]),&(binwidth[8]));
	break;
      case 's':
	constrain[7]=true;
	error = extract_constraints(argv[i]+12,&(bin_min[7]),&(bin_max[7]),&(smplpnts[7]),&(binwidth[7]));
	break;
      case 'L':
	constrain[6]=true;
	error = extract_constraints(argv[i]+12,&(bin_min[6]),&(bin_max[6]),&(smplpnts[6]),&(binwidth[6]));
	break;
      case 'm':
	constrain[5]=true;
	error = extract_constraints(argv[i]+12,&(bin_min[5]),&(bin_max[5]),&(smplpnts[5]),&(binwidth[5]));
	break;
      case 'E':
	constrain[4]=true;
	error = extract_constraints(argv[i]+12,&(bin_min[4]),&(bin_max[4]),&(smplpnts[4]),&(binwidth[4]));
	break;
      case 'P':
	constrain[3]=true;
	error = extract_constraints(argv[i]+12,&(bin_min[3]),&(bin_max[3]),&(smplpnts[3]),&(binwidth[3]));
	break;
      case 'a':
	constrain[2]=true;
	error = extract_constraints(argv[i]+12,&(bin_min[2]),&(bin_max[2]),&(smplpnts[2]),&(binwidth[2]));
	break;
      case 'q':
	constrain[1]=true;
	error = extract_constraints(argv[i]+12,&(bin_min[1]),&(bin_max[1]),&(smplpnts[1]),&(binwidth[1]));
	break;
      case 'e':
	constrain[0]=true;
	error = extract_constraints(argv[i]+12,&(bin_min[0]),&(bin_max[0]),&(smplpnts[0]),&(binwidth[0]));
	break;
      default:
	fprintf(stderr,"\n'%c' is not a valid OPD!\n\n",argv[i][10]);
	return(true);
      }

    }

  }

  if(error)
    return(true);

  if(!(*lib) && !(*clust) && !(*field) && !help)
    help = true;

  i=0;
  if(*lib) i++;
  if(*field) i++;
  if(*clust) i++;
  if(i>1) {

    fprintf(stderr,"Only one task at a time.\nDo not combine 'genlib', 'clust' and 'field'.\n");
    fflush(stderr);
    return(true);

  }

  if(help) {

    show_help(*lib,*clust,*field,*slidebinfrac,*init,*fin,*eigen,*ordered,spt);
    return(true);

  }

  if(*lib) {

    if(*mmin < 0.08) {

      fprintf(stderr,"ERROR: minimum mass below hydrogen-burning mass-limit (mmin<0.08 Msun) -> increase 'mmin'.\n");
      fflush(stderr);
      return(true);

    }

    if(*mmax > 150.0) {

      fprintf(stderr,"ERROR: maximum primary mass exceeds 150 Msun -> decrease 'mmax'.\n");
      fflush(stderr);
      return(true);

    }

    if(*mmax <= *mmin) {

      fprintf(stderr,"ERROR: 'mmax' needs to be larger than 'mmin'\n");
      fflush(stderr);
      return(true);

    }

    fprintf(stdout,"Task: Generate Library.\n");
    fprintf(stdout,"Your input:\n");
    fprintf(stdout,"\tOutput to = %s\n\t# of binaries  = %li\n\tmin. prim. mass / Msun = %g\n\tmax. prim. mass / Msun = %g\n\teigenevolution = %s",
	    libname,*Nlib,*mmin,*mmax,*eigen?("ON"):("OFF"));
    fprintf(stdout,"\n--------------------\n");
    fflush(stdout);

  }

  if(*clust || *field) {

    for(i=0;i<9;i++) {

      if(SpT[8-i]) break;
      if(i==8) {
	SpT[7] = true; // all masses
	//	fprintf(stderr,"ERROR: No spectral-type specified. (SpT=...)\n");
	//fflush(stderr);
	//return(true);
      }

    }

    j=0;
    for(i=0;i<NBDF;i++)
      if(bdf[NBDF-1-i]) j++;

    if(j==0) {
      bdf[7] = true; // apparent separation
      //      fprintf(stderr,"ERROR: No orbital-parameter distribution specified. (OPD=...)\n");
      //fflush(stderr);
      //return(true);
    }

    if(SpT[8] && *mmin ==  STD_MMIN && *mmax == STD_MMAX) {
      fprintf(stderr,"ERROR: When choosing user-defined mass range ('SpT=u') define mmin=... and  mmax=...\n");
      fflush(stderr);
      return(true);
    }

  }

  if(*clust) {

    if(*use_mk12) *rh=0.1*pow(*mecl,0.13);
    *rho=3*(*mecl)/(8*PI*(*rh)*(*rh)*(*rh));

    if(*use_mk12) sprintf(string,"%g pc (using rh=0.1xMecl^0.13 from Marks & Kroupa 2012)",*rh);
    else sprintf(string,"%g pc",*rh);

    fprintf(stdout,"Task: Synthesize Cluster Population\n");
    fprintf(stdout,"Your input:\n");
    fprintf(stdout,"\tLibrary: %s\n",libname);
    fprintf(stdout,"\tMecl = %g Msun\n",*mecl);
    fprintf(stdout,"\trh = %s\n",string);
    fprintf(stdout,"\t--> rho = %g Msun / pc^3\n",*rho);
    fprintf(stdout,"\tEvolve cluster for %.1lf Myr\n",*t);
    fprintf(stdout,"Requested output:\n");
    fprintf(stdout,"\t%s",(*init)?"Initial OPDs are ON.\n":"Initial OPDs are OFF.\n");
    fprintf(stdout,"\tThe ");
    for(i=0;i<NBDF;i++)
      if(bdf[NBDF-i-1]) fprintf(stdout,"%s- ",bdf_names[NBDF-1-i]);
    fprintf(stdout,"OPD(s) is/are output for ");
    for(i=0;i<9;i++)
      if(SpT[8-i]) fprintf(stdout,"%s- ",spectral_names[8-i]);
    fprintf(stdout,"type binaries.\n");
    fprintf(stdout,"\tOutput scaled to %i targets searched for companions.\n",*num_sys_obs);
    if(*slidebinfrac) fprintf(stdout,"Output binary fraction as a continous function of primary mass.\n");
    if(SpT[8]) fprintf(stdout,"\tUser-defined mass range: %.2lf-%.2lf Msun\n",*mmin,*mmax);
    fprintf(stdout,"Constraints:");
    j = 0;
    for(i=0;i<NBDF;i++) {
      if(constrain[i]) {
	fprintf(stdout,"\n\t%s in range [%g:%g] and %i bins",bdf_names[i],bin_min[i],bin_max[i],smplpnts[i]);
	j++;
      }
    }
    if(j == 0) fprintf(stdout,"\n\tNone.");
    fprintf(stdout,"\n--------------------\n");
    fflush(stdout);

  }

  if(*field) {

    /* check if meclmax<meclmin */
    s=1.34;k_ML=0.0144;
    meclmax=pow(pow(k_ML,s)*(*sfr)/pow(10.,-9.07),1./s);
    if(meclmax>M_ECL_MAX_star) meclmax=M_ECL_MAX_star;

    if(meclmax < (*meclmin)) {

      fprintf(stderr,"ERROR: max. cluster mass (%g) < min. cluster mass (%g) --> possible solution: increase SFR or reduce MECLMIN.\n",meclmax,*meclmin);
      fflush(stderr);
      return(true);

    }

    sprintf(string,"%g pc",*rh);

    fprintf(stdout,"Task: Synthesize Field Population\n");
    fprintf(stdout,"Your input:\n");
    fprintf(stdout,"\tLibrary: %s\n",libname);
    fprintf(stdout,"\tECMF slope beta = %g\n",*beta);
    fprintf(stdout,"\tMin. cluster mass / Msun  = %g\n",*meclmin);
    fprintf(stdout,"\tTyp. half-mass radius / pc  = %s\n",*use_mk12?"Mecl-dependent (using rh=0.1xMecl^0.13 from Marks & Kroupa 2012)":string);
    fprintf(stdout,"\tStar formation rate / Msun yr^-1 = %g (-> Max. cluster mass / Msun = %g)\n",*sfr,meclmax);
    fprintf(stdout,"Requested output:\n");
    fprintf(stdout,"\t%s",(*init)?"Initial OPDs are ON.\n":"Initial OPDs are OFF.\n");
    fprintf(stdout,"\t%s",(*fin)?"Field OPDs are ON.\n":"Field OPDs are OFF.\n");
    fprintf(stdout,"\tThe ");
    for(i=0;i<NBDF;i++)
      if(bdf[NBDF-1-i]) fprintf(stdout,"%s- ",bdf_names[NBDF-1-i]);
    fprintf(stdout,"OPD(s) is/are computed for ");
    for(i=0;i<9;i++)
      if(SpT[8-i]) fprintf(stdout,"%s- ",spectral_names[8-i]);
    fprintf(stdout,"type binaries.\n");
    fprintf(stdout,"\tOutput scaled to %i targets searched for companions.\n",*num_sys_obs);
    if(*slidebinfrac) fprintf(stdout,"Output binary fraction as a continous function of primary mass.\n");
    if(SpT[8]) fprintf(stdout,"\tUser-defined mass range: %.2lf-%.2lf Msun\n",*mmin,*mmax);
    fprintf(stdout,"Constraints:");
    j = 0;
    for(i=0;i<NBDF;i++) {
      if(constrain[i]) {
	fprintf(stdout,"\n\t%s in range [%g:%g] and %i bins",bdf_names[i],bin_min[i],bin_max[i],smplpnts[i]);
	j++;
      }
    }
    if(j == 0) fprintf(stdout,"\n\tNone.");
    fprintf(stdout,"\n--------------------\n");
    fflush(stdout);

  }

  return(false);

}

void show_help(bool lib,bool clust,bool field,bool slidebinfrac,bool init,bool fin,bool eigen,bool ordered,bool spt) {

  if(!clust && !field && !lib && !spt) {

    fprintf(stdout,"\tGenerate a library of binaries:\n\t\tsee ./BiPoS genlib help\n");
    fprintf(stdout,"--------------------\n");
    fprintf(stdout,"\tSynthesize star cluster population:\n\t\tsee ./BiPoS clust help\n");
    fprintf(stdout,"\tSynthesize galaxy-wide population:\n\t\tsee ./BiPoS field help\n");
    fprintf(stdout,"--------------------\n");
    fprintf(stdout,"\tShow spectral type ranges:\n\t\tsee ./BiPoS SpT help\n");
    fprintf(stdout,"--------------------\n");
    fprintf(stdout,"First-time users: Generate a library of binaries first!\n");
    fflush(stdout);

  }

  if(lib) {

    fprintf(stdout,"Generating a library of binaries\n");
    fprintf(stdout,"--------------------\n");
    fprintf(stdout,"In order to generate a library of binaries use the following Syntax:\n\n");
    fprintf(stdout,"./BiPoS genlib [libname=x] [Nlib=x] [mmin=x] [mmax=x]\n\n");
    fprintf(stdout,"\tExplanation ('*'=required, STD=standard value)\n");
    fprintf(stdout,"\tlibname=x\treplace x with filename of library in folder Lib/ (STD=%s)\n",STD_LIBNAME);
    fprintf(stdout,"\tNlib=x\t\treplace x with number of binaries to generate (STD=%g)\n",STD_NLIB);
    fprintf(stdout,"\tmmin=x\t\treplace x with minimum primary mass in Library in Msun (STD=%g)\n",STD_MMIN);
    fprintf(stdout,"\tmmax=x\t\treplace x with maximum primary mass in Library in Msun (STD=%g)\n",STD_MMAX);
    fprintf(stdout,"\t-eigen\t\tswitch eigenevolution OFF (STD=%s)\n",eigen?("ON"):("OFF"));
    fprintf(stdout,"\t-op\t\tswitch ordered pairing above 5 Msun OFF; pure random pairing (STD=%s)\n",ordered?("ON"):("OFF"));
    fprintf(stdout,"--------------------\n");
    fprintf(stdout,"Examples (distinguish upper- and lower-case characters!):\n");
    fprintf(stdout,"- Generate library of binaries with primary masses between 2 and 3 Msun:\n");
    fprintf(stdout,"  ./BiPoS genlib mmin=2.0 mmax=3.0\n");
    fprintf(stdout,"- Generate library of 1000 binaries into Lib/MyLib.dat with upper primary mass of 10 Msun:\n");
    fprintf(stdout,"  ./BiPoS genlib libname=MyLib.dat Nlib=1000 mmax=10.0\n");

  }

  if(clust) {

    fprintf(stdout,"Synthesizing a star cluster population\n");
    fprintf(stdout,"--------------------\n");
    fprintf(stdout,"In order to synthesize the binary population of a star cluster use the following Syntax:\n\n");
    fprintf(stdout,"./BiPoS clust [OPD=EPLaqemtihs] [SpT=OBAFGKMc] [+init] [-evolve] [mecl=x rh=x] [constrain=x,min,max,bin] [libname=x] [t=1|3|5]\n\n");
    fprintf(stdout,"\tExplanation ('*'=required, STD=standard value)\n");
    fprintf(stdout,"\tclust(*)\ttells the program to synthesize a star cluster population\n");
    fprintf(stdout,"\tOPD=EPLaqemstf\tspecify which orbital-parameter distributions are to be computed. (STD=s)\n");
    fprintf(stdout,"\t\t\t(E-Energy,P-Period,L-angular momentum,a-semimajor axis,q-Mass-ratio,e-eccentrity,m-reduced mass,\n");
    fprintf(stdout,"\t\t\t(s-apparent separation,t-orbital time,f-orbital phase)\n");
    fprintf(stdout,"\tSpT=OBAFGKMcu\tspecify for which spectral types the OPDs are to be computed. (STD=c)\n");
    fprintf(stdout,"\t\t\t(OBAFGKM-according to Havard sequence,c-all masses,u-user defined mass-range)\n");
    fprintf(stdout,"\tmmin=x\t\tdefines minimum mass if SpT=u is set. (STD=%g)\n",STD_MMIN);
    fprintf(stdout,"\tmmax=x\t\tdefines maximum mass if SpT=u is set. (STD=%g)\n",STD_MMAX);
    fprintf(stdout,"\t+init\t\toutput initial OPDs (STD=%s)\n",init?("ON"):("OFF"));
    fprintf(stdout,"\t-evolve\t\tDON'T output evolved OPDs (STD=%s)\n",fin?("ON"):("OFF"));
    fprintf(stdout,"\t+fbspt\t\toutput binary fraction as a continous function of SpT (STD=%s)\n",slidebinfrac?("ON"):("OFF"));
    fprintf(stdout,"\tmecl=x\t\treplace x with the star cluster mass in Msun (STD=%g)\n",STD_MECL);;
    fprintf(stdout,"\trh=x\t\treplace x with the half-mass radius in pc\n\t\t\tOR set rh=0 to use rh=0.1xMecl^0.13 from Marks & Kroupa 2012 (STD)\n");
    fprintf(stdout,"\tlibname=x\treplace x with filename of library in folder Lib/ (STD=%s)\n",STD_LIBNAME);
    fprintf(stdout,"\tt=1|3|5\t\tspecify the time-span for binary evolution in Myr (1,3 and 5 only, STD=%i)\n",STD_TIME);
    fprintf(stdout,"\tscale=x\t\treplace x with number of targets in survey searched for companions (STD=%i)\n",STD_TARGETS);
    fprintf(stdout,"--------------------\n");
    fprintf(stdout,"To place constraints on orbital-parameters use the following syntax:\n\n");
    fprintf(stdout,"\tconstrain=x,min,max,bin\t\treplace x with orbital-parameter identifier,\n\t\t\t\t\treplace 'min' with the minimum value for this parameter (log values for EPLa),\n\t\t\t\t\treplace 'max' with the maximum value for this parameter (log values for EPLa),\n\t\t\t\t\treplace 'bins' with the number of equal-size bins between min and max\n\n");
    fprintf(stdout,"Use 'constrain=...' repeatedly to place constraints on more than one parameter simultaneously. See examples below.\n");
    fprintf(stdout,"--------------------\n");
    fprintf(stdout,"Examples (distinguish upper- and lower-case characters!):\n");
    fprintf(stdout,"- Synthesize the period- and mass-ratio-OPD for G- and M-type binaries:\n");
    fprintf(stdout,"  ./BiPoS clust OPD=Pq SpT=GM mecl=1000 rh=0.5\n");
    fprintf(stdout,"- Synthesize the semi-major-axis-, mass-ratio- and eccentrity-OPD for F-,K-type and all binaries:\n");
    fprintf(stdout,"  ./BiPoS clust OPD=aqe SpT=FKc mecl=300 rh=0.2\n");
    fprintf(stdout,"- Synthesize the initial energy-OPD for G-type and all binaries from Library 'Lib/MyLib.dat',\n  but not the evolved population, and use standard values for Mecl and rh:\n");
    fprintf(stdout,"  ./BiPoS field OPD=E SpT=Gc +init -evolve libname=MyLib.dat\n");

  }

  if(field) {

    fprintf(stdout,"--------------------\n");
    fprintf(stdout,"BiPoS - The Binary Population Synthesizer\n");
    fprintf(stdout,"--------------------\n");
    fprintf(stdout,"Synthesizing galactic field populations\n");
    fprintf(stdout,"--------------------\n");
    fprintf(stdout,"In order to synthesize the binary population of a galactic field use the following Syntax:\n\n");
    fprintf(stdout,"./BiPoS field [OPD=EPLaqem] [SpT=OBAFGKMc] [+init] [-field] [beta=x sfr=x meclmin=x rh=x] [constrain=x,min,max,bin] [libname=x]\n\n");
    fprintf(stdout,"\tExplanation ('*'=required, STD=standard value)\n");
    fprintf(stdout,"\tfield(*)\ttells the program to synthesize a field population\n");
    fprintf(stdout,"\tOPD=EPLaqemstf\tspecify which orbital-parameter distributions are to be computed. (STD=s)\n");
    fprintf(stdout,"\t\t\t(E-Energy,P-Period,L-angular momentum,a-semimajor axis,q-Mass-ratio,e-eccentrity,m-reduced mass\n");
    fprintf(stdout,"\t\t\t(s-apparent separation,t-orbital time,f-orbital phase)\n");
    fprintf(stdout,"\tSpT=OBAFGKMcu\tspecify for which spectral types the OPDs are to be computed. (STD=c)\n");
    fprintf(stdout,"\t\t\t(OBAFGKM-according to Havard sequence,c-all masses,u-user defined mass-range)\n");
    fprintf(stdout,"\tmmin=x\t\tdefines minimum mass if SpT=u is set. (STD=%g)\n",STD_MMIN);
    fprintf(stdout,"\tmmax=x\t\tdefines maximum mass if SpT=u is set. (STD=%g)\n",STD_MMAX);
    fprintf(stdout,"\t+init\t\toutput initial OPDs (STD=%s)\n",init?("ON"):("OFF"));
    fprintf(stdout,"\t-field\t\tDON'T output field OPDs (STD=%s)\n",fin?("ON"):("OFF"));
    fprintf(stdout,"\t+fbspt\t\toutput binary fraction as a continous function of SpT (STD=%s)\n",slidebinfrac?("ON"):("OFF"));
    fprintf(stdout,"\tbeta=x\t\treplace x with the index (negative slope) of the embedded cluster mass-function (ECMF, STD=%g)\n",STD_BETA);
    fprintf(stdout,"\tsfr=x\t\treplace x with the star formation rate in Msun yr^-1 (STD=%g)\n",STD_SFR);
    fprintf(stdout,"\tmeclmin=x\treplace x with the minimum cluster mass in Msun (STD=%g)\n",STD_MECLMIN);
    fprintf(stdout,"\trh=x\t\treplace x with the typical cluster half-mass radius in pc\n\t\t\tOR set rh=0 to use rh=0.1xMecl^0.13 from Marks & Kroupa 2012 (STD)\n");
    fprintf(stdout,"\tlibname=x\treplace x with filename of library in folder Lib/ (STD=%s)\n",STD_LIBNAME);
    fprintf(stdout,"\tscale=x\t\treplace x with number of targets in survey searched for companions (STD=%i)\n",STD_TARGETS);
    fprintf(stdout,"--------------------\n");
    fprintf(stdout,"To place constraints on orbital-parameters use the following syntax:\n\n");
    fprintf(stdout,"\tconstrain=x,min,max,bin\t\treplace x with orbital-parameter identifier,\n\t\t\t\t\treplace 'min' with the minimum value for this parameter (log values for EPLa),\n\t\t\t\t\treplace 'max' with the maximum value for this parameter (log values for EPLa),\n\t\t\t\t\treplace 'bins' with the number of equal-size bins between min and max\n\n");
    fprintf(stdout,"Use 'constrain=...' repeatedly to place constraints on more than one parameter simultaneously. See examples below.\n");
    fprintf(stdout,"--------------------\n");
    fprintf(stdout,"Examples (distinguish upper- and lower-case characters!):\n");
    fprintf(stdout,"- Synthesize field period- and mass-ratio-OPD for G- and M-type binaries:\n");
    fprintf(stdout,"  ./BiPoS field OPD=Pq SpT=GM beta=2.0 sfr=0.01 meclmin=1.0 rh=0.5\n");
    fprintf(stdout,"- Synthesize field semi-major-axis-, mass-ratio- and eccentrity-OPD for F-,K-type and all binaries:\n");
    fprintf(stdout,"  ./BiPoS field OPD=aqe SpT=FKc beta=2.4 sfr=10.0 meclmin=5.0 rh=0.2\n");
    fprintf(stdout,"- Synthesize initial energy-OPD for G-type and all binaries from Library 'Lib/MyLib.dat',\n  but not for the field and use standard values for beta, sfr, meclmin and rh:\n");
    fprintf(stdout,"  ./BiPoS field OPD=E SpT=Gc +init -field libname=MyLib.dat\n");
    fprintf(stdout,"- Synthesize initial and final field OPD for apparent separations and mass-ratios for A-type stars from Library 'Lib/MyLib.dat';\n  binaries with mass-ratios between 0.15 and 1.0 only; distribute mass-ratios into four equal-size bins;\n  use standard values for beta, sfr, meclmin and rh:\n");
    fprintf(stdout,"  ./BiPoS field OPD=sq SpT=A +init libname=MyLib.dat constrain=q,0.15,1.0,4\n");
    fprintf(stdout,"- Synthesize final field eccentriciy OPD for solar-type stars which have mass-ratios q in [0.6:0.9]\n  and semi-major axes log10(a) in [1.5:3.5]; use standard values for remaining parameters\n");
    fprintf(stdout,"  ./BiPoS field OPD=e SpT=G constrain=q,0.6,0.9,3 constrain=a,1.5,3.5,4\n");

  }

  if(spt) {

    fprintf(stdout,"The code adopts the following spectral type (SpT) mass ranges\n");
    fprintf(stdout,"--------------------\n");
    fprintf(stdout,"SpT    |      O      |    B     |    A    |    F    |    G    |     K   |     M\n");
    fprintf(stdout,"m/Msun | %.1lf- %.1lf | %.1lf-%.1lf | %.1lf-%.1lf | %.1lf-%.1lf | %.1lf-%.1lf | %.1lf-%.1lf | %.2lf-%.1lf\n",
	    LOWMASS_O,MHIGH,LOWMASS_B,LOWMASS_O,LOWMASS_A,LOWMASS_B,LOWMASS_F,LOWMASS_A,LOWMASS_G,LOWMASS_F,LOWMASS_K,LOWMASS_G,LOWMASS_M,LOWMASS_K);
    fprintf(stdout,"--------------------\n");
    fprintf(stdout,"These values can be modified in 'Settings.h'; requires re-compilation via 'make clean' and 'make'.\nA user-defined mass range can be defined via the command-line.\n");

  }

}

void set_std_values(long int *Nlib,char *libname,double *mmin,double *mmax,bool *eigen,bool *ordered,double *rho,double *mecl,int *num_sys_obs,double *t,double *beta,double *sfr,double *meclmin,double *rh,bool *init,bool *fin,double *bin_min,double *bin_max,int *smplpnts,double *binwidth) {

  *Nlib = STD_NLIB;
  sprintf(libname,"Lib/%s",STD_LIBNAME);
  *mmin = STD_MMIN;
  *mmax = STD_MMAX;
  *eigen = STD_EIGEN;
  *ordered = STD_ORDERED;

  *num_sys_obs = STD_TARGETS;
  *mecl = STD_MECL;
  *t = STD_TIME;

  *beta = STD_BETA;
  *meclmin = STD_MECLMIN;
  *sfr = STD_SFR;
  *rh = STD_RH;

  *init = STD_INIT;
  *fin = STD_FIN;

  bin_min[0] = e_MIN;
  bin_min[1] = q_MIN;
  bin_min[2] = la_MIN;
  bin_min[3] = lP_MIN;
  bin_min[4] = lE_MIN;
  bin_min[5] = mu_MIN;
  bin_min[6] = lL_MIN;
  bin_min[7] = ls_MIN;
  bin_min[8] = t_MIN;
  bin_min[9] = phi_MIN;

  bin_max[0] = e_MAX;
  bin_max[1] = q_MAX;
  bin_max[2] = la_MAX;
  bin_max[3] = lP_MAX;
  bin_max[4] = lE_MAX;
  bin_max[5] = mu_MAX;
  bin_max[6] = lL_MAX;
  bin_max[7] = ls_MAX;
  bin_max[8] = t_MAX;
  bin_max[9] = phi_MAX;

  smplpnts[0] = SMPLPNTS_e;
  smplpnts[1] = SMPLPNTS_q;
  smplpnts[2] = SMPLPNTS_la;
  smplpnts[3] = SMPLPNTS_lP;
  smplpnts[4] = SMPLPNTS_lE;
  smplpnts[5] = SMPLPNTS_mu;
  smplpnts[6] = SMPLPNTS_lL;
  smplpnts[7] = SMPLPNTS_ls;
  smplpnts[8] = SMPLPNTS_t;
  smplpnts[9] = SMPLPNTS_phi;

  binwidth[0] = BINWIDTH_e;
  binwidth[1] = BINWIDTH_q;
  binwidth[2] = BINWIDTH_la;
  binwidth[3] = BINWIDTH_lP;
  binwidth[4] = BINWIDTH_lE;
  binwidth[5] = BINWIDTH_mu;
  binwidth[6] = BINWIDTH_lL;
  binwidth[7] = BINWIDTH_ls;
  binwidth[8] = BINWIDTH_t;
  binwidth[9] = BINWIDTH_phi;

}

bool extract_constraints(char argv[],double *bin_min,double *bin_max,int *smplpnts,double *binwidth) {

  int i;
  char *ptr,*arg[3];
  bool error = false;

  i = -1;
  ptr = strtok(argv,",");
  while(ptr != NULL) {
    arg[++i] = ptr;
    ptr = strtok(NULL,",");
  }

  if(i < 2) {
    fprintf(stderr,"ERROR: Provide 3 parameters for 'constrain=...'\n");
    error = true;
  }

  if(atof(arg[0]) > atof(arg[1])) {
    fprintf(stderr,"ERROR: min value > max value in 'constrain=...'\n");
    error = true;
  }

  if(atof(arg[0]) < *bin_min) {
    fprintf(stderr,"ERROR: min value too small in 'constrain=...', choose >= %lf\n",*bin_min);
    error = true;
  }

  if(atof(arg[1]) > *bin_max) {
    fprintf(stderr,"ERROR: max value too large in 'constrain=...', choose <= %lf\n",*bin_max);
    error = true;
  }

  if(!error) {
    *bin_min = atof(arg[0]);
    *bin_max = atof(arg[1]);
    *smplpnts = atoi(arg[2]);
    *binwidth = (*bin_max - *bin_min) / *smplpnts;
  }

  return error;

}

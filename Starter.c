/* program starter */

#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "InitFinDistr.h"
#include "Settings.h"

bool check_user_input(int argc,char *argv[],char *libname,long int *Nlib,double *mmin,double *mmax,int *num_sys_obs,bool *eigen,bool *ordered,double *rho,double *mecl,double *t,double *beta,double *sfr,double *meclmin,double *rh,bool *lib,bool *use_mk12,bool *synth,bool *clust,bool *slidebinfrac,bool *init,bool *fin,bool *SpT,bool *bdf,bool *constrain,double *bin_min,double *bin_max,int *smplpnts,double *binwidth);
void Library(FILE *dat,long int Nlib,double mmin,double mmax,bool eigen,bool ordered);
bool Synthesize(char *libname,double rho,double mecl,int t,double beta,double sfr,double meclmin,double rh,double mmin,double mmax,int num_sys_obs,bool use_mk12,bool field,bool clust,bool slidebinfrac,bool init,bool fin,bool SpT[],bool bdf[],bool constrain[],double bin_min[],double bin_max[],int smplpnts[],double binwidth[]);

int main(int argc,char *argv[]) {

  unsigned int i,n_SpT=9;
  long int Nlib;
  int num_sys_obs,smplpnts[NBDF];
  double mmin,mmax,rho,mecl,t,beta,meclmin,rh,sfr,bin_min[NBDF],bin_max[NBDF],binwidth[NBDF];
  bool lib=false,use_mk12=true,field=false,clust=false,slidebinfrac=false,init,fin,error=false,eigen,ordered,SpT[n_SpT],bdf[NBDF],constrain[NBDF];
  char libname[100];
  FILE *dat;

  fprintf(stdout,"\n--------------------\n");
  fprintf(stdout,"BiPoS1 - The Binary Population Synthesizer, version 1\n");
  fprintf(stdout,"by M. Marks @AIfA Bonn, Germany\n");
  fprintf(stdout,"Contact: astro.michi@yahoo.com or joerg@sirrah.troja.mff.cuni.cz\n");
  fprintf(stdout,"--------------------\n");

  for(i=0;i<n_SpT;i++)
    SpT[i] = false;
  for(i=0;i<NBDF;i++) {
    bdf[i] = false;
    constrain[i] = false;
  } 

  /* check command line arguments */
  error=check_user_input(argc,argv,libname,&Nlib,&mmin,&mmax,&num_sys_obs,&eigen,&ordered,&rho,&mecl,&t,&beta,&sfr,&meclmin,&rh,&lib,&use_mk12,&field,&clust,&slidebinfrac,&init,&fin,SpT,bdf,constrain,bin_min,bin_max,smplpnts,binwidth);
  /*
  for(i=0;i<NBDF;i++) {
    fprintf(stderr,"%s\t",constrain[i]?"TRUE":"FALSE");
    fprintf(stderr,"%lf\t",bin_min[i]);
    fprintf(stderr,"%lf\t",bin_max[i]);
    fprintf(stderr,"%i\t",smplpnts[i]);
    fprintf(stderr,"%lf\n",binwidth[i]);
  }
  */
  if(error) {
    fprintf(stdout,"--------------------\n\n");
    return(-1);
  }

  if(lib) {
    dat=fopen(libname,"w");
    Library(dat,Nlib,mmin,mmax,eigen,ordered);
    fclose(dat);
  }

  if(clust || field)
    error = Synthesize(libname,rho,mecl,t,beta,sfr,meclmin,rh,mmin,mmax,num_sys_obs,use_mk12,field,clust,slidebinfrac,init,fin,SpT,bdf,constrain,bin_min,bin_max,smplpnts,binwidth);

  fprintf(stdout,"--------------------\n\n");

  return 0;

}

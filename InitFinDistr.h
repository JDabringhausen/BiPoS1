#include <stdbool.h>

double edf(FILE *dat,double x);
double sigmoid(double x,double slope,double Ecut,double asymptote);
double *mult(double *x,double *y,unsigned int smplpnts);
void distribute_nrgs(double *res,double *nrgspect_fin,double num_stars,double *num_sys_fin,double *num_bin_tot,double num_clust_bin);
bool populate_initial_energy_distribution(char *libname,double *lib_entries,double *fE_in);
bool orbital_parameter_distributions(char *libname,double lib_entries,double *nrgspect_in,double *nrgspect_fin,double *fE_in,bool types[],bool bdf[],bool constrain[],double min[],double max[],int smplpnts[],double binwidth[],double mmin, double mmax,int num_sys_obs,bool slidebinfrac,bool init,bool fin);
double write_output(char filename[],double *f_x,double *f_y,double num_sys,int num_sys_obs,double binwidth,unsigned int smplpnts);
void write_moving(double num_bin[],double num_sys[],int smplpnts);
bool passes_constraints(double param[],bool constrain[],double min[],double max[]);
void bin_data(double *f,double sort_param,double min,double binwidth,double smplpnts);//double **E_f,unsigned int nrgbin
int check_spectral_type(double m);
void write_peqe_distribution(double *fE,double *fP,double *fq,double *fq1,double *fq2,double *fq_G,double *fe,double *fe1,double *fe2,double num_sys_tot_2,double num_bin_tot_2);

#define LOWMASS_M 0.079
#define LOWMASS_K 0.45
#define LOWMASS_G 0.8
#define LOWMASS_F 1.04
#define LOWMASS_A 1.4
#define LOWMASS_B 2.1
#define LOWMASS_O 16.0

#define NBDF 10

#define SMPLPNTS_lE 40
#define SMPLPNTS_lP 12
#define SMPLPNTS_q 10
#define SMPLPNTS_e 10
#define SMPLPNTS_la 30
#define SMPLPNTS_mu 20
#define SMPLPNTS_lL 10
#define SMPLPNTS_ls 30
#define SMPLPNTS_t 10
#define SMPLPNTS_phi 36

#define lE_MIN -10.0
#define lE_MAX 10.0
#define lP_MIN -1.0
#define lP_MAX 11.0
#define q_MIN 0.0
#define q_MAX 1.0
#define e_MIN 0.0
#define e_MAX 1.0
#define la_MIN -5.0
#define la_MAX 5.0
#define mu_MIN 0.0
#define mu_MAX 1.0
#define lL_MIN 15.0
#define lL_MAX 25.0
#define ls_MIN -5.0
#define ls_MAX 5.0
#define t_MIN 0.0
#define t_MAX 1.0
#define phi_MIN -1.0
#define phi_MAX 8.0

#define BINWIDTH_lE ((lE_MAX-lE_MIN)/SMPLPNTS_lE)
#define BINWIDTH_lP ((lP_MAX-lP_MIN)/SMPLPNTS_lP)
#define BINWIDTH_q ((q_MAX-q_MIN)/SMPLPNTS_q)
#define BINWIDTH_e ((e_MAX-e_MIN)/SMPLPNTS_e)
#define BINWIDTH_la ((la_MAX-la_MIN)/SMPLPNTS_la)
#define BINWIDTH_mu ((mu_MAX-mu_MIN)/SMPLPNTS_mu)
#define BINWIDTH_lL ((lL_MAX-lL_MIN)/SMPLPNTS_lL)
#define BINWIDTH_ls ((ls_MAX-ls_MIN)/SMPLPNTS_ls)
#define BINWIDTH_t ((t_MAX-t_MIN)/SMPLPNTS_t)
#define BINWIDTH_phi ((phi_MAX-phi_MIN)/SMPLPNTS_phi)

#define BIN_LIB_NAME "Lib/Lib.dat"
#define LIB_SCALE_FACTOR 0.9

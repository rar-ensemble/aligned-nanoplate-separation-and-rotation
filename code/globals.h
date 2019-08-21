/**
 * File              : globals.h
 * Author            : Christian Tabedzki <tabedzki@seas.upenn.edu>
 * Date              : 13.08.2018
 * Last Modified Date: 15.08.2018
 * Last Modified By  : Christian Tabedzki <tabedzki@seas.upenn.edu>
 */
#ifndef GLOBALS_H
#define GLOBALS_H

#include <complex>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>


#ifdef PAR
#include "mpi.h"
#include "fftw3-mpi.h"
#else
#include "fftw3.h"
#endif

#define min(A,B) ((A)<(B) ? (A) : (B) ) ///< Determine the min of two options
#define max(A,B) ((A)>(B) ? (A) : (B) ) ///< Determine the max of two options

#define PI 3.14159265358979323846263383
using namespace std;

#define Dim 2 ///< Dimension of my system
#define SHOW(a) std::cout << #a << " = " << (a) << std::endl

// Global variables
#ifndef MAIN
extern
#endif
double dx[Dim], ///< The resolution of the system in each of the Dimensions
       L[Dim], ///< The length of the box in each dimension
       V, Vf, V_1_fld_np, V_1_exp_np,
       a_smear, a_squared, rho0, 
       lam_real,  ///< Step size for the negative fields
       lam_imag, ///< Step size in the positve fields
       kappa, 
       chi_ab, ///< Can be positive or negative. multiplied by the length of the diblock chain
       chi_ac, ///< Can be positive or negative. Which chain length is multiplied is still TBD
       chi_bc, ///< Can be positive or negative. Which chain length is multiplied is still TBD 
       fD, Rg2, nD, nAH, phiH, nP, nFP, phiD, phiG,
       n_samples, ///< Number of times the state has been sampled
       ic_pre[2], ///< The prefactor of the initial condoition
       ic_dir[2], ///< The direction of the initial conditions
       ic_period[2], ///< The period for sine and cosine intial conditions
       top_wall_lamA, top_wall_lamB, bot_wall_lamA, bot_wall_lamB,
       wallT, wallXi, exp_nr_c[2][Dim],
       exp_nr_u[2][Dim], L_nr, R_nr, xi_nr,
       *theta, ///< used for integrating field based nanorods
       *phi, ///< used for integrating field based nanorods
       *theta_weights, ///< used for integrating field based nanorods
       *phi_weights, ///< used for integrating field based nanorods
       np_frac,
       error_tol, 
       sigma, ///< The units of sigma are chains/b^2
       ng_per_np;

#ifndef MAIN
extern
#endif
FILE *dot;

#ifndef MAIN
extern
#endif
complex<double> I, *wpl, ///< The w<SUB>+</SUB> field
      	 	*wa, 
		*wb,
	       	*wp, 
		*wabp, ///< The field of W<SUB>AB</SUB><SUP>(+)</SUP>
	       	*wabm, ///< The field of W<SUB>AB</SUB><SUP>(-)</SUP>
		*smwa, ///< 
		*smwb, ///<
		smwp_min,
                *wacp, ///< The field of W<SUB>AC</SUB><SUP>(+)</SUP>
		*wacm, ///< The field of W<SUB>AC</SUB><SUP>(-)</SUP>
		*smwg, ///< The smeareed field of grafted chain
		*wg, ///< The field of grafted chains
		*wbcp, 
		*wbcm,
                **qd, 
		**qddag, 
		**qha, 
		**qg, 
		**qgdag, 
		**qgdag_exp,
                Qd, 
		*rhoha, ///< Density vector for the homopolymer A
		*rhoda, ///< Density vector for the A part of the diblock
		*rhodb, ///< Density vector for the B part of the diblock 
		*rhog, ///< Density vector for the field based grafted chains
		*rhog_exp, ///< Density vector for the explicit nanoparticles chains
		Qha, 
		Qp,
                Qga_exp, *tmp, *tmp2, 
		*gd, 
		*gaa, ///< The Debye function for the A part of the diblock
	       	*gab, ///< The Debye function for the AB part of the diblock
		*gbb, ///< The Debye function for the B part of the diblock
                *gc, ///< The Debye function for the explicit grafted chains of the system
                *etap, *etam,
                *avg_rhoda, *avg_rhodb, *avg_rhoha, *avg_rhog,
                *avg_rhog_exp, *avg_grafts, *avg_expl_grafts,
                *avg_rho_fld_np, *avg_rho_fld_np_c,
                *rho_surf, *surfH, *rho_exp_nr, *exp_nrH,
                *rho_fld_np_c, *rho_fld_np, *fld_npH,
                Hcur, *poly_bond_fft,
                shift_wp , ///< Not referenced in the code as is.
	       	*hhat, **tmp_sph, ***Gamma_aniso, *Gamma_iso,
                ***smwp_aniso, *smwp_iso, ***tmp_aniso, *tmp_iso,
                *exp_neg_smwp_iso, 
		***exp_neg_smwp_aniso,
	       	*expl_grafts, *grafts;

#ifndef MAIN
extern
#endif
int Nx[Dim], ///< Length of the box in grid points
    NxT[Dim], ///< Not referenced in the code as is.
    M, do_CL, ///< Flag for whether to do Complex Langavan 
    iter, ///< Determine which iteration the system is on.
    print_freq, itermax, sample_freq, sample_wait, update_scheme,
    ML, NxL[Dim], zstart, size, 
    myrank, ///< Used to determine the head processor to insure output is only included once.
    nprocs, ///< Number of processors in this job
    N, Ng, Nda, Ndb, Nah, 
    negative_chi_ab_flag, negative_chi_bc_flag, negative_chi_ac_flag,
    ic_flag[3], ///< The flag to determine what type of initial condition is selected.
    do_film, np_type, np_chem, n_exp_nr, ///< The number of explicit nanoparticles 
    do_fld_np, Nu;

#ifndef MAIN
extern
#endif
long long total_alloced;

#ifdef PAR
#ifndef MAIN
extern
#endif
MPI_Status status;

#endif

#ifndef MAIN
extern
#endif
#ifdef PAR
fftw_complex *fmin0, *fmot0 ;
#else
fftw_complex *fin, *fout;
#endif

#ifndef MAIN
extern
#endif
fftw_plan fwd0, fwd1, fbk0, fbk1;


complex<double> integ_trapPBC(complex<double>*);
double gasdev2(void);
double ran2(void);
double get_k( int , double* );
double get_k_alias( int , double* );
int stack(int*);
int stack2(int);
int stack3(int, int);
void unstack(int, int*);
void unstack_local(int, int*);
int unstack_stack(int);
void fft_fwd_wrapper(complex<double>*, complex<double>*);
void fft_bck_wrapper(complex<double>*, complex<double>*);

void write_avg_rho(char* nm, complex<double>*);
void write_avg_dat(char* nm, complex<double>*);


double get_k(int, double*, int);

void calc_gd( complex<double>* , double );
complex<double> calc_wpart(void);

void initialize_1(void);
void initialize_2(void);
void allocate(void);
void calc_poly_density( void );

void zero_average(complex<double>*, int);
void read_input(void);
double pbc_mdr2( double*, double*, double*) ;
void write_data_bin(const char* , complex<double>* ) ;
void field_gradient( complex<double>* , complex<double>* , int );
complex<double> calc_H( void ) ;
void write_outputs( void ) ;

#endif // GLOBALS_H

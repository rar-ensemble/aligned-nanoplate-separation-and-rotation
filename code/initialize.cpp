/**
 * File              : initialize.cpp
 * Author            : Christian Tabedzki <tabedzki@seas.upenn.edu>
 * Date              : 14.08.2018
 * Last Modified Date: 14.08.2018
 * Last Modified By  : Christian Tabedzki <tabedzki@seas.upenn.edu>
 */
/**
 * @brief Initializes the variables and the fields of the System
 * 
 * @file initialize.cpp
 * @author Robert Riggleman
 * @author Jason Koski
 * @author Ben Lindsay
 * @author Christian Tabedzki
 * @date 2018-08-03
 */


#include "globals.h"
#include "cavity.hpp"
void read_resume_files(void);
void sphere_init(void);
complex<double> integ_sphere(complex<double> **);
void calc_gaa(complex<double> *, double);
void calc_gbb(complex<double> *, double);
void calc_gab(complex<double> *, double);
double get_r(int, double *);
double dot_prod(double[Dim], double[Dim]);
double cross_prod(double[Dim], double[Dim]);
void write_data(char *, complex<double> *); // Added for debugging
void init_fields(void);
void init_Gamma_sphere(void);
void init_Gamma_rod(void);
void explicit_nanorod(double, double, double, double[Dim], double[Dim], complex<double> *);
void explicit_nanosphere(double, double, double[Dim]);

 /// Initializes the Rg2, N, fD as well as M, V, and dx.

/**
 * This function has no parameters. 
 * 
 */

void initialize_1()
{

  if (myrank == 0){
    printf("\n---Initialization 1---\n");
    cout << endl;
    }
  I = complex<double>(0.0, 1.0);
  N = Nda + Ndb;
  fD = (N > 0 ? double(Nda) / double(N) : 0.0);
  Rg2 = double(N - 1) / 6.0;
  a_squared = a_smear * a_smear;


  /// Initializes M, V, dx
  M = 1;
  V = 1.0;
  for (int i = 0; i < Dim; i++)
  {
    M *= Nx[i];
    V *= L[i];
    dx[i] = L[i] / double(Nx[i]);
  }

  if (myrank == 0)
  {
    printf("V = %lf\n", V);
    printf("N = %d\n", N);
    printf("fD = %lf\n", fD);
    printf("a^2 = %lf\n", a_squared);
    printf("---Initialization 1 complete---\n\n");
    fflush(stdout);
  }

  n_samples = 0.0;
  total_alloced = 0;
}

/// Initializes fields and other variables


void initialize_2() 
{
  int i, j, k;
  double k2, kv[Dim];
  double mdr2, mdr, r1[Dim], r2[Dim], dr[Dim], x[Dim];

  if (myrank == 0){
    printf("---Initialization 2---");
    cout << endl;
    }

  // Initialize hhat
  for (i = 0; i < ML; i++) 
  {
    k2 = get_k(i, kv);
    hhat[i] = exp(-k2 * a_squared / 2.0);
  }

  // Zero all densities //
  for (i = 0; i < ML; i++)
  {
    rhoda[i] = rhodb[i] = rhoha[i] = rhog_exp[i] = rho_fld_np_c[i] = rho_fld_np[i] = 0.0;
    rho_surf[i] = surfH[i] = rho_exp_nr[i] = exp_nrH[i] = 0.0;
    grafts[i] = expl_grafts[i] = rhog[i] = 0.0;
  }

  // Initialize confining wall //
  if (do_film)
  {
    for (i = 0; i < ML; i++)
    {

      // Distance from current position to z=0 plane
      r1[Dim - 1] = 0.0;
      get_r(i, x);
      for (int j = 0; j < Dim - 1; j++)
        r1[j] = x[j];
      mdr2 = pbc_mdr2(x, r1, dr);
      mdr = sqrt(mdr2);

      rho_surf[i] = 0.5 * (1.0 - erf((fabs(mdr) - wallT) / wallXi));
    }
    write_data_bin("rho_surf", rho_surf);
    fft_fwd_wrapper(rho_surf, surfH);
  }

  // Initialize channel wall if it exists
  extern Cavity *channel;
  if (channel != NULL)
  {
    channel->init_rho();
  }

  // Initialize the densities of multiple particles
  complex<double> **particle_density          = (complex<double> **) fftw_malloc(n_exp_nr * sizeof(complex<double> *));
  complex<double> **particle_density_grad_mag = (complex<double> **) fftw_malloc(n_exp_nr * sizeof(complex<double> *));
  for (i = 0; i < n_exp_nr; i++)
  {
    particle_density[i]          = (complex<double> *) fftw_malloc(ML * sizeof(complex<double>));
    particle_density_grad_mag[i] = (complex<double> *) fftw_malloc(ML * sizeof(complex<double>));
  } 
  
  // Initializing the storage of the gradients
  complex<double> ***density_gradient = (complex<double> ***) fftw_malloc (n_exp_nr * sizeof(complex<double> **));
  for( i = 0; i < n_exp_nr; i++) {
    density_gradient[i]               = (complex<double> **) fftw_malloc (Dim * sizeof(complex<double> *));
    for(j = 0;j < Dim;j++) {
      density_gradient[i][j]          = (complex<double> *) fftw_malloc (ML* sizeof(complex<double> ));
    }
  }

  complex<double> temp_density;

  // Loop over nanorods and add nanorod density to rho_exp_nr
  for (i = 0; i < n_exp_nr; i++)
  {
    if (np_type == 1)
      explicit_nanosphere(R_nr, xi_nr, exp_nr_c[i]);
    else if (np_type == 2)
    {
      
      explicit_nanorod(L_nr*1.05, R_nr*1.05, xi_nr, exp_nr_c[i], exp_nr_u[i], particle_density[i]);
      for (j = 0; j < Dim; j++)
      {
        

          // fft_fwd_wrapper( particle_density[i], density_gradient[i][j] );
        field_gradient(particle_density[i], density_gradient[i][j], j);
      }
      for (j = 0; j < ML; j++)
      {
        temp_density = complex<double> (0.0,0.0);
        for (int k = 0; k < Dim; k++)
        { 
          // This is out of order to because of the way things were allocated earlier [k][j]
          temp_density += density_gradient[i][k][j] * density_gradient[i][k][j];
        }
        particle_density_grad_mag[i][j] = sqrt(temp_density);
        // rho_exp_nr[j] += particle_density[i][j];;
        expl_grafts[j] += particle_density_grad_mag[i][j];
      }
      if (n_exp_nr == (i + 1))
      {
        complex<double> total_int_grafts = integ_trapPBC(expl_grafts);
        for (j = 0; j < ML; j++)
          expl_grafts[j] /= total_int_grafts;
      }
    }
    else
    {
      printf("np_type must be 1 (spheres) or 2 (rods) if n_exp_nr>0 and/or"
             " do_fld_np is set to 1 (true). Exiting...\n");
      exit(1);
    }
  }
  // If there are nanorods (and/or nanospheres) write data and take FT
  if (n_exp_nr>0) {
    write_data_bin("rho_exp_nr",  rho_exp_nr);
    write_data_bin("expl_grafts", expl_grafts);
    fft_fwd_wrapper(rho_exp_nr, exp_nrH);
  }


  // Loop over nanorods and add nanorod density to rho_exp_nr
  for (i = 0; i < n_exp_nr; i++)
  {
    if (np_type == 2)
    {
      
      explicit_nanorod(L_nr, R_nr, xi_nr, exp_nr_c[i], exp_nr_u[i], particle_density[i]);

      for (j = 0; j < ML; j++)
      {
        // temp_density = complex<double> (0.0,0.0);
        // for (int k = 0; k < Dim; k++)
        // { 
        //   // This is out of order to because of the way things were allocated earlier [k][j]
        //   temp_density += density_gradient[i][k][j] * density_gradient[i][k][j];
        // }
        rho_exp_nr[j] += particle_density[i][j];;
      }
    }
    else
    {
      printf("np_type must be 1 (spheres) or 2 (rods) if n_exp_nr>0 and/or"
             " do_fld_np is set to 1 (true). Exiting...\n");
      exit(1);
    }
  }
  // If there are nanorods (and/or nanospheres) write data and take FT
  if (n_exp_nr>0) {
    write_data_bin("rho_exp_nr",  rho_exp_nr);
    fft_fwd_wrapper(rho_exp_nr, exp_nrH);
  }

  // Initialize fields
  init_fields();
  if (myrank == 0) printf("Initialized fields\n");

  ////////////////////////
  // Define the volumes //
  ////////////////////////

  double V_nps, V_exp_nps, V_fld_nps, V_poly;
  // Total volume of all explicit nanoparticles
  V_exp_nps = real(integ_trapPBC(rho_exp_nr));
  // "Free" Volume (everything excluding walls)
  Vf = V - real(integ_trapPBC(rho_surf));
  if (channel != NULL) Vf -= real(integ_trapPBC(channel->rho));
  if (do_fld_np) {
    // Error-check np_frac
    if (np_frac < 0.0) {
      if (myrank == 0) printf("Negative np_frac is not allowed\n");
      exit(1);
    }
    else if (np_frac == 0.0 && myrank == 0)
      printf("WARNING: You turned on field-based nanoparticles but set "
             "np_frac to 0. That's just plain silly.\n");
    // Total volume of all explicit and field-based nanoparticles
    V_nps = np_frac * Vf;
    // Total volume of all field-based nanoparticles (total np vol minus
    // explicit np vol)
    V_fld_nps = V_nps - V_exp_nps;
  }
  else {
    // If not doing field-based nps, total volume of all nanoparticles is just
    // the total volume of all explicit nanoparticles
    V_nps = V_exp_nps;
    // Set np_frac = volume of explicit nanoparticles / free volume
    np_frac = V_nps / Vf;
  }
  // Total volume taken up by polymer chains
  V_poly = Vf * (1.0 - np_frac);
  // Volume of just one explicit nanoparticle if applicable
  if (n_exp_nr > 0) {
    V_1_exp_np = V_exp_nps / double(n_exp_nr);
  }
  else
    V_1_exp_np = 0.0;
  // Initialize Gamma (and V_1_fld_np) if doing field-based nps
  if (do_fld_np) {
    if (np_type == 1) {
      init_Gamma_sphere();
    }
    else if (np_type == 2) {
      init_Gamma_rod();
    }
    else {
      printf("np_type=%d is not an acceptable value. Use 1 for spheres and 2 "
             "for rods.\n", np_type);
      exit(1);
    }
  }

  // Normalize expl_grafts and grafts
  complex<double> exp_norm, fld_norm;
  if (sigma > 0.0) {

    if (do_fld_np) {
      fld_norm = integ_trapPBC(grafts);
    }

    if (n_exp_nr > 0) {
      exp_norm = integ_trapPBC(expl_grafts);
    }

    fft_fwd_wrapper(grafts, grafts);
  }

  for (i=0; i<ML; i++) {
    if (sigma > 0.0) {
      if (do_fld_np) {
        grafts[i] *= V / fld_norm;
      }
      if (n_exp_nr > 0) {
        expl_grafts[i] *= 1.0 / exp_norm;
      }
    }
  }

  // Number of nanoparticles
  if (n_exp_nr == 0 && !do_fld_np)
    nP = 0.0;
  else if (do_fld_np)
    // Define nP based on nP*Vp/V = rho0*np_frac
    nP = rho0 * np_frac * Vf / V_1_fld_np;
  else
    nP = V_nps / V_1_exp_np;

  // Set ng_per_np and find phiG (grafted chain volume fraction)
  if (sigma > 0.0) {
    if (Dim == 2)
      ng_per_np = sigma * 2.0 * PI * R_nr;
    else if (Dim == 3)
      ng_per_np = sigma * 4.0 * PI * R_nr * R_nr;

    phiG = ng_per_np * nP * Ng / rho0 / Vf;
  }
  else {
    ng_per_np = 0.0;
    phiG = 0.0;
  }

  // Number of molecules of diblock and A homopolymer
  phiD = 1.0 - phiH - np_frac - phiG;
  nD = phiD * rho0 * Vf / N; // # of diblock chains
  if (Nah > 0.0)
    nAH = phiH * rho0 * Vf / double(Nah);
  else
    nAH = 0.0;

  if (do_fld_np) {
    // Number of field-based nanoparticles = total - explicit
    nFP = nP - n_exp_nr;
    if (nFP < 0.0) {
      if (myrank == 0)
        printf("The nanoparticle volume fraction you're using is lower than "
               "the volume fraction taken up by the explicit nanoparticles\n");
      exit(1);
    }
  }
  else
    nFP = 0.0;

  if (myrank == 0) {
    printf("Total V_segment actual: %lf\n", nD*N + nAH*Nah + nP*ng_per_np*Ng);
    printf("Total V_segment theoretical: %lf\n", rho0 * V_poly);
    printf("V - Vf: %lf\n", V - Vf);
    printf("rho0 = %lf\n", rho0);
    printf("V_1_exp_np = %lf\n", V_1_exp_np);
    printf("V_1_fld_np = %lf\n", V_1_fld_np);
    printf("nD = %lf\n", nD);
    printf("nP = %lf\n", nP);
    printf("n_exp_nr = %d\n", n_exp_nr);
    printf("nFP = %lf\n", nFP);
    printf("Volume fraction of particles excluding grafts = %lf\n", np_frac);
    printf("Volume fraction of grafts                     = %lf\n", phiG);
    printf("Volume fraction of particles including grafts = %lf\n",
            np_frac + phiG);
    printf("# of grafts per np = %lf\n", ng_per_np);
  }

  // Initialize Debye functions
  calc_gaa(gaa, fD);
  calc_gbb(gbb, fD);
  calc_gab(gab, fD);

  if (Nah > 0) calc_gd( gd, double(Nah-1)/double(N-1) );
  if (sigma > 0) calc_gd( gc, double(Ng-1)/double(N-1) );

  ///////////////////////////////
  // Set up bonding potentials //
  ///////////////////////////////
  if (myrank == 0)
    printf("Setting up Gaussian bonds\n");
  for (i=0; i<ML; i++) {
    k2 = get_k(i, kv);
    poly_bond_fft[i] = exp( -k2 / 6.0 ) ;
  }

  iter = 0;

  if (myrank == 0) {
    printf("---Initialization 2 complete---\n\n");
    fflush(stdout);
  }
} // initialize_2

// Initializes fields based on pattern specified in input file, then
// overwrites with restart file values if present.

// The following chunk of code was removed since we were not sure of the best way to initialize the code.

void init_fields()
{
  int sincos_dir;
  double x[Dim];
  complex<double> *tmp_rho_wall;
  extern Cavity* channel;
  if (channel == NULL) tmp_rho_wall = rho_surf;
  else tmp_rho_wall = channel->rho;
  // // Initialize fields based on custom pattern specified in input file
  for (int i=0; i<ML; i++) {
    wpl[i] = wabp[i] = wabm[i] = 0.0;
    wacp[i] = wacm[i] = 0.0;
    wbcp[i] = wbcm[i] = 0.0;

    get_r(i, x);
    // Sine and cosine pattern
    if (ic_flag[0] < -1) {
      if (real(tmp_rho_wall[i] + rho_exp_nr[i]) > 0.5)
        // -1 if there's some wall or nanorod present
        wabm[i] = 0.0;
      else {
        wabm[i] = 1.0;
        for (int j=0; j<2; j++) { // initial condition flag sets
          sincos_dir = int(ic_dir[j]);
          if (ic_flag[j] == -2) // cosine flag
            wabm[i] *= ic_pre[j]*tanh( cos(2.0 * PI * ic_period[j]
                  * x[sincos_dir] / L[sincos_dir])/0.2 ) 
                  * (chi_ab != 0 ? (negative_chi_ab_flag ? I : 1) : 1);
          else if (ic_flag[j] == -3) // sin flag
            wabm[i] *= ic_pre[j]*tanh( sin(2.0 * PI * ic_period[j]
                  * x[sincos_dir] / L[sincos_dir])/0.2 )
                  * (chi_ab != 0 ? (negative_chi_ab_flag ? I : 1) : 1);
          else if (ic_flag[j] < -4) { // -4 flag is just a factor of 1
            printf("Initial condition flags smaller than -4 not supported\n");
            exit(1);
          }
        } // j (initial condition flag sets)
      }
    }
    // Random fields
    else if (ic_flag[0] == -1) {
            if (chi_ab != 0)
      {
        wabp[i] = ic_pre[0] * (ran2() - 0.5) * (negative_chi_ab_flag ? 1 : I);
        wabm[i] = ic_pre[0] * (ran2() - 0.5) * (negative_chi_ab_flag ? I : 1);
      }

      if (chi_bc != 0)
      {
        wbcp[i] = ic_pre[0] * (ran2() - 0.5) * (negative_chi_bc_flag ? 1 : I);
        wbcm[i] = ic_pre[0] * (ran2() - 0.5) * (negative_chi_bc_flag ? I : 1);
      }

      if (chi_ac != 0)
      {
        wacp[i] = ic_pre[0] * (ran2() - 0.5) * (negative_chi_ac_flag ? 1 : I);
        wacm[i] = ic_pre[0] * (ran2() - 0.5) * (negative_chi_ac_flag ? I : 1);
      }

      // if (all_chi_0)
      {
        wpl[i] = ic_pre[0] * (ran2() - 0.5) * I;
      }
    }
    // Push polymer out of walls
    // if (real(tmp_rho_wall[i] + rho_exp_nr[i]) > 0.5) {
    //   if (chi_ab != 0) {  // This needs to be changed
    //     wabp[i] = -0.5;
    //   }
    //   else {
    //     wpl[i] = -0.5 * I;
    //   }
    }
  
  // bool all_chi_0 = !(bool(chi_ab || chi_bc || chi_ac));

  // for (int i = 0; i < ML; i++)
  // {
  //   wpl[i] = wabp[i] = wabm[i] = 0.0;
    
  //   get_r(i, x);

  //   for 


  //   switch (ic_flag[0])
  //   {
  //   case 0:
  //   if (ic_flag)
  //     wpl[i] =
  //         wabp[i] = wabm[i] =
  //             wbcp[i] = wbcm[i] =
  //                 wacp[i] = wacm[i] = 0.0;
  //     break;

  //   case -1:
  //     if (chi_ab != 0)
  //     {
  //       wabp[i] = ic_pre[0] * (ran2() - 0.5) * (negative_chi_ab_flag ? 1 : I);
  //       wabm[i] = ic_pre[0] * (ran2() - 0.5) * (negative_chi_ab_flag ? I : 1);
  //     }

  //     if (chi_bc != 0)
  //     {
  //       wbcp[i] = ic_pre[0] * (ran2() - 0.5) * (negative_chi_bc_flag ? 1 : I);
  //       wbcm[i] = ic_pre[0] * (ran2() - 0.5) * (negative_chi_bc_flag ? I : 1);
  //     }

  //     if (chi_ac != 0)
  //     {
  //       wacp[i] = ic_pre[0] * (ran2() - 0.5) * (negative_chi_ac_flag ? 1 : I);
  //       wacm[i] = ic_pre[0] * (ran2() - 0.5) * (negative_chi_ac_flag ? I : 1);
  //     }

  //     // if (all_chi_0)
  //     {
  //       wpl[i] = ic_pre[0] * (ran2() - 0.5) * I;
  //     }
  //     break;
  //   case -2:
  //     for (int j = 0; j < 2; j++)
  //     { // initial condition flag sets
  //       sincos_dir = int(ic_dir[j]);
  //       if (ic_flag[j] == -2) // cosine flag
  //         wabp[i] *= ic_pre[j] * tanh(cos(2.0 * PI * ic_period[j] * x[sincos_dir] / L[sincos_dir]) / 0.2);
  //     }
  //       break;
      
  //     case -3:

  //       break;

  //     case -4:
  //       wpl[i] =
  //           wabp[i] = wabm[i] =
  //               wbcp[i] = wbcm[i] =
  //                   wacp[i] = wacm[i] = 1.0;
  //       break;

  //     case default:
  //       cout << "Invalid value for field init parameter. Exiting." << endl;
  //     }
  //   }

    // overwrite fields with restart files if available
    read_resume_files();
  }

  // Initialize 1 explicit nanorod by adding appropriate densities to the
  // rho_exp_nr array
  //
  // INPUTS:
  //    len         Nanorod length
  //    rad         Nanorod radius
  //    xi          Thickness of interface between nanorod and surroundings
  //    rel_center  Position vector for the nanosphere center. Relative to the
  //                  length of each dimension (values from 0 to 1)
  //    u           Nanorod orientation vector in cartesian coordinates. Points
  //                  along long axis of nanorod. Must be unit vector.
  void explicit_nanorod(double len, double rad, double xi,
                        double rel_center[Dim], double u[Dim], complex<double> *nanoparticle_density_func_only)
  {
    int i_global;
    int nn[Dim];
    double u_dot_r, u_cross_r;
    double center[Dim], dr[Dim], x[Dim];

    // Multiply center (values between 0 and 1) by L to get the absolute
    // center point of the particle
    for (int i = 0; i < Dim; i++)
      center[i] = rel_center[i] * L[i];

    for (int i = 0; i < ML; i++)
    {
      // Get nn, the int array showing which point in the each dimension we're at
      i_global = unstack_stack(i);
      unstack(i_global, nn);

      // Get x, the current position
      for (int j = 0; j < Dim; j++)
        x[j] = double(nn[j]) * dx[j];

      // Get dr, the distance vector from nanorod center to current position
      pbc_mdr2(x, center, dr);

      // Compute explicit nanorod density
      u_dot_r = dot_prod(u, dr);
      u_cross_r = cross_prod(u, dr);

      //////////////// CHANGE THE VARIABLE BELOW SINCE IT IS GLOBAL

      nanoparticle_density_func_only[i] = 0.25 * erfc((fabs(u_dot_r) - 0.5 * len) / xi) * erfc((fabs(u_cross_r) - rad) / xi);
    }
} // explicit_nanorod

// Initialize 1 explicit nanosphere by adding appropriate densities to the
// rho_exp_nr array
//
// INPUTS:
//    rad         Nanosphere radius
//    xi          Thickness of interface between nanosphere and surroundings
//    rel_center  Position vector for the nanosphere center. Relative to the
//                  length of each dimension (values from 0 to 1)
void explicit_nanosphere(double rad, double xi, double rel_center[Dim]) {
  int i_global;
  int nn[Dim];
  double dr2, dr_abs;
  double center[Dim], dr[Dim], x[Dim];

  // Multiply rel_center (values between 0 and 1) by L to get the absolute
  // center point of the particle
  for (int i=0; i<Dim; i++)
    center[i] = rel_center[i] * L[i];

  for (int i=0; i<ML; i++) {
    // Get nn, the int array showing which point in the each dimension we're at
    i_global = unstack_stack(i);
    unstack(i_global, nn);

    // Get x, the current position
    for (int j=0; j<Dim; j++)
      x[j] = double(nn[j]) * dx[j];

    // Get dr_abs, the magnitude of the distance from nanorod center to current
    // position
    dr2 = pbc_mdr2(x, center, dr);
    dr_abs = sqrt(dr2);

    // Compute explicit nanosphere density (excluding rho0)
    rho_exp_nr[i] += 0.5 * erfc( (dr_abs-rad)/xi );

    if (sigma > 0.0) {
      double exp_arg = ( dr_abs - R_nr - xi_nr ) / xi_nr ;
      expl_grafts[i] += exp( -exp_arg * exp_arg ) ;
    }

  }
} // explicit_nanosphere

/////////////////////////////////////////////////////////////////////////////
// Initialize Gamma_iso for spherical (isotropic) particles
void init_Gamma_sphere() {
  if (R_nr <= 0.0 && myrank == 0) {
    printf("R_nr=%lf is invalid. Try again\n", R_nr);
    exit(1);
  }

  int nn[Dim];
  double dr[Dim], x[Dim], origin[Dim] = { 0.0 };

  for (int i=0; i<ML; i++) {
    // Get distance from nanorod center to current position (dr)
    int i_global = unstack_stack(i);
    unstack(i_global, nn);
    for (int j=0; j<Dim; j++)
      x[j] = double(nn[j])*dx[j];
    double dr2 = pbc_mdr2(x, origin, dr);
    double dr_abs = sqrt(dr2);

    // Compute nanosphere density (including rho0)
    Gamma_iso[i] = 0.5 * rho0 * erfc( (dr_abs-R_nr)/xi_nr );

    double exp_arg = ( dr_abs - R_nr - xi_nr ) / xi_nr ;

    grafts[i] = exp( -exp_arg * exp_arg ) ;

    // Multiply in factor of V because this will speed up convolution
    // calculations later, which all need a factor of V
    Gamma_iso[i] *= V;
  }

  // Calculate volume of 1 nanosphere for use later. The 1/V factor is to
  // cancel out the extra factor of V multiplied in above.
  V_1_fld_np = real(integ_trapPBC(Gamma_iso)) / V;

  // Fourier transform Gamma_iso and leave it that way. It's only used for
  // convolutions which are all done in k-space anyway.
  fft_fwd_wrapper(Gamma_iso, Gamma_iso);
} // END init_Gamma_sphere

// Initialize Gamma_aniso for nanorods (anisotropic)
void init_Gamma_rod() {
  if (Dim < 3) {
    if (myrank == 0)
      printf("Field-based anisotropic nanoparticles currently only supported"
             " in 3D\n");
    exit(1);
  }

  if (L_nr <= 0.0 || R_nr <= 0.0) {
    if (myrank == 0) {
      printf("L_nr=%lf, R_nr=%lf is invalid. Try again\n", L_nr, R_nr);
    }
    exit(1);
  }

  // Read legendre abscissa and weights once with processor 0
  if (myrank == 0) {
    sphere_init();
  }

#ifdef PAR
  // Broadcast the values to all other processors
  MPI_Bcast(theta,           Nu, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(theta_weights,   Nu, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(phi,           2*Nu, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(phi_weights,   2*Nu, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // Compute nanorod density (including rho0)
  for (int i=0; i<Nu; i++) {
    for (int j=0; j<2*Nu; j++) {
      for (int k=0; k<ML; k++) {
        double r[Dim], u[Dim], u_dot_r, u_cross_r;
        get_r(k, r);
        // Convert from spherical coords to x, y, and z
        u[0] = sin(theta[i]) * cos(phi[j]);
        u[1] = sin(theta[i]) * sin(phi[j]);
        u[2] = cos(theta[i]);
        // Get absolute values of dot and cross products of u and r
        u_dot_r = abs(dot_prod(u, r));
        u_cross_r = abs(cross_prod(u, r));
        // Nanorod density equation
        Gamma_aniso[i][j][k] = 0.25 * rho0
          * erfc( (u_dot_r-0.5*L_nr) / xi_nr )
          * erfc( (u_cross_r-R_nr) / xi_nr );
        // Multiply by V so we don't have to do it more expensively during
        // convolution
        Gamma_aniso[i][j][k] *= V;
      } // k
      tmp_sph[i][j] = integ_trapPBC(Gamma_aniso[i][j]);
      // Fourier transform Gamma_aniso and leave it that way. It's only used
      // for convolutions which are all done in k-space anyway.
      fft_fwd_wrapper(Gamma_aniso[i][j], Gamma_aniso[i][j]);
    } // j
  } // i

  // Calculate average volume of 1 nanorod based on
  // Vp = integral ( dr * 1/(4PI) * integral( du * Gamma ) )
  // The extra V in the denominator is to cancel out the extra V included in
  // Gamma to speed up convolution later
  V_1_fld_np = real( integ_sphere(tmp_sph) ) / (4 * PI * V);

} // init_Gamma_rod

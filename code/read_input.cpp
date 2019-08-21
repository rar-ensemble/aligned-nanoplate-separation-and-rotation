#include "globals.h"
#include "cavity.hpp"
#define BUF_LEN 256
void random_particles( void ) ;
void allocate_particles( void ) ;

void read_input() {

  FILE *inp;
  inp = fopen("bcp.input", "r");

  if (inp==NULL ) {
    if (myrank ==0 )
    printf("Failed to open bcp.input!\nTrying bcp2\n");
    inp = fopen("bcp.input2","r");
    if (inp==NULL) {
      printf("Failed to open bcp.input2!");
      exit(1);
    }

    
  }

  char tt[BUF_LEN];
  double dm1, dm2;
  int di1, di2, i, j;

  if (myrank == 0) printf("---Reading bcp.input---\n");

  // Main Parameters
  fgets(tt, BUF_LEN, inp);
  fscanf(inp, "%d %d", &Nda, &Ndb);    fgets(tt, BUF_LEN, inp);
  // if (Ndb > 0) {
  //   if (myrank == 0) {
  //     printf("Only A segments are allowed!!!!\n");
  //   }
  //   exit(1);
  // }
  fscanf(inp, "%lf %d", &phiH, &Nah);  fgets(tt, BUF_LEN, inp);
  if (phiH > 0) {
    if (myrank == 0) {
      printf("No homopolymer allowed unless it's a diblock homopolymer!!!!\n");
    }
    exit(1);
  }
  fscanf(inp, "%lf", &a_smear);        fgets(tt, BUF_LEN, inp);
  fscanf(inp, "%lf", &rho0);              fgets(tt, BUF_LEN, inp);

  fgets(tt, BUF_LEN, inp);
  fgets(tt, BUF_LEN, inp);

  // Interaction Parameters
  fscanf(inp, "%lf", &chi_ab);   fgets(tt, BUF_LEN, inp);
  if (chi_ab < 0.0) {
    if (myrank == 0) {
      printf("Absolute valuing chi_ab\n");
    }
    negative_chi_ab_flag = 1;
    chi_ab = - chi_ab;
  }
  else {
    negative_chi_ab_flag = 0;
  }
  fscanf(inp, "%lf", &chi_bc);   fgets(tt, BUF_LEN, inp);
  if (chi_bc < 0.0) {
    if (myrank == 0) {
      printf("Absolute valuing chi_bc\n");
    }
    negative_chi_bc_flag = 1;
    chi_bc = - chi_bc;
  }
  else {
    negative_chi_bc_flag = 0;
  }
  fscanf(inp, "%lf", &chi_ac);   fgets(tt, BUF_LEN, inp);
  if (chi_ac < 0.0) {
    if (myrank == 0) {
      printf("Absolute valuing chi_ac\n");
    }
    negative_chi_ac_flag = 1;
    chi_ac = - chi_ac;
  }
  else {
    negative_chi_ac_flag = 0;
  }

  fscanf(inp, "%lf", &kappa); fgets(tt, BUF_LEN, inp);
  if (myrank == 0){
  SHOW(chi_ab);
  SHOW(chi_bc);
  SHOW(chi_ac);

  SHOW(negative_chi_ab_flag);
  SHOW(negative_chi_bc_flag);
  SHOW(negative_chi_ac_flag);
  SHOW(kappa);
  }

  fgets(tt, BUF_LEN, inp);
  fgets(tt, BUF_LEN, inp);

  // Initial condition flags
  fscanf(inp, "%d %lf %lf %lf",
         &ic_flag[0], &ic_pre[0], &ic_dir[0], &ic_period[0]);
  fgets(tt, BUF_LEN, inp);
  fscanf(inp, "%d %lf %lf %lf", &ic_flag[1],
          &ic_pre[1], &ic_dir[1], &ic_period[1]);
  fgets(tt, BUF_LEN, inp);

  fgets(tt, BUF_LEN, inp);
  fgets(tt, BUF_LEN, inp);

  // Numerical parameters
  for (i=0; i<Dim; i++)
    fscanf(inp, "%d", &Nx[i]);
  fgets(tt, BUF_LEN, inp);
  for (i=0; i<Dim; i++)
    fscanf(inp, "%lf", &L[i]);
  fgets(tt, BUF_LEN, inp);
  if (myrank == 0)
    for (i=0; i<Dim; i++)
      printf("Nx%d: %d Lx%d: %lf\n", i, Nx[i], i, L[i]);
  fflush(stdout);

  // Lambda/time step parameters
  fscanf(inp, "%lf", &lam_imag); fgets(tt,BUF_LEN,inp);
  fscanf(inp, "%lf", &lam_real); fgets(tt,BUF_LEN,inp);

  // Iteration number parameters
  fscanf(inp, "%d", &itermax);                       fgets(tt, BUF_LEN, inp);
  fscanf(inp, "%d", &print_freq);                    fgets(tt, BUF_LEN, inp);
  fscanf(inp, "%d %d", &sample_freq, &sample_wait);  fgets(tt, BUF_LEN, inp);
  fscanf(inp, "%le", &error_tol);                    fgets(tt, BUF_LEN, inp);
  fscanf(inp, "%d", &update_scheme);                 fgets(tt, BUF_LEN, inp);
  fscanf(inp, "%d", &do_CL);                         fgets(tt, BUF_LEN, inp);

  fgets(tt, BUF_LEN, inp);
  fgets(tt, BUF_LEN, inp);

  // Film parameters
  fscanf(inp, "%d", &do_film);              fgets(tt, BUF_LEN, inp);
  fscanf(inp, "%lf %lf", &wallT, &wallXi);  fgets(tt, BUF_LEN, inp);
  fscanf(inp, "%lf %lf", &top_wall_lamA, &top_wall_lamB);
  fgets(tt, BUF_LEN, inp);
  fscanf(inp, "%lf %lf", &bot_wall_lamA, &bot_wall_lamB);
  fgets(tt, BUF_LEN, inp);

  fgets(tt, BUF_LEN, inp);
  fgets(tt, BUF_LEN, inp);

  // Channel wall parameters
  int do_channel, hor_dir, vert_dir;
  double wall_width, channel_width, xi;
  extern Cavity *channel;
  fscanf(inp, "%d %d %d", &do_channel, &hor_dir, &vert_dir);
  fgets(tt, BUF_LEN, inp);
  fscanf(inp, "%lf %lf %lf", &wall_width, &channel_width, &xi);
  fgets(tt, BUF_LEN, inp);
  if (do_film && do_channel)
  {
    cout << "Can't do film and channel in same simulation" << endl;
    exit(1);
  }
  else if (do_channel)
  {
    if (hor_dir == vert_dir || hor_dir >= Dim || vert_dir >= Dim ||
        hor_dir < 0 || vert_dir < 0) {
      cout << "vert_dir = " << vert_dir << ", hor_dir = " << hor_dir << endl;
      cout << "That's wrong!" << endl;
      exit(1);
    }
    channel = new Channel(channel_width, wall_width, xi, hor_dir, vert_dir);
  }
  else {
    channel = NULL;
  }

  fgets(tt, BUF_LEN, inp);
  fgets(tt, BUF_LEN, inp);

  // Nanoparticle parameters
  fscanf(inp, "%d", &n_exp_nr);                            fgets(tt, BUF_LEN, inp);
  fscanf(inp, "%d", &do_fld_np);                           fgets(tt, BUF_LEN, inp);
  if (do_fld_np) {
    if (myrank == 0) {
      printf("No field-based nanoparticles!\n");
    }
    exit(1);
  }
  fscanf(inp, "%d", &np_type);                             fgets(tt, BUF_LEN, inp);
  fscanf(inp, "%d", &np_chem);                             fgets(tt, BUF_LEN, inp);
  fscanf(inp, "%lf", &sigma);                              fgets(tt, BUF_LEN, inp);
  // if (sigma > 0) {
  //   if (myrank == 0) {
  //     printf("No grafting!\n");
  //   }
  //   exit(1);
  // }
  fscanf(inp, "%d", &Ng);                                  fgets(tt, BUF_LEN, inp);
  if (Ng == 0){sigma = 0;}
  fscanf(inp, "%d", &Nu);                                  fgets(tt, BUF_LEN, inp);
  fscanf(inp, "%lf", &np_frac);                            fgets(tt, BUF_LEN, inp);
  if (!do_fld_np) np_frac = 0.0;
  fscanf(inp, "%lf %lf %lf", &L_nr, &R_nr, &xi_nr);        fgets(tt, BUF_LEN, inp);
  for (j=0; j<2; j++) {
    for (i=0; i<Dim; i++) {
      fscanf(inp, "%lf", &exp_nr_c[j][i]);
      if ( (n_exp_nr > 0) &&
           (exp_nr_c[j][i] < 0.0 || exp_nr_c[j][i] > 1.0) ) {
        if (myrank == 0) {
          printf("Nanorod center[%d]=%lf isn't valid. Must be between 0 and 1 "
                 "since it's a fraction relative to L[%d]\n", i,
                 exp_nr_c[j][i], i);
        }
        exit(1);
      }
    }
    fgets(tt, BUF_LEN, inp);

    dm1=0;
    for (i=0; i<Dim; i++) {
      fscanf(inp, "%lf", &exp_nr_u[j][i]);
      // Get vector magnitude to normalize as it's scanned in
      dm1 += exp_nr_u[j][i] * exp_nr_u[j][i];
    }
    fgets(tt, BUF_LEN, inp);
    dm2 = sqrt(dm1);
    if (n_exp_nr+do_fld_np > 0 && dm2 == 0 && myrank == 0) {
      printf("Nanorod orientation vector must be nonzero\n");
      exit(1);
    }

    // Normalize nanorod orientation vector
    for (i=0; i<Dim; i++)
      exp_nr_u[j][i] /= dm2;
  }
  fclose(inp);
  if (myrank == 0) {
    printf("---Reading of bcp.input complete---\n");
    fflush(stdout);
  }
}



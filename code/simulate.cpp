#include "globals.h"
#include "r2.h"
void write_outputs(void);
void accumulate_all_averages(void);
void update_1s(void);
void update_Euler(void);
void write_data(char*, complex<double>*); // for debugging
complex<double> wpl_part(void);
complex<double> wab_part(void);
complex<double> wac_part(void) ;
complex<double> wbc_part(void) ;

// This is the routine that is essentially the main routine in a code that
// doesn't use Brent's method. Int calculates the equilibrium structure and
// energy and returns the real part of the hamiltonian.
double simulate() {
  complex<double> Hcur, Ho, H;
  double error;
  int close_to_tol_count = 0;
  FILE *otp;
  otp = fopen("data.dat", "w");
  if (otp == NULL) {
    printf("Failed to open data.dat!\n");
    exit(1);
  }

  if (myrank == 0) printf("------Starting simulation------\n");

  // Initialize variables and fields
  initialize_1();
  initialize_2();

#ifdef PAR
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  calc_poly_density();
  if (do_fld_np) {
    // Field nanoparticle sanity check
    complex<double> np_check = integ_trapPBC(rho_fld_np_c)
                               + double(n_exp_nr);
    complex<double> npVp_check = integ_trapPBC(rho_fld_np)
                               + integ_trapPBC(rho_exp_nr) * rho0;
    complex<double> rho0_check_easy = ( nD*N + nAH*Nah
                                        + ng_per_np * nP * Ng
                                        + nP*V_1_fld_np ) / Vf;
    complex<double> rho0_check_hard = ( nD*N + nAH*Nah
                                        + ng_per_np * nP * Ng
                                        + npVp_check ) / Vf;
    if (myrank==0) {
      printf("Field-based nanoparticle calculations sanity check:\n");
      printf("nP=%lf, np_check=%lf\n", nP, real(np_check));
      printf("nFP = %lf, n_exp_nr = %d, nFP + n_exp_nr = %lf =? nP\n",
              nFP,       n_exp_nr,      nFP + n_exp_nr );
      printf("nP*Vp = %lf, npVp_check=%lf\n",
              nP*V_1_fld_np , real(npVp_check) );
      printf("rho0 = %lf, rho0_check_easy=%lf, rho0_check_hard=%lf\n\n",
              rho0,  real(rho0_check_easy), real(rho0_check_hard) );
    }
  }

  double smrhodadb = real( integ_trapPBC(rhoda)+integ_trapPBC(rhodb) );
  double smrhoha = real( integ_trapPBC(rhoha) );

  if (myrank == 0) {
    printf("Initial densities calculated!\n");
    printf("Segment counts:\n");
    printf("nD * N = %lf integ(rhoda + rhodb) = %lf\n",
           nD * N, smrhodadb);
    printf("nAH * Nah = %lf integ(rhoha) = %lf\n",
           nAH * Nah, smrhoha);
    fflush(stdout);
  }

  Ho = calc_H();

  if (myrank == 0) printf("Starting H: %lf\n\n", real(Ho));

  write_outputs();

  if (myrank == 0) printf("---Entering main loop---\n");

  ///////////////
  // MAIN LOOP //
  ///////////////
  for (iter; iter<=itermax; iter++) {
    if (update_scheme == 0)
      update_Euler();
    else
      update_1s();

    if (do_CL && iter >= sample_wait && iter % sample_freq == 0)
      accumulate_all_averages();

    ////////////
    // OUTPUT //
    ////////////
    if (iter % print_freq == 0) {
      Ho = Hcur;
      H = Hcur = calc_H();
      if (Hcur != Hcur) { // If Hcur is NaN
        printf("Crashed! iteration %d\n", iter);
        printf("H: %lf Qd: %lf\n", real(H), real(Qd));
        write_data_bin((char *)"crashed.wpl", wpl);
        write_data_bin((char *)"crashed.wabp", wabp);
        write_data_bin((char *)"crashed.wabm", wabm);
        write_data_bin((char *)"crashed.wacm", wacm);
        write_data_bin((char *)"crashed.wacp", wacp);
        write_data_bin((char *)"crashed.wbcm", wbcm);
        write_data_bin((char *)"crashed.wbcp", wbcp);
        exit(1);
      }
      error = abs(H - Ho) / V / double(print_freq);
      if (error < error_tol * 10) {
        close_to_tol_count++;
      }
      if (myrank == 0) {
        printf("Iter: %d, H=%lf", iter, real(H) );
        if (do_CL)
          printf(" + i%lf", imag(H) );
        if (nD > 0.0)
          printf(", -log(Qd)=%lf", real(-log(Qd)) );
        if (nAH > 0.0)
          printf(", -log(Qha)=%lf", real(-log(Qha)) );
        if (do_fld_np)
          printf(", -log(Qp)=%lf", real(-log(Qp)+smwp_min) );
	if (sigma > 0)
		printf(", -log(Qg)=%lf", - double(n_exp_nr) * ng_per_np * real(Qga_exp));
        printf(", err=%1.1e\n", error);
        fflush(stdout);
      }
      if (myrank == 0) {
        fprintf(otp, "%d %5.6lf %1.3e %5.6lf %1.3e %5.6lf %1.3e %5.6lf %1.3e %1.3e",
                iter, real(H), imag(H), real(-log(Qd)), imag(-log(Qd)),
                real(-log(Qha)), imag(-log(Qha)),
                real(-log(Qp)+smwp_min), imag(-log(Qp)+smwp_min), error);
        // fprintf(otp, "%d %5.6lf %1.3e %5.6lf %1.3e %5.6lf %1.3e %5.6lf %1.3e %1.3e",
        //         iter, real(H), imag(H), real(wpl_part()+wab_part()++wbc_part()+wac_part()),
        //         imag(wpl_part()+wab_part()+wbc_part()+wac_part()),
        //         real(-nD*log(Qd)), imag(-nD*log(Qd)),
        //         real(n_exp_nr * ng_per_np * Qga_exp),
        //         imag(n_exp_nr * ng_per_np * Qga_exp),
        //         real(-log(Qp)+smwp_min), imag(-log(Qp)+smwp_min), error);
        fprintf(otp, "\n");
        fflush(otp);
      }
      write_outputs();
    } // Output

    if (!do_CL && iter > 25 && error < error_tol && close_to_tol_count > 10) {
      if (myrank == 0) {
        printf("Tolerance reached. Error = %.4e\n", error);
        printf("---Main loop complete---\n\n");
      }
      break;
    }

  }// Main Loop

  // Close output stream
  fclose(otp);

  if (myrank == 0) {
    printf("------Completed simulation------\n\n");
    fflush(stdout);
  }

}

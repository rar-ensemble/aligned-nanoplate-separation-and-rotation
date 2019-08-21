#include "globals.h"
#include "cavity.hpp"
void generate_1s_noise( complex<double>* , double ) ;

void update_Euler( ) {

  int i;
  complex<double> evar , F, numer;

  // printf("Euler method not supported yet!\n");
  // exit(1);
  
  // Here we make all the polymers into fourier transforms of themselves.
  if (nD > 0.0) {
    fft_fwd_wrapper(rhoda, rhoda);
    fft_fwd_wrapper(rhodb, rhodb);
  }

  if (nAH > 0.0)
    fft_fwd_wrapper(rhoha, rhoha);

  if (sigma > 0.0 && nFP > 0.0)
    fft_fwd_wrapper(rhog, rhog);

  if (sigma > 0.0 && n_exp_nr > 0.0)
    fft_fwd_wrapper(rhog_exp, rhog_exp);

 if (do_fld_np)
    fft_fwd_wrapper(rho_fld_np, rho_fld_np);

  fft_fwd_wrapper(wpl, wpl);


  if (chi_ab != 0.0) {
    fft_fwd_wrapper(wabm, wabm);
    fft_fwd_wrapper(wabp, wabp);
  }

  if (chi_bc != 0.0) {
    fft_fwd_wrapper(wbcm, wbcm);
    fft_fwd_wrapper(wbcp, wbcp);
  }
  
  if (chi_ac != 0.0) {
    fft_fwd_wrapper(wacm, wacm);
    fft_fwd_wrapper(wacp, wacp);
  }

  if (do_CL)
    generate_1s_noise(etap, lam_imag);

  for (i=0; i<ML; i++) {
    // Update w+ field //
    if (i == 0 && myrank == 0)
      evar = 1.0 ;
    else
      evar = 0.0 ;
    F = (kappa <= 0.0 ? 0.0 : rho0 * wpl[i] / kappa)
        + I * rho0 * (surfH[i] + exp_nrH[i] - evar)
        + I * rho_fld_np[i]  
        + I * hhat[i] * (rhoha[i] + rhoda[i] + rhodb[i] + rhog[i] + rhog_exp[i]);
    extern Cavity *channel;
    if (channel != NULL) F += I * rho0 * channel->rho_hat[i];
    numer = wpl[i] - lam_imag * F;
    if (do_CL)
      numer += etap[i];
    wpl[i] = numer;
  }

  // Update AB fields //
  if (chi_ab != 0.0) {
    if (do_CL) {
      generate_1s_noise(etap, lam_real);
      generate_1s_noise(etam, lam_imag);
    }

    for (i=0; i<ML; i++) {
      // AB+ //
      F = 2.0 * rho0 * wabp[i] / chi_ab
          + hhat[i] * ( 
            + (rhoha[i] + rhoda[i])  * (negative_chi_ab_flag ? -1 : I  )
            + rhodb[i]               * (negative_chi_ab_flag ? -1 : I  )
              ); 
      numer = wabp[i] - lam_imag * F;
      if (do_CL)
        numer += etap[i];
      wabp[i] = numer;

      // AB- //
      F = 2.0 * rho0 * wabm[i] / chi_ab
          + hhat[i] * ( 
            + (rhoha[i] + rhoda[i])  * (negative_chi_ab_flag ? I : -1  )
            + rhodb[i]               * (negative_chi_ab_flag ? -I : 1  )
              ); 
      numer = wabm[i] - lam_real * F;
      if (do_CL)
        numer += etam[i];
      wabm[i] = numer;
    }
  }

  else {
    for (i=0; i<ML; i++)
      wabp[i] = wabm[i] = 0.0 ;
  }  // End of Update AB Fields


  // Update BC fields //
  if (chi_bc != 0.0) {
    if (do_CL) {
      generate_1s_noise(etap, lam_imag);
      generate_1s_noise(etam, lam_real);
    }

    for (i=0; i<ML; i++) {
      // BC+ //
      F = 2.0 * rho0 * wbcp[i] / chi_bc
          + hhat[i] * ( 
            + (rhodb[i])  * (negative_chi_bc_flag ? -1 : I  )
            + rhog_exp[i] * (negative_chi_bc_flag ? -1 : I  )
              );
        numer = wbcp[i] - lam_imag * F;
      if (do_CL)
        numer += etap[i];
      wbcp[i] = numer;

      // BC- //
      F = 2.0 * rho0 * wbcm[i] / chi_bc
          + hhat[i] * ( 
            + (rhodb[i])  * (negative_chi_bc_flag ? I : -1  )
            + rhog_exp[i] * (negative_chi_bc_flag ? -I : 1  )
              ); 
      numer = wbcm[i] - lam_real * F;
      if (do_CL)
        numer += etam[i];
      wbcm[i] = numer;
    }
  }

  else {
    for (i=0; i<ML; i++)
      wbcp[i] = wbcm[i] = 0.0 ;
  }  // End of Update BC Fields


  // Update AC fields //
  if (chi_ac != 0.0) {
    if (do_CL) {
      generate_1s_noise(etap, lam_real);
      generate_1s_noise(etam, lam_imag);
    }

    for (i=0; i<ML; i++) {
      // AC+ //
      F = 2.0 * rho0 * wacp[i] / chi_ac
          + hhat[i] * ( 
            + (rhoha[i] + rhoda[i])  * (negative_chi_ac_flag ? -1 : I  )
            + rhog_exp[i]            * (negative_chi_ac_flag ? -1 : I  )
              );       
        numer = wacp[i] - lam_imag * F;
      if (do_CL)
        numer += etap[i];
      wacp[i] = numer;

      // AC- //
      F = 2.0 * rho0 * wacm[i] / chi_ac
          + hhat[i] * ( 
            + (rhoha[i] + rhoda[i])  * (negative_chi_ac_flag ? I : -1  )
            + rhog_exp[i]            * (negative_chi_ac_flag ? -I : 1  )
              );       
        numer = wacm[i] - lam_real * F;
      if (do_CL)
        numer += etam[i];
      wacm[i] = numer;
    }
  }

  else {
    for (i=0; i<ML; i++)
      wacp[i] = wacm[i] = 0.0 ;
  }  // End of Update AC Fields


  // Here we de-transform the variables so that they are normal variables. 
  fft_bck_wrapper(wpl, wpl);

  if (chi_ab != 0.0) {
    fft_bck_wrapper(wabm, wabm);
    fft_bck_wrapper(wabp, wabp);
  }

  if (chi_bc != 0.0) {
    fft_bck_wrapper(wbcm, wbcm);
    fft_bck_wrapper(wbcp, wbcp);
  }
  
  if (chi_ac != 0.0) {
    fft_bck_wrapper(wacm, wacm);
    fft_bck_wrapper(wacp, wacp);
  }

//     if (sigma > 0.0 && nFP > 0.0)
//     fft_bck_wrapper(rhog, rhog);

//   if (sigma > 0.0 && n_exp_nr > 0.0)
//     fft_bck_wrapper(rhog_exp, rhog_exp);

//  if (do_fld_np)
//     fft_bck_wrapper(rho_fld_np, rho_fld_np);


  calc_poly_density() ;
}

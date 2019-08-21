//This file has not been updated to include the diblock system.

#include "globals.h"
#include "cavity.hpp"
void generate_1s_noise(complex<double> *, double);

void update_1s()
{

  int i;
  complex<double> evar, F, A, B, numer, denom;

  // printf("1s update method not supported yet!\n");
  // exit(5);

  if (nD > 0.0)
  {
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

  if (chi_ab != 0.0)
  {
    fft_fwd_wrapper(wabm, wabm);
    fft_fwd_wrapper(wabp, wabp);
  }

  if (chi_bc != 0.0)
  {
    fft_fwd_wrapper(wbcm, wbcm);
    fft_fwd_wrapper(wbcp, wbcp);
  }

  if (chi_ac != 0.0)
  {
    fft_fwd_wrapper(wacm, wacm);
    fft_fwd_wrapper(wacp, wacp);
  }

  if (do_CL)
    generate_1s_noise(etap, lam_imag);

  // Update w+ field //
  for (i = 0; i < ML; i++)
  {
    // The 1 at the k=0 mode takes care of the delta function arising from
    // taking a fourier transform of a constant
    if (i == 0 && myrank == 0)
      evar = 1.0;
    else
      evar = 0.0;

    // F = (kappa <= 0.0 ? 0.0 : rho0*wpl[i]/kappa)
    F = I * rho0 * (surfH[i] + exp_nrH[i] - evar) + I * rho_fld_np[i] + I * hhat[i] * (rhoha[i] + rhoda[i] + rhog[i] + rhog_exp[i] + rhodb[i]);
    extern Cavity *channel;
    if (channel != NULL)
      F += I * rho0 * channel->rho_hat[i];
    A = (kappa <= 0.0 ? 0.0 : rho0 / kappa) 
    + phiD * rho0 * N * hhat[i] * hhat[i] * (gaa[i] + 2.0 * gab[i] + gbb[i]) 
    + (phiH * Nah* gd[i] + phiG * Ng * gc[i]) * rho0 * hhat[i] * hhat[i] ;
    if (nFP > 0.0 && np_type == 1)
      A += np_frac * rho0 * Gamma_iso[i] * Gamma_iso[i];
    if (kappa > 0)
      B = A - rho0 / kappa;
    numer = wpl[i] - lam_imag * (F - B * wpl[i]);
    if (do_CL)
      numer += etap[i];
    denom = 1.0 + lam_imag * A;
    wpl[i] = numer / denom;
  }

  // Update AB fields //
  if (chi_ab != 0.0)
  {
    if (do_CL)
    {
      generate_1s_noise(etap, lam_real);
      generate_1s_noise(etam, lam_imag);
    }

    for (i= 0; i < ML; i++)
    {
      if (!negative_chi_ab_flag) // Positive χ 
      {  
        // AB+ //
        // F = 2.0 * rho0 * wabp[i] / chi
        F = I * hhat[i] * (rhoda[i] + rhodb[i] + rhoha[i]);
        // If NP is A chemistry instead of neutral, add np contributions
        if (np_chem== 1)
          F += -rho_fld_np[i] - rho0 * exp_nrH[i];

        A  = 2.0 * rho0 / chi_ab 
          + phiD * rho0 * N * hhat[i] * hhat[i] * 
            (gaa[i] + 2.0 * gab[i] + gbb[i]) ;
          // (gaa[i] + gab[i]) ;
        if (nFP > 0 && np_chem == 1 && np_type == 1)
          A += np_frac * rho0 * Gamma_iso[i] * Gamma_iso[i];
        B  = A - 2.0 * rho0 / chi_ab;
        numer = wabp[i] - lam_imag * (F - B * wabp[i]);
        if (do_CL)
          numer += etap[i];
        denom  = 1.0 + lam_imag * A;
        wabp[i]   = numer / denom;

        // AB- //
        // F = 2.0 * rho0 * wabm[i] / chi
         F =  hhat[i] * (rhodb[i] - rhoda[i] - rhoha[i] );
        // If NP is A chemistry instead of neutral, add contributions
        if (np_chem == 1)
          F += I * rho_fld_np[i] + I * rho0 * exp_nrH[i];
        numer = wabm[i] - lam_real * F;
        if (do_CL)
          numer += etam[i];
        denom = 1.0 + lam_real * 2.0 * rho0 / chi_ab;
        wabm[i] = numer / denom;

      } 
      else
      {  
        // The case for when there is attraction between A and B. 
        // AB+ //
       F  = -hhat[i] * (rhoda[i] + rhodb[i] + rhoha[i]);
        numer = wabp[i] - lam_real * F;
        denom = 1.0 + lam_real * 2.0 * rho0 / chi_ab;
        wabp[i] = numer / denom;

         // AB- //
        // F = 2.0 * rho0 * wabp[i] / chi
        F  = I * hhat[i] * (rhoda[i] - rhodb[i] + rhoha[i]);
        // If NP is A chemistry instead of neutral, add NP contributions
        if (np_chem == 1)
          F += -rho_fld_np[i] - rho0 * exp_nrH[i];

        A  = 2.0 * rho0 / chi_ab + phiD * rho0 * N * hhat[i] * hhat[i] * 
        // (gaa[i] + 2.0 * gab[i] + gbb[i]) ;
        (gaa[i] - gbb[i] ) ;
        if (nFP > 0 && np_chem == 1 && np_type == 1)
          A += np_frac * rho0 * Gamma_iso[i] * Gamma_iso[i];
        B = A - 2.0 * rho0 / chi_ab;
        numer = wabm[i] - lam_imag * (F - B * wabm[i]);
        if (do_CL)
          numer += etap[i];
        denom = 1.0 + lam_imag * A;
        wabm[i] = numer / denom;       
      }
    }
  }
  else
  {
    for (i = 0; i < ML; i++)
      wabp[i] = wabm[i] = 0.0;
  }

  if (chi_bc != 0.0)
  {
    if (do_CL)
    {
      generate_1s_noise(etap, lam_real);
      generate_1s_noise(etam, lam_imag);
    }

    for (i = 0; i < ML; i++)
    {
      if (!negative_chi_bc_flag)   //This is the case for when the system does not have attractive interactions between A and B. Positive chi
      {
        // BC+ //
        F = I * hhat[i] * (rhog_exp[i] + rhodb[i] );
        // If NP is A chemistry instead of neutral, add NP contributions
        // if (np_chem == 1)
        //   F += -rho_fld_np[i] - rho0 * exp_nrH[i];

        A = 2.0 * rho0 / chi_bc + phiD * rho0 * N * hhat[i] * hhat[i] * ( 
        // + 2.0 * gab[i]   
        // + gbb[i]) ;
        + gbb[i]) 
        + phiG * rho0 * Ng * hhat[i] * hhat[i] * gc[i]  ;
        if (nFP > 0 && np_chem == 1 && np_type == 1)
          A += np_frac * rho0 * Gamma_iso[i] * Gamma_iso[i];
        B = A - 2.0 * rho0 / chi_bc;
        numer = wbcp[i] - lam_imag * (F - B * wbcp[i]);
        if (do_CL)
          numer += etap[i];
        denom = 1.0 + lam_imag * A;
        wbcp[i] = numer / denom;

        // BC- //
        // F = 2.0 * rho0 * wabm[i] / chi
        F =  hhat[i] * (rhog_exp[i] - rhodb[i] );
        // If NP is A chemistry instead of neutral, add NP contributions
        if (np_chem == 1)
          F += I * rho_fld_np[i] + I * rho0 * exp_nrH[i];
        numer = wbcm[i] - lam_real * F;
        if (do_CL)
          numer += etam[i];
        denom = 1.0 + lam_real * 2.0 * rho0 / chi_bc;
        wbcm[i] = numer / denom;
      }
      else{  // The case for when there is attraction between B and the graft.
        // BC+ //
       F = -hhat[i] * (rhodb[i] +  rhog_exp[i]);
        numer = wbcp[i] - lam_real * F;
        if (do_CL)
          numer += etap[i];
        denom = 1.0 + lam_real * 2.0 * rho0 / chi_bc;
        wbcp[i] = numer / denom;
        
        // BC- //
        // F = I * hhat[i] * (rhoda[i] + rhoha[i] - rhodb[i] );
        F =  -I*hhat[i] * (rhog_exp[i] - rhodb[i] );
        // If NP is A chemistry instead of neutral, add NP contributions
        if (np_chem == 1)
          F += -rho_fld_np[i] - rho0 * exp_nrH[i];

        A = 2.0 * rho0 / chi_bc 
        + phiD * rho0 * N * hhat[i] * hhat[i] * (gbb[i] )
         - phiG * Ng * hhat[i] * hhat[i] * rho0 * gc[i]
         ;
        if (nFP > 0 && np_chem == 1 && np_type == 1)
          A += np_frac * rho0 * Gamma_iso[i] * Gamma_iso[i];
        B = A - 2.0 * rho0 / chi_bc;
        numer = wbcm[i] - lam_imag * (F - B * wbcm[i]);
        if (do_CL)
          numer += etap[i];
        denom = 1.0 + lam_imag * A;
        wbcm[i] = numer / denom;
      }
    }
  }
  else
  {
    for (i = 0; i < ML; i++)
      wbcp[i] = wbcm[i] = 0.0;
  }


  if (chi_ac != 0.0)
  {
    if (do_CL)
    {
      generate_1s_noise(etap, lam_real);
      generate_1s_noise(etam, lam_imag);
    }

    for (i = 0; i < ML; i++)
    {
      if (!negative_chi_ac_flag)   //This is the case for when the system does not have attractive interactions between A and B.
      {
        // AC+ //
        F = I * hhat[i] * (rhog_exp[i] + rhoda[i] );
        // If NP is A chemistry instead of neutral, add NP contributions
        // if (np_chem == 1)
        //   F += -rho_fld_np[i] - rho0 * exp_nrH[i];

        A = 2.0 * rho0 / chi_ac + phiD * rho0 * N * hhat[i] * hhat[i] * (gaa[i] 
        /* + 2.0 * gab[i]   // This was removed because Rob said it did not have the middle part of the code */
        /* + gbb[i]) ; */
      	 ) 
        +  phiG * Ng * hhat[i] * hhat[i] * rho0 * gc[i] 
        ;
        if (nFP > 0 && np_chem == 1 && np_type == 1)
          A += np_frac * rho0 * Gamma_iso[i] * Gamma_iso[i];
        B = A - 2.0 * rho0 / chi_ac;
        numer = wacp[i] - lam_imag * (F - B * wacp[i]);
        if (do_CL)
          numer += etap[i];
        denom = 1.0 + lam_imag * A;
        wacp[i] = numer / denom;

        // AC- //
        // F = 2.0 * rho0 * wabm[i] / chi
        F =  hhat[i] * (rhog_exp[i] - rhoda[i] );
        // If NP is A chemistry instead of neutral, add NP contributions
        if (np_chem == 1)
          F += I * rho_fld_np[i] + I * rho0 * exp_nrH[i];
        numer = wacm[i] - lam_real * F;
        if (do_CL)
          numer += etam[i];
        denom = 1.0 + lam_real * 2.0 * rho0 / chi_ac;
        wacm[i] = numer / denom;
      }
      else{  // The case for when there is attraction between B and the graft.
        // AC+ //
       F = -hhat[i] * (rhoda[i] +  rhog_exp[i]);
        numer = wacp[i] - lam_real * F;
        if (do_CL)
          numer += etap[i];
        denom = 1.0 + lam_real * 2.0 * rho0 / chi_ac;
        wacp[i] = numer / denom;
        
        // AC- //
        // F = I * hhat[i] * (rhoda[i] + rhoha[i] - rhodb[i] );
        F =  -I*hhat[i] * (rhog_exp[i] - rhoda[i] );
        // If NP is A chemistry instead of neutral, add NP contributions
        if (np_chem == 1)
          F += -rho_fld_np[i] - rho0 * exp_nrH[i];

        A = 2.0 * rho0 / chi_ac + phiD * rho0 * N * hhat[i] * hhat[i] * (gaa[i] )
         - phiG * Ng * hhat[i] * hhat[i] * rho0  * gc[i]
         ;
        if (nFP > 0 && np_chem == 1 && np_type == 1)
          A += np_frac * rho0 * Gamma_iso[i] * Gamma_iso[i];
        B = A - 2.0 * rho0 / chi_ac;
        numer = wacm[i] - lam_imag * (F - B * wacm[i]);
        if (do_CL)
          numer += etap[i];
        denom = 1.0 + lam_imag * A;
        wacm[i] = numer / denom;
      }
    }
  }
  else
  {
    for (i = 0; i < ML; i++)
      wacp[i] = wacm[i] = 0.0;
  }

fft_bck_wrapper(wpl, wpl);

if (chi_ab != 0.0)
{
  fft_bck_wrapper(wabm, wabm);
  fft_bck_wrapper(wabp, wabp);
}

if (chi_bc != 0.0)
{
  fft_bck_wrapper(wbcm, wbcm);
  fft_bck_wrapper(wbcp, wbcp);
}

if (chi_ac != 0.0)
{
  fft_bck_wrapper(wacm, wacm);
  fft_bck_wrapper(wacp, wacp);
}

calc_poly_density();
}

// Generates Gaussian noise in k-space with appropriate statistics //
// for the 1s updating scheme.
void generate_1s_noise(complex<double> *et, double lambda)
{

  int i;
  double scale = sqrt(2.0 * lambda);
  for (i = 0; i < Dim; i++)
    scale /= sqrt(dx[i]);

  for (i = 0; i < ML; i++)
    et[i] = scale * gasdev2();

  fft_fwd_wrapper(et, et);
}

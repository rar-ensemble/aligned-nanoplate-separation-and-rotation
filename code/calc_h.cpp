#include "globals.h"
#include "cavity.hpp"

complex<double> wpl_part() ;
complex<double> wab_part() ;
complex<double> wac_part() ;
complex<double> wbc_part() ;

complex<double> calc_H() {

  int i;

  Hcur = wpl_part() + wbc_part() + wac_part()  + wab_part()  ;

  // Add diblock part if applicable
  if (nD > 0.0)
    Hcur += - nD*log(Qd) ;

  // Add homopolymer part if applicable
  if (nAH > 0.0)
    Hcur += - nAH * log(Qha);

  // Add field-based nanoparticle part if applicable. Note that smwp_min,
  // which was originally the minimum value of the smeared wp field (smwp),
  // is subtracted off here because when we subtracted smwp_min from smwp,
  // that made the value stored in Qp actually Qp * exp(smwp_min).
  if (do_fld_np)
    Hcur += - nP * (log(Qp) - smwp_min);

  // Add contribution of chains grafted to explicit nanoparticles
  if (n_exp_nr > 0) {
    if (sigma > 0.0) {
      Hcur += - double(n_exp_nr) * ng_per_np * real(Qga_exp);
    }
  }

  // Exit if H is NaN
  if ( Hcur != Hcur) {
    double wpl_part_tmp = real(wpl_part());
    double wab_part_tmp = real(wab_part());
    double wbc_part_tmp = real(wbc_part());
    double wac_part_tmp = real(wac_part());
    if (myrank == 0) {
      printf("Hcur is NaN!\n");
      printf("real(Hcur) = %lf\n", real(Hcur));
      printf("real(wpl_part) = %lf\n", wpl_part_tmp );
      printf("real(wab_part)=%lf\n", wab_part_tmp );
      printf("real(wac_part)=%lf\n", wac_part_tmp );
      printf("real(wbc_part)=%lf\n", wbc_part_tmp );
      printf("real(-nD*log(Qd))=%lf\n", real(-nD*log(Qd)) );
      printf("real(-nAH*log(Qha))=%lf\n", real(-nAH*log(Qha)) );
      printf("real(-nP*[log(Qp)-smwp_min])=%lf\n",
              real(-nP*(log(Qp)-smwp_min)) );
      printf("real(-n_exp_grafts*Qga_exp) = %lf\n",
              - double(n_exp_nr) * ng_per_np * real(Qga_exp) );
      printf("real(Qga_exp) = %lf\n", real(Qga_exp));
      printf("iter = %d\n", iter);
      printf("Qd = %lf\n", real(Qd));
    }
    exit(1);
  }

  return Hcur ;
}

complex<double> wpl_part() {

  int i;

  for (i=0; i<ML; i++) {
    tmp2[i] = ( kappa<=0.0 ? 0.0 : wpl[i]*wpl[i]*rho0/kappa/2.0 )
      - I * rho0 * ( 1.0 - rho_surf[i] -rho_exp_nr[i]  ) * wpl[i];
  }

  extern Cavity *channel;

  if (channel != NULL)
  {
    for (i=0; i<ML; i++)
    {
      tmp2[i] += I * rho0 * channel->rho[i] * wpl[i];
    }
  }

  return integ_trapPBC(tmp2);

}

complex<double> wab_part() {

  if ( chi_ab == 0.0 ) {
     return 0.0 ;
  }
  else {
    int i ;

    for ( i=0 ; i<ML ; i++ )
      tmp[i] = ( wabp[i] * wabp[i] + wabm[i] * wabm[i] ) * rho0 / chi_ab * (negative_chi_ab_flag ? -1.0 : 1.0);

    return integ_trapPBC( tmp ) ;
  }
}

complex<double> wac_part() {

  if ( chi_ac == 0.0 ) {
     return 0.0 ;
  }
  else {
    int i ;

    for ( i=0 ; i<ML ; i++ )
      tmp[i] = ( wacp[i] * wacp[i] + wacm[i] * wacm[i] ) * rho0 / chi_ac * (negative_chi_ac_flag ? -1.0 : 1.0);

    return integ_trapPBC( tmp ) ;
  }
}

complex<double> wbc_part() {

  if ( chi_bc == 0.0 ) {
     return 0.0 ;
  }
  else {
    int i ;

    for ( i=0 ; i<ML ; i++ )
      tmp[i] = ( wbcp[i] * wbcp[i] + wbcm[i] * wbcm[i] ) * rho0 / chi_bc * (negative_chi_bc_flag ? -1.0 : 1.0);

    return integ_trapPBC( tmp ) ;
  }
}

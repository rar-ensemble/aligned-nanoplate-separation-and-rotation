#include "globals.h"
/*
 * Note: "grafts" is expected to be the *single particle*, *normalized*
 * grafted density *in k space*. Its real-space integral should be unity. 
 * The grafted particle field will be constructed as:
 * 
 * WP_tot(r) = WP(r) - n_g / n_P * ( grafts * log( q[Ng-1] )(r)
 * 
 * where grafts is convolved with log(q).
 *
 */
complex<double> shift_field( complex<double>* ) ;
void write_data_bin(const char*, complex<double>*);

void graft_homopoly_free_ends(complex<double>* smwg_here, int Ng,
                              complex<double>** q) {

  // Equation 18 in H. Chao 2014 Soft Matter Paper
  // wg is either wa, wb, or wc in real space depending on whether grafts are A, B or C

  int i, j;

  // Initial condition
  for (i=0; i<ML; i++) 
    q[0][i] = exp(-smwg_here[i]) ;

  // Loop over all Ng segments in grafted chains
  for (j=1; j<Ng; j++) {

    // Convolve q with bond potential
    fft_fwd_wrapper(q[j-1], q[j]);
    for (i=0; i<ML; i++) 
      q[j][i] *= poly_bond_fft[i];
    fft_bck_wrapper(q[j], q[j]);

    // Multiply by exp(-wa) or exp(-wb)
    for (i=0; i<ML; i++)
      q[j][i] *= exp(-smwg_here[i]);

  }

}

complex<double> grafted_exp_nps(
    double ngrafts_per_np,
    int n_exp_np,
    int Ng,
    complex<double> *smwg_here,
    complex<double> *expl_grafts,
    complex<double> **q_loc,
    complex<double> **qdag_loc,
    complex<double> *rhog_here) {

  int i, j;

  // Initial condition for the complimentary graft propagator
  for ( i=0 ; i<ML ; i++ ) 
    qdag_loc[0][i] = expl_grafts[i] * exp( -smwg_here[i] ) / q_loc[Ng-1][i] ;

  // Rest of the graft propagator
  for ( j=1 ; j<Ng ; j++ ) {

    fft_fwd_wrapper( qdag_loc[j-1] , qdag_loc[j] ) ;

    for ( i=0 ; i<ML ; i++ ) 
      qdag_loc[j][i] *= poly_bond_fft[i] ;

    fft_bck_wrapper( qdag_loc[j] , qdag_loc[j] ) ;

    for ( i=0 ; i<ML ; i++ )
      qdag_loc[j][i] *= exp( -smwg_here[i] ) ;

  }

  // Calculate the grafted chain density
  for ( i=0 ; i<ML ; i++ ) {
    rhog_here[i] = 0.0 ;
    for ( j=0 ; j<Ng ; j++ ) 
      rhog_here[i] += q_loc[j][i] * qdag_loc[Ng-j-1][i] ;

    rhog_here[i] *= ngrafts_per_np * n_exp_np * exp( +smwg_here[i] ) ;
  }
  int counter=0;
  // Factor for the partition function
  for ( i=0 ; i<ML ; i++ ) {
    tmp[i] =  (real(expl_grafts[i]) <= 0.0 ? 0.0 : expl_grafts[i] * log( q_loc[Ng-1][i] ));
    if (real(q_loc[Ng-1][i]) <= 0.0){
      printf("q_loc[Ng-1][%d] = %g\n", i, real(q_loc[Ng-1][i]));
      cout << "wpl" << wpl[i] << endl;
      cout << "wabm" << wabm[i] << endl;
      cout << "wabp" << wabp[i] << endl;
      cout << "wacm" << wacm[i] << endl;
      cout << "wacp" << wacp[i] << endl;
      cout << "wbcm" << wbcm[i] << endl;
      cout << "wbcp" << wbcp[i] << endl;
      cout << "smwg" << smwg_here[i] << endl;

      if (++counter >= 5){
          write_data_bin("rhoda", rhoda);
          write_data_bin("rhodb", rhodb);
          
          write_data_bin("wa", wa);
          write_data_bin("wb", wb);
          write_data_bin("wpl", wpl);
          
          exit(4);
      }
      }
      
      
  }

  Qga_exp = integ_trapPBC(tmp);

  return Qga_exp ;

}

complex<double> grafted_fld_nps ( 
    complex<double> *WP,  // Potential felt by particle cores
    complex<double> *WG,  // Potential felt by grafts
    complex<double> *graft_pts,  // Graft distribution for one particle
    complex<double> *gammaNP,  // nanoparticle shape function
    complex<double> **q, complex<double> **qdag, // Propagators for grafts
    complex<double> *rp , complex<double> *smrp, // Nanoparticle densities
    complex<double> &shift_WP, // Shift of chx pot field of the NPs
    complex<double> *rhg, // Grafted chain densities
    double npar,    // number of particles
    double ngrafts_per_np, // Number of grafted chains 
    int Ng // Length of grafts
    )

{

  int i,j ;

  complex<double> Qgp ;

  // Distribution of bare particles //
  if ( sigma == 0.0 ) {

    shift_WP = shift_field( WP ) ;

    for ( i=0 ; i<ML ; i++ )
      tmp[i] = exp( -WP[i] ) ;

    Qgp = integ_trapPBC( tmp ) / V ;
 
    for ( i=0 ; i<ML ; i++ ) 
      rp[i] = tmp[i] * npar / ( Qgp * V ) ;

    fft_fwd_wrapper( rp , smrp ) ;
    for ( i=0 ; i<ML ; i++ ) 
      smrp[i] = gammaNP[i] * smrp[i] ;

    fft_bck_wrapper( smrp , smrp  ) ;

    return Qgp ;

  }

  // Grafted nanoparticles //
  else {

    // Construct the full particle field 
    for ( i=0 ; i<ML ; i++ ) 
      tmp[i] = ngrafts_per_np * log( q[Ng-1][i] ) ;

    fft_fwd_wrapper( tmp, tmp ) ;
    for ( i=0 ; i<ML ; i++ ) 
      tmp[i] *= graft_pts[i] ;
    fft_bck_wrapper( tmp, tmp ) ;

    for ( i=0 ; i<ML ; i++ ) 
      WP[i] = WP[i] - tmp[i] ;

    shift_WP = shift_field( WP ) ;

    // Calculate the particle partition function
    for ( i=0 ; i<ML ; i++ )
      tmp[i] = exp( - WP[i] ) ;

    Qgp = integ_trapPBC( tmp ) / V ;

    // Calculate the particle core densities //
    for ( i=0 ; i<ML ; i++ ) 
      rp[i] = tmp[i] * npar / ( Qgp * V ) ;

    fft_fwd_wrapper( rp , tmp ) ;
    for ( i=0 ; i<ML ; i++ ) 
      smrp[i] = gammaNP[i] * tmp[i] ;

    fft_bck_wrapper( smrp , smrp ) ;

    // Calculate the total graft density
    for ( i=0 ; i<ML ; i++ ) 
      tmp[i] = graft_pts[i] * tmp[i] / npar;

    fft_bck_wrapper( tmp , tmp ) ;

    // Initial condition for the complimentary graft propagator
    for ( i=0 ; i<ML ; i++ ) 
      qdag[0][i] = tmp[i] * exp( -WG[i] ) / q[Ng-1][i] ;

    // Rest of the graft propagator
    for ( j=1 ; j<Ng ; j++ ) {

      fft_fwd_wrapper( qdag[j-1] , qdag[j] ) ;

      for ( i=0 ; i<ML ; i++ ) 
        qdag[j][i] *= poly_bond_fft[i] ;

      fft_bck_wrapper( qdag[j] , qdag[j] ) ;

      for ( i=0 ; i<ML ; i++ )
        qdag[j][i] *= exp( -WG[i] ) ;

    }

    // Finally, calculate the grafted chain density //
    for ( i=0 ; i<ML ; i++ ) {
      rhg[i] = 0.0 ;
      for ( j=0 ; j<Ng ; j++ ) 
        rhg[i] += q[j][i] * qdag[Ng-j-1][i] ;

      rhg[i] *= ngrafts_per_np * npar * exp( +WG[i] ) ;
    }

    return Qgp ;

  }// if (sigma > 0)

}


complex<double> shift_field( complex<double> *w  ) {

  double wmin, wtmp ;
  int i;

  wtmp = 12342.0 ;

  for ( i=0 ; i<ML ; i++ )
    if ( real( w[i] ) < wtmp )
      wtmp = real( w[i] ) ;

#ifdef PAR
  MPI_Allreduce( &wtmp, &wmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD ) ;
#else
  wmin = wtmp;
#endif

  for ( i=0 ; i<ML ; i++ )
    w[i] -= wmin ;

  return complex<double>(wmin, 0.0) ;

}

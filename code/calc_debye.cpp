#include "globals.h"

// Check these routines - not working in Rg coordinates anymore //
// These routines need to receive factors of N or Rg to scale k //
void calc_gaa(complex<double> *g , double fD ) {
  double x, k2, k[Dim];
  int i;

  for (i=0; i<ML; i++) {

    k2 = get_k( i , k ) ;
    k2 *= Rg2;
    g[i] = 2.0 * ( exp( -k2*fD ) + fD*k2 - 1.0 ) / k2 / k2 ;

  }
  if ( myrank == 0 )
    g[0] = fD * fD ;

}

void calc_gbb(complex<double> *g , double fD ) {
  double x, k2, k[Dim];
  int i;

  for (i=0; i<ML; i++) {

    k2 = get_k( i , k ) ;
    k2 *= Rg2;
    g[i] = 2.0 * ( exp( -k2*(1.0-fD) ) + (1.0-fD)*k2 - 1.0 ) / k2 / k2 ;

  }
  if ( myrank == 0 )
    g[0] = (1.0-fD) * ( 1.0 - fD ) ;

}

void calc_gab( complex<double> *g , double fD ) {

  double k2, k[Dim] ;
  int i;

  for ( i=0 ; i<ML ; i++ ) {
    k2 = get_k ( i , k );
    k2 *= Rg2;
    g[i] = ( 1.0 - exp( -k2*fD ) ) * ( 1.0 - exp( -k2*(1.0-fD)) ) / k2 / k2 ;
  }

  if ( myrank == 0 )
    g[0] = fD * ( 1.0 - fD );

}

void calc_gd(complex<double> *g, double alpha) {
  double dm1, k2, k[Dim];
  int i;

  for (i=0; i<ML; i++) {
    
    k2 = get_k(i,k);
    dm1 = k2 * Rg2 * alpha;
    if ( dm1 == 0.0 )
      g[i] = 1.0 ;
    else 
      g[i] = 2.0*(exp(-dm1) + dm1 - 1.0) / dm1 / dm1;

  }

}



void calc_discrete_debye( complex<double> *g , int N ) {

  double dm1, k2, k[Dim], pref = 1.0 / double( N-1 ) , arg ;
  int i, j , m ;

  double Rg2 = double(N) / 6.0 ;

  for (i=0; i<ML; i++) {
    
    k2 = get_k(i,k);
    dm1 = k2 * Rg2 / double( N-1 ) ;

    g[i] = 0.0 ;
    for ( j=0 ; j<N ; j++ ) 
      for ( m=0 ; m<N ; m++ ) {
        arg = -dm1 * double( abs(j-m) ) / double(N-1) ;
        g[i] += exp( arg ) ;
      }

    g[i] *= 1.0 / double( (N) * (N) ) ;

  }

}

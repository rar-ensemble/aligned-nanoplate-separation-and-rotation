#include "globals.h"


complex<double> diblock_discrete( complex<double> *WA, complex<double> *WB, 
    complex<double> **q , complex<double> **qdag, int Ns , int Na ) {

  int i,n,j,k, nn[Dim];

  complex<double> *qf, *h, *a, *w ;

  double temp;
  double arg, indx, indy, indz, kvec[Dim];

  qf = (complex<double>*) fftw_malloc(sizeof(complex<double>)*size);
  h =  (complex<double>*) fftw_malloc(sizeof(complex<double>)*size);
  a =  (complex<double>*) fftw_malloc(sizeof(complex<double>)*size);

  // Initial condition on q(s=0, x)
  for (i=0; i<ML; i++) {
    if ( Na != 0 )
      q[0][i] = exp( -WA[i] );
    else
      q[0][i] = exp( -WB[i] ) ;
    
    if ( Na != N )
      qdag[0][i] = exp( -WB[i] );
    else
      qdag[0][i] = exp( -WA[i] );
  }

  for (n=1; n<Ns; n++) {

    ///////////////////////////
    // Forward propagator, q //
    ///////////////////////////

    // Switching function for the field
    if ( Na > 0 && n < Na ) 
      w = WA ;
    else 
      w = WB ;
    
    // Go to k-space to evaluate the convolution //
    fft_fwd_wrapper(q[n-1], a);

    for (j=0; j<ML; j++) 
      h[j] = a[j] * poly_bond_fft[j] ;

    // IFFT h to get q^{n+1/2}
    fft_bck_wrapper(h, qf);

    // Final transformation: q^{n+1} = q^{n+1/2}*exp(-w(x))
    for (j=0; j<ML; j++)
      q[n][j] = qf[j] * exp( -w[j] );

    ////////////////////////////////////////////////////
    // Section for the complimentary propagator, qdag //
    ////////////////////////////////////////////////////

    if ( Na != Ns && Na != 0 ) {

      if ( n <= ( Ns - Na - 1 ) ) 
        w = WB ;
      else 
        w = WA ;
      
      // To k-space
      fft_fwd_wrapper(qdag[n-1], a);
 
      // Convolve
      for (j=0; j<ML; j++) 
        h[j] = a[j] * poly_bond_fft[j] ;
 
      // to real space
      fft_bck_wrapper(h, qf);
 
      // Apply field
      for (j=0; j<ML; j++)
        qdag[n][j] = qf[j] * exp( -w[j] );
    }
    else
      for ( j=0 ; j<ML ; j++ )
        qdag[n][j] = q[n][j] ;

  }//for (n=0; n<N...

  complex<double> Q = integ_trapPBC(q[N-1]) / V;
  
  free(qf);
  free(h);
  free(a);

  return Q;

}


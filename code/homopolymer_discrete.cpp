#include "globals.h"

complex<double> homopolymer_discrete( complex<double> *WA,  
    complex<double> **q , 
    int N   ) {

  int i,n,j,k, nn[Dim];
  complex<double> *qf, *h, *a, *w ;
  double temp;
  double arg, indx, indy, indz, kvec[Dim];
  int Nflip[Dim];

  for (i=0; i<Dim; i++) 
    Nflip[Dim-i-1] = Nx[i];

  fftw_plan qplan, hplan;

  qf = (complex<double>*) fftw_malloc(sizeof(complex<double>)*size);
  h =  (complex<double>*) fftw_malloc(sizeof(complex<double>)*size);
  a =  (complex<double>*) fftw_malloc(sizeof(complex<double>)*size);

  // Initial condition on q(s=0, x)
  for (i=0; i<ML; i++) 
      q[0][i] = exp( -WA[i] );

  for (n=1; n<N; n++) {

    ///////////////////////////
    // Forward propagator, q //
    ///////////////////////////
    
    // Go to k-space to evaluate the convolution //
    fft_fwd_wrapper(q[n-1], a);

    for (j=0; j<ML; j++)
      h[j] = a[j] * poly_bond_fft[j] ;

    // IFFT h to get q^{n+1/2}
    fft_bck_wrapper(h, qf);

    // Final transformation: q^{n+1} = q^{n+1/2}*exp(-w(x))
    for (j=0; j<ML; j++)
      q[n][j] = qf[j] * exp( -WA[j] );

  } //for (n=0; n<N...

  complex<double> Q = integ_trapPBC(q[N-1]) / V;
  
  free(qf);
  free(h);
  free(a);

  return Q;

}


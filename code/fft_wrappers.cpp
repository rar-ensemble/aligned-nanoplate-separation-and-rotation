#include "globals.h"

// Forward transform
void fft_fwd_wrapper(complex<double>* in, complex<double>* out) {

  int i;

  // Store fft input
#ifdef PAR
  for (i=0; i<ML; i++) {
    fmin0[i][0] = real(in[i]);
    fmin0[i][1] = imag(in[i]);
  }
#else
  for (i=0; i<M; i++) {
    fin[i][0] = real(in[i]);
    fin[i][1] = imag(in[i]);
  }
#endif

  fftw_execute(fwd0);

  double norm = 1.0 / double(M);

  // Store fft output
#ifdef PAR
  for (i=0; i<ML; i++) {
    out[i] =( fmot0[i][0] + I * fmot0[i][1] ) * norm ;
  }
#else
  for (i=0; i<M; i++)
    out[i] = ( fout[i][0] + I * fout[i][1] ) * norm ;
#endif

}

// Backwards transform and normalization
void fft_bck_wrapper(complex<double>* in, complex<double>* out) {

  int i;
  // Store input
#ifdef PAR
  for (i=0; i<ML; i++) {
    fmin0[i][0] = real(in[i]);
    fmin0[i][1] = imag(in[i]);
  }
#else
  for (i=0; i<M; i++) {
    fin[i][0] = real(in[i]);
    fin[i][1] = imag(in[i]);
  }
#endif

  // Perform fft
  fftw_execute(fbk0);

  // Store output
#ifdef PAR
  for (i=0; i<ML; i++) {
    out[i] = (fmot0[i][0] + I * fmot0[i][1]);
  }
#else
  for (i=0; i<M; i++)
    out[i] = (fout[i][0] + I * fout[i][1]);
#endif

}

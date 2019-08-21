#include "globals.h"
using namespace std;



complex<double> integ_simpson(int np, complex<double>* dat, double h) {

  int i,j;
  complex<double> sum;

  sum = 3.0*dat[0]/8.0 + 7.0*dat[1]/6.0 + 23.0*dat[2]/24.0;

  for (i=3; i<np-3; i++)
    sum += dat[i];

  sum += 23.0*dat[np-3]/24.0 + 7.0*dat[np-2]/6.0 + 3.0*dat[np-1]/8.0;

  return sum*h;

}



complex<double> integ_trapPBC(complex<double>* dat) {
  int i;
  complex<double> sum = complex<double>(0.0,0.0);

  for (i=0; i<ML; i++)
    sum += dat[i];

  for (i=0; i<Dim; i++)
    sum *= dx[i];

#ifdef PAR
  double smr = real(sum), smi = imag(sum), rr, ri;

  MPI_Allreduce(&smr, &rr, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&smi, &ri, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  sum = rr + I*ri;
#endif

  return sum;

}




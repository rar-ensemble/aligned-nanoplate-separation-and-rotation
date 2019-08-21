#define MAIN
#include "main.hpp"

double simulate(void);

// The main routine is essentially a wrapper for the simulation routine
int main(int argc, char** argv) {

#ifdef PAR
  // MPI initialization stuff
  MPI_Init( &argc , &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank);
  MPI_Comm_size( MPI_COMM_WORLD, &nprocs);
  fftw_mpi_init();
#else
  myrank = 0;
#endif

  // Other initialization stuff
  idum = -long(time(0)) * (myrank + 1);
  read_input();
  initialize_1();
  allocate();

  // Single simulation if not doing Brent's method
  simulate();



#ifdef PAR
if (myrank == 0)
{
time_t current_time;
    char* c_time_string;

    current_time = time(NULL);

    /* Convert to local time format. */
    c_time_string = ctime(&current_time);

    printf("The program finished at %s.", c_time_string);

}
  MPI_Barrier( MPI_COMM_WORLD );
#else
time_t current_time;
    char* c_time_string;

    current_time = time(NULL);

    /* Convert to local time format. */
    c_time_string = ctime(&current_time);

    printf("The program finished at %s.", c_time_string);
#endif

#ifdef PAR
  MPI_Finalize();
#endif



  return 0;
}

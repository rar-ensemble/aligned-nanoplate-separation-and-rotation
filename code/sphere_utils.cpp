#include "globals.h"

complex<double> integ_sphere(complex<double>**);
void gauss_legendre_init( double*  , double* , double , double , int ) ;

complex<double> integ_sphere_trapPBC(complex<double> ***dat) {

  int i,j;

  for (i=0; i<Nu; i++) {
    for (j=0; j<2*Nu; j++) {
      tmp_sph[i][j] = integ_trapPBC(dat[i][j]);
//      cout<<utmp[i][j]<<endl;
    }
  }
  return integ_sphere(tmp_sph);

}

void integ_sphere_posits( complex<double> ***dat , complex<double> *out ) {

  int i,j, k;
  complex<double> sum = 0.0 ;

  for ( k=0 ; k<ML ; k++ ) {
    out[k] = 0.0 ;

    for ( i=0 ; i<Nu ; i++ ){ 
      for ( j=0 ; j<2 * Nu ; j++ ){ 

        out[k] += dat[i][j][k] * theta_weights[i] * phi_weights[j] ;
      }
    }
    //cout<<out[k]<<endl;
  }

}

complex<double> integ_sphere( complex<double> **dat ) {

  int i,j;
  complex<double> sum = 0.0 ;

  for ( i=0 ; i<Nu ; i++ ) 
    for ( j=0 ; j<2*Nu ; j++ ){
      sum += dat[i][j] * theta_weights[i] * phi_weights[j] ;
      //cout<<"integ "<<dat[i][j]<<" "<<theta_weights[i]<<" "<<phi_weights[j]<<" "<<sum<<endl; 
  }
  return sum ;
}

void sphere_init( ) {
  int i ;

  gauss_legendre_init( theta , theta_weights , 0.0 , PI ,Nu ) ;
  gauss_legendre_init( phi , phi_weights , 0.0 , 2.0*PI , 2*Nu ) ;

  for ( i=0 ; i<Nu ; i++ ){ 
    theta_weights[i] *= sin( theta[i] ) ;
   //cout<<sph_theta[i]<<endl;
 }
}

void gauss_legendre_init( double* x , double *w , double a , double b , int N) {

  FILE *inp ;
  char tt[350] ;
  char *home;
  int i,j ;

  home = getenv("HOME");

  // Make sure 'lgvalues-abscissa.txt' and 'lgvalues-weights.txt' are in ~/bin
  strcpy(tt, home);
  strcat(tt, "/bin/lgvalues-abscissa.txt");
  inp = fopen( tt , "r" );
  if ( inp == NULL ) {
    printf("Failed to find gauss-legendre abscissa file!\n");
    exit(1) ;
  }
  int line=1;
  for ( i=1 ; i < N-1 ; i++ ) {
    fgets( tt , 340 , inp );

    for ( j=0 ; j<i+1 ; j++ ) {
      fgets( tt , 340 , inp );
    }

    fgets( tt , 340 , inp );
  }

  fgets( tt , 340 , inp );

  for ( i=0 ; i<N ; i++ ) {
    fscanf( inp , "%lf" , &x[i] ) ;
    fgets( tt , 340 , inp  );
    // printf("Read %lf" , x[i] ) ;

    x[i] = ( b - a ) * x[i] / 2.0 + ( a + b ) / 2.0 ;
    if(x[i] != x[i]) {
      printf(" converted to %lf  %d radians\n" , x[i],i );
      exit(1);
    }
  }

  fclose( inp ) ;

  strcpy(tt, home);
  strcat(tt, "/bin/lgvalues-weights.txt");
  inp = fopen( tt , "r" );
  if ( inp == NULL ) {
    printf("Failed to find gauss-legendre weights file!\n");
    exit(1) ;
  }

  for ( i=1 ; i < N-1 ; i++ ) {
    fgets( tt , 340 , inp ) ;

    for ( j=0 ; j<i+1 ; j++ ) 
      fgets( tt , 340 , inp ) ;

    fgets( tt , 340 , inp ) ;
  }

  fgets( tt , 340 , inp ) ;

  for ( i=0 ; i<N ; i++ ) {
    fscanf( inp , "%lf" , &w[i] ) ;
    fgets( tt , 340 , inp  );
//    printf("Read %lf" , w[i] ) ;

    w[i] = ( b - a ) / 2.0 * w[i] ;
//    printf(" converted to %lf\n" , w[i] ) ;
  }

  fclose( inp ) ;
}

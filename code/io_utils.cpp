#include "globals.h"
void accumulate_all_averages(void);
void read_one_resume_file(FILE*, complex<double>*);
void write_avg_data(char*, complex<double>*);
void write_avg_kdata(char* , complex<double>*);
void write_kdata(char*,complex<double>*);
void write_avg_data_bin(const char*, complex<double>*);
void write_data_bin(const char*, complex<double>*);
void write_data(char*, complex<double>*);

void write_outputs() {
  int i , j;
  char nm[50];

  if (nD > 0.0) {
    fft_fwd_wrapper(rhoda, tmp);
    for (i=0; i<ML; i++)
      tmp[i] *= hhat[i];
    fft_bck_wrapper(tmp, tmp);
    write_data_bin("rhoda", tmp);
    write_data_bin("rhoda_c", rhoda);

    fft_fwd_wrapper(rhodb, tmp);
    for (i=0; i<ML; i++)
      tmp[i] *= hhat[i];
    fft_bck_wrapper(tmp, tmp);
    write_data_bin("rhodb", tmp);
    write_data_bin("rhodb_c", rhodb);
    // sprintf(nm, "rhodb_%d", iter);
    // write_data_bin(nm, tmp);
  }

  if (nAH > 0.0) {
    fft_fwd_wrapper(rhoha, tmp);
    for (i=0; i<ML; i++)
      tmp[i] *= hhat[i];
    fft_bck_wrapper(tmp, tmp);
    write_data_bin("rhoha", tmp);
    write_data_bin("rhoha_c", rhoha);
  }

  if (sigma > 0.0) {
    fft_fwd_wrapper(rhog, tmp);
    for (i=0; i<ML; i++)
      tmp[i] *= hhat[i];
    fft_bck_wrapper(tmp, tmp);
    write_data_bin("rhog", tmp);
    write_data_bin("rhog_c", rhog);
    fft_fwd_wrapper(rhog_exp, tmp);
    for (i=0; i<ML; i++)
      tmp[i] *= hhat[i];
    fft_bck_wrapper(tmp, tmp);
    write_data_bin("rhog_exp", tmp);
    write_data_bin("rhog_exp_c", rhog_exp);
  }

  if (n_exp_nr > 0 && iter == 1) {
    write_data_bin("rho_exp_nr", rho_exp_nr);
    write_data_bin("expl_grafts", expl_grafts);
  }

  if (do_fld_np) {
    write_data_bin("rho_fld_np", rho_fld_np);
    write_data_bin("rho_fld_np_c", rho_fld_np_c);
    write_data_bin("grafts", grafts);
  }

  if (do_CL && iter >= sample_wait) {
    int frame;
    frame = iter / print_freq;
    if (nD > 0.0) {
      write_avg_data_bin("avg_rhoda", avg_rhoda);
      write_avg_data_bin("avg_rhodb", avg_rhodb);
    }
    if (nAH > 0.0) {
      write_avg_data_bin("avg_rhoha", avg_rhoha);
    }
    if (do_fld_np) {
      write_avg_data_bin("avg_rho_fld_np", avg_rho_fld_np);
      write_avg_data_bin("avg_rho_fld_np_c", avg_rho_fld_np_c);
    }
    if (do_CL) {
      write_avg_data_bin("avg_rhog", avg_rhog);
      write_avg_data_bin("avg_rhog_exp", avg_rhog_exp);
    }
  }

  write_data_bin("wpl", wpl);
  if (chi_ab != 0){
    write_data_bin("wabp", wabp);
    write_data_bin("wabm", wabm);
  }
  
  if (chi_ac != 0){
    write_data_bin("wacp", wacp);
    write_data_bin("wacm", wacm);
  }

  if (chi_bc != 0){
    write_data_bin("wbcp", wbcp);
    write_data_bin("wbcm", wbcm);
  }

} // end write_outputs

// Saves average densities with the iteration number in the name
void save_averages() {
  char nm[50];
  if (nD > 0.0) {
    sprintf(nm, "avg_rhoda_%d", iter);
    write_avg_data_bin(nm, avg_rhoda);

    sprintf(nm, "avg_rhodb_%d", iter);
    write_avg_data_bin(nm, avg_rhodb);
  }
  if (nAH > 0.0) {
    sprintf(nm, "avg_rhoha_%d", iter);
    write_avg_data_bin(nm, avg_rhoha);
  }
  if (sigma > 0.0) {
    sprintf(nm, "avg_rhog_%d", iter);
    write_avg_data_bin(nm, avg_rhog);
    sprintf(nm, "avg_rhog_exp_%d", iter);
    write_avg_data_bin(nm, avg_rhog_exp);
  }
  if (do_fld_np) {
    sprintf(nm, "avg_rho_fld_np_%d", iter);
    write_avg_data_bin(nm, avg_rho_fld_np);

    sprintf(nm, "avg_rho_fld_np_c_%d", iter);
    write_avg_data_bin(nm, avg_rho_fld_np_c);
  }
}

// If running in parallel, this routine searches each file for
// the starting spatial position for this specific processor,
// then it reads in ML data points.
void read_one_resume_file(FILE *inp, complex<double> *w ) {

  int i,j, nn[Dim];
  double dr, di, bgn[Dim], dm[Dim];
  char tt[256];

#ifdef PAR
  unstack(unstack_stack(0), nn);
  for (i=0; i<Dim; i++)
    bgn[i] = dx[i] * double(nn[i]);
#endif

  // Set starting line as global index
  int startline = unstack_stack(0);
  // Since an extra line is added in the 2D output files for each distinct
  // value of y for gnuplot compatibility, account for those if Dim=2
  if (Dim == 2) {
    startline += nn[1];
  }

  // Skip to starting line
  for (i=0; i<startline; i++)
    fgets(tt, 256, inp);

  // Loop over all ML points and fill w
  for (i=0; i<ML; i++) {
    for (j=0; j<Dim; j++)
      fscanf(inp,"%lf ", &dm[j]);
    fscanf(inp, "%lf %lf\n", &dr, &di);
    if (i == 0) {
      // Print line and position info for starting line
      printf("Processor %d reading from line %d. x = [ ", myrank, i+1);
      for (j=0; j<Dim; j++)
        printf("%lf ", dm[j]);
      printf("]\n");
    }
    w[i] = dr + I * di;
  }

} // end read_one_resume_file

/**
 * @brief Reads the resume field data of other programs. The file must end in ".res"
 * 
 */
void read_resume_files() {
  FILE *inp;
  char nm[50];
  int iter_file_flag = 0, i , j;

  inp = fopen("wpl.res","r");
  if (inp!=NULL) {
    read_one_resume_file(inp, wpl);
    fclose(inp);
    printf("Read wpl.res\n");
  }

  inp = fopen("wabp.res","r");
  if (inp!=NULL) {
    read_one_resume_file(inp, wabp);
    fclose(inp);
    printf("Read wabp.res\n");
  }

  inp = fopen("wabm.res","r");
  if (inp!=NULL) {
    read_one_resume_file(inp, wabm);
    fclose(inp);
    printf("Read wabm.res\n");
  }
  
  
  inp = fopen("wbcp.res","r");
  if (inp!=NULL) {
    read_one_resume_file(inp, wbcp);
    fclose(inp);
    printf("Read wabp.res\n");
  }

  inp = fopen("wbcm.res","r");
  if (inp!=NULL) {
    read_one_resume_file(inp, wbcm);
    fclose(inp);
    printf("Read wabm.res\n");
  }
  
  
  inp = fopen("wacp.res","r");
  if (inp!=NULL) {
    read_one_resume_file(inp, wacp);
    fclose(inp);
    printf("Read wabp.res\n");
  }

  inp = fopen("wacm.res","r");
  if (inp!=NULL) {
    read_one_resume_file(inp, wacm);
    fclose(inp);
    printf("Read wabm.res\n");
  }
  
  
} // end read_resume_files

void write_fft_cpx(char* nm, fftw_complex *dt) {

  int i,j, nn[Dim];
  FILE *otp;
  otp = fopen(nm, "w");

  for (i=0; i<ML; i++) {
    unstack( unstack_stack(i) , nn );
    fprintf(otp,"%lf ", double(nn[0])*dx[0] )  ;
    for (j=1; j<Dim; j++)
      fprintf(otp,"%lf ", double(nn[j])*dx[j]);
    fprintf(otp,"%1.16e %1.16e\n", dt[i][0], dt[i][1] );

    if (Dim==2 && nn[0]==Nx[0]-1)
      fprintf(otp,"\n");

  }

  fclose(otp);
} // end write_fft_cpx

void write_kdata(char* nmi, complex<double> *dt) {

  int i,j, nn[Dim] ;
  double k2, kv[Dim];
  FILE *otp;
  char nm[50];
#ifndef PAR
  sprintf( nm , "%s.dat", nmi);
#else
  sprintf( nm , "%s.p%d.dat" , nmi, myrank);
#endif
  otp = fopen(nm, "w");

  for (i=0; i<ML; i++) {
    k2 = get_k(i, kv);
    fprintf(otp, "%lf ", kv[0]);
    for (j=1; j<Dim; j++)
      fprintf(otp, "%lf ", kv[j]);
    fprintf(otp,"%1.16e %1.16e %1.16e %1.16e\n", abs(dt[i]), sqrt(k2),
            real(dt[i]), imag(dt[i]) );

    if (Dim==2 && nn[0]==Nx[0]-1)
      fprintf(otp, "\n");
  }

  fclose(otp);
} // end write_kdata


void write_avg_kdata(char* nmi, complex<double> *dt) {

  int i,j, nn[Dim] ;
  double k2, kv[Dim];
  FILE *otp;
  char nm[50];
#ifndef PAR
  sprintf( nm , "%s.dat", nmi);
#else
  sprintf( nm , "%s.p%d.dat" , nmi, myrank);
#endif
  otp = fopen(nm, "w");

  for (i=0; i<ML; i++) {
    k2 = get_k( i , kv ) ;
    fprintf(otp,"%lf ", kv[0] ) ;
    for (j=1; j<Dim; j++)
      fprintf(otp,"%lf ", kv[j] ) ;
    fprintf(otp,"%1.16e %1.16e %1.16e %1.16e\n", abs(dt[i])/n_samples, sqrt(k2),
        real(dt[i])/n_samples, imag(dt[i])/n_samples );

    if (Dim==2 && nn[0]==Nx[0]-1)
      fprintf(otp,"\n");

  }

  fclose(otp);
} // end write_avg_kdata


void write_avg_data(char* nmi, complex<double> *dt) {

  int i,j, nn[Dim];
  FILE *otp;
  char nm[50];
#ifndef PAR
  sprintf( nm , "%s.dat", nmi);
#else
  sprintf( nm , "%s.p%d.dat" , nmi, myrank);
#endif
  otp = fopen(nm, "w");

  for (i=0; i<ML; i++) {
    unstack( unstack_stack(i) , nn );
    fprintf(otp,"%lf ", double(nn[0])*dx[0]  ) ;
    for (j=1; j<Dim; j++)
      fprintf(otp,"%lf ", double(nn[j])*dx[j]);
    fprintf(otp,"%1.16e %1.16e\n", real(dt[i])/n_samples, imag(dt[i])/n_samples );

    if (Dim==2 && nn[0]==Nx[0]-1)
      fprintf(otp,"\n");

  }

  fclose(otp);
} // end write_avg_data

void write_avg_data_bin(const char* nmi, complex<double> *dt) {

  int i,j, nn[Dim];
  FILE *otp;
  char nm[50];

  sprintf( nm , "%s.p%d.bin" , nmi, myrank);
  otp = fopen(nm, "wb");

  // The main cpu writes the number of processors, grid spacing, Nx, Ny, Nz
  if ( myrank == 0 ) {
    int dm = Dim ;
    fwrite( &dm , sizeof( int ) , 1 , otp ) ;
    fwrite( Nx , sizeof( int ) , dm , otp ) ;
    fwrite( L , sizeof( double ) , dm , otp ) ;
    fwrite( &nprocs , sizeof( int ) , 1 , otp ) ;
  }

  fwrite( &ML , sizeof( int ) , 1 , otp ) ;
  for ( i=0 ; i<ML ; i++ )
    tmp[i] = dt[i] / n_samples ;
  fwrite( tmp , sizeof( complex<double> ) , ML , otp ) ;

  fclose(otp);
} // end write_avg_data_bin

void write_data_bin(const char* nmi, complex<double> *dt) {

  int i,j, nn[Dim];
  FILE *otp;
  char nm[50];

  sprintf(nm, "%s.p%d.bin", nmi, myrank);
  otp = fopen(nm, "wb");

  // The main cpu writes the number of processors, grid spacing, Nx, Ny, Nz
  if (myrank == 0) {
    int dm = Dim;
    fwrite(&dm, sizeof(int), 1, otp);
    fwrite(Nx, sizeof(int), dm, otp);
    fwrite(L, sizeof(double), dm, otp);
    fwrite(&nprocs, sizeof(int), 1, otp);
  }

  fwrite(&ML, sizeof(int), 1, otp);
  fwrite(dt, sizeof(complex<double>), ML, otp);

  fclose(otp);
} // end write_data_bin

void write_data(char* nmi, complex<double> *dt) {

  int i,j, nn[Dim];
  FILE *otp;
  char nm[50];
#ifndef PAR
  sprintf(nm, "%s.dat", nmi);
#else
  sprintf(nm, "%s.p%d.dat", nmi, myrank);
#endif
  otp = fopen(nm, "w");

  for (i=0; i<ML; i++) {
    unstack(unstack_stack(i), nn);
    fprintf(otp, "%lf ", double(nn[0])*dx[0] );
    for (j=1; j<Dim; j++)
      fprintf(otp, "%lf ", double(nn[j])*dx[j] );
    fprintf(otp, "%1.16e %1.16e\n", real(dt[i]), imag(dt[i]) );

    if (Dim==2 && nn[0]==Nx[0]-1)
      fprintf(otp,"\n");
  }

  fclose(otp);
} // end write_data


void append_data(char* nmi, complex<double> *dt) {

  int i,j, nn[Dim];
  FILE *otp;
  char nm[50];
#ifndef PAR
  sprintf(nm, "%s.dat", nmi);
#else
  sprintf(nm, "%s.p%d.dat", nmi, myrank);
#endif

  otp = fopen(nm, "a");

  for (i=0; i<ML; i++) {
    unstack(unstack_stack(i), nn);
    for (j=0; j<Dim; j++)
      fprintf(otp, "%lf ", double(nn[j])*dx[j] );
    fprintf(otp, "%1.12e %1.12e\n", real(dt[i]), imag(dt[i]) );

    if (Dim==2 && nn[0]==Nx[0]-1)
      fprintf(otp,"\n");
  }

  fclose(otp);
} // end append_data

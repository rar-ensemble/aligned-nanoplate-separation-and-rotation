#ifndef CAVITY_HPP
#define CAVITY_HPP

#include "globals.h"
#include "array_utils.hpp"
#include "fft_wrappers.hpp"
#include <cmath>

class Cavity
{

  public:
    complex<double> *rho;
    complex<double> *rho_hat;

    // Cavity()
    // {
    //   array_utils::allocate_1d(&rho);
    //   array_utils::allocate_1d(&rho_hat);
    // };

    ~Cavity()
    {
      delete [] rho;
      delete [] rho_hat;
    };

    void allocate()
    {
      array_utils::allocate_1d(&rho);
      array_utils::allocate_1d(&rho_hat);
    };

    virtual void init_rho(void) = 0;
    virtual void write_bin(void) = 0;

};

class Channel: public Cavity
{

  private:
    int vert_dir;
    int hor_dir;
    double vert_center;
    double hor_center;
    double channel_width;
    double wall_width;
    double xi;

  public:
    Channel(double cw, double ww, double xi, int hd, int vd)
      : channel_width(cw), wall_width(ww), xi(xi), hor_dir(hd), vert_dir(vd)
    {
      // Initialize center positions
      vert_center = L[vert_dir] / 2.0;
      hor_center = L[hor_dir] / 2.0;
      if (myrank==0)
      {
        cout << "Initializing channel wall with channel width " << channel_width
          << ", wall width " << wall_width << ", and xi " << xi << endl;
        fflush(stdout);
      }
    };

    void init_rho(void);
    void write_bin(void) { write_data_bin("channel", rho); };

};

#endif // CAVITY_HPP

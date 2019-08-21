#include "cavity.hpp"

void Channel::init_rho()
{

  for (int i=0; i<ML; i++)
  {
    int i_global = unstack_stack(i);
    int nn[Dim];

    // Fill nn with global x, y [, and z] indices
    unstack(i_global, nn);

    // Get position in horizontal and vertical dimensions
    double x_hor = dx[hor_dir] * double(nn[hor_dir]);
    double x_vert = dx[vert_dir] * double(nn[vert_dir]);

    // Multiply horizontal and vertical functions to get wall/channel shape
    // The maximum value will be 1
    double hor_arg = (abs(x_hor-hor_center) - wall_width/2.0);
    double vert_arg = (abs(x_vert-vert_center) - channel_width/2.0);
    rho[i] = 0.5 * erfc(hor_arg / xi);
    rho[i] *= 0.5 * (1 + erf(vert_arg / xi));
  }

  // Initialize rho_hat (fourier transform of rho)
  fft_fwd_wrapper(rho, rho_hat);

  // Write data to file
  write_bin();

}

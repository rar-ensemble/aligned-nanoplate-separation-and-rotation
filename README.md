# Aligned Nanoplate-separation-and-rotation
Data repository for "Experiments and Simulations Probing Local Domain Bulge and
String Assembly of Aligned Nanoplates in a Lamellar Diblock Copolymer"

The non-code files are covered by the Creative Commons Attribution Share Alike
4.0 International license, and the code files are covered by  GNU General Public
License v3.0.

The raw results for the different simulation types (angle rotation, particle
separation) in the paper as well as the initialization files are provided. The
code is included in the `code` directory. To compile the code, `cd` into the
`code` directory and use `make`.

In order to compile, you must have a MPI compiler installed. For our code, we use
gcc compiler with the Open MPI standard. Additionally, our code requires the use
of a Fourier Transform library. We used FFTW 3.3.4 installed locally (in the
user's home directory) instead of a global installation. `${FFTWHOME}` refers to
the directory containing the local installation. Do not forget when installing
FFTW to include the `--enable-mpi` flag during configuration.

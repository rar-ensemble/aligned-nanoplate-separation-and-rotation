My Main Page                         {#mainpage}
============

This software is the cumulative work of Robert Riggleman, Jason Koski, Ben
Lindsay and Christian Tabedzki. This is by no means the most
complete documentation of the code but it is the most complete that exists.

This code has remnants of other code in it but it attempts to model the
interactions between an A/B diblock polymer and a C-type grafted nanoparticle,
allowing for any of the three possible \f$\chi\f$ combinations (AB/BC/AC) to be
either negative, positive or zero. 

# Grafted-NPS

Discrete Polymer Field Theory code that simulates grafted or bare, free or fixed
nanoparticles in a diblock copolymer or homopolymer matrix.

# Tips 
## Units/Smearing
This code is a _b_-scaled unit, which is the monomeric unit of the chain.
Generally, the suggested maximum smearing length for particles is 0.2 *Rg*
units. Therefore, you will need to convert the units to *b* units before
inputting the value into the input file. If you know the length of *Rg* and *b*
in terms of nanometers, then obtaining the maximum smearing is:

\f[0.2R_g\left(\frac{X \, \mbox{nm}}{1 R_g}\right)\left(\frac{1 b}{Y \,
R_g\mbox{nm}}\right)=Z\,\, b \f]

## ChiN
In this code, there is no need to specify the value of the \f$\chi\f$N. In order
to determine the value of chi for the bcp.input:
1. Calculate \f$\chi\f$<SUB>Balsara</SUB>, the statistical value of chi from the
   textbook.
2. Calculate N<SUB>Balsara</SUB> and N<SUB>Simulation</SUB>. The ratio between
   these two will be used to determine \f$\chi\f$<SUB>Simulation</SUB> from
   \f$\chi\f$<SUB>Balsara</SUB>. In future editions of the code, it might be
   necessary to specify which value of N should be used for \f$\chi\f$N.

## Xi Calculation
To determine the size of your resolution, you should look at the size of
\f$\xi\f$ (the value of how smeared out the particle is) and also the size of
smearing length. These three should be very similar to each other in magnitude.

As a refresher, the equation used to model a nanoparticle from 
<a href="http://dx.doi.org/10.1063/1.4853755">Jason's paper</a> 
\f[\Gamma ( \vert \mathbf{r} - \mathbf{r}'\vert, \mathbf{u} ) =
 \frac{\rho_0}{4}
 \mbox{erfc}
 \left[
 \frac{\vert \mathbf{u} \cdot (\mathbf{r} - \mathbf{r}'r) \vert - L/2 }
 {\xi}
 \right]
 \times
 \mbox{erfc}
 \left[\frac{\vert \mathbf{u} \times (\mathbf{r}
- \mathbf{r}'r) \vert - L/2 }
 {\xi}
 \right]
 \f]
  In this equation, we want to determine what effect \f$\xi\f$ has on the
  particle. The value of \f$\xi\f$ should be   **at least 1/6** of the particle
  length. Larger values (such as 1/2) will cause the system to have a cusp in
  the center of the particle. Although a smaller \f$\xi\f$ will cause the
  particle to have sharper edges, the tradeoff is that this will increase the
  runtime of the simulation since the resolution of the system will have to be
  increased to match \f$\xi\f$.

  Compare the following values of \f$\xi\f$ and the way they affect the particle
  shape:
1. <a href="http://www.wolframalpha.com/input/?i=y%3Derfc((norm((-1,-2).(x,0))+-+6%2F2)%2F(1%2F2*6))++from+-6+to+6">1/2</a>
2. <a href="http://www.wolframalpha.com/input/?i=y%3Derfc((norm((-1,-2).(x,0))+-+6%2F2)%2F(1%2F4*6))++from+-6+to+6">1/4</a>
3. <a href="http://www.wolframalpha.com/input/?i=y%3Derfc((norm((-1,-2).(x,0))+-+6%2F2)%2F(1%2F6*6))++from+-6+to+6">1/6</a>
4. <a href="http://www.wolframalpha.com/input/?i=y%3Derfc((norm((-1,-2).(x,0))+-+6%2F2)%2F(1%2F12*6))++from+-6+to+6">1/12</a>



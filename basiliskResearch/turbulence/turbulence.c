/**
# Decaying two-dimensional turbulence

We solve the two-dimensional incompressible Euler equations using a
vorticity--streamfunction formulation. */

#include "grid/multigrid.h"
#include "navier-stokes/stream.h"



/**
The domain is square of size unity by default. The resolution is
constant at $256^2$. */

int main() {
	size(4);
  init_grid (256);
  origin(-2,-2);
  run();
}

/**
The initial condition for vorticity is just a white noise in the range
$[-1:1]$ .*/

event init (i = 0) {
  foreach()
    omega[] = (x*x + y*y <1)* noise();
}

/**
We generate images of the vorticity field every 4 timesteps up to
$t=1000$. We fix the colorscale to $[-0.3:0.3]$.

![Evolution of the vorticity](turbulence/omega.mp4)(autoplay loop) */

event output (i += 4; t <= 1000) {
  output_ppm (omega, min = -0.3, max = 0.3, file = "omega.mp4");
}

event end(t=1000)
{
	printf(" i = %d , t = %g\n", &i,&t);
}

event logfile(i++)
{
	if(i==0){
		 fprintf (ferr,
	     "t dt perf.t\n");
  }
  fprintf (ferr,
	   "%g %g %.2e\n", 
	   t, dt, perf.t);
  fflush (ferr);
}
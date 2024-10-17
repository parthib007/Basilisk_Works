/**
# Spherical droplet blown by stream

Axisymetrical case, where a single drop is drag by an uniform stream of another.
After some time-step, a non physical velocity arise in the upstream apex of the droplet.
*/

#include "axi.h"

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"

/**
Fluids properties are defined next.
*/

//Drop properties
#define RHOL 824.
#define MUL 2.17e-3
#define SIGMA 2.e-2

//Stream properties
#define RHOG 1.1856
#define MUG 1.84e-5
#define DIAMETER 2.4e-3 
#define USTREAM 9.9198

#define MAXTIME 15.e-3


//Default refinements
int minlevel = 6;
int MAXLEVEL = 12;

double maxruntime = HUGE;

/**
Inlet and outlet boundary conditions. */

u.n[left]  = dirichlet(USTREAM);
u.t[left]  = dirichlet(0);
u.n[right] = neumann(0);
p[left]    = neumann(0);
p[right]   = dirichlet(0);

/**
MAXLEVEL can be overriden by input parameter. */

int main (int argc, char * argv[]) {
  if (argc > 1)
    MAXLEVEL = atoi (argv[1]);
  
  L0 = 10.667*DIAMETER;
  size (L0);
  //origin(L0/5,0); // origin was negative
  init_grid (pow(2.0,minlevel));

  // CFL number
  CFL = 0.4;

  f.sigma = SIGMA;
  rho1 = RHOL;
  rho2 = RHOG;

  mu1 = MUL;
  mu2 = MUG;
 
  TOLERANCE = 1e-5; //default is 1e-3.

  run();
}

/**
Initial condition is given by drop position only, there is no flow at
the begining.*/

event init (t = 0) {

	/* Doubt 
  foreach()
    f[] = 0;
	*/
    

  if (!restore (file = "dump")){

    refine (sq(x-DIAMETER) + sq(y) - sq(0.6*DIAMETER) < 0 &&
       sq(x-DIAMETER) + sq(y) - sq(0.4*DIAMETER) > 0 && level < MAXLEVEL);

    fraction (f, sq(0.5*DIAMETER) - (sq(x-DIAMETER) + sq(y)));
	
  }
}

event end (t = MAXTIME) {
  printf ("i = %d t = %g\n", i, t);
}

/**
Adaptation tolerances are defined as a function of the RMS velocity. */

event adapt (i++) {
  double uRMS = 0;
  foreach()
    foreach_dimension()
      uRMS += sq(u.x[]);

  double uemax = 2e-4*sqrt(uRMS);
  adapt_wavelet ({f,u}, (double[]){1e-4,uemax,uemax,uemax}, MAXLEVEL);
}

/**

## Post-processing

Log file reports the droplet position and velocity, besides some numerical
parameters. gfsview snapshots are each 100 time-steps.*/

event logfile (i += 10,first) {
  double xd = 0., yd = 0., sd = 0.;
  double vdx = 0., vdy = 0.;
  if (i == 0){
    fprintf (ferr,
	     "i dt t sd xd/sd yd/sd vdx/sd, vdy/sd"
	     " mgp.i mgpf.i mgu.i grid->tn perf.t perf.speed\n");
  }

  foreach(reduction(+:xd) reduction(+:yd)
	  reduction(+:vdx) reduction(+:vdy) reduction(+:sd)) {
    double dv = f[]*dv();
    xd += x*dv;
    yd += y*dv;
    vdx += u.x[]*dv;
    vdy += u.y[]*dv;
    sd += dv;
  }

  fprintf (ferr,
	   "%d %.4e %.4e %.4e %.2e %.2e %.2e %.2e %d %d %d %ld %.2e %.2e\n", 
	   i, dt, t, sd, xd/sd, yd/sd,
	   vdx/sd, vdy/sd, mgp.i, mgpf.i, mgu.i,
	   grid->tn, perf.t, perf.speed);

  fflush (ferr);
}

event snapshot (i = 0; i+=500; t <= MAXTIME)
{
  if(i > 0)
    dump (file = "dump");
  
  char name[80];
  sprintf (name, "snapshot-%d.png", i);
  output_ppm (f, file = name, spread = -1, linear = true, map = cool_warm, n = 512, box = {{0,0},{L0,L0/10}});
}
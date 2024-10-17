/**
# Droplet blow dy stream, bag mode secondary atomisation

We wish to study the behaviour of a single drop in free fall
affected dy a uniform flow of air, as in the [Opfer 2014](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.89.013023).

We use the centered Navier--Stokes solver. There are two phases, air and the liquid. The grid adaptation has been
modyfied to support region-depending maximum refinement (maxlevel).*/


#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "view.h"


/**
Physical properties are defined based on [Opfer 2014](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.89.013023). */

/**Liquid properties*/
#define RHOL 824.
#define MUL 2.17e-3
#define SIGMA 2.e-2

/**Air properties*/
#define RHOG 1.1856
#define MUG 1.84e-5

/**Drop diameter and stream velocity*/
#define DIAMETER 1.98e-4
#define USTREAM 67.8282

#define xc 3*DIAMETER
#define yc L0/2
#define zc L0/2

#define ACCEL 9.81      //On x axis

#define tsh 7.9657e-5    // shear time depending on We

#define MAXTIME 15.e-3

/**
Minimum refinement level will be used directly, maxlevel will only apply
on a certain region of the domain. */

int minlevel = 6;   // init grid of 32 points
int maxlevel = 12;  //11

double maxruntime = HUGE;


/**
Boundary conditions are set only at left (inlet) and right (outlet)
patches. Inlet is supposed to be uniform */

u.n[left]  = dirichlet(USTREAM);
u.t[left]  = dirichlet(0);
u.n[right] = neumann(0);

p[left]    = neumann(0);
p[right]   = dirichlet(0);
/**
The main function can take the maximum refinement level as argument.
The main program will run ten steps with both versions of AMR.*/

int main (int argc, char * argv[]) {
  if (argc > 1)
    maxlevel = atoi (argv[1]);
  
  /**
  We set the domain geometry and initial refinement.
  The drop will be centered at the origin, slightly
  moved to the left side*/
  L0 = 16* DIAMETER;
  size (L0);
  origin(0,0,0);
  init_grid (pow(2.0,minlevel));



  // CFL number
  CFL = 0.25;

  /**
  Physical parameters are set by the macros at the beggining.*/
  
  f.sigma = SIGMA;
  rho1 = RHOL;
  rho2 = RHOG;

  mu1 = MUL;
  mu2 = MUG;

  /**
  Tolerance of the Poisson problem should never be reduced,
  this could lead to non divergence free velocity fields
  and mass conservation errors. */
  
  TOLERANCE = 1e-4; 

 run();

  

 
}

/**
Initial condition is given by drop position only, there is no flow at the begining.
First step should give a velocity field close to potential flow solution around the drop,
plus small velocities inside the droplet.*/

event init (t = 0) {
  if (!restore (file = "restart")){
    // mask(y > 5*DIAMETER? top:none);
    //refine (sq(x-DIAMETER) + sq(y) - sq(0.6*DIAMETER) < 0 && sq(x-DIAMETER) + sq(y) - sq(0.4*DIAMETER) > 0 && level < maxlevel);
	#if TREE
	refine (sq(x-xc) + sq(y-yc) + sq(z-zc) - sq(0.6*DIAMETER) < 0 && sq(x-xc) + sq(y-yc) + sq(z-zc) - sq(0.4*DIAMETER) > 0 && level < maxlevel);
	#endif
        fraction (f, sq(0.5*DIAMETER) - (sq(x-xc) + sq(y-yc) + sq(z-zc)));
//    output_ppm(f,linear = true, n=512, box = {{0,0},{L0,L0/2}}, file = "vofinit.png");
	}
  }

/**
Only ten steps are performed, just to check mesh adaptation. */

event end (i=5) {
  printf ("i = %d t = %g\n", i, t);
}


/**
Every ten timesteps, we output the time, timestep, multigrid iterations, CPU time, speed and amount of cells. */

/*
event diameter_data (i++) {
  scalar pos_x[], pos_y[];
  position (f, pos_y, {0,1}); //y-axis
  position (f, pos_x, {1,0}); // along stream
  FILE *data = fopen( "diameter" , "a");
  fprintf (data, "%g %g %g\n", t, 2.*statsf(pos_y).max,(statsf(pos_x).max - statsf(pos_x).min));
  fclose(data);
}
*/

event logfile (i ++,first) {
   if (i == 0){
    fprintf (ferr, "i t dt mgp.i mgpf.i mgu.i perf.t perf.speed grid->tn\n");
  }
  fprintf (ferr, "%d %g %g %d %d %d %.2e %.2e %ld \n", i, t, dt, mgp.i, mgpf.i, mgu.i, perf.t, perf.speed, grid->tn);
  fflush (ferr);
}

/*
event interface (i += 100; t <= MAXTIME)
{
  char name[80];
  sprintf(name, "interface-%d", i);
  FILE *file = fopen (name , "a");
  output_facets(f, file);
  fclose(file);
}
*/


/*
event snapshot(i += 5; t <= MAXTIME)
{
  
  char name[80];
  sprintf (name, "snapshot-f-%i_ts.ppm", i);
  output_ppm ( f ,linear = true,n = 800, box = {{0,0},{L0,L0}}, file = name );

}
*/

/*

event movie (i += 50; t <= MAXTIME) {
  output_ppm (f, linear = true, n = 800, box = {{0,0},{L0,L0}}, file = "vof.mp4", spread = -1, map = cool_warm);
  output_ppm ( u.x, box = {{0,0},{L0,L0/2}}, n = 800, file = "speed.mp4" , min = 0, max = USTREAM, map = cool_warm);
  //scalar l[];
  //foreach()
    //l[] = level;
//  output_ppm (l, min = minlevel, max = maxlevel,box = {{0,0},{L0,L0/2}}, file = "level.mp4");
  //output_ppm (l, n = 800, box = {{0,0},{L0,L0/2}}, file = "level.mp4", min = 0, max = 1, linear = false);
}
*/

event viewing (i ++) {
  view (camera = "iso", fov = 44, tx = -0.418, ty = 0.288, width = 1600, height = 1200);
  clear();
  draw_vof ("f");
  char name[80];
  sprintf(name , "grid-%d.png", i);
  save (name);
}

#if TREE
event adapt (i++) {
  double uemax = 1e-2;
  
    adapt_wavelet ({f,u}, (double[]){1e-3,uemax,uemax,uemax}, maxlevel);
}
#endif

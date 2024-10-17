#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"


/*
#define RHOL 785.
#define MUL 16.1e-3
#define SIGMA 48.3e-3

#define RHOG 1.
#define MUG 1.81e-5

#define DIAMETER 2.7e-3 
#define USTREAM 15.

#define MAXTIME 30e-3


int minlevel = 6;   // init grid of 32 points
int MAXLEVEL = 11;  //12

double maxruntime = HUGE;


u.n[left]  = dirichlet(USTREAM);
u.t[left]  = dirichlet(0);
u.n[right] = neumann(0);

u.r[left] = dirichlet(0);

p[left]    = neumann(0);
p[right]   = dirichlet(0);
*/
#define DIAMETER 2.7e-3 

int main () {

  L0 = 5*DIAMETER;
  //size (L0);
  //origin(0,0,0);
  //init_grid (pow(2.0,minlevel));

  N = 64;
  
  /*
  CFL = 0.4;

  
  f.sigma = SIGMA;
  rho1 = RHOL;
  rho2 = RHOG;

  mu1 = MUL;
  mu2 = MUG;

 
  
  TOLERANCE = 1e-4; 
  */

  run();

 
}


event init (t = 0) {

    refine (sq(x) + sq(y) + sq(z) - sq(0.6*DIAMETER) < 0 && sq(x) + sq(y) + sq(z) - sq(0.4*DIAMETER) > 0 && level < 11);
    fraction (f, sq(0.5*DIAMETER) - (sq(x) + sq(y) + sq(z)));

}



event end (i = 10) {
  printf ("i = %d t = %g\n", i, t);
}


event adapt (i++) {
  double uemax = 5e-2;
    adapt_wavelet ({f,u}, (double[]){1e-3,uemax,uemax,uemax}, 11, 6);
}


event logfile (i ++,first) {
  if (i == 0){
    fprintf (ferr,
	     "dt t mgp.i mgpf.i mgu.i perf.t perf.speed grid->tn\n");
  }

  fprintf (ferr,
	   "%g %g %d %d %d %.2e %.2e %ld \n", 
	   dt, t, mgp.i, mgpf.i, mgu.i, perf.t, perf.speed, grid->tn);
  fflush (ferr);
}


event movies (i =0; i++; i<=10)
{  
  view (camera = "iso",
	fov = 14.5, tx = -0.418, ty = 0.288,
	width = 1600, height = 1200);
  clear();
  draw_vof("f");
  save("movie.mp4");
}
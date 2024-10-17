#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "adapt_wavelet_limited.h"
#include "view.h"


/* Parameters of the problem */

/* density ratio */
#define RHO_R 695.007

/* viscosity ratio */
#define MU_R  117.935

/* Reynolds number */
#define Re 1,534.032

/* Weber number */
#define We 14.

/* Surface tension */
#define SIGMA 2.e-2

/* radius of the drop */
//#define R 0.0012
#define DIAMETER 2.4e-3

/* We define a square domain of length 5*D */
#define L 10*DIAMETER

/* We define speed of air-flow */
#define USTREAM 9.9198

/* Initial location of the drop */
#define xc DIAMETER
#define yc 0.

#define TEND 15.e-3

/* adaptive mesh refinement levels */
#define LEVEL_MAX 11
#define LEVEL_MIN 5

/* error criteria for amr */
#define femax 1e-3
//#define uemax 1e-2

/* Tolerance for the poisson solver */
#define TOL 1e-4

/* inlet */  
u.n[left] = dirichlet(USTREAM); 

/* outflow bc at the outlet */
u.n[right] = neumann(0);
p[right]   = dirichlet(0);
pf[right]  = dirichlet(0.);

int main(){
  
  L0 = L;
  N = 1<<LEVEL_MIN;
  size (L0);
  init_grid (N);
  
  
  /* outer fluid */
  rho2 = 1.;
  mu2 = 1/Re;

  /* drop */
  rho1 = RHO_R;
  mu1 = MU_R/Re;
  
  f.sigma = 1/We;

  CFL = 0.2;
  
  TOLERANCE = TOL;
  run();

}

event init (t = 0) {
    refine (sq(x-xc) + sq(y) - sq(0.6*DIAMETER) < 0 && sq(x-xc) + sq(y) - sq(0.4*DIAMETER) > 0 && level < LEVEL_MAX);
    fraction (f, sq(0.5*DIAMETER) - (sq(x-xc) + sq(y)));
    //output_ppm(f,linear = true, n=512, box = {{0,0},{L0,L0/2}}, file = "vofinit.png");
  }
  
event end (i=3) {
  printf ("i = %d t = %g\n", i, t);
}

int refRegion(double x,double y, double z){
    int lev;
    if( sq(x-DIAMETER)+sq(y) < sq(DIAMETER*0.7) )
      lev = LEVEL_MAX;
    else
      lev = LEVEL_MIN;
    return lev;
}

event adapt (i++) {
	double uemax = 0.01*normf(u.x).avg;
	adapt_wavelet ({f,u}, (double[]){1e-4,uemax,uemax}, maxlevel = LEVEL_MAX);
}

event logfile (i ++,first) {
   if (i == 0){
    fprintf (ferr, "i dt t mgp.i mgpf.i mgu.i perf.t perf.speed grid->tn\n");
  }
  fprintf (ferr, "%d %g %g %d %d %d %.2e %.2e %ld \n", i, dt, t, mgp.i, mgpf.i, mgu.i, perf.t, perf.speed, grid->tn);
  fflush (ferr);
}

event snapshot(i++)
{
  
  char name[80];
  sprintf (name, "snapshot_f-%i.png", i);
  output_ppm ( f ,linear = true,n = 800, box = {{0,0},{L0,L0/2}}, file = name );

}

event viewing (i++) {
  view (width = 1024, height = 1024, tx = -0.5, ty = -0.5, fov =20);
  clear();
  draw_vof ("f");
  box (notics = true);
  cells();
  char name[80];
  sprintf(name , "grid-%d.png", i);
  save (name);
}



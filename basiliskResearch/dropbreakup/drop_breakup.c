/* Dynamics of drop-breakup in a cross flow */
//#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "tag.h"
#include "view.h"

/* Parameters of the problem */

/* density ratio */
#define RHO_R 656.7

/* viscosity ratio */
#define MU_R  69.9

/* Reynolds number */
#define Re 1420.2

/* Weber number */
#define We 11.5

/* radius of the drop */
#define R 0.5

/* We define a square domain of length 5*D */
#define L 8.

/* Initial location of the drop */
#define xc 1.
#define yc 0.

#define TEND 5.

/* adaptive mesh refinement levels */
#define LEVEL_MAX 12
#define LEVEL_MIN 8

/* error criteria for amr */
#define femax 1.e-3
#define uemax 1.e-2

/* Tolerance for the poisson solver */
#define TOL 1.e-6

/* inlet */  
u.n[left] = dirichlet(1.);     // shouldn't it be negative?

/* outflow bc at the outlet */
u.n[right] = neumann(0);
p[right]   = dirichlet(0);
pf[right]  = dirichlet(0.);

int main(){
  
  L0 = L;
  N = 1<<LEVEL_MIN;
  
  
  /* outer fluid */
  rho2 = 1.;
  mu2 = 1./Re;

  /* drop */
  rho1 = RHO_R;
  mu1 = MU_R/Re;
  
  f.sigma = 1./We;

  TOLERANCE = TOL;
  run();

}

//event init (t = 0) {
//  restore (file = "dump");


event init_drop (t = 0) {


  double  dR = 0.1*R;
  refine (level < LEVEL_MAX &&
	  sq(R + dR) - sq(x-xc) - sq(y-yc) > 0.0);

  fraction (f, sq(R)-sq(x-xc)-sq(y-yc));

}

event logfile(i++){
	if (i == 0){
    fprintf (ferr,
	     "i dt t mgp.i mgpf.i mgu.i perf.t perf.speed grid->tn\n");
  }
  
  fprintf (ferr,
	   "%d %g %g %d %d %d %.2e %.2e %ld \n", 
	    i,dt, t, mgp.i, mgpf.i, mgu.i, perf.t, perf.speed, grid->tn);
  fflush (ferr);
}

event adapt (i++) {

  scalar ff[];
  foreach()
    ff[] = f[];
  
  boundary({ff});
  adapt_wavelet ({f,ff, u},
		 (double[]){femax,femax,uemax,uemax},LEVEL_MAX);
  
}

/**
event log_stats (i +=10) {

  char fname[100];

  double vol_i, xi, yi, ui, vi; 
  vol_i = xi = yi = ui = vi = 0.;

  sprintf (fname, "drop_stats_R-%d.dat", LEVEL_MAX);
  static FILE * fp = fopen (fname, "a");
   
   foreach(reduction(+:vol_i) reduction(+:xi) reduction(+:yi)
	   reduction(+:ui) reduction(+:vi)){
     double c = f[]*dv();
     vol_i += c;
     xi += x*c;
     yi += y*c;
     ui += u.x[]*c;
     vi += u.y[]*c;
   }
   
   xi = xi/vol_i; yi = yi/vol_i; ui = ui/vol_i; vi = vi/vol_i;
   fprintf(fp, "%.10f %.10f %.10f %.10f %.10f %.10f \n", t,
	   vol_i, xi, yi, ui, vi);
   
}
*/


event movies (i += 4; t <= 15.)
{
  output_ppm (f, file = "f.mp4", box = {{0,0},{8.0,8.0}},
	      linear = true,n = 512, min = 0, max = 1);
}


event viewing (t += 0.01) {
  view (tx = -0.5, ty = -0.5);
  box(notics=true);
  draw_vof ("f", lw = 2, filled = 1, fc = {0.5, 0., 0.15});
  //cells();
  save ("movie.mp4");
} 

/**
event plotInterface (t++; t<=TEND) {

  char name[80];
  sprintf (name, "interface-%f.txt", t);
  FILE* fp = fopen (name, "w");

  output_facets (f, fp);
}
*/

/**
event dump_output (t = 0; t ++; t <= TEND)
{

  char name[100];
  sprintf (name, "snapshot_t-%03d", (int) t);
  scalar pid[], ff[];
  foreach(){
     pid[] = fmod(pid()*(npe() + 37), npe());
     ff[] = f[] < 1.0e-6 ? 0 : f[] > 1. - 1.0e-6 ? 1. : f[];
  }
  boundary ({pid, ff});
  dump(file="dump");
  dump (file = name);

  sprintf (name, "image_t-%03d.png", (int) t);
//view (fov = 6.43549, quat = {0,0,0,1}, tx = -0.497361, ty = -0.135413, bg = {1,1,1}, width = 844, height = 278,);
	view (tx = -0.5, ty = -0.5);
  box(notics=true);
  draw_vof ("f", lw = 2, filled = 1, fc = {0.5, 0., 0.15});
  //cells();
  save(name);
  
  
  sprintf (name, "interface_t-%f.txt", t);
 

  FILE* fp = fopen (name,"w");

  output_facets (f, fp);

}*/


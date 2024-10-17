/**
# Droplet blow dy stream, bag mode secondary atomisation

We wish to study the behaviour of a single drop in free fall
affected dy a uniform flow of air, as in the [Opfer 2014](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.89.013023).

We use the centered Navier--Stokes solver. There are two phases, air and the liquid. The grid adaptation has been
modyfied to support region-depending maximum refinement (maxlevel).*/

#include "axi.h"
#include "navier-stokes/centered.h"
//#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
//#include "adapt_wavelet_limited.h"
#include "view.h"

/**
Physical properties are defined based on [Opfer 2014](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.89.013023). */

/* Parameters of the problem */

/* density ratio */
#define RHO_R 695.007

/* viscosity ratio */
#define MU_R  117.935

/* Reynolds number */
#define Re 9040

/* Weber number */
#define We 14

/* radius of the drop */
#define R 1.2e-3 
#define DIAMETER 2.4e-3 
#define USTREAM 9.9198

#define ACCEL 9.81      //On x axis

#define tsh 6.378e-3    // shear time depending on We

#define MAXTIME 15.e-3
//#define MAXITER 5000

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
//u.t[left]  = dirichlet(0);
u.n[right] = neumann(0);

//p[left]    = neumann(0);
p[right]   = dirichlet(0);
pf[right]  = dirichlet(0);
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
  L0 = 10* DIAMETER;
  size (L0);
  origin(-L0/5,0);
  init_grid (pow(2.0,minlevel));
  
    // CFL number
  CFL = 0.25;
  
  /**
  Physical parameters are set by the macros at the beggining.*/
  
  /* outer fluid */
  rho2 = 1.;
  mu2 = 1./Re;

  /* drop */
  rho1 = RHO_R;
  mu1 = MU_R/Re;
  
  f.sigma = 1./We;
  TOLERANCE = 1e-4; 

  run();


 
}

/**
Initial condition is given by drop position only, there is no flow at the begining.
First step should give a velocity field close to potential flow solution around the drop,
plus small velocities inside the droplet.*/

event init (t = 0) {
 foreach()
    f[] = 0;
  if (!restore (file = "restart")){
//    mask(y > 10*DIAMETER? top:none);
    refine (sq(x) + sq(y) - sq(0.6*DIAMETER) < 0 && sq(x) + sq(y) - sq(0.4*DIAMETER) > 0 && level < maxlevel);
    fraction (f, sq(0.5*DIAMETER) - (sq(x) + sq(y)));
//    output_ppm(f,linear = true, n=800, box = {{0,0},{L0,L0/2}}, file = "vofinit.png");
	}
  }

/**
Only ten steps are performed, just to check mesh adaptation. */

event end (t=MAXTIME) {
  printf ("i = %d t = %g\n", i, t);
}

/*
event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] = ACCEL;
}
*/

/**
## Region limited mesh adaptation

Any function with input(x,y,z) returning an integer can be used to
define maxlevel locally for certain regions. Next, a simple example
using a circle (only for 2D case) slightly bigger than the drop is used.
Inside that region, AMR can refine until maxlevel, outside
only maxlevel-2 cells are allowed. Then, the alternative function
adapt_wavelet_limited is employed. */

int refRegion(double x,double y, double z){
    int lev;
    if( sq(x)+sq(y) < sq(DIAMETER*0.7) )
      lev = maxlevel;
    else
      lev = maxlevel-2;

    return lev;
}

/*
event adapt (i++) {
  double uemax = 1e-2;
  
//  Choose between region limited adaptation or standard adaptation:

 if(limitedAdaptation)
    adapt_wavelet_limited ({f,u}, (double[]){1e-3,uemax,uemax,uemax}, refRegion, minlevel);
  else
    adapt_wavelet ({f,u}, (double[]){1e-3,uemax,uemax,uemax}, maxlevel, minlevel);
}

*/


event adapt (i++){
  double uRMS = 0;
  foreach()
    foreach_dimension()
      uRMS += sq(u.x[]);
//  double uemax = 2e-4*sqrt(uRMS);
  double uemax = 1e-2;
  
 
    scalar ff[];
    foreach()
      ff[] = f[];
    boundary({ff});

    adapt_wavelet((scalar *){f,ff,u}, (double []){1e-2,1e-2,uemax,uemax,uemax}, maxlevel);
//    adapt_wavelet((scalar *){f,u}, (double []){0.001, 0.001, 0.001}, maxlevel);
}




event diameter_data (i++) {
  scalar pos_x[], pos_y[];
  position (f, pos_y, {0,1}); //y-axis
  position (f, pos_x, {1,0}); // along stream
  FILE *data = fopen( "diameter" , "a");
  fprintf (data, "%g %g %g %g\n", t/tsh, (2.*statsf(pos_y).max)/DIAMETER,((statsf(pos_x).max - statsf(pos_x).min))/DIAMETER,0.5*((statsf(pos_x).max - statsf(pos_x).min)) );
  fclose(data);
}

/**
Every ten timesteps, we output the time, timestep, multigrid iterations, CPU time, speed and amount of cells. */
event logfile (i ++,first) {
   if (i == 0){
    fprintf (ferr, "i dt t mgp.i mgpf.i mgu.i perf.t perf.speed grid->tn\n");
  }
  fprintf (ferr, "%d %g %g %d %d %d %.2e %.2e %ld \n", i, dt, t, mgp.i, mgpf.i, mgu.i, perf.t, perf.speed, grid->tn);
  fflush (ferr);
}

event interface (i += 1000; t <= MAXTIME)
{
  char name[80];
  sprintf(name, "interface-%d", i);
  FILE *file = fopen (name , "a");
  output_facets(f, file);
  fclose(file);
}

event urms_data(i++)
{
	FILE *file = fopen("urms" , "a");
	double uRMS = 0;
	foreach()
		foreach_dimension()
      uRMS += sq(u.x[]);
	double urms = sqrt(uRMS/N);
	fprintf(file, "%d %g\n",i, urms);
	fclose(file);
}

event snapshot(i += 1000; t <= MAXTIME)
{
  
  char name[80];
  sprintf (name, "snapshot-f-%i_ts.ppm", i);
//  output_ppm ( f ,linear = true,n = 800, box = {{-L0/5,0},{L0,L0/2}}, file = name );
  output_ppm ( f ,linear = true,n = 800, file = name );

}


event movie (i += 10; t <= MAXTIME) {
  output_ppm (f, linear = true, n = 800, box = {{0,0},{L0,L0/2}}, file = "vof.mp4", spread = -1, map = cool_warm);
 
}


event viewing (i += 500; t <= MAXTIME) {
  //view (width = 400, height = 400, fov = 20, ty = -0.5, theta =1.57079632679, phi =0, psi =0);
//  view (tx = 0., ty = 0.,fov=20);
//  view (width = 1024, height = 1024, tx = -0.5, ty = -0.5, fov =20);
  view (width = 1024, height = 1024, tx = -L0/5 - 0.5, ty = -0.5);
  clear();
  draw_vof ("f");
  box (notics = true);
  cells();
  mirror ({0,1}) {
   squares ("u.x", spread = -1, linear = true, map = cool_warm);
   draw_vof ("f");
  }
  char name[80];
  sprintf(name , "grid-%d.png", i);
  save (name);
}


/*
event pictures (t = end)
{
  char name[80];
#if 0  
  sprintf (name, "dump-test");
  dump (name);
#endif
	scalar l[];
	foreach()
	 	l[] = level;  
//  view (width = 600, height = 600);  
  squares ("l", linear = true, map = cool_warm);
  draw_vof ("f");
  
  sprintf (name,"final.png");
  save (name);
}
*/

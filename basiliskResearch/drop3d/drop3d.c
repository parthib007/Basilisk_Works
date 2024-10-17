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

#define ACCEL 9.81      //On x axis

#define tsh 7.9657e-5    // shear time depending on We

#define MAXTIME 15.e-3

/**
Minimum refinement level will be used directly, maxlevel will only apply
on a certain region of the domain. */

int minlevel = 5;//6   // init grid of 64 points
int maxlevel = 8;//10  //11

double maxruntime = HUGE;


/**
Boundary conditions are set only at left (inlet) and right (outlet)
patches. Inlet is supposed to be uniform */

u.n[left]  = dirichlet(USTREAM);
u.t[left]  = dirichlet(0);
u.n[right] = neumann(0);

p[left]    = neumann(0);
p[right]   = dirichlet(0);

void doublerainbow (double cmap[NCMAP][3]){
  for (int i = 0; i < NCMAP; i++){
    cmap[i][0]=sq(sin((double)i*M_PI/64.));
    cmap[i][1]=sq(sin(((double)i+20.)*M_PI/64.));
    cmap[i][2]=sq(sin(((double)i+40.)*M_PI/64.));           
  }//Black saturation
  cmap[NCMAP-1][0]=0;
  cmap[NCMAP-1][1]=0;
  cmap[NCMAP-1][2]=0;   
}

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
  origin(-3*DIAMETER,-L0/2, -L0/2);
  init_grid (pow(2.0,minlevel));



  // CFL number
  CFL = 0.5;

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
  
  TOLERANCE = 1e-6; 

 run();

  

 
}

/**
Initial condition is given by drop position only, there is no flow at the begining.
First step should give a velocity field close to potential flow solution around the drop,
plus small velocities inside the droplet.*/

event init (t = 0) {
  if (!restore (file = "restart")){
	refine (sq(x) + sq(y) + sq(z) - sq(0.6*DIAMETER) < 0 && sq(x) + sq(y) + sq(z) - sq(0.4*DIAMETER) > 0 && level < maxlevel);
        fraction (f, sq(0.5*DIAMETER) - (sq(x) + sq(y) + sq(z)));
	}
  }

/**
Only ten steps are performed, just to check mesh adaptation. */

event end (t=MAXTIME) {
  printf ("i = %d t = %g\n", i, t);
}


#if TREE
event adapt (i++) {
  double uemax = 1e-2;
  
    adapt_wavelet ({f,u}, (double[]){1e-3,uemax,uemax,uemax}, maxlevel);
}
#endif

/**
Every ten timesteps, we output the time, timestep, multigrid iterations, CPU time, speed and amount of cells. */

event diameter_data (i++) {
  scalar pos_x[], pos_y[];
  position (f, pos_y, {0,1}); //y-axis
  position (f, pos_x, {1,0}); // along stream
  FILE *data = fopen( "diameter" , "a");
  fprintf (data, "%g %g %g\n", t, 2.*statsf(pos_y).max,(statsf(pos_x).max - statsf(pos_x).min));
  fclose(data);
}


event logfile (i ++,first) {
   if (i == 0){
    fprintf (ferr, "i dt t mgp.i mgpf.i mgu.i perf.t perf.speed grid->tn\n");
  }
  fprintf (ferr, "%d %g %g %d %d %d %.2e %.2e %ld \n", i, dt, t, mgp.i, mgpf.i, mgu.i, perf.t, perf.speed, grid->tn);
  fflush (ferr);
}

event interface (i += 100; t <= MAXTIME)
{
  char name[80];
  sprintf(name, "interface-%d", i);
  FILE *file = fopen (name , "a");
  output_facets(f, file);
  fclose(file);
}



//event viewing (i += 250; t <= MAXTIME) {
event viewing (i ++) {

view (quat = {-0.207, -0.189, -0.010, 0.960},fov = 30, near = 0.01, far = 1000,
      tx = -0.290, ty = -0.030, tz = -1.483,width = 1548, height = 970,bg = {1,1,1});
draw_vof("f", color="sqrt(u.x[]*u.x[] + u.y[]*u.y[] + u.z[]*u.z[])", min=0, max=USTREAM, linear=true, map = doublerainbow);
//draw_vof("f");
//squares("sqrt(u.x[]*u.x[] + u.y[]*u.y[] + u.z[]*u.z[])*f[]", min=0, max=USTREAM, linear=true);
char name[80];
sprintf(name , "vof-%d.png", i);
save (name);
}


/*
event viewing (i += 500; t <= MAXTIME) {
view (quat = {-0.207, -0.189, -0.010, 0.960},
      fov = 30, near = 0.01, far = 1000,
      tx = -0.290, ty = -0.030, tz = -1.483,
      width = 1548, height = 970,bg = {1,1,1});
  box (notics = true);

  draw_vof ("f");
  
  cells();
  char name[80];
  sprintf(name , "grid-%d.png", i);
  save (name);
}


event snapshot (i += 1000; t <= MAXTIME) {  
  char name[80];
  sprintf (name, "dump-%d",  i);
  dump (file = name);
}
*/
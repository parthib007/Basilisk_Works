/**
# Impact of a viscoelastic drop on a solid

We solve the axisymmetric, incompressible, variable-density,
Navier--Stokes equations with two phases and use the log-conformation
method to include viscoelastic stresses. The curvature module is used
to compute interface properties. */

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "view.h"

// excluding the visoelastic properties

// #include "log-conform.h"
// #include "curvature.h"

#include "tension.h"

/**
The density and viscosity ratios are 1000. The Reynolds number based
on the droplet diameter, velocity and viscosity is 5 and the Froude
number is 2.26. */


#define RHO_r 0.001
#define MU_r 0.001
#define RE 5.
#define FR 1.63 //originally 2.26
#define LEVEL 7 // originally 7
#define DIAMETER 1.


/**Liquid properties
#define RHOL 785.
#define MUL 16.1e-3
#define SIGMA 48.3e-3

/**Air properties
#define RHOG 1.
#define MUG 1.81e-5

/**Drop diameter and stream velocity
#define DIAMETER 2.7e-3 
#define USTREAM 15.

#define MAXTIME 30e-3
#define MAXITER 10
*/

/**
The dimensionless viscoelastic properties used are the ratio of
the solvent to the total viscoelastic viscosity (polymeric
plus solvent) and the Weissenberg number. */

/**
#define BETA 0.1
#define WI 1.0

scalar lambdav[], mupv[];
*/

/**
The drop comes from the right. We allow the fluid to get through that
boundary. */

/**
u.n[right] = neumann(0);
p[right]   = dirichlet(0);
*/

/**
The wall is at the left side. We apply a no-slip boundary condition and a 
non-wetting condition for the VOF tracer. */

/**
u.t[left] = dirichlet(0);
// tau_qq[left] = dirichlet(0);
f[left]   = 0.;
*/

// writing my own boundary conditions where the droplet will fall from top except having a wall in the bottom
// here right -> inlet (top) , left -> outlet (bottom) 

#define USTREAM 11.
u.n[left] = dirichlet(USTREAM);
u.n[right] = neumann(0.);
p[left] = neumann(0.);
p[right] = dirichlet(0.);


int main() {

  /**
  The domain spans $[0:2.6]\times[0:2.6]$. */
  L0 = 10*DIAMETER;
  size (L0); // originally it was 2.6
 // origin(-L0/4,-L0/2);
  init_grid (1 << LEVEL); // 128 grids

  
  /**
  The densities and viscosities are defined by the parameters above. */

  rho1 = 1.;
  rho2 = RHO_r;
  mu1 = 0.1/RE;
  mu2 = MU_r/RE;

  /**
  The viscoelastic fields will be set below. */
  
/**
  mup = mupv;
  lambda = lambdav;
*/

  /**
  We set a maximum timestep. This is necessary for proper temporal
  resolution of the viscoelastic stresses. */
  
  DT = 2e-3;
  
  //DT = 1e-6;
  run();
}

event init (t = 0) {

  /**
  At a wall of normal $\mathbf{n}$ the component of the viscoelastic
  stress tensor $tau_p_{nn}$ is zero. Since the left boundary is a wall, we
  set $tau_p_{xx}$ equal to zero at that boundary. */

/**
  scalar s = tau_p.x.x;
  s[left] = dirichlet(0.);
*/

  /**
  The drop is centered on (2,0) and has a radius of 0.5. */

// for masking
//	mask (y > 1.5 ? top : none);
	
	if(!restore(file="restart")){
		mask (y>5*DIAMETER ? top : none);
		refine(sq(x-1.) + sq(y) - sq(0.6*DIAMETER) && sq(x-1.) + sq(y) - sq(0.4*DIAMETER) && level<LEVEL);
		fraction (f, - sq(x- 1.) - sq(y) + sq(DIAMETER*0.5));
	}

  /**
  The initial velocity of the droplet is -1. */

/**
  foreach()
    u.x[] =  f[];
*/
}

/**
We add the acceleration of gravity. */

/**
event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] -= 1./sq(FR);
	
}
*/

/**
We update the viscoelastic properties. Only the droplet is
viscoelastic. */

/**
event properties (i++) {
  foreach() {
    mupv[] = (1. - BETA)*clamp(f[],0,1)/RE;
    lambdav[] = WI*clamp(f[],0,1);
  }
}
*/

/**
We adapt the solution at every timestep based on the interface and
velocity errors. */

// #if TREE
event adapt (i++) {
  adapt_wavelet ({f, u.x, u.y}, (double[]){1e-2, 5e-3, 5e-3}, LEVEL, LEVEL - 2);

}
// #endif

// adding my own logfile for drop project

/**
We track the spreading diameter of the droplet. */

/*
event diameter_data (i += 20; t <= 4) {
  scalar pos_x[], pos_y[];
  position (f, pos_y, {0,1}); //y-axis
  position (f, pos_x, {1,0}); // along stream
  FILE *data = fopen( "diameter" , "a");
  fprintf (data, "%g %g %g %g %ld\n", t, dt, 2.*statsf(pos_y).max,(statsf(pos_x).max - statsf(pos_x).min), grid->tn);
  fclose(data);
}


event logfile (i += 20; t <= 4) {
   if (i == 0){
    fprintf (ferr, "t dt mgp.i mgpf.i mgu.i perf.t perf.speed grid->tn\n");
  }
  fprintf (ferr, "%g %g %d %d %d %.2e %.2e %ld \n", t, dt, mgp.i, mgpf.i, mgu.i, perf.t, perf.speed, grid->tn);
  fflush (ferr);
}
*/

/*
event interface (t += 1.)
{
  char name[80];
  sprintf(name, "interface-%d", (int)t);
  FILE *file = fopen (name , "a");
  output_facets(f, file);
  fclose(file);
}
*/

/**
We generate a movie of the interface shape. */

//trying to use output_ppm for movie

/**
event movies(i += 20; t <= 5)
{
	
//	output_ppm (f, file = "f-movie.mp4", box = {{-L0/4,-L0/2},{3*L0/4,L0/2}},linear = true);
	output_ppm (f, file = "f-movie.mp4", linear = true);

}
*/


/*
event end(t=end){
	char name[80];
	sprintf(name, "restart");
	dump(file=name);
}
*/

// graph plotting

/**
event graph(i++)
{
	char name[80];
	sprintf(name,"data");
	scalar pos[];
	position (f, pos, {0,1});
	stats s = statsf (pos);
	sprintf(name, "%g %g\n",t, 2.*s.max);
	//fclose(name);
}
*/

event viewing (i += 20; t <= 4) {
  //view (width = 400, height = 400, fov = 20, ty = -0.5,	theta =1.57079632679, phi =0, psi =0);
   view (tx = -0.5, ty = -0.5, theta = 0, fov = 30);
  clear();
  draw_vof ("f", lw = 2, linear = true, spread = 10);
  // squares ("f", spread = -1, linear = true, map = cool_warm);
  // squares ("u.x", linear = true);
  box (notics = true);
  cells();
//  mirror ({1,0}) {
//    draw_vof ("f", lw = 2);
//   squares ("omegaz", linear = true);
//    box (notics = true);
//  }
  save ("movie.mp4");
}


  
  /**
  We can optionally visualise the results while we run. */

#if 0
  static FILE * fp = popen ("bppm","w");
  save (fp = fp);
#endif


/**
## Results

![Animation of the interface shape. The color field on the right-hand-side 
(resp. l.h.s.) is the radial (resp. axial) velocity component.](fall/movie.mp4)

The time evolution of the maximum diameter is close to that reported
by [Figueiredo et al. (2016)](#figueiredo2016).

~~~gnuplot Time evolution of the maximum diameter
reset
set ylabel 'Maximum diameter'
set xlabel 't'
plot 'fall.figueiredo' lt 3 pt 5 t 'Figueiredo et al. (2016)', \
     'log' w l lw 2 t 'Basilisk'
~~~

## References

~~~bib
@article{figueiredo2016,
  title={A two-phase solver for complex fluids: Studies of 
  the {W}eissenberg effect},
  author={Figueiredo, RA and Oishi, CM and Afonso, AM and Tasso, 
  IVM and Cuminato, JA},
  journal={International Journal of Multiphase Flow},
  volume={84},
  pages={98--115},
  year={2016},
  publisher={Elsevier}
}
~~~
*/
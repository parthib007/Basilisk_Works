#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
//~ #include "reduced.h"
//#include "output_vtu_foreach.h" 
#include "navier-stokes/perfs.h"
//~ #include "maxruntime.h"

int MAXLEVEL = 11;
double  radius = 1.0;

//~ No slip, no penetration wall
u.n[left] = dirichlet(0.);
u.t[left] = dirichlet(0.);
p[left] = neumann(0.);

//~ No slip, no penetration wall
u.n[top] = dirichlet(0.);
u.t[top] = dirichlet(0.);
p[top] = neumann(0.);

//~ No slip, no penetration wall
u.n[right] = dirichlet(0.);
u.t[right] = dirichlet(0.);
p[right] = neumann(0.);

//~ Axis of symmetry
uf.n[bottom] = 0.;

int main(int argc, char* argv[]) {
  
	size (15.0);

  MAXLEVEL=atoi(argv[1]);
  N = 1 << MAXLEVEL;
  origin (0.0,0.0);

  rho1 = 1., mu1 = 0.001;
  rho2 = 0.001, mu2 = 0.00089;
  
	f.sigma = 72.0;
	
  TOLERANCE = 1e-4;
  run();
}

event init (t = 0) {
	//~ if (!restore(file = "restart")){
	//~ Initial shape of droplet
		fraction (f, sqrt(sq(x-7.5)+sq(y))<radius );
	//~ }
}

//~ Define gravity (All units are CGS)
event acceleration (i++) {
	face vector av = a;
	foreach_face(x)
	av.x[] = 981.0;
}

event logfile (i++) {
  if (i == 0)
   fprintf (ferr,"t dt \n");
   fprintf (ferr, "%g %g \n", t, dt);
}

//~ Output for post-processing
event snapshot (t = 0.; t += 0.001; t <= 5.0) {
	char name2[80];
	sprintf(name2, "bessel-%g.png",t*1000);
	output_ppm(f, file = name2, spread = -1, linear = true, map = cool_warm);
	
}
 
//~ #define IN_REGION (y < 5.0*radius)
//~ event adapt (i++) {
  //~ scalar f1[];
  //~ scalar region[];
  //~ foreach() 
  //~ f1[] = f[]; 	
  //~ foreach()
  //~ region[] = IN_REGION*noise();
  //~ adapt_wavelet ({region, f}, (double[]){0.01,1e-3}, minlevel = 6, maxlevel = MAXLEVEL);
//~ }

//~ Adaptivity
event adapt (i++) {
	double uemax = 0.01*normf(u.x).avg;
	adapt_wavelet ({f,u}, (double[]){1e-4,uemax,uemax}, minlevel = MAXLEVEL-2, maxlevel = MAXLEVEL);
}

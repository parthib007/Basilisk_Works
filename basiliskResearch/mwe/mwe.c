#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"

#define maxlevel 11
#define minlevel 6

#define RHOL 824.
#define MUL 2.17e-3
#define SIGMA 2.e-2

#define RHOG 1.1856
#define MUG 1.84e-5

#define USTREAM 9.9198

u.n[left] = dirichlet(USTREAM); 
u.t[left] = dirichlet(0);
p[left] = neumann(0);
p[right] = dirichlet(0);
u.n[right] = neumann(0);
u.t[right] = dirichlet(0);

int main()
{
	N = 1<<minlevel;
	L0 = 5;
	run();
	f.sigma = SIGMA;
	rho1 = RHOL;
	rho2 = RHOG;
	mu1 = MUL;
	mu2 = MUG;
	
}

event init(t=0)
{
	#if TREE
		refine (sq(x-1) + sq(y) - sq(0.6*1) < 0 && sq(x-1) + sq(y) - sq(0.4*1) > 0 && level < maxlevel);
	#endif
	fraction(f, sq(0.5) - sq(x-1) -sq(y));
}

event end(i=10)
{
	printf("%d %g", i ,t);
}

#if TREE
event adapt (i++){
    adapt_wavelet((scalar *){f,u}, (double []){1e-3,1e-2,1e-2}, maxlevel);
}
#endif

event logdata(i++)
{
	if(i==0)
		fprintf(ferr, "i dt t");
	fprintf(ferr, "%d %g %g", i, dt, t);
	fflush(ferr);
}

event snapshot(i++)
{
	char name[80];
	view( tx = -0.5, ty =-0.5, fov = 20);
	draw_vof("f");
	squares("f" , spread = -1, linear = true, map = randomap);
	box(notics = true);
	sprintf(name, "vof-%d.png",i);
	save(name);
}
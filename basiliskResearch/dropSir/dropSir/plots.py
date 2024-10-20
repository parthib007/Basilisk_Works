# -*- coding: utf-8 -*-
# /**
# # Droplet blow dy stream, bag mode secondary atomisation
# 
# We wish to study the behaviour of a single drop in free fall
# affected dy a uniform flow of air, as in the [Opfer 2014](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.89.013023).
# 
# We use the centered Navier--Stokes solver. There are two phases, air and the liquid. The grid adaptation has been
# modyfied to support region-depending maximum refinement (MAXLEVEL).*/
# 
# //#include "axi.h"
# #include "grid/octree.h"
# #include "navier-stokes/centered.h"
# #include "two-phase.h"
# #include "navier-stokes/conserving.h"
# #include "tension.h"
# #include "reduced.h"
# #include "view.h"
# 
# 
# /**Liquid properties*/
# #define RHOL 824.
# #define MUL 2.17e-3
# #define SIGMA 2.e-2
# 
# /**Air properties*/
# #define RHOG 1.196
# #define MUG 1.822e-5
# 
# /**Drop diameter and stream velocity*/
# #define DIAMETER 2.4e-3 
# #define USTREAM 9.9198
# 
# #define MAXTIME 15.e-3
# //#define MAXITER 5000
# 
# /**
# Minimum refinement level will be used directly, MAXLEVEL will only apply
# on a certain region of the domain. */
# 
# int minlevel = 6;   // init grid of 32 points
# int MAXLEVEL = 11;  //11
# 
# double maxruntime = HUGE;
# 
# 
# /**
# Boundary conditions are set only at left (inlet) and right (outlet)
# patches. Inlet is supposed to be uniform */
# 
# u.n[left]  = dirichlet(USTREAM);
# u.t[left]  = dirichlet(0);
# u.n[right] = neumann(0);
# 
# p[left]    = neumann(0);
# p[right]   = dirichlet(0);
# /**
# The main function can take the maximum refinement level as argument.
# The main program will run ten steps with both versions of AMR.*/
# 
# int main (int argc, char * argv[]) {
#   if (argc > 1)
#     MAXLEVEL = atoi (argv[1]);
#   
#   /**
#   We set the domain geometry and initial refinement.
#   The drop will be centered at the origin, slightly
#   moved to the left side*/
#   L0 = 5* DIAMETER;
#   size (L0);
#   origin(-L0/4,-L0/2, -L0/2);
#   init_grid (pow(2.0,minlevel));
#   
#   G.y = -9.8;
# 
# 
#   // CFL number
#   CFL = 0.5;
# 
#   /**
#   Physical parameters are set by the macros at the beggining.*/
#   
#   f.sigma = SIGMA;
#   rho1 = RHOL;
#   rho2 = RHOG;
# 
#   mu1 = MUL;
#   mu2 = MUG;
# 
#   /**
#   Tolerance of the Poisson problem should never be reduced,
#   this could lead to non divergence free velocity fields
#   and mass conservation errors. */
#   
#   TOLERANCE = 1e-4; 
# 
#  run();
# 
#  //limitedAdaptation = 1;
#   
#  //run();
#   
# 
#  
# }
# 
# /**
# Initial condition is given by drop position only, there is no flow at the begining.
# First step should give a velocity field close to potential flow solution around the drop,
# plus small velocities inside the droplet.*/
# 
# event init (t = 0) {
#   if (!restore (file = "restart")){
#     // mask(y > 5*DIAMETER? top:none);
#     refine (sq(x-DIAMETER) + sq(y) + sq(z) - sq(0.6*DIAMETER) < 0 && sq(x-DIAMETER) + sq(y) + sq(z) - sq(0.4*DIAMETER) > 0 && level < MAXLEVEL);
#     fraction (f, sq(0.5*DIAMETER) - (sq(x-DIAMETER) + sq(y) + sq(z)));
#     output_ppm(f,linear = true, n=512, box = {{0,0},{L0,L0/2}}, file = "vofinit.png");
# 	}
#   }
# 
# /**
# Only ten steps are performed, just to check mesh adaptation. */
# 
# event end (i=3) {
#   printf ("i = %d t = %g\n", i, t);
# }
# 
# /**
# ## Region limited mesh adaptation
# 
# Any function with input(x,y,z) returning an integer can be used to
# define MAXLEVEL locally for certain regions. Next, a simple example
# using a circle (only for 2D case) slightly bigger than the drop is used.
# Inside that region, AMR can refine until MAXLEVEL, outside
# only MAXLEVEL-2 cells are allowed. Then, the alternative function
# adapt_wavelet_limited is employed. */
# 
# event adapt (i++) {
#   double uemax = 5e-2;
#   
#   /**
#   Choose between region limited adaptation or standard adaptation: */
#     adapt_wavelet ({f,u}, (double[]){1e-3,uemax,uemax,uemax}, MAXLEVEL);
# }
# 
# /**
# Every ten timesteps, we output the time, timestep, multigrid iterations, CPU time, speed and amount of cells. */
# 
# event diameter_data (i++) {
#   scalar pos_x[], pos_y[];
#   position (f, pos_y, {0,1}); //y-axis
#   position (f, pos_x, {1,0}); // along stream
#   FILE *data = fopen( "diameter" , "a");
#   fprintf (data, "%g %g %g\n", t, 2.*statsf(pos_y).max,(statsf(pos_x).max - statsf(pos_x).min));
#   fclose(data);
# }
# 
# 
# event logfile (i ++,first) {
#    if (i == 0){
#     fprintf (ferr, "i t dt mgp.i mgpf.i mgu.i perf.t perf.speed grid->tn\n");
#   }
#   fprintf (ferr, "%d %g %g %d %d %d %.2e %.2e %ld \n", i, t, dt, mgp.i, mgpf.i, mgu.i, perf.t, perf.speed, grid->tn);
#   fflush (ferr);
# }
# 
# event interface (i ++)
# {
#   char name[80];
#   sprintf(name, "interface-%d", i);
#   FILE *file = fopen (name , "a");
#   output_facets(f, file);
#   fclose(file);
# }
# 
# /**
# First ten steps are saved in gfsview to compare initial flow field and mesh. */
# 
# /*event snapshot (i=0;i++; i <= 10)
# {
#    
#   scalar omegaz[];
#   foreach()
#     omegaz[] = (u.y[1] - u.y[-1] - u.x[0,1] + u.x[0,-1])/(2.*Delta);
#   boundary ({omegaz});
#   
#   char name[80];
#   if(limitedAdaptation) 
#     sprintf (name, "snapshot-limited-%i_ts.ppm", i);
#   else
#     sprintf (name, "snapshot-standard-%i_ts.ppm", i);
#   //output_gfs (file = name, t = t, list = {f,u,p, omegaz});
#   output_ppm ( omegaz , file = name );
# } */
# 
# 
# event snapshot(i ++)
# {
#   
#   char name[80];
#   sprintf (name, "snapshot-f-%i_ts.ppm", i);
#   output_ppm ( f ,linear = true,n = 800, box = {{0,0},{L0,L0/2}}, file = name );
# 
# }
# 
# /*
# event restarting(i += 1000)
# {
# 	char name[80];
# 	sprintf(name, "restart");
# 	dump(file = name);
# }
# */
# 
# /*
# //event images(i += 5; t <= 15.)
# event images(i += 5; i <= MAXITER)
# {
# 	scalar l[];
# 	foreach()
# 	 	l[] = level;
# //		char name[80];
# //		sprintf (name, "snapshot-l-%i_ts.ppm", i);
# //		output_ppm (l,linear = true,min = minlevel, max = MAXLEVEL,file = name);
# 		static FILE * fp = fopen ("grid.ppm", "w");
#                output_ppm (l, fp, min = 0, max = 8); // min max color levels
#                char name[80];
# 
# } 
# */
# 
# /*
# event movies(i += 250; t <= MAXTIME)
# {
# 
# //output_ppm (f, file = "f-movie.mp4", box = {{-L0/4,-L0/2},{3*L0/4,L0/2}},min =0, max = 10, linear = true);
#  output_ppm (f, file = "f-movie.mp4",box = {{0,0},{L0,L0/2}},linear = true;
# }
# */
# 
# /*
# event movies(i += 250; t <= MAXTIME)
# {
#   view (width = 1024, height = 768);
#   draw_vof("f", color = "Z", linear = true);
#   save ("movie.mp4");
# }
# */
# 
# 
# 
# event movie (i ++) {
#   output_ppm (f, linear = true, n = 800, box = {{0,0},{L0,L0/2}}, file = "vof.mp4", spread = -1, map = cool_warm);
#   output_ppm ( u.x, box = {{0,0},{L0,L0/2}}, n = 800, file = "speed.mp4" , min = 0, max = USTREAM, map = cool_warm);
#   //scalar l[];
#   //foreach()
#     //l[] = level;
# //  output_ppm (l, min = minlevel, max = MAXLEVEL,box = {{0,0},{L0,L0/2}}, file = "level.mp4");
#   //output_ppm (l, n = 800, box = {{0,0},{L0,L0/2}}, file = "level.mp4", min = 0, max = 1, linear = false);
# }
# 
# 
# event viewing (i ++) {
#   view (camera = "iso", fov = 14.5, tx = -0.418, ty = 0.288, width = 1600, height = 1200);
#   
#   draw_vof ("f");
#   char name[80];
#   sprintf(name , "grid-%d.png", i);
#   save (name);
#   
#   box();
#   squares ("f", linear = true);
#   cells();
#   sprintf (name, "cross-%d.png", MAXLEVEL);
#   save (name);  
# }
# 
# 
# /*
# event pictures (t = end)
# {
#   char name[80];
# #if 0  
#   sprintf (name, "dump-test");
#   dump (name);
# #endif
# 	scalar l[];
# 	foreach()
# 	 	l[] = level;  
# //  view (width = 600, height = 600);  
#   squares ("l", linear = true, map = cool_warm);
#   draw_vof ("f");
#   
#   sprintf (name,"final.png");
#   save (name);
# }
# 
# */
# /**event movie (t+=DT/100)
# {
# 	scalar omegaz[];
# 	vorticity(u,omegaz);
# 	foreach()
# 	omegaz[] = (u.y[1] - u.y[-1] - u.x[0,1] + u.x[0,-1])/(2.*Delta);
#   boundary ({omegaz}); 
# 	view(tx=-3*L0/4);
# 	clear();
# 	draw_vof("f");
# 	squares("omegaz", linear = true, spread = 10);
# 	box();
# 	save("movie.mp4");
# } */
# /**
# ##Mesh evolution
# 
# From the log of this example, a graph as the following can be obtained:
# 
# [Mesh size during runtime](/bagMode_meshInTime.pdf)
# 
# */

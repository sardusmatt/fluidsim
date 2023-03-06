**2D/3D Navier-Stokes Fluid Simulator**

In computer graphics, gaseous phenomena are often modelled as fluid flows, whose animation requires to numerically solve the Navier-Stokes equations that describe the motion
of such fluids.

The simulator implements an unconditionally stable physical model first presented by Stam, which can produce fluid-like (as in, visually plausible)
complex flows by solving the Navier-Stokes equations in (semi) real-time.This approach applies both to two and three dimensions, and is particularly suitable for smoke simulations. 

To produce even more realistic motion of smoke-like fluids, the simulator also implements an additional technique called _vorticity confinement_, which counters the numerical dissipation (i.e., visual dampening of the flow) introduced by the original solving algorithm.

Last but not least, the simulator makes use of several different shading techniques, including volumetric shading, to produce better renderings of the simulated fluids.
It also support offline rendering, to be able to simulate and render larger 3D fluid grids.
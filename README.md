# Fluid Simulator

Fluid simulator based on the 2D incompressible Navier-Stokes equations. You can find the theoritical framework for the simulation in [1].

## What is implemented

* Complete simulation of the 2D imcompressible Navier-Stokes equations, including obstacles.
* Add domain obstacles in the domain from PNG images.
* Parallel solver (SOR) for the Poisson pressure equation using OpenMP, which is constitutes the major computation time of the simulation.

## Compilation

You'll need the have a C/C++ compiler with OpenMP enabled and GNU make installed. There are two modes of compilation. One disables the debug messages, asserts and enables compiler optimizations (release mode) and the other enables them (debug mode).

To compile in release mode, type in a terminal:
> make

To compile in debug mode, type in a terminal:
> make TAG=GCCdebug

In this case, the GCC compiler toolchain is used to compile the project. However, one can also use ICC or LLVM by changing the 

## Simulations

The following (classical) simulation scenarios where tested with different Reynolds numbers and are present in the _simulations_ directory:
* Lid-Driven Cavity
* Flow through a channel
* Flow through a channel with a backstep

## References
[1] Griebel, Michael, Thomas Dornseifer, and Tilman Neunhoeffer. Numerical simulation in fluid dynamics: a practical introduction. Vol. 3. Siam, 1997.

# SPH Code - Project Overview

## Purpose
This is a Smoothed Particle Hydrodynamics (SPH) simulation code for compressible fluid dynamics. Written in C++ with OpenMP parallelization, it implements various SPH methods for astrophysical and computational fluid dynamics simulations.

## Key Features
- **Multiple SPH formulations**: Standard SPH (SSPH), Density Independent SPH (DISPH), Godunov SPH (GSPH)
- **Kernel functions**: Cubic spline, Wendland C4
- **Artificial viscosity**: Signal velocity formulation, Balsara switch, time-dependent AV
- **Self-gravity**: Tree-based gravity solver (Barnes-Hut)
- **Lane-Emden relaxation**: Polytropic sphere relaxation with checkpoint/resume capability
- **Multiple test cases**: Shock tube, KHI, Evrard collapse, Gresho-Chan vortex, etc.

## Tech Stack
- **Language**: C++14
- **Build system**: CMake (primary), legacy Makefile (deprecated)
- **Parallelization**: OpenMP
- **Dependencies**: 
  - Boost (for property tree JSON parsing)
  - OpenMP (for parallel execution)
- **Compiler**: GCC 7.4.0+ or compatible C++14 compiler

## Project Type
Scientific simulation code for astrophysics/computational fluid dynamics research.

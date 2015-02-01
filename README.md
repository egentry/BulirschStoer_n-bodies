# BulirschStoer_n-bodies
Bulirsch-Stoer integrator, with sample applications, as used in UCSC Astro 212
-------

Author: Eric Gentry   (gentry.e@gmail.com; egentry@ucsc.edu)   

Copyright (c) 2015
Licensed under the MIT License

-------

##Main Objectives
  - Given N-bodies, initial conditions, and a force term, numerically integrate
    - integration should accurately resolve numerically-difficult close interactions (e.g. within the Pythagorean 3-4-5 problem)
  - Generalize easily to different potentials and number of objects
  - Include no routines from Numerical Recipes

 

##Requires
  - Modern Fortran compiler
    - Must be at least approximately Fortran 2008 compliant
    - `gfortran` required for  `makefile`

##Getting Started
 - Choose an example
 - Switch to that directory
 - Run: `ipython notebook visualize_integration.ipynb`
 - s
  
##Optional Packages
  - `python` (2 or 3)
    - `ipython` (version > 2.0)
    - `numpy`
    - `matplotlib`
  - LaTeX (for use with matplotlib; not essential)

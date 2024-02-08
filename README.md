# fmm3dbie software demo

## Software demo for Workshop on Computational Tools for PDEs with Complicated Geometries and Interfaces

External dependencies
-----------------------
- [FMM3D](https://fmm3d.readthedocs.io/en/latest): Fast multipole libraries in
  three dimensions
- [fmm3dbie](https://fmm3dbie.readthedocs.io/en/latest): Fast multipole
  accelerated iterative solvers for integral equations
- [chebfun](https://www.chebfun.org): Numerical computing with functions
- [Surface smoother](https://github.com/fastalgorithms/smooth-surface): Converts
  flat triangulated meshes into high order meshes

Demos
---------------
- Demo 1: Illustrate loading of meshes, verify green's identity at an exterior
  point using smooth quadrature

- Demo 2: Compare application of layer potential on a sphere with its
  analytically known solution on and off surface

- Demo 3: Illustrate the use of the Surface Smoother GUI to convert a 
  flat triangulated boxed torus to a high-order smoothed boxed torus

- Demo 4: Illustrate Neumann solver on boxed torus geometry

Problem sets
----------------------

- Problem set 1: Illustrate order of convergence for a stellarator like
  geometry

- Problem set 2: Compute bistatic radar cross-section of a smoothed boxed torus


Note that solutions to the problem sets are included in the `solutions/` folder

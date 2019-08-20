# MPI_Spectral_Navier_Stokes
Parallel incompressible Navier-Stokes solver in a torus (periodic box), using Fourier pseudospectral method.


### Requirements: 1) OpenMPI, 2) GNU Fortran compiler, 3) A legacy .vtk file reader (e.g., Paraview of VisIt).

### How to run the code: in /BIN, type "make". Modify the input file as desired. Results will be saved in the /RESULTS/FLUID_FIELD/ directory, for separated processors. You can unify the outputs later using filters in Paraview or VisIt (e.g., "group datasets" in Paraview).

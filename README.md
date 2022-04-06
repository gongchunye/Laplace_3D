# Laplace_3D
Uses MPI parallel programming to code laplace equation in 3D


Initial Conditions:

100x100x100 Cube
Faces of the cube: top, bottom, right, left, front, back

Bottom, Right, Front face: T=100K

Everywhere else: T = 0K

(No gradient edges)


Instructions to run:

mpicc laplace_MPI_3D.c

mpirun -n 4 a.out

Maximum iterations [100-4000]?

>>4000 (converges at 3424)


The code will output files depending on PE number, each of which can be run through VisIt
producing sections of the cube.


1)-------------------------------------------------------------
old readme
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

2)-------------------------------------------------------------

how to change max_iterations?
max_iterations= 1000; 

3)-------------------------------------------------------------
how to change MPI process number:
#define NPES           4       // number of processors

4)-------------------------------------------------------------
how to change problem size:
#define ROWS_GLOBAL  100         // this is a "global" row count


5)-------------------------------------------------------------
may be trouble for big simulation with double Temperature[ROWS+2][COLUMNS+2][LAYERS+2];
static allocate memory

6) not test for too large problem size

7) results see: https://www.ydma.com/article-22532-1.html



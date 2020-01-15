Schrodinger Eigenvalue Solver:

Author: Sam Altier

Date: 12/13/2019
******************************************************************************************
This code set is designed to use the Jacobi rotation algorithm to numerically evaluate energy eigenvalues of the Schrodinger equation.
It uses the Jacobi algorithm for determination of eigenvalues for energy states.

******************************************************************************************
CODE:

Compiling:

-> ..._schro_solve.cpp, ..._harmonic_solve.cpp, fit_poly.cpp, jacobi.cpp with dependencies:

	g++ [file].cpp matrixmath.cpp

-> All others can be compiled as is.


Running:




******************************************************************************************
----------------------
-> ..._schro_solve.cpp, ..._harmonic_solve.cpp : performs Jacobi rotations to determine eigenvalues, outputs file with eigenvalues, file with all eigenvalues (in case of non-harmonic files)

./solve.exe filename M tol omegafactor

filename: potential file to be solved; harmonic_solve.cpp should only be provided quadratic potential files; schro_solve.cpp can take any potential
   Format:
   numpts
   x1	V1
   x2	V2
   ...

M: maximum number of iterations to perform

tol: tolerance to which to perform rotations

omegafactor: multiple of the space discretization to use as omega0

[all]: the harmonics have an optional parameter to include all of the eigenvalues in their output file 

----------------------
-> fit_poly.cpp : fits polynomial of order M to a dataset; used for analysis of generated eigenvalues

./fp.exe M filename

M: number of fit parameters

filename: file to be fit
   Format:
   numpts
   x1	V1
   x2	V2
   ...
   
----------------------
-> generate_potentials.cpp

./gp.exe N [dx]

N: number of points to generate, is forced to be odd so potential is symmetrical

[dx]: optional command to change discretization, default is 0.1 because is most effective for the method

----------------------
-> generatePureQuartic.cpp : generates the "analytic" solution to the pure quartic oscillator

./gpq.exe

No CMA

----------------------
-> compare_results.cpp : takes difference between two data sets and plots it

./cr.exe file1 file2

file1, file2: two files with the same number of points, takes the difference of file1 and file2 (no absolute value)



	

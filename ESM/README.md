# Readme file for the homogenization of the plate-array crystal

%%%%%%%%%%%%%
The Mathematica script "Maxwell_Garnett_equations" computes the determinant of the 8 by 8 dipole matrix and its expansion. 
alpha_0 is computed in the case of nr=1 and nr!=1. It also produces the surface and contour plots found in the paper. 
Note that this script takes several hours to run due to the expansion of the matrix determinant and the simplification of each term in the expansion.
%%%%%%%%%%%%%
The Matlab script "Dipole_Comparison" computes the numerical multipole results and compares these to the quasistatic approximations. 
It calls on the functions "LatticeSumMat_faster" and "getlatticesums_faster" to generate the lattice sums and "TwoDMetaMatS_nr" to generate 
the scattering matrix for the crystal.
%%%%%%%%%%%%%
The FEM data can be found in the "FEM_data" folder, and has to be moved to the main Effective Medium file in order for the "Dipole_comparison" 
script to run without errror. 
%%%%%%%%%%%%%
The colormaps can be found at https://www.fabiocrameri.ch/colourmaps/.
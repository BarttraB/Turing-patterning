# Turing-patterning

Supplementary Materials for analysis and simulation of model Eqs. (1-4) in Borek,et al (2016) "Turing patterning using gene circuits with gas-induced degradation of quorum sensing molecules"

The movies show 2D numerical simulation of the equations (using the nominal parameter set): 

1. Turing2Dmovie_periodicBC.mov 	has periodic boundary conditions.

2. Turing2Dmovie_nofluxBC.mov 		has no-flux boundary conditions.


The simulation scripts are:

1. modelE4_9_261.f		numerically solves the equations in 1D

2. modelE4_2D_noflux.f   	numerically solves the equations in 2D with no-flux boundary conditions

These scripts can be compiled using gfortran (ex.: in terminal type " gfortran  modelE4_9_261.f " and press enter, then type " a.out " and press enter to execute the program and generate the output .txt files)


The analysis scripts are:

1. evalWavenum_modE4_a3_a1is4p5.m	numerically solves for the eigenvalues for different wavenumbers. 

2. codim2_evalFP_modE4_a1_Pmax_iter4.m 	numerically solves for the parameter pairs around the zero discriminant from Eq. 5

These scripts are executed from inside Matlab by clicking on them or typing the script name in the command window. 


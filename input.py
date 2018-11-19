#!/usr/bin/env python

# input file for the stator geometry calculation.
# T4, p4, M4, rho4, 



# Define the trailing edge thickness 
D.ttmin = 1.0e-4  # m 
D.ttmax = 4.0e-4  # m
D.ttstep = 0.1e-4 # m

D.r3devr4min = 1.01
D.r3devr4max = 1.05


# stator properties
D.Z3min = 20
D.Z3max = 40
#D.b3    = 0.004251


D.fluidType = 'CO2'




#

D.r3       = 0.0485928   # m  0.0425187[1.05]   0.0485928[1.2]
D.delta    = 5.
# rotor properties/target properties
D.r4       = 0.040494   # m
D.T4       = 787.290547 # K
D.p4       = 14365321.162128 # pa
D.gamma    = 1.239186  
D.C4       = 322.402051   # m s-1
D.A4       = 0.001026   # m2
D.b4       = 0.004251   # m
D.alpha4   = 81.528855  # degree


D.b3       =  D.b4

#

#D.r3       = 0.020246   # m
#D.delta    = 1.
# rotor properties/target properties
#D.r4       = 0.020146   # m
#D.T4       = 795.485899 # K
#D.p4       = 15258730.196 # pa
#D.gamma    = 1.239186  
#D.C4       = 292.48     # m s-1
#D.A4       = 0.000114   # m2
#D.b4       = 0.000969   # m
#D.alpha4   = 71.146841  # degree


M.maxiter  = 100
M.gasModel = 'Ideal'#'Real'# 'Ideal'  # Ideal or Real


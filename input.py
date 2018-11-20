#!/usr/bin/env python

# input file for the stator geometry calculation.
# T4, p4, M4, rho4, 



# Define the trailing edge thickness 
D.ttmin = 1.0e-3  # m 
D.ttmax = 3.0e-3  # m
D.ttstep = 0.1e-3 # m


# Define the radius ratio range
D.r3devr4min  = 1.05
D.r3devr4max  = 1.35
D.r3devr4step = 0.02


# stator properties
D.Z3min = 14
D.Z3max = 30
D.Z3step = 1
#D.b3    = 0.004251


D.fluidType = 'CO2'




#




D.delta    = 2.
# rotor properties/target properties
D.r4       = 0.042016   # m
D.T4       = 795.091676 # K
D.p4       = 15214696.02081 # pa
D.gamma    = 1.239186  
D.C4       = 293.992329   # m s-1
D.A4       = 0.000602   # m2
D.b4       = 0.002398   # m
D.alpha4   = 73.300756  # degree


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

M.optimiser = 'Root' # or 'Nelder-Mead'
M.maxiter  = 500
M.gasModel = 'Real'#'Ideal'#'Real'# 'Ideal'  # Ideal or Real


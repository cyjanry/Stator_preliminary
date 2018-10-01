#!/usr/bin/env python

# input file for the stator geometry calculation.
# T4, p4, M4, rho4, 
input = []



V.tt = 0.0001 # m 



D.fluidType = 'CO2'
# stator properties
D.Z3 =    10
D.r3 =    0.020772   # m
D.delta = 1.
# rotor properties/target properties
D.r4 =    0.020146   # m
D.T4 =    795.485899 # K
D.p4 =    15258730.196 # pa
D.gamma = 1.239186  
D.C4 =    292.48     # m s-1
D.A4 =    0.000114   # m2
D.b4 =    0.000969   # m
D.alpha4= 71.146841  # degree


M.maxiter = 100
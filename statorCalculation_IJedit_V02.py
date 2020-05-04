#!/usr/bin/env python3
# This Python file uses the following encoding: utf-8



import os         as         os
import numpy      as         np
import shutil     as         sh
from   getopt            import getopt
import sys        as         sys 
import scipy      as sci
from   scipy.optimize    import minimize
import matplotlib.pyplot as plt 
import CoolProp.CoolProp as CP
import datetime



def pFunc(p3,M3,M4,gamma):
    gammaMinusOne = gamma - 1.
    p4 = p3 * pow((1. + 0.5 * gammaMinusOne* M4**2 )/( 1. + 0.5 * gammaMinusOne * M3**2) , -( gamma/gammaMinusOne ) )
    return p4 

def tFunc(T3,M3,M4,gamma):
    gammaMinusOne = gamma - 1.    
    T4 = T3 * (1. + 0.5 * gammaMinusOne * M3**2) / (1. + 0.5 * gammaMinusOne * M4**2)
    return T4 

def rhoFunc(rho3, M3, M4,gamma):
    gammaMinusOne = gamma -1.
    rho4 = rho3 * pow( (1. + 0.5 * gammaMinusOne * M4**2)/ (1. +0.5 * gammaMinusOne * M3**2)  , (-1./ gammaMinusOne) )
    return rho4




def CalInitial(M,V,D):
    # With given rotor initial conditions, calculate the rest rotor properties.

    D.gammaMinusOne = D.gamma - 1. 

    # Gas constant
    Ru =  CP.PropsSI(D.fluidType, "gas_constant")
    W  =  CP.PropsSI(D.fluidType, "molar_mass")
    D.R  = Ru/W
    D.Cp = D.R * D.gamma / (D.gamma-1)
    D.Cv = D.Cp - D.R
    # Calculate corresponding rotor fluid properties.

    # In the following section, we apply the real gas model
    #a4  = CP.PropsSI('A', 'T', D.T4, 'P', D.p4 , D.fluidType) # real gas model

    if M.gasModel == 'Ideal':

        D.a4    = np.sqrt(D.gamma * D.R * D.T4) # real gas model, acoustic speed.
        D.M4    = D.C4/D.a4
        D.T04   = D.T4 * (1. + 0.5* D.gammaMinusOne * D.M4**2 )
        D.p04   = D.p4 * pow(  (1. + 0.5* D.gammaMinusOne * D.M4**2), (D.gamma / D.gammaMinusOne)  )
        D.rho4  = D.p4/ (D.R*D.T4)
#        D.rho04 = D.rho4 * pow(  (1. + 0.5* D.gammaMinusOne * D.M4**2), (1 / D.gammaMinusOne)  )
        D.h4    = D.Cp * D.T4 
        D.h04   = D.Cp * D.T04
        # Reference conditions for entropy
        D.Tref = 298.15  #K 
        D.Pref = 1.0135e5 #Pa
        D.S0    = D.Cp * np.log(D.T4/D.Tref) - D.R * np.log(D.p4/D.Pref)
    else:
        D.a4    = CP.PropsSI('A','T',D.T4,'P',D.p4, D.fluidType)
        D.M4    = D.C4/D.a4

        # the following property is used in isentropic expansion
        D.S0    = CP.PropsSI('S', 'T',D.T4, 'P',D.p4, D.fluidType)
        D.h4    = CP.PropsSI('H', 'T',D.T4, 'P',D.p4, D.fluidType)
        D.h04   = D.h4 + D.C4**2/2
        D.T04   = CP.PropsSI('T', 'H', D.h04, 'S', D.S0, D.fluidType)
        D.p04   = CP.PropsSI('P',  'H', D.h04, 'S', D.S0, D.fluidType)

        D.rho4  = CP.PropsSI('D', 'T',D.T4, 'P',D.p4, D.fluidType) 
#        D.rho04 = CP.PropsSI('D',  'H', D.h04, 'S', D.S0, D.fluidType)

    #D.A4devb4 = D.A4/D.b4
    #A4cal = 2* np.pi*D.r4 - D.Z4*0.001 
    D.C4T   = D.C4 * np.sin( np.radians( D.alpha4) )
    D.C4R   = D.C4 * np.cos( np.radians( D.alpha4) )


    if M.verbosity >=1:
        print("The following are the rotor properties:")
        print("T04 is:",D.T04)
        print("p04 is:",D.p04)
        print("a4 is:", D.a4)
        print("M4 is:", D.M4)
        print("rho4 is:", D.rho4)
        print("A4devb4 is:",D.A4devb4)
        print("C4R is:",D.C4R)
        #print("A4cal", A4cal)
    
    # Assuming there is no loss

    D.T03 = D.T04
    D.p03 = D.p04
#    D.rho03 = D.rho04

    # Assuming there is no loss in momentum:
    D.C3T = D.C4T * ( 1/D.r3devr4min) # use the minimum r3/r4 ratio to initilize the calculation

    V.M3 = D.M4
    V.p3 = D.p4
    V.T3 = D.T4
    V.rho3 = D.rho4
    V.C3   = D.C4
    V.a3   = D.a4
    V.alpha3b = D.alpha4 # use alpha4 value as the intial value of stator blade

    return 0

def initSimplex(x):

    simplex = []
    simplex.append(x)
    for i in range(len(x)):
        temp = x[:]
        temp[i] = temp[i] + 0.1*temp[i]
        simplex.append(temp)
    print(simplex)
    print(np.array(simplex))

    return np.array(simplex)



def M4Func(M4temp,M,V,D,C4T,T04,p04):
    #
    if M.gasModel == 'Ideal':

        left  = ( V.C3R * V.rho3 *V.A3devb3 / (D.A4devb4 * rhoFunc(V.rho3, V.M3, M4temp,D.gamma) ) ) **2 +  C4T**2
        right =  M4temp**2 * D.gamma * D.R * tFunc(V.T3,V.M3,M4temp,D.gamma)
        cost = abs(left - right)
    else:

        H04 = CP.PropsSI('H', 'T', T04, 'P', p04, D.fluidType)
        S04 = CP.PropsSI('S', 'T', T04, 'P', p04, D.fluidType)
        A04 = CP.PropsSI('A', 'T', T04, 'P', p04, D.fluidType)

        H4temp = []
        H4temp.append( H04 - (M4temp[0] * A04)**2/2.)
        H4_simp = initSimplex(H4temp)
        result = minimize(HFuncReal, H4temp, args=(M,V,D,H04,S04,M4temp), method='Nelder-Mead',options={'initial_simplex': H4_simp, 'maxiter':M.maxiter})
        H4 = result.x[0]
        a4   = CP.PropsSI('A', 'H', H4, 'S', S04, D.fluidType)
        rho4 = CP.PropsSI('D', 'H', H4, 'S', S04, D.fluidType)
        #T4   = CP.PropsSI('T', 'H', H4, 'S', S04, D.fluidType)
        #p4   = CP.PropsSI('P', 'H', H4, 'S', S04, D.fluidType)

        left = ( V.C3R * V.rho3 *V.A3devb3 / (D.A4devb4 * rho4 ) ) **2 +  C4T**2
        right = M4temp**2 * a4**2
        cost = abs(left - right)

    return cost






def HFuncReal(Htemp, M,V,D, H0, S0, Mach):
    # function to calculate static enthalpy based on the given Mach number, stagnation enthalpy and stagnation entropy
    left = Htemp[0]
    acoustic   = CP.PropsSI('A','H', Htemp[0], 'S', S0, D.fluidType)
    right = H0 - ( Mach * acoustic)**2/2.
    cost = abs(left - right)
    return cost



    
def Calculation2(x,M,V,D):

    # New function called by the root optimiser, use constant S0 

    # unpack the state vector
    #x = [V.alpha3b, V.M3, V.p3]
    V.alpha3b = x[0]
    V.M3 = x[1]
    V.p3 = x[2]*1.e6  # here change mega pascal to pascal
    
    #print("COME here and alpha3b={0}, M3 = {1}, p3 ={2}".format(x[0],x[1],x[2]))
        
    # get geometry constraints, Z and tt
    tt = V.tt
    Z3 = V.Z3
    r3 = V.r3


    if M.verbosity > 0:
        print("\n")
        print("---------------------------")
        print("New Iteration Loop")
        print("X",x)
        print("alpha3b = {0:.3f} [deg]; M3 = {1:.3f} [-]; p3 = {2:.1f} [kPa]; T3 = {3:.4} [K]".format(x[0],x[1],x[2]/1.e3,V.T3))
        print("tt = {0:.5f} [m]; Z3 = {1} [-]".format(tt,Z3))


    #  calculate other state properties at point 3
    if M.gasModel == 'Ideal':
        V.T3 = D.Tref * np.exp((D.S0 + D.R * np.log(V.p3/D.Pref)) / D.Cp) 
        V.rho3 = V.p3/ (V.T3 * D.R)        
        V.h3  = D.Cp * V.T3
        V.a3  = np.sqrt(D.gamma * D.R * V.T3)
        V.p03 = V.p3 *   pow(  ( 1. + D.gammaMinusOne/2. * V.M3**2) , (D.gamma/ D.gammaMinusOne) )
        print("V.p03=",V.p03,"D.p04=",D.p04)

    else:
        V.T3 = CP.PropsSI('T', 'P', V.p3, 'S', D.S0, D.fluidType) 
        V.rho3 = CP.PropsSI('D', 'P', V.p3, 'T',V.T3, D.fluidType)
        V.h3 = CP.PropsSI('HMASS', 'P', V.p3, 'T',V.T3, D.fluidType)
        V.a3 = CP.PropsSI('A', 'P', V.p3, 'T',V.T3, D.fluidType)
    V.C3 = V.M3 * V.a3
    V.h3t = V.h3 + 1./2. * V.C3 * V.C3 

    # Calculate the geometry properties based on the fluid properties and actual 
    # flow angle at point 3
    #if V.alpha3b > 90:
    #    print('WARNING: alpha3b is larger than 90 degree! The value is: {}'.format(V.alpha3b))
    V.A3devb3 = 2.* np.pi * r3 - Z3 * 0.5* tt * ( 1. + 1./np.cos( np.radians(V.alpha3b) ) )

    #print(1./np.cos( np.radians(V.alpha3b))), 1/0

    #if V.A3devb3 < 0:
    #    print('WARNING! Negative A3devb3 value.')
    #    print( 2.* np.pi * r3,    Z3* 0.5 *tt , ( 1. + 1./np.cos( np.radians(V.alpha3b) )), V.alpha3b )


    # Calculate the flow angle based on the cosine rule
    AB = 2. * r3 * np.sin( np.pi / Z3)
    BC = np.cos( np.radians(V.alpha3b) - np.pi/Z3) *AB
    AC = np.sin( np.radians(V.alpha3b - 0.5*np.pi)) * 2 *r3 * np.sin(np.pi / Z3)
    DC = 0.5 * tt + AC * np.tan( np.radians(D.delta) )
    BE = 0.5 * tt
    V.ds = BC - BE - DC    
    V.s3 = 2* np.pi * r3/ Z3    
    V.alpha3 = np.degrees(np.arccos(V.ds/(V.s3+tt)))  # flow angle is respect to the radial  
    #V.alpha3 = np.degrees( np.radians(V.alpha3b) -   np.pi / Z3  )  
    
   

  
    # calculate radial and tangential velocities at stator exit
    V.C3R = V.C3 * np.cos( np.radians (V.alpha3))
    V.C3T = V.C3 * np.sin( np.radians (V.alpha3))
    # print("C3R, C3T, alpha3",C3R,C3T,V.alpha3)

    # Start to calculate the rotor properties.
    #!!! All the calculated rotor properties are constrained into local value.

    # calculate conserved quantities at rotor inlet:
    C4R_target = D.C4R 
    C4T_target = D.C4T 
    T4_target  = D.T4
    #p4_target = D.p4

    # do forward calculation from state 3 to predict conditions at 4:    
    # Momentum: mdot * V3T * r3 = mdot * V4T * r4
    #C4T = V.C3T * r3/D.r4  
    #Here is considering the momentume loss faction. 
    C4T = V.C3T * r3/ D.r4


    # this can be replaced by momentum equation incorporating losses, but this 
    # also needs modifications to gap_function, which assumes isentropic process.
    C4T_actual = C4T
    


    # calculate actual mass flow rate
    D.b3 = D.b4 # making the blade height is equal to the rotor inlet blade height
    mdot_actual = V.C3R * V.A3devb3*D.b3 * V.rho3 
    
    # need to solve iteratively for other conditions. 
    initial = D.C4R
    maxiter = M.maxiter 
    
    fun = lambda C4R:  gap_function_noloss(C4R,C4T,M,V,D)-mdot_actual
    
    sol = sci.optimize.root(fun, initial, method='hybr',options={'xtol':1.e-8, 'maxfev':maxiter})
    status = sol.status
    C4R = sol.x
    mesg = sol.message
    C4R_actual = C4R

    # calculate rho4 from continuity: C3R * A3 * rho3 = C4R * A4 * rho4 = mdot_actual
    rho4 = mdot_actual / (C4R * D.A4)

    # calcualte velocity
    C4 = np.sqrt(C4R*C4R + C4T*C4T)

    # calculate enthalpy at point 4
    h4 = V.h3t - 1./2. * C4 * C4

    # Use equation of state to get pressure and temperature
    if M.gasModel == 'Ideal':       
        T4 = h4 / D.Cp
        #p4 = T4 * D.R * rho4
    else:
        T4 = CP.PropsSI('T', 'HMASS', h4, 'D',rho4, D.fluidType)
        #p4 = CP.PropsSI('P', 'HMASS', h4, 'D',rho4, D.fluidType)   
    T4_actual = T4
    #p4_actual = p4

    errors = []
    # Calculate Errors
    errors.append(C4T_target - C4T_actual)  # mass flown rate error
    errors.append(float(C4R_target - C4R_actual))  # momentum error
    errors.append(float(T4_target - T4_actual))  # total energy error    
    #errors.append(float(p4_target - p4_actual))  # total energy error  

    if M.verbosity > 0:
        print("Errors:")
        print("C4T_target: {0}; C4T_actual: {1}; delta = {2}".format(C4T_target,C4T_actual,errors[0]))
        print("C4R_target: {0}; C4R_actual: {1}; delta = {2}".format(C4R_target,C4R_actual,errors[1]))
        print("T4_target: {0}; T4_actual: {1}; delta = {2}".format(T4_target,T4_actual,errors[2]))
        #print("P4_target: {0}; P4_actual: {1}; delta = n/a".format(p4_target/1.e3,p4_actual/1.e3))

    return errors

def gap_function_noloss(C4R_guess,C4T,M,V,D):
    # function to caluclate total mass flow at exit of gap for a guessed value of C4R

    C4 = np.sqrt( C4R_guess*C4R_guess + C4T*C4T )
    h4t = V.h3t
    h4 = h4t - 1./2. * C4 * C4
    s4 = D.S0 # assume an isentropic process
    if M.gasModel == 'Ideal':  
        T4 = h4 / D.Cp
        # S = Cp * log (T/D.Tref) - R log(P/D.Pref)  # definition of enthalpy
        P4 = np.exp((D.Cp * np.log(T4/D.Tref) - s4) / D.R) * D.Pref
        rho4 = P4 / (D.R * T4)
    else:
        T4 = CP.PropsSI('T', 'HMASS', h4, 'S', s4, D.fluidType)
        rho4 = CP.PropsSI('D', 'HMASS', h4, 'S', s4, D.fluidType)
    
    mdot4 = C4R_guess * rho4 * D.A4
    return mdot4


def Calculation3(x,M,V,D):
    # Function called by the root optimiser, use constant h0 to connect the stator properties and rotor properties
    # unpack the state vector
    #x = [V.alpha3b, V.M3, V.p3]

    V.alpha3b = x[0]
    V.M3 = x[1]
    V.p3 = x[2]*1.e6
    #print("COME here and alpha3b={0}, M3 = {1}, p3 ={2}".format(x[0],x[1],x[2]))

    if V.p3 < 0:
        print("WARNING!!!!!!,V.p3<0",1/0),
        print("V.p3=",V.p3)
        V.p3 = D.p4
            
    # get geometry constraints, Z and tt
    tt = V.tt
    Z3 = V.Z3
    r3 = V.r3

    if M.verbosity > 0:
        print("\n")
        print("---------------------------")
        print("New Iteration Loop")
        print("X",x)
        print("alpha3b = {0:.3f} [deg]; M3 = {1:.3f} [-]; p3 = {2:.2f} [MPa]; T3 = {3:.4} [K]".format(x[0],x[1],x[2],V.T3))
        print("tt = {0:.5f} [m]; Z3 = {1} [-]".format(tt,Z3))


    # The key step is to assuming there is no energy loss.
    V.h3t = D.h04  # assuming no energy loss

    #  calculate other state properties at point 3/stator outlet
    if M.gasModel == 'Ideal':
        V.T3  = V.h3t/(D.Cp + 0.5 * V.M3**2*D.gamma*D.R)
        #V.T3 = D.Tref * np.exp((D.S0 + D.R * np.log(V.p3/D.Pref)) / D.Cp) 
        V.rho3 = V.p3/ (V.T3 * D.R)        
        V.h3  = D.Cp * V.T3
        V.a3  = np.sqrt(D.gamma * D.R * V.T3)
        V.p03 = V.p3 *   pow(  ( 1. + D.gammaMinusOne/2. * V.M3**2) , (D.gamma/ D.gammaMinusOne) )
        print("V.p03=",V.p03,"D.p04=",D.p04)
        

        V.S3  = D.Cp * np.log(V.T3/D.Tref) - D.R * np.log(V.p3/D.Pref)

    else:
        #V.S3  = CP.PropsSI('S', 'P', V.p3, 'H', V.h3t, D.fluidType)
        #print('p3=',V.p3,'M3=',V.M3)
        initial = V.h3t
        maxiter = M.maxiter 
        fun = lambda h3:  h3_function(h3,V.p3,V.M3,M,V,D) - V.h3t
        
        sol = sci.optimize.root(fun, initial, method='hybr',options={'xtol':1.e-8, 'maxfev':maxiter})
        status = sol.status
        V.h3 = sol.x
        mesg = sol.message

        #print(sol),1/0

        V.T3 = float(CP.PropsSI('T', 'P', V.p3, 'H', V.h3, D.fluidType)) 
        V.rho3 = CP.PropsSI('D', 'P', V.p3, 'T', V.T3, D.fluidType)
        #h30 = CP.PropsSI('H', 'P', V.p3, 'T',V.T3, D.fluidType)
        #print('---', h30 - V.h3),1/0
        #V.h3c = V.h3t - 0.5 * 
        V.a3 = CP.PropsSI('A', 'P', V.p3, 'T',V.T3, D.fluidType)
        V.S3 = CP.PropsSI('S', 'P', V.p3, 'T',V.T3, D.fluidType) 
        V.p03 = CP.PropsSI('P', 'H', V.h3t, 'S', D.S0, D.fluidType )
    V.C3 = V.M3 * V.a3
    #V.h3t = V.h3 + 1./2. * V.C3 * V.C3 

    # Calculate the geometry properties based on the fluid properties and actual 
    # flow angle at point 3
    #if V.alpha3b > 90:
    #    print('WARNING: alpha3b is larger than 90 degree! The value is: {}'.format(V.alpha3b))
    V.A3devb3 = 2.* np.pi * r3 - Z3 * 0.5* tt * ( 1. + 1./np.cos( np.radians(V.alpha3b) ) )

    #print(1./np.cos( np.radians(V.alpha3b))), 1/0

    #if V.A3devb3 < 0:
    #    print('WARNING! Negative A3devb3 value.')
    #    print( 2.* np.pi * r3,    Z3* 0.5 *tt , ( 1. + 1./np.cos( np.radians(V.alpha3b) )), V.alpha3b )


    # Calculate the flow angle based on the cosine rule
    AB = 2. * r3 * np.sin( np.pi / Z3)
    BC = np.cos( np.radians(V.alpha3b) - np.pi/Z3) *AB
    AC = np.sin( np.radians(V.alpha3b - 0.5*np.pi)) * 2 *r3 * np.sin(np.pi / Z3)
    DC = 0.5 * tt + AC * np.tan( np.radians(D.delta) )
    BE = 0.5 * tt
    V.ds = BC - BE - DC    
    V.s3 = 2* np.pi * r3/ Z3    
    V.alpha3 = np.degrees(np.arccos(V.ds/(V.s3+tt)))  # flow angle is respect to the radial  
    #V.alpha3 = np.degrees( np.radians(V.alpha3b) -   np.pi / Z3  )  
    
    # calculate radial and tangential velocities at stator exit
    V.C3R = V.C3 * np.cos( np.radians (V.alpha3))
    V.C3T = V.C3 * np.sin( np.radians (V.alpha3))
    # print("C3R, C3T, alpha3",C3R,C3T,V.alpha3)

    # Start to calculate the rotor properties.
    #!!! All the calculated rotor properties are constrained into local value.

    # calculate conserved quantities at rotor inlet:
    C4R_target = D.C4R 
    C4T_target = D.C4T 
    T4_target  = D.T4
    p4_target  = D.p4  # Use P4 value as an error value to make sure the optimiser covergence.

    # do forward calculation from state 3 to predict conditions at 4:    
    C4T = V.C3T * r3/ D.r4


    # this can be replaced by momentum equation incorporating losses, but this 
    # also needs modifications to gap_function, which assumes isentropic process.
    C4T_actual = C4T
    


    # calculate actual mass flow rate
    D.b3 = D.b4 # making the blade height is equal to the rotor inlet blade height
    mdot_actual = V.C3R * V.A3devb3*D.b3 * V.rho3 

    
    # need to solve iteratively for other conditions. 
    initial = D.C4R
    maxiter = M.maxiter 
    fun = lambda C4R: ( gap_function_withloss(C4R,C4T,M,V,D)-mdot_actual)
    sol = sci.optimize.root(fun, initial, method='hybr',options={'xtol':1.e-8, 'maxfev':maxiter})
    status = sol.status
    C4R = sol.x
    mesg = sol.message
    C4R_actual = C4R
    #print("what is C4R:",C4R,"and D.C4R:",D.C4R,"mdot_actual=",mdot_actual)

    # calculate rho4 from continuity: C3R * A3 * rho3 = C4R * A4 * rho4 = mdot_actual
    rho4 = mdot_actual / (C4R * D.A4)

    # calcualte velocity
    C4 = np.sqrt(C4R*C4R + C4T*C4T)

    # calculate enthalpy at point 4
    h4 = V.h3t - 1./2. * C4 * C4

    # Use equation of state to get pressure and temperature
    if M.gasModel == 'Ideal':       
        T4 = h4 / D.Cp
        p4 = T4 * D.R * rho4
    else:
        T4 = CP.PropsSI('T', 'HMASS', h4, 'D',rho4, D.fluidType)
        p4 = CP.PropsSI('P', 'HMASS', h4, 'D',rho4, D.fluidType)   
    T4_actual = T4
    p4_actual = p4

    errors = []
    # Calculate Errors
    errors.append(C4T_target - C4T_actual)  # mass flown rate error
    errors.append(float(C4R_target - C4R_actual))  # momentum error
    #errors.append(float(T4_target - T4_actual))  # total energy error    
    errors.append(float(p4_target - p4_actual))  # total pressure error

    #print("there is error in error",errors,'T4=',T4)

    if M.verbosity > 0:
        print("Errors:")
        print("C4T_target: {0}; C4T_actual: {1}; delta = {2}".format(C4T_target,C4T_actual,errors[0]))
        print("C4R_target: {0}; C4R_actual: {1}; delta = {2}".format(C4R_target,C4R_actual,errors[1]))
        #print("T4_target: {0}; T4_actual: {1}; delta = {2}".format(T4_target,T4_actual,errors[2]))
        print("P4_target: {0}; P4_actual: {1}; delta = n/a".format(p4_target/1.e3,p4_actual/1.e3))

    return errors

def h3_function(h3_guess,p3,M3,M,V,D):
    T3 = CP.PropsSI('T',  'P', p3, 'H', h3_guess, D.fluidType)
    a3 = CP.PropsSI('A',  'T', float(T3), 'P', p3, D.fluidType)
    C3 = a3 * M3
    h3t = h3_guess + 0.5*C3**2 
    return h3t


def gap_function_withloss(C4R_guess,C4T,M,V,D):
    # function to caluclate total mass flow at exit of gap for a guessed value of C4R

    C4 = np.sqrt( C4R_guess*C4R_guess + C4T*C4T )
    h4t = V.h3t
    p4t = V.p03 *D.totalPressureLossFaction
    h4 = h4t - 0.5* C4**2
    
    s4 = D.S0 # assume an isentropic process
    if M.gasModel == 'Ideal':  

        T4 = h4 / D.Cp
        a4 = np.sqrt( D.gamma * D.R * T4)
        M4 = C4 / a4 
        # S = Cp * log (T/D.Tref) - R log(P/D.Pref)  # definition of enthalpy
        # With losses, we can only considering the effect of total pressure loss here.
        P4 = p4t * pow( ( 1. + 0.5*D.gammaMinusOne *M4**2 )  , (- D.gamma / D.gammaMinusOne) )
        #P4 = np.exp((D.Cp * np.log(T4/D.Tref) - s4) / D.R) * D.Pref
        #print("HER",P4c-P4,"P4c=",P4c,"P4=",P4)
        rho4 = P4 / (D.R * T4)

    else:

        T04 =  CP.PropsSI('T', 'H', h4t, 'P', p4t, D.fluidType)
        S4 =   CP.PropsSI('S', 'T', T04, 'P', p4t, D.fluidType)
        rho4 = CP.PropsSI('D', 'H', h4,  'S', S4, D.fluidType)
    
    mdot4 = C4R_guess * rho4 * D.A4
    #print("In gap_function_withloss, P4=",P4,"rho4=",rho4,'mdot4=',mdot4,"T4=",T4)


    return mdot4


def CalcStatorProps(M,V,D):
    # Functionc to calculate the rest properties for the stator


    # Make the gloable properties locally
    temp = []
    temp.append(V.r3)      # Add the radius
    temp.append(V.Z3)      # Add the stator blade number
    temp.append(V.tt)      # Add the trailing edge thickness
    temp.append(V.M3)      # Add the stator outlet radius 
    temp.append(V.A3devb3) # Add the unified pitch angle length
    temp.append(V.C3)      # Add the absolute velocity
    temp.append(V.alpha3)  # Add the outlet flow angle
    temp.append(V.alpha3b) # Add the stator blade angle
    temp.append(V.C3R)     # Add the radial velocity
    temp.append(V.C3T)     # Add the tangential velocity
    temp.append(V.p3)      # Add the pressure
    temp.append(V.T3)      # Add the temperature
    temp.append(V.rho3)    # Add the density
    temp.append(V.p03)     # Add the total pressure at the stator outlet
    temp.append(D.p04)     # Add the total pressure at the rotor inlet
    temp.append(D.T03)     # Add the total temperature at the stator outlet
    temp.append(D.T04)     # Add the total temperature at the rotor inlet
    return temp

def writeOutput(V,D):

    f = open('output','w')
    f.write('# '+str(datetime.datetime.now())+'\n')
    f.write('# The following are the data list for the stator design.\n')
    f.write('# The design parameters are:\n')
    f.write('# r3devr4min = ' + str(D.r3devr4min) + '\n')
    f.write('# r3devr4max = ' + str(D.r3devr4max) + '\n')
    f.write('# r3devr4step = ' + str(D.r3devr4step) + '\n')    
    f.write('# Z3min = ' + str(D.Z3min) + '\n')
    f.write('# Z3max = ' + str(D.Z3max) + '\n')
    f.write('# Z3step = ' + str(D.Z3step) + '\n')
    f.write('# ttmin = ' + str(D.ttmin) + '\n')
    f.write('# ttmax = ' + str(D.ttmax) + '\n')
    f.write('# ttstep = ' + str(D.ttstep) + '\n')
    f.write('#  [r3]   [Z3]  [tt]     [M3]         [A3devb3]       [C3]        [alpha3]       [alpha3b]         [C3R]         [C3T]        [p3]       [T3]       [rho3]        [p03]        [p04]       [T03]       [T04]\n')
    
    for i in range(len(V.Properties)):
        for j in range(len(V.Properties[i])):
            f.write(str(V.Properties[i][j]) + ' ')
        f.write('\n')
   

    print("I'm writing the lovely resutls ...")
    f.close()
    return 0




###
###
class Model:
    def __init__(self):
        self.test = 0
        self.iteration = 0
    ##
    def check(self):
        if not self.test == 0:
            raise MyError('M.testincorrect not specified')
        #if self.maxiter == 0:
        if not (self.gasModel == 'Ideal' or self.gasModel == 'Real'):
            raise MyError('You do not correctly set the gas model type.')
        if not (self.optimiser == 'Root' or self.optimiser == 'Nelder-Mead'):
            raise MyError('You do not correctly set the optimiser solver type.')
        # add additional checks to assess that basics have been set-up correct;
###
###
class Variables:
    def __init__(self):
        self.MasterCase_folder = []
        self.Stator_number = []
        self.Properties = []
        #self.constant = []
        #self.variables = []
        #self.ChildCase_folder_base = []

    ##
    def check(self):
        pass
        #if not self.MasterCase_folder:
        #    raise MyError('C.MasterCase_folder not specified')
        #if not self.ChildCase_folder_base:
        #    raise MyError('C.ChildCase_folder_base not specified')
###
###
class Data:
    def __init__(self):
        pass
        #self.Counter = 0
        #self.Use_initial_simplex = []
        #self.initial_simplex_filename = []

    ##
    def check(self):
        pass
        #if not (self.Use_initial_simplex == True or self.Use_initial_simplex == False):
        #    raise MyError('D.Use_initial_simplex has to be True or False')


###
###
def main(uoDict):
    """
    main function
    """
    # create string to collect warning messages
    warn_str = "\n"


    # main file to be executed
    jobFileName = uoDict.get("--job", "test")

    # strip .py extension from jobName
    jobName = jobFileName.split('.')[0]

    # create classes to store input data
    M = Model() # data definig optimisation run.
    V = Variables()  # Data that can be changed during calculation
    D = Data()      # Data that unchanged during calculation


    # set verbosity (can be overwritten from jobfile)
    M.verbosity = 0
    if "--verbosity" in uoDict:
        M.verbosity = int(uoDict.get("--verbosity", 0))

    # set ROOT 
    M.ROOT_DIR = os.getcwd()
    print(jobFileName)
    # Execute jobFile, this creates all the variables
    exec(open(M.ROOT_DIR+'/'+jobFileName).read(),globals(),locals())  

    # check that input data has been set correctly.
    # Currently empty
    M.check()
    V.check()
    D.check()    

    # Use the rotor target value as initial value
    CalInitial(M,V,D)


    #------------------------------------------------------------------------------------------------------------------#
    # test case to analyse operation of solver
    if False:
        # set trailing edge thickness and blade number
        V.tt = 1e-3
        V.Z3 = 15
        V.r3 = 0.0040 
        
        # Do single evaluation
        A0 = (D.alpha4,D.M4,D.p4) # use conditions at 4 as estimates to get conditions at 3.
        args = (M, V, D)
        
        max_fev = 100

        sol = sci.optimize.root(Calculation2,A0,args=args,method='lm',options={'eps':1.e-3, 'xtol':1.e-12, 'maxiter':max_fev})
        status = sol.status
        X = sol.x
        fun = sol.fun
        mesg = sol.message
        
        P3 = X[2]    
        if M.gasModel == 'Ideal':  
            T3 = D.Tref * np.exp((D.S0 + D.R * np.log(P3/D.Pref)) / D.Cp) 
        else:
            V.T3 = CP.PropsSI('T', 'P', V.p3, 'S',D.S0, D.fluidType) 
            
        print("\n")
        print("++++++++++++++++++++++++++++")
        print("Solver Status:",mesg)
        print("ERRORS: C4T={0:1.2e} [m/s]; C4R={1:1.2e} [m/s]; T4={2:1.2e} [K]; P4=n/a [kPa]".format(fun[0],fun[1],fun[2]))
        print("SOLUTION:",X)          #S = [V.alpha3b, V.M3, V.p3]
        print("alpha3b = {0:.3f} [deg]; M3 = {1:.3f} [-]; p3 = {2:.1f} [kPa]; T3 = {3:.4} [K]".format(X[0],X[1],X[2]/1.e3,T3))
        print("tt = {0:.5f} [m]; Z3 = {1} [-]".format(V.tt,V.Z3))
        print("++++++++++++++++++++++++++++")   
        print("Conditions at rotor")
        print("alpha4 = {0:.3f} [deg]; M4 = {1:.3f} [-]; p4 = {2:.1f} [kPa]; T4 = {3:.4} [K]".format(D.alpha4,D.M4,D.p4/1.e3,D.T4))
            
        P3 = X[2]; M3 = X[1]
        print("Checks:")
        T03 = T3 * (1. + 0.5* D.gammaMinusOne * M3**2)
        P03 = P3 * (1. + 0.5* D.gammaMinusOne * M3**2)**(D.gamma / D.gammaMinusOne) 
        print("Total Temperature: T03={0:1.2e} [K]; T04={1:1.2e} [kPa]".format(T03,D.T04))
        print("Total Pressure: P03={0:1.2e} [kPa]; P04={1:1.2e} [kPa]".format(P03,D.p04))
        print("Mach No: M3={0:1.2e} [-]; M4={1:1.2e} [-]".format(M3,D.M4))
        print("\n")
        #a=a+5
    #------------------------------------------------------------------------------------------------------------------#


    #Initialize the trailing edge thickness range, stator blade number range and the radius ratio range
    ttlist = np.arange(D.ttmin, D.ttmax, D.ttstep)

    Zlist  = np.arange(D.Z3min, D.Z3max, D.Z3step)
    r3devr4list = np.arange(D.r3devr4min, D.r3devr4max , D.r3devr4step)
    #TODO: why add 1 step here. DONE


    if M.verbosity >= 1:
        print("ttlist = {}".format(ttlist))
        print("Zlist = {}".format(Zlist))
        print("r3devr4list={}".format(r3devr4list))

    
    for i in range(len(r3devr4list)):
        for j in range(len(Zlist)):
            for k in range(len(ttlist)):
            
                # initialise the temp list to hold the properties
                tempProps = []

                V.r3 = r3devr4list[i] * D.r4
                V.Z3 = Zlist[j]
                V.tt = ttlist[k]
                

                if 'Root' == M.optimiser:
                    

                    print("\n \n")
                    print("+++++++++++++++++++++")
                    print("+++++++++++++++++++++")
                    print("Solving for tt = {0:.3f} [mm]; Z3 = {1} [-]; r3 = {2} [mm]".format(V.tt/1.e-3,V.Z3, V.r3/1.e-3))
                    print("Creating the simplex, alpha4 = {0}; M4 = {1}; p4 = {2} [MPa]".format(D.alpha4, D.M4, D.p4/1.e6))
                    A0 = (D.alpha4,D.M4,D.p4/1.e6) # use conditions at 4 as estimates to get conditions at 3.
                    args = (M, V, D)
                    max_fev = M.maxiter
                    print("p3=",V.p3)
                    if D.totalPressureLossFaction == 1.:
                        sol = sci.optimize.root(Calculation2,A0,args=args,method='lm',options={'eps':1.e-3, 'xtol':1.e-12, 'maxiter':max_fev})
                        status = sol.status
                        X = sol.x
                        fun = sol.fun
                        mesg = sol.message

                        P3 = X[2] *1.e6    
                        print("---",P3)
                        if M.gasModel == 'Ideal':  
                            T3 = D.Tref * np.exp((D.S0 + D.R * np.log(P3/D.Pref)) / D.Cp) 
                        else:
                            T3 = CP.PropsSI('T', 'P', V.p3, 'S',D.S0, D.fluidType) 


                    else:
                        sol = sci.optimize.root(Calculation3,A0,args=args,method='hybr',options={'eps':1.e-3, 'xtol':1.e-12, 'maxiter':max_fev})
                        status = sol.status
                        X = sol.x
                        fun = sol.fun
                        mesg = sol.message

                        P3 = X[2]*1.e6 


                        # TODO: to correct the T3 calculations along with the new one
                        #if M.gasModel == 'Ideal':  
                        #    T3 = D.Tref * np.exp((D.S0 + D.R * np.log(P3/D.Pref)) / D.Cp) 
                        #else:
                        #    T3 = CP.PropsSI('T', 'P', V.p3, 'S',D.S0, D.fluidType) 
                        T3 = V.T3


            
                    print("\n")
                    print("++++++++++++++++++++++++++++")
                    print("Solver Status:",mesg)
                    print("ERRORS: C4T={0:1.2e} [m/s]; C4R={1:1.2e} [m/s]; T4={2:1.2e} [K]; P4=n/a [kPa]".format(fun[0],fun[1],fun[2]))
                    print("SOLUTION:",X)          #S = [V.alpha3b, V.M3, V.p3]
                    print("alpha3b = {0:.3f} [deg]; M3 = {1:.3f} [-]; p3 = {2:.3f} [MPa]; T3 = {3:.4} [K]".format(X[0],X[1],X[2],T3))
                    print("tt = {0:.5f} [m]; Z3 = {1} [-]".format(V.tt,V.Z3))
                    print("++++++++++++++++++++++++++++")   

                    if V.A3devb3 < 0:
                        print('WARNING!!!!!!',V.A3devb3, X[0],), 1/0

                    V.Properties.append(CalcStatorProps(M,V,D))


                elif "Nelder-Mead" == M.optimiser:
                    #TODO: make Nelder-Mead optimiser available

                    # Here we use Nelder-Mead method to grab the gas and geometric properties.

                    # Two variables are optimised.
                    x0 = [D.alpha4,D.M4,D.p4/1.e6]

                    V.initial_simplex = initSimplex(x0)                
                               

                    ############################################################
                    # THE OPTIMISER!
                    # Here we use the minimize module from scipy.py
                    # The method is Nelder-Mead.
                    # For a given case, the initial simplex need to be adjust 
                    # from the input.py file, which located in the example folder
                    res =  minimize(Calculation3,x0,args=(M,V,D),method='Nelder-Mead',options={'initial_simplex': V.initial_simplex, 'maxiter':M.maxiter})
                    #
                    #
                    ############################################################

                    # TODO: need a function to calculate the rest geometry and fluid properties for stator
                    
                    print(res)
                    
                    
                    #tempProps.append(CalcStatorProps(M,V,D))
                    V.Properties.append(CalcStatorProps(M,V,D))

                    #CalcStatorProps(M,V,D)

                    
                    print('OPTIMISATION COMPLETE')
                    print('    Final Residual: {}'.format(res))
                    print('    alpha3b is:{}'.format(res.x[0]))
                    print('    M3 is:{}'.format(res.x[1]))
            
                else:
                    print("ERROR! PLEASE ENTER A CORRECT OPTIMISATION METHOD.",1/0)

    # write the output
    writeOutput(V,D)
    
    return 0
###
###
class MyError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
###
###
shortOptions = ""
longOptions = ["help", "job=", "verbosity="]
###
###
def printUsage():
    print("")
    print("Usage: Optimizer.py [--help] [--job=<jobFileName>] [--verbosity=<0,1,2>]")
    print("\n")
    print(" --help      Display help.")
    print("\n")
    print(" --job=      Use this to specify the job file.")
    print("\n")
    print(" --verbosity   Set level of screen output 0-none; 1-some; 2-all.")
    return
###
###
if __name__ == "__main__":

    userOptions = getopt(sys.argv[1:], shortOptions, longOptions)
    uoDict = dict(userOptions[0])

    if len(userOptions[0]) == 0 or "--help" in uoDict:
        printUsage()
        sys.exit(1)

    # execute the code
    try:
        main(uoDict)
        print("\n \n")
        print("SUCCESS.")
        print("\n \n")

    except MyError as e:
        print("This run has gone bad.")
        print(e.value)
        sys.exit(1)

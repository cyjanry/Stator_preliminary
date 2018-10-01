#!/usr/bin/env python
# This Python file uses the following encoding: utf-8



import os         as         os
import numpy      as         np
import shutil     as         sh
from   getopt            import getopt
import sys        as         sys 
from   scipy.interpolate import griddata
import scipy      as sci
from   scipy.optimize    import minimize
import matplotlib.pyplot as plt 
import CoolProp.CoolProp as CP




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

    D.gammaMinusOne = D.gamma - 1. 

    # Gas constant
    Ru =  CP.PropsSI(D.fluidType, "gas_constant")
    W  =  CP.PropsSI(D.fluidType, "molar_mass")
    D.R  = Ru/W
    # Calculate corresponding rotor fluid properties.

    # In the following section, we apply the real gas model
    #a4  = CP.PropsSI('A', 'T', D.T4, 'P', D.p4 , D.fluidType) # real gas model
    D.a4    = np.sqrt(D.gamma * D.R * D.T4) # real gas model, acoustic speed.
    D.M4    = D.C4/D.a4
    D.T04   = D.T4 * (1. + 0.5* D.gammaMinusOne * D.M4**2 )
    D.p04   = D.p4 * pow(  (1. + 0.5* D.gammaMinusOne * D.M4**2), (D.gamma / D.gammaMinusOne)  )
    D.rho4  = D.p4/ (D.R*D.T4)
    D.rho04 = D.rho4 * pow(  (1. + 0.5* D.gammaMinusOne * D.M4**2), (1 / D.gammaMinusOne)  )
    D.A4devb4 = D.A4/D.b4
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
    D.rho03 = D.rho04

    # Assuming there is no loss in momentum:
    D.C3T = D.C4T * D.r4 / D.r3

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
        temp[i] = temp[i] + 0.01*temp[i]
        simplex.append(temp)

    return np.array(simplex)



def M4Func(M4temp,M,V,D,C3R,rho3,A3devb3,M3,C4T,T3):
    #
    left  = ( C3R * rho3 *A3devb3 / (D.A4devb4 * rhoFunc(rho3, M3, M4temp,D.gamma) ) ) **2 +  C4T**2
    right =  M4temp**2 * D.gamma * D.R * tFunc(T3,M3,M4temp,D.gamma)
    
    cost = abs(left - right)

    return cost

def Calculation(x,M,V,D):

    x0 = x

    #x0 = [V.alpha3b, V.tt, V.M3]

    alpha3b = x[0]
    tt = x[1]
    M3 = x[2]

    T3 = D.T03 / (1. + 0.5* D.gammaMinusOne * M3**2 )
    p3 = D.p03 *  pow(  (1. + 0.5* D.gammaMinusOne * D.M4**2), (-D.gamma / D.gammaMinusOne)  )
    a3 = np.sqrt( D.gamma * D.R * T3 )
    rho3 = p3/ (T3 * D.R)

    C3 = M3 * a3

    A3devb3 = 2.* np.pi * D.r3 - D.Z3 * 0.5* tt * ( 1. - 1./np.cos( np.radians(alpha3b) ) )

    AB = 2. * D.r3 * np.sin( np.pi / D.Z3)
    BC = np.cos( np.radians(alpha3b) - np.pi/D.Z3) *AB
    AC = np.sin( np.radians(alpha3b - 0.5*np.pi)) * 2 *D.r3 * np.sin(np.pi / D.Z3)
    DC = 0.5 * tt + AC * np.tan( np.radians(D.delta) )
    BE = 0.5 * tt
    ds = BC - BE - DC    
    s3 = 2* np.pi * D.r3/ D.Z3    
    alpha3 = np.degrees(np.arccos(ds/s3))  # flow angle is respect to the radial


    C3R = C3 * np.cos( np.radians (alpha3))
    C3T = C3 * np.sin( np.radians (alpha3))

    # Then use the calculated stator properties to calculate the rotor properties.
    C4T = D.r4/ D.r3 * C3T
    T04 = D.T03
    P04 = D.p03 




    # first we calculate M4 
    # Based on the continuity:    
    #C4R  = C3R * rho3 *A3devb3 / (D.A4devb4 *rho4)
    #C4R**2 + C4T**2 = M4**2 * D.gamma * D.R * T4 
    #assuming M4 to make left = right
    # need to apply an optimiser to find the correct solution

    #M4temp = np.linspace(0.,1.0, 10000)
    #for i in range(len(M4temp)):
    #    left  = ( C3R * rho3 *A3devb3 / (D.A4devb4 * rhoFunc(rho3, M3, M4temp[i],D.gamma) ) ) **2 +  C4T**2
    #    right =  M4temp[i]**2 * D.gamma * D.R * tFunc(T3,M3,M4temp[i],D.gamma)
    #    if abs(left -right) < 10:
    #        M4inter =  M4temp[i]
    #        #print M4inter
    #        break
    #M4sectemp = np.linspace( M4inter-0.01, M4inter+0.01, 100000)
    #for i in range(len(M4sectemp)):
    #    left  = ( C3R * rho3 *A3devb3 / (D.A4devb4 * rhoFunc(rho3, M3, M4sectemp[i],D.gamma) ) ) **2 +  C4T**2
    #    right =  M4sectemp[i]**2 * D.gamma * D.R * tFunc(T3,M3,M4sectemp[i],D.gamma)
    #    if abs(left -right) <0.1:
    #        M4 = M4sectemp[i]
    #        break
    #Hence
    #T4 = T04 / (1. + 0.5* D.gammaMinusOne * D.M4**2 )
    #p4 = p04 * pow(  (1. + 0.5* D.gammaMinusOne * D.M4**2), (-D.gamma / D.gammaMinusOne)  )
    #x0 = [V.alpha3b, V.tt, V.M3]
    #V.initial_simplex = initSimplex(x0)
    #res =  minimize(Calculation,x0,args=(M,V,D),method='Nelder-Mead',options={'initial_simplex': V.initial_simplex, 'maxiter':M.maxiter})
    M4temp0 = [0.8] # initial value for Mach4
    init_simp = initSimplex(M4temp0)
    result = minimize(M4Func, M4temp0, args=(M,V,D,C3R,rho3,A3devb3,M3,C4T,T3), method='Nelder-Mead',options={'initial_simplex': init_simp, 'maxiter':M.maxiter})
    M4 = result.x[0]
 
    if M.verbosity >=1 :
        print("M3 is:", M3)
        print("alpha3b is:", alpha3b)
        #print("AB is:", AB)
        #print("BC is:", BC)
        #print("BE is:", BE)
        #print("DC is:", DC)
        #print("ds is:", ds)
        #print("s3 is:", s3)
        #print("alpha3 is:",alpha3)
        print("M4 is:", M4)
        print("----------")



    cost = abs(M4 - D.M4)

    return cost 
    


###
###
class Model:
    def __init__(self):
        self.test = 0
    ##
    def check(self):
        if not self.test == 0:
            raise MyError('M.testincorrect not specified')
        # add additional checks to assess that basics have been set-up correct;
###
###
class Variables:
    def __init__(self):
        self.MasterCase_folder = []
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
    M.verbosity = 1
    if "--verbosity" in uoDict:
        M.verbosity = int(uoDict.get("--verbosity", 1))

    # set ROOT 
    M.ROOT_DIR = os.getcwd()
    print(jobFileName)
    # Execute jobFile, this creates all the variables
    exec(open(M.ROOT_DIR+'/'+jobFileName).read(),globals(),locals())  

    # check that input data has been set correctly.
    M.check()
    V.check()
    D.check()    

    # Use the rotor target value as initial value
    CalInitial(M,V,D)
    
    x0 = [V.alpha3b, V.tt, V.M3]

    V.initial_simplex = initSimplex(x0)
 

    

    

    ############################################################
    # THE OPTIMISER!
    # Here we use the minimize module from scipy.py
    # The method is Nelder-Mead.
    # For a given case, the initial simplex need to be adjust 
    # from the input.txt file, which located in the example folder
    res =  minimize(Calculation,x0,args=(M,V,D),method='Nelder-Mead',options={'initial_simplex': V.initial_simplex, 'maxiter':M.maxiter})
    #
    #
    ############################################################

    print('OPTIMISATION COMPLETE')
    print('    Final Residual: {}'.format(res))
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

"""
reservoir simulation project 1(2022)
2D Multiphase reservoir simulation: Fluid properties file
Author: Promise O. Longe
Email: longe.promise@ku.edu
Date modified: 04/18/2022
"""
import numpy as np
from scipy import interpolate as interp
import matplotlib.pyplot as plt

class fluid:
    def __init__(self):
        self.mu = []
class reservoir:
    def __init__(self):
        self.dt = [] 

#fluid, reservoir parameters

# Pressure (psia)| Water FVF (RB/bbl)| Oil FVF (RB/bbl)| Gas FVF (RB/scf)| Water(lbm/ft3)|
# Oil(lbm/ft3)| Gas(lbm/ft3)| Water(cp)| Oil(cp)| Gas(cp)| Rs(scf/STB)| z-factor
pvt = np.loadtxt('pvt_prop.dat')
'''
def fluid_properties(reservoir, fluid, P, Pw):
    """fluid properties estimation using interpolation function"""
    
    f_Bo = interp.interp1d(pvt[:,0], pvt[:,2], kind = 'cubic')
    f_Bw = interp.interp1d(pvt[:,0], pvt[:,1], kind = 'cubic')
    f_rho = interp.interp1d(pvt[:,0], pvt[:,5], kind = 'cubic')
    f_rhw = interp.interp1d(pvt[:,0], pvt[:,4], kind = 'cubic')
    f_muo = interp.interp1d(pvt[:,0], pvt[:,8], kind = 'cubic')
    f_muw = interp.interp1d(pvt[:,0], pvt[:,7], kind = 'linear')
        
    fluid.Bw = f_Bw(Pw)
    fluid.Bo = f_Bo(P)
    fluid.rhoo = f_rho(P)
    fluid.rhow = f_rhw(P)
    fluid.muo = f_muo(P)
    fluid.muw = f_muw(Pw)

    return 
'''
def fluid_properties(reservoir, fluid, P, Pw):

    fluid.Bw = 1.00 / (1.0 + fluid.cw * (P - 14.7))
    fluid.Bo = 1.00 / (1.0 + fluid.co * (P - 14.7))
    fluid.rhgsc = (28.9586 * fluid.sg / 379.4)
    fluid.rho = (fluid.rhosc + fluid.Rs * fluid.rhgsc / 5.615) / fluid.Bo
    fluid.rhw = (fluid.rhwsc + fluid.Rs * fluid.rhgsc / 5.615) / fluid.Bw
    fluid.muo = 0.981 - 3.73E-05 * P + 5.0E-09 * P ** 2
    #fluid.rho = fluid.rho_ref * (1.0 + fluid.co * (P - fluid.BP))
    #fluid.rhw = fluid.rhw_ref * (1.0 + fluid.cw * (P - fluid.BP))
    fluid.muw = 0.52 * P/P
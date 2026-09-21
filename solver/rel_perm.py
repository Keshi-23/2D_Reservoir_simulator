"""
reservoir simulation project 1(2022)
2D Multiphase reservoir simulation - Capillary pressure and relative permeability: Rel perm
Author: Promise O. Longe
Email: longepromise@ku.edu
Date modified: 04/18/2022
"""
import numpy as np
import scipy.interpolate as interp
import matplotlib.pyplot as plt

#-----------------------------------------------------------------------------------------------------------------------------------------------------

class petro:
    def __init__(self):
        self.Swr = []

rel = np.loadtxt('rel_perm_data.txt')

def rel_perm(petro, Sw, rel=rel):
    """rel perm estimation using interpolation function"""

    f_krw = interp.interp1d(rel[:,0], rel[:,1], kind = 'cubic')
    f_kro = interp.interp1d(rel[:,0], rel[:,2], kind = 'cubic')

    S = (Sw - petro.Swr) / (1.0 - petro.Swr)  # Normalized saturation

    if Sw >= 1-petro.Sor:
        krw = petro.krw0
        kro = 0.0
    elif Sw <= petro.Swr:
        krw = 0.0
        kro = petro.kro0
    else:
        krw = f_krw(Sw)
        kro = f_kro(Sw)

    return np.array([krw, kro], dtype=np.float64)
#.----------------------------------------------------------------------------------------------------------------------------
'''

def rel_perm(petro, Sw):
    """s: saturation, s_wp: wetting percolation threshold, s_nwp: non-wetting percolation threshold,
        mu: viscosity, k_0: relative permeability threshold, n: power law parameter"""

    if Sw >= 1 - petro.Swr:
        krw = petro.krw0
        kro = 0.0
    elif Sw <= petro.Swr:
        krw = 0.0
        kro = petro.kro0
    else:
        S = (Sw - petro.Swr) / (1.0 - petro.Sor - petro.Swr)  # Normalized saturation

        krw = (petro.krw0 * S ** petro.nw)  # Corey-Brooks model
        kro = (petro.kro0 * (1.0 - S) ** petro.no)  # Corey-Brooks model

    return np.array([krw, kro], dtype=np.float64)
'''
'''
reservoir simulation project 1(2022)
2D Multiphase reservoir simulation: Capillary Pressure file
Author: Promise O. Longe
Email: longe.promise@ku.edu
Date modified: 04/18/2022
'''

import numpy as np
from scipy import interpolate as interp
import matplotlib.pyplot as plt

class petro:
    def __init__(self):
        self.Sor = []

cap = np.loadtxt("cap.txt")
def cap_press(petro, Sw, cap=cap):

    f_Pc = interp.interp1d(cap[:,0], cap[:,1], kind = 'cubic')

    S = (Sw - petro.Swr) / (1.0 - petro.Swr)  # Normalized saturation

    if np.where(Sw >= 1-petro.Sor):
        Pc = 0.0
        Pcprime = 0.0

    elif np.where(Sw <= petro.Swr):
        Pc = 0.0
        Pcprime = 0.0
    else:
        Pc = f_Pc(Sw)
        S1 = Sw - 0.1
        S2 = Sw + 0.1

        if np.where(S2 >= 1 - petro.Sor) or np.where(S1 <= petro.Swr):
            Pcprime = 0.0
        else:
            P1 = f_Pc(S1)
            P2 = f_Pc(S2)
            Pcprime = P2 - P1 / 0.1

    return Pc, Pcprime

def cap_press_init(petro, Sw, cap=cap):

    f_Pc = interp.interp1d(cap[:,0], cap[:,1], kind = 'cubic')

    S = (Sw - petro.Swr) / (1.0 - petro.Swr)  # Normalized saturation

    if Sw[np.argwhere(Sw >= 1-petro.Sor)]:
        Pc = 0.0
        Pcprime = 0.0

    elif Sw[np.argwhere(Sw <= petro.Swr)]:
        Pc = 0.0
        Pcprime = 0.0
    else:
        Pc = f_Pc(Sw)
        S1 = Sw - 0.1
        S2 = Sw + 0.1

        if S2[np.argwhere(S2 >= 1-petro.Sor)] or S1[np.argwhere(S1 <= petro.Swr)]:
            Pcprime = 0.0
        else:
            P1 = f_Pc(S1)
            P2 = f_Pc(S2)
            Pcprime = P2 - P1 / 0.1

    return Pc, Pcprime

'''
def cap_press(petro, Sw):

    # For initialization
    Se = (Sw - petro.Swr) / (1.0 - petro.Swr)  # {Eq. 1.28a}

    # Corey-Brooks model
    Pc = petro.Pe * (Se ** (-1.0 / petro.lamda))  # Drainage capillary pressure (used for initialization) {Eq. 1.28a}
    Pcprime = (-petro.Pe / petro.lamda) * Se ** ((-1 / petro.lamda)-1)

    return Pc, Pcprime
'''
# -------------------------------------------------------------------------------------------------------------------------------

'''
def cap_press(petro, Sw, Sw_hyst):

    # For initialization
    S = (Sw - petro.Swr) / (1.0 - petro.Sor - petro.Swr)    # Normalized saturation {Eq. 1.29}
    Se = (Sw - petro.Swr) / (1.0 - petro.Swr)               # {Eq. 1.28a}

    # Corey-Brooks model
    Pcd = petro.Pe * (Se ** (-1.0 / petro.lamda))           # Drainage capillary pressure (used for initialization) {Eq. 1.28a}
    Pci = petro.Pe * (S ** (-1.0 / petro.lamda) - 1.0)      # Imbibition capillary pressure (used for initialization) {Eq. 1.28b}

    # Capillary pressure scanning curve
    epspc = 1E-01
    Sw_max = 1.0 - petro.Sor
    f = ((Sw_max + epspc) / (Sw_max)) * ((Sw) / (Sw + epspc))
    # f = 1.0
    Pc = f * Pci + (1.0 - f) * Pcd

    # Calculate derivative
    S2 = (Sw + 0.001 - petro.Swr) / (1.0 - petro.Swr - petro.Sor)
    Se2 = (Sw + 0.001 - petro.Swr) / (1.0 - petro.Swr)

    Pcd2 = petro.Pe * Se2 ** (-1.0 / petro.lamda)
    Pci2 = petro.Pe * (S2 ** (-1.0 / petro.lamda) - 1.0)
    f2 = ((Sw_max - Sw_hyst + epspc) / (Sw_max - Sw_hyst)) * ((Sw + 0.001 - Sw_hyst) / (Sw + 0.001 - Sw_hyst + epspc))
    # f2 = 1.0
    Pc2 = f2 * Pci + (1.0 - f2) * Pcd

    Pcprime = (Pc2 - Pc) / 0.001

    return Pc, Pcprime
'''
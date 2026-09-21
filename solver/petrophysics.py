"""
reservoir simulation project 1(2022)
2D Multiphase reservoir simulation - Capillary pressure and relative permeability: Petrophysics main function
Author: Promise O. Longe
Email: longepromise@ku.edu
Date modified: 04/18/2022
"""
#import inbuilt libraries
import numpy as np
from scipy import interpolate as interp
import matplotlib.pyplot as plt

#importing personal libraries
from input_file_2D import inputfile
from rel_perm import rel_perm
from petroplots import petroplots
from cap_press import cap_press

# making petrophysics class
class fluid:
    def __init__(self):
        self.mu = []
class numerical:
    def __init__(self):
        self.Bw  = []
class reservoir:
    def __init__(self):
        self.dt = []
class grid:
    def __init__(self):
        self.xmin = []
class BC:
    def __init__(self):
        self.xmin = []
class IC:
    def __init__(self):
        self.xmin = []
class well:
    def __init__(self):
        self.xmin = []
        self.xblock = []
class petro:
    def __init__(self):
        self.Swr = []

inputfile(fluid, reservoir, petro, numerical, BC, IC, well)     #uploading petrophysical properties

Sw = np.transpose([np.linspace(petro.Swr,1-petro.Sor,10000)]) #initializing water saturation as a column vector
#Sw_actual = np.transpose([np.linspace(petro.Swr,1-petro.Sor,10000)]) #initializing water saturation

krw = np.zeros((len(Sw),1)) #allocating column vector for relative permeability (water)
kro = np.zeros((len(Sw),1)) #allocating column vector for relative permeability (oil)
Pci = np.zeros((len(Sw),1)) #allocating column vector for imbibition capillary pressure [psi]
Pcd = np.zeros((len(Sw),1)) #allocating column vector for drainage capillary pressure  [psi]


#evaluating relative permeabilities using Corey-Brooks model
for j in range(0, len(Sw)):
    Sw_hyst = Sw[j, 0]
    [krw[j], kro[j]] = rel_perm(petro, Sw[j, 0])
    [Pci[j], Pcd[j]] = cap_press(petro,Sw[j,0], Sw_hyst)

#plotting
petroplots(Sw,krw,kro,Pci, Pcd)
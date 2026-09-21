"""
reservoir simulation project 1(2022)
2D Multiphase reservoir simulation: Input file
Author: Promise O. Longe
Email: longepromise@ku.edu
Date modified: 04/18/2022
"""
class fluid:
    def __init__(self):
        self.mu = []

class numerical:
    def __init__(self):
        self.Bw = []

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

class petro:
    def __init__(self):
        self.Swr = []

class well:
    def __init__(self):
        self.xmin = []
        self.xblock = []

import numpy as np
import matplotlib.pyplot as plt
from Thalf import Thalf
from cap_press import cap_press
from myarrays_2 import myarrays
from rel_perm import rel_perm
from fluid_properties import fluid_properties
from rock_properties import rock_properties
import warnings

warnings.filterwarnings("ignore")

#reading files

dx = np.loadtxt('block_x.dat')
dy = np.loadtxt('block_y.dat')
pvt = np.loadtxt('pvt_prop.dat')
rel = np.loadtxt('rel_perm_data.txt')
cap = np.loadtxt("cap.txt")
perm_x = np.loadtxt('perm_x.dat')
perm_y = np.loadtxt('perm_y.dat')
poro = np.loadtxt('poro.dat')
depth = np.loadtxt('depth.dat')
thickness = np.loadtxt('thickness.dat')

#reading the grid flag
flag = np.loadtxt('flag.txt')
dx_old = np.copy(dx)
dy_old = np.copy(dy)

# Adding ghost grid
dx = np.pad(dx, [(1,1),(1,1)], mode='constant', constant_values = [(500,500), (500,500)])
dy = np.pad(dy, [(1,1),(1,1)], mode='constant', constant_values = [(500,500), (500,500)])
perm_x = np.pad(perm_x, [(1,1),(1,1)], mode='constant', constant_values = [(0,0), (0,0)])
perm_y = np.pad(perm_y, [(1,1),(1,1)], mode='constant', constant_values = [(0,0), (0,0)])
poro = np.pad(poro, [(1,1),(1,1)], mode='constant', constant_values = [(0,0), (0,0)])
depth = np.pad(depth, [(1,1),(1,1)], mode='constant', constant_values = [(0,0), (0,0)])
thickness = np.pad(thickness, [(1,1),(1,1)], mode='constant', constant_values = [(0,0), (0,0)])

#fluid, reservoir and simulation parameters   
def inputfile(fluid,reservoir,petro,numerical,BC,IC,well):
    # Numerical simulation parameters
    numerical.dt     = 0.2 #float(input('Enter timestep: ', ))  #time step (days)
    numerical.tfinal = 10 #final time [days]
    numerical.PV_final = 2  #final pore volume
    numerical.Nx     = len(dx[0,:])  #number of grid blocks in x-direction
    numerical.Ny     = len(dy[:,0])   #number of grid blocks in y-direction
    numerical.N  = numerical.Nx * numerical.Ny #Total number of grid blocks
    numerical.method = 'IMPES' #Implicit pressure explicit saturation solver
    
    # Fluid parameters
    fluid.muw = np.ones((numerical.N, 1))   #fluid viscosity [centipoise]
    fluid.Bw  = np.ones((numerical.N, 1))   #formation volume factor of water [rb/stb]
    fluid.cw  = 1.0E-6                      #total compressibility: rock [1/psi] 

    fluid.muo = np.ones((numerical.N, 1)) #fluid viscosity [centipoise]
    fluid.Bo  = np.ones((numerical.N, 1)) #formation volume factor of oil [rb/stb]
    fluid.co  = 5.0E-6     #oil compressibility: oil fluid [1/psi]
    fluid.ct  = 11.9E-6     #total compressibility: rock + fluid [1/psi]

    fluid.sg = 0.7              # SG of the gas
    fluid.Rs = 500              # GOR of the gas
    fluid.rhwsc   = 62.4        # density of the water [lbm/ft^3]
    fluid.rhosc = 52.40          # density of the stock tank oil [lbm/ft^3]
    fluid.BP     = 5000         # bubble point pressure [psi]
    fluid.Bo_ref   = 1.0        # formation volume factor at bubble point pressure [rb/stb]
    fluid.Bw_ref   = 1.0        # formation volume factor at bubble point pressure [rb/stb]
    fluid.muo_ref = 0.92        # viscosity constant
    fluid.Pref = 14.7           # Reference pressure
    
    #fluid.relperm = 1.0*np.ones((numerical.N, 1))  #Relative permeability of the fluid

    # Multiphase/Relative permeability values and Capillary pressure
    petro.Swr  = 0.18            #residual water saturation
    petro.Sor  = 0.1            #residual oil saturation
    petro.nw   = 3.866771       #Corey-Brooks exponent (water)
    petro.no   = 1.10423        #Corey-Brooks exponent (oil)
    petro.krw0 = 0.59439        #Corey-Brooks endpoint (water).
    petro.kro0 = 1.0000         #Corey-Brooks endpoint (oil)
    petro.lamda= 1.5202         #fitting parameter for Corey-Brooks model
    petro.Pe   = 0.7171         #capillary entry pressure [psi]

    # Reservoir parameters
    reservoir.L  = np.sum(dx[0,:])      #length of the reservoir [feet]
    reservoir.h  = 22                   # np.reshape(thickness, (numerical.N,1)) #height of the reservoir [feet]
    reservoir.W  = np.sum(dy[0,:])      #width of the reservoir [feet]
    reservoir.T  = 150.0                #Temperature of the reservoir [F]
    reservoir.phi   = np.reshape(poro, (numerical.N,1))         #porosity of the reservior vector [unitless]
    reservoir.permx = np.reshape(perm_x, (numerical.N,1))       #fluid permeability in x direction vector [mDarcy]
    #reservoir.permy = np.reshape(perm_y, (numerical.N,1))      #fluid permeability in x direction vector [mDarcy]
    reservoir.permx[reservoir.permx <= 1E-6] = 1E-16
    reservoir.permy  = 0.8 * reservoir.permx                    #fluid permeability in y direction vector [mDarcy]
    reservoir.permz  = 1.0 * reservoir.permx                    #fluid permeability in z direction vector [mDarcy]
    reservoir.Dref   = 9290.5                #reference depth [feet]
    reservoir.poro_ref   = 0.214            #reference porosity [-]
    reservoir.perm_ref   = 279.6            #reference perm [mD]
    reservoir.alpha  = 0.0*np.pi/6.0        #dip angle [in radians]
    reservoir.cfr    = 5.90E-6                 #formation of rock compressibility
    reservoir.Pref   = 7000                 #pressure at reference depth [psi]
    reservoir.Pwf    = 5300                 #reservoir BHP [psi]
    reservoir.Psc    = 14.7                 # Standard condition [psi]
    reservoir.phi[reservoir.phi==0] = 1E-16
    reservoir.flag   = np.reshape(flag,(numerical.N,1))           # Flag

    #Defining numerical parameters for discretized solution
    numerical.dx1 = np.reshape(dx[0,:], (numerical.Nx, 1))      #block thickness in x vector
    numerical.dy1 = np.reshape(dy[:,0], (numerical.Ny, 1))      #block thickness in y vector
    [numerical.dX,numerical.dY] = np.meshgrid(numerical.dx1,numerical.dy1)

    numerical.dx  = np.reshape(numerical.dX, (numerical.N,1))    #building the single dx column vector
    numerical.dy  = np.reshape(numerical.dY, (numerical.N,1))    #building the single dy column vector

    #position of the block centres x-direction
    numerical.xc = np.empty((numerical.Nx, 1))
    numerical.xc[0,0] = 0.5 * numerical.dx[0,0]
    for i in range(1,numerical.Nx):
        numerical.xc[i,0] = numerical.xc[i-1,0] + 0.5 * (numerical.dx1[i-1,0] + numerical.dx1[i,0])

    #position of the block centres y-direction
    numerical.yc = np.empty((numerical.Ny, 1))
    numerical.yc[0,0] = 0.5 * numerical.dy[0,0]
    for i in range(1,numerical.Ny):
        numerical.yc[i,0] = numerical.yc[i-1,0] + 0.5*(numerical.dy1[i-1,0] + numerical.dy1[i,0])
    
    [numerical.Xc,numerical.Yc] = np.meshgrid(numerical.xc,numerical.yc) 
    
    numerical.x1  = np.reshape(numerical.Xc, (numerical.N,1))    #building the single X column vector
    numerical.y1  = np.reshape(numerical.Yc, (numerical.N,1))    #building the single Y column vector

    #depth vector
    numerical.D = np.reshape(depth, (numerical.N,1))

    # Well parameters
    well.x = [numerical.Xc[2,8], numerical.Xc[3,3], numerical.Xc[1,2]]
    well.y = [numerical.Yc[2,8], numerical.Yc[3,3], numerical.Yc[1,2]]

    well.x = [[well.x[0]], [well.x[1]], [well.x[2]]]
    well.y = [[well.y[0]], [well.y[1]], [well.y[2]]]
    well.type = [[2], [1], [2]]                  # 1 for rate, 2 for BHP
    well.constraint = [[5300], [250 * 5.615], [5300]]      # rate = scf/day (+ for injector); BHP = psi (always +ve)
    well.rw = [[0.25], [0.25], [0.5]]              # well radius, ft
    well.skin = [[0], [0], [0]]                  # well skin friction factor, dimensionless
    well.direction = [['v'], ['v'], ['v']]         # direction of the well: vertical v or horizontal hx or hy
    print(well.x)
    print(well.y)
    # Boundary Conditions
    BC.type  = [['Dirichlet'],['Neumann'],['Dirichlet'],['Neumann']]    #type of BC: left, right, bottom, top
    BC.value = [[0],[0],[0],[0]]      #value of the boundary condition: psi, psi/ft or ft^3/day
    #BC.depth = [[reservoir.Dref],[reservoir.Dref-reservoir.L*np.sin(reservoir.alpha)],[reservoir.Dref],[reservoir.Dref-reservoir.L*np.sin(reservoir.alpha)]]
    
    IC.P     = reservoir.Pref * np.ones((numerical.N,1))      #Initial oil Pressure
    IC.Pw    = reservoir.Pref * np.ones((numerical.N,1))      #Initial water Pressure
    IC.Sw    = 0.50 * np.ones((numerical.N,1))                 #Initial saturation
    
    error = 1
    tol   = 1E-2

    while error > tol:
        IC.P_old = np.copy(IC.P)
        fluid_properties(reservoir, fluid, IC.P, IC.Pw)
        IC.P = reservoir.Pref + fluid.rho / 144 * (numerical.D - reservoir.Dref)
        IC.Pc, Pcprime = cap_press(petro, IC.Sw)
        IC.Pw = IC.P - IC.Pc
        fluid_properties(reservoir, fluid, IC.P, IC.Pw)
        error = abs(IC.P[0, 0] - IC.P_old[0, 0])

    IC.Sw[np.argwhere(IC.Sw <= petro.Swr)] = 0.18
    IC.Sw[np.argwhere(numerical.D == 0.0)] = 1.0
    IC.Pc[np.argwhere(numerical.D == 0.0)] = 0.0

inputfile(fluid,reservoir,petro,numerical,BC,IC,well)
from prodindex import prodindex
from updatewells import updatewells

#print(IC.Sw)
#print(IC.Pc)
'''
print(IC.Sw)
print(IC.Pw)
print(IC.P)
for i in range(0, numerical.Nx):
    for j in range(1,numerical.Nx):
        Tw, To = Thalf(i, j, 'x', fluid, reservoir, petro, numerical, IC.P, IC.Pw, IC.Pc, IC.Sw)
        print(Tw,To)

for i in range(0, numerical.Ny):
    for j in range(1,numerical.Ny):
        Tw, To = Thalf(i, j, 'y', fluid, reservoir, petro, numerical, IC.P, IC.Pw, IC.Pc, IC.Sw)
        print(Tw,To)

print(reservoir.phi, reservoir.permx, reservoir.permy)

reservoir.phi, reservoir.permx, reservoir.permy = rock_properties(reservoir, fluid, IC.P)

print(reservoir.phi, reservoir.permx, reservoir.permy)

print(flag)
Sw = IC.Sw
P = IC.P
Sw_hyst = IC.Sw
Tw, To, T, d11, d12, d21, d22, D, G, Q, Pc, Pw = myarrays(fluid,reservoir,petro,numerical,BC,P,Sw,Sw_hyst)
print(Pc)
i=1

Jwell_w,Jwell_o = prodindex(i,well,reservoir,fluid,petro,numerical,Sw)
'''


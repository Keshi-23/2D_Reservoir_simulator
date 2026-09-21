'''
Reservoir Simulation Assignment (2022)
Initializing a reservoir gridblock
Author: Promise O. Longe
Email: longepromise@ku.edu
Date modified: 04/18/2022
'''

#Importing required libraries

import numpy as np              #import numpy library
import matplotlib.pyplot as plt #library for plotting
plt.rcParams.update({'font.size': 22})
from matplotlib import cm

#Making grid class
class grid:
    def __init__(self):
        self.xmin = []
        self.xmax = []
        self.Nx = []
        self.N = []

#Importing data file
depth = np.loadtxt("depth.dat")
depth[depth==0.0] = np.nan                             #removing the zero depth for nice plotting
dx = np.loadtxt('block_x.dat')
dy = np.loadtxt('block_y.dat')
Nx = len(dx[0, :])  # number of grid blocks in x-direction
Ny = len(dy[:, 0])  # number of grid blocks in y-direction
N = Nx * Ny


#Input parameters
length    = np.sum(dx[0,:])       # length of the reservoir [feet]
breadth   = np.sum(dy[:,0])       # breadth of the reservoir [feet]
Nx        = Nx          # cells in x-direction
Ny        = Ny          # cells in y-direction
thic      = 22          # thickness of the reservoir [feet]
D_woc     = 9290     # depth of water oil contact line [feet]
P_w_woc   = 7000        # pressure at water oil contact line [psi]
rho_w_sc  = 62.4        # density of water at standard conditions [lbm/ft^3] 
rho_o_sc  = 53.0        # density of oil at standard conditions [lbm/ft^3]
rho_g_sc  = 0.0458171   # density of gas at standard conditions [lbm/ft^3]
c_w       = 1e-6     # compressibility of water [1/psi]
c_o       = 5.0e-5        # compressibility of oil [1/psi]
B_w       = 1.0        # formation volume factor of water
B_o       = 1.04567     # formation volume factor of oil 
Rs        = 923     # solution gas oil ratio [ft^3/bbl]
Pe        = 0.7171         # capillary entry pressure at water oil contact line [psi]
s_wr      = 0.18         # residual water saturation
s_or      = 0.1         # residual oil saturation
lam       = 1.52           # model parameter for capillary pressure

#Calculate water and oleic phase densities
rho_w     = rho_w_sc / B_w
rho_o     = (rho_o_sc + rho_g_sc * Rs /5.61)/B_o    #5.61 ft^3 in a barrel (bbl)

#Building grid
dx1 = np.reshape(dx[1,:], (Nx, 1))      #block thickness in x vector
dy1 = np.reshape(dy[:,1], (Ny, 1))      #block thickness in y vector

# position of the block centres x-direction
xc = np.empty((Nx, 1))
xc[0, 0] = 0.5 * dx[0, 0]
for i in range(1, Nx):
    xc[i, 0] = xc[i - 1, 0] + 0.5 * (dx1[i - 1, 0] + dx1[i, 0])

# position of the block centres y-direction
yc = np.empty((Ny, 1))
yc[0, 0] = 0.5 * dy[0, 0]
for i in range(1, Ny):
    yc[i, 0] = yc[i - 1, 0] + 0.5 * (dy1[i - 1, 0] + dy1[i, 0])

[Xc, Yc] = np.meshgrid(xc, yc)                  # centered- gridsblock

# Variation of pressure and saturation with depth (ignoring compressibility)
P_w = P_w_woc + rho_w / 144.0 * (depth - D_woc) #144 in^2 in 1 ft^2 
P_o = (P_w_woc + Pe) + rho_o / 144.0 * (depth - D_woc)  
s_w = s_wr + (1-s_wr)*((P_o-P_w)/Pe) **(-lam)  #From Corey-Brooks draining curve
s_w[depth>=D_woc] = 1.0
s_o = 1-s_w

#converting to a column vector
depth_col= np.reshape(np.transpose(depth), (N,-1))  #building the single depth vector
P_o_col  = np.reshape(np.transpose(P_o), (N,-1))    #building the single P_o vector
P_w_col  = np.reshape(np.transpose(P_w), (N,-1))    #building the single P_w vector
s_w_col  = np.reshape(np.transpose(s_w), (N,-1))    #building the single s_w vector
s_o_col  = np.reshape(np.transpose(s_o), (N,-1))    #building the single s_o vector

#Plotting starts here

#(a) Plotting the pressure variation with depth
fig = plt.figure(figsize=(15,7.5) , dpi=100)
plot = plt.plot(P_w_col,depth_col,'b.',label=r'$P_{water}$ [psi]')
plot = plt.plot(P_o_col,depth_col,'r.',label=r'$P_{oil}$ [psi]')
manager = plt.get_current_fig_manager()
manager.window.showMaximized()
plt.xlabel(r'$Pressure$ [psi]')
plt.ylabel(r'$Depth, D$ [m]')
plt.gca().invert_yaxis()
legend = plt.legend(loc='best', shadow=False, fontsize='x-large')
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig(f'Pressures.png',bbox_inches='tight', dpi = 600)

#(b) Plotting the water saturation variation with depth
fig = plt.figure(figsize=(15,7.5) , dpi=100)
plot = plt.plot(s_w_col,depth_col,'b.')
manager = plt.get_current_fig_manager()
manager.window.showMaximized()
plt.xlabel(r'$s_{w}$')
plt.ylabel(r'$Depth, D$ [m]')
plt.gca().invert_yaxis()
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig(f'Saturation.png',bbox_inches='tight', dpi = 600)

#(c) Plotting the 2D contour of depth
fig = plt.figure(figsize=(15,7.5) , dpi=100)
ax1= plt.contourf(Xc, Yc,depth,cmap=cm.coolwarm, antialiased=True)
plt.title('Depth of the reservoir')
plt.axis('scaled')
plt.xlabel(r'$x [feet]$')
plt.ylabel(r'$y [feet] $')
clb = plt.colorbar(ax1)
clb.set_label(r'$Depth,D [feet]$', labelpad=-40, y=1.1, rotation=0)
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig(f'reservoir_depth.png',bbox_inches='tight', dpi = 600)

#(d) Plotting the 2D contour of water saturation
fig = plt.figure(figsize=(15,7.5) , dpi=100)
ax1= plt.contourf(Xc, Yc,s_w,cmap=cm.coolwarm, antialiased=True)
plt.title('Water saturation in the reservoir')
plt.axis('scaled')
plt.xlabel(r'$ x [feet]$')
plt.ylabel(r'$y [feet] $')
clb = plt.colorbar(ax1)
clb.set_label(r'$s_w$', labelpad=-40, y=1.1, rotation=0)
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig(f'reservoir_water_saturation.png',bbox_inches='tight', dpi = 600)

#(e) Plotting the 2D contour of oil pressure
fig = plt.figure(figsize=(15,7.5) , dpi=100)
ax1= plt.contourf(Xc, Yc,P_o,cmap=cm.coolwarm, antialiased=True)
plt.title('Oil pressure in the reservoir')
plt.axis('scaled')
plt.xlabel(r'$x [feet]$')
plt.ylabel(r'$y [feet] $')
clb = plt.colorbar(ax1)
clb.set_label(r'$P_{oil} [psi]$', labelpad=-40, y=1.1, rotation=0)
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig(f'reservoir_oil_pressure.png',bbox_inches='tight', dpi = 600)
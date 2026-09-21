"""
reservoir simulation project 1(2022)
2D Multiphase reservoir simulation: Making arrays
Author: Promise O. Longe
Email: longe.promise@ku.edu
Date modified: 04/18/2022
"""
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import inv
from scipy.sparse.linalg import spsolve
from Thalf import Thalf   #for calculating transmissibility
from cap_press import cap_press
from rock_properties import rock_properties
from spdiaginv import spdiaginv
from numpy import zeros

class numerical:
    def __init__(self):
        self.Bw  = []
class reservoir:
    def __init__(self):
        self.dt = []
class fluid:
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
        self.xmin = []
class well:
    def __init__(self):
        self.xmin = []

cap = np.loadtxt("cap.txt")

#fluid, reservoir and simulation parameters
def myarrays(fluid,reservoir,petro,numerical,IC,BC,P,Sw):

    #Setting up matrix T, B, and Q
    T   = lil_matrix((numerical.N, numerical.N))
    Tw  = lil_matrix((numerical.N, numerical.N))
    To  = lil_matrix((numerical.N, numerical.N))
    B   = lil_matrix((numerical.N, numerical.N))
    d11 = lil_matrix((numerical.N, numerical.N))
    d12 = lil_matrix((numerical.N, numerical.N))
    d21 = lil_matrix((numerical.N, numerical.N))
    d22 = lil_matrix((numerical.N, numerical.N))
    D   = lil_matrix((numerical.N, numerical.N))
    G   = lil_matrix((numerical.N, 1))
    Q   = lil_matrix((numerical.N, 1))
    Qbw   = lil_matrix((numerical.N, 1))
    Qbo   = lil_matrix((numerical.N, 1))
    Pc  = zeros((numerical.N, 1))
    Pw  = zeros((numerical.N, 1))


    for l in range(0, numerical.N):
        Pc[l, 0], Pcprime = cap_press(petro, Sw[l, 0], cap=cap)
        Pw[l, 0] = P[l, 0] - Pc[l, 0]
        if (l + 1) % numerical.Nx != 1:  # not on left boundary
            Twhalf, Tohalf = Thalf(l, l - 1, 'x', fluid, reservoir, petro, numerical, P, Pw, Pc[l, 0], Sw)

            Tw[l, l - 1] = -Twhalf
            Tw[l, l] = Tw[l, l] - Tw[l, l - 1]

            To[l, l - 1] = -Tohalf
            To[l, l] = To[l, l] - To[l, l - 1]

        else:  # left boundary
            if 'Neumann' in BC.type[0]:
                None
                #Twblock, Toblock = Thalf(l, l, 'x', fluid, reservoir, petro, numerical, P, Pw, Pc[l, 0], Sw)

                #To[l, l] = To[l, l] + Toblock
                #Tw[l, l] = Tw[l, l] + Twblock

                #Qbw[l, 0] = Qbw[l, 0] - Twblock * ((BC.value[0][0]-fluid.rhw[l,0]/144) * (numerical.D[l,0]-BC.depth[0][0]) - Pc[l, 0])
                #Qbo[l, 0] = Qbo[l, 0] - Toblock * ((BC.value[0][0]-fluid.rho[l,0]/144) * (numerical.D[l,0]-BC.depth[0][0]))


            elif 'Dirichlet' in BC.type[0]:
                #None
                Twblock, Toblock = Thalf(l, l, 'x', fluid, reservoir, petro, numerical, P, Pw, Pc[l, 0], Sw)

                To[l, l] = To[l, l] + 2 * Toblock
                Tw[l, l] = Tw[l, l] + 2 * Twblock
                Qbw[l, 0] = Qbw[l, 0] + Twblock * (2 * BC.value[0][0])
                Qbo[l, 0] = Qbo[l, 0] + Toblock * (2 * BC.value[0][0])


        if (l + 1) % numerical.Nx != 0:  # not on right boundary

            Twhalf, Tohalf = Thalf(l, l + 1, 'x', fluid, reservoir, petro, numerical, P, Pw, Pc[l, 0], Sw)

            Tw[l, l + 1] = -Twhalf
            Tw[l, l] = Tw[l, l] - Tw[l, l + 1]

            To[l, l + 1] = -Tohalf
            To[l, l] = To[l, l] - To[l, l + 1]

        else:  # right  boundary
            if 'Neumann' in BC.type[1]:
                None
                #Twblock, Toblock = Thalf(l, l, 'x', fluid, reservoir, petro, numerical, P, Pw, Pc[l, 0], Sw)
                #To[l, l] = To[l, l] + Toblock
                #Tw[l, l] = Tw[l, l] + Twblock

                #Qbw[l, 0] = Qbw[l, 0] - Twblock * ((BC.value[1][0]-fluid.rhw[l,0]/144) * (numerical.D[l,0]-BC.depth[1][0]) - Pc[l, 0])
                #Qbo[l, 0] = Qbo[l, 0] - Toblock * ((BC.value[1][0]-fluid.rho[l,0]/144) * (numerical.D[l,0]-BC.depth[1][0]))

            elif 'Dirichlet' in BC.type[1]:
                #None
                Twblock, Toblock = Thalf(l, l, 'x', fluid, reservoir, petro, numerical, P, Pw, Pc[l, 0], Sw)

                To[l, l] = To[l, l] + 2 * Toblock
                Tw[l, l] = Tw[l, l] + 2 * Twblock
                Qbw[l, 0] = Qbw[l, 0] + Twblock * (2 * BC.value[1][0])
                Qbo[l, 0] = Qbo[l, 0] + Toblock * (2 * BC.value[1][0])

        if int(l / numerical.Nx) > 0:  # not bottom boundary
            Twhalf, Tohalf = Thalf(l, l - numerical.Nx, 'y', fluid, reservoir, petro, numerical, P, Pw, Pc[l, 0], Sw)

            Tw[l, l - numerical.Nx] = -Twhalf
            Tw[l, l] = Tw[l, l] - Tw[l, l - numerical.Nx]

            To[l, l - numerical.Nx] = -Tohalf
            To[l, l] = To[l, l] - To[l, l - numerical.Nx]

        else:  # bottom boundary
            if 'Neumann' in BC.type[2]:
                None
                #Twblock, Toblock = Thalf(l, l, 'y', fluid, reservoir, petro, numerical, P, Pw, Pc[l, 0], Sw)
                #To[l, l] = To[l, l] + Toblock
                #Tw[l, l] = Tw[l, l] + Twblock

                ##Qbw[l, 0] = Qbw[l, 0] - Twblock * ((BC.value[2][0]-fluid.rhw[l,0]/144) * (numerical.D[l,0]-BC.depth[2][0]) - Pc[l, 0])
                #Qbo[l, 0] = Qbo[l, 0] - Toblock * ((BC.value[2][0]-fluid.rho[l,0]/144) * (numerical.D[l,0]-BC.depth[2][0]))

            elif 'Dirichlet' in BC.type[2]:
                #None
                Twblock, Toblock = Thalf(l, l, 'y', fluid, reservoir, petro, numerical, P, Pw, Pc[l, 0], Sw)

                To[l, l] = To[l, l] + 2 * Toblock
                Tw[l, l] = Tw[l, l] + 2 * Twblock
                Qbw[l, 0] = Qbw[l, 0] + Twblock * (2 * BC.value[2][0])
                Qbo[l, 0] = Qbo[l, 0] + Toblock * (2 * BC.value[2][0])


        if int(l / numerical.Nx) < numerical.Ny - 1:  # not top boundary
            Twhalf, Tohalf = Thalf(l, l + numerical.Nx, 'y', fluid, reservoir, petro, numerical, P, Pw, Pc[l, 0], Sw)

            Tw[l, l + numerical.Nx] = -Twhalf
            Tw[l, l] = Tw[l, l] - Tw[l, l + numerical.Nx]

            To[l, l + numerical.Nx] = -Tohalf
            To[l, l] = To[l, l] - To[l, l + numerical.Nx]

        else:  # top boundary
            if 'Neumann' in BC.type[3]:
                None
                #Twblock, Toblock = Thalf(l, l, 'y', fluid, reservoir, petro, numerical, P, Pw, Pc[l, 0], Sw)
                #To[l, l] = To[l, l] + Toblock
                #Tw[l, l] = Tw[l, l] + Twblock

                #Qbw[l, 0] = Qbw[l, 0] - Twblock * ((BC.value[3][0]-fluid.rhw[l,0]/144) * (numerical.D[l,0]-BC.depth[3][0]) - Pc[l, 0])
                #Qbo[l, 0] = Qbo[l, 0] - Toblock * ((BC.value[3][0]-fluid.rho[l,0]/144) * (numerical.D[l,0]-BC.depth[3][0]))

            elif 'Dirichlet' in BC.type[3]:
                #None
                Twblock, Toblock = Thalf(l, l, 'y', fluid, reservoir, petro, numerical, P, Pw, Pc[l, 0], Sw)

                To[l, l] = To[l, l] + 2 * Toblock
                Tw[l, l] = Tw[l, l] + 2 * Twblock
                Qbw[l, 0] = Qbw[l, 0] + Twblock * (2 * BC.value[3][0])
                Qbo[l, 0] = Qbo[l, 0] + Toblock * (2 * BC.value[3][0])

        # B[l,l] = numerical.dx[l,0] * numerical.dy[l,0] * reservoir.h * reservoir.phi[l,0] * fluid.ct / fluid.Bw[l,0] #accumulation
        Vp = numerical.dx[l, 0] * numerical.dy[l, 0] * reservoir.h[l,0] * reservoir.phi[l, 0] / 5.615
        d11[l, l] = Vp * Sw[l, 0] * (fluid.cw + reservoir.cfr) / (fluid.Bw[l, 0] * numerical.dt)
        d12[l, l] = Vp / (fluid.Bw[l, 0] * numerical.dt) * (1.0 - Sw[l, 0] * reservoir.phi[l, 0] * fluid.cw * 0)
        d21[l, l] = Vp * (1 - Sw[l, 0]) * (fluid.co + reservoir.cfr) / (fluid.Bo[l, 0] * numerical.dt)
        d22[l, l] = -Vp / (fluid.Bo[l, 0] * numerical.dt)
        D[l, l] = -(d22[l, l] * d11[l, l] / d12[l, l]) + d21[l, l]

    d22 = d22.tocsr()
    d12 = d12.tocsr()
    d21 = d21.tocsr()
    d11 = d11.tocsr()

    Tw = (Tw).tocsr()  # multiplying with the conversion factor
    To = (To).tocsr()  # multiplying with the conversion factor
    T = (-d22 @ (spdiaginv(d12)) @ Tw) + To  # Weighing using the formula given in the sheet
    Q = (-d22 @ (spdiaginv(d12)) @ Qbw) + Qbo
    G = -d22 @ (spdiaginv(d12) @ (Tw @ (P - Pw))) - d22 @ (spdiaginv(d12) @ (Tw @ numerical.D)) * fluid.rhw / 144.0 + (fluid.rho[0,0] / 144.0 * To) @ numerical.D

    return Tw, To, T, d11, d12, d21, d22, D, Q, G, Pc, Pw;

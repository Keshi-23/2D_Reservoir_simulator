"""
reservoir simulation project 1(2022)
2D Multiphase reservoir simulation: Interblock transmissibility
Author: Promise O. Longe
Email: longe.promise@ku.edu
Date modified: 04/18/2022
"""

import numpy as np
import matplotlib.pyplot as plt
from rel_perm import rel_perm
from fluid_properties import fluid_properties
from rock_properties import rock_properties
from scipy import interpolate as interp
from scipy.special import logsumexp

#fluid,reservoir and simulation parameters

from rel_perm import rel_perm


# fluid, reservoir and simulation parameters
# fluid, reservoir and simulation parameters
def Thalf(i, j, direction, fluid, reservoir, petro, numerical, P, Pw, Pc, Sw):

    if direction == 'x':
        kAd = 0.001127 * (2 * reservoir.permx[i, 0] * numerical.dy[i, 0] * reservoir.h[i,0] * reservoir.permx[j, 0] * numerical.dy[j, 0] * reservoir.h[j,0]) / \
              (reservoir.permx[i, 0] * numerical.dy[i, 0] * reservoir.h[i,0] * numerical.dx[j, 0] + \
               reservoir.permx[j, 0] * numerical.dy[j, 0] * reservoir.h[j,0] * numerical.dx[i, 0])
    elif direction == 'y':
        kAd = 0.001127 * (2 * reservoir.permy[i, 0] * numerical.dx[i, 0] * reservoir.h[i,0] * reservoir.permy[j, 0] * numerical.dx[j, 0] * reservoir.h[j,0]) / \
              (reservoir.permy[i, 0] * numerical.dx[i, 0] * reservoir.h[i,0] * numerical.dy[j, 0] + \
               reservoir.permy[j, 0] * numerical.dx[j, 0] * reservoir.h[j,0] * numerical.dy[i, 0])

    # Calculating the potential for upstream weighting

    fluid_properties(reservoir, fluid, P, Pw)
    POT_i = P[i, 0] - (fluid.rho[i, 0] / (144.0 * fluid.Bo[i, 0])) * numerical.D[i, 0] - Pc
    POT_j = P[j, 0] - (fluid.rho[j, 0] / (144.0 * fluid.Bo[j, 0])) * numerical.D[j, 0] - Pc
    POTw_i = Pw[i, 0] - (fluid.rhw[i, 0] / (144.0 * fluid.Bw[i, 0])) * numerical.D[i, 0] - Pc
    POTw_j = Pw[j, 0] - (fluid.rhw[j, 0] / (144.0 * fluid.Bw[j, 0])) * numerical.D[j, 0] - Pc

    if POT_i >= POT_j:
        krw, kro = rel_perm(petro, Sw[i, 0])
    else:
        krw, kro = rel_perm(petro, Sw[j, 0])

    # Averaging fluid properties
    omega_i = reservoir.phi[i,0] * numerical.dy[i, 0] * reservoir.h[i,0] * numerical.dx[i, 0]
    omega_j = reservoir.phi[j,0] * numerical.dy[j, 0] * reservoir.h[j,0] * numerical.dx[j, 0]

    fluidw_avg = (omega_j / (fluid.muw[j, 0] * fluid.Bw[j, 0]) + omega_i / (fluid.muw[i, 0] * fluid.Bw[i, 0])) / (omega_i + omega_j)
    fluido_avg = (omega_j / (fluid.muo[j, 0] * fluid.Bo[j, 0]) + omega_i / (fluid.muo[i, 0] * fluid.Bo[i, 0])) / (omega_i + omega_j)

    fluidhalfw = krw * fluidw_avg
    fluidhalfo = kro * fluido_avg

    Tw = (fluidhalfw * kAd)
    To = (fluidhalfo * kAd)

    return Tw, To;
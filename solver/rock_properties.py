"""
reservoir simulation project 1(2022)
2D Multiphase reservoir simulation - Rock properties
Author: Promise O. Longe
Email: longe.promise@ku.edu
Date modified: 04/18/2022
"""
import numpy as np
from scipy import interpolate as interp
import matplotlib.pyplot as plt

class reservoir:
    def __init__(self):
        self.dt = [] 
        
def rock_properties(reservoir, P):

    reservoir.phi    = reservoir.phi
    reservoir.permx  = reservoir.permx
    reservoir.permy  = 0.8 * reservoir.permx
    reservoir.permz  = 1.0 * reservoir.permx

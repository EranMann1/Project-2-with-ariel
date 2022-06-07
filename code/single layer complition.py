# -*- coding: utf-8 -*-
"""
Created on Sat Feb 19 19:35:52 2022

@author: user
"""

import  numpy as np
import  matplotlib.pyplot as plt
import  scipy 

class single_layer:
    def __init__(self, wave_length, Lambda, D, H, Z, r_eff, y_grid, z_grid, sumaation_order):
        self.D = np.array(D)
        self.H = np.array(H)
        self.Z = np.array(Z)
        self.r_eff = np.array(r_eff)
        self.y_grid = np.array(y_grid)
        self.z_grid = np.array(z_grid)
        self.wave_length = wave_length
        self.k = np.divide(2 * np.pi, self.wave_length)
        self.Lambda = Lambda
        
        
    def calculate_alphas_and_betas(self):
        pass        
        
    
    
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 12 09:20:01 2022

@author: user
"""

import numpy as np


class layer:
    def __    def __init__(self, **kwargs):
        done = Falce
        if normalized in kwargs.keys:
            if kwargs.normalized:
                self.wave_length = 1
                self.eta = 1 
                self.c = 1
                self.k = 2 * np.pi/wave_length
                self.f = np.divide(self.c, cself.wave_length)
                done = True
        if not done:
            if 'c' in kwargs.keys():
                self.c = c 
            elif 
            
    # def __init__(self, wave_length, eta, r_eff, Hs, Ds, Zs, Lambda = 0):
        # self.wave_length = wave_length
        # self.eta = eta 
        # self.Lambda = Lambda 
        # self.r_eff = r_eff
        # self.Hs = np.array(Hs) + 0j
        # self.Ds = np.array(Ds) + 0j
        # self.Zs = np.array(Zs) + 0j
        # self.k = 2*np.pi/wave_length
     
        
    def calculate_alphas_betas_ns(self, phase, zeros = True, order = 100):
        if zeros:
            nums = np.linspace(-order, order, 2*order + 1) + 0j
        else:
            nums = np.concatenate((np.arange(- order, 0), np.arange(0, order))) + 0j
        alphas =  (2*np.pi*nums + phase)/self.Lambda
        betas = np.conj(np.sqrt((self.k)**2 - alphas**2))
        return alphas, betas, nums 
      
        
    def calculate_phase(self, theta_in):
        return np.mod(self.k * self.Lambda * np.sin(theta_in), 2 * np.pi)
      
        
    def calculate_theta(self, phase, Degree = False):
        theta = np.arcsin(np.divide(phase, self.k * self.Lambda))
        if Degree:
            return np.divide(theta * 180, np.pi)
        else:
            return theta
    
        
    def calculate_Lambda(self, theta_in, theta_out, order = 1, change = False):
        Lambda = np.abs(np.divide(order, np.sin(theta_in)-np.sin(theta_out)))
        if change:
            self.Lambda = Lambda
        return Lambda
    
    
      
        

    
            
        
    
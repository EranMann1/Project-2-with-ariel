# -*- coding: utf-8 -*-
"""
Created on Sat Feb 12 09:20:01 2022

@author: user
"""

import numpy as np


class layer:
    """
    tamplate to methagradings object with some layers, wavelentgh, and vacume parameters
    """
    def __init__(self, **kwargs):
        """
    
        Parameters
        ----------
        c \ speed : speed of light. defult 3e8
        freq\wave_lentgh\k : control wavelentgh in some way defult freq 10G
        Lambda\theta_in + (optional) theta_out, Lambda_order: control the period of the layers. defult is wave_lentgth
        r_eff : the effective radius of cabel. defult Lambda/100
        eta: impedance of vacume. defult:377
        
        normalized = True : makes eta, wavelength, c eqauls to 1


        """
        done = False
        # wave and vacume properties are 1
        if 'normalized' in kwargs.keys():
            if kwargs['normalized']:
                self.wave_length = 1
                self.eta = 1 
                self.c = 1
                self.k = np.divide(2 * np.pi, self.wave_length)
                self.freq = np.divide(self.c, self.wave_length)
                done = True
                
        # wave and vaccume properties are as stated or defult
        if not done:
            if 'c' in kwargs.keys():
                self.c = kwargs['c'] 
            elif 'speed' in kwargs.keys():
                self.c = kwargs['speed']
            else:
                self.c = 3e8
                
            if 'wave_length' in kwargs.keys():
                self.wave_length = kwargs['wave_length']
                self.freq = np.dividw(self.c, self.wave_length)
                self.k = np.divide(2 * np.pi, self.wave_length)
            elif 'freq' in kwargs.keys():
                self.freq = kwargs['freq']
                self.wave_length = np.divide(self.c, self.freq)
                self.k = np.divide(2 * np.pi, self.wave_length)
            elif 'k' in kwargs.keys():
                self.freq = kwargs['k']
                self.wave_length = np.divide(2 * np.pi, self.k)
                self.freq = np.dividw(self.c, self.wave_length)
            else:
                self.freq = 1e10
                self.wave_length = np.divide(self.c, self.freq)
                self.k = np.divide(2 * np.pi, self.wave_length)
                
            
            if 'eta' in kwargs.keys():
                self.eta = kwargs['eta']
            else:
                self.eta = 377
            
        # methagradings parameters    
        if 'Ds' in kwargs.keys():
            self.Ds = np.array(kwargs['Ds'])
        else:
            self.Ds = 0
        if 'Hs' in kwargs.keys():
            self.Hs = np.array(kwargs['Hs'])
        else:
            self.Hs = 0
        if 'Zs' in kwargs.keys():
            self.Zs = np.array(kwargs['Zs']) + 0j
        else:
            self.Zs = 0j
             
        if 'Lambda' in kwargs.keys():
            self.Lambda = kwargs['Lambda']
        elif 'theta_in' in kwargs.keys():
            if 'theta_out' in kwargs.keys():
                if 'Lambda_order' in kwargs.keys():
                    Lambda = np.abs(np.divide(kwargs['Lambda_order'],\
                                              np.sin(kwargs['theta_in']) - \
                                              np.sin(kwargs['theta_out'])))
                    self.Lambda = np.divide(Lambda, self.k)
                else:
                    Lambda = np.abs(np.divide(1,np.sin(kwargs['theta_in']) - \
                                              np.sin(kwargs['theta_out'])))
                    self.Lambda = np.divide(Lambda, self.k)
            elif 'Lambda_order' in kwargs.keys():
                Lambda = np.abs(np.divide(kwargs['Lambda_order'],\
                                          np.sin(kwargs['theta_in'])))
                self.Lambda = np.divide(Lambda, self.k)
            else:
                Lambda = np.abs(np.divide(1,np.sin(kwargs['theta_in'])))
                self.Lambda = np.divide(Lambda, self.k)
        else:
            self.Lambda = self.wave_length
        
        if 'r_eff' in kwargs.keys():
            self.r_eff = kwargs['r_eff']
        else:
            self.r_eff = np.divide(self.Lambda, 100)
            
     
        
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
    
    
      
        

    
            
        
    
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 12 09:20:01 2022

@author: user
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
plt.rc('font', size=10)          # controls default text sizes
plt.rc('axes', titlesize=15)     # fontsize of the axes title
plt.rc('axes', labelsize=15)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
plt.rc('legend', fontsize=10)    # legend fontsize
plt.rc('figure', titlesize=14)  # fontsize of the figure title


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
        Ds,Hs,Zs: vectors. defult:0
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
            self.Ds = np.array([0])
        if 'Hs' in kwargs.keys():
            self.Hs = np.array(kwargs['Hs'])
        else:
            self.Hs = np.array([0])
        if 'Zs' in kwargs.keys():
            self.Zs = np.array(kwargs['Zs']) + 0j
        else:
            self.Zs = np.array([0j])
             
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
    
    
class source:
    def __init__(self, y_vec, z_vec, I_vec):
        self.y_vec = y_vec
        self.z_vec = z_vec
        self.I_vec = I_vec
    
    
    def calc_field(self, y, z, grid = True):
        #either y_vec,z_vec or y_grid, z_frid
        pass
    
    def calc_phasor_field(self):
        pass
    
    
class space:
    def __init__(self, layer, y_vec, z_vec, source):
        self.layer = layer
        self.y_vec = np.array(y_vec)
        self.z_vec = np.array(z_vec)
        self.source = source
        
    
    
    def calc_phasor_field(self, phase, summation_order,  currents = 0, save = False):
        alphas, betas, nums = self.layer.calculate_alphas_betas_ns(phase, True, summation_order)
        # all vectors here are in directions(y,z,num_of_layers,summation order)
        y_vec = self.y_vec.reshape(len(self.y_vec), 1, 1, 1)
        z_vec = self.z_vec.reshape(1, len(self.z_vec), 1, 1)
        Alphas = alphas.reshape(1, 1, 1, len(nums))
        Betas = betas.reshape(1, 1, 1, len(nums))
        Ds = self.layer.Ds.reshape(1, 1, len(self.layer.Ds), 1)
        Hs = self.layer.Hs.reshape(1, 1, len(self.layer.Ds), 1)
        if currents == 0:
            Currents = np.ones(np.shape(Ds))
        else:
            Currents = np.array(currents).reshape(1, 1, len(self.layer.Ds), 1)
            
        sum_element = Currents * np.divide(np.exp(-1j * Alphas * (y_vec - Ds)) \
                                           * np.exp(-1j * Betas * \
                                                    np.abs(z_vec - Hs)), Betas)
        if save:
            self.field = np.sum(sum_element, axis = (-1, -2))
            return self.field
        else:
            return np.sum(sum_element, axis = (-1, -2))
    
    
    def calc_total_field(self, phasor_vec, summation_order):
        pass
    
    
    def plot(self, saved = True, total = True, plot_source = True, type = 'real', **kwargs,):
        # data aquasition
        if saved:
            y_vec, z_vec = np.meshgrid(self.y_vec, self.z_vec)
            field = self.field
        else:
            field = kwargs['field']
            if 'y_vec' in kwargs.keys():
                y_vec, z_vec = np.meshgrid(kwargs['y_vec'], kwargs['z_vec'])
            else:
                y_vec, z_vec = np.meshgrid(self.y_vec, self.z_vec)
        if total:
            field += self.source.calc_field(y_vec, z_vec)
        
        # taking layers place
        Ds = self.layer.Ds
        Hs = self.layer.Hs
        z_max = np.max(z_vec)
        y_max = np.max(y_vec)
        z_min = np.min(z_vec)
        y_min = np.min(y_vec)
        Lambda = self.layer.Lambda
        layer_list = []
        for i in range(len(Ds)):
            if Hs[i] < z_max and  Hs[i] > z_min:
                num = np.ceil(np.divide(y_min - Ds[i], Lambda))
                while Ds[i] + num * Lambda < y_max:
                    layer_list.append(np.array([Ds[i] + num * Lambda, Hs[i]]))
                    num += 1
        if plot_source:
            source_y_vec = self.source.y_vec
            source_z_vec = self.source.z_vec
            source_list = []
            for i in range(len(source_y_vec)):
                if source_y_vec[i] < y_max  and  source_y_vec[i] > y_min and \
                    source_z_vec[i] < z_max and source_z_vec[i] > z_min:
                        source_list.append(np.array([source_y_vec[i], source_z_vec[i]]))
                    
        
        # plotting:
        plt.figure()
        if type == 'real':
            fig = plt.pcolormesh(y_vec, z_vec, np.transpose(np.real(field)))
        if type == 'abs':
            fig = plt.pcolormesh(y_vec, z_vec, np.transpose(np.real(field)))
        if type == 'abslog':
            fig = plt.pcolormesh(y_vec, z_vec, np.transpose(np.abs(field)), norm=colors.LogNorm())
        if type == 'reallog':
            fig = plt.pcolormesh(y_vec, z_vec, np.transpose(np.abs(np.real(field))), norm=colors.LogNorm())
        if type == 'powerlog':
            fig = plt.pcolormesh(y_vec, z_vec, np.transpose(np.abs(np.abs(field) ** 2)), norm=colors.LogNorm())
        if type == 'power':
            fig = plt.pcolormesh(y_vec, z_vec, np.transpose(np.abs(np.abs(field) ** 2)))
        plt.colorbar(fig)
        
        for item in layer_list:
            plt.plot(item[0], item[1], 'b.')
        if plot_source:
            for item in source_list:
                plt.plot(item[0], item[1], 'k.')
            
        if 'title' in kwargs.keys():
            plt.title(kwargs['title'])
        if 'xlabel' in kwargs.keys():
            plt.xlabel(kwargs['xlabel'])
        if 'ylabel' in kwargs.keys():
            plt.ylabel(kwargs['ylabel'])
        
        
    
    
    
    
        

    
            
        
    
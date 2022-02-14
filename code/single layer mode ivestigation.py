# -*- coding: utf-8 -*-
"""
Created on Sat Feb 12 09:20:01 2022

@author: user
"""

import numpy as np
import matplotlib.pyplot as plt
from methagrading import layer
plt.rc('font', size=10)          # controls default text sizes
plt.rc('axes', titlesize=15)     # fontsize of the axes title
plt.rc('axes', labelsize=15)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
plt.rc('legend', fontsize=10)    # legend fontsize
plt.rc('figure', titlesize=14)  # fontsize of the figure title

class single_layer(layer):
    def __init__(self, **kwargs):
        layer.__init__(self, **kwargs)
        self.Lambda = np.array(kwargs['Lambda_vector'])
        self.order = kwargs['summation_order']
        self.phase_number = kwargs['phase_number']
            
    
    def calculate_phase_matrix(self):
        Lambdas = self.Lambda
        k = self.k
        phase_matrix = np.zeros([len(Lambdas), self.phase_number])
        theta_matrix = phase_matrix
        for Lambda_index in range(len(Lambdas)):
            Lambda = Lambdas[Lambda_index]
            start = k * Lambda
            stop = 2 * np.pi -k * Lambdas[Lambda_index]
            phase_matrix[Lambda_index, :] = np.linspace(start, stop, self.phase_number, endpoint = False)
            #theta_matrix[Lambda_index, :] = np.arcsin(np.divide(phase_matrix[Lambda_index, :], k * Lambda))
        self.phase_matrix = phase_matrix
        #self.theta_matrix = np.divide(theta_matrix * 180, np.pi)
        
        
    def calculate_imaginary_valur(self):
        Lambda_vector = self.Lambda
        phase_matrix = self.phase_matrix
        imaginary_value = np.zeros(np.shape(phase_matrix))
        k = self.k
        num_vector = np.concatenate([np.linspace(-self.order, 0, endpoint = False),\
                                     np.linspace(1, self.order)])
        for Lambda_index in range(len(Lambda_vector)):
            Lambda = Lambda_vector[Lambda_index]
            phase, nums = np.meshgrid(phase_matrix[Lambda_index, :], num_vector)
            sum_element_0 = -np.divide(1, np.sqrt(phase[0, :]**2 - (Lambda * k)**2))
            sum_elements = np.divide(1, 2 * np.pi * np.abs(nums)) - \
                np.divide(1, np.sqrt((phase + 2 * np.pi * nums)**2 - \
                                     (Lambda * k)**2))
            imaginary_value[Lambda_index, :] = np.divide(k * self.eta, 2) * \
                (sum_element_0 + np.sum(sum_elements,axis = 0))
        print(imaginary_value)
        print('and now')
        print(phase_matrix)
        self.imaginary_value = imaginary_value
        
        
    def plot(self):
        phase_matrix = self.phase_matrix
        imaginary_value = self.imaginary_value
        Lambda_vector = self.Lambda
        plt.figure()
        for Lambda_index in range(len(Lambda_vector)):
            plt.plot(np.divide(phase_matrix[Lambda_index, :], np.pi), \
                     - imaginary_value[Lambda_index, :],\
                     label = r'$\Lambda=' + str(np.around(Lambda_vector[Lambda_index],2)) + '$')
            plt.yscale('log')
        plt.legend(loc = 'upper right')
        plt.ylabel(r'$\frac{k}{2\pi}\ln{\frac{2\pi r_{eff}}{\Lambda}}-\frac{\Im(Z)}{\eta}$ ')
        plt.xlabel(r"$\frac{\Delta\phi}{2\pi}$")
        
        plt.grid(True, which = 'both')
        plt.show()
        
                     
             
def main():
    Lambda_vector = np.linspace(0.01, 0.49, 7)
    phase_number = 10000
    summation_order = 100000
    calculator = single_layer(normalized = True, Lambda_vector = Lambda_vector\
                              , phase_number = phase_number, \
                              summation_order = summation_order)
    calculator.calculate_phase_matrix()
    calculator.calculate_imaginary_valur()
    calculator.plot()
if __name__ == '__main__':
    phase = main()
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 27 20:20:12 2022

@author: user
"""
import matplotlib.pyplot as plt
import  numpy as np
from matplotlib import colors
plt.rc('font', size=10)          # controls default text sizes
plt.rc('axes', titlesize=15)     # fontsize of the axes title
plt.rc('axes', labelsize=15)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
plt.rc('legend', fontsize=10)    # legend fontsize
plt.rc('figure', titlesize=14)  # fontsize of the figure title
 

def calculate_alphas_betas_ns(Lambda_vec, phase_vec, order = 100):
        nums = np.linspace(-order, order, 2 * order + 1).reshape(1, 1, 1, -1) + 0j
        alphas =  np.divide(2 * np.pi * nums + phase_vec, Lambda_vec)
        betas = np.conj(np.sqrt((2 * np.pi) ** 2 - alphas ** 2))
        # first (0) direction phase, 4th (3) direction num
        return alphas, betas 
    

def plot(B_vec, Lambda_ky_vec, y_vector, z_vector): 
    Lambda_vec = Lambda_ky_vec[:, 0]
    num = np.max(Lambda_vec.shape)
    Lambda_vec = Lambda_vec.reshape(-1)
    ky_vec = Lambda_ky_vec[:, 1].reshape(-1)
    for index in range(num):
        plt.figure()
        fig = plt.pcolormesh(y_vector.reshape(-1), z_vector.reshape(-1),\
                             np.transpose(np.real(B_vec[index, :, :])),\
                             norm = colors.LogNorm())
        y_max = np.max(y_vector)
        z_max = np.max(z_vector)
        y_min = np.min(y_vector)
        z_min = np.min(z_vector)
        if z_min < 0 and z_max > 0:
            num = np.ceil(np.divide(y_min, Lambda_vec[index]))
            layer_list = []
            while num * Lambda_vec[index] < y_max:
                layer_list.append(np.array([num * Lambda_vec[index], 0]))
                num += 1
        plt.plot(layer_list[0], layer_list[1], 'b.')


def calc_B_vec(z_vector, y_vector, order, Lambda_ky_vec):
    y_vector = y_vector.reshape(1, -1 ,1, 1)
    z_vector = z_vector.reshape(1, 1, -1, 1)
    Lambda_vec = Lambda_ky_vec[:, 0].reshape(-1, 1, 1, 1)
    ky_vec = Lambda_ky_vec[:, 1].reshape(-1, 1, 1, 1)
    phase_vec = Lambda_vec * ky_vec 
    alphas, betas = calculate_alphas_betas_ns(Lambda_vec, phase_vec, order)
    B_vec = np.sum(np.divide(np.exp(-1j * (alphas * y_vector + betas * z_vector)), betas) , axis = -1)
    return B_vec # directions are (ky\Lambda pairs, y, z)

def main():
    z_vector = np.linspace(-0.1, 2.1, 100).reshape(1, -1, 1)
    y_vector = np.linspace(-0.1, 2.1, 100).reshape(-1, 1, 1)
    order = 1000
    Lambda_ky_vec = np.array([[0.17, 1.3], [0.17, 2.5], [0.17,3.5], [2.1 ,1], [2.1 ,3]])
    B_vec = calc_B_vec(z_vector, y_vector, order, Lambda_ky_vec)
    plot(B_vec, Lambda_ky_vec, y_vector, z_vector)
    
    
    
if __name__ == '__main__':
    main()
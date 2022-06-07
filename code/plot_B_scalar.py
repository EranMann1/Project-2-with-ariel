# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 09:08:46 2022

@author: user
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import pdb

def calc_alphas_betas(k, ky, Lambda, summation_order=100):
    phase = ky * Lambda
    nums = np.arange(-summation_order, summation_order)
    alphas = (2 * np.pi * nums + phase) / Lambda + 0j
    betas = np.conj(np.sqrt(k ** 2 - alphas ** 2 + 0j))
    #pdb.set_trace()
    return alphas, betas 


def calc_B_scalar(y_vec, z_vec, k, ky, Lambda, summation_order=100):
    alphas, betas = calc_alphas_betas(k, ky, Lambda, summation_order + 1)
    y_vec = y_vec.reshape(-1, 1, 1)
    z_vec = z_vec.reshape(1, -1, 1)
    alphas = alphas.reshape(1, 1, -1)
    betas = betas.reshape(1, 1, -1)
    pdb.set_trace()
    B = k * np.sum(np.divide(np.exp(-1j * alphas * y_vec) * np.exp(-1j * betas * np.abs(z_vec)), betas), axis=-1)
    return B


def plot_B_scalar(B, y_vec, z_vec, Lambda):
    plt.figure()
    y_vec, z_vec = np.meshgrid(y_vec, z_vec)
    fig = plt.pcolormesh(y_vec, z_vec, np.transpose(np.abs(B)))
    plt.colorbar(fig)
    ## plotting the source
    if np.max(z_vec) > 0 and np.min(z_vec) < 0:
        start_index = np.ceil(np.min(y_vec) / Lambda)
        index = start_index
        while Lambda * index < np.max(y_vec):
            plt.plot(Lambda * index, 0, 'b.')
            index += 1 

def main():
    summation_order = 10
    k = 2 * np.pi
    Lambdas = [0.17, 0.5, 0.6, 1.1, 1.5, 1.7]
    Kys = [0.3, 0.7, 1.3, 1.7]
    y_vec = np.linspace(0, 2, 101)
    z_vec = np.linspace(-0.1, 10, 100)
    for Lambda in Lambdas:
        for ky in Kys:
            B = calc_B_scalar(y_vec, z_vec, k, ky ,Lambda, summation_order)
            plot_B_scalar(B, y_vec, z_vec, Lambda)
            plt.xlabel(r'$y \left[\lambda\right]$')
            plt.ylabel(r'$z \left[\lambda\right]$')
            plt.title(r'$\left|B\left(\Lambda={},k_y={}\right)\right|$'.format(Lambda,ky))


if __name__ == '__main__':
    main()
    
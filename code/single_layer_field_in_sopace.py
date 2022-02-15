# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 12:58:05 2022

@author: user
"""

import matplotlib.pyplot as plt
import  numpy as np
from methagrading import layer, space


def main():
    single_layer = layer(normalized = True, Lambda = 0.17, Ds =[0], Hs = [0], Zs = [0])
    y_vec = np.linspace(- 0.1 * single_layer.Lambda, 10.1 * single_layer.Lambda, 1000)
    z_vec = np.linspace(-0.1, 2.1, 100)
    phase = 1.3 * single_layer.k * single_layer.Lambda
    my_space = space(single_layer, y_vec, z_vec, 0)
    summation_order = 1000
    my_space.calc_phasor_field(phase, summation_order, save = True, currents = [10])
    my_space.plot(total = False, plot_source = False, type = 'reallog', \
                  xlabel = r'$\frac{y}{\lambda}$', ylabel = r'$\frac{z}{\lambda}$'\
                  , title = r'Electric field $\left[\frac{I_0 \eta}{\lambda}\right]$'    )
    my_space.plot(total = False, plot_source = False, type = 'powerlog', \
                  xlabel = r'$\frac{y}{\lambda}$', ylabel = r'$\frac{z}{\lambda}$'\
                  , title = r'Energy per volume $\left[\frac{I_0^2 \eta^3}{2 \lambda^2}\right]$')
    plt.figure()
    plt.plot(my_space.z_vec, np.abs(my_space.field[1, :]) ** 2)
    plt.xlim(0, 2)
    plt.grid(which = 'both')
    plt.yscale('log')
    plt.xlabel(r'$\frac{z}{\lambda}$')
    plt.title(r'Energy per volume $\left[\frac{I_0^2 \eta^3}{2 \lambda^2}\right]$')

    plt.figure()
    plt.plot(my_space.y_vec, np.real(my_space.field[:, 1]) )
    plt.grid(which = 'both')
    plt.xlabel(r'$\frac{y}{\lambda}$')
    plt.title(r'Electric field $\left[\frac{I_0 \eta}{\lambda}\right]$')

if __name__ == '__main__':
    main()

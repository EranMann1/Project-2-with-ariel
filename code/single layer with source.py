# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 08:48:53 2022

@author: user
"""

from methagrading import single_layer, layer
import numpy as np
import  matplotlib.pyplot as plt


def main():
    # calculating approx and non approx STR and plotting
    # triel 1
    
    Lambda = np.array([0.17]) 
    k = 2 * np.pi 
    eta = 1
    summation_order = 10000000
    
    norm_k_y_vec = np.linspace(1, 2 * np.pi , 100000)
    norm_k_z_vec = np.sqrt(norm_k_y_vec ** 2 - 1 + 0j)
    approx_func = single_layer.Approx_func(norm_k_z_vec, Lambda)
    
    
    phase_vec = norm_k_y_vec * Lambda * k
    non_approx_func = single_layer.non_approx_func(Lambda, phase_vec, k , summation_order)
    
    wanted_phase = 2.5 * k * Lambda
    norm_k_y = np.divide(1, k * Lambda) * wanted_phase
    norm_k_z = np.array([np.sqrt(norm_k_y ** 2 - 1 + 0j)])
    wanted_func = single_layer.Approx_func(norm_k_z,Lambda)
    approx_STR_25 = np.divide(1, eta * np.abs(wanted_func - approx_func))
    non_approx_STR_25 = np.divide(1, eta * np.abs(1j * wanted_func - non_approx_func)).reshape(-1)
    
    wanted_phase = 3.5 * k * Lambda
    norm_k_y = np.divide(1, k * Lambda) * wanted_phase
    norm_k_z = np.array([np.sqrt(norm_k_y ** 2 - 1 + 0j)])
    wanted_func = single_layer.Approx_func(norm_k_z,Lambda)
    approx_STR_35 = np.divide(1, eta * np.abs(wanted_func - approx_func))
    non_approx_STR_35 = np.divide(1, eta * np.abs(1j * wanted_func - non_approx_func)).reshape(-1)
    
    wanted_phase = 1.5 * k * Lambda
    norm_k_y = np.divide(1, k * Lambda) * wanted_phase
    norm_k_z = np.array([np.sqrt(norm_k_y ** 2 - 1 + 0j)])
    wanted_func = single_layer.Approx_func(norm_k_z, Lambda)
    approx_STR_15 = np.divide(1, eta * np.abs(wanted_func - approx_func))
    non_approx_STR_15 = np.divide(1, eta * np.abs(1j * wanted_func - non_approx_func)).reshape(-1)
    
    plt.figure()
    plt.plot(norm_k_y_vec, approx_STR_15, 'r--', label = r'approximation $k_0=1.5$')
    plt.plot(norm_k_y_vec, non_approx_STR_15, 'r-', label = r'summation $k_0=1.5$')
    plt.plot(norm_k_y_vec, approx_STR_25, 'b--', label = r'approximation $k_0=2.5$')
    plt.plot(norm_k_y_vec, non_approx_STR_25, 'b-', label = r'summation $k_0=2.5$')
    plt.plot(norm_k_y_vec, approx_STR_35, 'g--', label = r'approximation $k_0=3.5$')
    plt.plot(norm_k_y_vec, non_approx_STR_35, 'g-', label = r'summation $k_0=3.5$')
    plt.ylim([1e-2, 1e6])
    plt.legend()
    plt.yscale('log')
    plt.xlabel(r'$\frac{k_y}{k}$')
    plt.ylabel('STR')
    plt.show()
    


    
    # wanted_func = single_layer.non_approx_func(Lambda, wanted_phase, k, summation_order)
    # phase_vector = np.linspace(0, 2 * np.pi , 1000)
    # norm_k_y_vec = np.divide(1, k * Lambda) * phase_vector
    # norm_k_z_vec = np.sqrt(norm_k_y_vec ** 2 - 1 + 0j)
    # Approx_func = single_layer.Approx_func(norm_k_z_vec, Lambda)
    # non_Approx_func = single_layer.non_approx_func(Lambda, phase_vector, k, summation_order)
    # approx_STR = np.divide(eta, np.abs(Approx_func - wanted_func))
    # STR = np.divide(eta, np.abs(non_Approx_func + wanted_func))
    # plt.figure()
    # plt.plot(norm_k_y_vec, STR.reshape(len(norm_k_z_vec)), 'b-', label = 'real')
    # plt.plot(norm_k_y_vec, approx_STR, 'r--', label = 'approximation')
    # plt.legend()
    # plt.yscale('log')
    # plt.xlabel(r'$\frac{k_y}{k}$')
    # plt.ylabel('STR')
    # plt.show()
    
    # plt.figure()
    # plt.plot(norm_k_y_vec, np.imag(non_Approx_func).reshape(len(norm_k_z_vec)), 'b-', label = 'real')
    # plt.plot(norm_k_y_vec, Approx_func, 'r--', label = 'approximation')
    # plt.legend()
    # plt.yscale('log')
    # plt.show()
    
    # #triel 2
    # phase = 1.3 * 2 * np.pi * 0.17
    # Zs = single_layer.calc_wanted_impedance(approximation = False, r_eff = 1/1000,\
    #                                         Lambda = np.array([0.17]), phase = phase)
    # print(Zs)
    # our_single_layer = layer(normalized = True, Lambda = 0.17, Zs = Zs, r_eff = 1/1000)
    # ky = np.linspace(0 * np.pi, 50 * np.pi, 10000)
    # summation_order = 1000
    # our_single_layer.calc_A_matricis(ky, summation_order)
    # print(our_single_layer.A_matricis)
    # our_single_layer.Plot_STR()
    
    
    
    
    
if __name__ == '__main__':
    main()
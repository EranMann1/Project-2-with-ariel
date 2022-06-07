# -*- coding: utf-8 -*-
"""
Created on Sat Feb 12 09:20:01 2022

@author: user
"""
import numpy as np
import matplotlib.pyplot as plt
from methagrading import layer, single_layer
from sklearn.linear_model import LinearRegression
plt.rc('font', size=10)          # controls default text sizes
plt.rc('axes', titlesize=15)     # fontsize of the axes title
plt.rc('axes', labelsize=15)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
plt.rc('legend', fontsize=10)    # legend fontsize
plt.rc('figure', titlesize=14)  # fontsize of the figure title




class single_layer_new(layer):
    def __init__(self, **kwargs):
        layer.__init__(self, **kwargs)
        self.Lambda = np.array(kwargs['Lambda_vector'])
        self.order = kwargs['summation_order']
        self.phase_number = kwargs['phase_number']
            
    
    def calculate_phase_matrix(self):
        Lambdas = self.Lambda
        k = self.k
        phase_matrix = np.zeros([len(Lambdas), self.phase_number])
        norm_k_y_matrix = np.zeros(np.shape(phase_matrix))
        for Lambda_index in range(len(Lambdas)):
            Lambda = Lambdas[Lambda_index]
            start = k * Lambda + np.divide(np.pi -k*Lambda, self.phase_number)
            stop = 2 * np.pi - k * Lambda
            phase_matrix[Lambda_index, :] = np.linspace(start, stop, self.phase_number, endpoint = False)
            norm_k_y_matrix[Lambda_index, :] = np.divide(phase_matrix[Lambda_index, :], Lambda * k)
        self.phase_matrix = phase_matrix
        self.norm_k_y_matrix = norm_k_y_matrix
        self.norm_k_z_matrix = np.sqrt(norm_k_y_matrix ** 2 - 1 ** 2)
        
        
    def calculate_imaginary_value(self):
        Lambda_vector = self.Lambda
        phase_matrix = self.phase_matrix
        imaginary_value = np.zeros(np.shape(phase_matrix))
        k = self.k
        num_vector = np.concatenate([np.linspace(-self.order, 0, endpoint = False),\
                                     np.linspace(1, self.order)])
        for Lambda_index in range(len(Lambda_vector)):
            Lambda = Lambda_vector[Lambda_index]
            phase = np.array(phase_matrix[Lambda_index, :]).reshape(-1, 1)
            nums = np.array(num_vector).reshape(1,-1)
            sum_element_0 = -np.divide(1, np.sqrt(phase ** 2 - (Lambda * k)**2))
            sum_elements = np.divide(1, 2 * np.pi * np.abs(nums)) - \
                np.divide(1, np.sqrt((phase + 2 * np.pi * nums)**2 - \
                                     (Lambda * k)**2))
            imaginary_value[Lambda_index, :] = np.divide(k * self.eta, 2) * \
                (sum_element_0[:,0] + np.sum(sum_elements,axis = -1))
        self.imaginary_value = imaginary_value
        
        
    def calculate_imaginary_value_with_alphas_betas(self):
        Lambda_vector = self.Lambda
        phase_matrix = self.phase_matrix
        imaginary_value = np.zeros(np.shape(phase_matrix)) + 0j
        k = self.k
        for Lambda_index in range(len(Lambda_vector)):
            Lambda = np.array(Lambda_vector[Lambda_index])
            phase = np.array(phase_matrix[Lambda_index, :])
            imaginary_value[Lambda_index, :] = single_layer.non_approx_func(Lambda, phase, k, self.order)
        self.imaginary_value = - np.divide(imaginary_value, 1j)
        
        
    def plot(self):
        phase_matrix = self.phase_matrix
        imaginary_value = self.imaginary_value
        Lambda_vector = self.Lambda
        plt.figure()
        for Lambda_index in range(len(Lambda_vector)):
            plt.plot(np.divide(phase_matrix[Lambda_index, :], np.pi), \
                     - np.real(imaginary_value[Lambda_index, :]),\
                     label = r'$\Lambda=' + str(np.around(Lambda_vector[Lambda_index],2)) + '\lambda$')
            plt.yscale('log')
        plt.legend(loc = 'upper right')
        plt.ylabel(r'$\frac{k}{2\pi}\ln{\frac{2\pi r_{eff}}{\Lambda}}-\frac{\Im(Z)}{\eta}$ ')
        plt.xlabel(r"$\frac{\Delta\phi}{2\pi}$")
        
        plt.grid(True, which = 'both')
        plt.show()
        
        plt.figure()
        for Lambda_index in range(len(Lambda_vector)):
            plt.plot(self.norm_k_y_matrix[Lambda_index, :], \
                     - np.real(imaginary_value[Lambda_index, :]),\
                     label = r'$\Lambda=' + str(np.around(Lambda_vector[Lambda_index],2)) + '\lambda$')
            plt.yscale('log')
            plt.xscale('log')
        plt.legend(loc = 'upper right')
        plt.ylabel(r'$\frac{k}{2\pi}\ln{\frac{2\pi r_{eff}}{\Lambda}}-\frac{\Im(Z)}{\eta}$ ')
        plt.xlabel(r"$\frac{k_y}{k}$")
        
        plt.grid(True, which = 'both')
        plt.show()
        
        plt.figure()
        for Lambda_index in range(len(Lambda_vector)):
            plt.plot(self.norm_k_z_matrix[Lambda_index, :], \
                     - np.real(imaginary_value[Lambda_index, :]),\
                     label = r'$\Lambda=' + str(np.around(Lambda_vector[Lambda_index],2)) + '\lambda$')
            plt.yscale('log')
            plt.xscale('log')
        plt.legend(loc = 'upper right')
        plt.ylabel(r'$\frac{k}{2\pi}\ln{\frac{2\pi r_{eff}}{\Lambda}}-\frac{\Im(Z)}{\eta}$ ')
        plt.xlabel(r"$\left| \Im{\left(\frac{k_z}{k}\right)}\right|$")
        
        plt.grid(True, which = 'both')
        plt.show()
     
        
    def linear_regration(self):
        Lambda_vec = self.Lambda
        A_vec = np.zeros(np.shape(Lambda_vec))
        B_vec = np.zeros(np.shape(Lambda_vec))
        score_vec = np.zeros(np.shape(Lambda_vec))
        X_matrix = np.log(self.norm_k_z_matrix)
        y_matrix = np.log(- self.imaginary_value)
        for Lambda_index in range(len(Lambda_vec)):
            X = X_matrix[Lambda_index, :].reshape(-1,1)
            y = y_matrix[Lambda_index, :]
            reg = LinearRegression().fit(X, y)
            A_vec[Lambda_index] = np.exp(reg.coef_)
            B_vec[Lambda_index] = reg.intercept_
            score_vec[Lambda_index] = reg.score(X, y)
        self.A_coefs = A_vec
        self.B_coefs = B_vec
        self.scores = score_vec


        
       
                     
             
def main():
    #graphs
    Lambda_vector = np.linspace(0.01, 0.49, 7)
    phase_number = 20000
    summation_order = 1000000
    calculator = single_layer_new(normalized = True, Lambda_vector = Lambda_vector\
                              , phase_number = phase_number, \
                              summation_order = summation_order)
    calculator.calculate_phase_matrix()
    calculator.calculate_imaginary_value()
    calculator.plot()
    
    calculator.calculate_imaginary_value_with_alphas_betas()
    calculator.plot()
    # regretion
    # Lambda_vector = np.linspace(0.01, 0.49, 1000)
    # phase_number = 100
    # summation_order = 10000
    # Regretor = single_layer_new(normalized = True, Lambda_vector = Lambda_vector\
    #                           , phase_number = phase_number, \
    #                           summation_order = summation_order)
    # Regretor.calculate_phase_matrix()
    # Regretor.calculate_imaginary_value()
    # Regretor.linear_regration()
    # fig, axs = plt.subplots(3)
    # fig.suptitle('coefitionts')
    # axs[0].plot(Regretor.Lambda, Regretor.A_coefs)
    # axs[0].set_ylabel(r'$A(\tilde{\Lambda})$')
    # axs[0].set_xlabel(r'$\tilde{\Lambda}$')
    # axs[1].plot(Regretor.Lambda, Regretor.B_coefs)
    # axs[1].set_ylabel(r'$B(\tilde{\Lambda})$')
    # axs[1].set_xlabel(r'$\tilde{\Lambda}$')
    # axs[2].plot(Regretor.Lambda, Regretor.scores)
    # axs[2].set_ylabel(r'$R^2$')
    # axs[2].set_xlabel(r'$\tilde{\Lambda}$')
    # plt.tight_layout()
    
    # plt.figure()
    # fig, axs = plt.subplots(2)
    # axs[0].plot(Regretor.Lambda, Regretor.A_coefs, label = 'real')
    # axs[0].plot(Regretor.Lambda, 0.327 * np.ones(np.shape(Regretor.Lambda)), label = 'zero order approximation')
    # axs[0].plot(Regretor.Lambda, 0.3045097 - 0.01790584 * Regretor.Lambda + 0.26870063 * Regretor.Lambda ** 2, label = 'second order approximation')
    # axs[0].set_ylabel(r'$A(\tilde{\Lambda})$')
    # axs[0].set_xlabel(r'$\tilde{\Lambda}$')
    # axs[0].legend()
    
    # axs[1].plot(Regretor.Lambda, Regretor.B_coefs, label = 'numeric')
    # axs[1].set_ylabel(r'$B(\tilde{\Lambda})$')
    # axs[1].set_xlabel(r'$\tilde{\Lambda}$')
    # axs[1].plot(Regretor.Lambda, - 1.15 * np.log( Regretor.Lambda) - 0.95, label = 'approximation')
    # axs[1].legend()
    



if __name__ == '__main__':
    phase = main()
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 18:55:46 2020

@author: Jonathan Browning
"""

import numpy as np
from scipy.stats import gaussian_kde as kdf
from scipy import special as sp

class ShadowedRice:
    numSamples = 2*(10**6)  # the number of samples used in the simulation
    r = np.linspace(0, 6, 6000) # theoretical envelope PDF x axes
    theta = np.linspace(-np.pi, np.pi, 6000)    # theoretical phase PDF x axes
    
    def __init__(self, K, m, r_hat, phi):
        
        # user input checks and assigns value
        self.K = self.input_Check(K, "K", 0.001, 50)
        self.m = self.input_Check(m, "m", 0.001, 50)
        self.r_hat = self.input_Check(r_hat, "\hat{r}", 0.5, 2.5)  
        self.phi = self.input_Check(phi, "\phi", -np.pi, np.pi)
        
        # simulating and theri densities 
        self.multipathFading = self.complex_Multipath_Fading()
        self.xdataEnv, self.ydataEnv = self.envelope_Density()
        self.xdataPh, self.ydataPh = self.phase_Density()
        
        # theoretical PDFs calculated
        self.envelopeProbability = self.envelope_PDF()
        self.phaseProbability = self.phase_PDF()

    def input_Check(self, data, inputName, lower, upper):
        # input_Check checks the user inputs
        
        # has a value been entered
        if data == "":
            raise ValueError(" ".join((inputName, "must have a numeric value")))
        
        # incase of an non-numeric input 
        try:
            data = float(data)
        except:    
            raise ValueError(" ".join((inputName, "must have a numeric value")))
    
        # data must be within the range
        if data < lower or data > upper:
            raise ValueError(" ".join((inputName, f"must be in the range [{lower:.2f}, {upper:.2f}]")))
        
        return data

    def calculate_Means(self):
        # calculate_means calculates the means of the complex Gaussians representing the
        # in-phase and quadrature components
        
        p = np.sqrt(self.K / (1+self.K)) * self.r_hat * np.cos(self.phi)
        q = np.sqrt(self.K / (1+self.K))  * self.r_hat * np.sin(self.phi)
        
        return p, q
    
    def scattered_Component(self):
        # scattered_Component calculates the power of the scattered signal component
        
        sigma = self.r_hat / np.sqrt( 2 * (1+self.K) )
        
        return sigma
    
    def generate_Gaussians(self, mean, sigma):
        # generate_Gaussians generates the Gaussian random variables
        
        gaussians = np.random.default_rng().normal(mean, sigma, self.numSamples)
        
        return gaussians
    
    def complex_Multipath_Fading(self):
        # complex_Multipath_Fading generates the complex fading random variables
        
        p, q = self.calculate_Means()
        sigma = self.scattered_Component()
        xi = np.sqrt(np.random.gamma(self.m, 1/self.m, self.numSamples))
        
        multipathFading = self.generate_Gaussians(xi*p, sigma) + (1j*self.generate_Gaussians(xi*q, sigma))
        
        return multipathFading
    
    def envelope_PDF(self):
        # envelope_PDF calculates the theoretical envelope PDF
        
        PDF = 2 * (1+self.K) * self.r *(self.m**(self.m)) / (self.r_hat**(2)*(self.m+self.K)**(self.m)) \
                * np.exp(- ((1+self.K) * self.r**(2)) / self.r_hat**(2)) \
                * sp.hyp1f1(self.m, 1, self.r**(2)*self.K*(self.K+1)/(self.r_hat**(2)*(self.K+self.m)))            
            
        return PDF

    def phase_PDF(self):
        # phase_PDF calculates the theoretical phase PDF
        
        PDF =  (self.m**self.m * np.sqrt(self.K)/(2 * np.sqrt(np.pi) * (self.K + self.m)**(self.m +1/2))) \
            * ( np.sqrt((self.K +self.m)/(np.pi*self.K)) * sp.hyp2f1(self.m, 1, 1/2,  (self.K*(np.cos(self.theta - self.phi))**(2))/(self.K +self.m)) \
            +  ((sp.gamma(self.m+1/2) / sp.gamma(self.m))*np.cos(self.theta-self.phi) \
            * (1-  (self.K*(np.cos(self.theta - self.phi))**(2)) / (self.K +self.m))**(-self.m-1/2)))
        
        return PDF
        
    
    def envelope_Density(self):
        # envelope_Density finds the envelope PDF of the simulated random variables
        
        R = np.sqrt((np.real(self.multipathFading))**2 + (np.imag(self.multipathFading))**2)
        kde = kdf(R)
        x = np.linspace(R.min(), R.max(), 100)
        p = kde(x)
        
        return x, p
    
    def phase_Density(self):
        # phase_Density finds the phase PDF of the simulated random variables
        
        R = np.angle(self.multipathFading)
        kde = kdf(R)
        x = np.linspace(R.min(), R.max(), 100)
        p = kde(x)
        
        return x, p
   




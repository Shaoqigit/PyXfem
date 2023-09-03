# This field is retrived from part of pymls, a software distributed under the MIT license.
# For any question, please contact one of the authors cited below.
#
# Copyright (c) 2017
# 	Olivier Dazel <olivier.dazel@univ-lemans.fr>
# 	Mathieu Gaborit <gaborit@kth.se>
# 	Peter Göransson <pege@kth.se>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#

import os
import numpy as np
import numpy.linalg as LA
import matplotlib.pyplot as plt
from pymls import from_yaml, Solver, Layer, backing

from pymls import Solver, Layer, backing
from fem.materials import Air, Fluid, EquivalentFluid, PoroElasticMaterial


phi          = 0.99  # porosity
sigma        = 1.0567e4 # resistivity
alpha        = 1.2  # Tortuosity
Lambda_prime = 490e-6  # Viscous characteristic length
Lambda      = 240e-6  # 
rho_1 = 9.2
nu = 0.285
E = 3.155e5
eta = 0.032
mat = PoroElasticMaterial('xfm', phi, sigma, alpha, Lambda_prime, Lambda, rho_1, E, nu, eta)

air = Air('classical air')

sigmaF = 0.
d = 1e-3
freq = 2000
omega = 2*np.pi*freq
theta_d = 0.
theta = theta_d*(np.pi/180)
mat.set_frequency(omega)

k_0 = omega/Air.c
ky = k_0*np.sin(theta_d*np.pi/180)
kx = k_0*np.cos(theta_d*np.pi/180)
M = np.zeros((2,2), dtype=complex)


M[0,0] = 1.
M[0,1] = 1j * omega * sigmaF * d
M[1,0] = 0.
M[1,1] = 1.

kx_1=np.sqrt(mat.delta_1**2-ky**2)
kx_2=np.sqrt(mat.delta_2**2-ky**2)
kx_3=np.sqrt(mat.delta_3**2-ky**2)

# # linear system matrix
SV =np.zeros((4,4),dtype=complex)
# u^t
SV[0,:] = np.array([1j*kx/(Air.rho*omega**2), -M[1,0]+M[1,1]*1j*kx_1/(mat.K_eq_til*mat.delta_1**2), 
                    -M[1,0]+M[1,1]*1j*kx_2/(mat.K_eq_til*mat.delta_2**2), M[1,1]*1j*ky*mat.mu_3])
# p
# SV[1,:] = np.array([1+sigmaF*d*kx/(omega*Air.rho), -1, -1, 0])
SV[1,:] = np.array([1, -M[0,0]+M[0,1]*1j*kx_1*mat.mu_1/(mat.mu_1*mat.K_eq_til*mat.delta_1**2), 
                    -M[0,0]+M[0,1]*1j*kx_2*mat.mu_2/(mat.mu_2*mat.K_eq_til*mat.delta_2**2), 
                    0+M[0,1]*1j*ky*mat.mu_3])
SV[2,:] = np.array([0, 2*mat.N*kx_1**2/(mat.mu_1*mat.K_eq_til*mat.delta_1**2)+mat.A_hat/(mat.mu_1*mat.K_eq_til),
                    2*mat.N*kx_2**2/(mat.mu_2*mat.K_eq_til*mat.delta_2**2)+mat.A_hat/(mat.mu_2*mat.K_eq_til),
                    2*mat.N*ky*kx_3])
SV[3,:] = np.array([0, -2*mat.N*kx_1*ky/(mat.mu_1*mat.K_eq_til*mat.delta_1**2), 
                    -2*mat.N*kx_2*ky/(mat.mu_2*mat.K_eq_til*mat.delta_2**2),
                    +(mat.N*kx_3**2-mat.N*ky**2)])
# vector
F = np.zeros(4,dtype=complex)
F[0] = 1j*kx/(Air.rho*omega**2)
# F[1] = -1+sigmaF*d*kx/(omega*Air.rho)
F[1] = -1

# unknown coefficient
C = LA.solve(SV,F)

P_i = lambda x, y: np.exp(-1j*kx*x-1j*ky*y)

P_r = lambda x, y: C[0]*np.exp(1j*kx*x-1j*ky*y)

P_a = lambda x, y: 1.*(np.exp(-1j*kx*x-1j*ky*y) + C[0]*np.exp(1j*kx*x-1j*ky*y))

u_a = lambda x, y: 1./(Air.rho*omega**2)*((-1j*kx*np.exp(-1j*kx*x-1j*ky*y) + 1j*kx*C[0]*np.exp(1j*kx*x-1j*ky*y)))

P_t = lambda x, y: 1*(C[1]*np.exp(-1j*kx_1*x-1j*ky*y) + C[2]*np.exp(-1j*kx_2*x-1j*ky*y))

us_x = lambda x, y: 1*(-1j*kx_1*C[1]/(mat.mu_1*mat.K_eq_til*mat.delta_1**2)*np.exp(-1j*kx_1*x-1j*ky*y) -
                      1j*kx_2*C[2]/(mat.mu_2*mat.K_eq_til*mat.delta_2**2)*np.exp(-1j*kx_2*x-1j*ky*y) -
                      1j*C[3]*ky*np.exp(-1j*kx_3*x-1j*ky*y))

u_t = lambda x, y: 1*(-1j*mat.mu_1*kx_1*C[1]/(mat.mu_1*mat.K_eq_til*mat.delta_1**2)*np.exp(-1j*kx_1*x-1j*ky*y) -
                      1j*mat.mu_2*kx_2*C[2]/(mat.mu_2*mat.K_eq_til*mat.delta_2**2)*np.exp(-1j*kx_2*x-1j*ky*y) -
                      1j*mat.mu_3*C[3]*ky*np.exp(-1j*kx_3*x-1j*ky*y)) - 1/(omega**2*mat.rho_f)*P_t(x,y)

us_y = lambda x, y: (-1j*ky*C[1]/(mat.mu_1*mat.K_eq_til*mat.delta_1**2)*np.exp(-1j*kx_1*x-1j*ky*y) -
                      1j*ky*C[2]/(mat.mu_2*mat.K_eq_til*mat.delta_2**2)*np.exp(-1j*kx_2*x-1j*ky*y) +
                      1j*C[3]*kx_3*np.exp(-1j*kx_3*x-1j*ky*y))

e_xx = lambda x, y: (-kx_1**2*C[1]/(mat.mu_1*mat.K_eq_til*mat.delta_1**2)*np.exp(-1j*kx_1*x-1j*ky*y) -
        kx_2**2*C[2]/(mat.mu_2*mat.K_eq_til*mat.delta_2**2)*np.exp(-1j*kx_2*x-1j*ky*y) -C[3]*
        ky*kx_3*np.exp(-1j*kx_3*x-1j*ky*y))

e_yy = lambda x, y: (-ky**2*C[1]/(mat.mu_1*mat.K_eq_til*mat.delta_1**2)*np.exp(-1j*kx_1*x-1j*ky*y) -
        ky**2*C[2]/(mat.mu_2*mat.K_eq_til*mat.delta_2**2)*np.exp(-1j*kx_2*x-1j*ky*y) +C[3]*
        ky*kx_3*np.exp(-1j*kx_3*x-1j*ky*y))

e_xy = lambda x, y: 0.5*(-2*C[1]*kx_1*ky/(mat.mu_1*mat.K_eq_til*mat.delta_1**2)*np.exp(-1j*kx_1*x-1j*ky*y)- 
                          2*C[2]*kx_2*ky/(mat.mu_2*mat.K_eq_til*mat.delta_2**2)*np.exp(-1j*kx_2*x-1j*ky*y)+
                          (kx_3**2-ky**2)*C[3]*np.exp(-1j*kx_3*x-1j*ky*y))

sigma_xx = lambda x, y: 2*mat.N*e_xx(x,y) + mat.A_hat*(e_xx(x,y)+e_yy(x,y))
sigma_yy = lambda x, y: 2*mat.N*e_yy(x,y) + mat.A_hat*(e_xx(x,y)+e_yy(x,y))
sigma_xy = lambda x, y: 2*mat.N*e_xy(x,y)



def Fluid_Biot_Pressure(X, Y):
    pre = np.zeros((len(X), len(Y)), dtype=complex)
    for j, y in enumerate(Y):
        for i, x in enumerate(X):    
            if x <= 0:
                pre[i,j] = P_a(x, y)
            
            else:
                pre[i,j] = P_t(x, y)
    return pre

def Fluid_Biot_Displacement(X, Y):
    Usx = np.zeros((len(X), len(Y)), dtype=complex)
    for j, y in enumerate(Y):
        for i, x in enumerate(X):    
            if x < 0:
                Usx[i,j] = 0
            
            else:
                Usx[i,j] = us_x(x, y)
    return Usx

def Fluid_Biot_Ut(X, Y):
    ux = np.zeros((len(X), len(Y)), dtype=complex)
    for j, y in enumerate(Y):
        for i, x in enumerate(X):    
            if x <= 0:
                ux[i,j] = u_a(x, y)
            
            else:
                ux[i,j] = u_t(x, y)
    return ux
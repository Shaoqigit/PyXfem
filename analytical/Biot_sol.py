#! /usr/bin/env python
# -*- coding:utf8 -*-
#
# FEM_1D_porous.py
#
# This file is part of pymls, a software distributed under the MIT license.
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
from mediapack import Air, PEM, EqFluidJCA

from scipy.special import legendre
Air = Air()
plot_color = ['r', 'm', 'b', 'k', 'g', 'c', 'y']
plot_color +=plot_color

def PEM_ks(mat,ky):
    ''' S={0:\hat{\sigma}_{xy}, 1:u_x^s, 2:u_x^t, 3:\hat{\sigma}_{xx}, 4:p, 5:u_y^s}'''
    kx_1=np.sqrt(mat.delta_1**2-ky**2)
    kx_2=np.sqrt(mat.delta_2**2-ky**2)
    kx_3=np.sqrt(mat.delta_3**2-ky**2)

    kx = np.array([kx_1, kx_2, kx_3])
    delta = np.array([mat.delta_1,mat.delta_2,mat.delta_3])
    return kx, delta



def PEM_SV(mat,ky):
    ''' S={0:\hat{\sigma}_{xy}, 1:u_x^s, 2:u_x^t, 3:\hat{\sigma}_{xx}, 4:p, 5:u_y^s}'''
    kx, delta = PEM_ks(mat,ky)

    alpha_1=-1j*mat.A_hat*mat.delta_1**2-1j*2*mat.N*kx[0]**2
    alpha_2=-1j*mat.A_hat*mat.delta_2**2-1j*2*mat.N*kx[1]**2
    alpha_3= -2*1j*mat.N*kx[2]*ky



    SV =np.zeros((6,6),dtype=complex)
    SV[0:6,0] = np.array([-2*1j*mat.N*kx[0]*ky, kx[0], mat.mu_1*kx[0], alpha_1, 1j*delta[0]**2*mat.K_eq_til*mat.mu_1,ky])
    SV[0:6,3] = np.array([ 2*1j*mat.N*kx[0]*ky,-kx[0],-mat.mu_1*kx[0], alpha_1, 1j*delta[0]**2*mat.K_eq_til*mat.mu_1,ky])

    SV[0:6,1] = np.array([-2*1j*mat.N*kx[1]*ky, kx[1], mat.mu_2*kx[1],alpha_2, 1j*delta[1]**2*mat.K_eq_til*mat.mu_2,ky])
    SV[0:6,4] = np.array([ 2*1j*mat.N*kx[1]*ky,-kx[1],-mat.mu_2*kx[1],alpha_2, 1j*delta[1]**2*mat.K_eq_til*mat.mu_2,ky])

    SV[0:6, 2] = np.array([1j*mat.N*(kx[2]**2-ky**2), ky, mat.mu_3*ky, alpha_3, 0., -kx[2]])
    SV[0:6, 5] = np.array([1j*mat.N*(kx[2]**2-ky**2), ky, mat.mu_3*ky, -alpha_3, 0., kx[2]])
    return SV

def solve_PW(mat,ky,x,type_CL):

    # S={0:\hat{\sigma}_{xy},1:u_x^s,2:u_x^t,3:\hat{\sigma}_{xx},4:p,5:u_y^s}
    if type_CL == 1 :
        BC = [2,0,3,1,2,5]
        BV = [-1.,0,0,0,0,0]
    if type_CL == 2 :
        BC = [3,0,2,1,2,5]
        BV = [1.,0,0,0,0,0]
    if type_CL == 3 :
        BC = [0,2,3,1,2,5]
        BV = [1.,0,0,0,0,0]

    L=x[-1]-x[0]
    x = x-x[0]-L # So that x = [–L;0]
    beta_1=np.sqrt(mat.delta_1**2-ky**2)
    beta_2=np.sqrt(mat.delta_2**2-ky**2)
    beta_3=np.sqrt(mat.delta_3**2-ky**2)

    k_z = np.array([beta_1,beta_2,beta_3])
    SV = PEM_SV(mat,ky)
    Mat_PW = np.zeros((6,6),dtype=complex)
    F_PW = np.zeros(6,dtype=complex)
    # S={\hat{\sigma}_{xy},u_x^s,u_x^t,\hat{\sigma}_{xx},p,u_y^s}
    Mat_PW[0,0]=SV[BC[0],0]
    Mat_PW[0,1]=SV[BC[0],1]
    Mat_PW[0,2]=SV[BC[0],2]
    Mat_PW[0,3]=SV[BC[0],3]*np.exp(-1j*k_z[0]*L)
    Mat_PW[0,4]=SV[BC[0],4]*np.exp(-1j*k_z[1]*L)
    Mat_PW[0,5]=SV[BC[0],5]*np.exp(-1j*k_z[2]*L)
    F_PW[0] = BV[0]

    Mat_PW[1,0]=SV[BC[1],0]
    Mat_PW[1,1]=SV[BC[1],1]
    Mat_PW[1,2]=SV[BC[1],2]
    Mat_PW[1,3]=SV[BC[1],3]*np.exp(-1j*k_z[0]*L)
    Mat_PW[1,4]=SV[BC[1],4]*np.exp(-1j*k_z[1]*L)
    Mat_PW[1,5]=SV[BC[1],5]*np.exp(-1j*k_z[2]*L)
    F_PW[1] = BV[1]

    Mat_PW[2,0]=SV[BC[2],0]
    Mat_PW[2,1]=SV[BC[2],1]
    Mat_PW[2,2]=SV[BC[2],2]
    Mat_PW[2,3]=SV[BC[2],3]*np.exp(-1j*k_z[0]*L)
    Mat_PW[2,4]=SV[BC[2],4]*np.exp(-1j*k_z[1]*L)
    Mat_PW[2,5]=SV[BC[2],5]*np.exp(-1j*k_z[2]*L)
    F_PW[2] = BV[2]

    Mat_PW[3,0]=SV[BC[3],0]*np.exp(-1j*k_z[0]*L)
    Mat_PW[3,1]=SV[BC[3],1]*np.exp(-1j*k_z[1]*L)
    Mat_PW[3,2]=SV[BC[3],2]*np.exp(-1j*k_z[2]*L)
    Mat_PW[3,3]=SV[BC[3],3]
    Mat_PW[3,4]=SV[BC[3],4]
    Mat_PW[3,5]=SV[BC[3],5]
    F_PW[3] = BV[3]

    Mat_PW[4,0]=SV[BC[4],0]*np.exp(-1j*k_z[0]*L)
    Mat_PW[4,1]=SV[BC[4],1]*np.exp(-1j*k_z[1]*L)
    Mat_PW[4,2]=SV[BC[4],2]*np.exp(-1j*k_z[2]*L)
    Mat_PW[4,3]=SV[BC[4],3]
    Mat_PW[4,4]=SV[BC[4],4]
    Mat_PW[4,5]=SV[BC[4],5]
    F_PW[4] = BV[4]

    Mat_PW[5,0]=SV[BC[5],0]*np.exp(-1j*k_z[0]*L)
    Mat_PW[5,1]=SV[BC[5],1]*np.exp(-1j*k_z[1]*L)
    Mat_PW[5,2]=SV[BC[5],2]*np.exp(-1j*k_z[2]*L)
    Mat_PW[5,3]=SV[BC[5],3]
    Mat_PW[5,4]=SV[BC[5],4]
    Mat_PW[5,5]=SV[BC[5],5]
    F_PW[5] = BV[5]

    X = LA.solve(Mat_PW,F_PW)

    sol_PW = np.zeros((6,len(x)),dtype=complex)
    for ix in range(len(x)):
            sol_PW[:,ix] =  SV[:,0]*np.exp(-1j*k_z[0]*(x[ix]+L))*X[0]
            sol_PW[:,ix] += SV[:,1]*np.exp(-1j*k_z[1]*(x[ix]+L))*X[1]
            sol_PW[:,ix] += SV[:,2]*np.exp(-1j*k_z[2]*(x[ix]+L))*X[2]
            sol_PW[:,ix] += SV[:,3]*np.exp( 1j*k_z[0]*(x[ix]  ))*X[3]
            sol_PW[:,ix] += SV[:,4]*np.exp( 1j*k_z[1]*(x[ix]  ))*X[4]
            sol_PW[:,ix] += SV[:,5]*np.exp( 1j*k_z[2]*(x[ix]  ))*X[5]
    return sol_PW

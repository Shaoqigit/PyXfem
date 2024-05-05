# This file is part of PyXfem, a software distributed under the MIT license.
# For any question, please contact the authors cited below.
#
# Copyright (c) 2023
# 	Shaoqi WU <shaoqiwu@outlook.com>
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

# precompte and store the elementary matrices for 1D Lobatto elements

import numpy as np
from SAcouS.acxfem.quadratures import GaussLegendre2DTri
from SAcouS.acxfem.polynomial import  Lagrange2DTri
# from numpy.polynomial.legendre import leggauss

order = 1
n_pts = 4
lagrange_o1 = Lagrange2DTri(order)
gl_intg_o1 = GaussLegendre2DTri(n_pts)
gl_pts, gl_wts = gl_intg_o1.points(), gl_intg_o1.weights()
N_o1 = lagrange_o1.get_shape_functions()
B_o1 = lagrange_o1.get_der_shape_functions()

import pdb; pdb.set_trace()

N_p1_1 = np.array([N_o1[0](gl_pts[i][0], gl_pts[i][1]) for i in range(n_pts)])
N_p1_2 = np.array([N_o1[1](gl_pts[i][0], gl_pts[i][1]) for i in range(n_pts)])
N_p1_3 = np.array([N_o1[2](gl_pts[i][0], gl_pts[i][1]) for i in range(n_pts)])
N_p1 = np.array([N_p1_1, N_p1_2, N_p1_3]).T
Ke_p1 = N_p1.T @ np.diag(gl_wts) @ N_p1

B_p1_11 = np.array([B_o1[0][0](gl_pts[i][0], gl_pts[i][1]) for i in range(n_pts)])
B_p1_12 = np.array([B_o1[0][1](gl_pts[i][0], gl_pts[i][1]) for i in range(n_pts)])
B_p1_21 = np.array([B_o1[1][0](gl_pts[i][0], gl_pts[i][1]) for i in range(n_pts)])
B_p1_22 = np.array([B_o1[1][1](gl_pts[i][0], gl_pts[i][1]) for i in range(n_pts)])
B_p1_31 = np.array([B_o1[2][0](gl_pts[i][0], gl_pts[i][1]) for i in range(n_pts)])
B_p1_32 = np.array([B_o1[2][1](gl_pts[i][0], gl_pts[i][1]) for i in range(n_pts)])
B_p1 = np.array([[B_p1_11.T, B_p1_12.T], [B_p1_21.T, B_p1_22.T], [B_p1_31.T, B_p1_32.T]])
Me_p1 = B_p1.T @ np.diag(gl_wts) @ B_p1

def compute_matrix(Ke, Me, Ce, order):
    # n_pts > (p+1)/2
    if order == 1:
        n_pts = 3
    elif order == 2:
        n_pts = 4
    else:
        raise ValueError("order should be 1 or 2")

    gl_integrate = GaussLegendre2DTri(n_pts)
    gl_pts, gl_wts = gl_integrate.points(), gl_integrate.weights()
    # gl_pts, gl_wts = leggauss(n_pts)
    l = Lagrange2DTri(order)
    B = l.get_der_shape_functions()
    N = l.get_shape_functions()
    import pdb; pdb.set_trace()
    len(Me[0])
    for i in range(len(Me[0])):
        for j in range(len(Me[0])):
            Ke[i, j] = sum(gl_wt*B[i](gl_pt)*B[j](gl_pt) for gl_pt, gl_wt in zip(gl_pts, gl_wts))
            Me[i, j] = sum(gl_wt*N[i](gl_pt)*N[j](gl_pt) for gl_pt, gl_wt in zip(gl_pts, gl_wts))
            Ce[i, j] = sum(gl_wt*N[i](gl_pt)*B[j](gl_pt) for gl_pt, gl_wt in zip(gl_pts, gl_wts))  #Coupling matrix that to be used in Biot equation
            if abs(Ke[i, j]) < 1e-10:
                Ke[i, j] = 0
            if abs(Me[i, j]) < 1e-10:
                Me[i, j] = 0
            if abs(Ce[i, j]) < 1e-10:
                Ce[i, j] = 0
            


# print(Ke1Do2)
# print(Me1Do2)
# print(Ke1Do3)
# print(Me1Do3)
# print(Ke1Do4)
# print(Me1Do4)
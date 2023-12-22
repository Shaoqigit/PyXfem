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
from acxfem.quadratures import GaussLegendreQuadrature
from acxfem.polynomial import Lobatto, Larange, Lagrange2DTri
# from numpy.polynomial.legendre import leggauss

def compute_matrix(Ke, Me, Ce, order):
    n_pts = order*2
    gl_q = GaussLegendreQuadrature(n_pts)
    gl_pts, gl_wts = gl_q.points(), gl_q.weights()
    # gl_pts, gl_wts = leggauss(n_pts)
    l = Lagrange2DTri(order)
    B = l.get_der_shape_functions()
    N = l.get_shape_functions()
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
            


# 2D Triangular lagrange element matrix: p=1
order = 1
Ke1Do1 = np.zeros((3,3))
Me1Do1 = np.zeros((3,3))
Ce1Do1 = np.zeros((3,3))
compute_matrix(Ke1Do1, Me1Do1, Ce1Do1, order)


# 1D lobatto element matrix: p=2
order = 2
Ke1Do2 = np.zeros((3,3))
Me1Do2 = np.zeros((3,3))
Ce1Do2 = np.zeros((3,3))
compute_matrix(Ke1Do2, Me1Do2, Ce1Do2, order)


Ke1D = [Ke1Do1, Ke1Do2]
Me1D = [Me1Do1, Me1Do2]
Ce1D = [Ce1Do1, Ce1Do2]

# print(Ke1Do1)
# print(Me1Do1)
# print(Ke1Do2)
# print(Me1Do2)
# print(Ke1Do3)
# print(Me1Do3)
# print(Ke1Do4)
# print(Me1Do4)
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

def compute_matrix(order):
    n_quad_pts = order*2 + 1
    lagrange = Lagrange2DTri(order)
    points, weights = GaussLegendre2DTri(n_quad_pts).points(), GaussLegendre2DTri(n_quad_pts).weights()
    N_w = np.array([weight*lagrange.get_shape_functions(*point) for point, weight in zip(points, weights)])
    N = np.array([lagrange.get_shape_functions(*point) for point in points])
    import pdb;pdb.set_trace()
    Me = N_w.T @ N

    B_w = np.array([weight*lagrange.get_der_shape_functions(*point) for point, weight in zip(points, weights)])
    B = np.array([lagrange.get_der_shape_functions(*point) for point in points])
    Ke = B_w[:,:,0].T @ B[:,:,0] *lagrange.determinant_jacobi + B_w[:,:,1].T @ B[:,:,1] *lagrange.determinant_jacobi

    return Ke, Me, N, B


# 2D lobatto element matrix: p=1
order = 1
Ke2Do1, Me2Do1 = compute_matrix(order)
Ke2DTri = [Ke2Do1]
Me2DTri = [Me2Do1]
# import numpy as np

# # Quadrature points and weights for the reference triangle
# quad_points = np.array([
#     [1/2, 1/2],
#     [0, 1/2],
#     [1/2, 0]
# ])
# quad_weights = np.array([1/6, 1/6, 1/6])

# # Derivatives of shape functions in the reference coordinates
# dN_dxi = np.array([
#     [-1, 1, 0],
#     [-1, 0, 1]
# ])

# # Vertices of the reference triangle
# vertices = np.array([
#     [0, 0],
#     [1, 0],
#     [0, 1]
# ])

# # Jacobian matrix and its determinant
# J = np.dot(dN_dxi, vertices)
# det_J = np.linalg.det(J)

# # Compute the inverse of the Jacobian matrix
# inv_J = np.linalg.inv(J)

# # Derivatives of shape functions in the physical coordinates
# dN_dx = np.dot(inv_J.T, dN_dxi)

# # Initialize the stiffness matrix
# stiffness_matrix = np.zeros((3, 3))

# # Perform numerical integration
# for i in range(3):
#     for j in range(3):
#         integral = 0
#         for k in range(len(quad_weights)):
#             weight = quad_weights[k]
#             # Evaluate B_i^T B_j
#             B_i = dN_dx[:, i]
#             B_j = dN_dx[:, j]
#             integral += np.dot(B_i, B_j) * weight
#         stiffness_matrix[i, j] = integral * det_J
# import pdb;pdb.set_trace()
# print("Stiffness matrix:\n", stiffness_matrix)

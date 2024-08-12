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
from .Quadratures import GaussLegendre2DTri, GaussLegendre3DTetra
from .Polynomial import Lagrange2DTri, Lagrange2DQuad, Lagrange3DTetra

lag2d_poly_o1 = Lagrange2DTri(1)
points_o1, weights_o1 = GaussLegendre2DTri(3).points(), GaussLegendre2DTri(
    3).weights()
# N_o1 = np.array(
#     [lag2d_poly_o1.get_shape_functions(*point) for point in points_o1])
# B_o1 = np.array(
#     [lag2d_poly_o1.get_der_shape_functions(*point) for point in points_o1])

#
N_tri_o1 = np.array([[0.66666667, 0.16666667, 0.16666667],
                     [0.16666667, 0.66666667, 0.16666667],
                     [0.16666667, 0.16666667, 0.66666667]])
B_tri_o1 = np.array([[[-1., -1.], [1., 0.], [0., 1.]],
                     [[-1., -1.], [1., 0.], [0., 1.]],
                     [[-1., -1.], [1., 0.], [0., 1.]]])

lag3d_poly_o1 = Lagrange3DTetra(1)
points_tetra_o1, weights_tetra_o1 = GaussLegendre3DTetra(
    4).points(), GaussLegendre3DTetra(4).weights()
N_tetra_o1 = np.array(
    [lag3d_poly_o1.get_shape_functions(*point) for point in points_tetra_o1])
B_tetra_o1 = np.array([
    lag3d_poly_o1.get_der_shape_functions(*point) for point in points_tetra_o1
])


def get_N_B_p1(dim: int) -> np.ndarray:
  if dim == 2:
    return N_tri_o1, B_tri_o1
  elif dim == 3:
    return N_tetra_o1, B_tetra_o1
  else:
    raise ValueError("Only support 2D and 3D")


if __name__ == "__main__":
  print(N_tetra_o1)
  print(B_tetra_o1)

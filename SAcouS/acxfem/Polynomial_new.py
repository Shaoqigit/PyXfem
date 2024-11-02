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

# polynomal definition, mainly developed for Lobatto shape functions
# partially support larange shape functions

import numpy as np
from abc import abstractmethod, ABCMeta


class BasePolynomial(metaclass=ABCMeta):
  """base abstract polynomial class
    parameters:
    order: int
        degree of the polynomial
    x: ndarray
        position of the Gauss points
        
    retruns:
    p_lobatto: ndarray
        value of the Lobatto polynomial at x
    d_lobatto: ndarray
        value of the Lobatto polynomial first derivative at x
    """
  __slots__ = ['order']

  def __init__(self, order):
    self.order = order

  @abstractmethod
  def polynomial(self):
    pass

  @abstractmethod
  def derivative(self):
    pass

  @abstractmethod
  def get_shape_functions(self):
    pass

  def get_der_shape_functions(self):
    pass

  def __call__(self):
    return self.polynomial()

  def __len__(self):
    return self.order

  def __getitem__(self, key):
    return self.polynomial()[key]


class Lobatto(BasePolynomial):
  """lobatto polynomial class"""

  def polynomial(self, u):
    if self.order == 1:
      lobatto = np.array([
          -1 / 2 * u + 1 / 2,
          1 / 2 * u + 1 / 2,
      ])
    elif self.order == 2:
      lobatto = np.array([
          -1 / 2 * u + 1 / 2, 1 / 2 * u + 1 / 2,
          1 / (6**0.5) * (1.5 * (u - 1) * (u + 1))
      ])
    elif self.order == 3:
      lobatto = np.array([
          -1 / 2 * u + 1 / 2, 1 / 2 * u + 1 / 2,
          1 / (6**0.5) * (1.5 * (u - 1) * (u + 1)),
          1 / (10**0.5) * (5 * u * (u - 1) * (u + 1) / 2)
      ])
    elif self.order == 4:
      lobatto = np.array([
          -1 / 2 * u + 1 / 2, 1 / 2 * u + 1 / 2,
          1 / (6**0.5) * (1.5 * (u - 1) * (u + 1)),
          1 / (10**0.5) * (5 * u * (u - 1) * (u + 1) / 2),
          1 / (14**0.5) * (7 * (u - 1) * (u + 1) * (5 * u**2 - 1) / 8)
      ])
    elif self.order == 5:
      lobatto = np.array([
          -1 / 2 * u + 1 / 2, 1 / 2 * u + 1 / 2,
          1 / (6**0.5) * (1.5 * (u - 1) * (u + 1)),
          1 / (10**0.5) * (5 * u * (u - 1) * (u + 1) / 2),
          1 / (14**0.5) * (7 * (u - 1) * (u + 1) * (5 * u**2 - 1) / 8),
          1 / (18**0.5) * (9 * u * (u - 1) * (u + 1) * (7 * u**2 - 3) / 8)
      ])
    elif self.order == 6:
      lobatto = np.array([
          -1 / 2 * u + 1 / 2, 1 / 2 * u + 1 / 2,
          1 / (6**0.5) * (1.5 * (u - 1) * (u + 1)),
          1 / (10**0.5) * (5 * u * (u - 1) * (u + 1) / 2),
          1 / (14**0.5) * (7 * (u - 1) * (u + 1) * (5 * u**2 - 1) / 8),
          1 / (18**0.5) * (9 * u * (u - 1) * (u + 1) * (7 * u**2 - 3) / 8),
          1 / (22**0.5) * (11 * (u - 1) * (u + 1) *
                           (21 * u**4 - 14 * u**2 + 1) / 16)
      ])
    elif self.order == 7:
      lobatto = np.array([
          -1 / 2 * u + 1 / 2, 1 / 2 * u + 1 / 2,
          1 / (6**0.5) * (1.5 * (u - 1) * (u + 1)),
          1 / (10**0.5) * (5 * u * (u - 1) * (u + 1) / 2),
          1 / (14**0.5) * (7 * (u - 1) * (u + 1) * (5 * u**2 - 1) / 8),
          1 / (18**0.5) * (9 * u * (u - 1) * (u + 1) * (7 * u**2 - 3) / 8),
          1 / (22**0.5) * (11 * (u - 1) * (u + 1) *
                           (21 * u**4 - 14 * u**2 + 1) / 16), 1 / (26**0.5) *
          (13 * u * (u - 1) * (u + 1) * (33 * u**4 - 30 * u**2 + 5) / 16)
      ])
    elif self.order == 8:
      lobatto = np.array([
          -1 / 2 * u + 1 / 2, 1 / 2 * u + 1 / 2,
          1 / (6**0.5) * (1.5 * (u - 1) * (u + 1)),
          1 / (10**0.5) * (5 * u * (u - 1) * (u + 1) / 2),
          1 / (14**0.5) * (7 * (u - 1) * (u + 1) * (5 * u**2 - 1) / 8),
          1 / (18**0.5) * (9 * u * (u - 1) * (u + 1) * (7 * u**2 - 3) / 8),
          1 / (22**0.5) * (11 * (u - 1) * (u + 1) *
                           (21 * u**4 - 14 * u**2 + 1) / 16), 1 / (26**0.5) *
          (13 * u * (u - 1) * (u + 1) * (33 * u**4 - 30 * u**2 + 5) / 16),
          1 / (30**0.5) * (15 * (u - 1) * (u + 1) *
                           (429 * u**6 - 495 * u**4 + 135 * u**2 - 5) / 128)
      ])
    elif self.order == 9:
      lobatto = np.array([
          -1 / 2 * u + 1 / 2, 1 / 2 * u + 1 / 2,
          1 / (6**0.5) * (1.5 * (u - 1) * (u + 1)),
          1 / (10**0.5) * (5 * u * (u - 1) * (u + 1) / 2),
          1 / (14**0.5) * (7 * (u - 1) * (u + 1) * (5 * u**2 - 1) / 8),
          1 / (18**0.5) * (9 * u * (u - 1) * (u + 1) * (7 * u**2 - 3) / 8),
          1 / (22**0.5) * (11 * (u - 1) * (u + 1) *
                           (21 * u**4 - 14 * u**2 + 1) / 16), 1 / (26**0.5) *
          (13 * u * (u - 1) * (u + 1) * (33 * u**4 - 30 * u**2 + 5) / 16),
          1 / (30**0.5) * (15 * (u - 1) * (u + 1) *
                           (429 * u**6 - 495 * u**4 + 135 * u**2 - 5) / 128),
          1 / (34**0.5) * (17 * u * (u - 1) * (u + 1) *
                           (715 * u**6 - 1001 * u**4 + 385 * u**2 - 35) / 128)
      ])
    elif self.order == 10:
      lobatto = np.array([
          -1 / 2 * u + 1 / 2, 1 / 2 * u + 1 / 2,
          1 / (6**0.5) * (1.5 * (u - 1) * (u + 1)),
          1 / (10**0.5) * (5 * u * (u - 1) * (u + 1) / 2),
          1 / (14**0.5) * (7 * (u - 1) * (u + 1) * (5 * u**2 - 1) / 8),
          1 / (18**0.5) * (9 * u * (u - 1) * (u + 1) * (7 * u**2 - 3) / 8),
          1 / (22**0.5) * (11 * (u - 1) * (u + 1) *
                           (21 * u**4 - 14 * u**2 + 1) / 16), 1 / (26**0.5) *
          (13 * u * (u - 1) * (u + 1) * (33 * u**4 - 30 * u**2 + 5) / 16),
          1 / (30**0.5) * (15 * (u - 1) * (u + 1) *
                           (429 * u**6 - 495 * u**4 + 135 * u**2 - 5) / 128),
          1 / (34**0.5) * (17 * u * (u - 1) * (u + 1) *
                           (715 * u**6 - 1001 * u**4 + 385 * u**2 - 35) / 128),
          1 / (38**0.5) *
          (19 * (u - 1) * (u + 1) *
           (2431 * u**8 - 4004 * u**6 + 2002 * u**4 - 308 * u**2 + 7) / 256)
      ])
    elif self.order == 11:
      lobatto = np.array([
          -1 / 2 * u + 1 / 2,
          1 / 2 * u + 1 / 2,
          1 / (6**0.5) * (1.5 * (u - 1) * (u + 1)),
          1 / (10**0.5) * (5 * u * (u - 1) * (u + 1) / 2),
          1 / (14**0.5) * (7 * (u - 1) * (u + 1) * (5 * u**2 - 1) / 8),
          1 / (18**0.5) * (9 * u * (u - 1) * (u + 1) * (7 * u**2 - 3) / 8),
          1 / (22**0.5) * (11 * (u - 1) * (u + 1) *
                           (21 * u**4 - 14 * u**2 + 1) / 16),
          1 / (26**0.5) * (13 * u * (u - 1) * (u + 1) *
                           (33 * u**4 - 30 * u**2 + 5) / 16),
          1 / (30**0.5) * (15 * (u - 1) * (u + 1) *
                           (429 * u**6 - 495 * u**4 + 135 * u**2 - 5) / 128),
          1 / (34**0.5) * (17 * u * (u - 1) * (u + 1) *
                           (715 * u**6 - 1001 * u**4 + 385 * u**2 - 35) / 128),
          1 / (38**0.5) *
          (19 * (u - 1) * (u + 1) *
           (2431 * u**8 - 4004 * u**6 + 2002 * u**4 - 308 * u**2 + 7) / 256),
          1 / (42**0.5) *
          (21 * u * (u - 1) * (u + 1) *
           (4199 * u**8 - 7315 * u**6 + 4004 * u**4 - 715 * u**2 + 21) / 256),
      ])
    return lobatto

  def get_shape_functions(self, u):
    return self.polynomial(u)

  def derivative(self, u):
    if self.order == 1:
      d_lobatto = np.array([-1 / 2, 1 / 2])
    elif self.order == 2:
      d_lobatto = np.array([-1 / 2, 1 / 2, 1 / (6**0.5) * 3 * u])
    elif self.order == 3:
      d_lobatto = np.array([
          -1 / 2, 1 / 2, 1 / (6**0.5) * 3 * u,
          1 / (10**0.5) * (15 * u**2 / 2 - 5 / 2)
      ])
    elif self.order == 4:
      d_lobatto = np.array([
          -1 / 2, 1 / 2, 1 / (6**0.5) * 3 * u,
          1 / (10**0.5) * (15 * u**2 / 2 - 5 / 2),
          1 / (14**0.5) * (7 * u * (5 * u**2 - 3) / 2)
      ])
    elif self.order == 5:
      d_lobatto = np.array([
          -1 / 2, 1 / 2, 1 / (6**0.5) * 3 * u,
          1 / (10**0.5) * (15 * u**2 / 2 - 5 / 2),
          1 / (14**0.5) * (7 * u * (5 * u**2 - 3) / 2),
          1 / (18**0.5) * (9 * (35 * u**4 - 30 * u**2 + 3) / 8)
      ])
    elif self.order == 6:
      d_lobatto = np.array([
          -1 / 2, 1 / 2, 1 / (6**0.5) * 3 * u,
          1 / (10**0.5) * (15 * u**2 / 2 - 5 / 2),
          1 / (14**0.5) * (7 * u * (5 * u**2 - 3) / 2),
          1 / (18**0.5) * (9 * (35 * u**4 - 30 * u**2 + 3) / 8),
          1 / (22**0.5) * (11 * u * (63 * u**4 - 70 * u**2 + 15) / 8)
      ])
    elif self.order == 7:
      d_lobatto = np.array([
          -1 / 2, 1 / 2, 1 / (6**0.5) * 3 * u,
          1 / (10**0.5) * (15 * u**2 / 2 - 5 / 2),
          1 / (14**0.5) * (7 * u * (5 * u**2 - 3) / 2),
          1 / (18**0.5) * (9 * (35 * u**4 - 30 * u**2 + 3) / 8),
          1 / (22**0.5) * (11 * u * (63 * u**4 - 70 * u**2 + 15) / 8), 1 /
          (26**0.5) * (13 * (231 * u**6 - 315 * u**4 + 105 * u**2 - 5) / 16)
      ])
    elif self.order == 8:
      d_lobatto = np.array([
          -1 / 2, 1 / 2, 1 / (6**0.5) * 3 * u,
          1 / (10**0.5) * (15 * u**2 / 2 - 5 / 2),
          1 / (14**0.5) * (7 * u * (5 * u**2 - 3) / 2),
          1 / (18**0.5) * (9 * (35 * u**4 - 30 * u**2 + 3) / 8),
          1 / (22**0.5) * (11 * u * (63 * u**4 - 70 * u**2 + 15) / 8), 1 /
          (26**0.5) * (13 * (231 * u**6 - 315 * u**4 + 105 * u**2 - 5) / 16),
          1 / (30**0.5) * (15 * u *
                           (429 * u**6 - 693 * u**4 + 315 * u**2 - 35) / 16)
      ])
    elif self.order == 9:
      d_lobatto = np.array([
          -1 / 2, 1 / 2, 1 / (6**0.5) * 3 * u,
          1 / (10**0.5) * (15 * u**2 / 2 - 5 / 2),
          1 / (14**0.5) * (7 * u * (5 * u**2 - 3) / 2),
          1 / (18**0.5) * (9 * (35 * u**4 - 30 * u**2 + 3) / 8),
          1 / (22**0.5) * (11 * u * (63 * u**4 - 70 * u**2 + 15) / 8), 1 /
          (26**0.5) * (13 * (231 * u**6 - 315 * u**4 + 105 * u**2 - 5) / 16),
          1 / (30**0.5) * (15 * u *
                           (429 * u**6 - 693 * u**4 + 315 * u**2 - 35) / 16),
          1 / (34**0.5) *
          (17 *
           (6435 * u**8 - 12012 * u**6 + 6930 * u**4 - 1260 * u**2 + 35) / 128)
      ])
    elif self.order == 10:
      d_lobatto = np.array([
          -1 / 2, 1 / 2, 1 / (6**0.5) * 3 * u,
          1 / (10**0.5) * (15 * u**2 / 2 - 5 / 2),
          1 / (14**0.5) * (7 * u * (5 * u**2 - 3) / 2),
          1 / (18**0.5) * (9 * (35 * u**4 - 30 * u**2 + 3) / 8),
          1 / (22**0.5) * (11 * u * (63 * u**4 - 70 * u**2 + 15) / 8), 1 /
          (26**0.5) * (13 * (231 * u**6 - 315 * u**4 + 105 * u**2 - 5) / 16),
          1 / (30**0.5) * (15 * u *
                           (429 * u**6 - 693 * u**4 + 315 * u**2 - 35) / 16),
          1 / (34**0.5) *
          (17 * (6435 * u**8 - 12012 * u**6 + 6930 * u**4 - 1260 * u**2 + 35) /
           128), 1 / (38**0.5) *
          (19 * u *
           (12155 * u**8 - 25740 * u**6 + 18018 * u**4 - 4620 * u**2 + 315) /
           128)
      ])
    elif self.order == 11:
      d_lobatto = np.array([
          -1 / 2, 1 / 2, 1 / (6**0.5) * 3 * u,
          1 / (10**0.5) * (15 * u**2 / 2 - 5 / 2),
          1 / (14**0.5) * (7 * u * (5 * u**2 - 3) / 2),
          1 / (18**0.5) * (9 * (35 * u**4 - 30 * u**2 + 3) / 8),
          1 / (22**0.5) * (11 * u * (63 * u**4 - 70 * u**2 + 15) / 8), 1 /
          (26**0.5) * (13 * (231 * u**6 - 315 * u**4 + 105 * u**2 - 5) / 16),
          1 / (30**0.5) * (15 * u *
                           (429 * u**6 - 693 * u**4 + 315 * u**2 - 35) / 16),
          1 / (34**0.5) *
          (17 * (6435 * u**8 - 12012 * u**6 + 6930 * u**4 - 1260 * u**2 + 35) /
           128), 1 / (38**0.5) *
          (19 * u *
           (12155 * u**8 - 25740 * u**6 + 18018 * u**4 - 4620 * u**2 + 315) /
           128),
          1 / (42**0.5) * (46189 * u**12 - 88179 * u**10 + 48450 * u**8 -
                           8280 * u**6 + 462 * u**4 - 7 * u**2)
      ])
    return d_lobatto

  def get_der_shape_functions(self, u):
    return self.derivative(u)


class Larange(BasePolynomial):

  def polynomial(self):
    larange = np.zeros((5))
    larange[0] = lambda x: (1 - x) / 2
    larange[1] = lambda x: (1 + x) / 2

    larange[2] = lambda x: (x - 0) * (x - 1) / (-1 - 0) * (-1 - 1)
    larange[3] = lambda x: (x + 1) * (x - 1) / (0 + 1) * (0 - 1)
    larange[4] = lambda x: (x + 1) * (x - 0) / (1 + 1) * (1 - 0)

    return larange

  def get_shape_functions(self):
    N_larange = []
    if self.order == 1:
      N_larange.append(self.polynomial()[0])
      N_larange.append(self.polynomial()[1])
    elif self.order == 2:
      N_larange.append(self.polynomial()[2])
      N_larange.append(self.polynomial()[3])
      N_larange.append(self.polynomial()[4])
    else:
      print("cubic larange not supported yet")

  def derivative(self):
    d_larange = np.zeros((5))
    d_larange[0] = lambda x: -1 / 2
    d_larange[1] = lambda x: 1 / 2

    d_larange[2] = lambda x: (2 * x - 1) / (-1 - 0) * (-1 - 1)
    d_larange[3] = lambda x: (2 * x) / (0 + 1) * (0 - 1)
    d_larange[4] = lambda x: (2 * x + 1) / (1 + 1) * (1 - 0)

    return

  def get_der_shape_functions(self):
    B_larange = np.zeros((self.order + 1))
    if self.order == 1:
      B_larange[0] = self.derivative()[0]
      B_larange[1] = self.derivative()[1]
    elif self.order == 2:
      B_larange[0] = self.derivative()[2]
      B_larange[1] = self.derivative()[3]
      B_larange[2] = self.derivative()[4]
    else:
      print("cubic larange not supported yet")

    return


# 2D triangular larange shape functions in reference coordiantates
class Lagrange2DTri:

  def __init__(self, order):
    self.order = order

  def polynomial(self, u, v):
    if self.order == 1:
      lagrange = np.array([
          1 - u - v,
          u,
          v,
      ])
    elif self.order == 2:
      lagrange = np.array([
          1 - 3 * u - 3 * v + 2 * u**2 + 4 * u * v + 2 * v**2,    # at (0,0)
          2 * u**2 - u,    # at (1,0)
          2 * v**2 - v,    # at (0,1)
          4 * u * (1 - u - v),    # between (0,1)
          4 * u * v,    # between (1,2)
          4 * v * (1 - u - v)    # between (0,2)
      ])
    else:
      print("cubic larange not supported yet")

    return lagrange

  def get_shape_functions(self, u, v):
    return self.polynomial(u, v)

  def derivative(self, u, v):
    if self.order == 1:
      d_lagrange = np.array([[-1, -1], [1, 0], [0, 1]])
    elif self.order == 2:
      d_lagrange = np.array([[-3 + 4 * u + 4 * v, -3 + 4 * u + 4 * v],
                             [4 * u - 1, 0], [0, 4 * v - 1],
                             [4 - 8 * u - 4 * v, -4 * u], [4 * v, 4 * u],
                             [-4 * v, 4 - 4 * u - 8 * v]])
    else:
      print("cubic larange not supported yet")

    return d_lagrange

  def get_der_shape_functions(self, u, v):
    return self.derivative(u, v)

  def jacobi(self):
    jacobi = np.array([[1, 0], [0, 1]])
    return jacobi

  def inverse_jacobi(self):
    inverse_jacobi = np.array([[1, 0], [0, 1]])
    return inverse_jacobi

  @property
  def determinant_jacobi(self):
    return 1


class Lagrange2DQuad:

  def __init__(self, order):
    self.order = order

  def polynomial(self, u, v):
    if self.order == 1:
      lagrange = np.array([
          1 / 4 * (1 - u) * (1 - v), 1 / 4 * (1 + u) * (1 - v),
          1 / 4 * (1 + u) * (1 + v), 1 / 4 * (1 - u) * (1 + v)
      ])
    elif self.order == 2:
      lagrange = np.array([
          1 / 4 * (u - 1) * (v - 1) * (u + v + 1),
          1 / 4 * (u + 1) * (v - 1) * (u - v - 1),
          1 / 4 * (u + 1) * (v + 1) * (u + v - 1),
          1 / 4 * (u - 1) * (v + 1) * (u - v + 1),
          1 / 2 * (1 - u**2) * (v - 1), 1 / 2 * (u + 1) * (1 - v**2),
          1 / 2 * (1 - u**2) * (v + 1), 1 / 2 * (u - 1) * (1 - v**2),
          1 / 2 * (1 - v**2) * (1 - u), 1 / 2 * (1 - u**2) * (1 + v),
          1 / 2 * (1 - v**2) * (1 + u), 1 / 2 * (1 + u) * (1 - v**2)
      ])
    else:
      print("cubic larange not supported yet")

    return lagrange

  def get_shape_functions(self, u, v):
    return self.polynomial(u, v)

  def derivative(self, u, v):
    if self.order == 1:
      d_lagrange = np.array([[-1 / 4 * (1 - v), -1 / 4 * (1 - u)],
                             [1 / 4 * (1 - v), -1 / 4 * (1 + u)],
                             [1 / 4 * (1 + v), 1 / 4 * (1 + u)],
                             [-1 / 4 * (1 + v), 1 / 4 * (1 - u)]])
    elif self.order == 2:
      d_lagrange = np.array([[1 / 4 * (2 * u + v), 1 / 4 * (u + 2 * v)],
                             [1 / 4 * (2 * u - v), 1 / 4 * (-u + 2 * v)],
                             [1 / 4 * (2 * u + v), 1 / 4 * (u + 2 * v)],
                             [1 / 4 * (2 * u - v), 1 / 4 * (-u + 2 * v)],
                             [-u * (v - 1), 1 / 2 * (u**2 - 1)],
                             [1 / 2 * (v**2 - 1), -v * (u + 1)],
                             [-u * (v + 1), 1 / 2 * (u**2 - 1)],
                             [1 / 2 * (v**2 - 1), -v * (u - 1)],
                             [-v * (1 - u), -1 / 2 * (v**2 - 1)],
                             [-u * (1 + v), -1 / 2 * (u**2 - 1)],
                             [-v * (1 + u), -1 / 2 * (v**2 - 1)],
                             [-1 / 2 * (1 - v**2), -u * (1 + v)]])
    else:
      print("cubic larange not supported yet")

    return d_lagrange

  def get_der_shape_functions(self, u, v):
    return self.derivative(u, v)


class Lagrange3DTetra:

  def __init__(self, order):
    self.order = order

  def polynomial(self, u, v, w):
    if self.order == 1:
      lagrange = np.array([1 - u - v - w, u, v, w])
    elif self.order == 2:
      lagrange = np.array([
          1 - 3 * u - 3 * v - 3 * w + 2 * u**2 + 4 * u * v + 4 * u * w +
          2 * v**2 + 4 * v * w + 2 * w**2,    # at (0,0,0)
          2 * u**2 - u,    # at (1,0,0)
          2 * v**2 - v,    # at (0,1,0)
          2 * w**2 - w,    # at (0,0,1)
          4 * u * (1 - u - v - w),    # between (0,1) 
          4 * u * v,    # between (1,2)
          4 * u * w,    # between (1,2) and (1,3)
          4 * v * (1 - u - v - w),    # between (0,2)
          4 * v * w,    # between (2,3)
          4 * w * (1 - u - v - w)    # between (0,3)
      ])
    else:
      print("cubic larange not supported yet")

    return lagrange

  def get_shape_functions(self, u, v, w):
    return self.polynomial(u, v, w)

  def derivative(self, u, v, w):
    if self.order == 1:
      d_lagrange = np.array([[-1, -1, -1], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    elif self.order == 2:
      d_lagrange = np.array([
          [
              -3 + 4 * u + 4 * v + 4 * w, -3 + 4 * u + 4 * v + 4 * w,
              -3 + 4 * u + 4 * v + 4 * w
          ],    # at (0,0,0)
          [4 * u - 1, 0, 0],    # at (1,0,0)
          [0, 4 * v - 1, 0],    # at (0,1,0)
          [0, 0, 4 * w - 1],    # at (0,0,1)
          [4 - 8 * v - 4 * w - 8 * u, -4 * u, -4 * u],    # between (0,1)
          [4 * v, 4 * u, 0],    # between (1,2)
          [0, 4 * w, 4 * u],    # between (1,2) and (1,3)
          [-3 + 4 * w + 4 * u + 4 * v, -3 + 4 * u, -3 + 4 * u],
          [0, 4 * w, 4 * v],
          [4 * w, 0, 4 * v]
      ])
    else:
      print("cubic larange not supported yet")

    return d_lagrange

  def get_der_shape_functions(self, u, v, w):
    return self.derivative(u, v, w)


class PolyBuilder:
  """build polynomial class"""

  def __init__(self, order):
    self.order = order

  def build(self, poly_type):
    if poly_type == 'lobatto':
      return Lobatto(self.order, )
    elif poly_type == 'larange':
      return Larange(self.order)
    else:
      raise ValueError('poly_type must be lobatto or larange')

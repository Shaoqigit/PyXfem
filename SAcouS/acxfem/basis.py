# This file is part of PyXfem, a software distributed under the MIT license.
# For any question, please contact the authors cited below.
#
# Copyright (c) 2023
# 	Shaoqi WU <shaoqiwu@outlook.com
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

# basis.py mutiple element types
# Lobatto element are recommended to use

import numpy as np
from abc import ABCMeta, abstractmethod
from functools import cached_property

from SAcouS.acxfem.precompute_matrices import Ke1D, Me1D, Ce1D
from SAcouS.acxfem.precompute_matrices_lag import Ke2DTri, Me2DTri


class Base1DElement(metaclass=ABCMeta):
  """base abstract FE elementary matrix class
    precompute the shape function and its derivative on gauss points
    parameters:
    order: int
        element order
    nodes: ndarray
        1d: [x1, x2]
        2d: [(x1, y1), (x2, y2)]
        3d: [(x1, y1, z1), (x2, y2, z2)]
    """

  def __init__(self, label, order, nodes):
    self.label = label
    self.order = order
    self.nodes = nodes
    self.is_discontinue = False

  @cached_property
  def Jacobian(self):
    """
    compute the Jacobian of the element
    returns:
    J: 
    J=dx/dxi"""

    return np.abs(self.nodes[0] - self.nodes[1]) / 2

  @cached_property
  def inverse_Jacobian(self):
    return 2 / np.abs(self.nodes[0] - self.nodes[1])


class Lobbato1DElement(Base1DElement):
  """FE lobatto 1D basis class
    parameters:
    order: int
        element order

    returns:
    """

  def __init__(self, label, order, nodes):
    super().__init__(label, order, nodes)

  @cached_property
  def ke(self):
    """compute the elementary stiffness matrix
    returns:
    K: ndarray
        elementary stiffness matrix
    """
    Ke = 0
    if self.order == 1:
      Ke = self.inverse_Jacobian * Ke1D[0]
    elif self.order == 2:
      Ke = self.inverse_Jacobian * Ke1D[1]
    elif self.order == 3:
      Ke = self.inverse_Jacobian * Ke1D[2]
    elif self.order == 4:
      Ke = self.inverse_Jacobian * Ke1D[3]
    else:
      print("quadrtic lobatto not supported yet")
    return Ke

  @cached_property
  def me(self):
    """compute the elementary stiffness matrix
    returns:
    m: ndarray
        elementary stiffness matrix
    """
    if self.order == 1:
      Me = self.Jacobian * Me1D[0]
    elif self.order == 2:
      Me = self.Jacobian * Me1D[1]
    elif self.order == 3:
      Me = self.Jacobian * Me1D[2]
    elif self.order == 4:
      Me = self.Jacobian * Me1D[3]
    else:
      print("quadrtic lobatto not supported yet")
    return Me

  @cached_property
  def ce(self):
    """compute the elementary coupling matrix, N(x)B(x)
    returns:
    c: ndarray
        elementary coupling matrix
    """
    if self.order == 1:
      Ce = Ce1D[0]
    elif self.order == 2:
      Ce = Ce1D[1]
    elif self.order == 3:
      Ce = Ce1D[2]
    elif self.order == 4:
      Ce = Ce1D[3]
    else:
      print("quadrtic lobatto not supported yet")
    return Ce

  def get_order(self):
    return self.order

  @property
  def nb_internal_dofs(self):
    return self.order - 1

  @property
  def local_dofs_index(self):
    return np.arange(self.order + 1)


class Base2DElement(metaclass=ABCMeta):
  """base abstract FE elementary matrix class
    precompute the shape function and its derivative on gauss points
    parameters:
    order: int
        element order
    nodes: ndarray
        2d: [(x1, y1), (x2, y2), (x3, y3)] for triangle
            [(x1, y1), (x2, y2), (x3, y3), (x4, y4)] for quad
    """

  def __init__(self, label, order, vertices):
    self.label = label
    self.order = order
    self.vertices = vertices
    self.is_discontinue = False

  @cached_property
  @abstractmethod
  def Jacobian(self):
    """
    compute the Jacobian of the element
    returns:
    J: 
    J=dx/dxi"""

    pass

  @cached_property
  @abstractmethod
  def inverse_Jacobian(self):
    pass


from SAcouS.acxfem.quadratures import GaussLegendre2DTri
from SAcouS.acxfem.polynomial import Lagrange2DTri


class Lagrange2DTriElement(Base2DElement):
  """FE lagrange 2D triangle basis class
    parameters:
    order: int
        element order
    vertices: ndarray
        [(x1, y1), (x2, y2), (x3, y3)] for triangle
    reference element: [(0, 0), (1, 0), (0, 1)]
    illustrated as below:
    2
    |\
    | \
    |  \
    |   \
    |    \
    0-----1
    """

  def __init__(self, label, order, vertices):
    super().__init__(label, order, vertices)
    n_quad_pts = order * 2 + 1
    lag_poly = Lagrange2DTri(order)

    self.points, self.weights = GaussLegendre2DTri(
        n_quad_pts).points(), GaussLegendre2DTri(n_quad_pts).weights()
    self.N = np.array(
        [lag_poly.get_shape_functions(*point) for point in self.points])
    self.B = np.array(
        [lag_poly.get_der_shape_functions(*point) for point in self.points])

    self.J = self.Jacobian
    self.inv_J = self.inverse_Jacobian
    self.det_J = self.determinant_Jacobian

  @cached_property
  def Jacobian(self):
    """
    compute the Jacobian of the element
    returns:
    J: 
    J=dx/dxi"""
    J = np.array([[
        self.vertices[1][0] - self.vertices[0][0],
        self.vertices[1][1] - self.vertices[0][1]
    ],
                  [
                      self.vertices[2][0] - self.vertices[0][0],
                      self.vertices[2][1] - self.vertices[0][1]
                  ]])

    return J

  @cached_property
  def inverse_Jacobian(self):
    return np.linalg.inv(self.J)

  @cached_property
  def determinant_Jacobian(self):
    return np.linalg.det(self.Jacobian)

  @cached_property
  def ke(self):
    """compute the elementary stiffness matrix
    returns:
    K: ndarray
        elementary stiffness matrix
    """
    if self.order == 1:
      Ke = sum(
          self.B[i, :, :] @ self.inv_J.T @ self.inv_J @ self.B[i, :, :].T *
          weight for i, weight in enumerate(self.weights)) * self.det_J

    else:
      print("quadrtic lagrange not implemented yet")

    return Ke

  @cached_property
  def me(self):
    """compute the elementary stiffness matrix
    returns:
    m: ndarray
        elementary stiffness matrix
    """
    if self.order == 1:
      #   import pdb
      #   pdb.set_trace()
      weight = np.diag(
          np.array([self.weights[0], self.weights[1], self.weights[2]]))
      Me = self.N[:, :].T @ weight @ self.N[:, :] * self.det_J
    else:
      print("quadrtic lagrange not implemented yet")
    return Me

  def get_order(self):
    return self.order

  @property
  def nb_internal_dofs(self):
    if order == 1 or order == 2:
      return 0
    elif order == 3:
      self.nb_internal_dofs = 1

  @property
  def nb_edge_dofs(self):
    if order == 1:
      return 0
    elif order == 2:
      return 3
    elif order == 3:
      return 6

  @property
  def local_dofs_index(self):
    return np.arange(self.order + 1)


if __name__ == "__main__":
  label = "fluid"
  order = 1
  nodes = np.array([[0, 0], [2, -1], [1, 0.5]])
  # plot the nodes
  import matplotlib.pyplot as plt
  plt.figure()
  plt.plot(nodes[:, 0], nodes[:, 1], 'o')
  plt.show()
  lag_2d_tri = Lagrange2DTriElement(1, order, nodes)
  print(lag_2d_tri.ke)
  print(lag_2d_tri.me)

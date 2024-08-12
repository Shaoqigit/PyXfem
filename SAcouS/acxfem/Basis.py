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

from .PrecomputeMatrices import Ke1D, Me1D, Ce1D
from .Quadratures import GaussLegendre2DTri


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
    if len(self.nodes.shape) == 1:
      return np.abs(self.nodes[0] - self.nodes[1]) / 2
    elif len(self.nodes.shape) == 2:
      return np.sqrt((self.nodes[0, 0] - self.nodes[1, 0])**2 +
                     (self.nodes[0, 1] - self.nodes[1, 1])**2) / 2

  @cached_property
  def inverse_Jacobian(self):
    if len(self.nodes.shape) == 1:
      return 2 / np.abs(self.nodes[0] - self.nodes[1])
    elif len(self.nodes.shape) == 2:
      return 2 / np.sqrt((self.nodes[0, 0] - self.nodes[1, 0])**2 +
                         (self.nodes[0, 1] - self.nodes[1, 1])**2)


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


class Helmholtz1DElement(Lobbato1DElement):

  def __init__(self, label, order, nodes, mat_coeffs=[]):
    super().__init__(label, order, nodes)
    self.mat_coeffs = mat_coeffs

  @cached_property
  def ke(self):
    """compute the elementary stiffness matrix
    returns:
    K: ndarray
        elementary stiffness matrix
    """
    Ke = 0
    if self.order == 1:
      Ke = self.inverse_Jacobian * self.mat_coeffs[0] * Ke1D[0]
    elif self.order == 2:
      Ke = self.inverse_Jacobian * self.mat_coeffs[0] * Ke1D[1]
    elif self.order == 3:
      Ke = self.inverse_Jacobian * self.mat_coeffs[0] * Ke1D[2]
    elif self.order == 4:
      Ke = self.inverse_Jacobian * self.mat_coeffs[0] * Ke1D[3]
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
      Me = self.Jacobian * self.mat_coeffs[1] * Me1D[0]
    elif self.order == 2:
      Me = self.Jacobian * self.mat_coeffs[1] * Me1D[1]
    elif self.order == 3:
      Me = self.Jacobian * self.mat_coeffs[1] * Me1D[2]
    elif self.order == 4:
      Me = self.Jacobian * self.mat_coeffs[1] * Me1D[3]
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
      Ce = self.mat_coeffs[2] * Ce1D[0]
    elif self.order == 2:
      Ce = self.mat_coeffs[2] * Ce1D[1]
    elif self.order == 3:
      Ce = self.mat_coeffs[2] * Ce1D[2]
    elif self.order == 4:
      Ce = self.mat_coeffs[2] * Ce1D[3]
    else:
      print("quadrtic lobatto not supported yet")
    return Ce


from SAcouS.acxfem.PrecomputeMatricesLag import points_o1, weights_o1
from SAcouS.acxfem.PrecomputeMatricesLag import get_N_B_p1


class BaseNDElement(metaclass=ABCMeta):
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

    self.Jacobian()
    self.determinant_Jacobian()
    self.inverse_Jacobian()
    self.inv_J_product = self.inv_J.T @ self.inv_J

  @abstractmethod
  def Jacobian(self):
    """
    compute the Jacobian of the element
    returns:
    J: 
    J=dx/dxi"""

    pass

  @abstractmethod
  def inverse_Jacobian(self):
    pass

  @abstractmethod
  def determinant_Jacobian(self):
    pass


class Lagrange2DTriElement(BaseNDElement):
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
  points = points_o1,
  weights = weights_o1

  def __init__(self, label, order, vertices):
    super().__init__(label, order, vertices)
    if order == 1:
      self.N, self.B = get_N_B_p1(2)
    else:
      print("quadrtic lagrange not supported yet")

  def Jacobian(self):
    """
    compute the Jacobian of the element
    returns:
    J: 
    J=dx/dxi"""
    self.J = np.array([[
        self.vertices[1][0] - self.vertices[0][0],
        self.vertices[1][1] - self.vertices[0][1]
    ],
                       [
                           self.vertices[2][0] - self.vertices[0][0],
                           self.vertices[2][1] - self.vertices[0][1]
                       ]])

  def inverse_Jacobian(self):
    self.inv_J = np.array([[self.J[1, 1], -self.J[0, 1]],
                           [-self.J[1, 0], self.J[0, 0]]]) * 1 / self.det_J

  def determinant_Jacobian(self):
    self.det_J = self.J[0, 0] * self.J[1, 1] - self.J[0, 1] * self.J[1, 0]

  @cached_property
  def ke(self):
    """compute the elementary stiffness matrix
    returns:
    K: ndarray
        elementary stiffness matrix
    """
    if self.order == 1:
      Ke = sum(
          self.B[i, :, :] @ self.inv_J_product @ self.B[i, :, :].T * weight
          for i, weight in enumerate(self.weights)) * self.det_J

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
      weight = np.diag(
          np.array([self.weights[0], self.weights[1], self.weights[2]]))
      Me = self.N[:, :].T @ weight @ self.N[:, :] * self.det_J
    else:
      print("quadrtic lagrange not implemented yet")
    return Me

  def weights_and_points(self, integ_order=1):
    if integ_order is None:
      integ_order = 2 * self.order + 1
    points, weights = GaussLegendre2DTri(
        integ_order).points(), GaussLegendre2DTri(3).weights()
    return points, weights

  def egde_basis(self, edge):
    """compute the edge basis function
    returns:
    N: ndarray
        edge basis function
    """
    N = np.zeros((self.order + 1, 2))
    if self.order == 1:
      if edge == 0:
        N[0, 0] = 1
      elif edge == 1:
        N[1, 0] = 1
      elif edge == 2:
        N[2, 0] = 1
    else:
      print("quadrtic lagrange not implemented yet")
    return N

  def edge_jacobian(self, edge):
    J = np.array([[
        self.vertices[edge[1]][0] - self.vertices[edge[0]][0],
        self.vertices[edge[1]][1] - self.vertices[edge[0]][1]
    ]])

  def get_order(self):
    return self.order

  @cached_property
  def nb_internal_dofs(self):
    if self.order == 1 or self.order == 2:
      return 0
    elif self.order == 3:
      return 1

  @cached_property
  def nb_edge_dofs(self):
    if self.order == 1:
      return 0
    elif self.order == 2:
      return 3
    elif self.order == 3:
      return 6

  @cached_property
  def local_dofs_index(self):
    return np.arange(self.order * 2 + 1)


# consider material properties: ease the matrix assembly
class Helmholtz2DElement(Lagrange2DTriElement):

  def __init__(self, label, order, vertices, mat_coeffs=[]):
    super().__init__(label, order, vertices)
    self.mat_coeffs = mat_coeffs

  @cached_property
  def ke(self):
    """compute the elementary stiffness matrix
    returns:
    K: ndarray
        elementary stiffness matrix
    """
    if self.order == 1:
      Ke = self.mat_coeffs[0] * sum(
          self.B[i, :, :] @ self.inv_J_product @ self.B[i, :, :].T * weight
          for i, weight in enumerate(self.weights)) * self.det_J

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
      weight = np.diag(
          np.array([self.weights[0], self.weights[1], self.weights[2]]))
      Me = self.mat_coeffs[
          1] * self.N[:, :].T @ weight @ self.N[:, :] * self.det_J
    else:
      print("quadrtic lagrange not implemented yet")
    return Me


class Lagrange2DQuadElement(BaseNDElement):
  """FE lagrange 2D quad basis class
    parameters:
    order: int
        element order
    vertices: ndarray
        [(x1, y1), (x2, y2), (x3, y3), (x4, y4)] for quad
    reference element: [(0, 0), (1, 0), (0, 1), (1, 1)]
    illustrated as below:
    3----2
    |    |
    |    |
    0----1
    """

  def __init__(self, label, order, vertices):
    super().__init__(label, order, vertices)

  def Jacobian(self):
    """
    compute the Jacobian of the element
    returns:
    J: 
    J=dx/dxi"""
    self.J = np.array([[
        self.vertices[1][0] - self.vertices[0][0],
        self.vertices[1][1] - self.vertices[0][1]
    ],
                       [
                           self.vertices[2][0] - self.vertices[0][0],
                           self.vertices[2][1] - self.vertices[0][1]
                       ]])

  def inverse_Jacobian(self):
    self.inv_J = np.array([[self.J[1, 1], -self.J[0, 1]],
                           [-self.J[1, 0], self.J[0, 0]]]) * 1 / self.det_J

  def determinant_Jacobian(self):
    self.det_J = self.J[0, 0] * self.J[1, 1] - self.J[0, 1] * self.J[1, 0]

  @cached_property
  def ke(self):
    """compute the elementary stiffness matrix
    returns:
    K: ndarray
        elementary stiffness matrix
    """
    if self.order == 1:
      Ke = sum(
          self.B[i, :, :] @ self.inv_J_product @ self.B[i, :, :].T * weight
          for i, weight in enumerate(self.weights)) * self.det_J

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
      weight = np.diag(
          np.array([self.weights[0], self.weights[1], self.weights[2]]))
      Me = self.N[:, :].T @ weight @ self.N[:, :] * self.det_J
    else:
      print("quadrtic lagrange not implemented yet")
    return Me


from SAcouS.acxfem.PrecomputeMatricesLag import points_tetra_o1, weights_tetra_o1


class Lagrange3DTetraElement(BaseNDElement):
  """FE lagrange 3D tetra basis class
    parameters:
    order: int
        element order
    vertices: ndarray
        [(x1, y1, z1), (x2, y2, z2), (x3, y3, z3), (x4, y4, z4)] for tetra
    reference element: [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
    illustrated as below:
         3
        /|\
       / | \
      /  |  \
     2---|---1
      \  |  /
       \ | /
        \|/
         0
    
    """

  points = points_tetra_o1,
  weights = weights_tetra_o1

  def __init__(self, label, order, vertices):
    super().__init__(label, order, vertices)
    if order == 1:
      self.N, self.B = get_N_B_p1(3)
    else:
      print("quadrtic lagrange not supported yet")

  def Jacobian(self):
    """
    compute the Jacobian of the element
    returns:
    J: 
    J=dx/dxi"""
    self.J = np.array([[
        self.vertices[1][0] - self.vertices[0][0],
        self.vertices[1][1] - self.vertices[0][1],
        self.vertices[1][2] - self.vertices[0][2]
    ],
                       [
                           self.vertices[2][0] - self.vertices[0][0],
                           self.vertices[2][1] - self.vertices[0][1],
                           self.vertices[2][2] - self.vertices[0][2]
                       ],
                       [
                           self.vertices[3][0] - self.vertices[0][0],
                           self.vertices[3][1] - self.vertices[0][1],
                           self.vertices[3][2] - self.vertices[0][2]
                       ]])

  def inverse_Jacobian(self):
    self.inv_J = np.linalg.inv(self.J)

  def determinant_Jacobian(self):
    self.det_J = np.linalg.det(self.J)

  @cached_property
  def ke(self):
    """compute the elementary stiffness matrix
    returns:
    K: ndarray
        elementary stiffness matrix
    """
    if self.order == 1:
      Ke = sum(
          self.B[i, :, :] @ self.inv_J_product @ self.B[i, :, :].T * weight
          for i, weight in enumerate(self.weights)) * self.det_J

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
      weight = np.diag(
          np.array([
              self.weights[0], self.weights[1], self.weights[2],
              self.weights[3]
          ]))
      Me = self.N[:, :].T @ weight @ self.N[:, :] * self.det_J
    else:
      print("quadrtic lagrange not implemented yet")
    return Me


if __name__ == "__main__":
  label = "fluid"
  order = 1
  # nodes = np.array([[0, 0], [2, -1], [1, 0.5]])
  nodes = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
  # plot the nodes
  import matplotlib.pyplot as plt
  # plt.figure()
  # plt.plot(nodes[:, 0], nodes[:, 1], 'o')
  # plot the nodes in 3d
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  ax.scatter(nodes[:, 0], nodes[:, 1], nodes[:, 2])
  plt.show()
  lag_3d_tetra = Lagrange3DTetraElement(1, order, nodes)
  print(lag_3d_tetra.ke)
  # lag_2d_tri = Lagrange2DTriElement(1, order, nodes)
  # print(lag_2d_tri.ke)
  # print(lag_2d_tri.me)

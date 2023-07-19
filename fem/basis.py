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

from fem.precompute_matrices import Ke1D, Me1D

            

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
    def __init__(self, order, nodes):
        self.order = order
        self.nodes = nodes
        self.is_discontinue = False

    @property
    def Jacobian(self):
        """
        compute the Jacobian of the element
        returns:
        J: 
        J=dx/dxi"""

        return np.abs(self.nodes[0]-self.nodes[1])/2

    @property
    def inverse_Jacobian(self):
        return 2/np.abs(self.nodes[0]-self.nodes[1])


class Lobbato1DElement(Base1DElement):
    """FE lobatto 1D basis class
    parameters:
    order: int
        element order

    returns:
    """
    def __init__(self, order, nodes):
        super().__init__(order, nodes)

    @property
    def ke(self):
        """compute the elementary stiffness matrix
        returns:
        K: ndarray
            elementary stiffness matrix
        """
        Ke = 0
        if self.order == 1:
            Ke =  self.inverse_Jacobian*Ke1D[0]
        elif self.order == 2:
            Ke =  self.inverse_Jacobian*Ke1D[1]
        elif self.order == 3:
            Ke =  self.inverse_Jacobian*Ke1D[2]
        elif self.order == 4:
            Ke =  self.inverse_Jacobian*Ke1D[3]
        else:
            print("quadrtic lobatto not supported yet")
        return Ke
    
    @property
    def me(self):
        """compute the elementary stiffness matrix
        returns:
        K: ndarray
            elementary stiffness matrix
        """
        if self.order == 1:
            Me =  self.Jacobian*Me1D[0]
        elif self.order == 2:
            Me =  self.Jacobian*Me1D[1]
        elif self.order == 3:
            Me =  self.Jacobian*Me1D[2]
        elif self.order == 4:
            Me =  self.Jacobian*Me1D[3]
        else:
            print("quadrtic lobatto not supported yet")
        return Me
    
    def get_order(self):
        return self.order
    
    def num_internal_dofs(self):
        return self.order-1
    
    def local_dofs_index(self):
        return np.arange(self.order+1)
    

# l_element = Lobbato1DElement(2, [0, 0.5])
# k, m = l_element.get_matrix()
# print(k)
# print(m)
# print(l_element.get_order())
# print(l_element.num_internal_dofs())
# print(l_element.local_dofs_index())
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

    def polynomial(self):
        lobatto = [
        lambda x: -1/2 * x + 1/2,
        lambda x: 1/2 * x + 1/2,
        lambda x: 1/(6**0.5) * (1.5*(x - 1)*(x + 1)),
        lambda x: 1/(10**0.5) * (5*x*(x - 1)*(x + 1)/2),
        lambda x: 1/(14**0.5) * (7*(x - 1)*(x + 1)*(5*x**2 - 1)/8),
        lambda x: 1/(18**0.5) * (9*x*(x - 1)*(x + 1)*(7*x**2 - 3)/8),
        lambda x: 1/(22**0.5) * (11*(x - 1)*(x + 1)*(21*x**4 - 14*x**2 + 1)/16),
        lambda x: 1/(26**0.5) * (13*x*(x - 1)*(x + 1)*(33*x**4 - 30*x**2 + 5)/16),
        lambda x: 1/(30**0.5) * (15*(x - 1)*(x + 1)*(429*x**6 - 495*x**4 + 135*x**2 - 5)/128),
        lambda x: 1/(34**0.5) * (17*x*(x - 1)*(x + 1)*(715*x**6 - 1001*x**4 + 385*x**2 - 35)/128),
        lambda x: 1/(38**0.5) * (19*(x - 1)*(x + 1)*(2431*x**8 - 4004*x**6 + 2002*x**4 - 308*x**2 + 7)/256),
        lambda x: 1/(42**0.5) * (21*x*(x - 1)*(x + 1)*(4199*x**8 - 7315*x**6 + 4004*x**4 - 715*x**2 + 21)/256),
        ]
        return lobatto

    def get_shape_functions(self):
        N_lobatto = []
        for i in range(self.order+1):
            N_lobatto.append(self.polynomial()[i])
            
        return N_lobatto
    
    def derivative(self):
        d_lobatto = [
        lambda x: -1/2,
        lambda x: 1/2,
        lambda x: 1/(6**0.5) * (3*x),
        lambda x: 1/(10**0.5) * (15*x**2/2 - 5/2),
        lambda x: 1/(14**0.5) * (7*x*(5*x**2 - 3)/2),
        lambda x: 1/(18**0.5) * (9*(35*x**4 - 30*x**2 + 3)/8),
        lambda x: 1/(22**0.5) * (11*x*(63*x**4 - 70*x**2 + 15)/8),
        lambda x: 1/(26**0.5) * (13*(231*x**6 - 315*x**4 + 105*x**2 - 5)/16),
        lambda x: 1/(30**0.5) * (15*x*(429*x**6 - 693*x**4 + 315*x**2 - 35)/16),
        lambda x: 1/(34**0.5) * (17*(6435*x**8 - 12012*x**6 + 6930*x**4 - 1260*x**2 + 35)/128),
        lambda x: 1/(38**0.5) * (19*x*(12155*x**8 - 25740*x**6 + 18018*x**4 - 4620*x**2 + 315)/128),
        lambda x: 1/(42**0.5) * (46189*x**12 - 88179*x**10 + 48450*x**8 - 8280*x**6 + 462*x**4 - 7*x**2)
        ]
        return d_lobatto
    
    def get_der_shape_functions(self):
        B_lobatto = []
        for i in range(self.order+1):
            B_lobatto.append(self.derivative()[i])
            
        return B_lobatto       
    
class Larange(BasePolynomial):
    def polynomial(self):
        larange = np.zeros((5))
        larange[0] = lambda x: (1-x)/2
        larange[1] = lambda x: (1+x)/2

        larange[2] = lambda x: (x-0)*(x-1)/(-1-0)*(-1-1)
        larange[3] = lambda x: (x+1)*(x-1)/(0+1)*(0-1)
        larange[4] = lambda x: (x+1)*(x-0)/(1+1)*(1-0)
        

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
        d_larange[0] = lambda x: -1/2
        d_larange[1] = lambda x: 1/2

        d_larange[2] = lambda x: (2*x-1)/(-1-0)*(-1-1)
        d_larange[3] = lambda x: (2*x)/(0+1)*(0-1)
        d_larange[4] = lambda x: (2*x+1)/(1+1)*(1-0)

            
        return 
    
    def get_der_shape_functions(self):
        B_larange = np.zeros((self.order+1))
        if self.order == 1:
            B_larange[0] = self.derivative()[0]
            B_larange[1] = self.derivative()[1]
        elif self.order == 2:
            B_larange[0] = self.derivative()[2]
            B_larange[1] = self.derivative()[3]
            B_larange[2] = self.derivative()[4]
        else:
            print("cubic larange not supported yet")

        return B_larange
class PolyBuilder:
    """build polynomial class"""
    def __init__(self, order):
        self.order = order

    def build(self, poly_type):
        if poly_type == 'lobatto':
            return Lobatto(self.order,)
        elif poly_type == 'larange':
            return Larange(self.order)
        else:
            raise ValueError('poly_type must be lobatto or larange')
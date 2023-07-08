from abc import abstractmethod, ABCMeta
from polynomial import Lobatto, Larange
import numpy as np

class BaseElementaryMatrix(metaclass=ABCMeta):
    """base abstract FE elementary matrix class
    precompute the shape function and its derivative on gauss points
    parameters:
    order: int
        element order
    """
    def Jacobian(self, coord):
        coord = np.array(coord)
        if len(coord) == 2:
            return 2/np.abs(coord[0]-coord[1])
        elif len(coord) == 3:
            print("not implemented yet")

    def inverse_Jacobian(self, coord):
        coord = np.array(coord)
        if len(coord) == 2:
            return np.abs(coord[0]-coord[1])/2
        elif len(coord) == 3:
            print("not implemented yet")


    @abstractmethod
    def matrix(self, order):
        pass

class Elementary1DStiffnessMatrix(BaseElementaryMatrix):
    """FE elementary stiffness matrix class
    parameters:
    order: int
        element order
    """
    def __init__(self, order):
        self.order = order
        self.gauss = Lobatto(order, np.linspace(-1, 1, order+1))
        self.gauss_d = self.gauss.derivative()

    def matrix(self):
        """compute the elementary stiffness matrix
        returns:
        K: ndarray
            elementary stiffness matrix
        """
        K = np.zeros((self.order+1, self.order+1))
        for i in range(self.order+1):
            for j in range(self.order+1):
                K[i, j] = np.sum(self.gauss_d[i]*self.gauss_d[j])
        return K
    
class Elementary1DMassMatrix(BaseElementaryMatrix):
    """FE elementary mass matrix class
    parameters:
    order: int
        element order
    """
    def __init__(self, order):
        self.order = order
        self.gauss = Lobatto(order, np.linspace(-1, 1, order+1))
        self.gauss_d = self.gauss.derivative()

    def matrix(self):
        """compute the elementary mass matrix
        returns:
        M: ndarray
            elementary mass matrix
        """
        M = np.zeros((self.order+1, self.order+1))
        for i in range(self.order+1):
            for j in range(self.order+1):
                M[i, j] = np.sum(self.gauss[i]*self.gauss[j])
        return M
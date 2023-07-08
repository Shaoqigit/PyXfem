from abc import abstractmethod, ABCMeta
from polynomial import Lobatto, Larange, PolyBuilder
from quadratures import GaussLegendreQuadrature
import numpy as np
from precompute_matrices import Ke1Do1, Me1Do1, Ke1Do2, Me1Do2, Ke1Do3, Me1Do3, Ke1Do4, Me1Do4

            

class Base1DElement(metaclass=ABCMeta):
    """base abstract FE elementary matrix class
    precompute the shape function and its derivative on gauss points
    parameters:
    order: int
        element order
    """
    def __init__(self, order):
        self.order = order
        self.dim = 0
        self.nodes = None

    @property
    def Jacobian(self):

        nodes = elem.nodes()
        if elem.dim() == 1:
            return 2/np.abs(nodes[0]-nodes[1])
        else:
            print("not implemented yet")

    @property
    def inverse_Jacobian(self, coord):
        coord = np.array(coord)
        if len(coord) == 2:
            return np.abs(coord[0]-coord[1])/2
        elif len(coord) == 3:
            print("not implemented yet")



class Elementary1DLobattoStiffnessMatrix(BaseElementaryMatrix):
    """FE elementary stiffness matrix class
    parameters:
    order: int
        element order
    """
    def __init__(self, order):
        super().__init__(order)

    def get_matrix(self):
        """compute the elementary stiffness matrix
        returns:
        K: ndarray
            elementary stiffness matrix
        """
        if self.order = 1:
            Ke =  self.Jacobian(elem)*Ke1Do1
        return Ke
    
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
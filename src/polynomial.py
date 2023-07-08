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
    def __init__(self, order, x):
        self.order = order
        self.x = x

    @abstractmethod
    def polynomial(self):
        pass

    @abstractmethod
    def derivative(self):
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
        if self.order == 0:
            p_lobatto = 1-self.x
            p_lobatto /= 2
        elif self.order == 1:
            p_lobatto = 1+self.x
            p_lobatto /= 2
        elif self.order == 2:
            p_lobatto = -1+self.x**2
            p_lobatto *= np.sqrt(3/2)/2
        elif self.order == 3:
            p_lobatto = self.x*(-1+self.x**2)
            p_lobatto *= np.sqrt(5/2)/2
        elif self.order == 4:
            p_lobatto = 1+self.x**2*(-6+5*self.x**2)
            p_lobatto *= np.sqrt(7/2)/8
        elif self.order == 5:
            p_lobatto = self.x*(3+self.x**2*(-10+7*self.x**2))
            p_lobatto *= (3/8)*np.sqrt(2)
        elif self.order == 6:
            p_lobatto = (-1)+self.x**2*(15+self.x**2*((-35)+21*self.x**2))
            p_lobatto *= (1/16)*np.sqrt(11/2)
        elif self.order == 7:
            p_lobatto = self.x*((-5)+self.x**2*(35+self.x**2*((-63)+33*self.x**2)))
            p_lobatto *= (1/16)*np.sqrt(13/2)
        
        return p_lobatto
    
    def derivative(self):
        if self.order == 0:
            d_lobatto = -self.x**0
            d_lobatto /= 2
        elif self.order == 1:
            d_lobatto = self.x**0
            d_lobatto /= 2
        elif self.order == 2:
            d_lobatto = 2*self.x
            d_lobatto *= np.sqrt(3/2)/2
        elif self.order == 3:
            d_lobatto = -1+3*self.x**2
            d_lobatto *= np.sqrt(5/2)/2
        elif self.order == 4:
            d_lobatto = 2*self.x*(-3+5*self.x**2)
            d_lobatto *= np.sqrt(7/2)/8
        elif self.order == 5:
            d_lobatto = -5+7*self.x**2*(3+self.x**2)
            d_lobatto *= (3/8)*np.sqrt(2)
        elif self.order == 6:
            d_lobatto = 2*self.x*(-15+7*self.x**2*(-5+3*self.x**2))
            d_lobatto *= (1/16)*np.sqrt(11/2)
        elif self.order == 7:
            d_lobatto = -35+7*self.x**2*(21+self.x**2*(-9+11*self.x**2))
            d_lobatto *= (1/16)*np.sqrt(13/2)

        return d_lobatto
    
class Larange(BasePolynomial):
    def polynomial(self):
        if self.order == 0:
            p_larange = 1-self.x
            p_larange /= 2
        elif self.order == 1:
            p_larange = 1+self.x
            p_larange /= 2
        elif self.order == 2:
            p_larange = 1+self.x*(-1+self.x)
        elif self.order == 3:

          

        return super().polynomial()

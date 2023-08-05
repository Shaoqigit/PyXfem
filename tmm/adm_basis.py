import numpy as np
from math import pi, sin, cos, exp


class admittance_element:

    def __init__(self, mat, omega, nodes):
        self.omega = omega
        self.nodes = nodes
        self.mat = mat
        self.thickness = np.abs(self.nodes[0]-self.nodes[1])
        self.k = self.omega/self.mat.c_f
        self.transfer_matrix()

    def transfer_matrix(self):
        self.tm = np.array([[np.cos(self.k*self.thickness), 1j*self.mat.Z_f*np.sin(self.k*self.thickness)], 
                            [1j*np.sin(self.k*self.thickness)/self.mat.Z_f, np.cos(self.k*self.thickness)]], dtype=np.complex128)
        return self.tm
    
    def admittance(self):
        tm_11 = self.tm[0,0]
        tm_12 = self.tm[0,1]
        tm_21 = self.tm[1,0]
        tm_22 = self.tm[1,1]
        self.adm = 1/tm_12*np.array([[-tm_22, -tm_21*tm_12+tm_11*tm_22],
                              [1, -tm_11]], dtype=np.complex128)
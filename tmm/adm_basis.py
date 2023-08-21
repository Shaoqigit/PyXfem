import numpy as np
from math import pi, sin, cos, exp


class AdmFluid:

    def __init__(self, mat, omega, theta, k_0, nodes, mode='continue'):
        self.omega = omega
        self.nodes = nodes
        self.mat = mat
        self.thickness = np.abs(self.nodes[0]-self.nodes[1])
        self.k = self.omega/self.mat.c_f
        ky = k_0*np.sin(theta*np.pi/180)
        self.kx = np.sqrt(self.k**2-ky**2)
        if mode=='continue':
            self.transfer_matrix()
        elif mode=='discrete':
            self.transfer_matrix_2()

    def transfer_matrix(self):
        self.tm = np.array([[np.cos(self.kx*self.thickness), -1*self.omega*self.mat.Z_f*np.sin(self.kx*self.thickness)], 
                            [np.sin(self.kx*self.thickness)/(self.mat.Z_f*self.omega), np.cos(self.kx*self.thickness)]], dtype=np.complex128)
        # import pdb; pdb.set_trace()
        return self.tm
    
    def transfer_matrix_2(self):
        alpha = np.array(
            [[0, self.mat.rho_f*self.omega**2],
             [-1/self.mat.K_f+self.kx**2/(self.mat.rho_f*self.omega**2), 0]])
        identity = np.ones((2,2), dtype=np.complex128)
        self.tm = identity-self.thickness*alpha
        # import pdb; pdb.set_trace()
        # self.tm = np.array([[np.exp(0*self.thickness), np.exp(-self.mat.rho_f*self.omega**2*self.thickness)], 
        #                     [np.exp(-1*(-1/self.mat.K_f+self.k**2/(self.mat.rho_f*self.omega**2))*self.thickness), np.exp(0*self.thickness)]], dtype=np.complex128)
        return self.tm
    
    def admittance(self):
        tm_11 = self.tm[0,0]
        tm_12 = self.tm[0,1]
        tm_21 = self.tm[1,0]
        tm_22 = self.tm[1,1]
        # import pdb; pdb.set_trace()
        adm_12 = -tm_21*tm_12+tm_11*tm_22
        self.adm = 1/tm_12*np.array([[-tm_22, 1.],
                              [1., -tm_11]], dtype=np.complex128)
        
    

class AdmElastic:

    def __init__(self, mat, omega, theta, k_0, nodes, mode='continue'):
        self.omega = omega
        self.nodes = nodes
        self.mat = mat
        self.thickness = np.abs(self.nodes[0]-self.nodes[1])
        self.k = self.omega/self.mat.c_f
        ky = k_0*np.sin(theta*np.pi/180)
        self.kx = np.sqrt(self.k**2-ky**2)
        if mode=='continue':
            self.transfer_matrix()
        elif mode=='discrete':
            self.transfer_matrix_2()

    def transfer_matrix(self):
        self.tm = np.array([[np.cos(self.kx*self.thickness), -1*self.omega*self.mat.Z_f*np.sin(self.kx*self.thickness)], 
                            [np.sin(self.kx*self.thickness)/(self.mat.Z_f*self.omega), np.cos(self.kx*self.thickness)]], dtype=np.complex128)
        # import pdb; pdb.set_trace()
        return self.tm
    
    def transfer_matrix_2(self):
        alpha = np.array(
            [[0, self.mat.rho_f*self.omega**2],
             [-1/self.mat.K_f+self.kx**2/(self.mat.rho_f*self.omega**2), 0]])
        identity = np.ones((2,2), dtype=np.complex128)
        self.tm = identity-self.thickness*alpha
        # import pdb; pdb.set_trace()
        # self.tm = np.array([[np.exp(0*self.thickness), np.exp(-self.mat.rho_f*self.omega**2*self.thickness)], 
        #                     [np.exp(-1*(-1/self.mat.K_f+self.k**2/(self.mat.rho_f*self.omega**2))*self.thickness), np.exp(0*self.thickness)]], dtype=np.complex128)
        return self.tm
    
    def admittance(self):
        tm_11 = self.tm[0,0]
        tm_12 = self.tm[0,1]
        tm_21 = self.tm[1,0]
        tm_22 = self.tm[1,1]
        # import pdb; pdb.set_trace()
        adm_12 = -tm_21*tm_12+tm_11*tm_22
        self.adm = 1/tm_12*np.array([[-tm_22, 1.],
                              [1., -tm_11]], dtype=np.complex128)
        

class AdmPoroElastic:

    def __init__(self, mat, omega, theta, k_0, nodes, mode='continue'):
        self.omega = omega
        self.nodes = nodes
        self.mat = mat
        self.thickness = np.abs(self.nodes[0]-self.nodes[1])
        self.k = self.omega/self.mat.c_f
        ky = k_0*np.sin(theta*np.pi/180)
        self.kx = np.sqrt(self.k**2-ky**2)
        if mode=='continue':
            self.transfer_matrix()
        elif mode=='discrete':
            self.transfer_matrix_2()

    def transfer_matrix(self):
        self.tm = np.array([[np.cos(self.kx*self.thickness), -1*self.omega*self.mat.Z_f*np.sin(self.kx*self.thickness)], 
                            [np.sin(self.kx*self.thickness)/(self.mat.Z_f*self.omega), np.cos(self.kx*self.thickness)]], dtype=np.complex128)
        # import pdb; pdb.set_trace()
        return self.tm
    
    def transfer_matrix_2(self):
        alpha = np.array(
            [[0, self.mat.rho_f*self.omega**2],
             [-1/self.mat.K_f+self.kx**2/(self.mat.rho_f*self.omega**2), 0]])
        identity = np.ones((2,2), dtype=np.complex128)
        self.tm = identity-self.thickness*alpha
        # import pdb; pdb.set_trace()
        # self.tm = np.array([[np.exp(0*self.thickness), np.exp(-self.mat.rho_f*self.omega**2*self.thickness)], 
        #                     [np.exp(-1*(-1/self.mat.K_f+self.k**2/(self.mat.rho_f*self.omega**2))*self.thickness), np.exp(0*self.thickness)]], dtype=np.complex128)
        return self.tm
    
    def admittance(self):
        tm_11 = self.tm[0,0]
        tm_12 = self.tm[0,1]
        tm_21 = self.tm[1,0]
        tm_22 = self.tm[1,1]
        # import pdb; pdb.set_trace()
        adm_12 = -tm_21*tm_12+tm_11*tm_22
        self.adm = 1/tm_12*np.array([[-tm_22, 1.],
                              [1., -tm_11]], dtype=np.complex128)
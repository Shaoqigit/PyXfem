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
        p = self.mat.phi
        R0 = self.mat.sigma
        rho_f = self.mat.rho_f
        rho_s = self.mat.rho_1
        t = self.mat.alpha
        K_f = self.mat.K_f
        biot=1
        lame2 = self.mat.E/(2*(1+self.mat.nu))
        lame1 = self.mat.E*self.mat.nu/((1+self.mat.nu)*(1-2*self.mat.nu))
        K_s = lame1+2*lame2

        b = p*R0
        rhot=p*rho_f*(1-t)
        rhoa=-rhot
        rho12=rhot+1j*b/self.omega
        rho22=p*rho_f-rho12
        rho1=(1-p)*rho_s
        rho11=rho1-rho12
        rho=rho11-rho12**2/rho22  # rho_til

        rg=p*rho12/rho22-(1-biot)  # gamma_til
        KKs=K_s
        KKf=p*p/rho22
        MMs=rho
        MMf=p/K_f

        G=KKs*MMf-KKf*MMs+rg**2
        D=G*G-4.0*KKs*KKf*MMs*MMf
        rac2D=np.sqrt(D)

        c1=np.sqrt((G+rac2D)/(2.0*MMs*MMf))
        c2=np.sqrt((G-rac2D)/(2.0*MMs*MMf))
        c3=np.sqrt(lame2/rho)

        fa = -1j*self.omega*(c1*rho-K_s/c1)/rg  #p/u factor for fast c-wave
        fb = -1j*self.omega*(c2*rho-K_s/c2)/rg  #p/u factor for slow c-wave
        va=KKf*(fa/c1-1j*self.omega*rho_f)  #v/u factor for fast c-wave
        vb=KKf*(fb/c2-1j*self.omega*rho_f)  #v/u factor for slow c-wave
        sa=-1j*self.omega/c1*K_s-biot*fa  #sigma/u factor for fast c-wave
        sb=-1j*self.omega/c2*K_s-biot*fb  #sigma/u factor for slow c-wave

        miw=-1j*self.omega

        h=self.thickness
        kx=self.k_0


        tm_11 = self.tm[0,0]
        tm_12 = self.tm[0,1]
        tm_21 = self.tm[1,0]
        tm_22 = self.tm[1,1]
        # import pdb; pdb.set_trace()
        adm_12 = -tm_21*tm_12+tm_11*tm_22
        self.adm = 1/tm_12*np.array([[-tm_22, 1.],
                              [1., -tm_11]], dtype=np.complex128)
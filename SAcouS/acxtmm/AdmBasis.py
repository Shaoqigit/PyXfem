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
        

class AdmPoroElastic2:

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
        alpha = np.zeros((6,6),dtype=complex)
        alpha[0,3] = 1j*self.kx*self.mat.A_hat/self.mat.P_hat
        alpha[0,4] = 1j*self.mat.gamma_til*self.kx
        alpha[0,5] = -(self.mat.A_hat**2-self.mat.P_hat**2)*self.kx**2/self.mat.P_hat-self.mat.rho_til*self.omega**2

        alpha[1,3] = 1/self.mat.P_hat
        alpha[1,5] = 1j*self.kx*self.mat.A_hat/self.mat.P_hat

        alpha[2,4] = -1./self.mat.K_eq_til+self.kx**2/(self.mat.rho_eq_til*self.omega**2)
        alpha[2,5] = -1j*self.kx*self.mat.gamma_til

        alpha[3,0] = 1j*self.kx
        alpha[3,1] = -self.mat.rho_s_til*self.omega**2
        alpha[3,2] = -self.mat.rho_eq_til*self.mat.gamma_til*self.omega**2

        alpha[4,1] = self.mat.rho_eq_til*self.mat.gamma_til*self.omega**2
        alpha[4,2] = self.mat.rho_eq_til*self.omega**2

        alpha[5,0] = 1./self.mat.N
        alpha[5,1] = 1j*self.kx

        import scipy.linalg as SA
        self.tm = SA.expm(-self.thickness*alpha)

        self.prim_M = np.zeros_like(alpha)
        self.prim_M[0,1] = -self.tm[0,5-1]
        self.prim_M[0,3] = -self.tm[0,2-1]
        self.prim_M[0,5] = -self.tm[0,6-1]

        self.prim_M[1,1] = -self.tm[2-1,5-1]
        self.prim_M[1,3] = -self.tm[2-1,2-1]
        self.prim_M[1,5] = -self.tm[2-1,6-1]
        self.prim_M[1,2] = 1

        self.prim_M[2,1] = -self.tm[3-1,5-1]
        self.prim_M[2,3] = -self.tm[3-1,2-1]
        self.prim_M[2,5] = -self.tm[3-1,6-1]

        self.prim_M[3,1] = -self.tm[4-1,5-1]
        self.prim_M[3,3] = -self.tm[4-1,2-1]
        self.prim_M[3,5] = -self.tm[4-1,6-1]

        self.prim_M[4,1] = -self.tm[5-1,5-1]
        self.prim_M[4,3] = -self.tm[5-1,2-1]
        self.prim_M[4,5] = -self.tm[5-1,6-1]
        self.prim_M[4,0] = 1

        self.prim_M[5,1] = -self.tm[6-1,5-1]
        self.prim_M[5,3] = -self.tm[6-1,2-1]
        self.prim_M[5,5] = -self.tm[6-1,6-1]
        self.prim_M[5,4] = 1

        self.dual_M = np.zeros_like(alpha)
        self.dual_M[0,1] = self.tm[1-1,3-1]
        self.dual_M[0,3] = self.tm[1-1,4-1]
        self.dual_M[0,5] = self.tm[1-1,1-1]
        self.dual_M[0,4] = -1

        self.dual_M[1,1] = self.tm[2-1,3-1]
        self.dual_M[1,3] = self.tm[2-1,4-1]
        self.dual_M[1,5] = self.tm[2-1,1-1]

        self.dual_M[2,1] = self.tm[3-1,3-1]
        self.dual_M[2,3] = self.tm[3-1,4-1]
        self.dual_M[2,5] = self.tm[3-1,1-1]
        self.dual_M[2,0] = -1

        self.dual_M[3,1] = self.tm[4-1,3-1]
        self.dual_M[3,3] = self.tm[4-1,4-1]
        self.dual_M[3,5] = self.tm[4-1,1-1]
        self.dual_M[3,2] = -1

        self.dual_M[4,1] = self.tm[5-1,3-1]
        self.dual_M[4,3] = self.tm[5-1,4-1]
        self.dual_M[4,5] = self.tm[5-1,1-1]

        self.dual_M[5,1] = self.tm[6-1,3-1]
        self.dual_M[5,3] = self.tm[6-1,4-1]
        self.dual_M[5,5] = self.tm[6-1,1-1]

        return self.tm
    
    def transfer_matrix_2(self):
        alpha = np.array(
            [[0, self.mat.rho_f*self.omega**2],
             [-1/self.mat.K_f+self.kx**2/(self.mat.rho_f*self.omega**2), 0]])
        identity = np.ones((2,2), dtype=np.complex128)
        self.tm = identity-self.thickness*alpha
        return self.tm
    
    def admittance(self):
        adm_full = np.zeros_like(self.tm)
        adm_full = np.linalg.inv(self.dual_M)@self.prim_M
        import pdb;pdb.set_trace()
        A11 = adm_full[:2,:2]
        A12 = adm_full[:2,2:]
        A21 = adm_full[2:,:2]
        A22 = adm_full[2:,2:]

        # shur complement of 6x6 matrix
        self.adm = A11 - A12@np.linalg.inv(A22)@A21

        return self.adm
import numpy as np
import scipy.linalg as SA

class TMMFluid:
    def __init__(self, mat, omega) -> None:
        self.mat = mat
        self.omega = omega
        self.mat.set_frequency(omega)

    def transfer_matrix(self, thickness):
        k = self.omega/self.mat.c_f
        tm = np.array([[np.cos(k*thickness), 1j*self.mat.Z_f*np.sin(k*thickness)], 
                        [1j*np.sin(k*thickness)/self.mat.Z_f, np.cos(k*thickness)]], dtype=np.complex128)
        return tm
    
class TMMElastic:
    def __init__(self, mat, omega) -> None:
        self.mat = mat
        self.omega = omega
        self.mat.set_frequency(omega)

    def transfer_matrix(self, thickness):
        k = self.omega/self.mat.c_f
        tm = np.array([[np.cos(k*thickness), 1j*self.mat.Z_f*np.sin(k*thickness)], 
                        [1j*np.sin(k*thickness)/self.mat.Z_f, np.cos(k*thickness)]], dtype=np.complex128)
        return tm
    

class TMMPoroElastic1:
    def __init__(self, mat, omega, k_0) -> None:
        self.mat = mat
        self.omega = omega
        self.k_0 = k_0
        self.mat.set_frequency(omega)

    def phi_matrix(self, x3):
        k13=np.sqrt(self.mat.delta_1**2-self.k_0**2)
        k23=np.sqrt(self.mat.delta_2**2-self.k_0**2)
        k33=np.sqrt(self.mat.delta_3**2-self.k_0**2)

        D1 = (self.mat.P_til+self.mat.Q_til*self.mat.mu_1)*(self.k_0**2+k13**2)-2*self.mat.N*self.k_0**2
        D2 = (self.mat.P_til+self.mat.Q_til*self.mat.mu_2)*(self.k_0**2+k23**2)-2*self.mat.N*self.k_0**2

        E1 = (self.mat.R_til*self.mat.mu_1+self.mat.Q_til)*(self.k_0**2+k13**2)
        E2 = (self.mat.R_til*self.mat.mu_2+self.mat.Q_til)*(self.k_0**2+k23**2)

        T1 = np.zeros((6,6),dtype=complex)
        T1[0,0]=self.omega*self.k_0*np.cos(k13*x3)
        T1[0,1]=-1j*self.omega*self.k_0*np.sin(k13*x3)
        T1[0,2]=self.omega*self.k_0*np.cos(k23*x3)
        T1[0,3]=-1j*self.omega*self.k_0*np.sin(k23*x3)
        T1[0,4]=1j*self.omega*k33*np.sin(k33*x3)
        T1[0,5]=-1*self.omega*k33*np.cos(k33*x3)

        T1[1,0]=-1j*self.omega*k13*np.sin(k13*x3)
        T1[1,1]=self.omega*k13*np.cos(k13*x3)
        T1[1,2]=-1j*self.omega*k23*np.sin(k23*x3)
        T1[1,3]=self.omega*k23*np.cos(k23*x3)
        T1[1,4]=self.omega*self.k_0*np.cos(k33*x3)
        T1[1,5]=-1j*self.omega*self.k_0*np.sin(k33*x3)

        T1[2,0]=-1j*self.omega*k13*self.mat.mu_1*np.sin(k13*x3)
        T1[2,1]=self.omega*k13*self.mat.mu_1*np.cos(k13*x3)
        T1[2,2]=-1j*self.omega*k23*self.mat.mu_2*np.sin(k23*x3)
        T1[2,3]=self.omega*k23*self.mat.mu_2*np.cos(k23*x3)
        T1[2,4]=self.omega*self.k_0*self.mat.mu_3*np.cos(k33*x3)
        T1[2,5]=-1j*self.omega*self.k_0*self.mat.mu_3*np.sin(k33*x3)

        T1[3,0]=-1*D1*np.cos(k13*x3)
        T1[3,1]=1j*D1*np.sin(k13*x3)
        T1[3,2]=-1*D2*np.cos(k23*x3)
        T1[3,3]=1j*D2*np.sin(k23*x3)
        T1[3,4]=2j*self.mat.N*k33*self.k_0*np.sin(k33*x3)
        T1[3,5]=-2*self.mat.N*k33*self.k_0*np.cos(k33*x3)

        T1[4,0]=2j*self.mat.N*self.k_0*k13*np.sin(k13*x3)
        T1[4,1]=-2*self.mat.N*self.k_0*k13*np.cos(k13*x3)
        T1[4,2]=2j*self.mat.N*self.k_0*k23*np.sin(k23*x3)
        T1[4,3]=-2*self.mat.N*self.k_0*k23*np.cos(k23*x3)
        T1[4,4]=self.mat.N*(k33**2-self.k_0**2)*np.cos(k33*x3)
        T1[4,5]=-1j*self.mat.N*(k33**2-self.k_0**2)*np.sin(k33*x3)

        T1[5,0]=-E1*np.cos(k13*x3)
        T1[5,1]=1j*E1*np.sin(k13*x3)
        T1[5,2]=-E2*np.cos(k23*x3)
        T1[5,3]=1j*E2*np.sin(k23*x3)
        T1[5,4]=0.
        T1[5,5]=0.

        return T1
    
    def transfer_matrix(self, thickness):
        T1 = self.phi_matrix(-thickness)
        T2 = self.phi_matrix(0)
        tm = T1@np.linalg.inv(T2)

        return tm


class TMMPoroElastic2:
    def __init__(self, mat, omega, k_0) -> None:
        self.mat = mat
        self.omega = omega
        self.k_0 = k_0
        self.mat.set_frequency(omega)

    def phi_matrix(self):
        k13=np.sqrt(self.mat.delta_1**2-self.k_0**2)
        k23=np.sqrt(self.mat.delta_2**2-self.k_0**2)
        k33=np.sqrt(self.mat.delta_3**2-self.k_0**2)

        D1 = (self.mat.P_hat+self.mat.Q_til*self.mat.mu_1)*(self.k_0**2+k13**2)-2*self.mat.N*self.k_0**2
        D2 = (self.mat.P_hat+self.mat.Q_til*self.mat.mu_2)*(self.k_0**2+k23**2)-2*self.mat.N*self.k_0**2

        E1 = (self.mat.R_til*self.mat.mu_1+self.mat.Q_til)*(self.k_0**2+k13**2)
        E2 = (self.mat.R_til*self.mat.mu_2+self.mat.Q_til)*(self.k_0**2+k23**2)

        T1 = np.zeros((6,6),dtype=complex)
        T1[0,0]=self.omega*self.k_0
        T1[0,1]=self.omega*self.k_0
        T1[0,2]=self.omega*self.k_0
        T1[0,3]=self.omega*self.k_0
        T1[0,4]=-1*self.omega*k33
        T1[0,5]=1*self.omega*k33

        T1[1,0]=self.omega*k13
        T1[1,1]=-self.omega*k13
        T1[1,2]=self.omega*k23
        T1[1,3]=-1*self.omega*k23
        T1[1,4]=self.omega*self.k_0
        T1[1,5]=self.omega*self.k_0

        T1[2,0]=self.omega*k13*self.mat.mu_1
        T1[2,1]=-1*self.omega*k13*self.mat.mu_1
        T1[2,2]=self.omega*k23*self.mat.mu_2
        T1[2,3]=-1*self.omega*k23*self.mat.mu_2
        T1[2,4]=self.omega*self.k_0*self.mat.mu_3
        T1[2,5]=self.omega*self.k_0*self.mat.mu_3

        T1[3,0]=-1*D1
        T1[3,1]=-1*D1
        T1[3,2]=-1*D2
        T1[3,3]=-1*D2
        T1[3,4]=-2*self.mat.N*k33*self.k_0
        T1[3,5]=2*self.mat.N*k33*self.k_0

        T1[4,0]=-2*self.mat.N*self.k_0*k13
        T1[4,1]=2*self.mat.N*self.k_0*k13
        T1[4,2]=-2*self.mat.N*self.k_0*k23
        T1[4,3]=2*self.mat.N*self.k_0*k23
        T1[4,4]=self.mat.N*(k33**2-self.k_0**2)
        T1[4,5]=self.mat.N*(k33**2-self.k_0**2)

        T1[5,0]=-1*E1
        T1[5,1]=-1*E1
        T1[5,2]=-1*E2
        T1[5,3]=-1*E2
        T1[5,4]=0.
        T1[5,5]=0.

        return T1
    
    def Lambda_matrix(self, x3):
        k13=np.sqrt(self.mat.delta_1**2-self.k_0**2)
        k23=np.sqrt(self.mat.delta_2**2-self.k_0**2)
        k33=np.sqrt(self.mat.delta_3**2-self.k_0**2)
        T2 = np.zeros((6,6),dtype=complex)
        T2[0,0] = np.exp(1j*k13*x3)
        T2[1,1] = np.exp(-1j*k13*x3)
        T2[2,2] = np.exp(1j*k23*x3)
        T2[3,3] = np.exp(-1j*k23*x3)
        T2[4,4] = np.exp(1j*k33*x3)
        T2[5,5] = np.exp(-1j*k33*x3)  
        return T2      
    
    def transfer_matrix(self, thickness):
        # import pdb;pdb.set_trace()
        phi = self.phi_matrix()
        Lambda = self.Lambda_matrix(thickness)
        tm = phi@Lambda@np.linalg.inv(phi)

        return tm

class TMMPoroElastic3:
    def __init__(self, mat, omega, k_0) -> None:
        self.mat = mat
        self.omega = omega
        self.k_0 = k_0
        self.mat.set_frequency(omega)

    def alpha_matrix(self):
        alpha = np.zeros((6,6),dtype=complex)
        alpha[0,3] = 1j*self.k_0*self.mat.A_hat/self.mat.P_hat
        alpha[0,4] = 1j*self.mat.gamma_til*self.k_0
        alpha[0,5] = -(self.mat.A_hat**2-self.mat.P_hat**2)*self.k_0**2/self.mat.P_hat-self.mat.rho_til*self.omega**2

        alpha[1,3] = 1/self.mat.P_hat
        alpha[1,5] = 1j*self.k_0*self.mat.A_hat/self.mat.P_hat

        alpha[2,4] = -1./self.mat.K_eq_til+self.k_0**2/(self.mat.rho_eq_til*self.omega**2)
        alpha[2,5] = -1j*self.k_0*self.mat.gamma_til

        alpha[3,0] = 1j*self.k_0
        alpha[3,1] = -self.mat.rho_s_til*self.omega**2
        alpha[3,2] = -self.mat.rho_eq_til*self.mat.gamma_til*self.omega**2

        alpha[4,1] = self.mat.rho_eq_til*self.mat.gamma_til*self.omega**2
        alpha[4,2] = self.mat.rho_eq_til*self.omega**2

        alpha[5,0] = 1./self.mat.N
        alpha[5,1] = 1j*self.k_0

        return alpha
    

    def transfer_matrix(self, thickness):
        # import pdb;pdb.set_trace()
        d = thickness
        alpha = self.alpha_matrix()
        tm = SA.expm(-d*alpha)
        # tm = np.identity(6, dtype=complex) - d*alpha

        return tm
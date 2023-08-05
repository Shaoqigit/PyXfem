import numpy as np
from numpy.lib.scimath import sqrt
from abc import abstractmethod, ABCMeta


class BaseMaterial(metaclass=ABCMeta):
    """base abstract material class
    parameters:
    """
    TYPE = 'Abstract Mtaterial'
    MODEL = 'Abstract Model'

    def __init__(self, name,  *args):
        self.name = name

    def set_frequency(self, omega):
        pass

    def __str__(self) -> str:
        return f"name: {self.name}, type: {self.__class__.TYPE}, model: {self.__class__.MODEL}"

class Air(BaseMaterial):
    """air material class
    parameters:
    """
    TYPE = 'Air'
    MODEL = 'Air Model'
    COMPATIBLE = ['Air', 'Fluid', 'Equivalent Fluid', 'Limp Equivalent Fluid']

    # atmospheric conditions
    T = 293.15  # reference temperature [K]
    P = 1.01325e5  # atmospheric Pressure [Pa]
    gamma = 1.400  # polytropic coefficient []
    lambda_ = 0.0262  # thermal conductivity [W.m^-1.K^-1]
    mu = 0.1839e-4  # dynamic viscosity [kg.m^-1.s^-1]
    Pr = 0.710  # Prandtl's number []
    molar_mass = 0.29e-1  # molar mass [kg.mol^-1]
    rho = 1.213  # density [kg.m^-3]
    C_p = 1006  # (mass) specific heat capacity as constant pressure [J.K^-1]

    K = gamma*P  # adiabatic bulk modulus
    c = np.sqrt(K/rho)  # adiabatic sound speed
    Z = rho*c  # characteristic impedance
    C_v = C_p/gamma  # (mass) specific heat capacity as constant volume [J.K^-1]
    nu = mu/rho  # kinematic viscosity [m.s^-2]
    nu_prime = nu/Pr  # viscothermal losses

    def __init__(self, name, *args):
        super().__init__(name, *args)
        assert(len(args) == 0)
        self.name = name
        self.K = Air.K
        self.c = Air.c
        self.Z = Air.Z
        self.rho = Air.rho

        # derived parameters
        self.rho_f = self.rho
        self.c_f = self.c
        self.Z_f = self.Z

    def set_frequency(self, omega):
        pass


class Fluid(BaseMaterial):
    """fluid material class
    parameters:
    """
    TYPE = 'Fluid'
    MODEL = 'Fluid Model'
    COMPATIBLE = ['Air', 'Fluid', 'Equivalent Fluid', 'Limp Equivalent Fluid']


    def __init__(self, name, *args):
        super().__init__(name, *args)
        assert(len(args) == 2)
        self.name = name
        self.rho = args[0]
        self.c = args[1]

        self.rho_f = self.rho
        self.c_f = self.c
        self.Z_f = self.rho*self.c

    def set_frequency(self, omega):
        return super().set_frequency(omega)


class EquivalentFluid(BaseMaterial):
    """equivalent fluid material class
    attributes:
    phi: porosity -> float
    sigma: flow resistivity -> float
    alpha: static tortuosity -> float
    Lambda_prime: thermal characteristic length -> float
    Lambda: viscous characteristic length -> float
    """
    TYPE = 'Equivalent Fluid'
    MODEL = 'JCAL Equivalent Fluid'
    COMPATIBLE = ['Air', 'Fluid', 'Equivalent Fluid', 'Limp Equivalent Fluid']


    def __init__(self, name, *args):
        super().__init__(name, *args)
        # assert(len(args) == 5)
        self.name = name
        self.phi = args[0]
        self.sigma = args[1]
        self.alpha = args[2]
        self.Lambda_prime = args[3]
        self.Lambda = args[4]

    def set_frequency(self, omega):
        """set frequency"""
        #  Johnson et al model for rho_eq_til
        self.omega_0 = self.sigma*self.phi/(Air.rho*self.alpha)
        self.omega_infty = (self.sigma*self.phi*self.Lambda)**2/(4*Air.mu*Air.rho*self.alpha**2)
        self.F_JKD = sqrt(1+1j*omega/self.omega_infty)
        self.rho_eq_til = (Air.rho*self.alpha/self.phi)*(1+(self.omega_0/(1j*omega))*self.F_JKD)
        self.alpha_til = self.phi*self.rho_eq_til/Air.rho

        #  Champoux-Allard model for K_eq_til
        self.omega_prime_infty = (16*Air.nu_prime)/(self.Lambda_prime**2)
        self.F_prime_CA = sqrt(1+1j*omega/self.omega_prime_infty)
        self.alpha_prime_til = 1+self.omega_prime_infty*self.F_prime_CA/(2*1j*omega)
        self.K_eq_til = (Air.gamma*Air.P/self.phi)/(Air.gamma-(Air.gamma-1)/self.alpha_prime_til)

        self.c_eq_til = sqrt(self.K_eq_til/self.rho_eq_til)

        self.rho_f = self.rho_eq_til
        self.c_f = self.c_eq_til 
        self.Z_f = self.rho_f*self.c_f


class LimpPorousMaterial(EquivalentFluid):
    """limp porous material class
    derived from equivalent fluid class
    take the solid density and inertia into account
    additional attributes:
    rho_1: solid density -> float
    """
    TYPE = 'Limp Fluid'
    MODEL = 'Limp Equivalent Fluid'
    COMPATIBLE = ['Air', 'Fluid', 'Equivalent Fluid', 'Limp Equivalent Fluid']


    def __init__(self, name, *args):
        super().__init__(name, *args)
        # assert(len(args) == 6)
        self.name = name
        self.rho_1 = args[5]  #

    def set_frequency(self, omega):
        super().set_frequency(omega)
        self.rho_12 = -self.phi*Air.rho*(self.alpha-1)
        self.rho_11 = self.rho_1-self.rho_12
        self.rho_2 = self.phi*Air.rho
        self.rho_22 = self.rho_2-self.rho_12

        self.rho_22_til = self.phi**2*self.rho_eq_til
        self.rho_12_til = self.rho_2-self.rho_22_til
        self.rho_11_til = self.rho_1-self.rho_12_til
        self.rho_til = self.rho_11_til-((self.rho_12_til**2)/self.rho_22_til)

        self.gamma_til = self.phi*(self.rho_12_til/self.rho_22_til-(1-self.phi)/self.phi)
        self.rho_s_til = self.rho_til+self.gamma_til**2*self.rho_eq_til

        self.rho_limp = self.rho_til*self.rho_eq_til/(self.rho_til+self.rho_eq_til*self.gamma_til**2)

        self.rho_f = self.rho_limp
        self.c_f= self.c_eq_til
        self.Z_f = self.rho_f*self.c_f

class ElasticMaterial(BaseMaterial):
    """elastic material class
    parameters:
    """
    TYPE = 'Elastic'
    MODEL = 'Elastic Model'
    COMPATIBLE = ['Elastic']


    def __init__(self, name, *args):
        super().__init__(name, *args)
        assert(len(args) == 4)
        self.name = name
        self.E = args[0]  # Young's modulus
        self.rho_1 = args[1]  # solid density
        self.nu = args[2]  # Poisson's ratio
        self.eta = args[3]  # loss factor

    def compute_missing(self):
        self.E *= (1+1j*self.eta)
        self.mu = (1+1j*self.eta)*(self.E)/(2*(1+self.nu))
        self.lambda_ = (self.E*self.nu)/((1+self.nu)*(1-2*self.nu))

    def set_frequency(self, omega):
        self.omega = omega
        P_mat = self.lambda_ + 2*self.mu
        self.delta_p = omega*sqrt(self.rho/P_mat)
        self.delta_s = omega*sqrt(self.rho/self.mu)

class PoroElasticMaterial(LimpPorousMaterial):
    TYPE = 'Poroelastic'
    MODEL = 'JCA-Biot Model'
    COMPATIBLE = ['Poroelastic', 'Elastic', 'Fluid', 'Equivalent Fluid', 'Limp Equivalent Fluid']

    def __init__(self, name, *args):
        super().__init__(name, *args)
        # assert(len(args) == 8)
        self.name = name
        self.E = args[6]
        self.nu = args[7]
        self.eta = args[8]
        
    def set_frequency(self, omega):
        super().set_frequency(omega)
        self.structural_loss = 1+1j*self.eta

        self.N = self.E/(2*(1+self.nu))*self.structural_loss
        self.A_hat = (self.E*self.nu)/((1+self.nu)*(1-2*self.nu))*self.structural_loss
        self.P_hat = self.A_hat+2*self.N

        # Biot 1956 elastic coefficients
        self.R_til = self.K_eq_til*self.phi**2
        self.Q_til = ((1-self.phi)/self.phi)*self.R_til
        self.P_til = self.P_hat+self.Q_til**2/self.R_til

        delta_eq = omega*sqrt(self.rho_eq_til/self.K_eq_til)
        delta_s_1 = omega*sqrt(self.rho_til/self.P_hat)
        delta_s_2 = omega*sqrt(self.rho_s_til/self.P_hat)

        Psi = ((delta_s_2**2+delta_eq**2)**2-4*delta_eq**2*delta_s_1**2)
        sdelta_total = sqrt(Psi)

        delta_1 = sqrt(0.5*(delta_s_2**2+delta_eq**2+sdelta_total))
        delta_2 = sqrt(0.5*(delta_s_2**2+delta_eq**2-sdelta_total))
        delta_3 = omega*sqrt(self.rho_til/self.N)

        mu_1 = self.gamma_til*delta_eq**2/(delta_1**2-delta_eq**2)
        mu_2 = self.gamma_til*delta_eq**2/(delta_2**2-delta_eq**2)
        mu_3 = -self.gamma_til

        self.delta_s_1 = delta_s_1
        self.delta_s_2 = delta_s_2
        self.delta_1 = delta_1
        self.delta_2 = delta_2
        self.delta_3 = delta_3
        self.delta_eq = delta_eq
        self.sdelta_total = sdelta_total
        self.mu_1 = mu_1
        self.mu_2 = mu_2
        self.mu_3 = mu_3

        self.rho_f = self.rho_eq_til
        self.c_f = self.c_eq_til



# TODO: correct Limp and BiotMaterial class
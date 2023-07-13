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
        assert(len(args) == 5)
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

class LimpPorousMaterial(EquivalentFluid):
    """limp porous material class
    derived from equivalent fluid class
    take the solid density and inertia into account
    additional attributes:
    rho_til: solid density -> float
    gamma_til: biot coupling coefficient -> float
    """
    TYPE = 'Limp Fluid'
    MODEL = 'Limp Equivalent Fluid'
    COMPATIBLE = ['Air', 'Fluid', 'Equivalent Fluid', 'Limp Equivalent Fluid']


    def __init__(self, name, *args):
        super().__init__(name, *args)
        assert(len(args) == 7)
        self.name = name
        self.rho_til = args[5]
        self.gamma_til = args[6]

    def set_frequency(self, omega):
        self.rho_limp = self.rho_til*self.rho_eq_til/(self.rho_til+self.rho_eq_til*self.gamma_til**2)



def check_material_compability(subdomains):
    mats = []
    for key in subdomains.keys():
        mats.append(key)

    compatibale_mats = mats[0].COMPATIBLE
    print(compatibale_mats)
    for mat in mats:
        print(mat.TYPE)
        if mat.TYPE not in compatibale_mats:
            raise ValueError("Material model is not compatible")
        else:
            print("Material models are compatible, computation continues ...")

# TODO: add ElasticMaterial class and BiotMaterial class
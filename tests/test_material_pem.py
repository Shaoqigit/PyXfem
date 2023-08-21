# This file is part of PyXfem, a software distributed under the MIT license.
# For any question, please contact the authors cited below.
#
# Copyright (c) 2023
# 	Shaoqi WU <shaoqiwu@outlook.com
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# Main Test case 
import sys
sys.path.append('/home/shaoqi/Devlop/PyXfem/PyAcoustiX/')
import numpy as np
import matplotlib.pyplot as plt
from pymls import *
from mediapack import Air, PEM, EqFluidJCA

from fem.materials import Air, Fluid, EquivalentFluid, PoroElasticMaterial


def test_case():
    # porous parameters
    phi = 0.4 # porosity
    sigma = 4e6 # flow resistivity
    alpha = 1.75 # static tortuosity
    Lambda_prime = 2.0e-5 # thermal characteristic length
    Lambda = 9.3e-6 # viscous characteristic length
    rho_1 = 120 # solid density
    nu = 0.3 # Poisson's ratio
    E = 4e4 # Young's modulus
    eta = 0.2 # dynamic viscosity
    loss_type = 'structural' # loss type

    freq = 1000
    omega = np.pi*2*freq

    # test equivalent fluid
    foam = EqFluidJCA()
    foam.phi = phi
    foam.sigma = sigma
    foam.alpha = alpha
    foam.Lambda_prime = Lambda_prime
    foam.Lambda = Lambda


    foam2 = EquivalentFluid('foam', phi, sigma, alpha, Lambda_prime, Lambda)


    foam.update_frequency(omega)
    Z = foam.rho_eq_til * foam.c_eq_til

    foam2.set_frequency(omega)
    Z2 = foam2.rho_f * foam2.c_f
    print("Impedance from reference:", Z)
    print("Impedance from pyxfem:", Z2)
    test_equivalent_fluid = False
    if Z == Z2:
        print("Equivalent fluid test passed!")
        test_equivalent_fluid = True


    foam = PEM()
    foam.phi = phi
    foam.sigma = sigma
    foam.alpha = alpha
    foam.Lambda_prime = Lambda_prime
    foam.Lambda = Lambda
    foam.rho_1 = rho_1
    foam.nu = nu
    foam.E = E
    foam.eta = eta
    foam.loss_type = 'structural'
    foam._compute_missing()
    foam.update_frequency(omega)
    Z = foam.rho_eq_til * foam.c_eq_til

    foam2 = PoroElasticMaterial('foam', phi, sigma, alpha, Lambda_prime, Lambda, rho_1, E, nu, eta)
    foam2.set_frequency(omega)
    Z2 = foam2.rho_f * foam2.c_f
    print("Impedance from reference:", foam.delta_1, foam.delta_2, foam.delta_3, foam.rho_til, foam.P_til)
    print("Impedance from reference:", foam2.delta_1, foam2.delta_2, foam2.delta_3, foam2.rho_til, foam2.P_til)
    delta_reference = [foam.delta_1, foam.delta_2, foam.delta_3]
    delta_lib = [foam2.delta_1, foam2.delta_2, foam2.delta_3]
    for i in range(3):
        if delta_reference[i] != delta_lib[i]:
            print("delta_{} is not equal!".format(i+1))
            test_poroelastic_material = False
            break
        else:
            test_poroelastic_material = True

    if test_equivalent_fluid and test_poroelastic_material:
        print("Test passed!")
        return True 
    else:
        print("Test failed!")
        return False


if __name__ == "__main__":
    result = test_case()
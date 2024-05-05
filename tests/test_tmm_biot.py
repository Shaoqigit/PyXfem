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
import os
current_dir = os.path.dirname(os.path.realpath(__file__))
working_dir = os.path.join(current_dir , "..")
import sys
sys.path.append(working_dir)

import numpy as np
import matplotlib.pyplot as plt

from SAcouS.acxtmm.adm_basis import AdmFluid
from SAcouS.acxtmm.tmm import TMMFluid, TMMPoroElastic1, TMMPoroElastic2, TMMPoroElastic3
from SAcouS.acxtmm.BC_matrix import bcm_poro_fluid, bcm_fluid_poro2, bcm_poro_rigid_wall, bcm_rigid_wall, bcm_rigid_wall2
from SAcouS.acxtmm.adm_assembler import AdmAssembler
from SAcouS.acxfem.mesh import Mesh1D
from SAcouS.acxfem.materials import Air, Fluid, EquivalentFluid, PoroElasticMaterial
from SAcouS.acxfem.utilities import check_material_compability, display_matrix_in_array, plot_matrix_partten
from SAcouS.acxfem.solver import AdmittanceSolver
from SAcouS.acxfem.postprocess import PostProcessField
from analytical.fluid_sol import DoubleleLayerKundltTube





# mesh.refine_mesh(1)


# ====================== Pysical Problem ======================
# define the materials
air = Air('classical air')

# given JCA porous material properties
phi          = 0.4  # porosity
sigma        = 4e6 # resistivity
alpha        = 1.75  # Tortuosity
Lambda_prime = 2.0e-5  # Viscous characteristic length
Lambda      = 9.3e-6  # 
rho_1 = 120
nu = 0.3
E = 4e4
eta = 0.2
xfm = PoroElasticMaterial('xfm', phi, sigma, alpha, Lambda_prime, Lambda, rho_1, E, nu, eta)
# import pdb; pdb.set_trace()
def discrete_in_frequency(freq):

    omega = 2 * np.pi * freq  # angular frequency
    lambda_ = air.c_f / freq  # wavelength
    h = lambda_/10  # size of elements
    num_elem = int(0.2 / h) +10 # number of elements
    # import pdb; pdb.set_trace()
    num_nodes = num_elem + 1  # number of nodes
    nodes = np.linspace(-0.01, 0.01, num_nodes)

    elem_connec1 = np.arange(0, num_elem)
    elem_connec2 = np.arange(1, num_nodes)
    connectivity = np.vstack((elem_connec1, elem_connec2)).T

    # read the mesh data structure
    mesh = Mesh1D(nodes, connectivity)
    # Harmonic Acoustic problem define the frequency

    k_0 = omega/Air.c
    theta = 50  #indidence angle
    ky = k_0*np.sin(theta*np.pi/180)
    nature_bcs = {'type': 'velocity', 'value': 1, 'position': 0}
    # right_hand_side = adm_assembler.assemble_nature_bc(nature_bcs)

    # ====================== Analytical Solution ======================
    # analytical solution
    # kundlt_tube = DoubleleLayerKundltTube(mesh, air, xfm, omega, nature_bcs)
    # ana_p = np.zeros(num_nodes, dtype=np.complex128)  #initialize the analytical solution vector
    # p_analy = kundlt_tube.sol_on_nodes(ana_p, sol_type='pressure')
    # ana_v = np.zeros(num_nodes, dtype=np.complex128)  #initialize the analytical solution vector
    # v_analy = kundlt_tube.sol_on_nodes(ana_v, sol_type='velocity')
    # import pdb; pdb.set_trace() 
    # Z_s = p_analy[0]/v_analy[0]
    # ref_analy = (Z_s - air.Z_f) / (Z_s + air.Z_f)
    # absop_analy = 1 - np.abs(ref_analy)**2
    thickness = 0.1
    tm_poro1 = TMMPoroElastic1(xfm, omega, k_0)
    tm_poro1_matrix = tm_poro1.transfer_matrix(thickness)
    # print(tm_poro1_matrix)
    tm_poro2 = TMMPoroElastic2(xfm, omega, k_0)
    tm_poro2_matrix = tm_poro2.transfer_matrix(thickness)
    # print(tm_poro2_matrix)

    diff = tm_poro1_matrix-tm_poro2_matrix
    diff[diff<1e-10] = 0

    # print("difference between two transfer matrix is:", diff)
    Tmm_a = TMMFluid(air, omega)
    T_a = Tmm_a.transfer_matrix(thickness)

    Tmm_xfm = TMMPoroElastic3(xfm, omega, ky)
    T_xfm = Tmm_xfm.transfer_matrix(thickness)

    I_fp, J_pf = bcm_fluid_poro2(xfm.phi)
    Y_p = bcm_rigid_wall2()

    D_dim = T_a.shape[1]+T_xfm.shape[1]
    D = np.zeros((D_dim,D_dim), dtype=np.complex128)
    D[1:I_fp.shape[0]+1, 0:I_fp.shape[1]] = I_fp
    D[1:I_fp.shape[0]+1, I_fp.shape[1]:] = J_pf@T_xfm
    D[I_fp.shape[0]+1:, I_fp.shape[1]:] = Y_p  # rigid wall
    D1 = D[1:,1:]
    D2 = D[1:,[0,2,3,4,5,6,7]]

    jomegaZs= -np.linalg.det(D1)/np.linalg.det(D2)
    Zs= jomegaZs/(1j*omega)
    # print(Zs)

    R= (Zs*np.cos(theta*np.pi/180)-air.Z_f)/(Zs*np.cos(theta*np.pi/180)+air.Z_f)
    # print(R)

    alpha = 1 - np.abs(R)**2


    return alpha


num_samp = int(1e5)
absops = np.zeros((1, num_samp))
freqs = np.linspace(10, 1e6, num=num_samp)

from mediapack import from_yaml
from pymls import Solver, Layer, backing
thickness_foam = 0.1
foam = from_yaml('tests/foam.yaml')
# instanciate the solver
# the layers are specified in order with the termination on the right
# Here:
#              |         |                         |\ 
#  Incident    |         |                         |\ 
#   Medium     | film    |      foam (2cm)         |\ rigid backing
#    (default) |  (.5mm) |                         |\ 
#              |         |                         |\ 
S = Solver(
    layers=[Layer(foam, thickness_foam)],
    backing=backing.rigid
)

result = S.solve(freqs.tolist(), 50)
# Result is a dict containing:
# - analysis name (str, default: auto)
# - enable_stochastic (bool, default: False)
# - f (frequency vector)
# - angle (angle vector)
# - R (list of complex reflexion coefficient)
# - T (list of complex transmission coefficient)
R = np.array(result['R'])
A = 1.-np.abs(R)**2
fig, axs = plt.subplots(nrows=1)
axs.plot(result['f'], A, label='pymls')



for i, freq in enumerate(freqs):
    absops[0,i] = discrete_in_frequency(freq)
    # plt.plot(freq, absop_analy, 'o')

# import pdb;pdb.set_trace()
if np.mean(np.abs(absops[0,:30] - A[:30])) < 1e-4:
    print("Test passed!")

# plt.plot(freqs, absops[0], '-', label='analytical')
plt.plot(freqs, absops[0], '--', label='tmm')
plt.ylim(0., 1.0)
plt.xscale('log')
plt.legend()
# plt.show()
plt.show(block=False)
plt.pause(1)
plt.close('all')
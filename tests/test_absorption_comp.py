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
working_dir = os.path.join(current_dir, "..")
import sys

sys.path.append(working_dir)

import numpy as np
import matplotlib.pyplot as plt

from SAcouS.acxtmm.adm_basis import AdmFluid
from SAcouS.acxtmm.tmm import TMMFluid, TMMPoroElastic1, TMMPoroElastic2
from SAcouS.acxtmm.adm_assembler import AdmAssembler
from SAcouS.acxfem.mesh import Mesh1D
from SAcouS.acxfem.materials import Air, Fluid, EquivalentFluid
from SAcouS.acxfem.utilities import check_material_compability, display_matrix_in_array, plot_matrix_partten
from SAcouS.acxfem.solver import AdmittanceSolver
from SAcouS.acxfem.postprocess import PostProcessField
from analytical.fluid_sol import DoubleleLayerKundltTube

# mesh.refine_mesh(1)

# ====================== Pysical Problem ======================
# define the materials
air = Air('classical air')

# given JCA porous material properties
phi = 0.4    # porosity
sigma = 4e6    # resistivity
alpha = 1.75    # Tortuosity
Lambda_prime = 9.3e-6    # Viscous characteristic length
Lambda = 2e-6    #
xfm = EquivalentFluid('xfm', phi, sigma, alpha, Lambda_prime, Lambda)


def discrete_in_frequency(freq):

  omega = 2 * np.pi * freq    # angular frequency
  lambda_ = air.c_f / freq    # wavelength
  h = lambda_ / 10    # size of elements
  num_elem = int(0.2 / h) + 10    # number of elements
  # import pdb; pdb.set_trace()
  num_nodes = num_elem + 1    # number of nodes
  nodes = np.linspace(-0.1, 0.1, num_nodes)

  elem_connec1 = np.arange(0, num_elem)
  elem_connec2 = np.arange(1, num_nodes)
  connectivity = np.vstack((elem_connec1, elem_connec2)).T

  # read the mesh data structure
  mesh = Mesh1D(nodes, connectivity)
  # Harmonic Acoustic problem define the frequency
  k_0 = omega / Air.c

  theta = 1e-8    #indidence angle

  # define the subdomains: domain name (material) and the elements in the domain
  air_elements = np.arange(0, num_elem / 2)
  xfm_elements = np.arange(num_elem / 2, num_nodes)
  subdomains = {air: air_elements, xfm: xfm_elements}

  # adm_assembler = AdmAssembler(mesh, subdomains, omega, dtype=np.complex128)
  # left_hand_side = adm_assembler.assemble_global_adm(theta, k_0, 'continue')

  nature_bcs = {'type': 'fluid_velocity', 'value': 1, 'position': 0}
  # right_hand_side = adm_assembler.assemble_nature_bc(nature_bcs)

  # adm_solver = AdmittanceSolver(left_hand_side, right_hand_side)
  # adm_solver.solve()
  # sol = adm_solver.sol

  # ====================== Analytical Solution ======================
  # analytical solution
  kundlt_tube = DoubleleLayerKundltTube(mesh, air, xfm, omega, nature_bcs)
  ana_p = np.zeros(
      num_nodes,
      dtype=np.complex128)    #initialize the analytical solution vector
  p_analy = kundlt_tube.sol_on_nodes(ana_p, sol_type='pressure')
  ana_v = np.zeros(
      num_nodes,
      dtype=np.complex128)    #initialize the analytical solution vector
  v_analy = kundlt_tube.sol_on_nodes(ana_v, sol_type='fluid_velocity')
  # import pdb; pdb.set_trace()
  Z_s = p_analy[0] / v_analy[0]
  ref_analy = (Z_s - air.Z_f) / (Z_s + air.Z_f)
  absop_analy = 1 - np.abs(ref_analy)**2

  tm_fluid = TMMFluid(xfm, omega)
  tm = tm_fluid.transfer_matrix(0.1)
  ref_tmm = (tm[0, 0] / (tm[1, 0] * air.Z_f) - 1) / (tm[0, 0] /
                                                     (tm[1, 0] * air.Z_f) + 1)
  absop_tmm = 1 - np.abs(ref_tmm)**2

  return absop_analy, absop_tmm


absops = np.zeros((2, 1000))
freqs = np.linspace(1, 1e5, 1000)
for i, freq in enumerate(freqs):
  absops[0, i], absops[1, i] = discrete_in_frequency(freq)
  # plt.plot(freq, absop_analy, 'o')
error = np.mean(np.abs(absops[0] - absops[1]))
if error < 1e-5:
  print("Test passed!")
print("error:", error)

plt.plot(freqs, absops[0], '-', label='analytical')
plt.plot(freqs, absops[1], '--', label='tmm')
plt.legend()
plt.show(block=False)
plt.pause(1)
plt.close('all')

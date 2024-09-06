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

from SAcouS.acxtmm.adm_assembler import AdmAssembler
from SAcouS.Mesh import Mesh1D
from SAcouS.Materials import Air, Fluid, EquivalentFluid
from SAcouS.acxfem import check_material_compability, display_matrix_in_array, plot_matrix_partten
from SAcouS.acxfem import AdmittanceSolver
from SAcouS.PostProcess import PostProcessField


def test_case_1():
  num_elem = 1000    # number of elements
  num_nodes = num_elem + 1    # number of nodes

  nodes = np.linspace(-1, 1e-2, num_nodes)

  elem_connec1 = np.arange(0, num_elem)
  elem_connec2 = np.arange(1, num_nodes)
  connectivity = np.vstack((elem_connec1, elem_connec2)).T
  # print(connectivity)

  # read the mesh data structure
  mesh = Mesh1D(nodes, connectivity)
  # mesh.refine_mesh(1)

  # ====================== Pysical Problem ======================
  # define the materials
  air = Air('classical air')
  air1 = Air('air gap')

  # given JCA porous material properties
  phi = 0.98    # porosity
  sigma = 3.75e3    # resistivity
  alpha = 1.17    # Tortuosity
  Lambda_prime = 742e-6    # Viscous characteristic length
  Lambda = 110e-6    #
  foam1 = EquivalentFluid('xfm', phi, sigma, alpha, Lambda_prime, Lambda)

  # given JCA porous material properties
  phi = 0.95    # porosity
  sigma = 2850    # resistivity
  alpha = 1.05    # Tortuosity
  Lambda_prime = 900e-6    # Viscous characteristic length
  Lambda = 300e-6    #
  foam2 = EquivalentFluid('xfm', phi, sigma, alpha, Lambda_prime, Lambda)

  # Harmonic Acoustic problem define the frequency
  freq = 2000
  omega = 2 * np.pi * freq    # angular frequency
  k_0 = omega / Air.c
  theta = 1e-8    #indidence angle

  # define the subdomains: domain name (material) and the elements in the domain
  air_elements = np.arange(0, num_elem / 2)
  foam1_elements = np.arange(num_elem / 2, 3 / 4 * num_elem)
  foam2_elements = np.arange(3 / 4 * num_elem, num_nodes)
  subdomains = {
      air: air_elements,
      foam1: foam1_elements,
      foam2: foam2_elements
  }
  check_material_compability(subdomains)

  adm_assembler = AdmAssembler(mesh, subdomains, omega, dtype=np.complex128)
  left_hand_side = adm_assembler.assemble_global_adm(theta, k_0, 'continue')

  nature_bcs = {
      'type': 'fluid_velocity',
      'value': 1 * np.exp(-1j * omega),
      'position': 0
  }
  right_hand_side = adm_assembler.assemble_nature_bc(nature_bcs)

  adm_solver = AdmittanceSolver(left_hand_side, right_hand_side)
  adm_solver.solve()
  sol = adm_solver.sol

  # # plot the solution
  post_processer = PostProcessField(mesh.nodes, r'1D Helmholtz (2000$Hz$)')
  post_processer.plot_sol((np.real(sol), f'FADM', 'solid'))

  plt.show(block=False)
  plt.pause(1)
  plt.close('all')


if __name__ == "__main__":
  result = test_case_1()

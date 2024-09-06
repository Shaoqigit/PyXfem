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
# import sys
# sys.path.append('/home/shaoqi/Devlop/PyXfem/PyAcoustiX/')
import os

current_dir = os.path.dirname(os.path.realpath(__file__))
working_dir = os.path.join(current_dir, "..")
import sys

sys.path.append(working_dir)

import numpy as np
import matplotlib.pyplot as plt

from SAcouS.Mesh import Mesh1D
from SAcouS.Materials import Air, Fluid, EquivalentFluid
from SAcouS.PostProcess import PostProcessField

from SAcouS.acxfem import Helmholtz1DElement
from SAcouS.acxfem import DofHandler1D, GeneralDofHandler1D, FESpace
from SAcouS.acxfem import HelmholtzAssembler, BiotAssembler, CouplingAssember
from SAcouS.acxfem import ApplyBoundaryConditions
from SAcouS.acxfem import check_material_compability, display_matrix_in_array, plot_matrix_partten
from SAcouS.acxfem import LinearSolver
from analytical.fluid_sol import DoubleleLayerKundltTube


def test_case_1():
  # ====================== Pysical Problem ======================
  # define the materials
  air = Air('classical air')

  # given JCA porous material properties
  phi = 0.98    # porosity
  sigma = 3.75e3    # resistivity
  alpha = 1.17    # Tortuosity
  Lambda_prime = 742e-6    # Viscous characteristic length
  Lambda = 110e-6    #
  xfm = EquivalentFluid('xfm', phi, sigma, alpha, Lambda_prime, Lambda)

  # Harmonic Acoustic problem define the frequency
  freq = 2000
  omega = 2 * np.pi * freq    # angular frequency
  xfm.set_frequency(omega)
  # ====================== Mesh and basis definition ======================
  num_elem = 1000    # number of elements
  num_nodes = num_elem + 1    # number of nodes
  nodes = np.linspace(-1, 1, num_nodes)
  elem_connec1 = np.arange(0, num_elem)
  elem_connec2 = np.arange(1, num_nodes)
  connectivity = np.vstack((elem_connec1, elem_connec2)).T
  # read the mesh data structure
  mesh = Mesh1D(nodes, connectivity)
  # mesh.plotmesh(withnode=True)
  # define the subdomains: domain name (material) and the elements in the domain
  air_elements = np.arange(0, int(num_elem / 2))
  xfm_elements = np.arange(int(num_elem / 2), num_elem)
  subdomains = {air: air_elements, xfm: xfm_elements}
  check_material_compability(subdomains)
  elements2node = mesh.get_mesh(
  )    # dict: elements number with nodes coodinates
  # print(elements_set)
  order = 2    # global order of the bases
  # applied the basis on each element
  Pf_bases = []
  for mat, elems in subdomains.items():
    if mat.TYPE == 'Fluid':
      Pf_bases += [
          Helmholtz1DElement('Pf', order, elements2node[elem],
                             (1 / mat.rho_f, 1 / mat.K_f)) for elem in elems
      ]
    # print(basis.ke)
    # print(basis.me)

  # handler the dofs: map the basis to mesh
  fe_space = FESpace(mesh, subdomains, Pf_bases)

  # initialize the assembler
  Helmholtz_assember = HelmholtzAssembler(fe_space, dtype=np.complex128)
  Helmholtz_assember.assembly_global_matrix(Pf_bases, 'Pf')
  left_hand_matrix = Helmholtz_assember.get_global_matrix(omega)

  right_hand_vec = np.zeros(Helmholtz_assember.nb_global_dofs,
                            dtype=np.complex128)
  #  natural boundary condition
  nature_bcs = {
      'type': 'fluid_velocity',
      'value': 1 * np.exp(-1j * omega),
      'position': -1.0
  }
  BCs_applier = ApplyBoundaryConditions(mesh, fe_space, left_hand_matrix,
                                        right_hand_vec, omega)
  BCs_applier.apply_nature_bc(nature_bcs, var='Pf')

  # solver the linear system
  linear_solver = LinearSolver(fe_space=fe_space)
  linear_solver.solve(left_hand_matrix, right_hand_vec)
  sol = linear_solver.u

  # ====================== Analytical Solution ======================
  # analytical solution
  kundlt_tube = DoubleleLayerKundltTube(1, 1, air, xfm, omega, nature_bcs)
  ana_sol = kundlt_tube.sol_on_mesh(mesh, sol_type='pressure')

  # plot the solution
  post_processer = PostProcessField(mesh.nodes, r'1D Helmholtz (2000$Hz$)')
  post_processer.plot_sol((np.real(sol), f'FEM ($p=3$)', 'solid'),
                          (np.real(ana_sol), 'Analytical', 'dashed'))
  plt.show(block=False)
  plt.pause(1)
  plt.close('all')

  error = post_processer.compute_error(sol, ana_sol)
  print("error:", error)
  if error < 1e-5:
    print("Test passed!")
    return True
  else:
    print("Test failed!")
    return False


if __name__ == "__main__":
  result = test_case_1()

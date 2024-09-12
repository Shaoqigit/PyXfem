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

from SAcouS.Materials import Air, Fluid, EquivalentFluid
from SAcouS.Mesh import Mesh2D, Mesh1D, MeshReader
from SAcouS.PostProcess import plot_field, save_plot, PostProcessField, read_solution

from SAcouS.acxfem import Helmholtz3DElement
from SAcouS.acxfem import FESpace
from SAcouS.acxfem import check_material_compability, display_matrix_in_array, plot_matrix_partten
from SAcouS.acxfem import HelmholtzAssembler
from SAcouS.acxfem import LinearSolver
from SAcouS.acxfem import ApplyBoundaryConditions
from analytical.fluid_sol import DoubleleLayerKundltTube


def test_case_2D():
  # ====================== Pysical Problem ======================
  # define the materials
  air = Air('classical air')
  freq = 10
  omega = 2 * np.pi * freq    # angular frequency
  # Harmonic Acoustic problem define the frequency
  current_dir = os.path.dirname(os.path.realpath(__file__))
  mesh_reader = MeshReader(current_dir + "/mesh/plane.msh", dim=3)
  mesh = mesh_reader.get_mesh()

  air_elements = np.arange(0, mesh.nb_elmes)
  elements2node = mesh.get_mesh()
  subdomains = {air: air_elements}
  Pf_bases = []
  order = 1
  for mat, elems in subdomains.items():
    if mat.TYPE == 'Fluid':
      Pf_bases += [
          Helmholtz3DElement('Pf', order, elements2node[elem],
                             (1 / mat.rho_f, 1 / mat.K_f)) for elem in elems
      ]
  # handler the dofs: map the basis to mesh
  fe_space = FESpace(mesh, subdomains, Pf_bases)
  # initialize the assembler
  import time
  start_time = time.time()
  Helmholtz_assember = HelmholtzAssembler(fe_space, dtype=float)
  Helmholtz_assember.assembly_global_matrix(Pf_bases, 'Pf')
  left_hand_matrix = Helmholtz_assember.get_global_matrix(omega)
  print("Time taken to assemble the matrix:", time.time() - start_time)
  right_hand_vec = np.zeros(Helmholtz_assember.nb_global_dofs,
                            dtype=np.complex128)

  # ====================== Boundary Conditions ======================
  # natural_edge = np.arange(64, 85)
  natural_facet = mesh_reader.get_facet_by_physical('inlet')
  natural_bcs = {
      'type':
      'fluid_velocity',
      'value':
      lambda x, y, z: np.array(
          [1e12 * np.exp(-1j * omega), 0, 1e12 * np.exp(-1j * omega)]),
      'position':
      natural_facet
  }    # position: number of facet number

  right_hand_vec = np.zeros(Helmholtz_assember.nb_global_dofs,
                            dtype=np.complex128)
  BCs_applier = ApplyBoundaryConditions(mesh, fe_space, left_hand_matrix,
                                        right_hand_vec, omega)
  BCs_applier.apply_nature_bc(natural_bcs, 'Pf', 13)
  linear_solver = LinearSolver(fe_space=fe_space)
  linear_solver.solve(left_hand_matrix, right_hand_vec, 'petsc')
  sol = linear_solver.u
  save_plot(mesh,
            np.abs(sol),
            'Pf_numerical',
            current_dir + "/Pf_3D_airplane.pos",
            engine='gmsh',
            binary=True)
  # #
  # num_elem = 200    # number of elements
  # num_nodes = num_elem + 1    # number of nodes

  # nodes = np.linspace(-0.5, 0.5, num_nodes)

  # elem_connec1 = np.arange(0, num_elem)
  # elem_connec2 = np.arange(1, num_nodes)
  # connectivity = np.vstack((elem_connec1, elem_connec2)).T
  # # print(connectivity)

  # # read the mesh data structure
  # mesh = Mesh1D(nodes, connectivity)
  # ana_sol2 = kundlt_tube.sol_on_mesh(mesh, sol_type='pressure')

  # # plot the solution
  # post_processer = PostProcessField(mesh.nodes, r'1D Helmholtz (2000$Hz$)')
  # post_processer.plot_sol((np.real(ana_sol2), 'Analytical', 'solid'))
  # plt.show()


if __name__ == "__main__":
  import time
  start = time.time()
  result = test_case_2D()
  print("Time taken:", time.time() - start)

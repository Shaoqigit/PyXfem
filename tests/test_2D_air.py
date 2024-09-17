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

from SAcouS.acxfem import Helmholtz2DElement
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
  freq = 1000
  omega = 2 * np.pi * freq    # angular frequency
  slice_points_1 = np.insert(np.arange(402, 797)[::-1], 0, 3)
  slice_points = np.append(slice_points_1, 2)
  # Harmonic Acoustic problem define the frequency
  current_dir = os.path.dirname(os.path.realpath(__file__))
  mesh_reader = MeshReader(current_dir + "/mesh/unit_tube_2.msh")
  mesh = mesh_reader.get_mesh()

  air_elements = np.arange(0, mesh.nb_elmes)
  elements2node = mesh.get_mesh_coordinates()
  mesh.set_subdomains({air: air_elements})
  Pf_bases = []
  order = 1
  for mat, elems in mesh.subdomains.items():
    if mat.TYPE == 'Fluid':
      Pf_bases += [
          Helmholtz2DElement('Pf', order, elements2node[elem],
                             (1 / mat.rho_f, 1 / mat.K_f)) for elem in elems
      ]
  # handler the dofs: map the basis to mesh
  fe_space = FESpace(mesh, Pf_bases)
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
  natural_edge = np.arange(796, 800)
  natural_bcs = {
      'type': 'fluid_velocity',
      'value': lambda x, y: np.array([1 * np.exp(-1j * omega), 0]),
      'position': natural_edge
  }    # position: number of facet number

  right_hand_vec = np.zeros(Helmholtz_assember.nb_global_dofs,
                            dtype=np.complex128)
  BCs_applier = ApplyBoundaryConditions(mesh, fe_space, left_hand_matrix,
                                        right_hand_vec, omega)
  BCs_applier.apply_nature_bc(natural_bcs, 'Pf')
  linear_solver = LinearSolver(fe_space=fe_space)
  linear_solver.solve(left_hand_matrix, right_hand_vec, 'petsc')
  sol = linear_solver.u
  mesh_io = mesh_reader.meshio_object
  save_plot(mesh_io,
            sol.real,
            'Numeral_solution',
            current_dir + "/Pressure_field.pos",
            engine='gmsh',
            binary=True)
  nodes = mesh.nodes[slice_points][:, 0]

  sol = sol[slice_points]

  elem_connec1 = np.arange(0, 396)
  elem_connec2 = np.arange(1, 397)
  connectivity = np.vstack((elem_connec1, elem_connec2)).T
  mesh_1d = Mesh1D(nodes, connectivity)
  natural_bcs_ana = {
      'type': 'fluid_velocity',
      'value': np.exp(-1j * omega),
      'position': -0.5
  }
  kundlt_tube = DoubleleLayerKundltTube(0.5, 0.5, air, air, omega,
                                        natural_bcs_ana)
  ana_sol = kundlt_tube.sol_on_mesh(mesh_1d, sol_type='pressure')

  post_processer = PostProcessField(mesh_1d.nodes, r'2D Helmholtz (2000$Hz$)')
  post_processer.plot_sol((np.real(sol), f'FEM ($p=3$)', 'solid'),
                          (np.real(ana_sol), 'Analytical', 'dashed'))
  # plt.show()
  error = np.mean(np.real(sol) - np.real(ana_sol)) / np.mean(np.real(ana_sol))
  print("error:", error)
  if error < 0.0005:
    print("Test passed!")
    return True
  else:
    print("Test failed!")
    return False


if __name__ == "__main__":
  import time
  start = time.time()
  result = test_case_2D()
  print("Time taken:", time.time() - start)

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

from SAcouS.Mesh import Mesh1D, MeshReader
from SAcouS.Materials import Air, EquivalentFluid
from SAcouS.PostProcess import plot_field, save_plot, PostProcessField, read_solution

from SAcouS.acxfem import Helmholtz3DElement
from SAcouS.acxfem import FESpace
from SAcouS.acxfem import check_material_compability, display_matrix_in_array, plot_matrix_partten
from SAcouS.acxfem import HelmholtzAssembler
from SAcouS.acxfem import LinearSolver
from SAcouS.acxfem import ApplyBoundaryConditions
from analytical.fluid_sol import DoubleleLayerKundltTube


def test_case_3D():
  # ====================== Pysical Problem ======================
  # define the materials
  air = Air('classical air')
  # given JCA porous material properties
  phi = 0.98    # porosity
  sigma = 3.75e3    # resistivity
  alpha = 1.17    # Tortuosity
  Lambda_prime = 742e-6    # Viscous characteristic length
  Lambda = 110e-6    #
  xfm = EquivalentFluid('foam', phi, sigma, alpha, Lambda_prime, Lambda)

  freq = 1000
  omega = 2 * np.pi * freq    # angular frequency
  # Harmonic Acoustic problem define the frequency
  current_dir = os.path.dirname(os.path.realpath(__file__))
  mesh_reader = MeshReader(current_dir + "/mesh/tube_3D_2f_refine.msh", dim=3)
  mesh = mesh_reader.get_mesh()
  air_elements = mesh_reader.get_elem_by_physical('air')
  foam_elements = mesh_reader.get_elem_by_physical('foam')

  elements2node = mesh.get_mesh()
  xfm.set_frequency(omega)
  subdomains = {air: air_elements, xfm: foam_elements}
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
  Helmholtz_assember = HelmholtzAssembler(fe_space, dtype=np.complex64)
  Helmholtz_assember.assembly_global_matrix(Pf_bases, 'Pf')
  left_hand_matrix = Helmholtz_assember.get_global_matrix(omega)
  print("Time taken to assemble the matrix:", time.time() - start_time)
  right_hand_vec = np.zeros(Helmholtz_assember.nb_global_dofs,
                            dtype=np.complex128)

  # ====================== Boundary Conditions ======================
  natural_facet = mesh_reader.get_facet_by_physical('inlet')
  natural_bcs = {
      'type': 'fluid_velocity',
      'value': lambda x, y, z: np.array([1 * np.exp(-1j * omega), 0, 0]),
      'position': natural_facet
  }    # position: number of facet number

  BCs_applier = ApplyBoundaryConditions(mesh, fe_space, left_hand_matrix,
                                        right_hand_vec, omega)
  BCs_applier.apply_nature_bc(natural_bcs, 'Pf')
  linear_solver = LinearSolver(fe_space=fe_space)
  linear_solver.solve(left_hand_matrix, right_hand_vec, solver='petsc')
  sol = linear_solver.u
  save_plot(mesh,
            sol.real,
            'Numeral_solution',
            current_dir + "/Pf_two_fluid_fem.pos",
            engine='gmsh',
            binary=True)
  # nodes = mesh.nodes[slice_points][:, 0]

  natural_bcs_ana = {
      'type': 'fluid_velocity',
      'value': np.exp(-1j * omega),
      'position': -0.5
  }
  kundlt_tube = DoubleleLayerKundltTube(0.5, 0.5, air, xfm, omega,
                                        natural_bcs_ana)
  ana_sol = kundlt_tube.sol_on_mesh(mesh, sol_type='pressure')
  save_plot(mesh,
            ana_sol.real,
            'Analytical_solution',
            current_dir + "/Pf_two_fluid_analy.pos",
            engine='gmsh',
            binary=True)
  error = np.mean(np.abs(sol - ana_sol)) / np.mean(np.abs(ana_sol))
  print("error:", error)
  if error < 0.02:
    print("Test passed!")
    return True
  else:
    print("Test failed!")
    return False


if __name__ == "__main__":
  import time
  start = time.time()
  result = test_case_3D()
  print("Time taken:", time.time() - start)

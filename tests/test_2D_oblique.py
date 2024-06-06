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

from SAcouS.acxfem.basis import Helmholtz2DElement
from SAcouS.acxfem.mesh import Mesh2D, MeshReader
from SAcouS.acxfem.dofhandler import DofHandler1D, GeneralDofHandler1D, FESpace
from SAcouS.acxfem.materials import Air, Fluid, EquivalentFluid
from SAcouS.acxfem.utilities import check_material_compability, display_matrix_in_array, plot_matrix_partten
from SAcouS.acxfem.physic_assembler import HelmholtzAssembler
from SAcouS.acxfem.solver import LinearSolver
from SAcouS.acxfem.postprocess import plot_field, save_plot, PostProcessField, read_solution
from SAcouS.acxfem.BCs_impose import ApplyBoundaryConditions
from analytical.fluid_sol import ObliquePlaneWave


def test_case_2D():
  # ====================== Pysical Problem ======================
  # define the materials
  air = Air('classical air')
  phi = 0.98    # porosity
  sigma = 3.75e3    # resistivity
  alpha = 1.17    # Tortuosity
  Lambda_prime = 742e-6    # Viscous characteristic length
  Lambda = 110e-6    #
  xfm = EquivalentFluid('foam', phi, sigma, alpha, Lambda_prime, Lambda)
  freq = 1000
  omega = 2 * np.pi * freq    # angular frequency

  # Harmonic Acoustic problem define the frequency
  import time
  start_time = time.time()
  current_dir = os.path.dirname(os.path.realpath(__file__))
  mesh_reader = MeshReader(current_dir + "/mesh/mat2_oblique_utra_fine.msh")
  # mesh_reader = MeshReader(current_dir + "/mesh/half_tube_2.msh")
  mesh = mesh_reader.get_mesh()

  # ====================== Analytical solution ======================
  analytical_solution = ObliquePlaneWave(mesh, air, xfm, omega, 45, 1.)
  nb_nodes = mesh.get_nb_nodes()
  analytical_solution_vec = np.zeros((nb_nodes), dtype=np.complex128)
  analytical_solution.sol_on_nodes(analytical_solution_vec)
  save_plot(mesh,
            analytical_solution_vec.real,
            current_dir + "/Pressure_field_oblique_ana.pos",
            engine='gmsh',
            binary=True)

  air_elements = mesh_reader.get_elem_by_physical('mat1')
  foam_elements = mesh_reader.get_elem_by_physical('mat2')
  left_boundary = mesh_reader.get_edge_by_physical('left')
  right_boundary = mesh_reader.get_edge_by_physical('right')
  left_top_boundary = mesh_reader.get_edge_by_physical('left_top')
  right_top_boundary = mesh_reader.get_edge_by_physical('right_top')
  left_bottom_boundary = mesh_reader.get_edge_by_physical('left_bot')
  right_bottom_boundary = mesh_reader.get_edge_by_physical('right_bot')
  elements2node = mesh.get_mesh()
  # breakpoint()
  # subdomains = {air: all_elements}
  subdomains = {air: air_elements, xfm: foam_elements}
  xfm.set_frequency(omega)
  Pf_bases = []
  order = 1
  for mat, elems in subdomains.items():
    if mat.TYPE == 'Fluid':
      Pf_bases += [
          Helmholtz2DElement('Pf', order, elements2node[elem],
                             (1 / mat.rho_f, 1 / mat.K_f)) for elem in elems
      ]
  fe_space = FESpace(mesh, subdomains, Pf_bases)
  print("Number of global dofs:", fe_space.get_nb_dofs())
  # initialize the assembler
  start_assembly_time = time.time()
  Helmholtz_assember = HelmholtzAssembler(fe_space, dtype=np.complex64)
  Helmholtz_assember.assembly_global_matrix(Pf_bases, 'Pf')
  left_hand_matrix = Helmholtz_assember.get_global_matrix(omega)
  print("Time taken to assemble the matrix:",
        time.time() - start_assembly_time)
  right_hand_vec = np.zeros(Helmholtz_assember.nb_global_dofs,
                            dtype=np.complex128)

  v1, v2 = analytical_solution.velocity_1, analytical_solution.velocity_2    #v_x(x<0), v_x
  BCs_applier = ApplyBoundaryConditions(mesh, fe_space, left_hand_matrix,
                                        right_hand_vec, omega)
  # (x>0), v_y(x<0), v_y(x>0)
  natural_bcs1 = {
      'type': 'analytical_gradient',
      'value': v1,
      'position': left_boundary
  }    # position

  BCs_applier.apply_nature_bc(natural_bcs1, 'Pf')

  natural_bcs2 = {
      'type': 'analytical_gradient',
      'value': v1,
      'position': left_top_boundary
  }
  # breakpoint()
  BCs_applier.apply_nature_bc(natural_bcs2, 'Pf')

  natural_bcs3 = {
      'type': 'analytical_gradient',
      'value': v1,
      'position': left_bottom_boundary
  }
  # breakpoint()
  BCs_applier.apply_nature_bc(natural_bcs3, 'Pf')

  natural_bcs4 = {
      'type': 'analytical_gradient',
      'value': v2,
      'position': right_boundary
  }
  BCs_applier.apply_nature_bc(natural_bcs4, 'Pf')

  natural_bcs5 = {
      'type': 'analytical_gradient',
      'value': v2,
      'position': right_top_boundary
  }
  BCs_applier.apply_nature_bc(natural_bcs5, 'Pf')

  natural_bcs6 = {
      'type': 'analytical_gradient',
      'value': v2,
      'position': right_bottom_boundary
  }
  # breakpoint()
  BCs_applier.apply_nature_bc(natural_bcs6, 'Pf')

  linear_solver = LinearSolver(fe_space=fe_space)
  linear_solver.solve(left_hand_matrix, right_hand_vec)
  sol = linear_solver.u

  print("Time taken of FEM process:", time.time() - start_time)

  save_plot(mesh,
            sol.real,
            current_dir + "/Pressure_field_oblique_succes.pos",
            engine='gmsh',
            binary=True)

  error = np.mean(np.abs(sol - analytical_solution_vec)) / np.mean(
      np.abs(analytical_solution_vec))
  print("error:", error)
  if error < 0.05 and (time.time() - start_time) < 80:
    print("Test passed!")
    return True
  else:
    print("Test failed!")
    return False


if __name__ == "__main__":
  # Python
  import cProfile

  # Profile the main function and save the results to 'profile_results.prof'
  #   cProfile.run('test_case_2D()', 'profile_results.prof')
  test_case_2D()
  # Python
  import pstats

  #   Create a pstats.Stats object from the profiling results
#   stats = pstats.Stats('profile_results.prof')

#   # Sort the statistics by the cumulative time spent in the function
#   stats.sort_stats('cumulative')

#   # Print the statistics
#   stats.print_stats()

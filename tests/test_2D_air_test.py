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

from SAcouS.acxfem.basis import Lagrange2DTriElement
from SAcouS.acxfem.mesh import Mesh2D, Mesh1D
from SAcouS.acxfem.dofhandler import DofHandler1D, GeneralDofHandler1D, FESpace
from SAcouS.acxfem.materials import Air, Fluid, EquivalentFluid
from SAcouS.acxfem.utilities import check_material_compability, display_matrix_in_array, plot_matrix_partten
from SAcouS.acxfem.physic_assembler import HelmholtzAssembler
from SAcouS.acxfem.solver import LinearSolver
from SAcouS.acxfem.postprocess import plot_field, save_plot, PostProcessField, read_solution
from SAcouS.acxfem.BCs_impose import ApplyBoundaryConditions
from analytical.fluid_sol import DoubleleLayerKundltTube


def test_case_2D():
  # ====================== Pysical Problem ======================
  # define the materials
  air = Air('classical air')
  freq = 1000
  omega = 2 * np.pi * freq    # angular frequency
  slice_points_1 = np.insert(np.arange(402, 797)[::-1], 0, 3)
  slice_points = np.append(slice_points_1, 2)
  if not os.path.exists("Pressure_field.pos"):
    # Harmonic Acoustic problem define the frequency
    current_dir = os.path.dirname(os.path.realpath(__file__))
    mesh = Mesh2D()
    mesh.read_mesh(current_dir + "/mesh/unit_tube_2.msh")
    # mesh.plotmesh(withedgeid=True)
    # import pdb
    # pdb.set_trace()
    air_elements = np.arange(0, mesh.nb_elmes)
    elements2node = mesh.get_mesh()
    subdomains = {air: air_elements}
    Pf_bases = []
    order = 1
    for mat, elems in subdomains.items():
      if mat.TYPE == 'Fluid':
        Pf_bases += [
            Lagrange2DTriElement('Pf', order, elements2node[elem])
            for elem in elems
        ]
    # handler the dofs: map the basis to mesh
    fe_space = FESpace(mesh, subdomains, Pf_bases)
    # initialize the assembler
    Helmholtz_assember = HelmholtzAssembler(fe_space, subdomains, dtype=float)
    Helmholtz_assember.assembly_global_matrix(Pf_bases, 'Pf', omega)
    left_hand_matrix = Helmholtz_assember.get_global_matrix()
    right_hand_vec = np.zeros(Helmholtz_assember.nb_global_dofs,
                              dtype=np.complex128)

    # ====================== Boundary Conditions ======================
    # natural_edge = np.arange(64, 85)
    natural_edge = np.arange(797, 801)
    natural_bcs = {
        'type': 'fluid_velocity',
        'value': lambda x, y: 1 * np.exp(-1j * omega),
        'position': natural_edge
    }    # position: number of facet number

    right_hand_vec = np.zeros(Helmholtz_assember.nb_global_dofs,
                              dtype=np.complex128)
    BCs_applier = ApplyBoundaryConditions(mesh, fe_space, left_hand_matrix,
                                          right_hand_vec, omega)
    BCs_applier.apply_nature_bc(natural_bcs, 'Pf')
    linear_solver = LinearSolver(fe_space=fe_space)
    linear_solver.solve(left_hand_matrix, right_hand_vec)
    sol = linear_solver.u
    save_plot(mesh,
              sol.real,
              current_dir + "/Pressure_field.pos",
              engine='gmsh',
              binary=True)
    nodes = mesh.nodes[slice_points][:, 0]
  else:
    mesh, sol = read_solution("Pressure_field.pos",
                              read_mesh=True,
                              engine='gmsh')
    nodes = mesh.points[slice_points][:, 0]

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
  kundlt_tube = DoubleleLayerKundltTube(mesh_1d, air, air, omega,
                                        natural_bcs_ana)
  ana_sol = np.zeros(
      397, dtype=np.complex128)    #initialize the analytical solution vector
  kundlt_tube.sol_on_nodes(ana_sol, sol_type='pressure')

  post_processer = PostProcessField(mesh_1d.nodes, r'2D Helmholtz (2000$Hz$)')
  post_processer.plot_sol((np.real(sol), f'FEM ($p=3$)', 'solid'),
                          (np.real(ana_sol), 'Analytical', 'dashed'))
  # plt.show()
  error = np.mean(sol - np.real(ana_sol)) / np.mean(np.real(ana_sol))
  print("error:", error)
  if error < 0.0005:
    print("Test passed!")
    return True
  else:
    print("Test failed!")
    return False


if __name__ == "__main__":
  # Python
  import cProfile

  # Profile the main function and save the results to 'profile_results.prof'
  # cProfile.run('test_case_2D()', 'profile_results.prof')
  # Python
  import pstats

  # Create a pstats.Stats object from the profiling results
  stats = pstats.Stats('profile_results.prof')

  # Sort the statistics by the cumulative time spent in the function
  stats.sort_stats('cumulative')

  # Print the statistics
  stats.print_stats()

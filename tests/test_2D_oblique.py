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
  current_dir = os.path.dirname(os.path.realpath(__file__))
  mesh_reader = MeshReader(current_dir + "/mesh/mat2_oblique_fine.msh")
  mesh = mesh_reader.get_mesh()
  # air_elements = mesh_reader.get_elem_by_physical('mat1')
  # foam_elements = mesh_reader.get_elem_by_physical('mat2')
  # elements2node = mesh.get_mesh()
  # subdomains = {air: air_elements, xfm: foam_elements}
  # Pf_bases = []
  # order = 1
  # for mat, elems in subdomains.items():
  #   if mat.TYPE == 'Fluid':
  #     Pf_bases += [
  #         Lagrange2DTriElement('Pf', order, elements2node[elem])
  #         for elem in elems
  #     ]
  # handler the dofs: map the basis to mesh

  # ====================== Analytical solution ======================
  analytical_solution = ObliquePlaneWave(mesh, air, xfm, omega, 45, 1.)
  nb_nodes = mesh.get_nb_nodes()
  analytical_solution_vec = np.zeros(nb_nodes)
  analytical_solution.sol_on_nodes(analytical_solution_vec)
  save_plot(mesh,
            analytical_solution_vec.real,
            current_dir + "/oblique_two_fluid.pos",
            engine='gmsh',
            binary=True)


if __name__ == "__main__":
  # Python
  import cProfile

  # Profile the main function and save the results to 'profile_results.prof'
  cProfile.run('test_case_2D()', 'profile_results.prof')
  # Python
  import pstats

  # Create a pstats.Stats object from the profiling results
  # stats = pstats.Stats('profile_results.prof')

  # Sort the statistics by the cumulative time spent in the function
  # stats.sort_stats('cumulative')

  # Print the statistics
  # stats.print_stats()

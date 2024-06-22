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

from SAcouS.Mesh import MeshReader
from SAcouS.Materials import Air, Fluid, EquivalentFluid
from SAcouS.PostProcess import plot_field, save_plot, PostProcessField, read_solution

from SAcouS.acxfem import FESpace
from SAcouS.acxfem import HelmholtzAssembler
from SAcouS.acxfem import Helmholtz2DElement
from SAcouS.acxfem import LinearSolver
from SAcouS.acxfem import ApplyBoundaryConditions
from analytical.fluid_sol import DoubleleLayerKundltTube


def test_case_2D():
  # ====================== Pysical Problem ======================
  # define the materials
  air = Air('classical air')
  freq = 1500
  omega = 2 * np.pi * freq    # angular frequency
  # Harmonic Acoustic problem define the frequency
  current_dir = os.path.dirname(os.path.realpath(__file__))
  mesh_reader = MeshReader(current_dir + "/mesh/square_air_imp.msh")
  mesh = mesh_reader.get_mesh()

  elements2node = mesh.get_mesh()
  air_elements = mesh_reader.get_elem_by_physical('air')
  imp_boundary = mesh_reader.get_edge_by_physical('impedance')
  subdomains = {air: air_elements}
  Pf_bases = []
  order = 1
  for mat, elems in subdomains.items():
    if mat.TYPE == 'Fluid':
      Pf_bases += [
          Helmholtz2DElement('Pf', order, elements2node[elem],
                             (1 / mat.rho_f, 1 / mat.K_f)) for elem in elems
      ]
  # handler the dofs: map the basis to mesh
  fe_space = FESpace(mesh, subdomains, Pf_bases)
  # initialize the assembler
  import time
  start_time = time.time()
  Helmholtz_assember = HelmholtzAssembler(fe_space, dtype=np.complex128)
  Helmholtz_assember.assembly_global_matrix(Pf_bases, 'Pf')
  left_hand_matrix = Helmholtz_assember.get_global_matrix(omega)
  print("Time taken to assemble the matrix:", time.time() - start_time)
  right_hand_vec = np.zeros((Helmholtz_assember.nb_global_dofs, 1),
                            dtype=np.complex128)

  # ====================== Boundary Conditions ======================
  # natural_edge = np.arange(64, 85)
  #source center
  Sx = 0.5
  Sy = 0.5
  source_amp = 0.001
  #Narrow normalized gauss distribution
  alfa = 0.015
  delta = lambda x, y: 1 / (np.abs(alfa) * np.sqrt(2 * np.pi)) * np.exp(-((
      (x - Sx)**2 + (y - Sy)**2) / (2 * alfa**2)))
  impedance_bcs = {
      'type':
      'impedance',
      'value':
      lambda x, y: air.rho_f * air.c_f /
      (1 + 1 / (1j * omega / air.c_f * np.sqrt((x - Sx)**2 + (y - Sy)**2))),
      'position':
      imp_boundary
  }    # position: number of facet number
  point_source = {
      'type': 'source',
      'value': lambda x, y: 1j * source_amp * delta(x, y) * omega * air.rho_f,
      'position': [Sx, Sy]
  }

  BCs_applier = ApplyBoundaryConditions(mesh, fe_space, left_hand_matrix,
                                        right_hand_vec, omega)
  left_hand_matrix = BCs_applier.apply_impedance_bc(impedance_bcs, 'Pf')
  right_hand_vec = BCs_applier.apply_source(point_source, Pf_bases, 'Pf')

  linear_solver = LinearSolver(fe_space=fe_space)
  linear_solver.solve(left_hand_matrix, right_hand_vec)
  sol = linear_solver.u
  save_plot(mesh,
            np.real(sol),
            current_dir + "/Pressure_field_real.pos",
            engine='gmsh',
            binary=True)
  save_plot(mesh,
            np.imag(sol),
            current_dir + "/Pressure_field_imag.pos",
            engine='gmsh',
            binary=True)


if __name__ == "__main__":
  import time
  start = time.time()
  result = test_case_2D()
  print("Time taken:", time.time() - start)

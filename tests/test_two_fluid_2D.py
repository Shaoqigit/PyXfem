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

from SAcouS.acxfem.basis import Lagrange2DTriElement
from SAcouS.acxfem.mesh import Mesh2D
from SAcouS.acxfem.dofhandler import DofHandler1D, GeneralDofHandler1D, FESpace
from SAcouS.acxfem.materials import Air, Fluid, EquivalentFluid
from SAcouS.acxfem.utilities import check_material_compability, display_matrix_in_array, plot_matrix_partten
from SAcouS.acxfem.physic_assembler import HelmholtzAssembler
from SAcouS.acxfem.solver import LinearSolver
from SAcouS.acxfem.postprocess import PostProcessField
from analytical.fluid_sol import ImpedenceKundltTube


def test_case_2D():
  # ====================== Pysical Problem ======================
  # define the materials
  air = Air('classical air')
  # Harmonic Acoustic problem define the frequency
  freq = 1000
  omega = 2 * np.pi * freq    # angular frequency

  mesh = Mesh2D()
  mesh.read_mesh("mesh/square_1.msh")
  # mesh.plotmesh()
  air_elements = np.arange(0, mesh.nb_elmes)
  elements2node = mesh.get_mesh()
  subdomains = {air: air_elements}
  import pdb
  pdb.set_trace()
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
  import pdb
  pdb.set_trace()
  # initialize the assembler
  Helmholtz_assember = HelmholtzAssembler(fe_space,
                                          subdomains,
                                          dtype=np.complex128)

  Helmholtz_assember.assembly_global_matrix(Pf_bases, 'Pf', omega)
  left_hand_matrix = Helmholtz_assember.get_global_matrix()


if __name__ == "__main__":
  result = test_case_2D()
